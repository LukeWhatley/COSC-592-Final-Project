import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd


GSE_REGEX = re.compile(r"(GSE\d+)", re.IGNORECASE)



@dataclass
class VitisSampleInfo:
    """Metadata about a single expression sample/column."""
    sample_id: str           # unique ID we create (tissue|GSE|colname)
    original_name: str       # original column name in the file
    tissue: str              # 'leaf', 'berry', etc.
    condition: str           # 'control', 'treatment', etc.
    gse_accession: str       # 'GSE97900', etc.
    cultivar: Optional[str]  # parsed from filename if possible
    file_path: str           # path to the source .txt file



class VitisExpressionData:
    """
    Loader for Vitis vinifera expression data in the described folder structure.

    Expects:
        data_root/
          vitis_vinifera/
            metadata/
            control/
              leaf/
                GSE*_control*.txt
              berry/
                ...
    """

    def __init__(
        self,
        data_root: str,
        subset: str = "control",
        intersect_genes: bool = True,
    ) -> None:
        """
        Parameters
        ----------
        data_root : str
            Path to the folder that contains 'vitis_vinifera'.
        subset : str
            Subfolder under vitis_vinifera (e.g. 'control').
        intersect_genes : bool
            If True, restrict to genes present in ALL loaded datasets.
        """
        self.data_root = os.path.abspath(data_root)
        self.subset = subset
        self.intersect_genes = intersect_genes

        self._counts: Optional[pd.DataFrame] = None
        self._metadata: Optional[pd.DataFrame] = None

        self._load_all()


    @property
    def counts(self) -> pd.DataFrame:
        """Raw counts matrix (genes x samples)."""
        return self._counts


    @property
    def metadata(self) -> pd.DataFrame:
        """Sample metadata table (one row per sample)."""
        return self._metadata


    @property
    def genes(self) -> List[str]:
        return list(self.counts.index)


    @property
    def samples(self) -> List[str]:
        return list(self.counts.columns)


    def _load_all(self) -> None:
        vitis_dir = os.path.join(self.data_root, "vitis_vinifera")
        subset_dir = os.path.join(vitis_dir, self.subset)

        if not os.path.isdir(subset_dir):
            raise FileNotFoundError(f"Subset dir not found: {subset_dir}")

        per_file_dfs: List[pd.DataFrame] = []
        sample_infos: List[VitisSampleInfo] = []

        # each tissue subfoler (like 'leaf', 'berry', etc.)
        for tissue in sorted(os.listdir(subset_dir)):
            tissue_dir = os.path.join(subset_dir, tissue)

            # check if it's a directory
            if not os.path.isdir(tissue_dir):
                continue

            # each txt file in the tissue folder
            for fname in sorted(os.listdir(tissue_dir)):
                if not fname.lower().endswith(".txt"):
                    continue

                fpath = os.path.join(tissue_dir, fname)

                df, infos = self._load_single_file(
                    fpath=fpath,
                    tissue=tissue,
                    condition=self.subset,
                )

                per_file_dfs.append(df)
                sample_infos.extend(infos)

        if not per_file_dfs:
            raise RuntimeError(f"No expression txt files in {subset_dir}")

        # restrict to intersecting genes
        # LW: I think this is the right default, but we can make it optional
        if self.intersect_genes:
            common_genes = set(per_file_dfs[0].index)

            for df in per_file_dfs[1:]:
                common_genes &= set(df.index)

            if not common_genes:
                raise RuntimeError("No intersecting genes across datasets.")
            
            common_genes = sorted(common_genes)
            per_file_dfs = [df.loc[common_genes] for df in per_file_dfs]

        # concat all columns together
        counts = pd.concat(per_file_dfs, axis=1)
        counts.index.name = "Gene"

        # build dataframe of sample metadata
        meta = pd.DataFrame([s.__dict__ for s in sample_infos])
        meta = meta.set_index("sample_id").loc[counts.columns]  # align ordering

        self._counts = counts
        self._metadata = meta


    def _load_single_file(
        self,
        fpath: str,
        tissue: str,
        condition: str,
    ) -> Tuple[pd.DataFrame, List[VitisSampleInfo]]:
        """
        Load a single .txt count file and return:

        - counts: DataFrame (genes x samples for this file)
        - sample_infos: metadata for each column
        """
        df = pd.read_csv(fpath, sep="\t", header=0)
        if df.shape[1] < 2:
            raise ValueError(f"File has no sample columns: {fpath}")

        # Ensure first column is 'Gene'
        gene_col = df.columns[0]
        if gene_col.lower() != "gene":
            df = df.rename(columns={gene_col: "Gene"})
        else:
            df = df.rename(columns={gene_col: "Gene"})
        df = df.set_index("Gene")

        if not df.index.is_unique:
            df = df.groupby(df.index).sum()

        # build sample IDs and metadata
        sample_infos: List[VitisSampleInfo] = []
        new_cols: Dict[str, str] = {}
        gse = self._parse_gse_from_filename(os.path.basename(fpath))
        cultivar = self._parse_cultivar_from_filename(os.path.basename(fpath))

        for col in df.columns:
            original_name = col
            sample_id = f"{tissue}|{gse}|{original_name}"

            info = VitisSampleInfo(
                sample_id=sample_id,
                original_name=original_name,
                tissue=tissue,
                condition=condition,
                gse_accession=gse,
                cultivar=cultivar,
                file_path=fpath,
            )
            sample_infos.append(info)
            new_cols[col] = sample_id

        df = df.rename(columns=new_cols)

        return df, sample_infos


    @staticmethod
    def _parse_gse_from_filename(fname: str) -> str:
        """
        'GSE97900_control.txt' -> 'GSE97900'
        """
        m = GSE_REGEX.search(fname)
        if m:
            return m.group(1).upper()
        return "UNKNOWN"


    @staticmethod
    def _parse_cultivar_from_filename(fname: str) -> Optional[str]:
        """
        Try to get cultivar from filename
        'GSE97900_control.txt' -> None
        'GSE12345_leaf_chardonnay.txt' -> 'chardonnay'
        """
        lower = fname.lower()
        if "mueller" in lower:
            return "mueller"
        
        if "regent" in lower:
            return "regent"
        
        return None