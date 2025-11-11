from typing import Optional

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset


class ExpressionDataset(Dataset):
    """
    PyTorch Dataset for expression data.

    Expects:
        counts_df: genes x samples
        metadata_df: samples x metadata_columns
    """

    def __init__(
        self,
        counts_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
        label_col: str = "tissue",
        dtype: torch.dtype = torch.float32,
    ) -> None:
        """
        Parameters
        ----------
        counts_df : pd.DataFrame
            Counts or normalized expression (genes x samples).
        metadata_df : pd.DataFrame
            Metadata with index = sample_id and a column for labels.
        label_col : str
            Column in metadata_df to use as the label (e.g. 'tissue').
        dtype : torch.dtype
            Tensor dtype for features.
        """
        # Align metadata to counts columns
        metadata_df = metadata_df.loc[counts_df.columns]

        # features are genes
        X = counts_df.T.values  # shape: samples x genes

        self.X = torch.tensor(X, dtype=dtype)
        self.sample_ids = list(counts_df.columns)
        self.metadata = metadata_df.copy()
        self.label_col = label_col

        # labels into integer indices
        labels = metadata_df[label_col].astype("category")
        self.label_categories = list(labels.cat.categories)
        self.y = torch.tensor(labels.cat.codes.values, dtype=torch.long)

    def __len__(self) -> int:
        return self.X.shape[0]

    def __getitem__(self, idx: int):
        return self.X[idx], self.y[idx]

    def label_to_name(self, y_idx: int) -> str:
        return self.label_categories[y_idx]

    def get_sample_metadata(self, idx: int) -> pd.Series:
        sample_id = self.sample_ids[idx]

        return self.metadata.loc[sample_id]
