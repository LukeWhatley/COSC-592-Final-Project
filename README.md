# COSC-592-Final-Project
Final project repository for COSC 592 - Bioinformatic Computing

## Setup & Installation
### Requirements

- pyenv
- pyenv-virtualenv
- Python 3.10.0

### Clone the Repository
```bash
git clone https://github.com/LukeWhatley/COSC-592-Final-Project.git
cd COSC-592-Final-Project
```

### Create environment
```bash
pyenv install 3.10.0
pyenv virtualenv 3.10.0 vitis
pyenv local vitis
pyenv activate vitis
python --version
```

### Install dependencies
```bash
python -m ensurepip --upgrade
python -m pip install --upgrade pip setuptools wheel

pip install --upgrade pip
pip install -e .
```