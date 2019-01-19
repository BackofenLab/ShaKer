# ShaKer

### Usage


```python

import ShaKer.rna_tools.rna_io as rio
import ShaKer.simushape as sim

#  Train a model 
data = rio.get_all_data("../data/RNA16.react","../data/RNA16.dbn") 
model  = sim.make_model(data,data.keys())

# Predict 
sim.predict(model,"AAAAAAGGGGCCCCCCCGGGGGUUUUUU")

```

### Installation

```fish
# get vienna RNA:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install viennarna

# this will take care of most of the python dependencies:
pip install graphlearn xgboost
pip install git+https://github.com/smautner/EDeN.git --user 

# also clone shaker and put it in your python path 
git clone https://github.com/smautner/ShaKer 
set -x PYTHONPATH .  # (this is the fish shell command to set a variable)
```

