# ShaKer

### Usage


```python

import ShaKer.rna_tools.rna_io as rio
import ShaKer.simushape as sim

#  Train a model 
data = rio.get_all_data("data/RNA16.react","data/RNA16.dbn") 
model  = sim.make_model(data,data.keys())

# Predict 
print (sim.predict(model,"AAAAAAGGGGCCCCCCCGGGGGUUUUUU"))

```

### Installation

```fish
# get vienna RNA binaries via conda or their ppa:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install viennarna

# python dependencies:
pip install graphlearn xgboost
pip install git+https://github.com/smautner/EDeN.git --user 

# clone shaker and put it in your python path:
git clone https://github.com/smautner/ShaKer 
set -x PYTHONPATH .  # (this is the fish shell command to set a variable)
# or directly:
git clone https://github.com/smautner/ShaKer ~/.local/lib/python2.7/site-packages/ShaKer

```

