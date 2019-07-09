# ShaKer

### Usage


```python

import shaker.rna_tools.rna_io as rio
import shaker.simushape as sim

#  Train a model 
data = rio.get_all_data("data/RNA16.react","data/RNA16.dbn") 
model  = sim.make_model(data,data.keys())

# Predict 
print (sim.predict(model,"AAAAAAGGGGCCCCCCCGGGGGUUUUUU"))

```

### Installation

```fish
wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
# get vienna RNA binaries via conda or their ppa:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install viennarna

# python dependencies:
pip install xgboost seaborn tabulate toolz
pip install git+https://github.com/smautner/EDeN.git --user 
pip install shaker-rna
```


### mustoe data set
```fish
www.chem.unc.edu/rna/data-files/mustoe_2018_DATA_SOFTWARE.zip 
```
