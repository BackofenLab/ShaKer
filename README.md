# ShaKer

### Usage


```python

import rna_tools.rna_io as rio
import simushape as sim

#  Train a model 
data = rio.get_all_data("../data/RNA16.react","../data/RNA16.dbn") 
model  = sim.make_model(data,data.keys())

# Predict 
sim.predict(model,"AAAAAAGGGGCCCCCCCGGGGGUUUUUU")

```
