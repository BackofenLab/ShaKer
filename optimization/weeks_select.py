import ShaKer.simushape as sim
import ShaKer.rna_tools.rna_io as rio

types = ["cellfree",
"incell",
"kasugamycin"]

def getdata(typ):
    return rio.get_all_data("../data/weeks194_orig/%s.react" % typ,"../data/weeks194_orig/%s.dbn" % typ)  # {key: rea, seq, stru}

typee = 'incell'
data= getdata(typee)

keys = data.keys()[:3] # 9 for testing
res = sim.crosspredict_nfold(data, keys,model=sim.make_xgbreg())


from scipy.stats import spearmanr as spear
import numpy as np

def compare(pred, dat):
    relevant = np.array(dat)!=None
    return spear(np.array(pred)[relevant], np.array(dat)[relevant]) 

for a,b in zip(keys, res):
    print a, compare(b,data[a][0])


print data[a][0]
