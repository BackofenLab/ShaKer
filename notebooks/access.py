import sys
sys.path.insert(0,'..')
import simushape as ss 
from scipy.stats import spearmanr as spear
import pandas
import numpy as np
from rna_tools.rnaplfold import rnaplfold

def run_access(data, keys):
    '''this is shit because we can do better: take the real structures...
    now that i think about it i should also use the real ones for training '''

    shaker_shape = ss.crosspredict_nfold(data,keys)
    seqs = [data[k][1] for k in keys] 
    acc_shaker = [rnaplfold(seq, rea) for seq,rea in zip(seqs,shaker_shape)]
    acc_nothing = [rnaplfold(seq) for seq in seqs]
    acc_real  = [rnaplfold(data[k][1],data[k][0]) for k in keys]
    corr = lambda other: [spear(accr, acco)[0] for accr,acco in zip(acc_real, other)]
    corr_shaker = corr(acc_shaker)
    corr_default = corr(acc_nothing)
    
    index= ['shaker', 'default']
    data = [corr_shaker,corr_default]
    map(lambda x:x.append(np.mean(x)),data)
    print pandas.DataFrame(data, columns=keys+['mean'], index=index).T
    

import rna_tools.rna_io as rio


# TODO: add 4th

types = ["cellfree",
"incell",
"kasugamycin"]

def getdata(typ):
    return rio.get_all_data("../data/weeks194_orig/%s.react" % typ,"../data/weeks194_orig/%s.dbn" % typ)  # {key: rea, seq, stru}

def getdata():
    return rio.get_all_data("../data/RNA16.react"  ,"../data/RNA16.dbn" )  # {key: rea, seq, stru}

#data=getdata(types[1])
data= getdata()
keys=data.keys()
run_access(data,keys)

