import sys
sys.path.insert(0,'..')
import simushape as ss 
from scipy.stats import spearmanr as spear
import pandas
import numpy as np
from rna_tools.rnaplfold import rnaplfold
import rna_tools.rna_io as rio
from rna_tools.sukosd import sukosd


def run_access_old(data, keys):

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
    


def run_access(data, keys):
    # get shaker predictions: 


    def get_str(seq):
	for _, data_seq, db in data.values():
	    if seq == data_seq:
		return [db]
	print "SHIT"

    predictions = { k:(data[k][1],v) for v,k  in zip(  ss.crosspredict_nfold(data,keys) , keys) }
    predictions_str = { k:(data[k][1],v) for v,k in zip(  ss.crosspredict_nfold(data,keys,seq_to_db_function=get_str) , keys)}
    

    acc_shaker = {k: rnaplfold(*predictions[k]) for  k in keys} # a predcition is seq,rea
    acc_shaker_str = {k: rnaplfold(*predictions_str[k]) for  k in keys} # a predcition is seq,rea

    acc_suko =   {k: rnaplfold(data[k][1], sukosd(data[k][2]) ) for  k in keys}
    acc_shape =  {k: rnaplfold(data[k][1],data[k][0]) for  k in keys}  
    acc_nodata = {k: rnaplfold(data[k][1]) for  k in keys}
    corr = lambda ac1, ac2:[spear(ac1[k], ac2[k])[0] for k in keys]
    corr_suko =  corr(acc_suko, acc_shape)
    corr_shaker =  corr(acc_shaker, acc_shape)
    corr_shaker_str =  corr(acc_shaker_str, acc_shape)
    corr_nodata =  corr(acc_nodata, acc_shape)


    index= ['suko',"shaker","shaker_str","nodata"]
    data = [corr_suko,corr_shaker, corr_shaker_str,corr_nodata]
    map(lambda x:x.append(np.mean(x)),data)
    asd = pandas.DataFrame(data, columns=keys+['mean'], index=index).T
    print asd 
    print asd.to_latex()

# TODO: add 4th

types = ["cellfree",
"incell",
"kasugamycin"]

def getdata(typ):
    return rio.get_all_data("../data/weeks194_orig/%s.react" % typ,"../data/weeks194_orig/%s.dbn" % typ)  # {key: rea, seq, stru}

#def getdata():
#    return rio.get_all_data("../data/RNA16.react"  ,"../data/RNA16.dbn" )  # {key: rea, seq, stru}

#data=getdata(types[1])
data= getdata("incell_nogenes")
keys=data.keys()[:9]
run_access(data,keys)

