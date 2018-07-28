import sys
sys.path.insert(0,'..')
import simushape as ss 
from scipy.stats import spearmanr as spear
import pandas
import numpy as np
from rna_tools.rnaplfold import rnaplfold
import rna_tools.rna_io as rio
from rna_tools.sukosd import sukosd

from rna_tools.rna_accuracy import get_structure_accuracy


def run_accuracy(data, keys):
    # get shaker predictions: 


    def get_str(seq):
	for _, data_seq, db in data.values():
	    if seq == data_seq:
		return [db]
	print "SHIT"

    predictions = { k:(data[k][1],data[k][2],v) for v,k  in zip(  ss.crosspredict_nfold(data,keys) , keys) }
    predictions_str = { k:(data[k][1],data[k][2],v) for v,k in zip(  ss.crosspredict_nfold(data,keys,seq_to_db_function=get_str) , keys)}

    acc_shaker = {k: get_structure_accuracy(*predictions[k])  for  k in keys} 
    acc_shaker_plain = {k: get_structure_accuracy(*predictions_str[k]) for  k in keys} 
    predictions_suko = { k: (data[k][1],data[k][2],sukosd(data[k][2])) for k in keys }
    acc_suko = {k: get_structure_accuracy(*predictions_suko[k])  for  k in keys} 
    acc_real = {k: get_structure_accuracy(data[k][1],data[k][2],data[k][0])  for  k in keys}
    acc_noshape = {k: get_structure_accuracy(data[k][1],data[k][2],None)  for  k in keys}

    index= ['suko',"shaker","shaker_plain",'real',"no_shape"]
    data = map(lambda x: [x[k] for k in keys] ,[acc_suko, acc_shaker,acc_shaker_plain, acc_real, acc_noshape])
    map(lambda x:x.append(np.mean(x)),data)
    df = pandas.DataFrame(data, columns=keys+['mean'], index=index).T
    print df
    print df.to_latex()





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
run_accuracy(data,keys)

