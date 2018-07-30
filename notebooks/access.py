import sys
sys.path.insert(0,'..')
import simushape as ss 
from scipy.stats import spearmanr as spear
from scipy.stats import pearsonr as pear
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
    


def run_access(data, keys,make_model=lambda: ss.make_xgbreg()):
    # get shaker predictions: 


    def get_str(seq):
	for _, data_seq, db in data.values():
	    if seq == data_seq:
		return [db]
	print "SHIT"

    predictions = { k:(data[k][1],v) for v,k  in zip(  ss.crosspredict_nfold(data,keys,model=make_model()) , keys) }
    predictions_str = { k:(data[k][1],v) for v,k in zip(  ss.crosspredict_nfold(data,keys,seq_to_db_function=get_str, model=make_model()) , keys)}
    

    acc_shaker = {k: rnaplfold(*predictions[k]) for  k in keys} # a predcition is seq,rea
    acc_shaker_str = {k: rnaplfold(*predictions_str[k]) for  k in keys} # a predcition is seq,rea

    acc_suko =   {k: rnaplfold(data[k][1], sukosd(data[k][2]) ) for  k in keys}
    acc_shape =  {k: rnaplfold(data[k][1],data[k][0]) for  k in keys}  
    acc_nodata = {k: rnaplfold(data[k][1]) for  k in keys}

    access_data = [acc_suko,acc_shaker, acc_shaker_str, acc_nodata]
    index= ['suko',"shaker","shaker_str","nodata"]

    corr_spear = lambda ac1:[spear(ac1[k], acc_shape[k])[0] for k in keys]
    corr_pear = lambda ac1:[pear(ac1[k], acc_shape[k])[0] for k in keys]

    data_spear = map(corr_spear, access_data)
    data_pear = map(corr_pear, access_data)


    ress = map(lambda x:np.mean(x),data_spear)
    resp=  map(lambda x:np.mean(x),data_pear)

    print ress, resp 
    return "{}\n{} spear\n{} pear".format(index,ress ,resp)

    #map(lambda x:x.append(np.mean(x)),data)
    #asd = pandas.DataFrame(data, columns=keys+['mean'], index=index).T
    #return asd



def getdata(typ):
    return rio.get_all_data("../data/weeks194_orig/%s.react" % typ,"../data/weeks194_orig/%s.dbn" % typ)  # {key: rea, seq, stru}

#def getdata():
#    return rio.get_all_data("../data/RNA16.react"  ,"../data/RNA16.dbn" )  # {key: rea, seq, stru}



types = ["cellfree_genes",
"incell_genes",
"kasugamycin_genes"]
models = [ss.make_forestregressor, ss.make_xgbreg, ss.RandomForestRegressor]


tasks = [(ty,mo) for ty in types for mo in models ]

#data=getdata(types[1])

def run(id):
    datum , mod = tasks[id]
    data= getdata(datum)
    keys=data.keys()[:9]
    ret = run_access(data,keys,make_model=mod)
    with open("%d.access.out" % id, "w") as f: f.write(ret)


if __name__ == "__main__":
    import sys
    run(int(sys.argv[1]))


