import ShaKer.simushape as sim
import ShaKer.rna_tools.rna_io as rio
import ShaKer.rna_tools.util as util

types = ["cellfree",
"incell",
"kasugamycin"]

def getdata(typ):
    return rio.get_all_data("../data/weeks194_orig/%s.react" % typ,"../data/weeks194_orig/%s.dbn" % typ)  # {key: rea, seq, stru}

typee = 'incell'
data= getdata(typee)

keys = data.keys()# 3 for testing
res = sim.crosspredict_nfold(data, keys,model=sim.make_xgbreg())

res = dict(zip(keys,res))
rio.dump_shape(res,"%s_predicted.react" % typee)

