

import rna_tools.rna_io as rio
import simushape as sim
DEBUG = False


################
## train on RNA16
################
'''
# LAOD DATA
data = rio.get_all_data("data/RNA16.react","data/RNA16.dbn")  # {key: rea, seq, stru}
if DEBUG: 
    for k,v in data.items():print k,"\t",  len(v[1]) 
# TRAIN MODEL 
if DEBUG:
    model  = sim.make_model(data,data.keys()[4:6]) 
else: 
    model  = sim.make_model(data,data.keys())
'''
##################
# train on mustoe incell
##################
# LAOD DATA
data = rio.get_all_data("data/weeks194_orig/incell.react","data/weeks194_orig/incell.dbn") 
# {key: rea, seq, stru}

from rna_tools.rnafold import fold
def preproc(data,keys):
    for k in keys:
        data[k]= (data[k][0],data[k][1], fold(data[k][1],data[k][0]))
        print ".",
    return data

# TRAIN MODEL 
if DEBUG:
    data = preproc(data,data.keys()[4:6])
    model  = sim.make_model(data,data.keys()[4:6]) 
else: 
    data = preproc(data,data.keys())
    model  = sim.make_model(data,data.keys())

print "training done"


######################
# ok hacking the fasta stuff 
########################

import os 
import rna_tools.rnaplfold as plf

# for all fasta files
for X in os.listdir("fasta"):
    fasta = rio.read_fasta("fasta/%s" % X) # [(name,seq)]

    # predict shape for all seq in fasta
    shapedata = [sim.predict(model,seq) for _,seq in fasta]

    # and run plfold 
    for  i,(rea,(_,seq)) in enumerate(zip(shapedata,fasta)):
        plf.call_vienna_plfold(seq,rea,L=100,W=150,u=150, seq_name= X[:-5]+"%d."%i)
        print ".",









