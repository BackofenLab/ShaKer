

import rna_tools.rna_io as rio
import simushape as sim
DEBUG = False


# LAOD DATA
data = rio.get_all_data("data/RNA16.react","data/RNA16.dbn")  # {key: rea, seq, stru}
if DEBUG: 
    for k,v in data.items():print k,"\t",  len(v[1]) 




# TRAIN MODEL 
if DEBUG:
    model  = sim.make_model(data,data.keys()[4:6]) 
else: 
    model  = sim.make_model(data,data.keys())


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









