import simushape
import numpy as np
import eden_rna
from scipy.sparse import vstack
import eden
from sklearn.ensemble import RandomForestRegressor
#####
# laod  data is same as simushape.py
####

#######
# wrap shapes
######

import subprocess
import re

def shexec(cmd):
    '''
    :param cmd:
    :return: (exit-code, stderr, stdout)

    the subprocess module is chogeum.. here is a workaround
    '''
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, stderr = process.communicate()
    retcode = process.poll()
    return (retcode, stderr, output)


def shape(sequence):
    retcode,err,out = shexec("RNAshapes -u -s -t 5 -c 10 %s | sort -n " % sequence)
    if retcode != 0:
        print "RNAshapes failed"
        return

    energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)
    shape =[ a.strip() for a in re.findall(r' [.()]+ ',out) ]
    return shape, energy


def getgraphs(sequences, react,maxstructs=3):
   for seq,r in zip(sequences,react):
       for stru in shape(seq)[0][:maxstructs]: # MAXIMUM structures looked at
           graph=  eden_rna.sequence_dotbracket_to_graph(seq,stru)
           graph.graph['rea']=r
           yield graph
#######
#  train
#######

def getXY(data,keys, maxstru=3):
    # data is  name -> (react,sequence,dotbacket)
    # we first make some graphs
    graphs = list(getgraphs([data[k][1] for k in keys], [data[k][0] for k in keys ], maxstructs=maxstru))

    # then we edenize
    x = vstack( eden.graph.vertex_vectorize(graphs,r=3,d=3) )

    y= [y for g in graphs for y in g.graph['rea']]
    y= np.array(y)
    # then done
    #print x,y
    return simushape.mask(x,y)

def make_model(data,sequence_names=[],
               model=RandomForestRegressor(**{'oob_score': False,
                                              'min_impurity_split': 0.01,
                                              'bootstrap': True,
                                              'min_samples_leaf': 1,
                                              'n_estimators': 16,
                                              'min_samples_split': 6,
                                              'min_weight_fraction_leaf': 0.02,
                                              'max_features': None}),
               maxstruct=3):

    x,y = getXY(data,sequence_names, maxstru=maxstruct)
    model.fit(x,y)
    return model

def predict(model, sequence, maxstruct=3):
    graphs = getgraphs([sequence],['None'], maxstructs=maxstruct)
    vecs = eden.graph.vertex_vectorize(graphs,r=3,d=3)
    res= [ model.predict(blob) for blob in vecs ]
    res = np.vstack(res)
    #return res.mean(axis=0)
    return np.median(res,axis=0)





def remove(li, it):
    li2 = list(li)
    li2.remove(it)
    return li2


def crossval(data,keys):

    for key in keys():
        trainset = remove(keys, key)
        mod = make_model(data,trainset)
        res= predict(mod,data[key][1])