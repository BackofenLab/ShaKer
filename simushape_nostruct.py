###
# this is hard cancer, it trains a model on some rna shapes structures,  also precition works on a set of structures
###


import numpy as np
from scipy.sparse import vstack
import eden
import eden_rna
from sklearn.ensemble import RandomForestRegressor
import simushape
from rna_tools import rnashapes



def getgraphs(sequences, react,maxstructs=3):
   for seq,r in zip(sequences,react):
       for stru in rnashapes(seq)[0][:maxstructs]: # MAXIMUM structures looked at
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




def predict_shreps_proba(model, sequence, cutoff=0.01):
    structs, probabilities = r
    graphs = getgraphs([sequence],['None'], maxstructs=maxstruct)
    vecs = eden.graph.vertex_vectorize(graphs,r=3,d=3)
    res= [ model.predict(blob) for blob in vecs ]
    res = np.vstack(res)
    #return res.mean(axis=0)
    return np.median(res,axis=0)

def predict(model, sequence, maxstruct=3):
    graphs = getgraphs([sequence],['None'], maxstructs=maxstruct)
    vecs = eden.graph.vertex_vectorize(graphs,r=3,d=3)
    res= [ model.predict(blob) for blob in vecs ]
    res = np.vstack(res)
    #return res.mean(axis=0)
    return np.median(res,axis=0)

def predict_and_scale(model, sequence, maxstruct=3):
    graphs = getgraphs([sequence],['None'], maxstructs=maxstruct)
    vecs = eden.graph.vertex_vectorize(graphs,r=3,d=3)
    res= [ model.predict(blob) for blob in vecs ]
    res = np.vstack(res)
    med= np.median(res,axis=0)

    scal =  res.max()/ med.max()
    return med, med*scal


def remove(li, it):
    li2 = list(li)
    li2.remove(it)
    return li2


def crossval(data,keys):
    for key in keys():
        trainset = remove(keys, key)
        mod = make_model(data,trainset)
        res= predict(mod,data[key][1])