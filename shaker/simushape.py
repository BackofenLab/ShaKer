import numpy as np
from scipy.sparse import vstack
import eden
import rna_tools.util as util
from eden import graph as eg
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import normalize
import xgboost
from scipy.stats import uniform as uni
from rna_tools.rnasubopt import  rnasubopt

from rna_tools.structureprobability import probabilities_of_structures


from sklearn.model_selection import KFold
def make_forestregressor():
    return RandomForestRegressor(**{'oob_score': False,
                                    'min_impurity_decrease': 0.01,
                                    'bootstrap': True,
                                    'min_samples_leaf': 1,
                                    'n_estimators': 16,
                                    'min_samples_split': 6,
                                    'min_weight_fraction_leaf': 0.02,
                                    'max_features': None})

def make_xgbreg():
    return xgboost.XGBRegressor(
            **{'reg_alpha': 0.81547748872761927, 'learning_rate': 0.03, 'max_delta_step': 1, 'min_child_weight': 3, 'n_estimators': 65, 'reg_lambda': 0.93307324674007364, 'max_depth': 14, 'gamma': 0, 'booster': 'gbtree'}
                )

def crosspredict_nfold(data, keys, seq_to_db_function=rnasubopt, n_splits=3, model=make_xgbreg()):
    '''
    data = {seqname:[shapearray, sequence, structure]}

    train on n-1 keys, predict on the last,
    yield result for each
    '''
    print "crosspredict:",
    res = []
    kf = KFold(n_splits=n_splits,shuffle=False)
    for train_n, test_n in kf.split(keys):
        train_keys = [keys[i] for i in train_n]
        mod = make_model(data,train_keys, model)
        for i in test_n:
            res.append( predict(mod, data[keys[i]][1], seq_to_db_function=seq_to_db_function))
            print ".",
    print "\n"
    return res



def crosspredict(data, keys, seq_to_db_function=rnasubopt):
    '''
    data = {seqname:[shapearray, sequence, structure]}

    train on n-1 keys, predict on the last,
    yield result for each
    '''
    print "crosspredict:",
    for key in keys:
        trainkeys = remove(keys, key)
        mod = make_model(data,trainkeys, model= make_forestregressor())
        print ".",
        yield predict(mod, data[key][1], seq_to_db_function=seq_to_db_function)




def remove(li, it):
    '''returns copy of li(st) without "it"'''
    li2 = list(li)
    li2.remove(it)
    return li2






def make_model( data,
                sequence_names=[],
                model= make_xgbreg()):
    x,y = getXY(data,sequence_names)
    model.fit(x,y)
    return model


def mask(x,y):
    mask = np.array([ i for i,e in enumerate(y) if e!=None])
    y= np.array(y)
    y=y[mask]
    x=x[mask]
    return x,y


def getXY(data,keys):
    '''takes entries in data that are in the list keys, returns X,Y for regression task'''
    # data is  name -> (react,sequence,dotbacket)
    # we first make some graphs
    react,sequence,stru = zip(*[ data[k] for k in keys ])
    graphs  = map( getgraph, sequence,stru)

    # then we edenize
    x = vstack( eg.vertex_vectorize(graphs,r=3,d=3))
    y= [y for reactlist in react for y in reactlist]
    y= np.array(y)
    # then done
    #print x,y
    return mask(x,y)



def getgraph(sequence,structure):
    """returns networkx graph"""
    return util.sequence_dotbracket_to_graph(sequence,structure)



def weighted_average(weights, react_arrays):
    '''
    generates l1 norm of weights (sum of absolutes =1)
    then multiplies with the react_arraus,

    t = [np.array(range(3)) for i in range(3)]
    g = [.5] * 3
    d = [t*g for t , g in zip(t,g)]
    print d
    print sum(d)
    '''
    weights = np.reshape(weights, (1,-1))
    weights = normalize(weights, norm='l1').tolist()[0]
    return sum([ array*weight for array, weight in zip(react_arrays,weights) ])


def predict(model, sequence,seq_to_db_function= rnasubopt):
    db_list = seq_to_db_function(sequence)

    if len(db_list)==1:
        graph = util.sequence_dotbracket_to_graph(sequence, db_list[0])
        return model.predict(eg.vertex_vectorize([graph])[0])

    # get probability for each structure
    struct_proba = probabilities_of_structures(sequence, db_list)
    structures, weights =  zip(*struct_proba)

    # edenize and predict reacticuty
    graphs = map(lambda x: getgraph(sequence,x), structures)
    vecs = list(eg.vertex_vectorize(graphs,r=3,d=3))
    predictions_all_structures = [ model.predict(blob) for blob in vecs ]

    # mix reactivity with probabilities
    return weighted_average(weights, predictions_all_structures)


def main():
    import sys
    import shaker.rna_tools.rna_io as rio
    import pickle
    
    helpstr =  '''
    I am displaying this help-message becuase your parameters are bad

    USAGE:
    1. train a model:
    shaker makemodel react-file dbn-file output-model-file

    2. use a model:
    shaker predict model-file sequence
    '''


    if len(sys.argv) < 3:
        print helpstr
    elif  sys.argv[1]  == "makemodel" and len(sys.argv)==5:
	data = rio.get_all_data(sys.argv[2],sys.argv[3])
	model  = make_model(data,data.keys())
	with open(sys.argv[4], 'wb') as file:
	    pickle.dump(model, file)

    elif  sys.argv[1]  == "predict" and len(sys.argv)==4:
	with open(sys.argv[2], 'rb') as file:
	    model = pickle.load(file)
	res = predict(model,sys.argv[3])
        print rio.format_shape("result", res)
    else:
        print helpstr

def test():
        import rna_tools.rna_io as rio
        #  Train a model 
        data = rio.get_all_data("../data/RNA16.react","../data/RNA16.dbn") 
        model  = make_model(data,data.keys())

        # Predict 
        print (predict(model,"AAAAAAGGGGCCCCCCCGGGGGUUUUUU"))


























