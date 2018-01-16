
import argparse
def get_argparser():
    parser = argparse.ArgumentParser(description='shape annotation')
    parser.add_argument('--react',  type=str, nargs=1, default='data.react', help='react (>seqname\\n1\\t shapevalue\\n2\\t shapevalue etc)')
    parser.add_argument('--dbn',  type=str, nargs=1, default='data.dbn', help=' dbn files (>seqname\\nsequence\\ndotbacked etc)')
    parser.add_argument('--target',  type=str, nargs=1, default='target.dbn', help='we predict on these...  dbn files (>seqname\\nsequence\\ndotbacked etc)')
    parser.add_argument('--output',  type=str, nargs=1, default='out.react', help='we predict on these...  dbn files (>seqname\\nsequence\\ndotbacked etc)')
    parser.add_argument('--kernelargs',  type=dict, nargs=1, default={}, help='dictionary for the eden constructor.. UNUSED ')
    return parser



#####
# READ DATA
#####
def read_dbn(path):
    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res=[]
        for e in text:
            a_thing = [thing for thing in e.split('\n') if len(thing) > 1]
            if len(a_thing)!=3:
                print "ERRER", a_thing ,e
                return
            if len(a_thing[1])!=len(a_thing[2]):
                print "ERRER", a_thing ,e
                return
            res.append(a_thing)
    return res

def read_react(path):
    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res={}
        for e in text:
            lines = [thing for thing in e.split('\n') if len(thing) > 1]
            header = lines[0]
            data = [ e.strip().split() for e in lines[1:]  ]
            data = [ float(d[1].strip()) for d in data if len(d)==2 ]
            res[header.strip()] = data
    return res


def combine_dbn_react(dbn,react):
    res = {}
    for name,seq,brack in dbn:
        re = react[name]
        if len(re) == len(seq) == len(brack):
            res[name]= (re,seq,brack)
        else:
            print "data for '%s' is corrupted, ignoring..." % name
    return res

def get_all_data(react, dbn):
    dbn = read_dbn(dbn)
    react = read_react(react)
    return combine_dbn_react(dbn,react)






#######
# generate data and learn model
########

import eden_rna
import eden.graph as eg
from scipy.sparse import vstack
from sklearn.linear_model import SGDRegressor
import numpy as np

from sklearn.ensemble import RandomForestRegressor
def mask(x,y):
    mask = np.array([ i for i,e in enumerate(y) if e!=-999.0])
    y=y[mask]
    x=x[mask]
    return x,y



def make_graphs(data, good_keys):
    return  [ eden_rna.sequence_dotbracket_to_graph(*data[key][1:]) for key in good_keys]

def getXY(data,good_keys,r,d,DEBUG=False):
    graphs =  make_graphs(data,good_keys)
    x = vstack( eg.vertex_vectorize(graphs,r=r,d=d) )

    if DEBUG:
        print "x matrix: ", x.shape
        print "graphs", len(graphs),
        nodecount = [len(g) for g in graphs]
        print nodecount, "sum", sum(nodecount)

    y = [f for key in good_keys for f in data[key][0]]
    if DEBUG:
        print "found this many y values:", len(y)
    # filter x and y
    y=np.array(y)
    x,y= mask(x,y)
    if DEBUG:
        print "x,y after filtering",x.shape, y.shape
    return x,y

def make_model(data,sequence_names=[], DEBUG=False,r=0,d=2,model=RandomForestRegressor()):
    x,y = getXY(data,sequence_names,r,d,DEBUG=DEBUG)
    model.fit(x,y)
    return model







##############
# PREDICT
#############
from graphlearn01.utils import draw as draw

def predict2(model,seq,stru):
    graph = eden_rna.sequence_dotbracket_to_graph(seq,stru)
    return model.predict(eg.vertex_vectorize([graph])[0])

def predict(model,graph):
    #graph = eden_rna.sequence_dotbracket_to_graph(seq,stru)
    return model.predict(eg.vertex_vectorize([graph])[0])




def dump_shape(result, fname):
    with open(fname,'w') as f:
        for k,v in result.items():
            vout = lambda x: "\n".join( [str(i + 1) + "\t" + str(e) for i, e in enumerate(v)])
            f.write(">%s\n%s\n\n" % (k,vout(v)))






##########
# OPTIMIZE
##########


def remove(li, it):
    li2 = list(li)
    li2.remove(it)
    return li2

from scipy.stats import spearmanr as corr
def calc(names,data, item,r,d,model=RandomForestRegressor()):
    model = make_model(data,names,False,r,d, model)
    graph = eden_rna.sequence_dotbracket_to_graph(data[item][1],data[item][2])
    res = np.array(predict(model,graph))
    other = np.array(data[item][0])

    res,other = mask(res,other)
    value =  corr(res,other)[0]
    #print '\t',len(data[item][1]),"\t", value
    return value

import multiprocessing as mp
def funmap(f,args):
    pool=mp.Pool(10)
    res=pool.map(f,args)
    pool.close()
    return res
#------
# EDEN


def opti_eden(data,r,d):
    names= data.keys()
    res= [ calc( remove(names,item), data , item,r,d) for item in names ]
    mederror= np.array(res).mean()
    print "r=%d d=%d res=%.2f\n\n" % (r,d, mederror)

def optimize_eden_multi(data, jobs=4):
    funmap(opti_eden,[(data,r,d) for r in range(5) for d in range(5)])


def optimize_eden_serial(data):
    # data is now  a dict name -> (react,sequence,dotbacket)
    for r in range(0,5):
        for d in range(0,5):
            opti_eden(data,r,d)



from sklearn.model_selection import RandomizedSearchCV as rsearch
from scipy.stats import randint as rint
from sklearn.ensemble import RandomForestRegressor

def quickladdata():
    ok = ['ADDRSW', 'ZHCV', 'Z-CIDGMP-1', '23sRNA', 'p564', 'srRNA', 'MDLOOP', 'R009', 'TRP5']
    data = get_all_data('data/RNA16.react','data/RNA16.dbn')
    for k in data.keys():
        if k not in ok:
            data.pop(k)
    return data

def opti_forest(data,r=3,d=3, n_jobs=1,n_iter=10):
    model = RandomForestRegressor()
    param_dist = {'n_estimators': rint(15, 30),
                  # 'criterion': ['mse','mae'],  # not in 0.18 but in 0.19
                  'min_samples_split': rint(2, 10),
                  'min_samples_leaf': rint(1, 5),
                  'min_weight_fraction_leaf': [0.02],# opti said this is good
                  'max_features': [None], # None is best
                  'min_impurity_split': [0.03, 0.02, 0.01, 0.04],  # min_impurity_decrease
                  "bootstrap": [True],  # false conflicts with oob score thing
                  "oob_score": [False, True]}

    X,y = getXY(data,data.keys(),r,d)
    blu = rsearch(model, param_distributions=param_dist, n_iter=n_iter,n_jobs=n_jobs)
    blu.fit(X, y)
    print blu.best_params_
    print blu.best_score_



if __name__=='__main__':

    parser= get_argparser()
    args = parser.parse_args()
    DEBUG = True
    if DEBUG:
        print 'argparse args: ', args

    dbn = read_dbn(args.dbn)
    react = read_react(args.react)
    data= combine_dbn_react(dbn,react)
    # data is now  a dict name -> (react,sequence,dotbacket)
    model=make_model(data,data.keys(), DEBUG)
    targets ={a:[b,c] for (a,b,c) in read_dbn(args.target)}
    result = {}
    dump_shape(result,args.output)
    #for target in targets:
    #    result[target] = predict2(model,*targets[target])
    #optimize(data)
