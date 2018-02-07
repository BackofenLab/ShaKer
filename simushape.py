


#####
# READ AND WRITE DATA
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

def read_fasta(path):
    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res=[]
        for e in text:
            a_thing = [thing for thing in e.split('\n') if len(thing) > 1]
            if len(a_thing)!=2:
                print "ERRER", a_thing ,e
                return
            res.append(a_thing)
    return res


def float_or_none(item):
    if item == 'NA':
        return None
    return float(item)


def read_react(path):
    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res={}
        for e in text:
            lines = [thing for thing in e.split('\n') if len(thing) > 1]
            header = lines[0]
            data = [ e.strip().split() for e in lines[1:]  ]
            data = [ float_or_none(d[1].strip()) for d in data if len(d)==2 ]
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


def dump_shape(result, fname):
    with open(fname,'w') as f:
        for k,v in result.items():
            vout = lambda x: "\n".join( [str(i + 1) + "\t" + str(e) for i, e in enumerate(v)])
            f.write(">%s\n%s\n\n" % (k,vout(v)))





#######
# generate data _  struct given
########

import eden_rna
import eden.graph as eg
from scipy.sparse import vstack
import numpy as np

from sklearn.ensemble import RandomForestRegressor
def mask(x,y):
    mask = np.array([ i for i,e in enumerate(y) if e!=None])
    y=y[mask]
    x=x[mask]
    return x,y


def make_graphs(data, good_keys):
    return  [ eden_rna.sequence_dotbracket_to_graph(*data[key][1:]) for key in good_keys]

def getXY(data,good_keys,r=3,d=3,DEBUG=False,n_bits=16):
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

# {'reg_alpha': 0.81547748872761927, 'learning_rate': 0.03, 'max_delta_step': 1, 'min_child_weight': 3, 'n_estimators': 65, 'reg_lambda': 0.93307324674007364, 'max_depth': 14, 'gamma': 0, 'booster': 'gbtree'}
def make_model(data,sequence_names=[],
               DEBUG=False,
               r=3,
               d=3,
               model=RandomForestRegressor(**{'oob_score': False,
                                              'min_impurity_split': 0.01,
                                              'bootstrap': True,
                                              'min_samples_leaf': 1,
                                              'n_estimators': 16,
                                              'min_samples_split': 6,
                                              'min_weight_fraction_leaf': 0.02,
                                              'max_features': None})):

    x,y = getXY(data,sequence_names,r,d,DEBUG=DEBUG)
    model.fit(x,y)
    return model







##############
# PREDICT          _ struct given
#############

def predict2(model,seq,stru):
    graph = eden_rna.sequence_dotbracket_to_graph(seq,stru)
    return model.predict(eg.vertex_vectorize([graph])[0])

def predict(model,graph):
    #graph = eden_rna.sequence_dotbracket_to_graph(seq,stru)
    return model.predict(eg.vertex_vectorize([graph])[0])




'''import argparse

def get_argparser():
    parser = argparse.ArgumentParser(description='shape annotation')
    parser.add_argument('--react',  type=str, nargs=1, default='data.react', help='react (>seqname\\n1\\t shapevalue\\n2\\t shapevalue etc)')
    parser.add_argument('--dbn',  type=str, nargs=1, default='data.dbn', help=' dbn files (>seqname\\nsequence\\ndotbacked etc)')
    parser.add_argument('--target',  type=str, nargs=1, default='target.dbn', help='we predict on these...  dbn files (>seqname\\nsequence\\ndotbacked etc)')
    parser.add_argument('--output',  type=str, nargs=1, default='out.react', help='we predict on these...  dbn files (>seqname\\nsequence\\ndotbacked etc)')
    parser.add_argument('--kernelargs',  type=dict, nargs=1, default={}, help='dictionary for the eden constructor.. UNUSED ')
    return parser
    
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
'''
