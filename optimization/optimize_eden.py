import multiprocessing as mp
import numpy as np
from scipy.stats import spearmanr as corr

import eden_rna
from sklearn.ensemble import RandomForestRegressor

from ShaKer.simushape import make_model, predict, mask


def remove(li, it):
    li2 = list(li)
    li2.remove(it)
    return li2


def calc(names,data, item,r,d,model=RandomForestRegressor()):
    model = make_model(data,names,False,r,d, model)
    graph = eden_rna.sequence_dotbracket_to_graph(data[item][1],data[item][2])
    res = np.array(predict(model,graph))
    other = np.array(data[item][0])

    res,other = mask(res,other)
    value =  corr(res,other)[0]
    #print '\t',len(data[item][1]),"\t", value
    return value


def funmap(f,args):
    pool=mp.Pool(10)
    res=pool.map(f,args)
    pool.close()
    return res


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