from scipy.stats import randint as rint
from scipy.stats import uniform as uni

#from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV as rsearch

from rna_tools import rna_io
import simushape as ss
import xgboost


def quickladdata():
    data = rna_io.get_all_data('../data/RNA16.react', '../data/RNA16.dbn')
    data2 = rna_io.get_all_data('../data/RNA20.react', '../data/RNA20.dbn')
    data.update(data2)
    for e in ['ZHCV', 'Lysine', 'GLYCFN']:
        data.pop(e)
    return data


def uniform(lower, upper):
    return uni(lower, upper-lower)

def optimize_xgb(data, n_jobs=1, n_iter=10):
    model = xgboost.XGBRegressor()
    param_dist = {
                     'max_depth': [14], # 13 and 14 are bestd
                     'learning_rate': [0.03], # rate id 0.022 and nesti of 94 is also ok
                        'n_estimators': [65],
                          'booster': ['gbtree',  'dart'],
                    'gamma':[0,0.005,0.001],
        'min_child_weight' : [3], # best
'max_delta_step' : [1], # best
'reg_alpha' : uniform(.8,1),
'reg_lambda' : uniform(0.6,1)
    }

    X,y = ss.getXY(data,data.keys())
    #X = xgboost.DMatrix(X)
    #y= xgboost.DMatrix(y)
    blu = rsearch(model, param_distributions=param_dist, n_iter=n_iter,n_jobs=n_jobs)
    blu.fit(X, y)
    print blu.best_params_
    print blu.best_score_


if __name__ == '__main__':
    data= quickladdata()
    optimize_xgb(data, n_jobs=24, n_iter=3500)
# BEST {'reg_alpha': 0.81547748872761927, 'learning_rate': 0.03, 'max_delta_step': 1, 'min_child_weight': 3, 'n_estimators': 65, 'reg_lambda': 0.93307324674007364, 'max_depth': 14, 'gamma': 0, 'booster': 'gbtree'}
