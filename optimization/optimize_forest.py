from scipy.stats import randint as rint

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV as rsearch

import simushape as ss
from simushape import get_all_data, getXY



def quickloaddata(dataset='36'):
    
    if dataset == '36':
        data = ss.get_all_data('../data/RNA16.react','../data/RNA16.dbn')
        data2 = ss.get_all_data('../data/RNA20.react','../data/RNA20.dbn')
        data.update(data2)
    else:
        data = ss.get_all_data('../data/RNA%s.react' % dataset,'../data/RNA%s.dbn' % dataset)

    for e in ['ZHCV', 'Lysine', 'GLYCFN']:
        data.pop(e,None)
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
                  "oob_score": [False]}

    X,y = getXY(data,data.keys(),r,d)
    blu = rsearch(model, param_distributions=param_dist, n_iter=n_iter,n_jobs=n_jobs)
    blu.fit(X, y)
    print blu.best_params_
    print blu.best_score_

if __name__ == "__main__":
    print "*"*80
    print "16"
    data= quickloaddata('16')
    opti_forest(data, n_jobs=24, n_iter=1000)
    print "20"
    data= quickloaddata('20')
    opti_forest(data, n_jobs=24, n_iter=1000)
    print "36"
    data= quickloaddata('36')
    opti_forest(data, n_jobs=24, n_iter=1000)
