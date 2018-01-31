from scipy.stats import randint as rint
from scipy.stats import uniform as uni
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV as rsearch
import simushape as ss


def quickladdata():
    data = ss.get_all_data('data/RNA16.react','data/RNA16.dbn')
    data2 = ss.get_all_data('data/RNA20.react','data/RNA20.dbn')
    data.update(data2)
    for e in ['ZHCV', 'Lysine', 'GLYCFN']:
        data.pop(e)
    return data

data= quickladdata()

def uniform(lower, upper):
    return uni(lower, upper-lower)

def opti_forest(data,r=3,d=3, n_jobs=1,n_iter=10):
    model = GradientBoostingRegressor()
    param_dist = {
                     'loss': ['ls', 'lad', 'huber', 'quantile'],
                     'learning_rate': uniform(0.01,0.2),
                     'max_depth': rint(1, 6),
                     'min_samples_split': rint(2, 5),
                     'min_samples_leaf': rint(1, 3),
                     'min_weight_fraction_leaf': uniform(0,.05),
                     'max_features': ["auto", 'sqrt' , 'log2'],
                     'alpha': uniform(0.8,0.97)
    }

    X,y = ss.getXY(data,data.keys(),r,d)
    X=X.toarray()
    blu = rsearch(model, param_distributions=param_dist, n_iter=n_iter,n_jobs=n_jobs)
    blu.fit(X, y)

    print blu.best_params_
    print blu.best_score_


opti_forest(data, n_jobs=4, n_iter=8)
