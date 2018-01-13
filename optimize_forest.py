import simushape as ss


data=ss.quickladdata()
ss.opti_forest(data,n_jobs=10,n_iter=10)
