import simushape as ss


data=ss.quickladdata()
ss.opti_forest(data,n_jobs=24,n_iter=5000)
