




import simushape_nostruct  as sn
import simushape as ss
import numpy as np


def get_test_data():
    data = ss.get_all_data('data/RNA16.react','data/RNA16.dbn')
    #data2 = ss.get_all_data('data/RNA20.react','data/RNA20.dbn')
    data.pop("GLYCFN",None) # data is bad
    data.pop("23sRNA",None) #2 long for RNAshapes
    data.pop("R009",None)   #2 long for shapes
    return data


def run(rs_model=False,rs_structure=True, filename='' , maxstruct_train=3, maxstruct_predict=3):
    data=get_test_data()
    res=[]
    for e in data:
        train = data.keys()
        train.remove(e)
        if rs_model:
            model = sn.make_model(data,train,maxstruct=maxstruct_train)
        else:
            model = ss.make_model(data,train)

        target_sequence = data[e][1]
        target_struct = data[e][2]

        if rs_structure:
            my_react = np.array(sn.predict(model,target_sequence,maxstruct=maxstruct_predict))
        else:
            my_react = np.array(ss.predict2(model,target_sequence,target_struct))

        print 'x',
        res.append(">%s"%e)
        res.append('\n'.join(["%s\t%.4f" % (i,e) for i,e in enumerate(my_react)]))
        res.append('')

    with open(filename,'w') as f: f.write('\n'.join(res))

