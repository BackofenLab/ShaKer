from rna_tools import rna_io
import rna_tools
import simushape_nostruct  as sn
import simushape as ss
import numpy as np
from graphlearn01.utils import draw

import pylab as plt
import matplotlib.colors as colors

def annotate(g, shap):
    """ sets the col argument in the nodes of g, accorgin to shape data shap"""
    n = g.nodes()
    n.sort()
    cm= plt.get_cmap("viridis")
    for e, i in zip(n, shap):
        g.node[e]["col"] =  "#AAAAAA"   if i == None else colors.rgb2hex(cm(i))
        g.node[e]["none"] = ''
    return g


def draw3(graph, shape_list):
    '''draws graph len(shape_list) times in a row, colored by shape data'''
    graphs = [annotate(graph.copy(), shape) for shape in shape_list]
    draw.graphlearn(graphs, size=15, layout="RNA",vertex_color='col', vertex_label='none', edge_alpha=0.05, vertex_size=150,
                    vertex_border=False, n_graphs_per_line=len(shape_list))


def draw_seq_rea(sequence, react, stru=None):
    ''' given sequence, reactivuty and maybe a structure: draw graph row showing reactivity '''
    if stru != None:
        brack = stru
    else:
        print "FIXME, i should use rnafold here"
        brack = rna_tools.shape(sequence)[0][0]
    graph = ss.eden_rna.sequence_dotbracket_to_graph(sequence,brack)
    graph.graph['structure']= brack
    draw3(graph,react)









def get_test_data():
    print "removed get_test_data"
    '''
    data = rna_io.get_all_data('../data/RNA16.react', '../data/RNA16.dbn')
    #data2 = ss.get_all_data('data/RNA20.react','data/RNA20.dbn')
    data.pop("GLYCFN",None) # data is bad
    data.pop("23sRNA",None) #2 long for RNAshapes
    data.pop("R009",None)   #2 long for shapes

    data2 = rna_io.get_all_data('../data/RNA20.react', '../data/RNA20.dbn')
    for e in ['P546', 'Adenine', 'tRNA-phe', 'M-Box', 'tRNA-asp','5srRNA' ]:
        data[e] = data2[e]
    return data
    '''

def get_train_data():
    print "removed get train data"
    return get_test_data()





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



def run2(rs_structure=True, filename='' ,  maxstruct_predict=3, train_data=None, test_data=None):
    '''
    :param rs_structure:
    :param filename:
    :param maxstruct_train:
    :param maxstruct_predict:
    :return:

    the idea is to split train and test data..
    '''

    if test_data==None:
        test_data=get_test_data()
    if train_data==None:
        train_data= get_train_data()

    res=[]
    for e in test_data:
        train = train_data.keys()
        if e in train: train.remove(e)
        model = ss.make_model(train_data,train)

        target_sequence = test_data[e][1]
        target_struct = test_data[e][2]

        if rs_structure:
            my_react = np.array(sn.predict(model,target_sequence,maxstruct=maxstruct_predict))
        else:
            my_react = np.array(ss.predict2(model,target_sequence,target_struct))

        print 'x',
        res.append(">%s"%e)
        res.append('\n'.join(["%s\t%.4f" % (i,e) for i,e in enumerate(my_react)]))
        res.append('')

    with open(filename,'w') as f: f.write('\n'.join(res))
