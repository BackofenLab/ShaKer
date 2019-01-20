import ShaKer.rna_tools
import ShaKer.simushape as ss
import numpy as np
from graphlearn01.utils import draw
import matplotlib.pyplot as plt
import pylab as plt
import matplotlib.colors as colors
import ShaKer.rna_tools.rnafold as rnafold
import ShaKer.data.weeks194_orig.remove_genes as d
def getgenedict():
    return d.read_genes()


def annotate(g, shap):
    """ sets the col argument in the nodes of g, accorgin to shape data shap"""
    n = g.nodes()
    n.sort()
    cm= plt.get_cmap("viridis")
    for e, i in zip(n, shap):
        g.node[e]["col"] =  "#AAAAAA"   if i == None else colors.rgb2hex(cm(i))
        g.node[e]["none"] = ''
    return g


import eden.display as ed
def draw_print(seq, shape, file_name= 'lol.svg',stru=None):
    '''seq+shape => image'''
    if stru != None:
        brack = stru
    else:
        brack = rnafold.fold(seq)

    graph = ss.eden_rna.sequence_dotbracket_to_graph(seq,brack)
    graph.graph['structure']= brack
    graph = annotate(graph.copy(), shape)

    ed.draw_graph(graph, size=5, layout="RNA",vertex_color='col',
                    vertex_label=None, edge_alpha=0.4, vertex_size=150,
                    vertex_border=False, file_name=file_name)


def draw3(graph, shape_list, **kwargs):
    '''draws graph len(shape_list) times in a row, colored by shape data'''
    graphs = [annotate(graph.copy(), shape) for shape in shape_list]
    draw.graphlearn(graphs, size=15, layout="RNA",vertex_color='col',
                    vertex_label='none', edge_alpha=0.05, vertex_size=150,
                    vertex_border=False, n_graphs_per_line=max(len(shape_list),5),**kwargs)


def draw_seq_rea(sequence, react_list, stru=None, **kwargs):
    ''' given sequence, reactivuty and maybe a structure: draw graph row showing reactivity '''
    if stru != None:
        brack = stru
    else:
        print "FIXME, i should use rnafold here"
        brack = rna_tools.shape(sequence)[0][0]
    graph = ss.eden_rna.sequence_dotbracket_to_graph(sequence,brack)
    graph.graph['structure']= brack
    draw3(graph, react_list,**kwargs)


def reactivitiesBoxplot(cross_predictions_list, data, keys, fig_title, xlabels, label_title, nr=5, nc=3):     
    '''i assume that this draws a grid of boxplots.'''
    j=0
    k=0
    fig, axes = plt.subplots(nrows=nr, ncols=nc, figsize=(20, 30))
    for i, key in enumerate(keys):
        if i==0: k=0
        elif (i!=0)&(i%3==0):
            j=j+1
            k=0
        else: k=k+1
        mydata = [x for x in data[key][0] if x != None]
        #if Sukosd !=None:
            #ddata=[mydata, cross_predictions_list[key]]
        #else: 
        ddata=[mydata, cross_predictions_list[i]]
        axes[j, k].boxplot(ddata, labels = xlabels)
        axes[j, k].set_title(label_title + '(' + fig_title + ') for ' + key + 'length of' + str(len(cross_predictions_list[i]))) #Comparison of reactivities of RealShape and SimuShape (computed by ' + fig_title + ') for ' + key
    plt.show()
    fig.savefig(fig_title + '-Reactivities.png', bbox_inches='tight')

def boxplotDraw(corrData, axis, plt_title, fig_title, filesave):
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(1, 1, 1)
        major_ticks = np.arange(0, 1, 0.05)
        minor_ticks = np.arange(0, 1, 0.05)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        # And a corresponding grid
        ax.grid(which='both')
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.2)
        plt.boxplot(corrData, 0, 'rs', 1)
        if len(axis) == 2:
            plt.scatter([1,2], [np.mean(x) for x in corrData] )
        elif len(axis) == 3: 
            plt.scatter([1,2,3], [np.mean(x) for x in corrData] )
        plt.xticks([y+1 for y in range(len(corrData))], axis)
        plt.xlabel(fig_title)
        plt.ylim([0, 1])
        t = plt.title(plt_title)
        fig.savefig(filesave + '.png', bbox_inches='tight')
        plt.show()







def get_genetrack(sequencename, data, genestartdict, low=0.1, high = 0.2):
    # name is not understood by the genedirectory
    index = int(sequencename.split("_")[-1])
    length = len(data[sequencename][1])
    genes = genestartdict[index]

    genetrack = [low]*length
    for start, end in genes:
        print start,end
        for e in range(start,end):
            genetrack[e] = high
    return genetrack



def get_genetrack_multigene(sequencename, data, genestartdict, drawindex=(0,9999999), low=0.1, high = 0.1):
    # returns genenames for the legend, value arrays for the line plot



    # get gene location list
    index = int(sequencename.split("_")[-1])
    length = len(data[sequencename][1])
    genes = genestartdict[index]

    # get  gene name ist
    names = sequencename.split("_")[:-1]
    relevant_names = []

    # background
    genetrack = []

    for (start,end),name in zip(genes,names):
        if end < drawindex[0] or start > drawindex[1]:
            continue
        nutrack = [np.nan]*length
        for e in range(start,end):
            nutrack[e]=high
        genetrack.append(nutrack[drawindex[0]:drawindex[1]])

        relevant_names.append("Gene: "+name)

    return relevant_names, genetrack























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
