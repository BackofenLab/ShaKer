import collections
import sys
sys.path.insert(0,'../..')
datapath="/home/ikea/Mustoe2018_data/"

##################
# phase 1 get data
###################


from collections import defaultdict
import rna_tools.rna_io as rio
import rna_tools.rnafold as rf

def read_genes():
    '''return {transcript_id:[(genestart,stop),...]}'''
    res=defaultdict(list)
    with open(datapath+"genes.txt", "r") as g:
        g.next()
        for e in g:
            assert "," in e
            name, Gstart, Gend, Tstart, Tend, idx =  e.strip().split(',')
            res[int(idx)].append((int(Tstart),int(Tend))) 
    return res

def read_dbn(fname):
    '''return sorted dbn ;sort because we will zip with the gene keys'''
    dbn = rio.read_dbn(fname)
    dbn.sort( key=lambda x:int(x[0].split("_")[-1]))
    return dbn


##############
# phase 2 mixing data
##############

def removegenes(gen, dbn, rea, lencutoff=10):
    ''' returns {seqname:react} and [seqname, seq]'''
    res_react = {}
    res_dbn = []
    k=gen.keys()
    k.sort()
    for kk, (seqname,seq,_) in zip (k,dbn): # for all transcripts
        genelocations = [(0,0)]+gen[kk] # fake gene for edge case
        genelocations.reverse()
        react = rea[seqname]
        for i, (start, stop )in enumerate(genelocations): # for all genes
            myreact = react[stop:]
            react = react[:start-1]
            myseq = seq[stop:]
            seq = seq[:start-1]
            if len(myseq) > lencutoff:
                sequencename = "%s_%d" % (seqname, i)
                res_react[sequencename] = myreact
                res_dbn.append([sequencename,myseq,rf.fold(myseq,myreact)])

    return res_react, res_dbn

##############################
# phase 3 dump data 
###########################


for cdata in ['incell','kasugamycin', 'cellfree']:
    gen = read_genes()
    dbn = read_dbn(cdata+".dbn")
    rea = rio.read_react(cdata+".react")
    react_dict, fasta_list = removegenes(gen,dbn,rea)
    rio.dump_shape(react_dict, cdata+"_nogenes.react")
    rio.dump_dbn(fasta_list, cdata+"_nogenes.dbn")



