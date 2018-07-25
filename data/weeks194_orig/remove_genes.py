import collections
import sys
sys.path.insert(0,'../../rna_tools')
datapath="/home/ikea/Mustoe2018_data/"

##################
# phase 1 get data
###################


from collections import defaultdict
import rna_io as rio

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
            


##############
# phase 2 mixing data
##############

def removegenes(gen, dbn, rea, lencutoff=20):
    ''' returns {seqname:react} and [seqname, seq]'''
    res_react = {}
    res_dbn = []
    k=gen.keys()
    k.sort()
    for kk, (seqname,seq,_) in zip (k,dbn): # for all transcripts
        genelocations = gen[kk] + [(0,0)] # fake gene for edge case
        genelocations.reverse()
        react = rea[seqname]
        
        for i, (start, stop )in enumerate(genelocations): # for all genes
            myreact = react[stop:]
            react = react[:start]
            myseq = seq[stop:]
            seq = seq[:start]
            if myseq > lencutoff:
                res_react[seqname+str(i)] = myreact
                res_dbn.append([seqname+str(i),myseq])

    return res_react, res_dbn

##############################
# now there should be a way to dump le data, but that seems trivial, also there might be off-by-1-errors here
###########################

cdata = 'incell'            

gen = read_genes()
dbn = rio.read_dbn(cdata+".dbn")
rea = rio.read_react(cdata+".react")
react_dict, fasta_list = removegenes(gen,dbn,rea)

