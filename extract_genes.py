import collections
path="/home/montaser/ShaKer/data/weeks194_orig/"
        
#for genes
def extractgenes():
    endIdxG = 0 
    i = 0
    filew = openfile("cellfreeGenes.dbn", "w")
    filer = openfile("genes.txt", "r")
    extractfilerInfo(filer, filew, endIdxG, i)# extract information from read file
                
                
def openfile(filename, type):
    if type == "w":
        file = open(path + filename, "w") # write into a file
    else: file = open(path + filename, "r") # read a file
    return file 

def extractfilerInfo(filer, filew, endIdxG, i, reactdict = None):
    with filer as g:  # read readable file
        g.next()
        for e in g:
            if "," in e:  # all the lines with "," are interesting
                name, Gstart, Gend, Tstart, Tend, idx =  e.strip().split(',')[0:6]
                assert (int(Gend) - int(Gstart)) == (int(Tend) - int(Tstart))
                l = extractLSeq(idx)
                if reactdict == None: endIdxG, i = extractTranscriptomes(idx, Tstart, Tend, filew, endIdxG, i)
                else: endIdxG, i, reactdict = computereact(idx, Tstart, Tend, l, filew, endIdxG, i, reactdict)

def extractTranscriptomes(idx, Tstart, Tend, filew, endIdxG, i):
    cutoff = 20
    structure = ""
    filer = openfile("cellfree.dbn", "r")
    with filer as D:
        for d in D:
            if ">" in d: 
                namesidx = d.strip().split('_')# CONTAINS NAME AND INDEX
                if namesidx[len(namesidx)-1] == idx:
                    i = i + 1
                    seq = D.next()
                    if endIdxG == 0: # write the first transcript
                        startIdxT = endIdxG 
                        if len(seq[startIdxT: int(Tstart) - 1]) >= cutoff: 
                            structure = makeStructure(len(seq[startIdxT: int(Tstart) - 1]))
                            filew.write(">%s\n%s\n%s\n" % ("T_" + idx + "_" + str(i), seq[startIdxT: int(Tstart) - 1], structure))
                    else: # write the rest transcript
                        startIdxT = endIdxG + 1 
                        if len(seq[startIdxT - 1: int(Tstart) - 1])>= cutoff: 
                            structure = makeStructure(len(seq[startIdxT - 1: int(Tstart) - 1]))
                            filew.write(">%s\n%s\n%s\n" % ("T_" + idx + "_" + str(i), seq[startIdxT - 1: int(Tstart) - 1], structure))
                    endIdxG = int(Tend)
                    if i == len(namesidx) - 1: # write the last transcript
                        if len(seq[endIdxG: ]) - 1 >= cutoff: 
                            structure = makeStructure(len(seq[endIdxG: ]) -1)
                            filew.write(">%s\n%s%s\n" % ("T_" + idx + "_" + str(i+1), seq[endIdxG: ], structure))
                        endIdxG = 0
                        i = 0
                    break
    return (endIdxG, i)
                    
def makeStructure(lseq):
    struct = ""
    for k in range(lseq):
        struct = struct + "."
    return struct
   
    
#for reactivities                    
def extractreacts():
    endIdxG = 0
    i = 0
    l = 0
    react = ""
    reactdict = {}
    filew = openfile("cellfreereacts.react", "w")
    filer = openfile("genes.txt", "r")
    extractfilerInfo(filer, filew, endIdxG, i, reactdict)# extract information from read file
                 
                    
def extractLSeq(idx):
    with open(path + "cellfree.dbn", "r") as C:
        for c in C:
            if ">" in c:
                namesidx = c.strip().split('_')
                indx1 = namesidx[len(namesidx)-1]
                if indx1 == idx:
                    seq = C.next()
                    l =  len(seq)
    return l  

def computereact(idx, Tstart, Tend, l, filew, endIdxG, i, reactdict):
    cutoff = 20
    filer = openfile("cellfree.react", "r")
    with filer as D:
        for d in D:
            if ">" in d:
                namesidx = d.strip().split('_')# contains name and index
                if namesidx[len(namesidx)-1] == idx:
                    i = i + 1
                    for d in D:
                        if ">" not in d:
                            index, react = d.strip().split(' ')[0:2] # contains index and reactivity
                            reactdict[int(index)] = react# contains reactivity
                        else: break
                    reactdict = collections.OrderedDict(sorted(reactdict.items()))# sort reactdict
                    text1= ""
                    if endIdxG == 0: # Index of the first transcript
                        for k in range (endIdxG - 1, int(Tstart) - 2):
                            text1 = text1 + str(k - endIdxG + 2) + " " + reactdict.values()[k + 1] + "\n"
                    else: # Index of the rest transcripts
                        for k in range (endIdxG, int(Tstart) - 1):
                            text1 = text1 + str(k - endIdxG + 1) + " " + reactdict.values()[k] + "\n" 
                    if int(Tstart) - endIdxG - 1 >= cutoff:
                        filew.write(">%s\n%s\n" % ("T_" + idx + "_" + str(i), text1)) # write the transcripts    
                    endIdxG = int(Tend)
                    if i == len(namesidx) - 1:  # write the last transcript
                        text1 = ""
                        for k in range (endIdxG, l - 1):
                            text1 = text1 + str(k - endIdxG + 1) + " " + reactdict.values()[k] + "\n"
                        if l - endIdxG - 1 >= cutoff:
                            filew.write(">%s\n%s\n" % ("T_" + idx + "_" + str(i+1), text1))
                        endIdxG = 0
                        i = 0
                    break

    return (endIdxG, i, reactdict)

    


