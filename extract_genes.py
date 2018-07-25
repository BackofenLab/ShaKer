import collections
path="/home/montaser/ShaKer/data/weeks194_orig/"
        
#for genes
def extractgenes():
    endIdxG = 0
    i = 0
    dbnfileGenes = open(path + "cellfreeGenes.dbn", "w") # we write s.th
    with open(path + "genes.txt", "r") as g:  # we read something
        g.next()
        for e in g:
            if "," in e:  # all the lines with "," are interesting
                name, Gstart, Gend, Tstart, Tend, idx =  e.strip().split(',')[0:6]
                assert (int(Gend) - int(Gstart)) == (int(Tend) - int(Tstart)) 
                endIdxG, i = extractTranscriptomes(idx, Tstart, Tend, dbnfileGenes, endIdxG, i)

def extractTranscriptomes(idx, Tstart, Tend, dbnfileGenes, endIdxG, i):
    with open(path + "cellfree.dbn", "r") as D:
        for d in D:
            if ">" in d:
                namesidx = d.strip().split('_')
                indx = namesidx[len(namesidx)-1]
                if indx == idx:
                    i = i + 1
                    seq = D.next()
                    if endIdxG == 0:
                        startIdxT = endIdxG
                    else: startIdxT = endIdxG + 1
                    startIdxG = int(Tstart)
                    endIdxG = int(Tend)
                    endIdxT = startIdxG
                    nameT = "T_" + idx + "_" + str(i)
                    if startIdxT == 0:
                        text =">%s\n%s\n" % (nameT, seq[startIdxT: endIdxT - 1])
                    else: text =">%s\n%s\n" % (nameT, seq[startIdxT - 1: endIdxT - 1])
                    dbnfileGenes.write(text)
                    if i == len(namesidx) - 1: 
                        text = ">%s\n%s\n" % ("T_" + idx + "_" + str(i+1), seq[endIdxG: ])
                        dbnfileGenes.write(text)
                        endIdxG = 0
                        i = 0
                    break

    return (endIdxG, i)

#for reactivities                    
def extractreacts():
    endIdxG = 0
    i = 0
    l = 0
    react = ""
    reactdict = {}
    dbnfileGenes = open(path + "cellfreereacts.react", "w")
    with open(path + "genes.txt", "r") as g:
        g.next()
        for e in g:
            if "," in e:
                data = e.strip().split(',')
                name, Gstart, Gend, Tstart, Tend, idx = data[0:6]
                assert (int(Gend) - int(Gstart)) == (int(Tend) - int(Tstart))
                l = extractLSeq(idx)
                endIdxG, i, reactdict = computereact(idx, Tstart, Tend, l, dbnfileGenes, endIdxG, i, reactdict)
                    
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

def computereact(idx, Tstart, Tend, l, dbnfileGenes, endIdxG, i, reactdict):
    with open(path + "cellfree.react", "r") as D:
        for d in D:
            if ">" in d:
                namesidx = d.strip().split('_')
                indx = namesidx[len(namesidx)-1]
                if indx == idx:
                    i = i + 1
                    for d in D:
                        indxreact = d.strip().split(' ')
                        if ">" not in d:
                            index, react = indxreact[0:2]
                            reactdict[int(index)] = react
                        else: break
                    reactdict = collections.OrderedDict(sorted(reactdict.items()))
                    if endIdxG == 0:
                        startIdxT = endIdxG
                    else: startIdxT = endIdxG + 1
                    startIdxG = int(Tstart)
                    endIdxG = int(Tend)
                    endIdxT = startIdxG
                    nameT = "T_" + idx + "_" + str(i)
                    text1= ""
                    if startIdxT == 0:
                        endT = endIdxT - 2
                    else: endT = endIdxT - 1

                    for k in range (startIdxT - 1, endT):
                        if startIdxT == 0: 
                            text1 = text1 + str(k - startIdxT + 2) + " " + reactdict.values()[k + 1] + "\n"
                        else: text1 = text1 + str(k - startIdxT + 2) + " " + reactdict.values()[k] + "\n"
                    text =">%s\n%s\n" % (nameT, text1)
                    dbnfileGenes.write(text)
                    if i == len(namesidx) - 1: 
                        text1 = ""
                        for k in range (endIdxG, l - 1):
                            text1 = text1 + str(k - endIdxG + 1) + " " + reactdict.values()[k] + "\n"
                        text = ">%s\n%s\n" % ("T_" + idx + "_" + str(i+1), text1)
                        dbnfileGenes.write(text)
                        endIdxG = 0
                        i = 0
                    break

    return (endIdxG, i, reactdict)


