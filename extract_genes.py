


def extractgenes():
    path="/home/montaser/ShaKer/data/weeks194_orig/"
    endIdxG = 0
    i = 0
    dbnfileGenes = open(path + "cellfreeGenes.dbn", "w")
    with open(path + "genes.txt", "r") as g:
        g.next()
        for e in g:
            print "line in genes"
            if "," in e:
                data = e.strip().split(',')
                #print "data ", data
                name, Gstart, Gend, Tstart, Tend, idx = data[0:6]
                print "name, Gstart, Gend, Tstart, Tend, idx", name, Gstart, Gend, Tstart, Tend, idx
                print "idx ", idx
                if (int(Gend) - int(Gstart)) == (int(Tend) - int(Tstart)):
                    with open(path + "cellfree.dbn", "r") as D:
                        for d in D:
                            if ">" in d:
                                namesidx = d.strip().split('_')
                                #print "len(namesidx) ", len(namesidx)
                                #print "namesidx ", namesidx
                                indx = namesidx[len(namesidx)-1]
                                #print "indx ", indx
                                if indx == idx:
                                    i = i + 1
                                    seq = D.next()
                                    #print "seq", seq
                                    if endIdxG == 0:
                                        startIdxT = endIdxG
                                    else: startIdxT = endIdxG + 1
                                    #print int(Gstart),int(Tstart)
                                    startIdxG = int(Tstart)
                                    #print "startIdxG ", startIdxG
                                    endIdxG = int(Tend)
                                    endIdxT = startIdxG
                                    print "startIdxG,  endIdxG ", startIdxG,  endIdxG
                                    print "startIdxT,  endIdxT ", startIdxT,  endIdxT
                                    nameT = "T_" + idx + "_" + str(i)
                                    text =">%s\n%s\n" % (nameT, seq[startIdxT: endIdxT])
                                    print "text ", text
                                    dbnfileGenes.write(text)
                                    if i == len(namesidx) - 1: 
                                        text = ">%s\n%s\n" % ("T_" + idx + "_" + str(i+1), seq[endIdxG + 1: ])
                                        print "text ", text
                                        dbnfileGenes.write(text)
                                        endIdxG = 0
                                        i = 0
                                    break
                                #else: 
                                    #endIdxG = 0
                                    #i = 0
                                   
                                    




'''
        array = []
        for line in D:
            array.append(line)
    l = len(array)
    str = []
    for i in range(l):
        if (i % lSeq) != 0:
            column = array[i].split()
            if column[4] == '0':
                str.append('.')
            elif int(column[0]) < int(column[4]):
                str.append('(')
            else:
                str.append(')')


    filename_to_genname = {}
    with open(prefix + "transcripts.txt", "r") as f:
        f.next()  # first line is a table header
        for e in f:
            if "," in e:
                data = e.split(',')
                idx, _, start, end = data[:4]
                name = "_".join(data[4:]).strip() + "_" + idx
                filename_to_genname[start + '-' + end] = name

    transcriptcfile = open(path + "transcript.txt", "r")
    genesfile = open(path + "genes.txt", "r")
    dbnfileGenes = open(path + "cellfreeGenes.dbn", "w")
    dbnfileGenes.write("asd")
    readallline
#extract sequences. the structure are computed using RNAfold

def rnastructure_wrap(sequence):
    filename = 'asdasd.seq'
    with open(filename, "w") as f:
        f.write(">asdasd\n%s\n" % sequence)
    res = rnastructure("asdasd.seq", filename, 'test3.ct')
    return res

'''

