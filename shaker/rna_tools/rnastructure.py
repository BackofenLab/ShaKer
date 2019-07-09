from util import shexec


# this wraps the tool: rnastructure by soheila montaseri
# lol what is this???

def rnastructure_wrap(sequence):
    filename = 'asdasd.seq'
    with open(filename, "w") as f:
        f.write(">asdasd\n%s\n" % sequence)
    res = rnastructure("asdasd.seq", filename, 'test3.ct')
    return res


def rnastructure(seqname, fastafile, ctFile, return_energy=False):
        """call rnasubopt, return (dotbracket, energy)"""
        print "Fold " + "%s" % seqname + " " + "%s" % ctFile
        retcode, err, out = shexec("Fold " + "%s" % seqname + " " + "%s" % ctFile)

        with open(fastafile, "r") as S:
            arraySeq = []
            for lineS in S:
                arraySeq.append(lineS)
        lSeq= len(arraySeq[1])

        with open(ctFile, "r") as C:
            array = []
            for line in C:
                array.append(line)

            #array = C.read().split("\n")

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
            elif i != 0:
                str.append(',')

        strSeq = ''.join(str)
        shape = strSeq.split(",")

        if return_energy:
            return shape
        else:
            return shape