import subprocess
import rna_io
import re
import math
import collections
import rna_io
import os
from rna_tools import rna_accuracy, shexec
from os.path import isfile
from sklearn.preprocessing import normalize
import rna_tools.rna_accuracy  # is this even in the repo?


#################
# I HAVE NO IDEA WHAT THIS IS
###################

def soheilas_workaround(sequence):
    #sequence = "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"
    print "  I AM HERE"
    filename = 'asdasd.seq'
    with open(filename, "w") as f:
        f.write(">asdasd\n%s\n" % sequence)

    res = rnastructure("asdasd.seq", filename, 'test3.ct')
    print "RNAstructureeeeeeee", res
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

        #read from the file ctFile and save it as output
        #F = open("test2.ct", "r")
        #E = F.read()
        #energy = re.findall(r"[-+][0-9]*\.[0-9]+", E)

        with open(ctFile, "r") as C:
            array = []
            for line in C:
                array.append(line)

        l = len(array)
        column = []
        str = []
        k = 0
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
        #convertCtToDotbracket(E)
        #shape = [a.strip() for a in re.findall(r'\n[.()]+', out)]
        if return_energy:
            return shape
        else:
            return shape