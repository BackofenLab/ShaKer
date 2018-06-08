import subprocess
import rna_io
import re
import math
import collections
import rna_io
import os #Milads
import rna_accuracy

from os.path import isfile #Milads
#from sklearn.preprocessing import normalize #soheila commented it

import rna_accuracy


def fold(sequence,react = None):
    """call rna fold, return dotbracket string"""
    
    print "RNAfold"
    if react==None:
        cmd = 'echo "%s" | RNAfold -p '  % sequence
    else:
        with open('shap.tmp','w') as f:
            f.write(rna_io.format_shape("", react, noheader=True))
        cmd = 'echo "%s" | RNAfold --shape shap.tmp -p'  % sequence
    res = shexec(cmd)[2]

    #rna_accuracy.get_structure_accuracy(./notebooks/dot.ps, res[res.find("\n")+1: res.find(" ")]) #to compute accuracy
    #print "accuracy : ", rna_accuracy.get_structure_accuracy(['./notebooks/dot.ps'],[res[res.find("\n")+1: res.find(" ")]])
    
    print "accuracy : ", rna_accuracy.get_structure_accuracy(['./notebooks/dot.ps'],[[a.strip() for a in re.findall(r'\n[.()]+', res)]])
    
    #print "dot.ps", dot.ps
    return res[res.find("\n")+1: res.find(" ")]





def rnashapes(sequence, return_energy=False):
    """call rnashapes, return (dotbracket, energy)"""
    retcode,err,out = shexec("RNAshapes %s" % sequence)

    #retcode,err,out = shexec("RNAshapes -u -s -t 5 -c 10 %s | sort -n " % sequence)
    if retcode != 0:
        print "RNAshapes failed"
        return
    energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)
    shape =[ a.strip() for a in re.findall(r' [.()]+ ',out) ]
    if return_energy:
        return shape, energy
    else:
        return shape

    #soheila implemented rnasubopt


def rnasubopt(sequence, return_energy=False):
        """call rnasubopt, return (dotbracket, energy)"""
        #print "echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150"
        retcode, err, out = shexec("echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150")
        energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+", out)
        shape = [a.strip() for a in re.findall(r'\n[.()]+', out)]
        if return_energy:
            return shape, energy
        else:
            return shape



#Milad's stuff

def soheilas_plfold(sequence, react=None):
    #sequence = "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"
    #print "  I AM in RNAplfold"
    #print "sequence ", sequence
    resultlist = []
    if react==None:
        dotplot_file, accessibility_file = call_vienna_plfold(sequence, 'seq8')
        #dotplot_file, accessibility_file = call_vienna_plfold('ACGCGGCGAAGAGACGTTGATGAGCCTGCTTAGACCGATAACAC', 'seq1')
        #print "result ", \
        result = get_unpaired_probs(accessibility_file)
    else:
        dotplot_file, accessibility_file = call_vienna_plfold(sequence, 'seq8', react)
        #print "result ", \
        result = get_unpaired_probs(accessibility_file)

    result = collections.OrderedDict(sorted(result.items()))
    for key, value in result.iteritems():
        resultlist.append(value)
        #print "RNAplfold"
    return resultlist
        #get_unpaired_probs(accessibility_file)


def call_vienna_plfold(sequence, seq_name, react=None, W=200, L=150):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file
    Input: sequence as string of nucelutides
           seq_name as string to be used for intermediate and output file names
           :rtype: object
    '''

    from subprocess import Popen, PIPE
    dp_file_name = "{}_dp.ps".format(seq_name)
    unp_file_name = "{}_lunp".format(seq_name)
    if isfile(dp_file_name):  # Caution Race condition may occur
        os.remove(dp_file_name)
    if isfile(unp_file_name):  # Caution Race condition may occur
        os.remove(unp_file_name)

    #shapefile = open("tRNAp.txt", "r")
    if react == None:
        RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1'.format(W, L)  # -u 1 for unpaired probablitiy
    else:
        rna_io.write_shape('tmp.txt', react)
        RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 --shape tmp.txt --shapeMethod="D"'.format(W, L)  # -u 1 for unpaired probablitiy
        #RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 --shape tRNAp.txt'.format(W, L)  # -u 1 for unpaired probablitiy
    #print (sequence)
    assert len(sequence.split()) == 1
    cmd = ('echo -e ">%s\\n%s\\n" | ' % (seq_name, sequence))
    cmd += RNAPLFOLD
    #print(cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    #print (out)
    #print (err)
    has_error = False
    if err:
        if b"warn" in err.lower():
            print ("Warning in calling call_RNAfold:\n {} {}\n".format(out, err))
        else:
            raise RuntimeError("Error in calling call_RNAfold: {} {}\n".format(out, err))
    if not isfile(dp_file_name):
        raise RuntimeError("Error: Expected dp file: {} is not created!".format(dp_file_name))
    if not isfile(unp_file_name):
        raise RuntimeError("Error: Expected lunp file: {} is not created!".format(unp_file_name))

    return dp_file_name, unp_file_name


def get_unpaired_probs(unp_file):
    '''
    Reads Vienna RNAplfold unpaired prob file (-u 1) into dict
    Returns: a dictionary with 1-based positions as key and
             unpaired-probabilities, aka accessibilities, as values
    '''

    with open(unp_file) as unp_in:
        line1 = unp_in.readline()
        if "#unpaired probabilities" not in line1:
            raise IOError('Unexpected header for lunp file: {}'.format(line1))
        line2 = unp_in.readline()
        if "#i$\tl=1" not in line2:
            raise IOError('Unexpected second header for lunp file: {}'.format(line2))
        up_dic = dict()
        for line in unp_in:
            splits = line.split()
            assert len(splits) >= 2
            pos = int(splits[0])
            up_prob = float(splits[1])
            assert pos >= 1
            assert up_prob >= 0 and up_prob <= 1
            assert pos not in up_dic
            up_dic[pos] = up_prob

    return up_dic


# soheila implemented rnastructure
#USAGE: Fold <seq file> <ct file> [options]

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


def get_ens_energy(seq,react=None):
    '''calculate  ensemble energy'''
    if react==None:
        retcode,err,out = shexec("echo %s | RNAfold -p0" % seq)
    else:
        rna_io.write_shape('tmp.react',react)
        retcode,err,out = shexec("echo %s | RNAfold --shape tmp.react -p0" % seq)
    # a float followed by kcal/mol
    return float(  re.findall( r"([-+]?[0-9]*\.?[0-9]+) kcal/mol",out )[0] )


def get_stru_energy(struct, sequence,react=None):
    """calculate energy of a structure"""
    if react == None:
        cmd = "echo \"%s\n%s\" | RNAeval" % (sequence,struct)
    else:
        rna_io.write_shape('tmp.react',react)
        cmd = "echo \"%s\n%s\" | RNAeval --shape tmp.react" % (sequence,struct)
    retcode,err,out = shexec(cmd)
    return float(re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)[0])



def energy_to_proba(ensemble,other):
    """use the obvious formula to calculate the probability of a structure given its energy and the energy of the ensemble"""
    RT= 0.61632
    return math.exp(-other/RT) / math.exp(-ensemble/RT)


def probability(structure,seq, react=None):
    """calc probabity of a structure given a sequence and optionaly reactivity data"""
    return energy_to_proba(get_ens_energy(seq,react),get_stru_energy(structure,seq,react))


def probabilities_of_structures(sequence, structure_list):
    """calculate probabilities of structures for a sequence
    returns [(dotbracket,probability),..]
    """
    ensemble_energy = get_ens_energy(sequence)
    energies = map(lambda x: get_stru_energy(x, sequence), structure_list)
    probabilities = map(lambda x:energy_to_proba(ensemble_energy, x), energies)
    #print "probabilities_of_structures", probabilities
    #probabilities = normalize(probabilities, norm='l1').tolist()[0]
    return [(stru,proba) for stru,proba in zip(structure_list,probabilities)]


def test():
    """a test :D"""
    testseq= "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"
    testseq= "ggaaauaaUCGGAUGAAGAUAUGAGGAGAGAUUUCAUUUUAAUGAAACACCGAAGAAGUAAAUCUUUCAGGUAAAAAGGACUCAUAUUGGACGAACCUCUGGAGAGCUUAUCUAAGAGAUAACACCGAAGGAGCAAAGCUAAUUUUAGCCUAAACUCUCAGGUAAAAGGACGGAGaaaacaaaacaaagaaacaacaacaacaac"
    seqname = "/home/montaser/miniconda2/bin/ASE_00159.fa"
    testseq= "GAAUACAGCAGUUCUCUCUUCGUAAGAAGAGAGAUUGUGUUACCCUUCGUAACACAUUACCCGUGAGGGAACUAUUAGGUGGAAACCUAGUUAGCUUUGUCGUCAUAAGUGAUCAACGUCUUUCCCCUCGAUUUAUUUUCGGCGGGCGAAUGACGGAUGGGAACGGCGAUAGUUACAAUGUUUGACGUAUCACGUGAUAUUAUCACAAGGUAUGGUCUGUUGAGUGCAAUCGUAGGACCCCGGUUGUUAUCUACAACCGGUUGAA"
    testseq = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"

    #print "asdasdasd", rnashapes(testseq)
    #print "sadasd",probabilities_of_structures(testseq, rnashapes(testseq))

    #print rnasubopt(testseq) #soheila added
    #print rnastructure(seqname, "test2.ct") #soheila added that works correctly
    #print "RNAstructure111 ", probabilities_of_structures(testseq, rnastructure(seqname, "test2.ct"))

    #print "rnaplfold", rnaplfold(testseq)
    print "rnaplfold", soheilas_plfold(testseq)




def shexec(cmd):
    '''
    takes cmd, the command
    returns (exit-code, stderr, stdout)
    '''
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, stderr = process.communicate()
    retcode = process.poll()
    return (retcode, stderr, output)

#testseq = "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"

#soheilas_workaround(testseq)

#test()


