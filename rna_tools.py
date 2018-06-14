####
#  strange code at bottom, plfold code looks fishy... should be for the accessibility.. maybe move access and accuracy stuff in their own file
##
import subprocess
import rna_io
import re
import math
from sklearn.preprocessing import normalize


def fold(sequence,react = None):
    """call rna fold, return dotbracket string"""
    if react==None:
        cmd = 'echo "%s" | RNAfold '  % sequence
    else:
        with open('shap.tmp','w') as f:
            f.write(rna_io.format_shape("", react, noheader=True))
        cmd = 'echo "%s" | RNAfold --shape shap.tmp'  % sequence
    res = shexec(cmd)[2]
    return res[res.find("\n")+1: res.find(" ")]

def rnasubopt(sequence, return_energy=False):
    """call rnasubopt, return (dotbracket, energy)"""
    print "echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150"
    retcode, err, out = shexec("echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150")
    energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+", out)
    shape = [a.strip() for a in re.findall(r'\n[.()]+', out)]
    if return_energy:
        return shape, energy
    else:
        return shape


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


def get_ens_energy(seq,react=None):
    '''calculate  ensemble energy'''
    if react==None:
        retcode,err,out = shexec("echo %s | RNAfold -p0" % seq)
    else:
        rna_io.write_shape('tmp.react',react)
        retcode,err,out = shexec("echo %s | RNAfold --shape tmp.react -p0" % seq)
    # a float followed by kcal/mol
    return float(  re.findall( r"([-+]?[0-9]*\.?[0-9]+) kcal/mol",out )[0] )

############
# this block is just proba of structures to omou
#####
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

def probabilities_of_structures(sequence, structure_list, cutoff = 0.01):
    """calculate probabilities of structures for a sequence
    returns [(dotbracket,probability),..]
    """
    ensemble_energy = get_ens_energy(sequence)
    energies = map(lambda x: get_stru_energy(x, sequence), structure_list)
    probabilities = map(lambda x:energy_to_proba(ensemble_energy, x), energies)
    #probabilities = normalize(probabilities, norm='l1').tolist()[0]
    return [(stru,proba) for stru,proba in zip(structure_list,probabilities) if proba >= cutoff ]
def test():
    """a test :D"""
    testseq= "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"
    testseq= "ggaaauaaUCGGAUGAAGAUAUGAGGAGAGAUUUCAUUUUAAUGAAACACCGAAGAAGUAAAUCUUUCAGGUAAAAAGGACUCAUAUUGGACGAACCUCUGGAGAGCUUAUCUAAGAGAUAACACCGAAGGAGCAAAGCUAAUUUUAGCCUAAACUCUCAGGUAAAAGGACGGAGaaaacaaaacaaagaaacaacaacaacaac"
    print probabilities_of_structures(testseq, rnashapes(testseq))

######################

def shexec(cmd):
    '''
    takes cmd, the command
    returns (exit-code, stderr, stdout)
    '''
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, stderr = process.communicate()
    retcode = process.poll()
    return (retcode, stderr, output)




####
## PLfold wrapper ,
####
import os
from os.path import isfile


def call_vienna_plfold(sequence, seq_name, W=200, L=150):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file
    Input: sequence as string of nucelutides
           seq_name as string to be used for intermediate and output file names
    '''

    from subprocess import Popen, PIPE
    dp_file_name = "{}_dp.ps".format(seq_name)
    unp_file_name = "{}_lunp".format(seq_name)
    if isfile(dp_file_name):  # Caution Race condition may occur
        os.remove(dp_file_name)
    if isfile(unp_file_name):  # Caution Race condition may occur
        os.remove(unp_file_name)

    RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 '.format(W, L)  # -u 1 for unpaired probablitiy
    print (sequence)
    assert len(sequence.split()) == 1
    cmd = ('echo ">%s\\n%s\\n" | ' % (seq_name, sequence))
    cmd += RNAPLFOLD
    print(cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    # print (out)
    # print (err)
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


dotplot_file, accessibility_file = call_vienna_plfold('ACGCGGCGAAGAGACGTTGATGAGCCTGCTTAGACCGATAACAC', 'seq1')
get_unpaired_probs(accessibility_file)
