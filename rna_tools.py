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




def rnashapes(sequence):
    """call rnashapes, return (dotbracket, energy)"""
    retcode,err,out = shexec("RNAshapes %s" % sequence)
    if retcode != 0:
        print "RNAshapes failed"
        return

    energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)
    shape =[ a.strip() for a in re.findall(r' [.()]+ ',out) ]
    return shape, energy


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

def get_struct_and_proba(seq, cutoff = 0.01):
    """uses rnashapes to determine structure representatives, for the representatives, probabilities are calculated,

    returns [(dotbracket,probability),..]
    """
    ensemble_energy = get_ens_energy(seq)
    structs, _ = rnashapes(seq)
    energies = map( lambda x: get_stru_energy(x,seq), structs)
    probabilities = map(lambda x:energy_to_proba(ensemble_energy, x), energies)
    #probabilities = normalize(probabilities, norm='l1').tolist()[0]
    return [(stru,proba) for stru,proba in zip(structs,probabilities) if proba >= cutoff ]
    #energy = get_energy(seq, structs[0])
    #print energy
    #cat ss cc | RNAfold -p0
    #cat ss cc | RNAeval

def test():
    """a test :D"""
    testseq= "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"
    testseq= "ggaaauaaUCGGAUGAAGAUAUGAGGAGAGAUUUCAUUUUAAUGAAACACCGAAGAAGUAAAUCUUUCAGGUAAAAAGGACUCAUAUUGGACGAACCUCUGGAGAGCUUAUCUAAGAGAUAACACCGAAGGAGCAAAGCUAAUUUUAGCCUAAACUCUCAGGUAAAAGGACGGAGaaaacaaaacaaagaaacaacaacaacaac"
    print get_struct_and_proba(testseq)


def shexec(cmd):
    '''
    takes cmd, the command
    returns (exit-code, stderr, stdout)
    '''
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, stderr = process.communicate()
    retcode = process.poll()
    return (retcode, stderr, output)


def shape(sequence):
    '''given a sequence,  call rnashapes with some hardcoded params and return [dotbracketstrings],[energies]'''
    retcode,err,out = shexec("RNAshapes -u -s -t 5 -c 10 %s | sort -n " % sequence)
    if retcode != 0:
        print "RNAshapes failed"
        return
    energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)
    shape =[ a.strip() for a in re.findall(r' [.()]+ ',out) ]
    return shape, energy

