import simushape_nostruct as sn
import simushape as ss
import re
import math

def fold(sequence,react):
    with open('shap.tmp','w') as f:
        f.write(ss.format_shape("",react,noheader=True))
    cmd = 'echo "%s" | RNAfold --shape shap.tmp'  % sequence
    res = sn.shexec(cmd)[2]
    return res[res.find("\n")+1: res.find(" ")]





'''now we have a sequence,
1. use rnashapes to get the representatives
2. use rnafold to get the energies for these
3. calculate probas
4. weight annotation by proba
'''


def rnashapes(sequence):
    retcode,err,out = sn.shexec("RNAshapes %s" % sequence)
    if retcode != 0:
        print "RNAshapes failed"
        return

    energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)
    shape =[ a.strip() for a in re.findall(r' [.()]+ ',out) ]
    return shape, energy

def get_ens_energy(seq):
    retcode,err,out = sn.shexec("echo %s | RNAfold -p0" % seq)
    # a float followed by kcal/mol
    return float(  re.findall( r"([-+]?[0-9]*\.?[0-9]+) kcal/mol",out )[0] )

def get_stru_energy(struct, sequence):
    cmd = "echo \"%s\n%s\" | RNAeval" % (sequence,struct)
    retcode,err,out = sn.shexec(cmd)
    return float(re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)[0])

def energy_to_proba(ensemble,other):
    RT= 0.61632
    return math.exp(-other/RT) / math.exp(-ensemble/RT)

def get_struct_and_proba(seq):
    ensemble_energy = get_ens_energy(seq)
    structs, _ = rnashapes(seq)
    energies = map( lambda x: get_stru_energy(x,seq), structs)
    probabilities = map(lambda x:energy_to_proba(ensemble_energy, x), energies)
    print probabilities
    #energy = get_energy(seq, structs[0])
    #print energy
    #cat ss cc | RNAfold -p0
    #cat ss cc | RNAeval
    #return structs, probas

def test():
    testseq= "CCAUGAAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACCCCC"
    print get_struct_and_proba(testseq)

test()



