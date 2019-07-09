import rna_io
from util import shexec
import re
import math

import tempfile
def get_ens_energy(seq,react=None):
    '''calculate  ensemble energy'''
    if type( react) == type(None):
        retcode,err,out = shexec("echo %s | RNAfold --noPS -p0" % seq)
    else:
        fname = tempfile._get_default_tempdir() + '/' + next(tempfile._get_candidate_names()) + "rea.tmp"
        rna_io.write_shape(fname, react)
        retcode,err,out = shexec("echo %s | RNAfold --noPS --shape %s -p0" % (seq,fname))
    # a float followed by kcal/mol
    return float(  re.findall( r"([-+]?[0-9]*\.?[0-9]+) kcal/mol",out )[0] )


def get_stru_energy(struct, sequence,react=None):
    """calculate energy of a structure"""
    if type( react) == type(None):
        cmd = "echo \"%s\n%s\" | RNAeval" % (sequence,struct)
    else:
        fname = tempfile._get_default_tempdir() + '/' + next(tempfile._get_candidate_names()) + "rea.tmp"
        rna_io.write_shape(fname, react)
        cmd = "echo \"%s\n%s\" | RNAeval --shape %s" % (sequence,struct,fname)
    retcode,err,out = shexec(cmd)
    return float(re.findall(r"[-+]?[0-9]*\.?[0-9]+",out)[0])


def energy_to_proba(ensemble,other):
    """use the obvious formula to calculate the probability of a structure given its energy and the energy of the ensemble"""
    RT= 0.61632
    #return math.exp(-other/RT) / math.exp(-ensemble/RT)
    try:
        #res = math.exp(-other/RT) / math.exp(-ensemble/RT)
        res = math.exp((-other+ensemble)/RT)
    except:
        print "energy to proba failed: other:%.100f ens: %.100f" % (other, ensemble)

    return res

def probability(structure,seq, react=None):
    """calc probabity of a structure given a sequence and optionaly reactivity data"""
    return energy_to_proba(get_ens_energy(seq,react),get_stru_energy(structure,seq,react))


def probabilities_of_structures(sequence, structure_list, react=None):
    """calculate probabilities of structures for a sequence , the readct option is never used because it distorts the ens energies
    returns [(dotbracket,probability),..]
    """
    ensemble_energy = get_ens_energy(sequence, react = react)
    energies = map(lambda x: get_stru_energy(x, sequence, react = react), structure_list)
    probabilities = map(lambda x:energy_to_proba(ensemble_energy, x), energies)
    #probabilities = normalize(probabilities, norm='l1').tolist()[0]
    return [(stru,proba) for stru,proba in zip(structure_list,probabilities)]
