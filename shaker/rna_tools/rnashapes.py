from util import shexec
import re


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