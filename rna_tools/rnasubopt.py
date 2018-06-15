from rna_tools import shexec
import re


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