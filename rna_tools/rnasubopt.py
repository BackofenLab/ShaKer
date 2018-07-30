from rna_tools import shexec
import os
import re


import rna_tools
def rnasubopt(sequence, return_energy=False):
        """call rnasubopt, return (dotbracket, energy)"""
        #print "echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150"

        d={}
        if False: # deactivated cache...
            if len(sequence) > 300:
                if os.path.isfile(".rnasubopt_cache"):
                    d = rna_tools.loadfile(".rnasubopt_cache")
                    if sequence in d:
                        return d[sequence] if return_energy else d[sequence][0]
                else:
                    rna_tools.dumpfile({}, '.rnasubopt_cache')

        retcode, err, out = shexec("echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150")
        energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+", out)
        shape = [a.strip() for a in re.findall(r'\n[.()]+', out)]

        if len(sequence) > 300:
            d[sequence]= (shape, energy)
            #rna_tools.dumpfile(d,".rnasubopt_cache")

        if return_energy:
            return shape, energy
        else:
            return shape
