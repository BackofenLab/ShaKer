import util
import os
import re



def rnasubopt(sequence, return_energy=False):
        """call rnasubopt, return (dotbracket, energy)"""
        #print "echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150"

        d={}
        if len(sequence) > 100:
            if os.path.isfile(".rnasubopt_cache"):
                d = util.loadfile(".rnasubopt_cache")
                if sequence in d:
                    return d[sequence] if return_energy else d[sequence][0]
            else:
                util.dumpfile({}, '.rnasubopt_cache')
        retcode, err, out = util.shexec("echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150")
        energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+", out)
        shape = [a.strip() for a in re.findall(r'\n[.()]+', out)]

        if len(sequence) > 100:
            d[sequence]= (shape, energy)
            util.dumpfile(d, ".rnasubopt_cache")

        if return_energy:
            return shape, energy
        else:
            return shape
