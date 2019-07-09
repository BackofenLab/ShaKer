import util
import os
import re


cachefile = '.subopt_cache'

def rnasubopt(sequence, return_energy=False, samples=60):
        """call rnasubopt, return (dotbracket, energy)"""
        #print "echo " + "%s" % sequence + " | RNAsubopt -p 60 --maxBPspan=150"

        d={}
        if len(sequence) > 100:
            if os.path.isfile(cachefile):
                d = util.jloadfile(cachefile)
                if sequence in d:
                        
                    # hack for my adhoc cache system
                    if "(" in d[sequence][0] and return_energy ==False:
                        return d[sequence]
                    if "(" in d[sequence][0] and return_energy:
                        print "continuing would yield a problem"
                        exit()
                        
                    return d[sequence] if return_energy else d[sequence][0]
                    
            else:
                util.jdumpfile({},cachefile)
        retcode, err, out = util.shexec("echo " + "%s" % sequence + " | RNAsubopt -p %d --maxBPspan=150" % samples)
        energy = re.findall(r"[-+]?[0-9]*\.?[0-9]+", out)
        shape = [a.strip() for a in re.findall(r'\n[.()]+', out)]

        if len(sequence) > 100:
            d[sequence]= (shape, energy)
            util.jdumpfile(d,cachefile)

        if return_energy:
            return shape, energy
        else:
            return shape
