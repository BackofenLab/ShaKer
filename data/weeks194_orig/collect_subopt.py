
import ShaKer.rna_tools.rna_io as rio

path ="/scratch/bi01/mautner/ShaKer/ShaKer" 

def getdata():
        return rio.get_all_data("%s/data/weeks194_orig/%s.react" % (path,"kasugamycin")
                ,"%s/data/weeks194_orig/%s.dbn" % (path,"kasugamycin"))
data = getdata() 


k = data.keys()

import sys
import ShaKer.rna_tools.rnasubopt as so
import ShaKer.rna_tools.util as util



path = '/home/mautner/JOBZ/asd'


res = {}
for i,e in enumerate(k[:5]):
    f = "/5297087.o_%d" % (i+1)
    with open(path+f,"r") as fi:
        res[data[e][1]]= eval(fi.read())

util.jdumpfile(res,'.subopt_cache') 



