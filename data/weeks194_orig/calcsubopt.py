
import ShaKer.rna_tools.rna_io as rio

path ="/scratch/bi01/mautner/ShaKer/ShaKer" 

def getdata():
        return rio.get_all_data("%s/data/weeks194_orig/%s.react" % (path,"kasugamycin")
                ,"%s/data/weeks194_orig/%s.dbn" % (path,"kasugamycin"))
data = getdata() 


k = data.keys()

import sys
import ShaKer.rna_tools.rnasubopt as so
if __name__=="__main__":
    task = int(sys.argv[1])-1
    print so.rnasubopt(data[k[task]][1])




