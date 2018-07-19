#http://www.chem.unc.edu/rna/data-files/mustoe_2018_DATA_SOFTWARE.zip   


import sys
sys.path.append("../..")
import rna_tools.rnafold
from rna_tools import rna_io
import collections

prefix ="/home/montaser/ikea/mustoe_2018_DATA_SOFTWARE/Mustoe2018_data/"

filename_to_genname = {}
with open(prefix+"transcripts.txt","r") as f:
    f.next() # first line is a table header
    for e in f:
        if "," in e:
            data = e.split(',')
            idx,_,start,end = data[:4]
            name = "_".join(data[4:]).strip()+"_"+idx
            filename_to_genname[start+'-'+end] = name
    
def process_folder(foldername):
    allreacts = []
    allsequences = []
    allheaders = []
    allstructures = []
    for location in filename_to_genname:
        with open(prefix+"%s_SHAPE/%s.shape" % (foldername,location),"r") as ff:
            react = []
            sequence = ""
            for e in ff:
                if " " in e:
                    idx, react_str, dontknow, base =  e.strip().split(' ')
                    sequence+=base
                    react.append(react_str)
            allreacts.append(react)
            allsequences.append(sequence)
            allheaders.append(filename_to_genname[location])
        
        #strutures    
        with open(prefix+"%s_structures/%s.ct" % (foldername,location),"r") as ss:
            structure = ""
            ss.next()
            for s in ss:
                if " " in s:
                    idx1, base, indx0, indx2, bp, indx11 =  s.strip().split()
                    if int(bp) == 0: structure+="." 
                    elif int(idx1) < int(bp): structure+="(" 
                    else: structure+=")" 
            allstructures.append(structure)
            
    
    with open(foldername+".react","w") as f:
        transform = lambda x: "NA" if x == "nan" else x
        text=''
        for header,values in zip(allheaders,allreacts):
            sss  = '\n'.join( [ "%d %s" % (i+1, transform(v)) for i,v in enumerate(values)   ]  )
            text+=">%s\n%s\n" % (header, sss)

        f.write(text)    

    
    with open(foldername+".dbn","w") as f:
        text = ''
        react = rna_io.read_react(foldername + ".react")
        for header, sequence, structure in zip(allheaders, allsequences, allstructures):
            #db = rna_tools.rnafold.fold(sequence, react[header])
            text += ">%s\n%s\n%s\n\n" % (header, sequence, structure)
        f.write(text)
        

process_folder("incell")
process_folder("cellfree")
process_folder("kasugamycin")
