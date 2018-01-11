#####
# ARGPARSE
#####
import argparse
parser = argparse.ArgumentParser(description='shape annotation')
parser.add_argument('--bpseq',  type=str, default='BPseq', help='path to folder containing .bpseq files')
parser.add_argument('--seq',  type=str,  default='seq', help='path to folder containing .seq files')
parser.add_argument('--shape',  type=str, default='shape', help='path to folder containing .shape files ')
args = parser.parse_args()

DEBUG = True
if DEBUG:
    print 'argparse args: ', args



########
# LOAD WEIRDLY FORMATED DATA
########
from os import listdir
from os.path import isfile, join

def getfiles(mypath):
    files =  [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files.sort()
    files = sorted(files, key=lambda s: s.lower())
    def load(f):
        with open(join(mypath,f),'r') as fi:
            return fi.read()
    return [load(f) for f in files], files

seqs, _= getfiles(args.seq)
bpseq,_ = getfiles(args.bpseq)
shape,_= getfiles(args.shape)


########
# dump data, sane format
##########

def dump(fname,content):
    with open(fname,'w') as f:
        f.write(content)


def dump_fasta(se,dp,stringname):
    stringname = [ e[:e.find(".")] for e in stringname ]
    textlist = [">%s\n%s\n%s\n\n" % (s.strip(),d.strip(),n.strip())  for (s,d,n) in zip(stringname,se,dp)]
    res =  ''.join(textlist)
    dump("data.dbn",res)
    print res

def dump_shape(shape,fname):
    print shape
    fnames = [ e[:e.find(".")] for e in fname ]
    fix_data = lambda x: "\n".join([ str(i+1)+"\t"+e.strip() for i,e in enumerate( x.split("\n") ) if e.strip()>0  ])
    dataz = map(fix_data, shape)
    textlist = [">%s\n%s\n\n" % (name,data)  for (name,data) in zip(fnames,dataz)]
    res = ''.join(textlist)
    dump("data.react",res)
    print res


dump_fasta(seqs,bpseq,_)
dump_shape(shape,_)