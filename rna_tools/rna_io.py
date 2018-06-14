# i am ok with this :)

def read_dbn(path):
    '''returns (name sequcence dotbrack tupples)'''
    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res=[]
        for e in text:
            a_thing = [thing for thing in e.split('\n') if len(thing) > 1]
            if len(a_thing)!=3:
                print "ERRER", a_thing ,e
                return
            if len(a_thing[1])!=len(a_thing[2]):
                print "ERRER", a_thing ,e
                return
            res.append(a_thing)
    return res


def read_fasta(path):

    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res=[]
        for e in text:
            a_thing = [thing for thing in e.split('\n') if len(thing) > 1]
            if len(a_thing)!=2:
                print "ERRER", a_thing ,e
                return
            res.append(a_thing)
    return res


def float_or_none(item):
    if item == 'NA':
        return None
    return float(item)


def read_react(path):
    '''returns {seqname:reactlist}'''
    with open(path,'r') as fi:
        text = fi.read()
        text = text.split(">")[1:]
        res={}
        for e in text:
            lines = [thing for thing in e.split('\n') if len(thing) > 1]
            header = lines[0]
            data = [ e.strip().split() for e in lines[1:]  ]
            data = [ float_or_none(d[1].strip()) for d in data if len(d)==2 ]
            res[header.strip()] = data
    return res


def combine_dbn_react(dbn,react):
    '''returns {seqname:[react, sequence, dotbracketstring]}'''
    res = {}
    for name,seq,brack in dbn:
        re = react[name]
        if len(re) == len(seq) == len(brack):
            res[name]= (re,seq,brack)
        else:
            print "data for '%s' is corrupted, ignoring..." % name
    return res


def get_all_data(react, dbn):
    '''takes paths to .react and .dbn  returns {seqname:[react, sequence, dotbracketstring]}'''
    dbn = read_dbn(dbn)
    react = read_react(react)
    return combine_dbn_react(dbn,react)


def format_shape(name,data, noheader=False):
    '''takes shape data, returns .react string e.i. >seqname\nVal\nVal etc'''
    vout = lambda x: "\n".join( [str(i + 1) + "\t" + str(e) for i, e in enumerate(x)])
    if noheader:
        return "%s\n\n" % vout(data)
    return ">%s\n%s\n\n" % (name,vout(data))



def write_shape(fname,react):
    '''writes react file'''
    with open(fname,'w') as f:
        f.write(format_shape("", react, noheader=True))

def dump_shape(result, fname):
    '''writes many .react files'''
    with open(fname,'w') as f:
        for k,v in result.items():
            f.write(format_shape(k, v))

