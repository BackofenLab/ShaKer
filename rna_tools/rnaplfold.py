from rna_tools import rna_io


def soheilas_plfold(sequence, react=None):

    resultlist = []

    if react==None:
        dotplot_file, accessibility_file = call_vienna_plfold(sequence, 'seq8')
        #dotplot_file, accessibility_file = call_vienna_plfold('ACGCGGCGAAGAGACGTTGATGAGCCTGCTTAGACCGATAACAC', 'seq1')
        #print "result ", \
        result = get_unpaired_probs(accessibility_file)
    else:
        dotplot_file, accessibility_file = call_vienna_plfold(sequence, 'seq8', react)
        #print "result ", \
        result = get_unpaired_probs(accessibility_file)

    result = collections.OrderedDict(sorted(result.items()))
    for key, value in result.iteritems():
        resultlist.append(value)
        #print "RNAplfold"
    return resultlist
        #get_unpaired_probs(accessibility_file)


def call_vienna_plfold(sequence, seq_name, react=None, W=200, L=150):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file
    Input: sequence as string of nucelutides
           seq_name as string to be used for intermediate and output file names
           :rtype: object
    '''

    from subprocess import Popen, PIPE
    dp_file_name = "{}_dp.ps".format(seq_name)
    unp_file_name = "{}_lunp".format(seq_name)
    if isfile(dp_file_name):  # Caution Race condition may occur
        os.remove(dp_file_name)
    if isfile(unp_file_name):  # Caution Race condition may occur
        os.remove(unp_file_name)

    #shapefile = open("tRNAp.txt", "r")
    if react == None:
        RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1'.format(W, L)  # -u 1 for unpaired probablitiy
    else:
        rna_io.write_shape('tmp.txt', react)
        RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 --shape tmp.txt --shapeMethod="D"'.format(W, L)  # -u 1 for unpaired probablitiy
        #RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 --shape tRNAp.txt'.format(W, L)  # -u 1 for unpaired probablitiy
    #print (sequence)
    assert len(sequence.split()) == 1
    cmd = ('echo -e ">%s\\n%s\\n" | ' % (seq_name, sequence))
    cmd += RNAPLFOLD
    #print(cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    #print (out)
    #print (err)
    has_error = False
    if err:
        if b"warn" in err.lower():
            print ("Warning in calling call_RNAfold:\n {} {}\n".format(out, err))
        else:
            raise RuntimeError("Error in calling call_RNAfold: {} {}\n".format(out, err))
    if not isfile(dp_file_name):
        raise RuntimeError("Error: Expected dp file: {} is not created!".format(dp_file_name))
    if not isfile(unp_file_name):
        raise RuntimeError("Error: Expected lunp file: {} is not created!".format(unp_file_name))

    return dp_file_name, unp_file_name


def get_unpaired_probs(unp_file):
    '''
    Reads Vienna RNAplfold unpaired prob file (-u 1) into dict
    Returns: a dictionary with 1-based positions as key and
             unpaired-probabilities, aka accessibilities, as values
    '''

    with open(unp_file) as unp_in:
        line1 = unp_in.readline()
        if "#unpaired probabilities" not in line1:
            raise IOError('Unexpected header for lunp file: {}'.format(line1))
        line2 = unp_in.readline()
        if "#i$\tl=1" not in line2:
            raise IOError('Unexpected second header for lunp file: {}'.format(line2))
        up_dic = dict()
        for line in unp_in:
            splits = line.split()
            assert len(splits) >= 2
            pos = int(splits[0])
            up_prob = float(splits[1])
            assert pos >= 1
            assert up_prob >= 0 and up_prob <= 1
            assert pos not in up_dic
            up_dic[pos] = up_prob

    return up_dic