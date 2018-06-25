import rna_io
import os
import collections




from rna_tools import shexec


def rnaplfold(sequence, react=None):

    dotplot_file, accessibility_file = call_vienna_plfold(sequence, react)
    result = get_unpaired_probs(accessibility_file)
    return result



def call_vienna_plfold(sequence, react=None, W=200, L=150):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file
    Input: sequence as string of nucelutides
           seq_name as string to be used for intermediate and output file names
           :rtype: object
    '''
    seq_name='plfold'
    assert len(sequence) == len(react), "len seq and len react are not the same %d == %d" %(len(sequence),len(react))
    dp_file_name = "{}_dp.ps".format(seq_name)
    unp_file_name = "{}_lunp".format(seq_name)
    if os.path.isfile(dp_file_name): os.remove(dp_file_name)
    if os.path.isfile(unp_file_name): os.remove(unp_file_name)

    if react == None:
        RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1'.format(W, L)  # -u 1 for unpaired probablitiy
    else:
        rna_io.write_shape('tmp.txt', react)
        RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 --shape tmp.txt --shapeMethod="D"'.format(W, L)  # -u 1 unpaired proba
    cmd = ('echo  ">%s\\n%s\\n" | ' % (seq_name, sequence))
    cmd += RNAPLFOLD

    ret, error, out =  shexec(cmd)
    if ret != 0:
        print "error:", error
    assert os.path.isfile(dp_file_name), "Error: Expected dp file: {} is not created!".format(dp_file_name)
    assert os.path.isfile(unp_file_name),"Error: Expected lunp file: {} is not created!".format(unp_file_name)

    return dp_file_name, unp_file_name


def get_unpaired_probs(unp_file):
    '''
    Reads Vienna RNAplfold unpaired prob file (-u 1) into dict
    Returns: access list
    '''

    with open(unp_file) as unp_in:
        line1 = unp_in.readline()
        if "#unpaired probabilities" not in line1:
            raise IOError('Unexpected header for lunp file: {}'.format(line1))
        line2 = unp_in.readline()
        if "#i$\tl=1" not in line2:
            raise IOError('Unexpected second header for lunp file: {}'.format(line2))

        res = []
        for i, line in enumerate(unp_in):

            # there are 2 parts
            splits = line.split()
            assert len(splits) >= 2
            # is the possition arg as expeceted?
            assert (   int(splits[0])  == i+1 )

            # nice
            res.append( float(splits[1]) )

    return res
