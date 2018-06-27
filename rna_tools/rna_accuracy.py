
import rnaplfold
import os, re
import numpy as np

def dotbracket_to_dict(struct):
    '''Returns a dictionary where basepairs are keys with !ONE! based indices joined by ":" ,
    e.g. dict {'0:10': 1, '2:8': 1} '''
    assert len(struct.replace('.', '').replace('(', '').replace(')', '')) == 0
    stack = list()
    pairs = dict()
    for pos, ch in enumerate(list(struct)):
        #         print pos+1, ch
        if ch == '(':
            stack.append(pos)
        elif ch == ')':
            left = stack.pop()
            key = "{}:{}".format(left+1, pos+1)
            pairs[key] = 1

    assert len(stack) == 0
    return pairs


def parse_dp_ps(ps_file):
    '''Extracts base pair probabliies from vienna ps file
    returns: Dictinary of form dict[i:j]=p(i,j) '''

    # Extract sequence from ps file
    myseq = ""
    read_seq = False
    with open(ps_file) as in_ps:
        for line in in_ps:
            if "/sequence" in line:
                read_seq = True
            elif read_seq and ") } def" in line:
                read_seq = False
            elif read_seq:
                myseq += line.rstrip().rstrip("\\")
    #     print ps_file.rstrip("_dp.ps") , myseq

    ureg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+ubox\s*')
    bp_prob_dict = dict()
    bp_prob_mat = np.zeros((len(myseq), len(myseq)))

    with open(ps_file) as in_ps:
        for line in in_ps:
            if "ubox" in line:
                um = ureg.match(line)
                if um:
                    i, j, sqrp = um.groups()

                    #                     print i, j, sqrp

                    # keys are pair of indexes as smaller:larger
                    key = ":".join(sorted([i, j], reverse=True))
                    assert (key not in bp_prob_dict)
                    bpprob = float(sqrp)*float(sqrp)
                    bp_prob_dict[key] = bpprob

                    i, j = int(i), int(j)
                    bp_prob_mat[i-1, j-1] = bpprob
    return bp_prob_mat


def get_expected_accuracy(reference_struct, dp_matrix):
    '''dp_matrix is a numpy matrix where base indeices are ZERO based'''
    assert dp_matrix.shape[0] == dp_matrix.shape[1]
    assert dp_matrix.shape[0] == len(reference_struct)
    reference_struct_dict = dotbracket_to_dict(reference_struct)
    sum_TP_prob = 0.0
    for bp_key in reference_struct_dict:
        i, j = bp_key.split(":")
        i, j = int(i), int(j)
        sum_TP_prob += dp_matrix[i-1, j-1]
#         print i,j, dp_matrix[i-1,j-1]

#     print "    TP_score: %.2f" % (sum_TP_prob/len(reference_struct_dict))
    if len(reference_struct_dict) == 0:
        return 1
    return (sum_TP_prob/len(reference_struct_dict))

def get_structure_accuracy_from_ps_file(dp_ps, ref_struct):
    assert(os.path.isfile(dp_ps))
    dp_matrix = parse_dp_ps(dp_ps)
    print 'matrix dp', sum(dp_matrix)
    return get_expected_accuracy(ref_struct, dp_matrix)

#get_structure_accuracy(['./notebooks/dot.ps'], ['(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'])

def get_structure_accuracy(sequence,ref_struct, react=None):
    dpfile,_ = rnaplfold.call_vienna_plfold(sequence, react=react)
    return get_structure_accuracy_from_ps_file(dpfile, ref_struct=ref_struct)