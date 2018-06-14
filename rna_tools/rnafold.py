from rna_tools import rna_io, shexec


def fold(sequence,react = None, return_dotplot= False):
    """call rna fold, return dotbracket string"""
    if react==None:
        cmd = 'echo "%s" | RNAfold '  % sequence
    else:
        with open('shap.tmp','w') as f:
            f.write(rna_io.format_shape("", react, noheader=True))
        cmd = 'echo "%s" | RNAfold --shape shap.tmp '  % sequence
    if return_dotplot:
        cmd+= " -p"
    res = shexec(cmd)[2]
    return res[res.find("\n")+1: res.find(" ")]