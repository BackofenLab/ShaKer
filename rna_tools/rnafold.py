import rna_io
from rna_tools import shexec


import tempfile


def fold(sequence,react = None, return_dotplot= False):
    """call rna fold, return dotbracket string"""
    if type(react)==type(None):
        cmd = 'echo "%s" | RNAfold --noPS '  % sequence
    else:
        tmpname = tempfile._get_default_tempdir()+'/'+next(tempfile._get_candidate_names()) + "shap.tmp"
        with open(tmpname,'w') as f:
            f.write(rna_io.format_shape("", react, noheader=True))
        cmd = 'echo "%s" | RNAfold --noPS --shape %s '  % (sequence,tmpname)
    if return_dotplot:
        cmd+= " -p"
    res = shexec(cmd)[2]
    return res[res.find("\n")+1: res.find(" ")]