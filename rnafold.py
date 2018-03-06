import simushape_nostruct as sn
import simushape as ss
def fold(sequence,react):
    with open('shap.tmp','w') as f:
        f.write(ss.format_shape("",react,noheader=True))
    cmd = 'echo "%s" | RNAfold --shape shap.tmp'  % sequence
    res = sn.shexec(cmd)[2]
    return res[res.find("\n")+1: res.find(" ")]



