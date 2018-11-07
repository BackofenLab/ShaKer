import subprocess

import dill

loadfile = lambda filename: dill.load(open(filename, "r"))
dumpfile = lambda thing, filename: dill.dump(thing, open(filename, "w"))


def shexec(cmd):
    '''
    takes cmd, the command
    returns (exit-code, stderr, stdout)
    '''
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, stderr = process.communicate()
    retcode = process.poll()
    return (retcode, stderr, output)