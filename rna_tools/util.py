import subprocess
from collections import defaultdict

import networkx as nx

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

def sequence_dotbracket_to_graph(seq_info=None, seq_struct=None):
    """Given a sequence and the dotbracket sequence make a graph.
    Parameters
    ----------
    seq_info string
        node labels eg a sequence string
    seq_struct  string
        dotbracket string
    Returns
    -------
        returns a nx.Graph
        secondary struct associated with seq_struct
    """
    graph = nx.Graph()

    lifo = defaultdict(list)
    open_brace_string={")":"(",
                "]":"[",
                ">":"<"}

    for i, (c, b) in enumerate(zip(seq_info, seq_struct)):
        graph.add_node(i, label=c, position=i)
        if i > 0:
            graph.add_edge(i, i - 1, label='-', type='backbone', len=1)
        if b in ['(','[','<']:
            lifo[b].append(i)
        if b in [')',']','>']:
            j = lifo[open_brace_string[b]].pop()
            graph.add_edge(i, j, label='=', type='basepair', len=1)


    return graph
