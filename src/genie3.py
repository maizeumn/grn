#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re
import _pickle

dir_genie3 = '/home/springer/zhoux379/source/git/GENIE3/GENIE3_python'
sys.path.insert(0, dir_genie3)
from GENIE3 import *
#dir_dyngenie3 = '/home/springer/zhoux379/source/git/dynGENIE3/dynGENIE3_python'
#sys.path.insert(0, dir_dyngenie3)
#from dynGENIE3 import *

def run_dynGENIE3(fi, fo, thread, tree_method, K, ntrees):
    fhi = open(fi, 'rb')
    (TS_data, time_points, decay_rates, gene_names) = _pickle.load(fhi)
    fhi.close()
    (VIM, alphas, prediction_score, stability_score, treeEstimators) = \
            dynGENIE3(TS_data, time_points, 
                    gene_names = gene_names, 
                    regulators = regulators, 
                    tree_method = tree_method, 
                    K = K, 
                    ntrees = ntrees, 
                    nthreads = thread)
    with open(fo, "wb") as fho:
        _pickle.dump(VIM, fho)
 
def run_GENIE3(fi, fo, thread, tree_method, K, ntrees):
    fhi = open(fi, 'rb')
    (exp_mat, tids, rids) = _pickle.load(fhi)
    fhi.close()
    VIM = GENIE3(exp_mat, gene_names = tids, regulators = rids,
                 tree_method = tree_method, 
                 K = K, 
                 ntrees = ntrees, 
                 nthreads = thread)
    with open(fo, "wb") as fho:
        _pickle.dump(VIM, fho, protocol = 4)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'run GENIE3 or dynGENIE3'
    )
    parser.add_argument(
            'fi', help = 'input file w. expression matrix and TF IDs'
    )
    parser.add_argument(
            'fo', help = 'output file'
    )
    parser.add_argument(
            '-p', '--thread', type = int, default = 1, help = 'threads'
    )
    parser.add_argument(
            '--tree_method', default = 'RF', help = 'tree method'
    )
    parser.add_argument(
            '--K', default = 'sqrt', help = 'K'
    )
    parser.add_argument(
            '--ntrees', type = int, default = 1000, help = 'num. trees to grow'
    )
    args = parser.parse_args()

    fi, fo = args.fi, args.fo
    thread = args.thread
    tree_method, K, ntrees = args.tree_method, args.K, args.ntrees
    run_GENIE3(fi, fo, thread, tree_method, K, ntrees)
