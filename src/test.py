#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re
import _pickle
from sklearn.tree.tree import BaseDecisionTree
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from numpy import *
import time
from operator import itemgetter
from multiprocessing import Pool

diri = '/home/springer/zhoux379/projects/grn/data/11_input'
tid = 'Zm00001d040002'

f1 = op.join(diri, 'n14a.pkl')
(exp_mat, rids, tids, tids0, exp_mat0) = _pickle.load(open(f1,'rb'))

exp_mat.shape
len(rids)
len(tids)

ridxs = [i for i, gene in enumerate(tids) if gene in rids]
len(ridxs)
assert(len(ridxs)==len(rids))

tidx = [i for i, gene in enumerate(tids) if gene == tid][0]
tidx in ridxs
tids[tidx]

x_sample = exp_mat[:,ridxs]
x_label = exp_mat[:,tidx]


ntrees = 1000
K = 'sqrt'
max_features = K
res = RandomForestRegressor(n_estimators=ntrees,max_features=max_features)
res.fit(x_sample, x_label)
res.score(x_sample, x_label)


f2 = op.join(diri, 'n16a.pkl')
(exp_mat, rids2, tids2, tids0, exp_mat0) = _pickle.load(open(f2,'rb'))

ridxs = [i for i, gene in enumerate(tids0) if gene in rids]
len(ridxs)
assert(len(ridxs)==len(rids))

tidx = [i for i, gene in enumerate(tids0) if gene == tid][0]
tidx in ridxs
tids0[tidx]

y_sample = exp_mat0[:,ridxs]
y_label = exp_mat0[:,tidx]

res.score(y_sample, y_label)
