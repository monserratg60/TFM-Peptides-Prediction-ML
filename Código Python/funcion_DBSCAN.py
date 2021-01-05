#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:36:02 2020

@author:
    Chen Z, Zhao P, Li F, Leier A, Marquez-Lago TT, Wang Y, Webb GI, Smith AI, Daly RJ, Chou KC, Song J. 
    iFeature: a Python package and web server for features extraction and selection from protein and peptide sequences. 
    Bioinformatics. 2018 Jul 15;34(14):2499-2502. doi: 10.1093/bioinformatics/bty140. PMID: 29528364; PMCID: PMC6658705. 

Enlace online: https://github.com/Superzchen/iFeature

"""

#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

def dbscan(encodings, **kw):

    encodings4 = np.array(encodings)

    if len(encodings4) < 5:
        return 0, 'sample number should be greater than 5.'

    data2 = encodings4[:, 1:]
    shape2 = data2.shape
    data2 = np.reshape(data2, shape2[0] * shape2[1])
    data2 = np.reshape([float(i) for i in data2], shape2)

    # e = ''
    dataNew = StandardScaler().fit_transform(data2)
    db = DBSCAN().fit(dataNew)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    res = []
    for i in range(len(data2)):
        res.append([encodings['#'][i], labels[i]])
    return res



import pandas as pd

encodings_Allpeptides = pd.read_csv(r'BBDD/encodings_Allpeptides.tsv', sep = '\t')

res=dbscan(encodings=encodings_Allpeptides)

DBSCAN_Allpeptides=pd.DataFrame(res)

DBSCAN_Allpeptides.to_csv("BBDD/DBSCAN_Allpeptides.csv", 
                sep='\t', index=False, header=False)






