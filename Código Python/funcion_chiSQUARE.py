#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:39:49 2020

@author:
    Chen Z, Zhao P, Li F, Leier A, Marquez-Lago TT, Wang Y, Webb GI, Smith AI, Daly RJ, Chou KC, Song J. 
    iFeature: a Python package and web server for features extraction and selection from protein and peptide sequences. 
    Bioinformatics. 2018 Jul 15;34(14):2499-2502. doi: 10.1093/bioinformatics/bty140. PMID: 29528364; PMCID: PMC6658705. 

Enlace online: https://github.com/Superzchen/iFeature

"""

#!/usr/bin/env python
#_*_coding:utf-8_*_

import numpy as np
import pandas as pd

binBox = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

def CHI2(encodings, labelFile):
    features = list(encodings.columns[1:])
    encodings = np.array(encodings)
    data = encodings[:, 1:]
    shape = data.shape
    data = np.reshape(data, shape[0] * shape[1])
    data = np.reshape([float(i) for i in data], shape)

    e = ''
    # if shape[0] < 5 or shape[1] < 2:
    #     return 0, e

    with open(labelFile) as f:
        records = f.readlines()
    myDict = {}
    try:
        for i in records:
            array = i.rstrip().split() if i.strip() != '' else None
            myDict[array[0]] = int(array[1])
    except IndexError as e:
        return 0, e

    labels = []
    for i in encodings:
        labels.append(myDict.get(i[0], 0))

    dataShape = data.shape

    if dataShape[1] != len(features):
        print('Error: inconsistent data shape with feature number.')
        return 0, 'Error: inconsistent data shape with feature number.'
    if dataShape[0] != len(labels):
        print('Error: inconsistent data shape with sample number.')
        return 0, 'Error: inconsistent data shape with sample number.'

    sampleNumber = len(data)
    labelClass = set(labels)
    myFea = {}
    for i in range(len(features)):
        array = data[:, i]
        newArray = list(pd.cut(array, len(binBox), labels= binBox))
        binBoxClass = set(newArray)
        myObservation = {}
        for  j in range(len(labels)):
		#print(labels[j], newArray[j])
        	myObservation[str(labels[j]) + str(newArray[j])] = myObservation.get(str(labels[j]) + str(newArray[j]), 0) + 1

        myExpect = {}
        for j in labelClass:
            for k in binBox:
                myExpect[str(j) + str(k)] = labels.count(j) * newArray.count(k) / sampleNumber

        chiValue = 0
        for j in labelClass:
            for k in binBoxClass:
                chiValue= chiValue + pow(((myObservation.get(str(j) + str(k), 0)) - myExpect.get(str(j) + str(k), 0)), 2) / myExpect[str(j) + str(k)]
                myFea[features[i]] = chiValue

    res = []
    res.append(['feature', 'CHI-value'])
    for key in sorted(myFea.items(), key=lambda item:item[1], reverse=True):
        res.append([key[0], '{0:.3f}'.format(myFea[key[0]])])
    return res





encodings_Allpeptides = pd.read_csv(r'BBDD/encodings_Allpeptides.tsv', sep = '\t')

chi2_allpeptides=CHI2(encodings=encodings_Allpeptides, labelFile='BBDD/ID_LABEL.csv')

chi2_allpeptides=pd.DataFrame(chi2_allpeptides)

chi2_allpeptides.to_csv("BBDD/chi2_allpeptides.csv", 
                sep='\t', index=False, header=False)





    
    
    
    
    
    
    

