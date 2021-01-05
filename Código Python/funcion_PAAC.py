#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:41:00 2020

@author:
    Chen Z, Zhao P, Li F, Leier A, Marquez-Lago TT, Wang Y, Webb GI, Smith AI, Daly RJ, Chou KC, Song J. 
    iFeature: a Python package and web server for features extraction and selection from protein and peptide sequences. 
    Bioinformatics. 2018 Jul 15;34(14):2499-2502. doi: 10.1093/bioinformatics/bty140. PMID: 29528364; PMCID: PMC6658705. 

Enlace online: https://github.com/Superzchen/iFeature

"""


# Hemos utilizado una función proporcionada por iFeature para obtener los descriptores PAAC
# a partir de las secuencias peptídicas en formato fasta.

from codes import *
import re, sys, os, platform
import math
pPath = os.path.dirname(os.path.abspath("__file__"))
import os
# pPath = os.getcwd()
# os.path.join(os.path.dirname(__file__))
# pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta
import readFasta
import saveCode

USAGE = """
USAGE:
	python PAAC.py input.fasta <lambda> <output>

	input.fasta:      the input protein sequence file in fasta format.
	lambda:           the lambda value, integer, defaule: 30
	output:           the encoding file, default: 'encodings.tsv'
"""
# lambdaValue = 30
# fastas = readFasta.readFasta('AyB.txt')
# w=0.05

def Rvalue(aa1, aa2, AADict, Matrix):
	return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

def PAAC(fastas, lambdaValue=30, w=0.05, **kw):
# 	if checkFasta.minSequenceLengthWithNormalAA(fastas) < lambdaValue + 1:
# 		print('Error: all the sequence length should be larger than the lambdaValue+1: ' + str(lambdaValue + 1) + '\n\n')
# 		return 0


	dataFile = 'data/PAAC.txt' 
	with open(dataFile) as f:
		records = f.readlines()
	AA = ''.join(records[0].rstrip().split()[1:])
	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i
	AAProperty = []
	AAPropertyNames = []
	for i in range(1, len(records)):
		array = records[i].rstrip().split() if records[i].rstrip() != '' else None
		AAProperty.append([float(j) for j in array[1:]])
		AAPropertyNames.append(array[0])

	AAProperty1 = []
	for i in AAProperty:
		meanI = sum(i) / 20
		fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
		AAProperty1.append([(j-meanI)/fenmu for j in i])

	encodings = []
	header = ['#']
	for aa in AA:
		header.append('Xc1.' + aa)
	for n in range(1, lambdaValue + 1):
		header.append('Xc2.lambda' + str(n))
	encodings.append(header)
    
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		theta = []
		for n in range(1, lambdaValue + 1):
			theta.append(
				sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (
				len(sequence) - n) if (len(sequence) - n) != 0 else 0)
		myDict = {}
		for aa in AA:
			myDict[aa] = sequence.count(aa)
		code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
		code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
		encodings.append(code)
	return encodings



fastas = readFasta.readFasta('datasets_fasta/dataset_balanceado.txt')
 
encodings=PAAC(fastas)

saveCode.savetsv(encodings, 'BBDD/encodings_dataset_balanceado.tsv')



fastas = readFasta.readFasta('datasets_fasta/All_peptides.txt')
 
encodings=PAAC(fastas)

saveCode.savetsv(encodings, 'BBDD/encodings_Allpeptides.tsv')



fastas = readFasta.readFasta('datasets_fasta/dataset_ejemplo.csv')
 
encodings=PAAC(fastas)

saveCode.savetsv(encodings, 'BBDD/encodings_dataset_ejemplo.tsv')






