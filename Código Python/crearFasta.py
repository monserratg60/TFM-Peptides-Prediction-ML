#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:35:00 2020

@author: marianomonserratgomez
"""

# Función para rear ficheros en formato fasta.



# Creamos el dataset balanceado en formato fasta:


# IMPORTANTE:

# Antes de crear el fichero fasta del datasets 'dataset_balanceado.txt '
# debemos ejecutar el script 'BBDD.py'


# Leemos el fichero.
import pandas as pd
dataset_balanceado = pd.read_csv(r'BBDD/dataset_balanceado.txt', 
                           sep = ',')

# En el script 'BBDD.py' también hemos exportado el fichero 'dataset_balanceado.txt'.
f = open('datasets_fasta/dataset_balanceado.txt', 'w')
    
for i in range(len(dataset_balanceado)):
    
    f.write('>'+str(dataset_balanceado.iloc[i].ID)+' '+dataset_balanceado.iloc[i].Label+'\n')
    f.write(dataset_balanceado.iloc[i].Sequence+'\n')

f.close()
   
    




# Creamos el dataset All_peptides en formato fasta:

 
# IMPORTANTE:
    
# Antes de crear el fichero fasta del datasets 'All_peptides.txt '
# debemos ejecutar el script 'BBDD.py'



# Leemos el fichero.
import pandas as pd
All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',')

# En el script 'BBDD.py' también hemos exportado el fichero 'All_peptides.txt'.
f = open('datasets_fasta/All_peptides.txt', 'w')
    
for i in range(len(All_peptides)):
    
    f.write('>'+str(All_peptides.iloc[i].ID)+' '+All_peptides.iloc[i].Label+'\n')
    f.write(All_peptides.iloc[i].Sequence+'\n')

f.close()




# Creamos el dataset externo en formato fasta:

    
# IMPORTANTE:
    
# Antes de crear el fichero fasta del datasets 'dataset_ejemplo.csv'
# debemos ejecutar el script 'dt_pepExternos.py'



# Leemos el fichero.
import pandas as pd
dataset_balanceado = pd.read_csv(r'PEPTIDOS EXTERNOS/dataset_ejemplo.csv', 
                           sep = '\t')

# Es necesario copiar el fichero 'dataset_ejemplo.csv' y pegarlo dentro de la 
# carpeta datasets_fasta antes de ejecutar las siguientes líneas y crear el fichero.


f = open('datasets_fasta/dataset_ejemplo.csv', 'w')
    
for i in range(len(dataset_balanceado)):
    
    f.write('>'+str(dataset_balanceado.iloc[i].ID)+' '+dataset_balanceado.iloc[i].Label+'\n')
    f.write(dataset_balanceado.iloc[i].Sequence+'\n')

f.close()











