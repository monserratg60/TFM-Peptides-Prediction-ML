# -*- coding: utf-8 -*-
"""
Spyder Editor

author: Mariano Monserrat Gómez
"""


# Script para tratar las bases de datos analizadas y agrupar todos los péptidos 
# tóxicos y no tóxicos en un dataset.

# Creamos un dataset llamado 'All_peptides.txt' que contendra un identificador (ID), 
# la secuencia peptídica y la etiqueta (Label) del péptido tóxicos o no tóxico correspondiente.

# También, realizamos el balanceo de 'All_peptides.txt' y creamos el fichero 'dataset_balanceado.txt'


# Por último, creamos un fichero 'ID_LABEL.CSV' QUE VA A CONTENER EL ID DE CADA PÉPTIDO
# Y SU ETIQUETA.
#  LA ETIQUETA LA TRANSFORMAREMOS EN UNA VARIABLE CUANTITATIVA, SIENDO:  
#    Toxic = 1
#    NonToxic = 0

# Este fichero lo useramos para calcular los mejores descriptores según CHI2.



# Creación del fichero 'All_peptides.txt'

#### TOXINPRED.



# Leemos los péptidos NO TÓXICOS de TOXINPRED.

with open("BBDD/ToxinPred/toxinPred_NoToxic.txt", "r") as archivo:
    NoToxic_ToxinPred_1 = list(map(str.rstrip, archivo))
import pandas as pd
NoToxic_ToxinPred_1 = pd.DataFrame(NoToxic_ToxinPred_1)
NoToxic_ToxinPred_1['Label'] = 'NonToxic'

with open("BBDD/ToxinPred/toxinPred_NoToxic2.txt", "r") as archivo:
    NoToxic_ToxinPred_2 = list(map(str.rstrip, archivo))
NoToxic_ToxinPred_2 = pd.DataFrame(NoToxic_ToxinPred_2)
NoToxic_ToxinPred_2['Label'] = 'NonToxic'


# Leemos los péptidos TÓXICOS de ToxinPred.

with open("BBDD/ToxinPred/toxinPred_Toxic.txt", "r") as archivo:
    Toxic_ToxinPred_1 = list(map(str.rstrip, archivo))
import pandas as pd
Toxic_ToxinPred_1 = pd.DataFrame(Toxic_ToxinPred_1)
Toxic_ToxinPred_1['Label'] = 'Toxic'

with open("BBDD/ToxinPred/toxinPred_toxic2.txt", "r") as archivo:
    Toxic_ToxinPred_2 = list(map(str.rstrip, archivo))
Toxic_ToxinPred_2 = pd.DataFrame(Toxic_ToxinPred_2)
Toxic_ToxinPred_2['Label'] = 'Toxic'


# Agrupamos y eliminamos duplicados de péptidos no tóxicos:
NoToxic_ToxinPred=pd.concat([NoToxic_ToxinPred_1, NoToxic_ToxinPred_2])
NoToxic_ToxinPred = NoToxic_ToxinPred.drop_duplicates()


# Agrupamos y eliminamos duplicados de péptidos tóxicos:
Toxic_ToxinPred=pd.concat([Toxic_ToxinPred_1, Toxic_ToxinPred_2])
Toxic_ToxinPred = Toxic_ToxinPred.drop_duplicates()


# Renombramos los nombres de las columnas:
Toxic_ToxinPred = Toxic_ToxinPred.rename(columns={0:'seq'})
NoToxic_ToxinPred = NoToxic_ToxinPred.rename(columns={0:'seq'})






#### ARACHNO SERVER.


# Creamos el dataset con los péptidos tóxicos de ARACHNOSERVER.

import pandas as pd
ToxicPeptides_ArachnoServer = pd.read_csv("BBDD/ArachnoServer/ToxicPeptides_ArachnoServer.csv")
ToxicPeptides_ArachnoServer=ToxicPeptides_ArachnoServer.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_ArachnoServer['Label'] = 'Toxic'





#### SWISS PROT 

# Creamos el dataset de péptidos tóxicos de SWISS PROT
ToxicPeptides_SwissProt = pd.read_csv("BBDD/Swiss Prot/PeptideToxic_SwissProt.csv")
ToxicPeptides_SwissProt=ToxicPeptides_SwissProt.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_SwissProt['Label'] = 'Toxic'


# Creamos el dataset de péptidos NO tóxicos de SWISS PROT
Non_ToxicPeptides_SwissProt = pd.read_csv("BBDD/Swiss Prot/Peptide_Non_Toxic_SwissProt.csv")
Non_ToxicPeptides_SwissProt=Non_ToxicPeptides_SwissProt.drop(['Unnamed: 0'], axis=1)
Non_ToxicPeptides_SwissProt['Label'] = 'NonToxic'




#### CELIAC DATABASE


# Creamos el dataset de péptidos tóxicos:
ToxicPeptides_CeliacDataBase = pd.read_csv("BBDD/celiacDataBase/Celiac_Database_ToxinPeptides.csv")

# Eliminamos las columnas que no nos serán necesarias y añadimos la Label.
ToxicPeptides_CeliacDataBase=ToxicPeptides_CeliacDataBase.drop(['ID', 'Type', 'Description', 'Toxicity', 'Form', 'HLADQ', 'Refs', 'SeqLen'], axis=1)
ToxicPeptides_CeliacDataBase['Label'] = 'Toxic'
ToxicPeptides_CeliacDataBase = ToxicPeptides_CeliacDataBase.rename(columns={'Sequence':'seq'})






#### NTXpred


# Creamos el dataset de péptidos tóxicos de NTXpred
ToxicPeptides_NTXpred = pd.read_csv("BBDD/NTXpred/Peptide_neuroToxin_NTXpred.csv")
ToxicPeptides_NTXpred=ToxicPeptides_NTXpred.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_NTXpred['Label'] = 'Toxic'


# Creamos el dataset de péptidos NO tóxicos de NTXpred
Non_ToxicPeptides_NTXpred = pd.read_csv("BBDD/NTXpred/Peptide_Non_Toxic_NTXpred.csv")
Non_ToxicPeptides_NTXpred=Non_ToxicPeptides_NTXpred.drop(['Unnamed: 0'], axis=1)
Non_ToxicPeptides_NTXpred['Label'] = 'NonToxic'





#### BTXpred


# Creamos el dataset de péptidos tóxicos de BTXpred
ToxicPeptides_BTXpred = pd.read_csv("BBDD/BTXpred/Peptide_Toxic_BTXpred.csv")
ToxicPeptides_BTXpred=ToxicPeptides_BTXpred.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_BTXpred['Label'] = 'Toxic'



# Creamos el dataset de péptidos NO tóxicos de BTXpred
Non_ToxicPeptides_BTXpred = pd.read_csv("BBDD/BTXpred/Peptide_Non_Toxic_BTXpred.csv")
Non_ToxicPeptides_BTXpred=Non_ToxicPeptides_BTXpred.drop(['Unnamed: 0'], axis=1)
Non_ToxicPeptides_BTXpred['Label'] = 'NonToxic'





#### Kalium


# Creamos el dataset de péptidos tóxicos de Kalium
ToxicPeptides_Kalium = pd.read_csv("BBDD/Kalium/Kalium_ToxinPeptides.csv")
ToxicPeptides_Kalium=ToxicPeptides_Kalium.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_Kalium['Label'] = 'Toxic'
ToxicPeptides_Kalium = ToxicPeptides_Kalium.rename(columns={'x':'seq'})






#### T3DB


# Creamos el dataset de péptidos tóxicos de T3DB
ToxicPeptides_T3DB = pd.read_csv("BBDD/T3DB/Peptide_Toxic_T3DB.csv")
ToxicPeptides_T3DB=ToxicPeptides_T3DB.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_T3DB['Label'] = 'Toxic'






#### DBETH


# Creamos el dataset de péptidos tóxicos de T3DB
ToxicPeptides_DBETH = pd.read_csv("BBDD/DBETH/Peptide_Toxic_DBETH.csv")
ToxicPeptides_DBETH=ToxicPeptides_DBETH.drop(['Unnamed: 0'], axis=1)
ToxicPeptides_DBETH['Label'] = 'Toxic'





# A continuación, crearemos un dataset con los péptidos tóxicos y 
# otro que contenga los péptidos no tóxicos.   
# Eliminaremos los duplicados y los agruparemos en un conjunto de datos final.



# Péptidos tóxicos.


peptidos_toxicos=Toxic_ToxinPred.append(ToxicPeptides_ArachnoServer)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_SwissProt)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_CeliacDataBase)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_NTXpred)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_BTXpred)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_Kalium)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_T3DB)
peptidos_toxicos=peptidos_toxicos.append(ToxicPeptides_DBETH)


peptidos_toxicos = peptidos_toxicos.drop_duplicates()





# Péptidos NO tóxicos.


peptidos_NO_toxicos=NoToxic_ToxinPred.append(Non_ToxicPeptides_SwissProt)
peptidos_NO_toxicos=peptidos_NO_toxicos.append(Non_ToxicPeptides_NTXpred)
peptidos_NO_toxicos=peptidos_NO_toxicos.append(Non_ToxicPeptides_BTXpred)


peptidos_NO_toxicos = peptidos_NO_toxicos.drop_duplicates()



# Agrupamos todos los péptidos en un set ('All_peptides')

# Vamos a crear un índice en una columna. Ya que luego lo usaremos 
# para crear un ID que concatene la letra 'P' con el índice para cada péptido.

All_peptides = peptidos_toxicos.append(peptidos_NO_toxicos)

indice = range(0,len(All_peptides))
indice=list(indice)

All_peptides['indice']=indice

# Reconfiguramos el índice para que vaya de 0 a 480656.
All_peptides = All_peptides.reset_index(drop=True)

# Veo si hay duplicados. Si los hay, estaría mal el dataset.
All_peptides_v2 = All_peptides.drop_duplicates()
# Presenta el mismo número de péptidos.


# Exportamos un fichero .txt con todos los péptidos y lo ordenamos.
All_peptides.to_csv("BBDD/All_peptides.txt",index=None)


f = open('BBDD/All_peptides.txt', 'w')

for i in range(len(All_peptides)):

    f.write('P'+str(All_peptides.iloc[i].indice)+','+All_peptides.iloc[i].seq+','+All_peptides.iloc[i].Label+'\n')
        
f.close()




All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',', header=None)


All_peptides = All_peptides.rename(columns={0:'ID'})
All_peptides = All_peptides.rename(columns={1:'Sequence'})
All_peptides = All_peptides.rename(columns={2:'Label'})


# Exportamos un fichero .txt con todos los péptidos.
# Será el 'Dataset inicial' con el que partimos para realizar modelos predictivos y
# el que tratamos mediante  técnicas de balanceo y clustering.
All_peptides.to_csv("BBDD/All_peptides.txt",index=None)


# Exportamos el fichero 'All_peptides.txt' a la carpeta 'datasets_fasta' del proyecto
# para utilizar el script 'crearFasta.py' y convertilo en un fichero fasta.
All_peptides.to_csv("datasets_fasta/All_peptides.txt",index=None)







# Creación del fichero 'dataset_balanceado.txt'




## Balanceamos el dataset 'All_peptides' siguiendo la estrategia Subsampling
# Por lo que vamos reduciendo la clase mayoritaria hasta alcanzar la minoritaria.

# Leemos el fichero.
import pandas as pd
All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',')


# Creamos el índice 'ID_Sequence' que luego nos servirá para mergear el dataset 'All_peptides'
# con el dataset balanceado que creemos.
All_peptides['ID_Sequence']=range(len(All_peptides))

A=All_peptides.loc[:, ['Label']]
B=All_peptides.loc[:, ['ID_Sequence']]


# Categorizamos la variable cualitativa 'Label' a numérica
# Lable: Toxic = 1
#        NonToxic = 0
A = pd.get_dummies(A, columns = ["Label"], drop_first = True)

# Realizamos el balanceo:
from imblearn.under_sampling import NearMiss
nr = NearMiss()
B_res, A_res = nr.fit_sample(B, A)

dataset_balanceado=pd.merge(B_res, A_res, right_index=True, left_index=True)

dataset_balanceado=pd.merge(dataset_balanceado, All_peptides, on='ID_Sequence')

dataset_balanceado=dataset_balanceado.drop(columns=['ID_Sequence','Label_Toxic'])


# Vemos que el balanceo se ha realizado correctamente, ya que tenemos
# 11364 péptidos de cada clase.
dataset_balanceado.groupby('Label').size()


# Exportamos el fichero 'dataset_balanceado.txt' a la carpeta BBDD del proyecto.
dataset_balanceado.to_csv("BBDD/dataset_balanceado.txt",index=None)


# Exportamos el fichero 'dataset_balanceado.txt' a la carpeta 'datasets_fasta' del proyecto
# para al utilizar el script 'crearFasta.py' convertilo en un fichero fasta.
dataset_balanceado.to_csv("datasets_fasta/dataset_balanceado.txt",index=None)











# CREACIÓN DE UN FICHERO 'ID_LABEL.CSV' QUE VA A CONTENER EL ID DE CADA PÉPTIDO Y SU ETIQUETA.
#  LA ETIQUETA LA TRANSFORMAREMOS EN UNA VARIABLE CUANTITATIVA, SIENDO:  
#    Toxic = 1
#    NonToxic = 0

# Este fichero lo useramos para calcular los mejores descriptores según CHI2.


import pandas as pd

# Leemos los péptidos
All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',')


# Creamos el dataframe 'ID_LABEL' con el ID y Label de cada péptido.
ID_LABEL=pd.DataFrame()
ID_LABEL['ID']=All_peptides.iloc[0:len(All_peptides)].ID
ID_LABEL['label']=All_peptides.iloc[0:len(All_peptides)].Label
# Transformamos la variable cualitativa 'Label' en cuantitativa.
ID_LABEL = pd.get_dummies(ID_LABEL, columns = ["label"], drop_first = True)

# Exportamos el archivo que lo usaremos en la función chi-cuadrado proporcionada por iFeature.
ID_LABEL.to_csv("BBDD/ID_LABEL.csv", 
                sep='\t', index=False, header=False)












