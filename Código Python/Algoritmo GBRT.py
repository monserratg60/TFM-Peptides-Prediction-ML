#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 11:39:57 2021

@author: marianomonserratgomez
"""

# ALGORITMO Gradient Boasted Tree (GBRT).



import pandas as pd

# Leemos los péptidos
All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',')

# Leemos los descriptores PAAC de cada péptido.
PAAC_descriptors = pd.read_csv(r'BBDD/encodings_Allpeptides.tsv', sep = '\t')
PAAC_descriptors.rename(columns={'#': 'ID'}, inplace=True)



# Unimos los dos datasets.
All_peptides=pd.merge(PAAC_descriptors, All_peptides, on='ID')






# MODELO PREDICTIVO SVM CON EL DATASET INICIAL.


# Preparación del dataset para la generación de los modelos predictivos.

prep_modelo=All_peptides.drop(columns=['Sequence'])


# Barajamos los datos.
from sklearn.utils import shuffle
prep_modelo = shuffle(prep_modelo, random_state = 0)


# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X=prep_modelo.drop(columns=['ID','Label'])
Y = prep_modelo.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)



# Entrenamiento del modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [10],
              'learning_rate' : [0.1]}


clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 

y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix


mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_datasetInicial=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]

accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_ds_inicial = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_ds_inicial.loc[len(resGenerales_GBRT_ds_inicial)]=['Modelo dataset inicial', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_ds_inicial = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_ds_inicial.loc[len(matrizConfusion_GBRT_ds_inicial)]=['Modelo dataset inicial', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_ds_inicial.loc[len(matrizConfusion_GBRT_ds_inicial)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_ds_inicial.loc[len(matrizConfusion_GBRT_ds_inicial)]=['Toxic', FN, TN]


















# MODELO PREDICTIVO SVM para el mejor descriptor.


# Dataset con la columna Xc2.lambda30

dataset_Xc2_lambda30 = All_peptides.loc[:, ['ID', 'Xc2.lambda30', 'Label']]


# Lo primero de todo va a ser barajar los datos.
from sklearn.utils import shuffle
dataset_Xc2_lambda30 = shuffle(dataset_Xc2_lambda30, random_state = 0)


X=dataset_Xc2_lambda30.drop(columns=['ID','Label'])
Y = dataset_Xc2_lambda30.iloc[:,-1]



# import numpy as np
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split, GridSearchCV

x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


# Entrenamiento del modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [50],
              'learning_rate' : [0.1]}


clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 
y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix


mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_1descriptor=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]
accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_1descrip = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_1descrip.loc[len(resGenerales_GBRT_1descrip)]=['Modelo mejor descriptor', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_1descrip = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_1descrip.loc[len(matrizConfusion_GBRT_1descrip)]=['Modelo mejor descriptor', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_1descrip.loc[len(matrizConfusion_GBRT_1descrip)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_1descrip.loc[len(matrizConfusion_GBRT_1descrip)]=['Toxic', FN, TN]

















# MODELO SVM PARA DATASET CON LOS DOS MEJORES DESCRIPTORES.


# Dataset con las columnas Xc2.lambda30 y Xc2.lambda29.

dataset_Xc2_lambda29y30 = All_peptides.loc[:, ['ID', 'Xc2.lambda30', 'Xc2.lambda29', 'Label']]


# Lo primero de todo va a ser barajar los datos.
from sklearn.utils import shuffle
dataset_Xc2_lambda29y30 = shuffle(dataset_Xc2_lambda29y30, random_state = 0)



X=dataset_Xc2_lambda29y30.drop(columns=['ID','Label'])
Y = dataset_Xc2_lambda29y30.iloc[:,-1]




# import numpy as np
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split, GridSearchCV

x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


# Entrenamiento del modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [50],
              'learning_rate' : [0.1]}

clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 

y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix


mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_2descriptores=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]

accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_2descrip = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_2descrip.loc[len(resGenerales_GBRT_2descrip)]=['Modelo 2 mejores descriptores', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_2descrip = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_2descrip.loc[len(matrizConfusion_GBRT_2descrip)]=['Modelo 2 mejores descriptores', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_2descrip.loc[len(matrizConfusion_GBRT_2descrip)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_2descrip.loc[len(matrizConfusion_GBRT_2descrip)]=['Toxic', FN, TN]















# MODELO SVM PARA DATASET BALANCEADO (ESTRATEGIA SUBSAMPLING).



# MODELO RF PARA DATASET BALANCEADO (ESTRATEGIA SUBSAMPLING).

# Modelos generados a partir del dataset balanceado con subsampling de la clase mayoritaria.

# Leemos el dataset balanceado y el fichero con los descriptores:
import pandas as pd


dataset_balanceado = pd.read_csv(r'BBDD/dataset_balanceado.txt', sep = ',')


encodings_dataset_balanceado = pd.read_csv(r'BBDD/encodings_dataset_balanceado.tsv', sep = '\t')
encodings_dataset_balanceado.rename(columns={'#': 'ID'}, inplace=True)


# Unimos los dos datasets.
dataset_balanceado=pd.merge(encodings_dataset_balanceado, dataset_balanceado, on='ID')




# Preparación del dataset para la generación de los modelos predictivos.
prep_modelo=dataset_balanceado.drop(columns=['Sequence'])


# Barajamos los datos.
from sklearn.utils import shuffle
prep_modelo = shuffle(prep_modelo, random_state = 0)


# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X=prep_modelo.drop(columns=['ID','Label'])
Y = prep_modelo.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)



# Entrenamiento del modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [50, 100, 200, 500],
              'learning_rate' : [0.1]}


clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 

# MEJOR PARÁMETRO.
# clf.best_params_
# {'learning_rate': 0.1, 'n_estimators': 500}

y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix


mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_Subsampling=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]
accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_SUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_SUBSAMPLING.loc[len(resGenerales_GBRT_SUBSAMPLING)]=['Modelo Subsampling', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_SUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_SUBSAMPLING.loc[len(matrizConfusion_GBRT_SUBSAMPLING)]=['Modelo Subsampling', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_SUBSAMPLING.loc[len(matrizConfusion_GBRT_SUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_SUBSAMPLING.loc[len(matrizConfusion_GBRT_SUBSAMPLING)]=['Toxic', FN, TN]













# MODELO GBRT PARA DATASET CLUSTERING DBSCAN.



# Utilizamos la función 'DBSCAN.py' que tendrá como argumento de entrada
# el fichero 'encodings_Allpeptides.tsv'. Leemos el fichero 'DBSCAN_Allpeptides' resultante
# con el clustering ya realizado.
DBSCAN_ifeature = pd.read_csv(r'BBDD/DBSCAN_Allpeptides.csv', 
                              sep="\t", header=None)

# Renombramos columnas
DBSCAN_ifeature.rename(columns={0: 'ID'}, inplace=True)
DBSCAN_ifeature.rename(columns={1: 'cluster'}, inplace=True)

# Unimos los ficheros
All_peptides=pd.merge(DBSCAN_ifeature, All_peptides, on='ID')


# Vamos a ver el número de outliers.
All_peptides.groupby('cluster').size()

# Eliminamos los 446726 outliers:
copy = All_peptides.drop(All_peptides[All_peptides['cluster'] == -1].index)



# MODELO GBRT del modelo DBSCAN.

# Preparación del dataset para la generación de los modelos predictivos.
# Vamos a quedarnos con las columnas que nos interesan del dataset creado.

# Modelo DBSCAN.
m_DBSCAN=copy.drop(columns=['Sequence', 'cluster'])


# Barajamos los datos.
from sklearn.utils import shuffle
m_DBSCAN = shuffle(m_DBSCAN, random_state = 0)

# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X=m_DBSCAN.drop(columns=['ID','Label'])
Y = m_DBSCAN.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)



# Entrenamos el modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [50,100,200,500],
              'learning_rate' : [0.1]}



clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 

# MEJOR PARÁMETRO.
# clf.best_params_
# {'learning_rate': 0.1, 'n_estimators': 50}

y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix


mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_DBSCAN=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]
accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_DBSCAN = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_DBSCAN.loc[len(resGenerales_GBRT_DBSCAN)]=['Modelo DBSCAN', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_DBSCAN = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_DBSCAN.loc[len(matrizConfusion_GBRT_DBSCAN)]=['Modelo DBSCAN', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_DBSCAN.loc[len(matrizConfusion_GBRT_DBSCAN)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_DBSCAN.loc[len(matrizConfusion_GBRT_DBSCAN)]=['Toxic', FN, TN]












# MODELO GBRT PARA EL DATASET DBSCAN + SUBSAMPLING.

DBSCAN_subsampling = pd.read_csv(r'BBDD/DBSCAN_subsampling.txt', sep = ',')



# Preparación del dataset para la generación de los modelos predictivos.
# Vamos a quedarnos con las columnas que nos interesan del dataset creado.

# Modelo DBSCAN con subsampling
m_DBSCAN_subsampling=DBSCAN_subsampling.drop(columns=['Sequence', 'ID_Sequence', 'Label_Toxic'
                                                      ,'cluster'])


# Barajamos los datos.
from sklearn.utils import shuffle
m_DBSCAN_subsampling = shuffle(m_DBSCAN_subsampling, random_state = 0)


# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X=m_DBSCAN_subsampling.drop(columns=['ID','Label'])
Y = m_DBSCAN_subsampling.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


# Entrenamos el modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [50, 100, 200, 500],
              'learning_rate' : [0.01, 0.1]}


clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 

# MEJOR PARÁMETRO.
# clf.best_params_
# {'learning_rate': 0.1, 'n_estimators': 200}

y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix

mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_DBSCANySubsampling=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]
accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_DBSCANySUBSAMPLING.loc[len(resGenerales_GBRT_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+Subsampling', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_DBSCANySUBSAMPLING.loc[len(matrizConfusion_GBRT_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+Subsampling', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_DBSCANySUBSAMPLING.loc[len(matrizConfusion_GBRT_DBSCANySUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_DBSCANySUBSAMPLING.loc[len(matrizConfusion_GBRT_DBSCANySUBSAMPLING)]=['Toxic', FN, TN]










# MODELO GBRT PARA DATASET CLUSTERING LINCLUST/MMseqs2.

# Leemos el fichero creado en el script del algoritmo SVM
peptidosLinclust = pd.read_csv(r'BBDD/peptidosLinclust.csv'
                                 , sep = '\t')


# Preparación del dataset para la generación de los modelos predictivos.
# Vamos a quedarnos con las columnas que nos interesan del dataset creado.


Modelo_linclust=peptidosLinclust.drop(columns=['Sequence'])


# Barajamos los datos.
from sklearn.utils import shuffle
Modelo_linclust = shuffle(Modelo_linclust, random_state = 0)


# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X = Modelo_linclust.drop(columns=['ID','Label'])
Y = Modelo_linclust.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


# Generamos el MODELO GBRT.


# Entrenamos el modelo.

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV


param_grid = {'n_estimators' : [50, 100, 200, 500],
              'learning_rate' : [0.01, 0.1]}


clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)

model=clf.fit(x_train,y_train) 

# MEJOR PARÁMETRO.
# clf.best_params_
# {'learning_rate': 0.1, 'n_estimators': 200}

y_pred_GBRT = model.predict(x_test)



from sklearn.metrics import classification_report, confusion_matrix

mat_Grad = confusion_matrix(y_test, y_pred_GBRT)
c_GBRT_Linclust=classification_report(y_test,y_pred_GBRT)
TP = mat_Grad[0][0]
TN = mat_Grad[1][1]
FP = mat_Grad[0][1]
FN = mat_Grad[1][0]
accuracy_Grad =(TP+TN)/(TP+TN+FN+FP)
precision_Grad = TP/(TP+FP)
sensitivity_Grad = TP/(TP+FN)
specifity_Grad = TN/(TN+FP)
F1_GBRT = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_Linclust = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_Linclust.loc[len(resGenerales_GBRT_Linclust)]=['Modelo Linclust', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_Linclust = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_Linclust.loc[len(matrizConfusion_GBRT_Linclust)]=['Modelo Linclust', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_Linclust.loc[len(matrizConfusion_GBRT_Linclust)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_Linclust.loc[len(matrizConfusion_GBRT_Linclust)]=['Toxic', FN, TN]











# Guardamos los resultados de los todos los modelos generados.

Resultados_GBRT_I = resGenerales_GBRT_ds_inicial.append(resGenerales_GBRT_1descrip)
Resultados_GBRT_I = Resultados_GBRT_I.append(resGenerales_GBRT_2descrip)


Resultados_GBRT_II = resGenerales_GBRT_SUBSAMPLING.append(resGenerales_GBRT_DBSCAN)
Resultados_GBRT_II = Resultados_GBRT_II.append(resGenerales_GBRT_DBSCANySUBSAMPLING)
Resultados_GBRT_II = Resultados_GBRT_II.append(resGenerales_GBRT_Linclust)





Resultados_GBRT= pd.ExcelWriter('RESULTADOS MODELOS/Resultados_GBRT_I.xlsx')
Resultados_GBRT_I.to_excel(Resultados_GBRT, sheet_name="Resultados_GBRT_I", index=False)

Resultados_GBRT.save()
Resultados_GBRT.close()

Resultados_GBRT= pd.ExcelWriter('RESULTADOS MODELOS/Resultados_GBRT_II.xlsx')
Resultados_GBRT_II.to_excel(Resultados_GBRT, sheet_name="Resultados_GBRT_II", index=False)


Resultados_GBRT.save()
Resultados_GBRT.close()


# Matriz de confusión.

matrizConfusion_GBRT_I = matrizConfusion_GBRT_ds_inicial.append(matrizConfusion_GBRT_1descrip)
matrizConfusion_GBRT_I = matrizConfusion_GBRT_I.append(matrizConfusion_GBRT_2descrip)


matrizConfusion_GBRT_II = matrizConfusion_GBRT_SUBSAMPLING.append(matrizConfusion_GBRT_DBSCAN)
matrizConfusion_GBRT_II = matrizConfusion_GBRT_II.append(matrizConfusion_GBRT_DBSCANySUBSAMPLING)
matrizConfusion_GBRT_II = matrizConfusion_GBRT_II.append(matrizConfusion_GBRT_Linclust)



matrizConfusion_GBRT= pd.ExcelWriter('RESULTADOS MODELOS/matrizConfusion_GBRT_I.xlsx')
matrizConfusion_GBRT_I.to_excel(matrizConfusion_GBRT, sheet_name="matrizConfusion_GBRT_I", index=False)

matrizConfusion_GBRT.save()
matrizConfusion_GBRT.close()

matrizConfusion_GBRT= pd.ExcelWriter('RESULTADOS MODELOS/matrizConfusion_GBRT_II.xlsx')
matrizConfusion_GBRT_II.to_excel(matrizConfusion_GBRT, sheet_name="matrizConfusion_GBRT_II", index=False)


matrizConfusion_GBRT.save()
matrizConfusion_GBRT.close()

