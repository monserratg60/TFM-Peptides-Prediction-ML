# -*- coding: utf-8 -*-
"""
Spyder Editor

author: Mariano Monserrat Gómez

"""




import pandas as pd

# Leemos los péptidos
All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',')



# Leemos los descriptores PAAC de cada péptido.
PAAC_descriptors = pd.read_csv(r'BBDD/encodings_Allpeptides.tsv', sep = '\t')
PAAC_descriptors.rename(columns={'#': 'ID'}, inplace=True)



# Unimos los dos datasets.
All_peptides=pd.merge(PAAC_descriptors, All_peptides, on='ID')






###### MODELO PREDICTIVO SVM CON EL DATASET INICIAL.


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


from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV


# Seleccionamos los parámetros del modelo.
parameters = {
    'kernel':['rbf'],
    'gamma':[1],
    'C':[1]
}


# Realizamos el entrenamiento del modelo con cross-validation (CV) = 10.
clf = GridSearchCV(SVC(), parameters, cv=10)
model_SVM=clf.fit(x_train,y_train)


# Una vez ajustado el modelo, vamos a realizar predicciones.
y_pred_SVM = model_SVM.predict(x_test)



# Parámetros para evaluar la calidad del modelo:
    
mat_SVM = confusion_matrix(y_test, y_pred_SVM)
c_SVM=classification_report(y_test,y_pred_SVM)

TP = mat_SVM[0][0]
TP
TN = mat_SVM[1][1]
TN
FP = mat_SVM[0][1]
FP
FN = mat_SVM[1][0]
FN

accuracy =(TP+TN)/(TP+TN+FN+FP)
accuracy
precision = TP/(TP+FP)
precision
sensitivity = TP/(TP+FN)
sensitivity
specifity = TN/(TN+FP)
specifity

F1 = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_SVM, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_ds_inicial = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_ds_inicial.loc[len(resGenerales_SVM_ds_inicial)]=['Modelo dataset inicial', 
                                        round(accuracy,2),
                                        round(precision,2), round(sensitivity,2), 
                                        round(specifity,2),
                                        round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_ds_inicial = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_ds_inicial.loc[len(matrizConfusion_SVM_ds_inicial)]=['Modelo dataset inicial', 'NonToxic', 'Toxic']
matrizConfusion_SVM_ds_inicial.loc[len(matrizConfusion_SVM_ds_inicial)]=['NonToxic', TP, FP]
matrizConfusion_SVM_ds_inicial.loc[len(matrizConfusion_SVM_ds_inicial)]=['Toxic', FN, TN]








###### MODELO PREDICTIVO SVM para el mejor descriptor.


# Dataset con la columna Xc2.lambda30

dataset_Xc2_lambda30 = All_peptides.loc[:, ['ID', 'Xc2.lambda30', 'Label']]


# Lo primero de todo va a ser barajar los datos.
from sklearn.utils import shuffle
dataset_Xc2_lambda30 = shuffle(dataset_Xc2_lambda30, random_state = 0)


X=dataset_Xc2_lambda30.drop(columns=['ID','Label'])
Y = dataset_Xc2_lambda30.iloc[:,-1]



# import numpy as np
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, GridSearchCV

x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


parameters = {
    'kernel':['rbf'],
    'gamma':[1],
    'C':[1]
}



clf = GridSearchCV(SVC(), parameters, cv=10)


model=clf.fit(x_train,y_train)



y_pred = model.predict(x_test)


mat = confusion_matrix(y_test, y_pred)
c_SVM_1descrip=classification_report(y_test,y_pred)
TP = mat[0][0]
TN = mat[1][1]
FP = mat[0][1]
FN = mat[1][0]
F1 = 2*TP/(2*TP+FP+FN)
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)



from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = False)
predictions_num = pd.get_dummies(y_pred, drop_first = False)

fpr, tpr, thresholds = roc_curve(y_test_num['NonToxic'], predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num['NonToxic'], predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_1descrip = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_1descrip.loc[len(resGenerales_SVM_1descrip)]=['Modelo mejor descriptor', 
                                           round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_1descrip = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_1descrip.loc[len(matrizConfusion_SVM_1descrip)]=['Modelo mejor descriptor', 'NonToxic', 'Toxic']
matrizConfusion_SVM_1descrip.loc[len(matrizConfusion_SVM_1descrip)]=['NonToxic', TP, FP]
matrizConfusion_SVM_1descrip.loc[len(matrizConfusion_SVM_1descrip)]=['Toxic', FN, TN]







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
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, GridSearchCV

x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


parameters = {
    'kernel':['rbf'],
    'gamma':[1],
    'C':[1]
}

# parameters = {
#     'kernel':['rbf'],
#     'gamma':[1e-3,1e-2,0.1,0.5,1],
#     'C':[0.1,0.5,1, 10, 100, 1000]
# }

clf = GridSearchCV(SVC(), parameters, cv=10)


model=clf.fit(x_train,y_train)



y_pred = model.predict(x_test)


mat = confusion_matrix(y_test, y_pred)
c_SVM_2descrip=classification_report(y_test,y_pred)
TP = mat[0][0]
TN = mat[1][1]
FP = mat[0][1]
FN = mat[1][0]
F1 = 2*TP/(2*TP+FP+FN)
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)





from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = False)
predictions_num = pd.get_dummies(y_pred, drop_first = False)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_2descrip = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_2descrip.loc[len(resGenerales_SVM_2descrip)]=['Modelo 2 mejores descriptores', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_2descrip = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_2descrip.loc[len(matrizConfusion_SVM_2descrip)]=['Modelo 2 mejores descriptores', 'NonToxic', 'Toxic']
matrizConfusion_SVM_2descrip.loc[len(matrizConfusion_SVM_2descrip)]=['NonToxic', TP, FP]
matrizConfusion_SVM_2descrip.loc[len(matrizConfusion_SVM_2descrip)]=['Toxic', FN, TN]
















# MODELO SVM PARA DATASET BALANCEADO (ESTRATEGIA SUBSAMPLING).



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




###### Generación del modelo SVM.


from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV


# Seleccionamos los parámetros del modelo.

parameters = {
    'kernel':['rbf'],
    'gamma':[1e-3,1e-2,0.1,0.5,1],
    'C':[0.1,0.5,1, 10, 100, 1000]
}


# Realizamos el entrenamiento del modelo con 10fold-cross-validation.
clf = GridSearchCV(SVC(), parameters, cv=10)
model_SVM=clf.fit(x_train,y_train)



# MEJORES PARÁMETROS.
# clf.best_params_
# {'C': 10, 'gamma': 0.01, 'kernel': 'rbf'}


# Una vez ajustado el modelo, vamos a realizar predicciones.
y_pred_SVM = model_SVM.predict(x_test)



# Parámetros para evaluar la calidad del modelo:
    
mat_SVM = confusion_matrix(y_test, y_pred_SVM)
c_SVM_subsampling=classification_report(y_test,y_pred_SVM)
TP = mat_SVM[0][0]
TN = mat_SVM[1][1]
FP = mat_SVM[0][1]
FN = mat_SVM[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)

from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_SVM, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_SUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_SUBSAMPLING.loc[len(resGenerales_SVM_SUBSAMPLING)]=['Modelo Subsampling', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_SUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_SUBSAMPLING.loc[len(matrizConfusion_SVM_SUBSAMPLING)]=['Modelo Subsampling', 'NonToxic', 'Toxic']
matrizConfusion_SVM_SUBSAMPLING.loc[len(matrizConfusion_SVM_SUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_SVM_SUBSAMPLING.loc[len(matrizConfusion_SVM_SUBSAMPLING)]=['Toxic', FN, TN]










# MODELO SVM PARA DATASET CLUSTERING DBSCAN.



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

# Y agrupamos por categoría para ver si nos queda un modelo desbalanceado.
# peptidos no tóxicos=33227 y tóxicos=703
copy.groupby('Label').size()


# Así que balanceamos siguiendo nuevamente la estrategia subsampling de la clase mayoritaria.

# Creamos un índice que luego nos servirá para mergear el dataset 'All_peptides'
# con el dataset balanceado que creemos.
copy['ID_Sequence']=range(len(copy))

A=copy.loc[:, ['Label']]
B=copy.loc[:, ['ID_Sequence']]


# Categorizamos la variable cualitativa 'Label' a numérica
# Lable: Toxic = 1
#        NonToxic = 0
A = pd.get_dummies(A, columns = ["Label"], drop_first = True)

# Realizamos el balanceo:
from imblearn.under_sampling import NearMiss
nr = NearMiss()
B_res, A_res = nr.fit_sample(B, A)

DBSCAN_subsampling=pd.merge(B_res, A_res, right_index=True, left_index=True)

DBSCAN_subsampling=pd.merge(DBSCAN_subsampling, copy, on='ID_Sequence')

# Vemos que el balanceo se ha realizado correctamente, ya que tenemos
# 703 péptidos de cada clase.
DBSCAN_subsampling.groupby('Label').size()



# Exportamos el fichero para poder utilizarlo en los demás algoritmos.
DBSCAN_subsampling.to_csv("BBDD/DBSCAN_subsampling.txt",index=None)




# MODELO SVM del modelo DBSCAN.

# Preparación del dataset para la generación de los modelos predictivos.
# Vamos a quedarnos con las columnas que nos interesan del dataset creado.

# Modelo DBSCAN.
m_DBSCAN=copy.drop(columns=['Sequence', 'cluster', 'ID_Sequence'])






# Barajamos los datos.
from sklearn.utils import shuffle
m_DBSCAN = shuffle(m_DBSCAN, random_state = 0)


# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X=m_DBSCAN.drop(columns=['ID','Label'])
Y = m_DBSCAN.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)








from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

parameters = {
    'kernel':['rbf'],
    'gamma':[1e-3,1e-2,0.1,0.5,1],
    'C':[0.1,0.5,1, 10, 100, 1000]
}




# Realizamos el entrenamiento del modelo con cross-validation (CV) = 10.
clf = GridSearchCV(SVC(), parameters, cv=10)
model=clf.fit(x_train,y_train)


# clf.best_params_
# BEST PARAMETERS PARA EL MODELO DBSCAN.
# {'C': 1, 'gamma': 0.1, 'kernel': 'rbf'}



y_pred = model.predict(x_test)
mat = confusion_matrix(y_test, y_pred)
c_SVM_DBSCAN=classification_report(y_test,y_pred)
TP = mat[0][0]
TN = mat[1][1]
FP = mat[0][1]
FN = mat[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)



from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_DBSCAN = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_DBSCAN.loc[len(resGenerales_SVM_DBSCAN)]=['Modelo DBSCAN', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_DBSCAN = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_DBSCAN.loc[len(matrizConfusion_SVM_DBSCAN)]=['Modelo DBSCAN', 'NonToxic', 'Toxic']
matrizConfusion_SVM_DBSCAN.loc[len(matrizConfusion_SVM_DBSCAN)]=['NonToxic', TP, FP]
matrizConfusion_SVM_DBSCAN.loc[len(matrizConfusion_SVM_DBSCAN)]=['Toxic', FN, TN]



















# MODELO SVM PARA EL DATASET DBSCAN + SUBSAMPLING.


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







# Entrenamiento del MODELO SVM del modelo DBSCAN + SUBSAMPLING.

from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

parameters = {
    'kernel':['rbf'],
    'gamma':[1e-3,1e-2,0.1,0.5,1],
    'C':[0.1,0.5,1, 10, 100, 1000]
}


# Realizamos el entrenamiento del modelo con cross-validation (CV) = 10.
clf = GridSearchCV(SVC(), parameters, cv=10)
model=clf.fit(x_train,y_train)

# clf.best_params_
# BEST PARAMETERS PARA EL MODELO DBSCAN_SUBSAMPLING.
# {'C': 100, 'gamma': 0.01, 'kernel': 'rbf'}


y_pred = model.predict(x_test)
mat = confusion_matrix(y_test, y_pred)
c_SVM_DBSCANySubsampling=classification_report(y_test,y_pred)
TP = mat[0][0]
TN = mat[1][1]
FP = mat[0][1]
FN = mat[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)



from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_DBSCANySUBSAMPLING.loc[len(resGenerales_SVM_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+SUBSAMPLING', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_DBSCANySUBSAMPLING.loc[len(matrizConfusion_SVM_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+SUBSAMPLING', 'NonToxic', 'Toxic']
matrizConfusion_SVM_DBSCANySUBSAMPLING.loc[len(matrizConfusion_SVM_DBSCANySUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_SVM_DBSCANySUBSAMPLING.loc[len(matrizConfusion_SVM_DBSCANySUBSAMPLING)]=['Toxic', FN, TN]









# MODELO SVM PARA DATASET CLUSTERING LINCLUST/MMseqs2.




# Leemos el fichero que contiene el clustering Linclust/MMSeqs2 realizado.
clusterRes_cluster = pd.read_csv(r'CLUSTERING_MMSEQS2/clusterRes_cluster.tsv'
                                 , sep = '\t', header = None)

clusterRes_cluster = clusterRes_cluster.rename(columns={0:'cluster-representative'})
clusterRes_cluster = clusterRes_cluster.rename(columns={1:'cluster-member'})
           

cantidadGrupo =  pd.DataFrame()
cantidadGrupo['cantidad']=clusterRes_cluster.groupby('cluster-representative').size()


# Eliminamos los 'single' péptidos, ya que son péptidos que no han formado clúster
# con ningún otro péptido.
df = cantidadGrupo.drop(cantidadGrupo[cantidadGrupo['cantidad']==1].index)
ID =  pd.DataFrame()
ID= df.index
df['ID'] = ID
clusterRes_cluster = clusterRes_cluster.rename(columns={'cluster-representative':'ID'})

res=pd.merge(df, clusterRes_cluster, on='ID')
res_2=pd.merge(res, All_peptides, on='ID')
res_2=res_2.loc[:,['ID','Sequence','Label']]


# Mergeamos los descriptores PAAC con cada secuencia peptídicia del dataset.
res_3=pd.merge(PAAC_descriptors, res_2, on='ID')


# Exportamos el fichero para reutilizarlo en los demás algoritmos.
res_3.to_csv("BBDD/peptidosLinclust.csv", 
                sep='\t', index=False, header=True)



# Preparación del dataset para la generación de los modelos predictivos.
# Vamos a quedarnos con las columnas que nos interesan del dataset creado.


Modelo_linclust=res_3.drop(columns=['Sequence'])


# Barajamos los datos.
from sklearn.utils import shuffle
Modelo_linclust = shuffle(Modelo_linclust, random_state = 0)


# Dividimos los datos en el conjunto train y test en un 80% y 20%, respectivamente.
X = Modelo_linclust.drop(columns=['ID','Label'])
Y = Modelo_linclust.iloc[:,-1]

from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test = train_test_split(X,Y, test_size = 0.20, random_state = 0)


# Generamos el MODELO SVM.
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

parameters = {
    'kernel':['rbf'],
    'gamma':[1e-3,1e-2,0.1,0.5,1],
    'C':[0.1,0.5,1, 10, 100, 1000]
}


# Realizamos el entrenamiento del modelo con cross-validation (CV) = 10.
clf = GridSearchCV(SVC(), parameters, cv=10)
model=clf.fit(x_train,y_train)

# clf.best_params_
# BEST PARAMETERS PARA EL MODELO DBSCAN_SUBSAMPLING.
# {'C': 1000, 'gamma': 0.1, 'kernel': 'rbf'}


y_pred = model.predict(x_test)
mat = confusion_matrix(y_test, y_pred)
c_SVM_Linclust=classification_report(y_test,y_pred)
TP = mat[0][0]
TN = mat[1][1]
FP = mat[0][1]
FN = mat[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)




from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(y_test, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_Linclust = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_Linclust.loc[len(resGenerales_SVM_Linclust)]=['Modelo Linclust', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_Linclust = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_Linclust.loc[len(matrizConfusion_SVM_Linclust)]=['Modelo Linclust', 'NonToxic', 'Toxic']
matrizConfusion_SVM_Linclust.loc[len(matrizConfusion_SVM_Linclust)]=['NonToxic', TP, FP]
matrizConfusion_SVM_Linclust.loc[len(matrizConfusion_SVM_Linclust)]=['Toxic', FN, TN]








# Guardamos los resultados de los todos los modelos generados.

Resultados_SVM_I=resGenerales_SVM_1descrip


Resultados_SVM= pd.ExcelWriter('RESULTADOS MODELOS/Resultados_SVM_I.xlsx')
Resultados_SVM_I.to_excel(Resultados_SVM, sheet_name="Resultados_SVM_I", index=False)


Resultados_SVM.save()
Resultados_SVM.close()



Resultados_SVM_II = resGenerales_SVM_SUBSAMPLING.append(resGenerales_SVM_DBSCAN)
Resultados_SVM_II = Resultados_SVM_II.append(resGenerales_SVM_DBSCANySUBSAMPLING)
Resultados_SVM_II = Resultados_SVM_II.append(resGenerales_SVM_Linclust)


Resultados_SVM= pd.ExcelWriter('RESULTADOS MODELOS/Resultados_SVM_II.xlsx')
Resultados_SVM_II.to_excel(Resultados_SVM, sheet_name="Resultados_SVM_II", index=False)


Resultados_SVM.save()
Resultados_SVM.close()


# Matriz de confusión.

matrizConfusion_SVM_I=matrizConfusion_SVM_1descrip

matrizConfusion_SVM= pd.ExcelWriter('RESULTADOS MODELOS/matrizConfusion_SVM_I.xlsx')
matrizConfusion_SVM_I.to_excel(matrizConfusion_SVM, sheet_name="matrizConfusion_SVM_I", index=False)

matrizConfusion_SVM.save()
matrizConfusion_SVM.close()



matrizConfusion_SVM_II = matrizConfusion_SVM_SUBSAMPLING.append(matrizConfusion_SVM_DBSCAN)
matrizConfusion_SVM_II = matrizConfusion_SVM_II.append(matrizConfusion_SVM_DBSCANySUBSAMPLING)
matrizConfusion_SVM_II = matrizConfusion_SVM_II.append(matrizConfusion_SVM_Linclust)


matrizConfusion_SVM= pd.ExcelWriter('RESULTADOS MODELOS/matrizConfusion_SVM_II.xlsx')
matrizConfusion_SVM_II.to_excel(matrizConfusion_SVM, sheet_name="matrizConfusion_SVM_II", index=False)


matrizConfusion_SVM.save()
matrizConfusion_SVM.close()


















