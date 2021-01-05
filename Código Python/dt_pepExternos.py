#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 00:41:56 2021

@author: marianomonserratgomez
"""

# DATASET DE PÉPTIDOS EXTERNOs.

# Dado un nuevo conjunto de péptidos externos, queremos ver cómo los mejores modelos
# generados clasifican dichos péptidos.


# Por ello, primero leemos el dataset con el nuevo conjunto de péptidos.
import pandas as pd
dataset_ejemplo = pd.read_csv("PEPTIDOS EXTERNOS/dataset_ejemplo.csv", sep='\t')


All_peptides = pd.read_csv(r'BBDD/All_peptides.txt', 
                           sep = ',')

# Comprobamos si el dataset inicial con todos los péptidos contiene alguno de los 
# péptidos del dataset externo. 
# Hacemos un merge, el dataframe 'merged_inner' creado, contendrá
# los péptidos que comparten los datasets, por lo que debería estar vacío.
merged_inner = pd.merge(left=dataset_ejemplo,right=All_peptides, left_on='Sequence', right_on='Sequence')




# Segundo, creamos un fichero fasta y obtenemos los descriptores PAAC de cada péptido.
# Por tanto, exportamos el fichero 'dataset_ejemplo.csv' a la carpeta 'datasets_fasta'.
# y por útlimo, utilizamos las funciones de los scripts 'crearFasta.py' y 'funcion_PAAC.py'
dataset_ejemplo.to_csv("datasets_fasta/dataset_ejemplo.csv",index=None)



# Leemos los descriptores PAAC.
PAAC_descriptors = pd.read_csv(r'BBDD/encodings_dataset_ejemplo.tsv', sep = '\t')
PAAC_descriptors.rename(columns={'#': 'ID'}, inplace=True)

# Unimos los dos datasets.
dataset_ejemplo=pd.merge(PAAC_descriptors, dataset_ejemplo, on='ID')

# Preparación del dataset externo.
prep_modelo=dataset_ejemplo.drop(columns=['Sequence'])
X_externo=prep_modelo.drop(columns=['ID','Label'])
Y_externo = prep_modelo.iloc[:,-1]



# Generamos de nuevo los mejores modelos predictivos seleccionados anteriormente.


# MODELO SUBSAMPLING CON ALGORITMO SVM.

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
# MEJORES PARÁMETROS.
# clf.best_params_
# {'C': 10, 'gamma': 0.01, 'kernel': 'rbf'}

parameters = {
    'kernel':['rbf'],
    'gamma':[0.01],
    'C':[10]
}
# Realizamos el entrenamiento del modelo con 10fold-cross-validation.
clf = GridSearchCV(SVC(), parameters, cv=10)
model_SVM_Subsampling=clf.fit(x_train,y_train)




# Realizamos el testeo con el dataset externo para ver si el modelo
# generado clasifica correctamente los péptidos de este dataset nuevo.


# Realizamos las predicciones.
y_pred_externo_SVM_Subsampling = model_SVM_Subsampling.predict(X_externo)




# Parámetros para evaluar la calidad del modelo:

mat_SVM_Subsampling = confusion_matrix(Y_externo, y_pred_externo_SVM_Subsampling)
c_SVM_Subsampling=classification_report(Y_externo,y_pred_externo_SVM_Subsampling)
TP = mat_SVM_Subsampling[0][0]
TN = mat_SVM_Subsampling[1][1]
FP = mat_SVM_Subsampling[0][1]
FN = mat_SVM_Subsampling[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)

from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(Y_externo, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_externo_SVM_Subsampling, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_SUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_SUBSAMPLING.loc[len(resGenerales_SVM_SUBSAMPLING)]=['Modelo Subsampling SVM', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_SUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_SUBSAMPLING.loc[len(matrizConfusion_SVM_SUBSAMPLING)]=['Modelo Subsampling SVM', 'NonToxic', 'Toxic']
matrizConfusion_SVM_SUBSAMPLING.loc[len(matrizConfusion_SVM_SUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_SVM_SUBSAMPLING.loc[len(matrizConfusion_SVM_SUBSAMPLING)]=['Toxic', FN, TN]







# MODELO SUBSAMPLING CON ALGORITMO RF.

from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

# model_RF.best_params_
# LOS MEJORES PARAMETROS SON:
# {'n_estimators': 200}
model_RF_Subsampling = GridSearchCV(RandomForestClassifier(),
                        param_grid={'n_estimators':[200]},
                        cv = 10)


model_RF_Subsampling=model_RF_Subsampling.fit(x_train,y_train)

# Realizamos las predicciones.
y_predictions_RF_Subsampling=model_RF_Subsampling.predict(X_externo)

# Parámetros para evaluar la calidad del modelo:
mat_RF_Subsampling = confusion_matrix(Y_externo, y_predictions_RF_Subsampling)
c_RF_Subsampling=classification_report(Y_externo,y_predictions_RF_Subsampling)
TP = mat_RF_Subsampling[0][0]
TN = mat_RF_Subsampling[1][1]
FP = mat_RF_Subsampling[0][1]
FN = mat_RF_Subsampling[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)


from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(Y_externo, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_predictions_RF_Subsampling, drop_first = True)

fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_RF = roc_auc_score(y_test_num, predictions_num)


# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_RF_SUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_RF_SUBSAMPLING.loc[len(resGenerales_RF_SUBSAMPLING)]=['Modelo Subsampling RF', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_RF,2)]


matrizConfusion_RF_SUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_RF_SUBSAMPLING.loc[len(matrizConfusion_RF_SUBSAMPLING)]=['Modelo Subsampling RF', 'NonToxic', 'Toxic']
matrizConfusion_RF_SUBSAMPLING.loc[len(matrizConfusion_RF_SUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_RF_SUBSAMPLING.loc[len(matrizConfusion_RF_SUBSAMPLING)]=['Toxic', FN, TN]









# MODELO DEBSCAN+SUBSAMPLING CON ALGORITMO SVM



DBSCAN_subsampling = pd.read_csv(r'BBDD/DBSCAN_subsampling.txt', sep = ',')
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
# clf.best_params_
# BEST PARAMETERS PARA EL MODELO DBSCAN_SUBSAMPLING.
# {'C': 100, 'gamma': 0.01, 'kernel': 'rbf'}

parameters = {
    'kernel':['rbf'],
    'gamma':[0.01],
    'C':[100]
}
# Realizamos el entrenamiento del modelo con cross-validation (CV) = 10.
clf = GridSearchCV(SVC(), parameters, cv=10)
model_DBSCANySubsampling_SVM=clf.fit(x_train,y_train)


# Realizamos las predicciones.
y_pred_DBSCANySubsampling_SVM = model_DBSCANySubsampling_SVM.predict(X_externo)

# Parámetros para evaluar la calidad.
mat_DBSCANySubsampling_SVM = confusion_matrix(Y_externo, y_pred_DBSCANySubsampling_SVM)
c_SVM_DBSCANySubsampling=classification_report(Y_externo,y_pred_DBSCANySubsampling_SVM)
TP = mat_DBSCANySubsampling_SVM[0][0]
TN = mat_DBSCANySubsampling_SVM[1][1]
FP = mat_DBSCANySubsampling_SVM[0][1]
FN = mat_DBSCANySubsampling_SVM[1][0]
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
F1 = 2*TP/(2*TP+FP+FN)
from sklearn.metrics import roc_curve, roc_auc_score
y_test_num = pd.get_dummies(Y_externo, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_DBSCANySubsampling_SVM, drop_first = True)
fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# calculate AUC
auc_SVM = roc_auc_score(y_test_num, predictions_num)



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_SVM_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_SVM_DBSCANySUBSAMPLING.loc[len(resGenerales_SVM_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+SUBSAMPLING SVM', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_SVM,2)]


matrizConfusion_SVM_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_SVM_DBSCANySUBSAMPLING.loc[len(matrizConfusion_SVM_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+SUBSAMPLING SVM', 'NonToxic', 'Toxic']
matrizConfusion_SVM_DBSCANySUBSAMPLING.loc[len(matrizConfusion_SVM_DBSCANySUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_SVM_DBSCANySUBSAMPLING.loc[len(matrizConfusion_SVM_DBSCANySUBSAMPLING)]=['Toxic', FN, TN]





# MODELO DEBSCAN+SUBSAMPLING CON ALGORITMO GBRT



# Entrenamos el modelo.
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
# MEJOR PARÁMETRO.
# clf.best_params_
# {'learning_rate': 0.1, 'n_estimators': 200}
param_grid = {'n_estimators' : [200],
              'learning_rate' : [0.1]}
clf = GridSearchCV(GradientBoostingClassifier(random_state=0), param_grid, cv=10)
model_GBRT_DBSCANySubsampling=clf.fit(x_train,y_train) 

# Predicciones
y_pred_GBRT_DBSCANySubsampling = model_GBRT_DBSCANySubsampling.predict(X_externo)


# Parámetros para evaluar el modelo.
from sklearn.metrics import classification_report, confusion_matrix
mat_Grad = confusion_matrix(Y_externo, y_pred_GBRT_DBSCANySubsampling)
c_GBRT_DBSCANySubsampling=classification_report(Y_externo,y_pred_GBRT_DBSCANySubsampling)
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
y_test_num = pd.get_dummies(Y_externo, columns = ["Label"], drop_first = True)
predictions_num = pd.get_dummies(y_pred_GBRT_DBSCANySubsampling, drop_first = True)
fpr, tpr, thresholds = roc_curve(y_test_num, predictions_num)
# Calculamos AUC
auc_GBRT = roc_auc_score(y_test_num, predictions_num)



resGenerales_GBRT_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_GBRT_DBSCANySUBSAMPLING.loc[len(resGenerales_GBRT_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+Subsampling GBRT', round(accuracy_Grad,2),
                                           round(precision_Grad,2), round(sensitivity_Grad,2), 
                                           round(specifity_Grad,2),
                                           round(F1_GBRT,2), round(auc_GBRT,2)]



matrizConfusion_GBRT_DBSCANySUBSAMPLING = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_GBRT_DBSCANySUBSAMPLING.loc[len(matrizConfusion_GBRT_DBSCANySUBSAMPLING)]=['Modelo DBSCAN+Subsampling GBRT', 'NonToxic', 'Toxic']
matrizConfusion_GBRT_DBSCANySUBSAMPLING.loc[len(matrizConfusion_GBRT_DBSCANySUBSAMPLING)]=['NonToxic', TP, FP]
matrizConfusion_GBRT_DBSCANySUBSAMPLING.loc[len(matrizConfusion_GBRT_DBSCANySUBSAMPLING)]=['Toxic', FN, TN]






# MODELO LINCLUST CON ALGORITMO RF


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


# Entrenamos el modelo.
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
# model_RF.best_params_
# LOS MEJORES PARAMETROS SON:
# {'n_estimators': 500}
model_RF_Linclust = GridSearchCV(RandomForestClassifier(),
                        param_grid={'n_estimators':[500]},
                        cv = 10)


model_RF_Linclust=model_RF_Linclust.fit(x_train,y_train)




# Predicciones.
predictions_RF_Linclust=model_RF_Linclust.predict(X_externo)
a=pd.DataFrame(predictions_RF_Linclust)
from sklearn.metrics import classification_report, confusion_matrix

mat_RF_Linclust = confusion_matrix(Y_externo, predictions_RF_Linclust)
c_RF_Linclust=classification_report(Y_externo,predictions_RF_Linclust)
TP = mat_RF_Linclust[0][0]
TN = mat_RF_Linclust[1][1]
FP = mat_RF_Linclust[0][1]
FN = mat_RF_Linclust[1][0]
F1 = 2*TP/(2*TP+FP+FN)
accuracy =(TP+TN)/(TP+TN+FN+FP)
precision = TP/(TP+FP)
sensitivity = TP/(TP+FN)
specifity = TN/(TN+FP)
from sklearn.metrics import roc_curve, roc_auc_score

y_test_num = pd.get_dummies(Y_externo, columns = ["Label"])
predictions_num = pd.get_dummies(predictions_RF_Linclust)
fpr, tpr, thresholds = roc_curve(y_test_num['Toxic'], predictions_num['Toxic'])
# calculate AUC
auc_RF = roc_auc_score(y_test_num['Toxic'], predictions_num['Toxic'])



# Creamos dataframes donde guardar los resultados que exportaremos en un archivo Excel.
resGenerales_RF_Linclust = pd.DataFrame(columns=('Modelo','Accuracy',
                                        'Precision','Sensitivity',
                                        'Sprecifity', 'F1', 'AUC'))

resGenerales_RF_Linclust.loc[len(resGenerales_RF_Linclust)]=['Modelo Linclust RF', round(accuracy,2),
                                           round(precision,2), round(sensitivity,2), 
                                           round(specifity,2),
                                           round(F1,2), round(auc_RF,2)]


matrizConfusion_RF_Linclust = pd.DataFrame(columns=('Modelo','',
                                        ''))

matrizConfusion_RF_Linclust.loc[len(matrizConfusion_RF_Linclust)]=['Modelo Linclust RF', 'NonToxic', 'Toxic']
matrizConfusion_RF_Linclust.loc[len(matrizConfusion_RF_Linclust)]=['NonToxic', TP, FP]
matrizConfusion_RF_Linclust.loc[len(matrizConfusion_RF_Linclust)]=['Toxic', FN, TN]






# Vamos a exportar los dataframe creados en dos ficheros Excel.

Resultados_datasetExterno_I = resGenerales_SVM_SUBSAMPLING.append(resGenerales_RF_SUBSAMPLING)
Resultados_datasetExterno_I = Resultados_datasetExterno_I.append(resGenerales_SVM_DBSCANySUBSAMPLING)
Resultados_datasetExterno_I = Resultados_datasetExterno_I.append(resGenerales_GBRT_DBSCANySUBSAMPLING)
Resultados_datasetExterno_I = Resultados_datasetExterno_I.append(resGenerales_RF_Linclust)

Resultados_datasetExterno= pd.ExcelWriter('RESULTADOS MODELOS/Resultados_datasetExterno.xlsx')
Resultados_datasetExterno_I.to_excel(Resultados_datasetExterno, sheet_name="Resultados_datasetExterno", index=False)

Resultados_datasetExterno.save()
Resultados_datasetExterno.close()


matrizConfusion_datasetExterno_I = matrizConfusion_SVM_SUBSAMPLING.append(matrizConfusion_RF_SUBSAMPLING)
matrizConfusion_datasetExterno_I = matrizConfusion_datasetExterno_I.append(matrizConfusion_SVM_DBSCANySUBSAMPLING)
matrizConfusion_datasetExterno_I = matrizConfusion_datasetExterno_I.append(matrizConfusion_GBRT_DBSCANySUBSAMPLING)
matrizConfusion_datasetExterno_I = matrizConfusion_datasetExterno_I.append(matrizConfusion_RF_Linclust)

matrizConfusion_datasetExterno= pd.ExcelWriter('RESULTADOS MODELOS/matrizConfusion_datasetExterno.xlsx')
matrizConfusion_datasetExterno_I.to_excel(matrizConfusion_datasetExterno, sheet_name="matrizConfusion_datasetExterno", index=False)

matrizConfusion_datasetExterno.save()
matrizConfusion_datasetExterno.close()






