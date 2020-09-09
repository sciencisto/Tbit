#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 13:53:53 2020

@author: kristine
"""
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import r2_score
from sklearn import linear_model
from scipy.special import expit

def import_filescdr12(path):   
    d = {"TM": [1], "RMSD": [1], "Blosum62": [1], "hydrophobicity": [1], "Bin": [0]}
    df = pd.DataFrame(data=d)    
    file = pd.read_csv(path, index_col=None, header = None)
    file.columns = df.columns
    df = pd.concat([df, file])
    return df


def import_filescdr3(path):   
    d = {"Hamming": [1], "Blosum62": [1], "Blosum45": [1], "Pam10": [1], "K-Mer_3": [1], "K-Mer_4": [1], 
         "K-Mer_5": [1], "Atchley_Factors": [1], "hydrophobicity": [1], "Bin": [0]}
    df = pd.DataFrame(data=d)    
    file = pd.read_csv(path, index_col=None, header = None)
    file.columns = df.columns
    df = pd.concat([df, file])
    return df
#sub_CDR3 = import_files('/home/kristine/Desktop/CDR3_subsampled.csv')
#sub_CDR1 = import_filescdr12('/home/kristine/Desktop/CDR1_subsampled.csv')
#sub_CDR2 = import_filescdr12('/home/kristine/Desktop/CDR2_subsampled.csv')

def graph_auc(df):
    AUC_list = []
    metho_list = list(df.columns[:4])
    y_true = np.array(df.Bin).flatten()
    
    for i in range(len(metho_list)):
        y_scores = np.array(df[metho_list[i]]).flatten()
        AUC = roc_auc_score(y_true, 1-y_scores)
        AUC_list.append(AUC)
        
        fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores, pos_label = 1)
        #plt.plot(thresholds, 1-tpr)
        #plt.plot(thresholds, fpr)

        plt.plot(1-fpr, 1-tpr, label = "%s"% (metho_list[i]))    

        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive rate')
        plt.title("Reciver Operator Characteristics curves for CDR2")
        plt.legend(loc='lower right')
        
        
    return AUC_list    
    
test = graph_auc(sub_CDR1)


"""
def log_Regress(total):
    y = total.Bin
    X = 1-total.drop("Bin", axis=1)
    
    clf= LogisticRegression().fit(X, y)
    
    plt.figure(1, figsize=(4, 3))
    plt.clf()
    plt.scatter(X.ravel(), y, color='black', zorder=20)
    X_test = np.linspace(-5, 10, 300)

    loss = expit(X_test * clf.coef_ + clf.intercept_).ravel()
    plt.plot(X_test, loss, color='red', linewidth=3)

    ols = linear_model.LinearRegression()
    ols.fit(X, y)
    plt.plot(X_test, ols.coef_ * X_test + ols.intercept_, linewidth=1)
    plt.axhline(.5, color='.5')

    plt.ylabel('y')
    plt.xlabel('X')
    plt.xticks(range(-5, 10))
    plt.yticks([0, 0.5, 1])
    plt.ylim(-.25, 1.25)
    plt.xlim(-4, 10)
    plt.legend(('Logistic Regression Model', 'Linear Regression Model'),
           loc="lower right", fontsize='small')
    plt.tight_layout()
    plt.show()
    
    
    pred_y_2 = clf.predict(X)
    
    print("Classes predicted")
    print( np.unique( pred_y_2 ) )
    print("R squared")
    print( r2_score(y, pred_y_2))
    print("accuray score:")
    print( accuracy_score(y, pred_y_2) )
    print("Coeficients:")
    print(clf.coef_)
    print("intercept")
    print(clf.intercept_)

    return clf


#test1 = log_Regress(sub_CDR3)


plt.scatter(X.Hamming,y, color='black', zorder=20)
X_test = np.linspace(0, 1, 300)
loss = expit(X_test * 2.73918869 + -2.34217803).ravel()
plt.plot(X_test, loss, color='darkred', linewidth=3)
ols = linear_model.LinearRegression()
ols.fit(X, y)
plt.plot(X_test, ols.coef_ * X_test + ols.intercept_, linewidth=1)
plt.axhline(.5, color='.5')



y_scores = np.array(df[metho_list[i]]).flatten()
fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores, pos_label = 1)
        






y_true = np.array(sub_CDR3.Bin).flatten()
y_scores = np.array(sub_CDR3.Pam10).flatten()

fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores, pos_label = 1)
plt.plot(thresholds, 1-tpr)
plt.plot(thresholds, fpr)
plt.plot(1-fpr, 1-tpr)

AUC = roc_auc_score(y_true, 1-y_scores)
print(AUC)


# General a toy dataset:s it's just a straight line with some Gaussian noise:


#plt.scatter(1-X["K-Mer_3"],y, color='black', zorder=20)
#plt.scatter(1-X["K-Mer_4"],y, color='black', zorder=20)
#plt.scatter(1-X["K-Mer_5"],y, color='black', zorder=20)
#plt.scatter(1-X["Hamming"],y, color='black', zorder=20)
#plt.scatter(1-X["Blosum62"],y, color='black', zorder=20)
#plt.scatter(1-X["Blosum45"],y, color='black', zorder=20)
plt.scatter(1-X["Pam10"],y, color='black', zorder=20)
#plt.scatter(1-X["Atchley_Factors"],y, color='black', zorder=20)
#plt.scatter(1-X["hydrophobicity"],y, color='black', zorder=20)

X_test = np.linspace(0, 1, 300)
loss = expit(X_test * -1.74 + 2.3).ravel()
plt.plot(X_test, loss, color='darkred', linewidth=3)


#2.73918869,  1.08830587, -0.69914219, -1.74655054,  0.17892066 
# 3.54059637, -4.46657442, -0.5527015 ,  0.13295723

new_line = np.array([1-X["Hamming"]*-2.74 + 1-X["Blosum62"]*-1.08 + 
                     1-X["Blosum45"]*0.70 + 1-X["Pam10"]*1.74 + 1-X["K-Mer_3"]*-0.17 
                     + 1-X["K-Mer_4"]*-3.5 + 1-X["K-Mer_5"]*4.46 +1-X["Atchley_Factors"]*0.55
                     + 1-X["hydrophobicity"]*-0.13])

"""