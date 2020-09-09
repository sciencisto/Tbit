#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 08:26:05 2020

@author: kristine
"""
import pandas as pd 
import numpy as np
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn import metrics
from sys import argv 
import time 
import seaborn as sns
sns.set_style("white")
sns.set_context("paper")

def parse_commandline():  
    """
    cmmandline call: python3 ROC.py path 
	python3 ROC.py "/work1/s111518/CDR2dist_TM/Train/TCR_subsampled.csv" 
	"""
    path_dist = argv[1] #path to the subsampled data

    return path_dist

def import_file(path):
    BetaDist_sub = pd.read_csv(path, header = None)
    BetaDist_sub.columns = ["BetaDist", "Bin"]
    return BetaDist_sub
    

def roc_data(df):
    y_true = np.array(df["Bin"]).flatten()
    y_scores = np.array(df["BetaDist"]).flatten()   
    AUC = roc_auc_score(y_true, 1-y_scores) 
    fpr, tpr, thresholds = metrics.roc_curve(y_true, 1-y_scores, pos_label = 1)
    print("BetaDist AUC is:")
    print(AUC)
    return fpr, tpr, thresholds


def ROC_graph(TCR_v1_ROC):    
    #graph
    plt.figure(figsize=(10,7))

    plt.plot(np.array([0,1]), np.array([0,1]), '--', c = "green", alpha = 0.6, label ="Random Guess")
    
    plt.plot(TCR_v1_ROC[0], TCR_v1_ROC[1], '-o', c = "red", alpha = 0.7, label = "BetaDist")

    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive rate')
    plt.title("Reciver Operator Characteristics curve")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5));
    plt.savefig("ROC.png", bbox_inches = 'tight')

def pairplot(file):
    g = sns.pairplot(file, hue="Bin", palette="dark")
    g.savefig("pairplot.png", bbox_inches = 'tight')


if __name__ == '__main__':
    start_overall_time = time.time()   
    path_dist = parse_commandline()
    file = import_file(path_dist)
    AUC = roc_data(file)
    graph = ROC_graph(AUC)
    pair = pairplot(file)
  
    print("--- %s seconds ---" % (time.time()- start_overall_time))
