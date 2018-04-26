# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:08:33 2018

@author: Sumit
"""



import numpy as np 
from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.model_selection import cross_val_score

def GetBestModel(X,y,pen='l2'):
    L_Cv = LogisticRegressionCV(Cs = 10, cv=10, penalty=pen, solver = 'saga')
    L_Cv.fit(X, y)
    C_best = L_Cv.C_    
    return C_best

def BalanceClasses(X,y):
    In = np.argwhere(y)
    In = In[:,0]
    X_new = X
    Reps = (len(y)-len(In))/len(In)
    Y_new = y
    for i in range(Reps):
        X_new = np.vstack((X_new,X[In]))
        Y_new = np.vstack((Y_new,np.ones((len(In),1))))
        
    return {'X':X_new, 'y':Y_new}
        

def IterativeTraining(X,y,C = 1.0,pen='l2',Iter = 10, 
                      Thresh = 0.9):
    
    prev = y + 0.0  
    
    NoPos = []
    for i in range(Iter):
        #creating training matrices at each iteration
        
        tmp = BalanceClasses(X + 0.0 ,y + 0.0)
        X_new, y_new = tmp['X'],tmp['y']
        
        LR = LogisticRegression(penalty = pen, solver = 'saga', C = C)
        LR.fit(X_new,y_new)
        
        y_pred = LR.predict_proba(X)
        y_pred = y_pred[:,1]
        
        y[y_pred >=Thresh] = 1.0
        
        if sum(abs(y - prev)) == 0: 
            break 
        
        prev = y + 0.0 
        NoPos += [sum(y)]       
               
    return {'LR':LR, 'y':y, 'NoPos':NoPos}  
               
    
def CrossValScore(Mdl,X,y):
    c, r = y.shape
    y = y.reshape(c,)    
    Prec = cross_val_score(Mdl, X, y, cv=10, scoring='precision')
    Rec = cross_val_score(Mdl, X, y, cv=10, scoring='recall')
    AUC = cross_val_score(Mdl, X, y, cv=10, scoring='roc_auc')
    
    return {'Prec':Prec, 'Rec':Rec, 'AUC':AUC}

    
