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


def GetBestModel(X,y,pen='l2',Cs = 10, cv=10, tol = 1e-2):
    L_Cv = LogisticRegressionCV(Cs = Cs, cv=cv, penalty=pen, solver = 'saga',
                                n_jobs = -1, tol = tol)
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
                      Thresh = 0.9, DT = False):

    prev = y + 0.0
    y_orig = y + 0.0

    NoPos = []
    for i in range(Iter):
        #creating training matrices at each iteration

        tmp = BalanceClasses(X + 0.0 ,y + 0.0)
        X_new, y_new = tmp['X'],tmp['y']


        LR = LogisticRegression(penalty = pen, solver = 'saga', C = C)
        LR.fit(X_new,y_new)

        y_pred = LR.predict_proba(X)
        y_pred = y_pred[:,1]

        y = y_orig + 0.0
        y[y_pred >=Thresh] = 1.0

        if sum(abs(y - prev)) == 0:
            break

        prev = y + 0.0
        NoPos += [sum(y)]

    return {'LR':LR, 'y':y, 'NoPos':NoPos}


def CrossValScore(Mdl,X,y):
    c, r = y.shape
    y = y.reshape(c,)
    Prec = cross_val_score(Mdl, X, y, cv=10, scoring='precision', n_jobs = -1)
    Rec = cross_val_score(Mdl, X, y, cv=10, scoring='recall', n_jobs = -1)
    AUC = cross_val_score(Mdl, X, y, cv=10, scoring='roc_auc', n_jobs = -1)

    return {'Prec':Prec, 'Rec':Rec, 'AUC':AUC}



def IterativeCotraining(X1,X2,X3,y,C1 = 1.0, C2 = 1.0, C3= 1.0, Iter=10,rho=0.1,pen='l2', Thresh = 0.9, BC = True):

    #balance feature classes

    if BC:
        tmp = BalanceClasses(X1 + 0.0 ,y + 0.0)
        X1_new, y1_new = tmp['X'],tmp['y']

        tmp = BalanceClasses(X2 + 0.0 ,y + 0.0)
        X2_new, y2_new = tmp['X'],tmp['y']

        tmp = BalanceClasses(X3 + 0.0 ,y + 0.0)
        X3_new, y3_new = tmp['X'],tmp['y']

    else:
        X1_new = X1 + 0.0
        y1_new = y + 0.0

        X2_new = X2 + 0.0
        y2_new = y + 0.0

        X3_new = X3 + 0.0
        y3_new = y + 0.0




    n = len(y1_new)
    n1 = len(y)
    m = n + 0.0


    y1_new = np.reshape(y1_new,(n,1))
    y2_new = np.reshape(y2_new,(n,1))
    y3_new = np.reshape(y3_new,(n,1))

    y_orig = y1_new + 0.0

    y1_rec = y + 0.0
    y2_rec = y + 0.0
    y3_rec = y + 0.0


    y1 = y + 0.0
    y2 = y + 0.0
    y3 = y + 0.0

    for i in range(Iter):

        LR1 = LogisticRegression(penalty = pen, solver = 'saga', C = C1)
        LR1.fit(X1_new,y1_new)

        LR2 = LogisticRegression(penalty = pen, solver = 'saga', C = C2)
        LR2.fit(X2_new,y2_new)

        LR3 = LogisticRegression(penalty = pen, solver = 'saga', C = C3)
        LR3.fit(X3_new,y3_new)


        a1 = (1.0/(rho))*(np.asscalar(LR1.intercept_) + X1.dot(np.transpose(LR1.coef_))) + 0.5*(y2 + y3 + 0.0)
        a2 = (1.0/(rho))*(np.asscalar(LR2.intercept_) + X2.dot(np.transpose(LR2.coef_))) + 0.5*(y1 + y3 + 0.0)
        a3 = (1.0/(rho))*(np.asscalar(LR3.intercept_) + X3.dot(np.transpose(LR3.coef_))) + 0.5*(y1 + y2 + 0.0)


        #print(a1.shape)

        y1 = np.maximum(np.zeros((n1,1)),np.minimum(np.ones((n1,1)),a1)) + 0.0
        y2 = np.maximum(np.zeros((n1,1)),np.minimum(np.ones((n1,1)),a2)) + 0.0
        y3 = np.maximum(np.zeros((n1,1)),np.minimum(np.ones((n1,1)),a3)) + 0.0

        y1[y1 >= Thresh] = 1.0
        y1[y==1.0] = 1.0
        y1[y1 < Thresh] = 0.0

        y2[y2 >= Thresh] = 1.0
        y2[y==1.0] = 1.0
        y2[y2 < Thresh] = 0.0

        y3[y3 >= Thresh] = 1.0
        y3[y==1.0] = 1.0
        y3[y3 < Thresh] = 0.0

        if BC:
            tmp = BalanceClasses(X1 + 0.0 ,y1 + 0.0)
            X1_new, y1_new = tmp['X'],tmp['y']

            tmp = BalanceClasses(X2 + 0.0 ,y2 + 0.0)
            X2_new, y2_new = tmp['X'],tmp['y']

            tmp = BalanceClasses(X3 + 0.0 ,y3 + 0.0)
            X3_new, y3_new = tmp['X'],tmp['y']

        else:
            X1_new = X1 + 0.0
            y1_new = y + 0.0

            X2_new = X2 + 0.0
            y2_new = y3 + 0.0

            X3_new = X3 + 0.0
            y3_new = y3 + 0.0



        y1_rec = np.hstack((y1_rec,np.reshape(LR1.predict_proba(X1)[:,1],(n1,1))))
        y2_rec = np.hstack((y2_rec,np.reshape(LR2.predict_proba(X2)[:,1],(n1,1))))
        y3_rec = np.hstack((y3_rec,np.reshape(LR3.predict_proba(X3)[:,1],(n1,1))))


    y_rec = [y1_rec,y2_rec,y3_rec]

    LR_all = [LR1,LR2,LR3]

    return {'LR':LR_all,'y_rec':y_rec}
