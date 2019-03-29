from __future__ import division
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import manifold
from scipy import sparse
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import euclidean_distances
import SemisupFns as SF
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.model_selection import cross_val_score
import MiscFcns as MF
import scipy.stats as ss
import mygene



FileOut = 'EC2_ICCT_subNet_de_agg.csv'
FileOut2 = 'EC2_ICCT_subNet_de_agg_probs.csv'

#Function to process datasets
def ProcessData(X):
    X = X.drop(['GeneID'], axis = 1)
    X.head()
    X = X.fillna(0)
    X2 = X.abs()
    X_normed = X.values / X2.values.max(axis=0)
    return(X_normed)

#Initially loading datasets
X0 = pd.read_csv('./amp_ad_de_feature_set.csv')
X1 = pd.read_csv('./amp_ad_agg_feature_set.csv')
#X1 = pd.read_csv('./amp_ad_de_feature_set.csv')
#X0 = pd.read_csv('./amp_ad_agg_feature_set.csv')
X2 = pd.read_csv('./amp_ad_deVal_feature_set.csv')
X3 = pd.read_csv('./amp_ad_subNet_feature_set.csv')

GeneId = X1['GeneID']

#Processing datasets
X0 = ProcessData(X0)
X1 = ProcessData(X1)
#X2 = ProcessData(X2)
X3 = ProcessData(X3)

#Reading response vector
Y = pd.read_excel('./ResponseVec_040318.xlsx')
Y = Y.drop(['GeneID'], axis = 1)
Y.head()

#Crossvalidation and training
def TrainingCrossVal(X_normed,Y,Pen=''):
    #Balancing classes for better classification performance
    tmp = SF.BalanceClasses(X_normed + 0.0 ,Y.values + 0.0)
    X_train, y_train = tmp['X'],tmp['y']

    print 'Crossvalidation'
    Cs = np.logspace(-4,4,10)
    Cs = Cs.tolist()
    C_best = SF.GetBestModel(X_train,y_train,pen=Pen,Cs = Cs)
    C_best = C_best[0]
    print C_best
    #print 'IterativeTraining'
    #D = SF.IterativeTraining(X_normed + 0.0, Y.values + 0.0 ,C = C_best,pen=Pen,
    #                        Iter = 10,Thresh = 0.9)
    return(C_best)


def HoldOneOut(X0,X1,X3,Y,C1,C2,C3,rho):
    n = sum(Y)
    print(Y.shape)

    In = np.argwhere(Y==1.0)

    IC = []
    ICCT = []

    for i in range(len(In)):
        X0_i = np.delete(X0, In[i][0], 0)
        X1_i = np.delete(X1, In[i][0], 0)
        X3_i = np.delete(X3, In[i][0], 0)
        Y_i = np.delete(Y, In[i][0], 0)
        print(X3.shape)
        print(X3_i.shape)
        print(Y_i.shape)

        d = SF.IterativeCotraining(X0_i,X1_i,X3_i,Y_i,C1 = C0, C2 = C1, C3= C3, rho=rho)
        d0 = SF.IterativeTraining(X0_i + 0.0, Y_i ,C = C0,Iter = 10,Thresh = 0.9)
        d1 = SF.IterativeTraining(X1_i + 0.0, Y_i ,C = C1,Iter = 10,Thresh = 0.9)
        d3 = SF.IterativeTraining(X3_i + 0.0, Y_i ,C = C3,Iter = 10,Thresh = 0.9)

        t_IC = (d0['LR'].predict_proba(X0)[In[i][0],1] +
        d1['LR'].predict_proba(X1)[In[i][0],1] +
        d3['LR'].predict_proba(X3)[In[i][0],1])/3.0

        t_ICCT = (d['LR'][0].predict_proba(X0)[In[i][0],1] +
        d['LR'][1].predict_proba(X1)[In[i][0],1] +
        d['LR'][2].predict_proba(X3)[In[i][0],1])/3.0

        IC += [t_IC]
        ICCT += [t_ICCT]


    return({'IC':IC,'ICCT':ICCT})

C0 = TrainingCrossVal(X0,Y,'l2')
C1 = TrainingCrossVal(X1,Y,'l2')
C3 = TrainingCrossVal(X3,Y,'l2')

d = HoldOneOut(X0,X1,X3,Y.values,C1 = C0, C2 = C1, C3= C3, rho=1)

df = pd.DataFrame.from_dict(d)
df.to_csv('ICCT_hold1out.csv')
