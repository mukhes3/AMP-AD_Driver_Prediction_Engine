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



FileOut = 'EC2_ConsUnion1_NoDe_1015.csv'
FileOut2 = 'EC2_ConsUnion1_NoDeProbs_DLPFC_PWB.csv'

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
X2 = pd.read_csv('./amp_ad_xQTL_CIT_feature_set_DLPFC.csv')
X3 = pd.read_csv('./amp_ad_subNet_feature_set_GW.csv')

GeneId = X1['GeneID']

#Processing datasets
X0 = ProcessData(X0)
X1 = ProcessData(X1)
X2 = ProcessData(X2)
X3 = ProcessData(X3)

#Reading response vector
Y1 = pd.read_excel('./ResponseVec_SWB1000G.xlsx')
Y1 = Y1.drop(['GeneID'], axis = 1)

Y2 = pd.read_excel('./ResponseVec_LS.xlsx')
Y2 = Y2.drop(['GeneID'], axis = 1)

Y3 = pd.read_excel('./ResponseVec_SWB.xlsx')
Y3 = Y3.drop(['GeneID'], axis = 1)

Y4 = pd.read_excel('./ResponseVec_PA.xlsx')
Y4 = Y4.drop(['GeneID'], axis = 1)

Y = (Y1 + Y2 + Y3 + Y4)>0 + 0.0

#Crossvalidation and training
def TrainingCrossVal(X_normed,Y,Pen):
    #Balancing classes for better classification performance
    tmp = SF.BalanceClasses(X_normed + 0.0 ,Y.values + 0.0)
    X_train, y_train = tmp['X'],tmp['y']

    print 'Crossvalidation'
    Cs = np.logspace(-4,4,10)
    Cs = Cs.tolist()
    C_best = SF.GetBestModel(X_train,y_train,pen=Pen,Cs = Cs)
    C_best = C_best[0]
    print C_best
    print 'IterativeTraining'
    D = SF.IterativeTraining(X_normed + 0.0, Y.values + 0.0 ,C = C_best,pen=Pen,
                            Iter = 10,Thresh = 0.9)
    return(D)


print 'training X0'
D0a = TrainingCrossVal(X0,Y,'l1')
D0b = TrainingCrossVal(X0,Y,'l2')
print 'training X1'
D1a = TrainingCrossVal(X1,Y,'l1')
D1b = TrainingCrossVal(X1,Y,'l2')
print 'training X2'
D2a = TrainingCrossVal(X2,Y,'l1')
D2b = TrainingCrossVal(X2,Y,'l2')
print 'training X3'
D3a = TrainingCrossVal(X3,Y,'l1')
D3b = TrainingCrossVal(X3,Y,'l2')


#Getting consensus predictions
def PredictProbs(D,X):
    Y = D['LR'].predict_proba(X)
    Y = Y[:,1]
    return(Y)

Y0a = PredictProbs(D0a,X0)
Y0b = PredictProbs(D0b,X0)

Y1a = PredictProbs(D1a,X1)
Y1b = PredictProbs(D1b,X1)

Y2a = PredictProbs(D2a,X2)
Y2b = PredictProbs(D2b,X2)

Y3a = PredictProbs(D3a,X3)
Y3b = PredictProbs(D3b,X3)

#getting weights for each model
#W0a = D0a['LR'].coef_
#W0b = D0b['LR'].coef_

#W1a = D1a['LR'].coef_
#W1b = D1b['LR'].coef_

#W2a = D2a['LR'].coef_
#W2b = D2b['LR'].coef_

#W3a = D3a['LR'].coef_
#W3b = D3b['LR'].coef_

Y1 = (Y0a + Y0b + Y1a + Y1b + Y2a + Y2b + Y3a + Y3b)/8.0
#Y1 = (Y3a + Y3b)/4.0

GP = pd.DataFrame(data = {'Gene':GeneId,
'Y0a':Y0a,'Y0b':Y0b,
'Y1a':Y1a,'Y1b':Y1b,
'Y2a':Y2a,'Y2b':Y2b,
'Y3a':Y3a,'Y3b':Y3b})

#GP = pd.DataFrame(data = {'Gene':GeneId,
#'Y3a':Y3a,'Y3b':Y3b})

#W0 = pd.DataFrame(data = {'L1':W0a,'L2':W0b})
#W1 = pd.DataFrame(data = {'L1':W1a,'L2':W1b})
#W2 = pd.DataFrame(data = {'L1':W2a,'L2':W2b})
#W3 = pd.DataFrame(data = {'L1':W3a,'L2':W3b})

#W0.to_csv('ConsUnion_Wts_0904_x0.csv')
#W1.to_csv('ConsUnion_Wts_0904_x1.csv')
#W2.to_csv('ConsUnion_Wts_0904_x2.csv')
#W3.to_csv('ConsUnion_Wts_0904_x3.csv')

GP.to_csv(FileOut2)
