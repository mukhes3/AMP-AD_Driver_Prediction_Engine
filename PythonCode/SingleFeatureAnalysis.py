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



FileOut = 'EC2_de_pred.csv'
FileOut2 = 'EC2_de_probs.csv'

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

GeneId = X0['GeneID']

#Processing datasets
X0 = ProcessData(X0)

#Reading response vector
Y = pd.read_excel('./ResponseVec_040318.xlsx')
Y = Y.drop(['GeneID'], axis = 1)
Y.head()

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


#Getting consensus predictions
def PredictProbs(D,X):
    Y = D['LR'].predict_proba(X)
    Y = Y[:,1]
    return(Y)

Y0a = PredictProbs(D0a,X0)
Y0b = PredictProbs(D0b,X0)


Y1 = (Y0a + Y0b)/2.0

GP = pd.DataFrame(data = {'Gene':GeneId,
'Y0a':Y0a,'Y0b':Y0b})


GP.to_csv(FileOut2)

#Getting list of genes which have SNPs according to IGAP_stage1
cnt = 0
GeneNames = []
MinPval = []
AvgLogPval = []
with open("./IGAP_geneAnalysis.genes.raw") as f:
    for line in f:
        if cnt <= 1:
            cnt+=1
            continue

        cnt += 1
        temp = line.split(' ')

        if len(temp)>9:
            GeneNames += [temp[0]]
            temp = map(float,temp[9:])
            MinPval += [min(temp)]
            AvgLogPval += [np.mean(np.log10(temp))]

#Getting the ENSEMBL ID's of the predicted regulator genes
#lgcl = (Y0a>0.5) + (Y0b>0.5)
lgcl = (Y0b>0.9)
bla = np.argwhere(lgcl)
bla2 = np.argwhere(np.logical_not(lgcl))
GenePred = list(GeneId[bla[:,0]])
GenePredVal = list(Y1[bla[:,0]])
#Converting to Gene Symbols
GenePred2 = MF.ConvertToSymb(GenePred)
GeneNotPred2 = MF.ConvertToSymb(list(GeneId[bla2[:,0]]))
GeneNotPredVal = list(Y1[bla2[:,0]])

IGAP_genes = pd.DataFrame(data = {'Genes':GeneNames,'Pval':MinPval,'MeanPval':AvgLogPval})


def PredGenesPval(IGAP_genes,GenePred, GenePredVal):
    Int = list(set(GeneNames).intersection(GenePred))
    GenePred = pd.DataFrame({'G':GenePred})

    G = []
    P = []
    Pmean = []
    Y = []

    for i in range(len(Int)):
        In = IGAP_genes['Genes'][IGAP_genes['Genes']==Int[i]].index[0]
        In2 = GenePred['G'][GenePred['G']==Int[i]].index[0]
        G += [IGAP_genes['Genes'][In]]
        P += [IGAP_genes['Pval'][In]]
        Pmean += [IGAP_genes['MeanPval'][In]]
        Y += [GenePredVal[In2]]


    OR = np.log10(np.array(Y)/(1-np.array(Y)))
    PredGenes_pval = pd.DataFrame(data = {'GeneSymb':G,'Pval':P,'MeanPval':Pmean,'Y':Y,'OR':OR})
    return PredGenes_pval



#Getting p-values of predicted genes
PredGenes_pval = PredGenesPval(IGAP_genes,GenePred2, GenePredVal)
NotPredGenes_pval = PredGenesPval(IGAP_genes,GeneNotPred2, GeneNotPredVal)

#Get t-test value (min)
temp = ss.ttest_ind(np.log10(PredGenes_pval['Pval']),np.log10(NotPredGenes_pval['Pval']))
print temp


#Get t-test value (mean)
temp = ss.ttest_ind((PredGenes_pval['MeanPval']),(NotPredGenes_pval['MeanPval']))
print temp

#saving to File
#PredGenes_pval.to_csv(FileOut)
