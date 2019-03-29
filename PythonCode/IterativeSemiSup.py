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

FileOut = 'EC2_L1_agg.csv'
Pen = 'l1'

#Reading feature vector
X = pd.read_csv('./amp_ad_agg_feature_set.csv')

#Dropping GeneID column
GeneId = X['GeneID']
X = X.drop(['GeneID'], axis = 1)
X.head()
X = X.fillna(0)

#Reading response vector
Y = pd.read_excel('./ResponseVec_040318.xlsx')
Y = Y.drop(['GeneID'], axis = 1)
Y.head()

#Performing dimensionality reduction using PCA to 2 dimensions
X2 = X.abs()
X_normed = X.values / X2.values.max(axis=0)

#Balancing classes for better classification performance
tmp = SF.BalanceClasses(X_normed + 0.0 ,Y.values + 0.0)
X_train, y_train = tmp['X'],tmp['y']

print 'Crossvalidation'
Cs = np.logspace(-4,4,10)
Cs = Cs.tolist()
print Cs
C_best = SF.GetBestModel(X_train,y_train,pen=Pen,Cs = Cs)
C_best = C_best[0]
print C_best

D = SF.IterativeTraining(X_normed + 0.0, Y.values + 0.0 ,C = C_best,pen=Pen,
                        Iter = 10,Thresh = 0.9)

Y_pred = D['LR'].predict(X_normed)

Y1 = D['LR'].predict(X_normed)

Y1 = D['LR'].predict_proba(X_normed)
Y1 = Y1[:,1]

print(sum(Y1))


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
bla = np.argwhere(Y1>0.5)
GenePred = list(GeneId[bla[:,0]])
GenePredVal = list(Y1[bla[:,0]])
#Converting to Gene Symbols
GenePred2 = MF.ConvertToSymb(GenePred)

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

#Get t-test value (min)
temp = ss.ttest_ind(np.log10(IGAP_genes['Pval']),np.log10(PredGenes_pval['Pval']))
print temp


#Get t-test value (mean)
temp = ss.ttest_ind((IGAP_genes['MeanPval']),(PredGenes_pval['MeanPval']))
print temp

#saving to File
PredGenes_pval.to_csv(FileOut)
