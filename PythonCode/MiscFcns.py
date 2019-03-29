import numpy as np 
import mygene
import scipy.stats as ss


def ConvertToEnsb(GenePred):
    #GenePred - 1-d Python List
    mg = mygene.MyGeneInfo()
    temp = mg.querymany(GenePred , scopes='symbol', fields='ensembl.gene', species='human')
    GenePredEnsb = []
    GenePredNew = []
    Ind = []
    for i in range(len(GenePred)):
        if 'ensembl.gene' in temp[i]:
            GenePredEnsb += [str(temp[i]['ensembl.gene'])]
            GenePredNew += [GenePred[i]]
            Ind += [i]
    
    D = {'Symbol':GenePredNew, 'Ensembl':GenePredEnsb, 'Ind':Ind}
    return D
    


def ConvertToSymb(GenePred):
    #GenePred - 1-d Python List
    mg = mygene.MyGeneInfo()
    temp = mg.querymany(GenePred , scopes='ensembl.gene', fields='symbol', species='human')
    GenePredSyms = []
    for i in range(len(GenePred)):
        if 'symbol' in temp[i]:
            GenePredSyms += [str(temp[i]['symbol'])]
        else:
            GenePredSyms += ['Not found']
        
        
    return GenePredSyms

def PredGenesPval(IGAP_genes,GenePred):
    Int = list(set(GeneNames).intersection(GenePred))
    
    G = [] 
    P = []
    
    for i in range(len(Int)):
        In = IGAP_genes['Genes'][IGAP_genes['Genes']==Int[i]].index[0]
        G += [IGAP_genes['Genes'][In]]
        P += [IGAP_genes['Pval'][In]]

        
    PredGenes_pval = pd.DataFrame(data = {'GeneSymb':G,'Pval':P})
    return PredGenes_pval


def GetInGwas(Gwas,GenePred):
    
    Int = list(set(Gwas['Names']).intersection(GenePred))
    
    Mn = []
    M = []
    
    for i in range(len(Int)):
        In = Gwas['Names'][Gwas['Names']==Int[i]].index[0]
        Mn += [Gwas['Min'][In]]
        M += [Gwas['Mean'][In]]
           
    return({'Min':Mn,'Mean':M})
            

def GetTTestResult(lgcl,GeneId,Gwas):
    
    bla = np.argwhere(lgcl)
    bla2 = np.argwhere(np.logical_not(lgcl))
    
    GenePred = list(GeneId[bla[:,0]])
    
    GenePred2 = ConvertToSymb(GenePred)
    GeneNotPred2 = ConvertToSymb(list(GeneId[bla2[:,0]]))
    
    P = GetInGwas(Gwas,GenePred2)
    NP = GetInGwas(Gwas,GeneNotPred2)
    
    #temp = ss.ttest_ind((P['Mean']),(NP['Mean']))
    temp = ss.mannwhitneyu((P['Mean']),(NP['Mean']))
    print temp
    
    #temp = ss.ttest_ind(np.log10(P['Min']),np.log10(NP['Min']))
    temp = ss.mannwhitneyu(np.log10(P['Min']),np.log10(NP['Min']))
    print temp
    

    
    

    
