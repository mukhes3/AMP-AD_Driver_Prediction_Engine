import numpy as np 
import mygene


def ConvertToSymb(GenePred):
    #GenePred - 1-d Python List
    mg = mygene.MyGeneInfo()
    temp = mg.querymany(GenePred , scopes='ensembl.gene', fields='symbol', species='human')
    GenePredSyms = []
    for i in range(len(GenePred)):
        GenePredSyms += [str(temp[i]['symbol'])]
        
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