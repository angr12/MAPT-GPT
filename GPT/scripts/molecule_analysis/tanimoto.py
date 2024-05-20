from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import pandas as pd
import numpy as np

def external_dviersity(smiles1, smiles2):
    """
    Function to calculate the average Tanimoto similarity between two sets of molecules
    Based on https://arxiv.org/pdf/1708.08227
    """
    
    mols_1 = [Chem.MolFromSmiles(mol) for mol in smiles1 if Chem.MolFromSmiles(mol) is not None] # convert SMILES strings to RDKit molecules
    mols_2 = [Chem.MolFromSmiles(mol) for mol in smiles2 if Chem.MolFromSmiles(mol) is not None] 
    
    fp1 = [FingerprintMols.FingerprintMol(mol) for mol in mols_1] # generate fingerprints for each molecule
    fp2 = [FingerprintMols.FingerprintMol(mol) for mol in mols_2]
    
    T = []
    A = len(fp1)*len(fp2)
    
    for i in fp1:
        for j in fp2:
            T.append(DataStructs.TanimotoSimilarity(i, j))
    
    diversity = np.mean(T)
    return diversity

def internal_diversity(smiles1):
    """
    Function to calculate the average Tanimoto similarity within a set of molecules
    Based on https://arxiv.org/pdf/1708.08227
    """
    A = len(smiles1)
    
    mols = [Chem.MolFromSmiles(mol) for mol in smiles1 if Chem.MolFromSmiles(mol) is not None]
    fp = [FingerprintMols.FingerprintMol(mol, minPath=1, maxPath=7, fpSize=2048,
                               bitsPerHash=2, useHs=True, tgtDensity=0.0,
                               minSize=128) for mol in mols] # args to handle mismatch in fpSize
    
    T = []
    
    for i in range(len(fp)-1):
        s = DataStructs.BulkTanimotoSimilarity(fp[i], fp[i+1:])
        T.extend(s)

    return(np.mean(T))

    
if __name__ == "__main__":
    generated = pd.read_csv('scripts/generated_molecules/top_k.csv')    
    training_set = pd.read_csv('scripts/data_preprocessing/actives.csv')
    print(training_set.head())
    
    ed = external_dviersity(generated['Molecules'], training_set['Smiles'].astype(str))
    print(f'External Diversity: {ed}')
    
    i_div = internal_diversity(generated['Molecules'])
    print(f'Internal Diversity: {i_div}')