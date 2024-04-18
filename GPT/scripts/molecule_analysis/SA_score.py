from rdkit import Chem
from rdkit.Chem import RDConfig
import pandas as pd
import os 
import sys

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))  # add the path to SA_Score
import sascorer

def append_SA(path):
    """Function to append SA score to the csv file of SMILES strings"""
    molecules = pd.read_csv(path)
    
    print(f'Number of molecules loaded: {len(molecules)}')
    synthetic_accessibility = []
    
    for mol in molecules['Valid_Molecules']:
        m = Chem.MolFromSmiles(mol)
        sa = sascorer.calculateScore(m)
        synthetic_accessibility.append(sa)
        
    molecules = molecules.assign(SA_Score = synthetic_accessibility)
    molecules.to_csv(path, index=False)
    
if __name__ == "__main__":
    path = 'scripts/generated_molecules/valid_molecules.csv'
    append_SA(path)