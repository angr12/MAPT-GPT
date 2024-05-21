from rdkit import Chem
from rdkit.Chem import RDConfig
import pandas as pd
import os 
import sys

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))  # add the path to SA_Score
import sascorer

def append_SA(input_df: pd.DataFrame):
    """Function to append SA score for a df of SMILES strings"""
    
    print(f'Number of molecules loaded: {len(input_df)}')
    synthetic_accessibility = []
    
    for mol in input_df:
        m = Chem.MolFromSmiles(mol)
        sa = sascorer.calculateScore(m)
        synthetic_accessibility.append(sa)
        
    synthetic_accessibility = pd.DataFrame(synthetic_accessibility, columns=['SA_Score'])
    return synthetic_accessibility
    
if __name__ == "__main__":
    path = 'scripts/generated_molecules/valid_molecules.csv'
    molecules = pd.read_csv(path)
    sa_df = append_SA(molecules['Smiles'])
    print(sa_df.head())
    print(f'Average SA score of dataset: {sa_df["SA_Score"].mean()}')
    
    # add to the csv
    molecules.insert(len(molecules.columns), 'SA_Score', sa_df)
    print(molecules.head())
    molecules.to_csv(path, index=False)