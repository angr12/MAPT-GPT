from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import os

def draw_molecule(input_df: pd.DataFrame, output_path):
    "Function to draw molecules from the csv file of SMILES strings and saves them as images"
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    print(f'Number of molecules loaded: {len(input_df)}')
    
    for c, molecule in enumerate(input_df):
        m = Chem.MolFromSmiles(molecule)
        img = Draw.MolToImage(m)
        img.save(f'{output_path}molecule_{c}.png')
    
    
if __name__ == '__main__':
    pass