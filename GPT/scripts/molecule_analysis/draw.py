from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import os

def draw_molecule(path, output_path):
    "Function to draw molecules from the csv file of SMILES strings and saves them as images"
    
    molecules = pd.read_csv(path)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    print(f'Number of molecules loaded: {len(molecules)}')
    
    for c, molecule in enumerate(molecules['Valid_Molecules']):
        m = Chem.MolFromSmiles(molecule)
        img = Draw.MolToImage(m)
        img.save(f'{output_path}molecule_{c}.png')
    
    
if __name__ == '__main__':
    path = 'scripts/generated_molecules/valid_molecules.csv'
    output_path = 'scripts/generated_molecules/molecule_images/'
    draw_molecule(path, output_path)