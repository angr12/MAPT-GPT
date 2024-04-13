from rdkit import Chem
import pandas as pd

path = 'scripts/generated_molecules/top_k.csv'
folder_path = 'scripts/generated_molecules/'
molecules = pd.read_csv(path)
valid_molecules = []
invalid_molecules = []

for c, mol in enumerate(molecules['Molecules']):
    m = Chem.MolFromSmiles(mol)
    
    if m is None:
        invalid_molecules.append(mol)
        print(f'Invalid molecule at index {c}: {mol}')
    else:
        valid_molecules.append(mol)
                
print(f'Valid molecules: {len(valid_molecules)}')
print(f'Invalid molecules: {len(invalid_molecules)}')

valid_molecules = pd.DataFrame(valid_molecules, columns=['Valid_Molecules'])
valid_molecules.to_csv(folder_path + 'valid_molecules.csv', index=False)
