from rdkit import Chem
import pandas as pd

def filter_valid(input_path, output_path, output_name='valid_molecules.csv'):
    """Function to filter out valid molecules from the generated molecules from path to csv file"""
    molecules = pd.read_csv(input_path)
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
    valid_molecules.to_csv(output_path + output_name, index=False)

if __name__ == '__main__':
    path = 'scripts/generated_molecules/top_k.csv'
    folder_path = 'scripts/generated_molecules/'
    filter_valid(path, folder_path)