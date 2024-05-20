from rdkit import Chem
import pandas as pd

def filter_valid(input_df: pd.DataFrame):
    """Function to filter out valid molecules from the generated molecules from dataframe
    Returns: df of valid molecules, validity out of 1
    """
    valid_molecules = []
    invalid_molecules = []

    for c, mol in enumerate(input_df['Molecules']):
        m = Chem.MolFromSmiles(mol)
        
        if m is None:
            invalid_molecules.append(mol)
            print(f'Invalid molecule at index {c}: {mol}')
        else:
            valid_molecules.append(mol)
                
    print(f'Valid molecules: {len(valid_molecules)}')
    print(f'Invalid molecules: {len(invalid_molecules)}')

    valid_molecules = pd.DataFrame(valid_molecules, columns=['Smiles'])
    
    return valid_molecules, len(valid_molecules)/len(input_df)

if __name__ == '__main__':
    input_path = 'scripts/generated_molecules/top_k.csv'
    output_path = 'scripts/generated_molecules/'
    input_df = pd.read_csv(input_path)
    
    valid_df, validity = filter_valid(input_df)
    print(f'Validity of dataset: {validity}')
    valid_df.to_csv(output_path + 'valid_molecules.csv', index=False)