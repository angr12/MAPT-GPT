from rdkit import Chem
import pandas as pd

def novelty(training, generated):
    """Function to calculate novelty of generated molecules, returns novelty (0-1)
    based on https://github.com/eyalmazuz/MolGen/blob/master/MolGen/src/utils/metrics.py"""
    training_set = set(training)
    generated_set = set(generated)
    
    new_molecules = generated_set - training_set
    for i in generated:
        if i in training_set:
            print(f'Repeated Molecule: {i}')
    return len(new_molecules)/len(generated_set)
    

if __name__ == '__main__':
    training_path = 'scripts/data_preprocessing/actives_list.csv'
    generated_path = 'scripts/generated_molecules/valid_molecules.csv'
    training_mols = pd.read_csv(training_path)['Molecules']
    generated_mols = pd.read_csv(generated_path)['Smiles']
    
    print(training_mols.head())
    print(generated_mols.head())
    
    novelty_score = novelty(training_mols.to_list(), generated_mols.to_list())
    print(f'Novelty of dataset: {novelty_score} or {novelty_score*100}%')