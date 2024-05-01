from rdkit import Chem
import pandas as pd

def novelty(training, generated):
    """Function to calculate novelty of generated molecules, returns novelty (0-1)
    based on https://github.com/eyalmazuz/MolGen/blob/master/MolGen/src/utils/metrics.py"""
    training_set = set(training)
    generated_set = set(generated)
    
    new_molecules = generated_set - training_set
    return len(new_molecules)/len(generated_set)
    
def calculate_novelty(training_path, generated_path):
    """Function to calculate novelty of generated molecules from a csv file"""
    training = pd.read_csv(training_path)
    generated = pd.read_csv(generated_path)
    generated = generated['Valid_Molecules']
    
    novelty_score = novelty(training, generated)
    print(f'training set: {training.head(1)}')
    print(f'generated set: {generated.head(1)}')
    
    return novelty_score

if __name__ == '__main__':
    training_path = 'scripts/data_preprocessing/actives_list.csv'
    generated_path = 'scripts/generated_molecules/valid_molecules.csv'
    novelty_score = calculate_novelty(training_path, generated_path)
    print(f'Novelty of dataset: {novelty_score} or {novelty_score*100}%')