import pandas as pd
from rdkit import Chem

def calculate_diversity(input_df: pd.DataFrame):
    """Function to calculate diversity of generated molecules from a df"""
    molecules = input_df.to_list()
    diversity = len(set(molecules))/len(molecules)
    
    return diversity 


if __name__ == "__main__":
    molecules_path = 'scripts/generated_molecules/valid_molecules.csv'
    molecules = pd.read_csv(molecules_path)
    diversity = calculate_diversity(molecules['Smiles'])
    print(f'Diversity of dataset: {diversity} or {diversity*100}%')