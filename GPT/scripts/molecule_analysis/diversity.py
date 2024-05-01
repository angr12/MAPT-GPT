import pandas as pd

def calculate_diversity(input_path: str):
    """Function to calculate diversity of generated molecules from a csv file"""
    molecules = pd.read_csv(input_path)
    molecules = molecules['Valid_Molecules']
    diversity = len(set(molecules))/len(molecules)
    
    return diversity 


if __name__ == "__main__":
    molecules_path = 'scripts/generated_molecules/valid_molecules.csv'
    diversity = calculate_diversity(molecules_path)
    print(f'Diversity of dataset: {diversity} or {diversity*100}%')