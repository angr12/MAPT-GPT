from rdkit import Chem
from rdkit.Chem import QED
import pandas as pd

def append_QED(input_df: pd.DataFrame):
    """Function to append QED values to the df of SMILES strings"""
    print(f'Number of molecules loaded: {len(input_df)}')
    druglikeness = []
    
    for c, molecule in enumerate(input_df):
        if Chem.MolFromSmiles(molecule) is not None:
            m = Chem.MolFromSmiles(molecule)
            qed = QED.qed(m) # calculate QED for each smiles string
            druglikeness.append(qed)
            # print(f'QED for molecule at index {c}: {qed}')
        
    druglikeness = pd.DataFrame(druglikeness, columns=['QED'])
    return druglikeness
    
if __name__ == '__main__':
    input_path = 'scripts/generated_molecules/valid_molecules.csv'
    molecules = pd.read_csv(input_path)
    qed_df = append_QED(molecules['Smiles'])
    # print(qed_df.head())
    
    av_qed = qed_df['QED'].mean()
    print(f'Average QED of dataset: {av_qed}')
    molecules.insert(len(molecules.columns), 'QED', qed_df)
    # print(molecules.head())
    
    molecules.to_csv(input_path, index=False)
    
    training = pd.read_csv('scripts/data_preprocessing/actives_list.csv')
    training_qed = append_QED(training['Molecules'].astype(str))
    # print(training_qed.head())
    av_qed = training_qed['QED'].mean()
    print(f'Average QED of training dataset: {av_qed}')