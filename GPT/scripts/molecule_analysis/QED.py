from rdkit import Chem
from rdkit.Chem import QED
import pandas as pd

def append_QED(file_path):
    """Function to append QED values to the csv file of SMILES strings"""
    molecules = pd.read_csv(file_path)    
    print(f'Number of molecules loaded: {len(molecules)}')
    druglikeness = []
    
    for c, molecule in enumerate(molecules['Valid_Molecules']):
        m = Chem.MolFromSmiles(molecule)
        qed = QED.qed(m) # calculate QED for each smiles string
        druglikeness.append(qed)
        # print(f'QED for molecule at index {c}: {qed}')
        
    molecules.insert(1, 'QED', druglikeness)
    molecules.to_csv(file_path, index=False)
    
if __name__ == '__main__':
    path = 'scripts/generated_molecules/valid_molecules.csv'
    append_QED(path)
    av_qed = pd.read_csv(path)['QED'].mean()
    print(f'Average QED of dataset: {av_qed}')
    
    training = pd.read_csv('scripts/data_preprocessing/actives.csv')
    training_qed = []
    for mol in training['Smiles'].astype(str):
        if Chem.MolFromSmiles(mol) is not None:
            m = Chem.MolFromSmiles(mol)
            qed = QED.qed(m)
            training_qed.append(qed)
    print(f'Average QED of training set: {sum(training_qed)/len(training_qed)}')
