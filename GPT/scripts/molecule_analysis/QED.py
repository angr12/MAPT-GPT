from rdkit import Chem
from rdkit.Chem import QED
import pandas as pd

path = 'scripts/generated_molecules/valid_molecules.csv'
molecules = pd.read_csv(path)
print(f'Number of molecules loaded: {len(molecules)}')

druglikeness = []
for c, molecule in enumerate(molecules['Valid_Molecules']):
    m = Chem.MolFromSmiles(molecule)
    qed = QED.qed(m)
    druglikeness.append(qed)
    # print(f'QED for molecule at index {c}: {qed}')
molecules.insert(1, 'QED', druglikeness)

molecules.to_csv(path)