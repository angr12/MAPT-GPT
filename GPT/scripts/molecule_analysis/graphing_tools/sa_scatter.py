import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('scripts/molecule_analysis/')
from utils.SA_score import *

# Load in the data
generated = pd.read_csv('scripts/generated_molecules/valid_molecules.csv')
training_mols = pd.read_csv('scripts/data_preprocessing/actives_list.csv', header=None, names=['Molecules'], dtype=str)
curated = pd.read_csv('scripts/generated_molecules/curated_molecules.csv')

# calculate SAscore
generated_sa = append_SA(generated['Smiles'])
generated.insert(len(generated.columns), 'SA_Score', generated_sa)

training_sa = append_SA(training_mols['Molecules'])
training_mols.insert(len(training_mols.columns), 'SA_Score', training_sa)

curated_sa = append_SA(curated['Smiles'])
curated.insert(len(curated.columns), 'SA_Score', curated_sa)
print(f'Avg SA Score of Curated: {curated_sa.mean()}')
print(f'Max SA Score of Curated: {curated_sa.max()}')
print(f'Min SA Score of Curated: {curated_sa.min()}')

print(f'Avg SA Score of Generated: {generated_sa.mean()}')
print(f'Max SA Score of Generated: {generated_sa.max()}')
print(f'Min SA Score of Generated: {generated_sa.min()}')

print(f'Avg SA Score of Training: {training_sa.mean()}')
print(f'Max SA Score of Training: {training_sa.max()}')
print(f'Min SA Score of Training: {training_sa.min()}')

# Box plot
data = [np.asarray(generated_sa['SA_Score']), np.asarray(training_sa['SA_Score']), np.asarray(curated_sa['SA_Score'])]
plt.boxplot(data)
plt.ylabel('SA Score (A.U.)')
plt.xlabel('Dataset of Molecules')
plt.title('SA Score of Generated vs Training Set of Molecules')
plt.xticks([1, 2, 3], ['Generated', 'Training', 'Generated (Curated)'])
plt.ylim(0)
plt.axhline(y=6, color='r', linestyle=':', label='Threshold for Easy Synthesis')
plt.axhline(y=3.2, color='b', linestyle=':', label='Modal Bioactive Molecule SAScore')
plt.legend()
plt.show()