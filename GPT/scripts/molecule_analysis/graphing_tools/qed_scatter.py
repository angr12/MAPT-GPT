import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('scripts/molecule_analysis/')

# Load in the data
generated = pd.read_csv('scripts/generated_molecules/valid_molecules.csv')
training_mols = pd.read_csv('scripts/data_preprocessing/actives_list.csv', header=None, names=['Molecules'], dtype=str)

# calculate druglikeness
from utils.QED import *
generated_qed = append_QED(generated['Smiles'].astype(str))
generated.insert(len(generated.columns), 'QED', generated_qed)

training_qed = append_QED(training_mols['Molecules'])
training_mols.insert(len(training_mols.columns), 'QED', training_qed)

# Filter out molecules with QED > 0.67
curated = generated[generated.QED > 0.67]
# print(curated.head())
# print(curated.shape)

av_qed = generated_qed['QED'].mean()

# box plot
generated_arr = np.asarray(generated_qed['QED'])
training_arr = np.asarray(training_qed['QED'])
curated_arr = np.asarray(curated['QED'])
data = [generated_arr, training_arr, curated_arr]

plt.boxplot(data)
plt.ylabel('QED')
plt.title('QED of Generated vs Training Set of Molecules')
plt.xticks([1, 2, 3], ['Generated', 'Training', 'Generated (Curated)'])
plt.show()