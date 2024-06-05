import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('scripts/molecule_analysis/')

from utils.tanimoto import *

generated = pd.read_csv('scripts/generated_molecules/valid_molecules.csv')
training = pd.read_csv('scripts/data_preprocessing/actives_list.csv', header=None, names=['Molecules'], dtype=str)
curated = pd.read_csv('scripts/generated_molecules/curated_molecules.csv')

# calculate external diversity
generated_ed = external_diversity(generated['Smiles'], training['Molecules'])
curated_ed = external_diversity(curated['Smiles'], training['Molecules'])

print(f'External Diversity of Generated: {generated_ed}')
print(f'External Diversity of Curated: {curated_ed}')

# Bar chart
data = [generated_ed, curated_ed]
plt.bar(['Generated', 'Curated'], data)
plt.ylabel('External Diversity')
plt.xlabel('Dataset')
plt.title('External Diversity of Generated and Curated Molecules')
plt.show()