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

# Define x coordinates
x = np.arange(len(['Generated', 'Curated'])*0.7)

# Define bar width
bar_width = 0.4

# Plot bars
plt.bar(x, data, width=bar_width)

# Replace numerical x-axis labels with your labels
plt.xticks(x, ['Generated', 'Curated'])

# Other plot settings
plt.ylim(0, 1.1)
plt.ylabel('External Diversity')
plt.xlabel('Dataset of Molecules', labelpad=10)
plt.title('External Diversity of Generated and Curated Molecules')
plt.show()