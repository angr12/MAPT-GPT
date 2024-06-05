import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('scripts/molecule_analysis/')

from utils.tanimoto import *

generated = pd.read_csv('scripts/generated_molecules/valid_molecules.csv')
training = pd.read_csv('scripts/data_preprocessing/actives_list.csv', header=None, names=['Molecules'], dtype=str)
curated = pd.read_csv('scripts/generated_molecules/curated_molecules.csv')

# calculate internal diversity 
training_div = internal_diversity(training['Molecules'])
generated_div = internal_diversity(generated['Smiles'])
curated_div = internal_diversity(curated['Smiles'])

print(f'Internal Diversity of Training: {training_div}')
print(f'Internal Diversity of Generated: {generated_div}')
print(f'Internal Diversity of Curated: {curated_div}')

# Bar chart
data = [training_div, generated_div, curated_div]
plt.bar(['Training', 'Generated', 'Curated'], data)
plt.ylabel('Internal Diversity')
plt.xlabel('Dataset')
plt.title('Internal Diversity of Training, Generated and Curated Molecules')
plt.show()