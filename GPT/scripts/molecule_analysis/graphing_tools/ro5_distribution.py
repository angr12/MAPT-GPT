import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('scripts/molecule_analysis/')

from utils.RO5 import *

# Load in the data
generated = pd.read_csv('scripts/generated_molecules/valid_molecules.csv')
training = pd.read_csv('scripts/data_preprocessing/actives_list.csv', header=None, names=['Molecules'], dtype=str)
curated = pd.read_csv('scripts/generated_molecules/curated_molecules.csv')

# calculate druglikeness
generated_ro5 = append_RO5(generated['Smiles'])
training_ro5 = append_RO5(training['Molecules'])
curated_ro5 = append_RO5(curated['Smiles'])

# metrics
av_vio_gen = generated_ro5['Num_violations'].mean()
av_vio_train = training_ro5['Num_violations'].mean()
av_vio_cur = curated_ro5['Num_violations'].mean()

print(f'Average number of RO5 violations of Generated: {av_vio_gen}')
print(f'Average number of RO5 violations of Training: {av_vio_train}')
print(f'Average number of RO5 violations of Curated: {av_vio_cur}')
print('\n')

print(f'Max/Min number of RO5 violations of Generated: {generated_ro5["Num_violations"].max()} {generated_ro5["Num_violations"].min()}')
print(f'Max/Min number of RO5 violations of Training: {training_ro5["Num_violations"].max()} {training_ro5["Num_violations"].min()}')
print(f'Max/Min number of RO5 violations of Curated: {curated_ro5["Num_violations"].max()} {curated_ro5["Num_violations"].min()}')
print('\n')

print(f'Percentage of molecules that pass RO5 of Generated: {generated_ro5["Pass"].value_counts()[True]/len(generated_ro5)*100}%')
print(f'Percentage of molecules that pass RO5 of Training: {training_ro5["Pass"].value_counts()[True]/len(training_ro5)*100}%')
print(f'Percentage of molecules that pass RO5 of Curated: {curated_ro5["Pass"].value_counts()[True]/len(curated_ro5)*100}%')

# count number of RO5 violations for the datasets (from 0 to 5)
generated_count = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
for i in generated_ro5['Num_violations']:
    generated_count[i] += 1

curated_count = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
for i in curated_ro5['Num_violations']:
    curated_count[i] += 1

plt.bar(generated_count.keys(), generated_count.values(), align='center', width=0.5, alpha=0.8)
plt.axvline(x=generated_ro5['Num_violations'].mean(), color='r', linestyle=':', label='Mean Number of Violations')
plt.ylabel('Number of Molecules')
plt.xlabel('Number of RO5 Violations')
plt.title('RO5 Violations of Generated Molecules')
plt.savefig('scripts/molecule_analysis/graphs/ro5_generated.png')

plt.clf()
plt.bar(np.arange(0,5), training_ro5['Num_violations'].value_counts().sort_index(), align='center', width=0.5, alpha=0.8)
plt.axvline(x=training_ro5['Num_violations'].mean(), color='r', linestyle=':', label='Mean Number of Violations')
plt.ylabel('Number of Molecules')
plt.xlabel('Number of RO5 Violations')
plt.title('RO5 Violations of Training Molecules')
plt.savefig('scripts/molecule_analysis/graphs/ro5_training.png')

plt.clf()
plt.bar(curated_count.keys(), curated_count.values(),  align='center', width=0.5, alpha=0.8)
plt.axvline(x=curated_ro5['Num_violations'].mean(), color='r', linestyle=':', label='Mean Number of Violations')
plt.ylabel('Number of Molecules')
plt.xlabel('Number of RO5 Violations')
plt.title('RO5 Violations of Curated Molecules')
plt.savefig('scripts/molecule_analysis/graphs/ro5_curated.png')