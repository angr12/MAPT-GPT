import pandas as pd

generated_mols = pd.read_csv('scripts/generated_molecules/generated200_p96_k50.csv')
training_mols = pd.read_csv('scripts/data_preprocessing/actives_list.csv', header=None, names=['Molecules'], dtype=str)

    
# calculate validity
from utils.validity import *
valid_df, validity = filter_valid(generated_mols)
print(f'Validity of dataset: {validity*100}%')

output_df = valid_df.copy()
# valid_df.to_csv('scripts/generated_molecules/valid_molecules.csv', index=False)


# RO5 Calculations
from utils.RO5 import *
RO5_df = append_RO5(valid_df['Smiles'])
av_vio = RO5_df['Num_violations'].mean()
print(f'Average number of RO5 violations: {av_vio}')
print(f'Max number of RO5 violations: {RO5_df["Num_violations"].max()}')
print(f'Number of molecules that pass RO5: {RO5_df["Pass"].value_counts()[True]}')
print(f'Percentage of molecules passing RO5: {RO5_df["Pass"].value_counts()[True]/len(RO5_df)*100}%')

#calculate diversity
print('\n')
from utils.tanimoto import *
from utils.diversity import *
from utils.novelty import *
ed = external_diversity(valid_df['Smiles'], training_mols['Molecules'])
i_div = internal_diversity(valid_df['Smiles'])

print(f'External Tanimoto Diversity: {ed}')
print(f'Internal Tanimoto Diversity: {i_div}')

novelty = novelty(training_mols, valid_df['Smiles'])
print(f'Simple Novelty of dataset: {novelty*100}%')

diversity = calculate_diversity(valid_df['Smiles'])
print(f'Simple Diversity of dataset: {diversity} or {diversity*100}%')

print('\n Training set values')
print(f'Internal Tanimoto Diversity: {internal_diversity(training_mols['Molecules'])}')
# valid_df.to_csv( 'scripts/generated_molecules/'+ 'valid_molecules.csv', index=False)
