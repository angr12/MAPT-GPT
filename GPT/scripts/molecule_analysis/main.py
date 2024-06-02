import pandas as pd

generated_mols = pd.read_csv('scripts/generated_molecules/generated.csv')
training_mols = pd.read_csv('scripts/data_preprocessing/actives_list.csv')


# calculate validity
from utils.validity import *
valid_df, validity = filter_valid(generated_mols)
print(f'Validity of dataset: {validity*100}%')

output_df = valid_df.copy()

# calculate druglikeness
from utils.QED import *
QED_df = append_QED(valid_df['Smiles'])

av_qed = QED_df['QED'].mean()
print(f'Average QED of Generated dataset: {av_qed}')
print(f'Max QED of Generated: {QED_df["QED"].max()}')
print(f'Min QED of Generated: {QED_df["QED"].min()}')
print(f'Number of molecules with QED > 0.67: {QED_df[QED_df["QED"] > 0.67].shape[0]}')

# RO5 Calculations
from utils.RO5 import *
RO5_df = append_RO5(valid_df['Smiles'])
av_vio = RO5_df['Num_violations'].mean()
print(f'Average number of RO5 violations: {av_vio}')
print(f'Max number of RO5 violations: {RO5_df["Num_violations"].max()}')
print(f'Number of molecules that pass RO5: {RO5_df["Pass"].value_counts()[True]}')
print(f'Percentage of molecules passing RO5: {RO5_df["Pass"].value_counts()[True]/len(RO5_df)*100}%')

# calculate SAscore
from utils.SA_score import *
SAscore_df = append_SA(valid_df['Smiles'])
av_sa = SAscore_df['SA_Score'].mean()
print(f'Average SA score of dataset: {av_sa}')
print(f'Max SA score: {SAscore_df["SA_Score"].max()}')
print(f'Min SA score: {SAscore_df["SA_Score"].min()}')


#calculate diversity
from utils.tanimoto import *
# ed = external_diversity(valid_df['Smiles'], training_mols['Molecules'])
i_div = internal_diversity(valid_df['Smiles'])

# print(f'External Diversity: {ed}')
print(f'Internal Diversity: {i_div}')

# valid_df.to_csv( 'scripts/generated_molecules/'+ 'valid_molecules.csv', index=False)
