import pandas as pd

generated_mols = pd.read_csv('scripts/generated_molecules/generated.csv')


# calculate validity
from utils.validity import *
valid_df, validity = filter_valid(generated_mols)
print(f'Validity of dataset: {validity*100}%')

output_df = valid_df.copy()

# calculate druglikeness
from utils.QED import *
QED_df = append_QED(valid_df['Smiles'])

av_qed = QED_df['QED'].mean()
print(f'Average QED of dataset: {av_qed}')
print(f'Max QED: {QED_df["QED"].max()}')
print(f'Min QED: {QED_df["QED"].min()}')

# calculate SAscore
from utils.SA_score import *
SAscore_df = append_SA(valid_df['Smiles'])
av_sa = SAscore_df['SA_Score'].mean()
print(f'Average SA score of dataset: {av_sa}')
print(f'Max SA score: {SAscore_df["SA_Score"].max()}')
print(f'Min SA score: {SAscore_df["SA_Score"].min()}')

# valid_df.to_csv( 'scripts/generated_molecules/'+ 'valid_molecules.csv', index=False)