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
print(f'Max QED: {QED_df["QED"].max()}')0

# valid_df.to_csv( 'scripts/generated_molecules/'+ 'valid_molecules.csv', index=False)
