import pandas as pd
import fcd
import numpy as np

# load molecule sets from chEMBL (split into two sets by chEMBL download)
chEMBL1_path = 'scripts/data_preprocessing/chEMBL1.csv'
chEMBL1 = pd.read_csv(chEMBL1_path, sep=';', usecols=['Smiles'])

chEMBL2_path = 'scripts/data_preprocessing/chEMBL2.csv'
chEMBL2 = pd.read_csv(chEMBL2_path, sep=';', usecols=['Smiles'])

frames = [chEMBL1, chEMBL2]
chEMBL = pd.concat(frames)

print(f'Total molecules loaded: {len(chEMBL)}')
# print(chEMBL.head())

# load generated molecules
generated_path = 'scripts/generated_molecules/top_k_10000.csv'
generated = pd.read_csv(generated_path)
print(f'Generated molecules: {len(generated)}')


# preprocess smiles strings to get canonical smiles and filter invalid smiles
chEMBL_sample = np.random.choice(chEMBL['Smiles'], 5000, replace=False)
generated_sample = np.random.choice(generated['Molecules'], 5000, replace=False)


can_chEMBL= [fcd.utils.canonical(s) for s in chEMBL_sample if fcd.utils.canonical(s) is not None]
can_generated= [fcd.utils.canonical(s) for s in generated_sample if fcd.utils.canonical(s) is not None]


print(f'Valid chEMBL molecules: {len(can_chEMBL)}')
print(f'Valid generated molecules: {len(can_generated)}')

print(can_chEMBL[0:2])
print(can_generated[0:2])
# calculate fcd score

fcd_score = fcd.get_fcd(can_chEMBL, can_generated)
print(f'FCD score: {fcd_score}')
