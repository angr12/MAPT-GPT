import sys
print(f'Initial sys.path: {sys.path}')
# sys.path.append('/rds/general/user/rla120/home/FYP/GPT/scripts/smiles-gpt/smiles_gpt')
sys.path.append('/Users/angr1/FYP/GPT/scripts/smiles-gpt/smiles_gpt')
# sys.path.insert(0, '/rds/general/user/rla120/home/FYP/GPT/scripts/smiles-gpt/smiles_gpt')

print(f'Updated sys.path: {sys.path}')
import smiles_gpt as gpt
import pandas as pd