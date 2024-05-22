from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def append_RO5(input_df: pd.DataFrame):
    """Function to append RO5 descriptors and number of violations to the df of SMILES strings"""
    print(f'Number of molecules loaded: {len(input_df)}')
    
    # list to store values of RO5 descriptors
    MW_list = list()
    HBA_list = list()
    HBD_list = list()
    LogP_list = list()
    num_v_list = list()
    pass_list = list()
    
    # check for RO5 violations in each molecule
    for molecule in input_df:
        m = Chem.MolFromSmiles(molecule)
        MW = Descriptors.MolWt(m)
        HBA = Descriptors.NOCount(m)
        HBD = Descriptors.NHOHCount(m)
        LogP = Descriptors.MolLogP(m)
        
        num_v = [MW > 500, HBA > 10, HBD > 5, LogP > 5]
        if sum(num_v) >= 2:
            pass_list.append(False)
        else:
            pass_list.append(True)
        
        MW_list.append(MW)
        HBA_list.append(HBA)
        HBD_list.append(HBD)
        LogP_list.append(LogP)
        num_v_list.append(sum(num_v))
        
    # append new rows to df
    d = {'MW': MW_list, 'HBA': HBA_list, 'HBD': HBD_list, 'LogP': LogP_list, 'Num_violations': num_v_list, 'Pass': pass_list}
    RO5 = pd.DataFrame(d)
    return RO5
    
if __name__ == "__main__":
    path = 'scripts/generated_molecules/valid_molecules.csv'
    molecules_df = pd.read_csv(path)
    RO5_df = append_RO5(molecules_df['Smiles'])
    print(RO5_df.head())
    
    # add to the csv
    molecules_df = pd.concat([molecules_df, RO5_df], axis=1)
    # print(molecules_df.head())
    # print(molecules_df.shape)
    
    molecules_df.to_csv(path, index=False)