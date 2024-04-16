from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def append_RO5(path):
    """Function to append RO5 descriptors and number of violations to the csv file of SMILES strings"""
    molecules = pd.read_csv(path)
    print(f'Number of molecules loaded: {len(molecules)}')
    
    # list to store values of RO5 descriptors
    MW_list = list()
    HBA_list = list()
    HBD_list = list()
    LogP_list = list()
    num_v_list = list()
    
    # check for RO5 violations in each molecule
    for molecule in molecules["Valid_Molecules"]:
        m = Chem.MolFromSmiles(molecule)
        MW = Descriptors.MolWt(m)
        HBA = Descriptors.NOCount(m)
        HBD = Descriptors.NHOHCount(m)
        LogP = Descriptors.MolLogP(m)
        
        num_v = [MW > 500, HBA > 10, HBD > 5, LogP > 5]
        
        MW_list.append(MW)
        HBA_list.append(HBA)
        HBD_list.append(HBD)
        LogP_list.append(LogP)
        num_v_list.append(sum(num_v))
        
    # append new rows to df
    molecules = molecules.assign(MW = MW_list, HBA = HBA_list, HBD = HBD_list, LogP = LogP_list, Num_violations = num_v_list)
    molecules.to_csv(path, index=False)
    
if __name__ == "__main__":
    path = 'scripts/generated_molecules/valid_molecules.csv'
    append_RO5(path)