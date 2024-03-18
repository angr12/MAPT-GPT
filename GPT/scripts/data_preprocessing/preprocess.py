import pandas as pd
import numpy as np

# Load data from chEMBL csv
data = pd.read_csv("data_preprocessing/chEMBL_potency.csv", on_bad_lines='skip') # line 513 is skipped

print(data.head(2))