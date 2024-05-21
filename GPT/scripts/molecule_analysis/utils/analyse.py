import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def histogram_from_column(path, column: str):
    """Function to plot histogram of a column from a csv file"""
    molecules = pd.read_csv(path)
    
    col = molecules[column].to_numpy()
    
    # calculate mean and plot
    mean = np.mean(col)
    plt.axvline(mean, color='r', linestyle='dashed', linewidth=1)
    
    # plot histogram
    plt.hist(col)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    plt.title(f'Histogram of {column} values')
    
    plt.show()
    
if __name__ == '__main__':
    path = 'scripts/generated_molecules/valid_molecules.csv'
    histogram_from_column(path, 'QED')