import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# load metrics.csv
path = 'scripts/checkpoints/trained_models/final_model/metrics.csv'
df = pd.read_csv(path)
filtered_df = df[df['ppl_epoch'].notnull()]

# plot ppl_epoch
# print(filtered_df.head())

epoch = np.asarray(filtered_df['epoch'])
ppl_epoch = np.asarray(filtered_df['ppl_epoch'])

min_x = np.argmin(ppl_epoch)
min_y = np.min(ppl_epoch)

plt.plot(epoch, ppl_epoch)
plt.scatter(min_x, min_y,c='r', label='Minimum Perplexity')
plt.text(min_x+1, min_y+1, f'({min_x}, {min_y:.{2}f})', ha='right')
plt.legend()
plt.title('Perplexity per Epoch vs Epochs')
plt.xlabel('Epoch')
plt.ylabel('Perplexity per Epoch')
plt.show()