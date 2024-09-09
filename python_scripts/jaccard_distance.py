import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy.spatial.distance import jaccard, pdist, squareform

project_path = "~/webber_group/Gregory_Wickham/Methylation_modification/"
methylation_predictions = pd.read_csv(project_path+"mapped_methylation_sites.csv")

methylation_predictions['presence'] = 1
presence_absence_matrix = methylation_predictions.pivot_table(
    index = 'genome', 
    columns = 'methylated_sequence',
    values = 'presence',
    )\
    .fillna(0)\
    .astype(int)

jaccard_distances = pdist(presence_absence_matrix.values, metric='jaccard')
jaccard_matrix = squareform(jaccard_distances)
jaccard_df = pd.DataFrame(jaccard_matrix, index=presence_absence_matrix.index, columns=presence_absence_matrix.index)

plt.figure(figsize=(10, 8))
sns.heatmap(
    jaccard_df, 
    annot=True, 
    cmap='viridis', 
    cbar=True, 
    linewidths=0.5,
    mask = np.triu(np.ones_like(jaccard_df, dtype=bool))
    )
plt.title("Jaccard Distance")
plt.show()
