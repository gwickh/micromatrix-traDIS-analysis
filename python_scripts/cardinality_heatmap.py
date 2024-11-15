import matplotlib.pyplot as plt

# create df with gene names and their corresponding sets
gene_list = []
for set_name, gene_set in sets.items():
    for gene in gene_set:
        gene_list.append((gene, set_name))
df = pd.DataFrame(gene_list, columns=['gene', 'set'])

# Create binary matrix indicating set membership for each gene
df_binary = pd.get_dummies(df['set'])
df_binary['gene'] = df['gene']

# Group by gene and count the membership in each set
df_binary = df_binary.groupby('gene').sum()
gene_counts = df_binary.sum(axis=1)

#Filter genes that appear in >1 set
df_binary_filtered = df_binary[gene_counts > 1].loc[gene_counts[gene_counts > 1].sort_values(ascending=False).index]

#Plot heatmap
fig, ax = plt.subplots(figsize=(21, 3))

cax = ax.imshow(df_binary_filtered.T, aspect='auto', cmap='Blues')

ax.set_xticks(range(len(df_binary_filtered.index)))
ax.set_xticklabels(df_binary_filtered.index, fontsize=8, rotation=90)  # Gene names on x-axis
ax.set_yticks(range(len(df_binary_filtered.columns)))
ax.set_yticklabels(df_binary_filtered.columns, fontsize=8)  # Set names on y-axis
ax.set_xlabel('Genes')
ax.set_title(r'$| \{ A_i : G \in A_i \} | > 1$')

plt.tight_layout()
plt.show()