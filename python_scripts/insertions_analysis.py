import pandas as pd
import numpy as np

# plotting_functions/__init__.py
from plotting_functions import plot_volcano

project_path = "~/webber_group/Gregory_Wickham/Micromatrix_TraDIS/TraDIS-Microbiome-240903/reports/"
tradis_output = pd.read_csv(project_path+"tradis_output.csv")

#populate a new column 'insertion site' with insertion location
tradis_output['insertion_site'] = tradis_output['gene_name']\
    .apply(lambda x: "5'" if '5prime' in str(x) 
        else ("3'" if '3prime' in str(x) else 'Coding')
    )

#remove all insertion locations in gene name column and '.csv' string from name column
tradis_output['gene_name'] = tradis_output['gene_name'].str.replace(r'__[^ ]*prime', '', regex=True)
tradis_output['name'] = tradis_output['name'].str.replace('.csv', '', regex=False)

#populate a new column based on antibiotic
tradis_output['antibiotic'] = tradis_output['name'].str.extract(r'.*_([^_]+)$')

#add neg_log_qval column
tradis_output['neg_log_qval'] = -np.log10(tradis_output['q.value'])


##plot volcano plot
#add column to determine conditional colouring (red = down, blue = up, grey = no change)
tradis_output['color_condition'] = \
    np.where(
        (tradis_output['logFC'] > 1) & (tradis_output['q.value'] < 0.05), 
        'Upregulated', 
    np.where(
        (tradis_output['logFC'] < -1) & (tradis_output['q.value'] < 0.05), 
        'Downregulated', 
    'Neutral'))

#subset based on comparison
naive_vs_abx = tradis_output[tradis_output['name'].isin(
    ['naive_vs_ciprofloxacin','naive_vs_azithromycin','naive_vs_meropenem']
    )
]
abx_vs_microbiome_abx = tradis_output[tradis_output['name'].isin(
    ['ciprofloxacin_vs_microbiome_ciprofloxacin',
    'azithromycin_vs_microbiome_azithromycin',
    'meropenem_vs_microbiome_meropenem']
    )
]
microbiome_vs_microbiome_abx = tradis_output[tradis_output['name'].isin(
    ['microbiome_vs_microbiome_ciprofloxacin',
    'microbiome_vs_microbiome_azithromycin',
    'microbiome_vs_microbiome_meropenem']
    )
]

#loop through dataframes
datasets = [naive_vs_abx, abx_vs_microbiome_abx, microbiome_vs_microbiome_abx]
titles = [
    "Change in Insertions Between Monoculture Control and Monoculture with Antibiotic",
    "Change in Insertions Between Monoculture with Antibiotic and Microbiome with Antibiotic",
    "Change in Insertions Between Microbiome Control and Microbiome with Antibiotic"
]

for df,title in zip(datasets, titles):
    plot_volcano(
        data=df,
        col_name='antibiotic',
        x_var='logFC',
        y_var='neg_log_qval',
        hue_var='color_condition',
        style_var='insertion_site',
        label_var='gene_name',
        title=title,
        top_n=10,
        min_dist=1,
        min_inter_annotation_dist=1
    )