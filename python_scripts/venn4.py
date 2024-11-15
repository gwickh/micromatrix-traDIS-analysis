from venny4py.venny4py import *
import pandas as pd

project_path = "~/webber_group/Gregory_Wickham/Micromatrix_TraDIS/TraDIS-Microbiome-240903/reports/"
tradis_output_filtered = pd.read_csv(project_path+"tradis_output_filtered.csv", encoding='ISO-8859-1')

#dict of sets
sets = {
    'Naive Monoculture vs Naive Microbiome': 
        set(
            tradis_output_filtered[
                tradis_output_filtered['name'].isin(['nave_vs_microbiome'])
            ]['gene_name']
        ),
    'Naive Microbiome vs Microbiome with Antibiotics': 
        set(
            tradis_output_filtered[
                tradis_output_filtered['name'].isin(
                    ['microbiome_vs_microbiome_ciprofloxacin',
                     'microbiome_vs_microbiome_azithromycin',
                     'microbiome_vs_microbiome_meropenem']
                )
            ]['gene_name']
        ),
    'Monoculture with Antibiotics vs Microbiome with Antibiotics': 
       set(
            tradis_output_filtered[
                tradis_output_filtered['name'].isin(
                    ['ciprofloxacin_vs_microbiome_ciprofloxacin',
                     'azithromycin_vs_microbiome_azithromycin',
                     'meropenem_vs_microbiome_meropenem']
                )
            ]['gene_name']
        ),
    'Naive Monoculture vs Monoculture with Antibiotics': 
       set(
            tradis_output_filtered[
                tradis_output_filtered['name'].isin(
                    ['naive_vs_ciprofloxacin',
                     'naive_vs_azithromycin',
                     'naive_vs_meropenem']
                )
            ]['gene_name']
        )
}
    
venny4py(sets=sets)
