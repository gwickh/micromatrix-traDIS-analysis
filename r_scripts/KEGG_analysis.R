library(clusterProfiler)
library(AnnotationHub)
library(AnnotationForge)
library(rtracklayer)
library(biomaRt)
library(KEGGREST)
library(pathview)
library(tidyverse)
library(egg)
library(ggbeeswarm)

# cleanup -----------------------------------------------------------------
setwd("/Volumes/Research-Groups/Mark-Webber/Gregory_Wickham/Micromatrix_TraDIS/Pathway_analysis/")

for (k in list("PA_combined", "EC_tagged_combined", "EC_tagless_combined")) {
  assign(
    k,
    read_csv(paste0(k, ".csv")) %>%
      separate(
        locus_tag, 
        c("locus_tag", "insertion"), 
        "__", 
        remove = FALSE) %>%
      separate(
        gene_name, 
        "gene_name", 
        "__", 
        remove = TRUE) %>%  
      mutate(
        insertion = replace_na(insertion, "coding"),
        across(
          all_of("gene_name"), 
          ~ ifelse(str_detect(., "^[0-9]"), NA, .)
        ),
        gene_name = coalesce(gene_name, locus_tag),
        neg_log_pval = -log10(PValue),
        ori_distance = sqrt(logFC^2 + neg_log_pval^2),
        Direction = ifelse(logFC > 0, "Positive", "Negative")
      )
  )
}

EC_tagged_combined <- left_join(EC_tagged_combined, sequence, by = c("gene_name"="gene"), relationship = "many-to-many")
EC_tagless_combined <- left_join(EC_tagless_combined, sequence, by = c("gene_name"="gene"), relationship = "many-to-many")


# KEGG --------------------------------------------------------------------
PA_combined_significant <- PA_combined %>%
  filter(abs(logFC) > 2 & PValue < 0.01) %>%
  pull(1) %>%
  unique() %>%
  enrichKEGG(
    gene = .,
    organism = 'pae',
    pvalueCutoff = 0.05)
PA_combined_pathways <- tibble(PA_combined_significant@result) %>%
  filter(p.adjust < 0.05) %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(
    Description = str_remove(Description, " - Pseudomonas aeruginosa PAO1")) %>%
  right_join(
    filter(PA_combined, abs(logFC) > 2 & PValue < 0.01),
    by = c("geneID" = "locus_tag"),
    relationship = "many-to-many",    
  ) %>%
  drop_na()

EC_tagged_combined_significant <- EC_tagged_combined %>%
  filter(abs(logFC) > 2 & PValue < 0.01) %>%
  pull(-1) %>%
  unique() %>%
  enrichKEGG(
    gene = .,
    organism = 'eco',
    pvalueCutoff = 0.05)
EC_tagged_combined_pathways <- tibble(EC_tagged_combined_significant@result) %>%
  filter(p.adjust < 0.05) %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(
    Description = str_remove(Description, " - Escherichia coli K-12 MG1655")) %>%
  right_join(
    filter(EC_tagged_combined, abs(logFC) > 2 & PValue < 0.01),
    by = c("geneID" = "locus_tag"),
    relationship = "many-to-many",    
  ) %>%
  drop_na()

EC_tagless_combined_significant <- EC_tagless_combined %>%
  filter(abs(logFC) > 2 & PValue < 0.01) %>%
  pull(-1) %>%
  unique() %>%
  enrichKEGG(
    gene = .,
    organism = 'eco',
    pvalueCutoff = 0.05)
EC_tagless_combined_pathways <- tibble(EC_tagless_combined_significant@result) %>%
  filter(p.adjust < 0.05) %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(
    Description = str_remove(Description, " - Escherichia coli K-12 MG1655")) %>%
  right_join(
    filter(EC_tagless_combined, abs(logFC) > 2 & PValue < 0.05),
    by = c("geneID" = "locus_tag"),
    relationship = "many-to-many",    
  ) %>%
  drop_na()



# plot KEGG scatter -------------------------------------------------------
plot_kegg_scatter <- function(data, species) {
  ggplot(
    data,
    aes(y = logFC, x = reorder(Description, -logFC, mean), colour = neg_log_pval)
  )+
    geom_beeswarm(
      size = 2,
      alpha = 0.75,
      cex = 2
    )+
    geom_hline(yintercept = 0, size = 0.5) +
    ggtitle(
      expression("Log"["2"]*" Fold Change of Insertions in Significantly Enriched KEGG Pathways During Co-culture of P. aeruginosa with Faecal Microbiome"))+
    labs(y = expression("-Log"["2"]*" Fold Change"),
         x = "Enriched KEGG Pathway",
         colour = expression("-log"["10 "]*"Adjusted p-value"))+
    theme_article()+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust =1),
      plot.title = element_text(size = 10),
      panel.grid.major.y = element_line(
        colour="black", 
        size=0.125),
      panel.grid.major.x = element_line(
        colour="grey", 
        size=0.125),
      plot.margin = margin(10, 10, 10, 100),
      strip.text =  element_text(
        size = 11, 
        family = "sans", 
        colour = "black"))+
    scale_colour_viridis_c()
}

map(
  list(PA_combined_pathways),
  ~ plot_kegg_scatter(.x, list("P. aeruginosa"))
)