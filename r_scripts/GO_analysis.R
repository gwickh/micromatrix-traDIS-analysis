library(tidyverse)
library(egg)
library(ggbeeswarm)
library(ggrepel)
library(rbioapi)

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
        neg_log_qval = -log10(q.value),
        ori_distance = sqrt(logFC^2 + neg_log_qval^2),
        Direction = ifelse(logFC > 0, "Positive", "Negative")
      )
  )
}

EC_tagged_combined <- left_join(EC_tagged_combined, sequence, by = c("gene_name"="gene"), relationship = "many-to-many")
EC_tagless_combined <- left_join(EC_tagless_combined, sequence, by = c("gene_name"="gene"), relationship = "many-to-many")

# GO --------------------------------------------------------------------
PA_combined_significant <- PA_combined %>%
  filter(abs(logFC) > 2 & q.value < 0.01) %>%
  pull(1) %>%
  unique()
  
PA_enrichment <- rba_panther_enrich(
  PA_combined_significant,
  208964,
  "GO:0008150",
  test_type = "FISHER",
  correction = "FDR",
  cutoff = 0.05
)

PA_combined_pathways <- GO_enriched %>%
  left_join(
    filter(PA_combined, abs(logFC) > 2 & q.value < 0.01),
    by = c("gene_name" = "locus_tag"),
    relationship = "many-to-many",    
  )

# plot GO scatter -------------------------------------------------------
plot_GO_scatter <- function(data, species) {
  ggplot(
    data,
    aes(y = logFC, x = term.label, colour = neg_log_qval)
  )+
    geom_hline(yintercept = 0, size = 0.5) +
    geom_text_repel(
      aes(label = gene_name.y),
      size = 2.5,
      position = position_jitter(seed = 1),
    )+
    geom_point(
      size = 1,
      alpha = 0.5,
      position = position_jitter(seed = 1),
    ) +
    facet_wrap(
      ~ reorder(term.label, -logFC, mean), 
      ncol = 5,
      scales = "free_x"
    ) +
    ggtitle(
      expression("Log"["2"]*" Fold Change of Insertions in Significantly Enriched GO Families During Co-culture of P. aeruginosa with Faecal Microbiome"))+
    labs(y = expression("Log"["2"]*" Fold Change"),
         x = "Enriched GO Family",
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
      panel.spacing = unit(0, "lines"),
      plot.margin = margin(10, 10, 10, 100),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )+
    scale_colour_viridis_c()
}
map(
  list(PA_combined_pathways),
  ~ plot_GO_scatter(.x, list("P. aeruginosa"))
)