library(tidyverse)
library(readxl)
library(egg)
library(ggrepel)
library(tidytext)

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

# volcano plot ------------------------------------------------------------
plot_volcano <- function(data) {
  ggplot(data)+
    geom_point(
      aes(x = logFC, y = neg_log_pval),
      colour = "grey"
    ) +
    geom_point(data = subset(
      data,
      (logFC > 2 & PValue < 0.05) | (logFC < -2 & PValue < 0.05)), 
      aes(x = logFC, y = neg_log_pval, colour = Direction),
      alpha = 0.6) +
    geom_hline(yintercept = 1.301) +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2) +
    geom_label_repel(
      data = subset(
        data,
        (logFC > 2 & PValue < 0.05) | (logFC < -2 & PValue < 0.05)) %>%
        group_by(insertion) %>%
        slice_max(abs(ori_distance), n = 10), 
      aes(x = logFC, y = neg_log_pval, label = gene_name),
      force = 1,
      max.overlaps = Inf,
      size = 2.5
    ) +
    facet_wrap(~ insertion) +
    ggtitle("Change in Insertions after Co-culture with Faecal Microbiome") +
    labs(y = expression("-log"["10 "]*"p-value"),
         x = expression("log"["2"]*" Fold Change")) +
    theme_article() +
    theme(
      legend.position="none",
      panel.grid.major.y = element_line(colour="black", size=0.125),
      panel.grid.minor.y = element_line(colour="black", size=0.125),
      panel.grid.major.x = element_line(colour="grey", size=0.125),
      strip.text =  element_text(
        size = 11, 
        family = "sans", 
        colour = "black"))+
    scale_colour_manual(values = c("Red", "Blue"))+
    scale_x_continuous(
      limits = c(-20, 20), 
      breaks = seq(-20, 20, by = 5))
}

map(
  list(PA_combined, EC_tagged_combined, EC_tagless_combined),
  ~ plot_volcano(.x)
)