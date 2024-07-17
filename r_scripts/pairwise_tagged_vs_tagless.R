library(tidyverse)
library(readxl)
library(egg)
library(see)

# cleanup -----------------------------------------------------------------
setwd("/Volumes/Research-Groups/Mark-Webber/Gregory_Wickham/Micromatrix_TraDIS/Pathway_analysis/")

for (k in list("EC_tagged_combined", "EC_tagless_combined")) {
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

EC_pairwise_comparisons <- EC_tagless_combined %>%
  full_join(
    EC_tagged_combined,
    by = c("locus_tag", "insertion", "function")
  ) %>%
  subset(
    select = c("locus_tag", "insertion", "function", "logFC.x", "logFC.y") 
  ) %>%
  rename(c('tagless'='logFC.x', 'tagged'='logFC.y', 'product'='function')) %>%
  mutate(
    delta_log2FC = tagless - tagged 
  )

EC_pairwise_comparisons_long <- EC_pairwise_comparisons %>%
  pivot_longer(
    cols = c('tagless', 'tagged'),
    names_to = 'library',
    values_to = 'log2FC'
  ) %>%
  mutate(
    log2FC = ifelse(is.na(log2FC), -17, log2FC),
    group = paste(locus_tag, product, insertion)
  )

mean_abs_delta_log2FC <- EC_pairwise_comparisons %>%
  subset(
    select = c("insertion", "delta_log2FC")) %>%
  group_by(insertion) %>%
  summarise(
    mean_abs_delta_log2FC = mean(
      abs(delta_log2FC), 
      na.rm = TRUE
    ),
    sd_abs_delta_log2FC = sd(
      abs(delta_log2FC), 
      na.rm = TRUE
    ),
    count_undef = sum(is.na(delta_log2FC))
  )

ANOVA <- aov(
  log2FC ~ library * insertion, 
  data = EC_pairwise_comparisons_long %>%
    subset(
      abs(log2FC) < 15
    )
)
summary(ANOVA)
Tukey <- as.data.frame(TukeyHSD(ANOVA)[3])[c(1,10,15), 4]

mean_abs_delta_log2FC <- merge(mean_abs_delta_log2FC, Tukey, by = "row.names", all = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

# plot tagged-tagless pairwise comparisons -------------------------------------------------------

ggplot(
  EC_pairwise_comparisons_long,
  aes(
    x = library,
    y = log2FC,
    colour = insertion,
    fill = insertion
  )
)+
  geom_line(
    aes(group = group),
    alpha = 0.025
  ) +
  geom_violinhalf(
    data = subset(
      EC_pairwise_comparisons_long,
      abs(log2FC) < 15
    ),
    colour = "black",
    fill = NA,
    flip = 1 
  ) +
  geom_point(
    alpha = 0.1,
    stroke = 0.75,
    shape = 21,
    color = "black"
  ) +
  geom_hline(yintercept = -16) +
  geom_text(
    data = mean_abs_delta_log2FC,
    x = 1.5,
    y = 5,
    size = 3,
    colour = "black",
    label = paste(
      "mean |Δlog2FC|", 
      mean_abs_delta_log2FC$mean_abs_delta_log2FC,
      "±",
      mean_abs_delta_log2FC$sd_abs_delta_log2FC,
      "\np-adj =",
      mean_abs_delta_log2FC$y
    )
  ) +
  geom_text(
    data = mean_abs_delta_log2FC,
    x = 1.5,
    y = -17,
    size = 3,
    colour = "black",
    label = paste(
      "undefined\nn =", 
      mean_abs_delta_log2FC$count_undef
    )
  ) +
  facet_wrap(~ insertion) +
  theme_article() +
  theme(legend.position = "none")
