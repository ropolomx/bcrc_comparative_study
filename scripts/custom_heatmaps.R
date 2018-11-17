# Script for generating custom heatmaps of the BCRC data

# Load packages -----------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(here)
#source('scripts/meg_utility_functions.R')

# Read data ---------------------------------------------------------------

metadata <- read.csv('BCRC_metadata.csv')

# annotations <- read.csv('megares_annotations_v1.01.csv')

amr_group_norm <- read.csv('results_normalized_all/amr_matrices_updated/normalized/AMR_Group_Normalized.csv')

tax_ranks <- c(
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

krakenNorm <-
  tax_ranks %>%
  map(~ read.csv(
    here(
      'results_normalized_all',
      'kraken_cladeReads_matrices_updated',
      'normalized',
      paste0('kraken_', .x, '_Normalized.csv')
    )
  )) %>%
  set_names(nm = tax_ranks)

krakenCladeParsed <- read.csv('aggregated_data_for_analysis/krakenAnalytical_cladeReads.csv') 

krakenCladeParsed <- 
  krakenCladeParsed %>%
  rename(id = Lineage)

krakenCladeParsed <- 
  krakenCladeParsed %>%
  left_join(.,kraken_tax_dt$cladeReads, by = "id")
  

# AMR heatmaps ------------------------------------------------------------

# Tidy the normalized data

amr_group_tidy <-
  tidyr::gather(amr_group_norm, ID, counts, 2:ncol(amr_group_norm))

amr_group_tidy <-
  amr_group_tidy %>%
  mutate(ID = str_replace(ID, "FC_Con_V055", "FC_V055"))

amr_group_merge <- merge(amr_group_tidy, annotations, by = "group")

amr_group_merge <- merge(amr_group_merge, metadata)

amr_group_merge <-
  amr_group_merge[amr_group_merge$Type != "Wetlands", ]

amr_group_merge <-
  amr_group_merge %>%
  group_by(class) %>%
  arrange(class)

names(amr_group_merge) <-
  str_replace(names(amr_group_merge), "group", "Group")
names(amr_group_merge) <-
  str_replace(names(amr_group_merge), "counts", "Normalized_Counts")

names(amr_group_merge) <-
  str_replace(names(amr_group_merge), "class", "Class")

amr_group_merge$Sample_Type <-
  interaction(amr_group_merge$ID, amr_group_merge$Type)

# Plot the tidy data

# Re-order factors of Type column

amr_group_merge$Type <-
    str_replace(amr_group_merge$Type, "\\_", " ")

amr_group_merge$Type <-
  str_replace(amr_group_merge$Type, "Sewage Treatment", "Sewage Influent")

amr_group_merge$Type <-
  factor(
    amr_group_merge$Type,
    levels = c('Fecal Composite', 'Catch Basin', 'Soil', 'Sewage Influent')
  )

# Generating heatmap of subset of top X groups by sum of counts across Types
# 1. Create another dataframe grouping data by AMR Group
# 2. Summarising groups by calculating the sum
# 3. Pick top X groups
# 3. Cross-reference top groups vs. merge database

# amr_group_merge <- 
#   amr_group_merge %>%
#   mutate(., ID = str_replace(ID, "FC_(006|007|008)_(.*)$", "FC_\\2")) %>%
#   mutate(., ID = str_replace(ID, "Soil_N_(.*)_(006|007)$", "Soil_N_\\1")) %>%
#   mutate(ID = str_replace(ID, "FC_Con_V055", "FC_V055"))


amr_group_top <-
  amr_group_merge %>%
  group_by(Group) %>%
  summarise(Count_sum = sum(Normalized_Counts)) %>%
  arrange(-Count_sum) %>%
  slice(1:100)

amr_group_subset <-
  amr_group_merge[amr_group_merge$Group %in% amr_group_top$Group, ]

# amr_group_subset$Type <- str_replace(amr_group_subset$Type, "\\.", " ")

amr_group_subset$Type <-
  factor(
    amr_group_subset$Type,
    levels = c(
      'Fecal Composite', 
      'Catch Basin', 
      'Soil', 
      'Sewage Influent'
      )
  )

amr_group_subset$ID <-
  factor(
    amr_group_subset$ID,
    levels = c(
      "FC_A062",
      "FC_A070",
      "FC_S034",
      "FC_S040",
      "FC_N003",
      "FC_N013",
      "FC_V046",
      "FC_V053",
      "FC_V055",
      "FC_Nat_V042",
      "FC_Nat_V045",
      "FC_Nat_V052",
      "CB_A018",
      "CB_A038",
      "CB_S008",
      "CB_S026",
      "CB_S043",
      "CB_N001",
      "CB_N002",
      "CB_N021",
      "CB_N022",
      "CB_V011",
      "CB_V027",
      "CB_V012",
      "CB_V028",
      "Soil_N_WF_02May14",
      "Soil_N_WF_20May15",
      "Soil_N_EF_02May14",
      "Soil_N_EF_20May15",
      "ST_C039",
      "ST_C053",
      "ST_C062",
      "ST_M029",
      "ST_M050",
      "ST_M060"          
    )
  )

# amr_group_subset$ID <-
#   factor(
#     amr_group_subset$,
#     levels = c('')
#   )


#meg_heatmap <- ggplot(amr80NormMerge, aes(x=ID, y=Group)) +

amr_group_by_class_hm <-
  ggplot(amr_group_subset, aes(x = ID, y = Group)) +
  geom_tile(aes(fill = log2(Normalized_Counts + 1))) +
  facet_grid(
    Class ~ Type,
    scales = 'free',
    switch = 'x',
    drop = TRUE,
    space = 'free_y'
  ) +
  #facet_wrap(~ Type, scales ='free_x', strip.position = 'bottom', nrow = 1) +
  #facet_grid(Class ~ ., scales ='free') +
  theme(
    panel.background = element_rect(fill = "black", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 17, face = "bold"),
    strip.text.y = element_text(size = 10, angle = 0, face = "bold", margin = margin(4.5,0,4.5,0, "cm")),
    #strip.text.y=element_blank(),
    #strip.background= element_rect(fill = amr80NormMerge$Class),
    #axis.text.y=element_text(size=3),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle = 90,hjust=1),
    #axis.text.y=element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = ),
    panel.spacing.x = unit(0.01, "lines"),
    panel.spacing.y = unit(0.06, "lines"),
    plot.title = element_text(size = 24, hjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  xlab('\nSample Matrix Type') +
  scale_fill_gradient(low = "black", high = "cyan") +
  labs(fill = 'Log2 Normalized Count') +
  ggtitle(paste(
    'Normalized AMR Group Counts',
    ' by Class',
    '\n',
    sep = '',
    collapse = ''
  ))

amr_group_by_class_hm

ggsave(
  here('graphs_updated', 'AMR', 'amr_top_group_by_class_hm_updated.png'),
  amr_group_by_class_hm,
  width = 18,
  height = 14.5,
  units = "in",
  dpi = 600
)

# Kraken heatmap ----------------------------------------------------------

krakenClade <- krakenCladeParsed[, -(2:36)]

krakenNormTidy <- 
  krakenNorm %>%
  map(
    ~ gather(.x, ID, counts, 2:36)
  )

krakenNormTidy <-
  krakenNormTidy %>%
  map(
    ~ mutate(.x, ID = str_replace(ID, "FC_Con_V055", "FC_V055"))
  )

krakenNormMerge <- 
  krakenNormTidy %>%
  imap(
    ~ merge(krakenClade, .x, by = .y)
  )

krakenNormMerge <- 
  krakenNormMerge %>%
  map(
    ~ merge(.x, metadata, by = "ID")
  )

krakenNormMerge <-
  krakenNormMerge %>% 
  map(
    ~ arrange(.x, Phylum) %>% 
      filter(., Domain != "Viruses") %>% 
      filter(., Phylum != "NA")
  )

krakenNormMerge <-
  krakenNormMerge %>%
  map(
    ~ rename(.x, Normalized_Counts = counts)
  )

# Re-order factors of Type column

krakenNormMerge <-
  krakenNormMerge %>%
  map(
    ~ mutate(.x, Matrix_Type = str_replace(Matrix_Type, "_", " "))
  )

krakenNormMerge <-
  krakenNormMerge %>%
  map(~ mutate(.x, Matrix_Type = factor(
    Matrix_Type,
    levels = c('Fecal Composite', 'Catch Basin', 'Soil', 'Wastewater')
  )))

krakenNormTop_Phylum <- 
    krakenNormMerge$Phylum %>%
    group_by(Phylum) %>%
    summarise(Count_sum = sum(Normalized_Counts)) %>%
    arrange(-Count_sum) %>%
    filter(Phylum != "NA") %>%
    slice(1:100)

krakenNormTop_Class <- 
    krakenNormMerge$Class %>%
    group_by(Class) %>%
    summarise(Count_sum = sum(Normalized_Counts)) %>%
    arrange(-Count_sum) %>%
    filter(Class != "NA") %>%
    slice(1:100)

krakenNormTop_Order <- 
    krakenNormMerge$Order %>%
    group_by(Order) %>%
    summarise(Count_sum = sum(Normalized_Counts)) %>%
    arrange(-Count_sum) %>%
    filter(Order != "NA") %>%
    slice(1:100)

krakenNormTop_Family <- 
    krakenNormMerge$Family %>%
    group_by(Family) %>%
    summarise(Count_sum = sum(Normalized_Counts)) %>%
    arrange(-Count_sum) %>%
    filter(Family != "NA") %>%
    slice(1:100)

krakenNormTop_Genus <- 
    krakenNormMerge$Genus %>%
    group_by(Genus) %>%
    summarise(Count_sum = sum(Normalized_Counts)) %>%
    arrange(-Count_sum) %>%
    filter(Genus != "NA") %>%
    slice(1:100)

krakenNormTop_Species <- 
    krakenNormMerge$Species %>%
    group_by(Species) %>%
    summarise(Count_sum = sum(Normalized_Counts)) %>%
    arrange(-Count_sum) %>%
    filter(Species != "NA") %>%
    slice(1:100)

krakenNormTop <- list(
  "Phylum" = krakenNormTop_Phylum,
  "Class" = krakenNormTop_Class,
  "Order" = krakenNormTop_Order,
  "Family" = krakenNormTop_Family,
  "Genus" = krakenNormTop_Genus,
  "Species" = krakenNormTop_Species
)

krakenNormSubset_Phylum <- krakenNormMerge$Phylum[krakenNormMerge$Phylum$Phylum %in% krakenNormTop$Phylum$Phylum,]
krakenNormSubset_Class <- krakenNormMerge$Class[krakenNormMerge$Class$Class %in% krakenNormTop$Class$Class,]
krakenNormSubset_Order <- krakenNormMerge$Order[krakenNormMerge$Order$Order %in% krakenNormTop$Order$Order,]
krakenNormSubset_Family <- krakenNormMerge$Family[krakenNormMerge$Family$Family %in% krakenNormTop$Family$Family,]
krakenNormSubset_Genus <- krakenNormMerge$Genus[krakenNormMerge$Genus$Genus %in% krakenNormTop$Genus$Genus,]
krakenNormSubset_Species <- krakenNormMerge$Species[krakenNormMerge$Species$Species %in% krakenNormTop$Species$Species,]

krakenNormSubset <- list(
    "Phylum" = krakenNormSubset_Phylum,
    "Class" = krakenNormSubset_Class,
    "Order" = krakenNormSubset_Order,
    "Family" = krakenNormSubset_Family,
    "Genus" = krakenNormSubset_Genus,
    "Species" = krakenNormSubset_Species
)
  # krakenNormMerge[krakenNormMerge$Family %in% krakenNormTop$Family, ]

krakenNormSubset <-
  krakenNormSubset %>%
  map(~ mutate(.x, Matrix_Type = factor(
    .x$Matrix_Type,
    levels = c('Fecal Composite',
      'Catch Basin',
      'Soil',
      'Wastewater')
  )))

krakenNormSubset <-
  krakenNormSubset %>%
  map(~ mutate(.x, ID = factor(
    .x$ID,
    levels = c(
      "FC_A062",
      "FC_A070",
      "FC_S034",
      "FC_S040",
      "FC_N003",
      "FC_N013",
      "FC_V046",
      "FC_V053",
      "FC_V055",
      "FC_Nat_V042",
      "FC_Nat_V045",
      "FC_Nat_V052",
      "CB_A018",
      "CB_A038",
      "CB_S008",
      "CB_S026",
      "CB_S043",
      "CB_N001",
      "CB_N002",
      "CB_N021",
      "CB_N022",
      "CB_V012",
      "CB_V028",
      "CB_V011",
      "CB_V027",
      "Soil_N_WF_02May14",
      "Soil_N_WF_20May15",
      "Soil_N_EF_02May14",
      "Soil_N_EF_20May15",
      "ST_C039",
      "ST_C053",
      "ST_C062",
      "ST_M029",
      "ST_M050",
      "ST_M060"          
    )
  )))

krakenNormSubset <-
  krakenNormSubset %>%
  map(
    ~ .x %>%
      mutate(.,Phylum = str_replace(Phylum, "^p_(.*)$", "\\1")) %>%
      mutate(.,Class = str_replace(Class, "^c_(.*)$", "\\1")) %>%
      mutate(.,Order = str_replace(Order, "^o_(.*)$", "\\1")) %>%
      mutate(.,Family = str_replace(Family, "^f_(.*)$", "\\1")) %>%
      mutate(.,Genus = str_replace(Genus, "^g_(.*)$", "\\1")) %>%
      mutate(.,Species = str_replace(Species, "^s_(.*)$", "\\1"))
  )

kraken_hm <- function(kraken_subset, tax){
  ggplot(kraken_subset, aes_string(x = 'ID', y = tax)) +
  geom_tile(aes(fill = log2(Normalized_Counts + 1))) +
  facet_grid(
    Phylum ~ Matrix_Type,
    scales = 'free',
    switch = 'x',
    drop = TRUE,
    space = 'free_y'
  ) +
  #facet_wrap(~ Type, scales ='free_x', strip.position = 'bottom', nrow = 1) +
  #facet_grid(Class ~ ., scales ='free') +
  theme(
    panel.background = element_rect(fill = "black", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 17, face = "bold"),
    strip.text.y = element_text(size = 10, angle = 0, face = "bold", margin = margin(4.5,0,4.5,0, "cm")),
    #strip.text.y=element_blank(),
    #strip.background= element_rect(fill = amr80NormMerge$Class),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(size = 10, angle = 90,hjust=1),
    #axis.text.y=element_blank(),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = ),
    panel.spacing.x = unit(0.00, "lines"),
    panel.spacing.y = unit(0.02, "lines"),
    plot.title = element_text(size = 24, hjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  xlab('\nSample Matrix Type') +
  scale_fill_gradient(low = "black", high = "cyan") +
  labs(fill = 'Log2 Normalized Count') +
    ggtitle(paste(
      'Normalized',
      tax,
      'Counts',
      'by',
      'Phylum',
      '\n',
      sep = ' ',
      collapse = ''
    ))
}

kraken_hm_species <- kraken_hm(krakenNormSubset$Species, "Species")

ggsave(
  here(
    'graphs_updated',
    'Microbiome_cladeReads',
    'kraken_top_Species_by_Phylum_hm.png'
  ),
  kraken_hm_species,
  width = 18,
  height = 14.5,
  units = "in",
  dpi = 600
)

all_kraken_hm <- 
  krakenNormSubset %>%
  imap(
    ~ kraken_hm(.x,.y)
  )

all_kraken_hm %>%
  iwalk( ~ ggsave(
    filename = here(
      'graphs_updated',
      'Microbiome_cladeReads',
      paste0('kraken_top_', .y, '_by_Phylum_hm_nolabs.png')
    ),
    plot = .x,
    width = 18,
    height = 14.5,
    units = "in",
    dpi = 600
  ))


# Fancier Kraken heatmaps -------------------------------------------------

meg_heatmap_kraken <- function(df,
  sample_var,
  taxon,
  facet1,
  facet2a,
  facet2b) {
  kraken_heatmap <- ggplot(df, aes_string(x = sample_var, y = taxon)) +
    geom_tile(aes(fill = log2(Normalized_Counts + 1))) +
    facet_grid(
      as.formula(paste0(facet1, "~", facet2a, "+", facet2b)),
      scales = 'free',
      switch = 'x',
      drop = TRUE,
      space = 'free_y'
    ) +
    #facet_wrap(~ Type, scales ='free_x', strip.position = 'bottom', nrow = 1) +
    #facet_grid(Class ~ ., scales ='free') +
    theme(
      panel.background = element_rect(fill = "black", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 10, angle = 0),
      #strip.text.y=element_blank(),
      #strip.background= element_rect(fill = amr80NormMerge$Class),
      #axis.text.y=element_text(size=3),
      axis.text.x = element_blank(),
      #axis.text.y=element_blank(),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 20),
      panel.spacing.x = unit(0.1, "lines"),
      panel.spacing.y = unit(0.02, "lines"),
      plot.title = element_text(size = 24, hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    xlab(paste('Samples by ', 'Type', sep = '', collapse = '')) +
    scale_fill_gradient(low = "black", high = "cyan") +
    labs(fill = 'Log2 Normalized Count') +
    ggtitle(
      paste(
        'Microbiome ',
        taxon,
        ' Normalized Counts by ',
        '\n',
        facet1,
        ', ',
        facet2a,
        ', and ',
        facet2b,
        '\n',
        sep = '',
        collapse = ''
      )
    )
  
  print(kraken_heatmap)
}

meg_barplot_kraken <-
  ggplot(krakenNatConvFCSubset,
    aes(x = NatType, y = Normalized_Counts, fill = Class)) +
  geom_bar(stat = 'identity') +
  #scale_fill_brewer(palette="Spectral") +
  theme(
    strip.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(
      size = 22,
      vjust = 1,
      hjust = 1,
      angle = 33
    ),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 26, hjust = 0.25)
  ) +
  xlab('Class') +
  ylab('Mean of Normalized Count\n') +
  ggtitle(
    paste(
      'Mean ',
      'Microbiome',
      ' ',
      'Class',
      ' Normalized Count by ',
      'NatType',
      ' in FC',
      '\n',
      sep = '',
      collapse = ''
    )
  )
meg_barplot_kraken

meg_heatmap_kraken(
  df = krakenNormSubset,
  sample_var = "ID",
  taxon = "Species",
  facet1 = "Phylum",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Genus",
  facet1 = "Phylum",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Species",
  facet1 = "Family",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Genus",
  facet1 = "Family",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Species",
  facet1 = "Class",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Genus",
  facet1 = "Class",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Species",
  facet1 = "Order",
  facet2a = "Type",
  facet2b = "NatType"
)

meg_heatmap_kraken(
  df = krakenNatConvSubset,
  sample_var = "ID",
  taxon = "Genus",
  facet1 = "Order",
  facet2a = "Type",
  facet2b = "NatType"
)

