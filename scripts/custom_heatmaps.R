# Script for generating custom heatmaps of the BCRC data

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
#source('scripts/meg_utility_functions.R')

metadata <- read.csv('BCRC_metadata.csv')

# annotations <- read.csv('megares_annotations_v1.01.csv')

amr_group_norm <- read_csv('amr_matrices/normalized/AMR_Group_Normalized.csv')

krakenNorm <- read.csv('kraken_matrices/normalized/kraken_Genus_Normalized.csv')

# krakenParsed <- read.csv('kraken_analytic_parsed.csv')

# Tidy the normalized data

amr_group_tidy <- tidyr::gather(amr_group_norm, ID, counts, 2:ncol(amr_group_norm))

amr_group_merge <- merge(amr_group_tidy, annotations, by ="group")

amr_group_merge <- merge(amr_group_merge, metadata)

amr_group_merge <- amr_group_merge[amr_group_merge$Type != "Wetlands",]

amr_group_merge <- 
  amr_group_merge %>% 
  group_by(class) %>% 
  arrange(class)

names(amr_group_merge) <- str_replace(names(amr_group_merge), "group", "Group")
names(amr_group_merge) <- str_replace(names(amr_group_merge), "counts", "Normalized_Counts")

names(amr_group_merge) <- str_replace(names(amr_group_merge), "class", "Class")

amr_group_merge$Sample_Type <- interaction(amr_group_merge$ID, amr_group_merge$Type)

# Plot the tidy data

# Re-order factors of Type column

amr_group_merge$Type <- str_replace(amr_group_merge$Type, "\\_", " ")

amr_group_merge$Type <- factor(amr_group_merge$Type, levels = c('Fecal Composite', 'Catch Basin', 'Soil', 'Sewage Treatment'))

# Generating heatmap of subset of top X groups by sum of counts across Types
# 1. Create another dataframe grouping data by AMR Group
# 2. Summarising groups by calculating the sum
# 3. Pick top X groups
# 3. Cross-reference top groups vs. merge database

amr_group_top <- 
  amr_group_merge %>% 
  group_by(Group) %>%
  summarise(Count_sum=sum(Normalized_Counts)) %>% 
  arrange(-Count_sum) %>%
  slice(1:100)

amr_group_subset <-amr_group_merge[amr_group_merge$Group %in% amr_group_top$Group,]

# amr_group_subset$Type <- str_replace(amr_group_subset$Type, "\\.", " ")

amr_group_subset$Type <- factor(amr_group_subset$Type, levels = c('Fecal Composite', 'Catch Basin', 'Soil', 'Sewage Treatment'))

#meg_heatmap <- ggplot(amr80NormMerge, aes(x=ID, y=Group)) + 
  
amr_group_by_class_hm <- ggplot(amr_group_subset, aes(x=ID, y=Group)) + 
  geom_tile(aes(fill=log2(Normalized_Counts+1))) +
  facet_grid(Class ~ Type, scales='free', switch = 'x', drop=TRUE, space = 'free_y') +
  #facet_wrap(~ Type, scales ='free_x', strip.position = 'bottom', nrow = 1) +
  #facet_grid(Class ~ ., scales ='free') +
  theme(panel.background=element_rect(fill="black", colour="black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.x=element_text(size=15),
              strip.text.y=element_text(size=11, angle = 0),
              #strip.text.y=element_blank(),
              #strip.background= element_rect(fill = amr80NormMerge$Class),
              #axis.text.y=element_text(size=3),
              axis.text.x=element_blank(),
              #axis.text.y=element_blank(),
              axis.text.y=element_text(size=6.5),
              axis.title.x=element_text(size=17),
              axis.title.y=element_blank(),
              legend.position="bottom",
              legend.title=element_text(size=),
              panel.spacing.x = unit(0.1, "lines"), 
              panel.spacing.y = unit(0.06, "lines"),
              plot.title=element_text(size=24, hjust=0.5),
              plot.margin=unit(c(0,0,0,0), "cm")) +
        xlab(paste('Samples by ', 'Type', sep='', collapse='')) +
        scale_fill_gradient(low="black", high="cyan") +
        labs(fill= 'Log2 Normalized Count') + 
        ggtitle(paste('Normalized AMR Group Counts', ' by Class', '\n',
                      sep='', collapse=''))

amr_group_by_class_hm

ggsave(
  here('graphs','AMR', 'amr_top_group_by_class_heatmap.png'),
  amr_group_by_class_hm,
  width = 18,
  height = 14.5,
  units = "in",
  dpi = 320
)

krakenTax <- krakenParsed[,-(2:49)]

krakenNormTidy <- gather(krakenNorm, ID, counts, 2:49)

krakenNormMerge <- merge(krakenTax, krakenNormTidy, by="Genus")

krakenNormMerge <- merge(krakenNormMerge, metadata, by = "ID")

krakenNormMerge <- krakenNormMerge[krakenNormMerge$Type != "Wetlands",]

krakenNormMerge <- krakenNormMerge %>% arrange(Phylum) %>% filter(!Domain %in% "Viruses") %>% filter(!Phylum %in% "") %>% drop_na(Phylum)

names(krakenNormMerge) <- str_replace(names(krakenNormMerge), "counts", "Normalized_Counts")

# Re-order factors of Type column

krakenNormMerge$Type <- str_replace(krakenNormMerge$Type, "\\.", " ")

krakenNormMerge$Type <- factor(krakenNormMerge$Type, levels = c('Fecal Composite', 'Catch Basin', 'Soil', 'Sewage Treatment'))

krakenNormTop <- amr80NormMerge %>% 
  group_by(Class) %>%
  summarise(Count_sum=sum(Normalized_Counts)) %>% 
  arrange(-Count_sum) %>%
  slice(1:100)

amr80NormSubset <- amr80NormMerge[amr80NormMerge$Group %in% amr80NormTop$Group,]



meg_heatmap_kraken <- function(df,
                               sample_var,
                               taxon,
                               facet1,
                               facet2a,
                               facet2b) {
  
  kraken_heatmap <- ggplot(df, aes_string(x=sample_var, y=taxon)) + 
  geom_tile(aes(fill=log2(Normalized_Counts+1))) +
  facet_grid(as.formula(paste0(facet1,"~", facet2a, "+", facet2b)), scales='free', switch = 'x', drop=TRUE, space = 'free_y') +
  #facet_wrap(~ Type, scales ='free_x', strip.position = 'bottom', nrow = 1) +
  #facet_grid(Class ~ ., scales ='free') +
  theme(panel.background=element_rect(fill="black", colour="black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.x=element_text(size=15),
              strip.text.y=element_text(size=10, angle = 0),
              #strip.text.y=element_blank(),
              #strip.background= element_rect(fill = amr80NormMerge$Class),
              #axis.text.y=element_text(size=3),
              axis.text.x=element_blank(),
              #axis.text.y=element_blank(),
              axis.title.x=element_text(size=20),
              axis.title.y=element_blank(),
              legend.position="bottom",
              legend.title=element_text(size=20),
              panel.spacing.x = unit(0.1, "lines"), 
              panel.spacing.y = unit(0.02, "lines"),
              plot.title=element_text(size=24, hjust=0.5),
              plot.margin=unit(c(0,0,0,0), "cm")) +
        xlab(paste('Samples by ', 'Type', sep='', collapse='')) +
        scale_fill_gradient(low="black", high="cyan") +
        labs(fill= 'Log2 Normalized Count') + 
        ggtitle(paste('Microbiome ', taxon, ' Normalized Counts by ', '\n', facet1, ', ', facet2a, ', and ', facet2b, '\n',
                      sep='', collapse=''))

print(kraken_heatmap)
  }

meg_barplot_kraken <- ggplot(krakenNatConvFCSubset, aes(x=NatType, y=Normalized_Counts, fill=Class)) +
        geom_bar(stat='identity') + 
        #scale_fill_brewer(palette="Spectral") +
        theme(strip.text.x=element_text(size=18),
              axis.text.y=element_text(size=20),
              axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
              axis.title.x=element_text(size=24),
              axis.title.y=element_text(size=24),
              legend.title=element_text(size=20),
              legend.text=element_text(size=18),
              plot.title=element_text(size=26, hjust=0.25)) +
        xlab('Class') +
        ylab('Mean of Normalized Count\n') +
        ggtitle(paste('Mean ', 'Microbiome', ' ', 'Class', ' Normalized Count by ', 'NatType', ' in FC', '\n',
                      sep='', collapse=''))
meg_barplot_kraken


meg_heatmap_kraken(krakenNatConvCBSubset, ID, Species, Family, Type, NatType)

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Species", facet1="Phylum", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Genus", facet1="Phylum", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Species", facet1="Family", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Genus", facet1="Family", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Species", facet1="Class", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Genus", facet1="Class", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Species", facet1="Order", facet2a="Type", facet2b="NatType")

meg_heatmap_kraken(df=krakenNatConvSubset, sample_var="ID", taxon="Genus", facet1="Order", facet2a="Type", facet2b="NatType")
