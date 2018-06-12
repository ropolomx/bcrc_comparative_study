## R script template for the production of basic qualitative
## statistics on the AMR and kraken output from AMR++

## Analysis of resistome and microbiome
## Analysis of microbiome with two approaches:
## taxReads: using taxon-specific read counts
## cladeReads: using clade read counts

## Original author: Steven Lakin (Colorado State University)
## Modified and adapted for BCRC comparative analysis by: Rodrigo Ortega Polo (University of Lethbridge)

## The files you want to use for input to this (for the MEG group analyses)
## are the AMR_analytic_matrix.csv and kraken_analytic_matrix.csv.  The AMR
## matrix is identical to the Gene.csv matrix, however the kraken analytic
## matrix is not due to the way that reads get classified using each of
## these methods.

## So you should have pulled these files from the output of the Galaxy pipeline
## and you are now performing this analysis on your local machine. We will assume
## that you've set your working directory to where these files are located on your
## local machine and that you have installed the metagenomeSeq package.

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.

# Loading packages --------------------------------------------------------

all_packages <- c(
  "tidyverse",
  "data.table",
  "vegan",
  "metagenomeSeq",
  "ggplot2",
  "here"
  # "drake"
)

lapply(all_packages, require, character.only = TRUE)

# User Controls -----------------------------------------------------------

## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# Set your working directory to the main folder for analysis:
# Already set in project
# setwd('/home/lakinsm/Documents/morleyBioinformatics/CanadaAnalyticData/20March2017Analysis/')


# Set the output directory for graphs:
graph_output_dir = 'graphs'


# Set the output directory for statistics:
stats_output_dir = 'stats'

# Where is the metadata file stored on your machine?
metadata_filepath = here('BCRC_metadata.csv')

# Name of the megares annotation file used for this project

megares_annotation_filename = here('aggregated_data_for_analysis', 'megaresMegabioUpdated.csv')


# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'

# The following is a list of analyses based on variables in 
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
exploratory_analyses = list(
  # Analysis 1
  # Description: Type comparison for all locations
  list(
    name = 'TypeOverall',
    subsets = list('Type != Wetlands', 'NatType != Natural'),
    exploratory_var = 'Type'
  ),
  
  # Analysis 2
  # Description: Location comparison for all types
  list(
    name = 'LocationOverall',
    subsets = list('Type != Wetlands', 'NatType != Natural'),
    exploratory_var = 'Location'
  ),
  
  # Analysis 3
  # Description: Location comparison within Fecal Composite type
  list(
    name = 'LocationFC',
    subsets = list('Type == Fecal_Composite', 'NatType != Natural'),
    exploratory_var = 'Location'
  ),
  
  # # Analysis 4
  # # Description: Location comparison within Catch Basin type
  # list(
  #     name = 'LocationCB',
  #     subsets = list('Type == Catch_Basin', 'NatType != Natural'),
  #     exploratory_var = 'Location'
  # ),

  # Analysis 5
  # Description: Location comparison within Waste Water sewage treatment type
  list(
    name = 'LocationST',
    subsets = list('Type == Sewage_Treatment'),
    exploratory_var = 'Location'
  ),
  
  # Analysis 6
  # Description: Natural vs conventional Fecal Composite
  list(
    name = 'NaturalConventionalFC',
    subsets = list('Type == Fecal_Composite', 'Location == Vegreville', 'NatType != None'),
    exploratory_var = 'NatType'
  ),
  
  # Analysis 7
  # Description: Natural vs conventional Catch Basin
  list(
    name = 'NaturalConventionalCB',
    subsets = list('Type == Catch_Basin', 'Location == Vegreville'),
    exploratory_var = 'NatType'
  ),
  
  # Analysis 8
  # Description: FieldType comparison within Soil type
  list(
    name = 'SoilFieldType',
    subsets = list('Type == Soil'),
    exploratory_var = 'FieldType'
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.
statistical_analyses = list(
  # Analysis 1
  # Description: Fixed effect for type, control for location using random effect
  list(
    name = 'TypeFixedLocationRandom',
    subsets = list('Type != Wetlands', 'NatType != Natural'),
    model_matrix = '~ 0 + Type',
    contrasts = list('TypeFecal_Composite - TypeCatch_Basin',
                     'TypeFecal_Composite - TypeSewage_Treatment',
                     'TypeCatch_Basin - TypeSewage_Treatment'),
    random_effect = 'Location'
  ),
  
  # Analysis 2
  # Description: Fixed effect for location, control for type using fixed effect
  list(
    name = 'LocationFixedTypeFixed',
    subsets = list('Type != Wetlands', 'NatType != Natural'),
    model_matrix = '~ 0 + Location + Type',
    contrasts = list('LocationAcme - LocationCalgary',
                     'LocationAcme - LocationIron_Springs',
                     'LocationAcme - LocationMedicine_Hat',
                     'LocationAcme - LocationNanton',
                     'LocationAcme - LocationVegreville',
                     'LocationCalgary - LocationIron_Springs',
                     'LocationCalgary - LocationMedicine_Hat',
                     'LocationCalgary - LocationNanton',
                     'LocationCalgary - LocationVegreville',
                     'LocationIron_Springs - LocationMedicine_Hat',
                     'LocationIron_Springs - LocationNanton',
                     'LocationIron_Springs - LocationVegreville',
                     'LocationMedicine_Hat - LocationNanton',
                     'LocationMedicine_Hat - LocationVegreville',
                     'LocationNanton - LocationVegreville'),
    random_effect = NA
  ),
  
  # Analysis 3
  # Description: Natural vs Conventional fixed effect within Vegreville, Fecal Composite
  list(
    name = 'NaturalConventionalFCVegreville',
    subsets = list('Type != Wetlands',
                   'Type != Sewage_Treatment',
                   'Type != Catch_Basin',
                   'NatType != None',
                   'Location == Vegreville'),
    model_matrix = '~ 0 + NatType',
    contrasts = list('NatTypeNatural - NatTypeConventional'),
    random_effect = NA
  )
)


# Automated Code ----------------------------------------------------------

## Modify this as necessary, though you shouldn't need to for basic use.

# Source the utility functions file, which should be in the scripts folder with this file
source(here('scripts','meg_utility_functions.R'))

set.seed(154)  # Seed the RNG, necessary for reproducibility


# We usually filter out genes with wild-type potential.  If you want to include these
# in your analysis, comment this vector out
snp_regex = c('ACRR',
              'CATB',
              'CLS',
              'DFRC',
              'DHFR',
              'DHFRIII',
              'DHFRIX',
              'EMBA',
              'embB',
              'EMBB',
              'EMBC',
              'EMBR',
              'ETHA',
              'FOLP',
              'GIDB',
              'GYRA',
              'gyrB',
              'GYRB',
              'INHA',
              'INIA',
              'INIC',
              'KASA',
              'LIAFSR',
              'LMRA',
              'MARR',
              'MEXR',
              'MEXZ',
              'mprF',
              'MPRF',
              'NDH',
              'omp36',
              'OMP36',
              'OMPF',
              'OPRD',
              'PARC',
              'parE',
              'PARE',
              'PGSA',
              'phoP',
              'PHOP',
              'PNCA',
              'POR',
              'PORB',
              'RAMR',
              'rpoB',
              'RPOB',
              'RPOC',
              'RPSL',
              'SOXS',
              'tetR',
              'TETR',
              'TLYA',
              'TUFAB'
  )


# Import & format Data ----------------------------------------------------

## These files should be standard for all analyses, as they are
## the output matrices from AMR++ nextflow. Additionally,
## you will need to obtain the most recent megares annotations file
## from megares.meglab.org

# If subdirs for stats and exploratory variables don't exist, create them
ifelse(!dir.exists(file.path(graph_output_dir)), dir.create(file.path(graph_output_dir), mode='777'), FALSE)
ifelse(!dir.exists(file.path(stats_output_dir)), dir.create(file.path(stats_output_dir), mode='777'), FALSE)

for( dtype in c('AMR', 'Microbiome_taxonReads', 'Microbiome_cladeReads') ) {
  ifelse(!dir.exists(file.path(graph_output_dir, dtype)),
         dir.create(file.path(graph_output_dir, dtype), mode='777'), FALSE)
  
  for( v in 1:length(exploratory_analyses) ) {
    ifelse(!dir.exists(file.path(graph_output_dir, dtype, exploratory_analyses[[v]]$name)),
           dir.create(file.path(graph_output_dir, dtype, exploratory_analyses[[v]]$name), mode='777'), FALSE)
  }
  
  ifelse(!dir.exists(file.path(stats_output_dir, dtype)),
         dir.create(file.path(stats_output_dir, dtype), mode='777'), FALSE)
  
  for( a in 1:length(statistical_analyses) ) {
    ifelse(!dir.exists(file.path(stats_output_dir, dtype, statistical_analyses[[a]]$name)),
           dir.create(file.path(stats_output_dir, dtype, statistical_analyses[[a]]$name), mode='777'), FALSE)
  }
}

ifelse(!dir.exists(file.path('amr_matrices')), dir.create(file.path('amr_matrices'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_matrices')), dir.create(file.path('kraken_matrices'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_taxonReads_matrices')), dir.create(file.path('kraken_taxonReads_matrices'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_cladeReads_matrices')), dir.create(file.path('kraken_cladeReads_matrices'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('amr_matrices/sparse_normalized')), dir.create(file.path('amr_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/normalized')), dir.create(file.path('amr_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/raw')), dir.create(file.path('amr_matrices/raw'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('kraken_taxonReads_matrices/sparse_normalized')), dir.create(file.path('kraken_taxonReads_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_taxonReads_matrices/normalized')), dir.create(file.path('kraken_taxonReads_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_taxonReads_matrices/raw')), dir.create(file.path('kraken_taxonReads_matrices/raw'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('kraken_cladeReads_matrices/sparse_normalized')), dir.create(file.path('kraken_cladeReads_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_cladeReads_matrices/normalized')), dir.create(file.path('kraken_cladeReads_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_cladeReads_matrices/raw')), dir.create(file.path('kraken_cladeReads_matrices/raw'), mode='777'), FALSE)

# Load the Kraken data, MEGARes annotations, and metadata

# Load post-PhiX filtered Kraken reports

kraken_analytical <- Sys.glob(here("aggregated_data_for_analysis", "krakenAnalytical_*.csv"))

kraken_names <- map_chr(
  kraken_analytical,
  ~ str_replace(.x, "^.*_(.*)\\.csv", "\\1")
)

temp_kraken_list <- map(
  kraken_analytical,
  ~ read.table(.x, header = T, row.names = 1, sep = ",")
) %>%
  set_names(nm = kraken_names)

kraken_new_mr <-
  temp_kraken_list %>%
  map( ~ newMRexperiment(.x[rowSums(.x) > 0, ]))

amr <- newMRexperiment(read.table(here('aggregated_data_for_analysis', 'amrBioAnalytical.csv'), 
                                  header=T, row.names=1, sep=','))

annotations <- data.table(read.csv(megares_annotation_filename, header=T))
annotations$class <- str_replace(annotations$class, "betalactams", "Betalactams")
  
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys

metadata <- read.csv(metadata_filepath, header=T)
metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])

# Normalizing unsplit Kraken and AMR

kraken_css <- 
  kraken_new_mr %>%
  map(~ cumNorm(.x))

cumNorm(amr)

# Normalize clade reads kraken

# Extract the normalized counts into data tables for aggregation

kraken_norm <- 
  kraken_css %>%
  map(~ data.table(MRcounts(.x, norm=T)))
  
kraken_raw <- 
  kraken_css %>%
  map(~ data.table(MRcounts(.x, norm=F)))

amr_norm <- data.table(MRcounts(amr, norm=T))
amr_raw <- data.table(MRcounts(amr, norm=F))

# Aggregate the normalized counts for AMR using the annotations data table, SQL
# outer join, and aggregation with vectorized lapply

amr_norm[, header :=( rownames(amr) ), ]
setkey(amr_norm, header)
amr_norm <- annotations[amr_norm]  # left outer join

amr_raw[, header :=( rownames(amr) ), ]
setkey(amr_raw, header)
amr_raw <- annotations[amr_raw]  # left outer join

# Remove groups that correspond to potentially wild-type genes
amr_raw <- amr_raw[!(group %in% snp_regex), ]
amr_norm<- amr_norm[!(group %in% snp_regex), ]


# Group the AMR data by level for analysis
amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class'])
rownames(amr_class_analytic) <- amr_class$class

amr_class_raw <- amr_raw[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_raw_analytic <- newMRexperiment(counts=amr_class_raw[, .SD, .SDcols=!'class'])
rownames(amr_class_raw_analytic) <- amr_class_raw$class

amr_mech <- amr_norm[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_analytic <- newMRexperiment(counts=amr_mech[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_analytic) <- amr_mech$mechanism

amr_mech_raw <- amr_raw[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_raw_analytic <- newMRexperiment(counts=amr_mech_raw[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_raw_analytic) <- amr_mech_raw$mechanism

amr_group <- amr_norm[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_analytic <- newMRexperiment(counts=amr_group[, .SD, .SDcols=!'group'])
rownames(amr_group_analytic) <- amr_group$group

amr_group_raw <- amr_raw[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_raw_analytic <- newMRexperiment(counts=amr_group_raw[, .SD, .SDcols=!'group'])
rownames(amr_group_raw_analytic) <- amr_group_raw$group

amr_gene_analytic <- newMRexperiment(
  counts=amr_norm[!(group %in% snp_regex),
                  .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])
amr_gene_raw_analytic <- newMRexperiment(
  counts=amr_raw[!(group %in% snp_regex),
                 .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])

rownames(amr_gene_analytic) <- amr_norm$header
rownames(amr_gene_raw_analytic) <- amr_raw$header


# Make long data frame for plotting with ggplot2
amr_melted_analytic <- rbind(melt_dt(MRcounts(amr_class_analytic), 'Class'),
                             melt_dt(MRcounts(amr_mech_analytic), 'Mechanism'),
                             melt_dt(MRcounts(amr_group_analytic), 'Group'),
                             melt_dt(MRcounts(amr_gene_analytic), 'Gene'))

amr_melted_raw_analytic <- rbind(melt_dt(MRcounts(amr_class_raw_analytic), 'Class'),
                                 melt_dt(MRcounts(amr_mech_raw_analytic), 'Mechanism'),
                                 melt_dt(MRcounts(amr_group_raw_analytic), 'Group'),
                                 melt_dt(MRcounts(amr_gene_raw_analytic), 'Gene'))



# Split Kraken taxonomic lineages -----------------------------------------

# Aggregate the kraken data using the rownames:
# this set of commands splits the rownames into their taxonomic levels and
# fills empty values with NA.  We then join that taxonomy data table with
# the actual data and aggregate using lapply as before.

# kraken_taxonomy <- data.table(id=rownames(kraken))

kraken_taxonomy <- 
  temp_kraken_list %>%
  map(~ data.table(id=rownames(.x)))

lineage <- 
  kraken_taxonomy %>%
  map(~ .x$id)

kraken_taxonomy_split <- 
  kraken_taxonomy %>%
  map(~ str_split(string = .x$id, pattern = "\\|"))

domain_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^d_.*")]) %>%
  modify_depth(., .depth=2, ~ if(length(.x) == 0){.x=NA} else{.x})

phylum_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^p_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

class_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^c_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

order_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^o_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

family_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^f_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

genus_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^g_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

species_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^s_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

kraken_tax_dt_clade <- data.table(
  Domain = domain_tax$cladeReads,
  Phylum = phylum_tax$cladeReads,
  Class = class_tax$cladeReads,
  Order = order_tax$cladeReads,
  Family = family_tax$cladeReads,
  Genus = genus_tax$cladeReads,
  Species = species_tax$cladeReads
)

kraken_tax_dt_taxon <- data.table(
  Domain = domain_tax$taxonReads,
  Phylum = phylum_tax$taxonReads,
  Class = class_tax$taxonReads,
  Order = order_tax$taxonReads,
  Family = family_tax$taxonReads,
  Genus = genus_tax$taxonReads,
  Species = species_tax$taxonReads
)

kraken_tax_dt <- list(
  "cladeReads" = kraken_tax_dt_clade,
  "taxonReads" = kraken_tax_dt_taxon
)

kraken_tax_dt <- map2(
  kraken_tax_dt,
  lineage,
  ~ .x[, id := .y]
)

kraken_tax_dt <- 
  kraken_tax_dt %>%
  map(
    ~ .x[, lowest := str_split(id, pattern = "\\|") %>% map_chr(~ tail(.x, n=1))]
  )

# Use tax patterns named vector to replace lowest taxon name
# for taxonomy level

tax_levels <- c(
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

tax_regex <- c(
  "^d_.*",
  "^p_.*", 
  "^c_.*", 
  "^o_.*", 
  "^f_.*", 
  "^g_.*", 
  "^s_.*" 
)

# Zip tax levels and tax regex into named vector
tax_patterns <- tax_levels
names(tax_patterns) <- tax_regex

kraken_tax_dt <-
  kraken_tax_dt %>%
  map(
    ~ .x[,lowest_level := str_replace_all(lowest, tax_patterns)]
  )

kraken_tax_dt <- map(kraken_tax_dt, ~ setkey(.x, id))

kraken_norm <- map2(
  kraken_norm,
  kraken_css,
  ~ .x[, id :=(rownames(.y)), ]
)

kraken_norm <- map(
  kraken_norm,
  ~ setkey(.x, id)
)

kraken_norm <- map2(
  kraken_norm,
  kraken_tax_dt,
  ~ .y[.x] # left outer join
)

kraken_raw <- map2(
  kraken_raw,
  kraken_css,
  ~ .x[, id :=(rownames(.y)), ]
)

kraken_raw <- map(
  kraken_raw,
  ~ setkey(.x, id)
)
 
kraken_raw <- map2(
  kraken_raw,
  kraken_tax_dt,
  ~ .y[.x] # left outer join
) 

# Kraken taxon analytic matrices ------------------------------------------

# Group the kraken taxonReads data by level for analysis, removing NA entries

kraken_taxon_norm_summarised <- 
  tax_levels %>%
  map(
    ~ kraken_norm$taxonReads[!is.na(eval(as.name(.x))) & eval(as.name(.x)) != 'NA' , lapply(.SD, sum), by=.x, .SDcols=!1:10]
  ) %>% 
  set_names(nm=tax_levels)

kraken_taxon_raw_summarised <- 
  tax_levels %>%
  map(
    ~ kraken_raw$taxonReads[!is.na(eval(as.name(.x))) & eval(as.name(.x)) != 'NA' , lapply(.SD, sum), by=.x, .SDcols=!1:10]
  ) %>% 
  set_names(nm=tax_levels)

make_analytic <- function(x,y){
  analytic <- newMRexperiment(counts=x[, .SD, .SDcols=!y])
  rownames(analytic) <- x[,eval(as.name(y))]
  analytic
}

kraken_taxon_norm_analytic <- map2(
  kraken_taxon_norm_summarised,
  tax_levels,
  ~ make_analytic(.x, .y)
)

kraken_taxon_raw_analytic <- map2(
  kraken_taxon_raw_summarised,
  tax_levels,
  ~ make_analytic(.x, .y)
)


# Examples for further reference  
# kraken_species_raw <- kraken_raw[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
# kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!'Species'])
# rownames(kraken_species_raw_analytic) <- kraken_species_raw$Species

# Kraken clade analytic matrices ------------------------------------------

reorder_tax_ranks <- function(level_id){
  level_id <- factor(level_id,
                     levels = c(
                       "Domain",
                       "Phylum",
                       "Class",
                       "Order",
                       "Family",
                       "Genus",
                       "Species"
                     )
  )
  level_id
}



kraken_norm$cladeReads$lowest_level <- reorder_tax_ranks(kraken_norm$cladeReads$lowest_level)
kraken_raw$cladeReads$lowest_level <- reorder_tax_ranks(kraken_raw$cladeReads$lowest_level)
kraken_norm$taxonReads$lowest_level <- reorder_tax_ranks(kraken_norm$taxonReads$lowest_level)
kraken_raw$taxonReads$lowest_level <- reorder_tax_ranks(kraken_raw$taxonReads$lowest_level)

kraken_clade_norm_list <-
  kraken_norm$cladeReads %>%
  split(.$lowest_level) %>%
  map(
    ~ .x %>%
      select(-c(id:Species,lowest_level))
  )

kraken_clade_raw_list <-
  kraken_raw$cladeReads %>%
  split(.$lowest_level) %>%
  map(
    ~ .x %>%
      select(-c(id:Species,lowest_level))
  )

make_analytic_clade <- function(x){
  analytic <- newMRexperiment(counts=x[, .SD, .SDcols=!'lowest'])
  rownames(analytic) <- x$lowest
  analytic
}

kraken_clade_norm_analytic <- map(
  kraken_clade_norm_list,
  ~ make_analytic_clade(.x)
)

kraken_clade_raw_analytic <- map(
  kraken_clade_raw_list,
  ~ make_analytic_clade(.x)
)

# Make long data frame for plotting with ggplot2

kraken_taxon_norm_melted <- imap_dfr(
  kraken_taxon_norm_analytic,
  ~ melt_dt(MRcounts(.x), .y)
) # getting warning: binding character and factor vector, coercing into character vector

kraken_taxon_raw_melted <- imap_dfr(
  kraken_taxon_raw_analytic,
  ~ melt_dt(MRcounts(.x), .y)
) # getting warning: binding character and factor vector, coercing into character vector

kraken_clade_norm_melted <- imap_dfr(
  kraken_clade_norm_analytic,
  ~ melt_dt(MRcounts(.x), .y)
)

kraken_clade_raw_melted <- imap_dfr(
  kraken_clade_raw_analytic,
  ~ melt_dt(MRcounts(.x), .y)
) # getting warning: binding character and factor vector, coercing into character vector


# kraken_taxon_norm_melted$Level_ID <- reorder_tax_ranks(kraken_taxon_norm_melted$Level_ID)
# kraken_taxon_raw_melted$Level_ID <- reorder_tax_ranks(kraken_taxon_raw_melted$Level_ID)
# kraken_clade_norm_melted$Level_ID <- reorder_tax_ranks(kraken_clade_norm_melted$Level_ID)
# kraken_clade_raw_melted$Level_ID <- reorder_tax_ranks(kraken_clade_raw_melted$Level_ID)


# Match metadata ----------------------------------------------------------

# Ensure that the metadata entries match the factor order of the MRexperiments
metadata <- data.table(metadata[match(colnames(MRcounts(amr_class_analytic)), metadata[, "ID"]), ])
setkeyv(metadata, sample_column_id)

reorder_environments <- function(env_column){
  env_column <- factor(env_column,
                     levels = c(
                       "Fecal_Composite",
                       "Catch_Basin",
                       "Soil",
                       "Sewage_Treatment"
                      )
  )
}

metadata$Type <- reorder_environments(metadata$Type)


# Vector of objects for iteration and their names

AMR_analytic_data <- c(amr_class_analytic,
                       amr_mech_analytic,
                       amr_group_analytic,
                       amr_gene_analytic)

AMR_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')
AMR_raw_analytic_data <- c(amr_class_raw_analytic,
                           amr_mech_raw_analytic,
                           amr_group_raw_analytic,
                           amr_gene_raw_analytic)
AMR_raw_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')

for( l in 1:length(AMR_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(AMR_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(AMR_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(AMR_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(AMR_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_analytic_data[[l]])))
    rownames(fData(AMR_analytic_data[[l]])) <- rownames(MRcounts(AMR_analytic_data[[l]]))
}

for( l in 1:length(AMR_raw_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(AMR_raw_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(AMR_raw_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(AMR_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(AMR_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_raw_analytic_data[[l]])))
    rownames(fData(AMR_raw_analytic_data[[l]])) <- rownames(MRcounts(AMR_raw_analytic_data[[l]]))
}

names(AMR_analytic_data) <- AMR_analytic_names
names(AMR_raw_analytic_data) <- AMR_raw_analytic_names

for( l in 1:length(kraken_taxon_norm_analytic) ) {
    sample_idx <- match(colnames(MRcounts(kraken_taxon_norm_analytic[[l]])), metadata[[sample_column_id]])
    pData(kraken_taxon_norm_analytic[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(kraken_taxon_norm_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(kraken_taxon_norm_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_taxon_norm_analytic[[l]])))
    rownames(fData(kraken_taxon_norm_analytic[[l]])) <- rownames(MRcounts(kraken_taxon_norm_analytic[[l]]))
}
# 
for( l in 1:length(kraken_taxon_raw_analytic) ) {
    sample_idx <- match(colnames(MRcounts(kraken_taxon_raw_analytic[[l]])), metadata[[sample_column_id]])
    pData(kraken_taxon_raw_analytic[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(kraken_taxon_raw_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(kraken_taxon_raw_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_taxon_raw_analytic[[l]])))
    rownames(fData(kraken_taxon_raw_analytic[[l]])) <- rownames(MRcounts(kraken_taxon_raw_analytic[[l]]))
}

kraken_taxon_norm_analytic <-
  kraken_taxon_norm_analytic %>%
  map(function(x){
    pData(x)[["Type"]] <- reorder_environments(pData(x)[["Type"]])
    x
  }
)

kraken_taxon_raw_analytic <-
  kraken_taxon_raw_analytic %>%
  map(function(x){
    pData(x)[["Type"]] <- reorder_environments(pData(x)[["Type"]])
    x
  }
)

kraken_clade_norm_analytic <-
  kraken_clade_norm_analytic %>%
  map(function(x){
    pData(x)[["Type"]] <- reorder_environments(pData(x)[["Type"]])
    x
  }
)

kraken_clade_raw_analytic <-
  kraken_clade_raw_analytic %>%
  map(function(x){
    pData(x)[["Type"]] <- reorder_environments(pData(x)[["Type"]])
    x
  }
)
# for( l in 1:length(kraken_clade_norm_analytic) ) {
#     sample_idx <- match(colnames(MRcounts(kraken_clade_norm_analytic[[l]])), metadata[[sample_column_id]])
#     pData(kraken_clade_norm_analytic[[l]]) <- data.frame(
#         metadata[sample_idx, .SD, .SDcols=!sample_column_id])
#     rownames(pData(kraken_clade_norm_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
#     fData(kraken_clade_norm_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_clade_norm_analytic[[l]])))
#     rownames(fData(kraken_clade_norm_analytic[[l]])) <- rownames(MRcounts(kraken_clade_norm_analytic[[l]]))
# }
# 
# for( l in 1:length(kraken_clade_raw_analytic) ) {
#     sample_idx <- match(colnames(MRcounts(kraken_clade_raw_analytic[[l]])), metadata[[sample_column_id]])
#     pData(kraken_clade_raw_analytic[[l]]) <- data.frame(
#         metadata[sample_idx, .SD, .SDcols=!sample_column_id])
#     rownames(pData(kraken_clade_raw_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
#     fData(kraken_clade_raw_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_clade_raw_analytic[[l]])))
#     rownames(fData(kraken_clade_raw_analytic[[l]])) <- rownames(MRcounts(kraken_clade_raw_analytic[[l]]))
# }

match_metadata <- function(x, meta){
 sample_idx <- match(colnames(MRcounts(x)), meta[[sample_column_id]])
 pData(x) <- data.frame(
   meta[sample_idx, .SD, .SDcols=!sample_column_id]
   )
 rownames(pData(x)) <- meta[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
 fData(x) <- data.frame(Feature=rownames(MRcounts(x)))
 rownames(fData(x)) <- rownames(MRcounts(x))
 x
}

# Normalized Kraken analytic matrices

kraken_taxon_norm_analytic <- 
  as.vector(kraken_taxon_norm_analytic) %>%
  map(
    ~ match_metadata(.x)
)

kraken_clade_norm_analytic <- 
as.vector(kraken_clade_norm_analytic) %>%
  map(
    ~ match_metadata(.x)
  )

# Raw Kraken analytic matrices

kraken_taxon_raw_analytic <- 
  as.vector(kraken_taxon_raw_analytic) %>%
  map(
    ~ match_metadata(.x)
  )

kraken_clade_raw_analytic <- 
  as.vector(kraken_clade_raw_analytic) %>%
  map(
    ~ match_metadata(.x)
  )

kraken_taxon_names <- names(kraken_taxon_raw_analytic)
kraken_clade_names <- names(kraken_clade_raw_analytic)



# Exploratory Analyses: Alpha Rarefaction ---------------------------------

for( v in 1:length(exploratory_analyses) ) {
    # AMR
    meg_alpha_rarefaction(data_list=AMR_raw_analytic_data,
                          data_names=AMR_raw_analytic_names,
                          metadata=metadata,
                          sample_var=sample_column_id,
                          group_var=exploratory_analyses[[v]]$exploratory_var,
                          analysis_subset=exploratory_analyses[[v]]$subsets,
                          outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                       sep='/', collapse=''),
                          data_type='AMR')
    
    # Microbiome (taxon reads)
    meg_alpha_rarefaction(data_list=kraken_taxon_raw_analytic,
                          data_names=kraken_taxon_names,
                          metadata=metadata,
                          sample_var=sample_column_id,
                          group_var=exploratory_analyses[[v]]$exploratory_var,
                          analysis_subset=exploratory_analyses[[v]]$subsets,
                          outdir=paste(graph_output_dir, 'Microbiome_taxonReads', exploratory_analyses[[v]]$name,
                                       sep='/', collapse=''),
                          data_type='Microbiome_taxonReads')
    
}    

exploratory_analyses %>%
  walk(safely(~ meg_alpha_rarefaction(
      data_list=kraken_clade_raw_analytic,
      data_names=kraken_clade_names,
      metadata=metadata,
      sample_var=sample_column_id,
      group_var=.x$exploratory_var,
      analysis_subset=.x$subsets,
      outdir=paste(graph_output_dir, 'Microbiome_cladeReads', .x$name,
        sep='/', collapse=''),
      data_type='Microbiome_cladeReads')
  ))


# Exploratory Analyses: Alpha Normalized ----------------------------------

# widen_amr <- function(x){
#   amr_norm_wide <- spread(x, key=categoryNames, value=normCountsSum, fill=0)
#   row.names(amr_norm_wide) <- amr_norm_wide$samples
#   amrNormWide <- amrNormWide %>%
#     select(2:ncol(amrNormWide))
#   return(amrNormWide)
# }

# amr_norm_div_mat <- map(AMR_analytic_data, ~ widen_amr(MRcounts(.x)))

calc_diversity_df <- function(x){
  observed_richness <- specnumber(x, MARGIN=2)
  invsimpson <- diversity(x, index="invsimpson", MARGIN=2)
  simpson <- diversity(x, index="simpson", MARGIN=2)
  shannon <- diversity(x, index="shannon", MARGIN=2)
  evenness <- shannon/log(observed_richness)
  div_df <- data.frame(
    ID = names(observed_richness),
    Observed_Richness = observed_richness,
    Inv_Simpson = invsimpson,
    Simpson = simpson,
    Shannon = shannon,
    Evenness = evenness
  )
  div_df
}

amr_norm_diversity <- 
  AMR_analytic_data %>%
  set_names(nm = AMR_analytic_names) %>%
  map_dfr(
    ~ calc_diversity_df(MRcounts(.x)) %>%
        left_join(., metadata, by = "ID"),
    .id = "Level"
  )

kraken_clade_norm_diversity <- 
  kraken_clade_norm_analytic %>%
  map(
    ~ calc_diversity_df(MRcounts(.x))
  )

kraken_taxon_norm_diversity <- 
  kraken_taxon_norm_analytic %>%
  map(
    ~ calc_diversity_df(MRcounts(.x))
  )

amr_observed_species <- function(amr_df) {
  alphaDivBoxPlot <-
    ggplot(amr_df, aes(Type, Observed_Richness, color = Type)) +
    geom_boxplot(lwd = 0.5, aes(size = 5)) +
    theme(
      strip.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 30),
      axis.text.x = element_text(size = 25, angle = 90),
      axis.title.x = element_text(size = 32),
      axis.title.y = element_text(size = 32),
      legend.position = "none",
      legend.title = element_text(size = 36),
      legend.text = element_text(size = 36, vjust = 0.5),
      plot.title = element_text(size = 50, hjust = 0.5)
    ) +
    xlab("Type") +
    ylab("Number of Unique\nAMR Categories\n") +
    #ggtitle('AMR Category Richness by Depth for Raw Data') +
    # scale_color_manual(values=rev(cbPalette)) +
    facet_wrap( ~ Level, nrow = 2, scales = "free_y")
}

amr_norm_rich_boxplots <- amr_observed_species(amr_norm_diversity)

amr_inv_simpson <- function(amr_df) {
  alphaDivBoxPlot <-
    ggplot(amr_df, aes(Type, Inv_Simpson, color = Type)) +
    geom_boxplot(lwd = 0.5, aes(size = 5)) +
    theme(
      strip.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 30),
      axis.text.x = element_text(size = 25, angle = 90),
      axis.title.x = element_text(size = 32),
      axis.title.y = element_text(size = 32),
      legend.position = "none",
      #legend.title=element_text(size=36),
      #legend.text=element_text(size=36, vjust=0.5),
      plot.title = element_text(size = 50, hjust = 0.5)
    ) +
    xlab("Type") +
    ylab("Inverse Simpson Index") +
    #ggtitle('AMR Category Richness by Depth for Raw Data') +
    # scale_color_manual(values=rev(cbPalette)) +
    facet_wrap( ~ Level, nrow = 2, scales = "free_y")
}

amr_norm_inv_simpson_box <- amr_inv_simpson(amr_norm_diversity)

amr_simpson <- function(amr_df) {
  alphaDivBoxPlot <- ggplot(amr_df, aes(Type, Simpson, color = Type)) +
    geom_boxplot(lwd = 0.5, aes(size = 5)) +
    theme(
      strip.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 30),
      axis.text.x = element_text(size = 25, angle = 90),
      axis.title.x = element_text(size = 32),
      axis.title.y = element_text(size = 32),
      legend.position = "none",
      #legend.title=element_text(size=36),
      #legend.text=element_text(size=36, vjust=0.5),
      plot.title = element_text(size = 50, hjust = 0.5)
    ) +
    xlab("Type") +
    ylab("Inverse Simpson Index") +
    #ggtitle('AMR Category Richness by Depth for Raw Data') +
    # scale_color_manual(values=rev(cbPalette)) +
    facet_wrap( ~ Level, nrow = 2, scales = "free_y")
}

amr_norm_simpson_box <- amr_simpson(amr_norm_diversity)

amr_shannon <- function(amr_df) {
  alphaDivBoxPlot <- ggplot(amr_df, aes(Type, Shannon, color = Type)) +
    geom_boxplot(lwd = 0.5, aes(size = 5)) +
    theme(
      strip.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 30),
      axis.text.x = element_text(size = 25, angle = 90),
      axis.title.x = element_text(size = 32),
      axis.title.y = element_text(size = 32),
      legend.position = "none",
      #legend.title=element_text(size=36),
      #legend.text=element_text(size=36, vjust=0.5),
      plot.title = element_text(size = 50, hjust = 0.5)
    ) +
    xlab("Type") +
    ylab("Shannon Index") +
    #ggtitle('AMR Category Richness by Depth for Raw Data') +
    # scale_color_manual(values=rev(cbPalette)) +
    facet_wrap( ~ Level, nrow = 2, scales = "free_y")
}

amr_norm_shannon_box <- amr_shannon(amr_norm_diversity)

exploratory_analyses %>%
  walk(safely(
    ~ meg_alpha_rarefaction(
      data_list = kraken_clade_raw_analytic,
      data_names = kraken_clade_names,
      metadata = metadata,
      sample_var = sample_column_id,
      group_var = .x$exploratory_var,
      analysis_subset = .x$subsets,
      outdir = paste(
        graph_output_dir,
        'Microbiome_cladeReads',
        .x$name,
        sep = '/',
        collapse = ''
      ),
      data_type = 'Microbiome_cladeReads'
    )
  ))


# Exploratory Analyses: Ordination ----------------------------------------

for( v in 1:length(exploratory_analyses) ) {
    # AMR NMDS
    meg_ordination(data_list = AMR_analytic_data,
                   data_names = AMR_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_analyses[[v]]$exploratory_var,
                   analysis_subset=exploratory_analyses[[v]]$subsets,
                   outdir = paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                  sep='/', collapse=''),
                   data_type = 'AMR',
                   method = 'NMDS')
    
    # AMR PCA
    meg_ordination(data_list = AMR_analytic_data,
                   data_names = AMR_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_analyses[[v]]$exploratory_var,
                   analysis_subset=exploratory_analyses[[v]]$subsets,
                   outdir = paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                  sep='/', collapse=''),
                   data_type = 'AMR',
                   method = 'PCA')
    
}
    
    # Microbiome NMDS (taxon reads)

# Mapping quietly allows to capture output on procrustes and stress
# It also allows to store warnings
# Consider using walk with capture.output to store the procrustes information

micro_taxon_nmds <- exploratory_analyses %>%
  map(
    quietly(
    ~ meg_ordination(
      data_list = kraken_taxon_norm_analytic,
      data_names = kraken_taxon_names,
      metadata = metadata,
      sample_var = sample_column_id,
      hull_var = .x$exploratory_var,
      analysis_subset=.x$subsets,
      outdir = paste(graph_output_dir, 'Microbiome_taxonReads', .x$name, sep='/', collapse=''),
      data_type = 'Microbiome',
      method = 'NMDS')
    ))

# dev.off() ?
    
    # Microbiome PCA
exploratory_analyses %>%
  walk(
    safely(
    ~ meg_ordination(data_list = kraken_taxon_norm_analytic,
                   data_names = kraken_taxon_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = .x$exploratory_var,
                   analysis_subset=.x$subsets,
                   outdir = paste(graph_output_dir, 'Microbiome_taxonReads', .x$name, sep='/', collapse=''),
                   data_type = 'Microbiome_taxonReads',
                   method = 'PCA')
  ))

# Microbiome NMDS (clade Reads)

exploratory_analyses %>%
  walk(
    safely(
    ~ meg_ordination(
      data_list = kraken_clade_norm_analytic,
      data_names = kraken_clade_names,
      metadata = metadata,
      sample_var = sample_column_id,
      hull_var = .x$exploratory_var,
      analysis_subset=.x$subsets,
      outdir = paste(graph_output_dir, 'Microbiome_cladeReads', .x$name, sep='/', collapse=''),
      data_type = 'Microbiome',
      method = 'NMDS')
    ))
    
    # Microbiome PCA
exploratory_analyses %>%
  walk(
    safely(
    ~ meg_ordination(data_list = kraken_clade_norm_analytic,
                   data_names = kraken_clade_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = .x$exploratory_var,
                   analysis_subset=.x$subsets,
                   outdir = paste(graph_output_dir, 'Microbiome_cladeReads', .x$name, sep='/', collapse=''),
                   data_type = 'Microbiome_cladeReads',
                   method = 'PCA')
  ))

# Exploratory Analyses: Heatmaps ------------------------------------------

# AMR Heatmaps for each level
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        meg_heatmap(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=AMR_analytic_names[[l]],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='AMR')
    }
}

# Microbiome (taxon reads)
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(kraken_taxon_names) ) {
        meg_heatmap(melted_data=kraken_taxon_norm_melted,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=kraken_taxon_names[[l]],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'Microbiome_taxonReads', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='Microbiome_taxonReads')
    }
}

# Microbiome (clade reads)
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(kraken_clade_names) ) {
        meg_heatmap(melted_data=kraken_clade_norm_melted,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=kraken_clade_names[[l]],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'Microbiome_cladeReads', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='Microbiome_cladeReads')
    }
}

# Exploratory Analyses: Barplots ------------------------------------------

# AMR
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=AMR_analytic_names[l],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='AMR')
        )
    }
}

# Microbiome (taxon reads)
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(kraken_taxon_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=kraken_taxon_norm_melted,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=kraken_taxon_names[[l]],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'Microbiome_taxonReads', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='Microbiome_taxonReads')
        )
    }
}

# Microbiome (clade reads)

for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(kraken_clade_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=kraken_clade_norm_melted,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=kraken_taxon_names[[l]],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'Microbiome_cladeReads', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='Microbiome_cladeReads')
        )
    }
}

# Statistical Analyses ----------------------------------------------------

for( a in 1:length(statistical_analyses) ) {
    meg_fitZig(data_list=AMR_analytic_data,
               data_names=AMR_analytic_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(amr))),
               data_mod=statistical_analyses[[a]]$model_matrix,
               filter_min_threshold=0.15,
               contrast_list=statistical_analyses[[a]]$contrasts,
               random_effect_var=statistical_analyses[[a]]$random_effect,
               outdir=paste(stats_output_dir, 'AMR', statistical_analyses[[a]]$name,
                            sep='/', collapse=''),
               analysis_name=statistical_analyses[[a]]$name,
               analysis_subset=statistical_analyses[[a]]$subsets,
               data_type='AMR',
               pval=0.1,
               top_hits=1000)
}

for (a in 1:length(statistical_analyses)){
    meg_fitZig(data_list=kraken_taxon_norm_analytic,
               data_names=kraken_taxon_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(kraken_css$taxonReads))),
               data_mod=statistical_analyses[[a]]$model_matrix,
               filter_min_threshold=0.15,
               contrast_list=statistical_analyses[[a]]$contrasts,
               random_effect_var=statistical_analyses[[a]]$random_effect,
               outdir=paste(stats_output_dir, 'Microbiome_taxonReads', statistical_analyses[[a]]$name,
                            sep='/', collapse=''),
               analysis_name=statistical_analyses[[a]]$name,
               analysis_subset=statistical_analyses[[a]]$subsets,
               data_type='Microbiome_taxonReads',
               pval=0.1,
               top_hits=1000)
}

for (a in 1:length(statistical_analyses)){
    meg_fitZig(data_list=kraken_clade_norm_analytic,
               data_names=kraken_clade_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(kraken_css$cladeReads))),
               data_mod=statistical_analyses[[a]]$model_matrix,
               filter_min_threshold=0.15,
               contrast_list=statistical_analyses[[a]]$contrasts,
               random_effect_var=statistical_analyses[[a]]$random_effect,
               outdir=paste(stats_output_dir, 'Microbiome_cladeReads', statistical_analyses[[a]]$name,
                            sep='/', collapse=''),
               analysis_name=statistical_analyses[[a]]$name,
               analysis_subset=statistical_analyses[[a]]$subsets,
               data_type='Microbiome_cladeReads',
               pval=0.1,
               top_hits=1000)
}


# Output of matrices ------------------------------------------------------

# Attempt to include purrr functional programming approach


write.csv(make_sparse(amr_class, 'class', c('class')), 'amr_matrices/sparse_normalized/AMR_Class_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_class, 'amr_matrices/normalized/AMR_Class_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_class_raw, 'amr_matrices/raw/AMR_Class_Raw.csv', sep=',', row.names = F, col.names = T)


write.csv(make_sparse(amr_mech, 'mechanism', c('mechanism')), 'amr_matrices/sparse_normalized/AMR_Mechanism_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_mech, 'amr_matrices/normalized/AMR_Mechanism_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Mechanism_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_group, 'group', c('group')), 'amr_matrices/sparse_normalized/AMR_Group_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_group, 'amr_matrices/normalized/AMR_Group_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Group_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_norm, 'header', c('header', 'class', 'mechanism', 'group')),
          'amr_matrices/sparse_normalized/AMR_Gene_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_norm, 'amr_matrices/normalized/AMR_Gene_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_raw, 'amr_matrices/raw/AMR_Gene_Raw.csv', sep=',', row.names = F, col.names = T)


kraken_taxon_norm_summarised %>%
  iwalk(
    ~ write.csv(
      make_sparse(.x, .y, c(.y)), 
      here('kraken_taxonReads_matrices', 'sparse_normalized', paste0('kraken_',.y,'_Sparse_Normalized.csv')),
      row.names = T)
    )

kraken_taxon_norm_summarised %>%
  iwalk(
    ~ write.csv(
      .x, 
      here('kraken_taxonReads_matrices', 'normalized', paste0('kraken_',.y,'_Normalized.csv')),
      row.names = T)
  )

kraken_taxon_raw_summarised %>%
  iwalk(
    ~ write.csv(
      .x,
      here('kraken_taxonReads_matrices', 'raw', paste0('kraken_',.y,'_Raw.csv')),
      row.names = T)
    )


# Kraken clade reads

kraken_clade_norm_list <-
  kraken_clade_norm_list %>%
  imap(
    ~ rename(.x, !!.y := lowest) # quosure notation; tip from SO https://stackoverflow.com/questions/46616591/rename-multiple-dataframe-columns-using-purrr
  )

kraken_clade_raw_list <-
  kraken_clade_raw_list %>%
  imap(
    ~ rename(.x, !!.y := lowest)
  )

kraken_clade_norm_list %>%
  iwalk(~ write.csv(
    make_sparse(.x, .y, c(.y)), 
    here('kraken_cladeReads_matrices', 'sparse_normalized', paste0('kraken_',.y,'_Sparse_Normalized.csv')),
    row.names = T)
  )

kraken_clade_norm_list %>%
  iwalk(
    ~ write.csv(
      .x,
      here('kraken_cladeReads_matrices', 'normalized', paste0('kraken_',.y,'_Normalized.csv')),
      row.names = T)
    )

kraken_clade_raw_list %>%
  iwalk(
    ~ write.csv(
      .x,
      here('kraken_cladeReads_matrices', 'raw', paste0('kraken_',.y,'_Raw_Normalized.csv')),
      row.names = T)
  )
