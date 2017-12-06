## R script template for the production of basic qualitative
## statistics on the AMR and kraken output from AMR++

## Author: Steven Lakin


## The files you want to use for input to this (for the MEG group analyses)
## are the AMR_analytic_matrix.csv and kraken_analytic_matrix.csv.  The AMR
## matrix is identical to the Gene.csv matrix, however the kraken analytic
## matrix is not due to the way that reads get classified using each of
## these methods.

## So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine.  We will assume
## that you've set your working directory to where these files are located on your
## local machine and that you have installed the metagenomeSeq package.

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.


###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# Set your working directory to the main folder for analysis:
setwd('/home/lakinsm/Documents/morleyBioinformatics/CanadaAnalyticData/20March2017Analysis/')

# Set the output directory for graphs:
graph_output_dir = 'graphs'


# Set the output directory for statistics:
stats_output_dir = 'stats'


# Where is the metadata file stored on your machine?
metadata_filepath = 'BCRC_metadata.csv'


# Name of the megares annotation file used for this project
megares_annotation_filename = 'megares_annotations_v1.01.csv'


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
        subsets = list('Type == Fecal.Composite', 'NatType != Natural'),
        exploratory_var = 'Location'
    ),
    
    # # Analysis 4
    # # Description: Location comparison within Catch Basin type
    # list(
    #     name = 'LocationCB',
    #     subsets = list('Type == Catch.Basin', 'NatType != Natural'),
    #     exploratory_var = 'Location'
    # ),
    
    # Analysis 5
    # Description: Location comparison within Waste Water sewage treatment type
    list(
        name = 'LocationST',
        subsets = list('Type == Sewage.Treatment'),
        exploratory_var = 'Location'
    ),
    
    # Analysis 6
    # Description: Natural vs conventional Fecal Composite
    list(
        name = 'NaturalConventionalFC',
        subsets = list('Type == Fecal.Composite', 'Location == Vegreville', 'NatType != None'),
        exploratory_var = 'NatType'
    ),
    
    # Analysis 7
    # Description: Natural vs conventional Catch Basin
    list(
        name = 'NaturalConventionalCB',
        subsets = list('Type == Fecal.Composite', 'Location == Vegreville'),
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
        contrasts = list('TypeFecal.Composite - TypeCatch.Basin',
                         'TypeFecal.Composite - TypeSewage.Treatment',
                         'TypeCatch.Basin - TypeSewage.Treatment'),
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
                       'Type != Sewage.Treatment',
                       'Type != Catch.Basin',
                       'NatType != None',
                       'Location == Vegreville'),
        model_matrix = '~ 0 + NatType',
        contrasts = list('NatTypeNatural - NatTypeConventional'),
        random_effect = NA
    )
)






####################
## Automated Code ##
####################
## Modify this as necessary, though you shouldn't need to for basic use.

# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')

require(metagenomeSeq)
require(data.table)
require(ggplot2)
require(vegan)

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
              'TUFAB')


##########################
## Import & Format Data ##
##########################
## These files should be standard for all analyses, as they are
## the output matrices from AMR++ nextflow.  Additionally,
## you will need to obtain the most recent megares annotations file
## from megares.meglab.org


# If subdirs for stats and exploratory variables don't exist, create them
ifelse(!dir.exists(file.path(graph_output_dir)), dir.create(file.path(graph_output_dir), mode='777'), FALSE)
ifelse(!dir.exists(file.path(stats_output_dir)), dir.create(file.path(stats_output_dir), mode='777'), FALSE)

for( dtype in c('AMR', 'Microbiome') ) {
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

ifelse(!dir.exists(file.path('50amr_matrices')), dir.create(file.path('50amr_matrices'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('50amr_matrices/sparse_normalized')), dir.create(file.path('50amr_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('50amr_matrices/normalized')), dir.create(file.path('50amr_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('50amr_matrices/raw')), dir.create(file.path('50amr_matrices/raw'), mode='777'), FALSE)



# Load the data, MEGARes annotations, and metadata
amr <- newMRexperiment(read.table('50AMR_biometal_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- data.table(read.csv(megares_annotation_filename, header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys

metadata <- read.csv(metadata_filepath, header=T)
metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])


# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore
cumNorm(amr)


# Extract the normalized counts into data tables for aggregation
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


# Ensure that the metadata entries match the factor order of the MRexperiments
metadata <- data.table(metadata[match(colnames(MRcounts(amr_class_analytic)), metadata[, sample_column_id]), ])
setkeyv(metadata, sample_column_id)


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



#############################################
## Exploratory Analyses: Alpha Rarefaction ##
#############################################
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
}


######################################
## Exploratory Analyses: Ordination ##
######################################
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


####################################
## Exploratory Analyses: Heatmaps ##
####################################

# AMR Heatmaps for each level
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        meg_heatmap(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=AMR_analytic_names[l],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='AMR')
    }
}


####################################
## Exploratory Analyses: Barplots ##
####################################

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


##########################
## Statistical Analyses ##
##########################
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


########################
## Output of matrices ##
########################
write.csv(make_sparse(amr_class, 'class', c('class')), '50amr_matrices/sparse_normalized/AMR_Class_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_class, '50amr_matrices/normalized/AMR_Class_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_class_raw, '50amr_matrices/raw/AMR_Class_Raw.csv', sep=',', row.names = F, col.names = T)


write.csv(make_sparse(amr_mech, 'mechanism', c('mechanism')), '50amr_matrices/sparse_normalized/AMR_Mechanism_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_mech, '50amr_matrices/normalized/AMR_Mechanism_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, '50amr_matrices/raw/AMR_Mechanism_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_group, 'group', c('group')), '50amr_matrices/sparse_normalized/AMR_Group_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_group, '50amr_matrices/normalized/AMR_Group_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, '50amr_matrices/raw/AMR_Group_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_norm, 'header', c('header', 'class', 'mechanism', 'group')),
          '50amr_matrices/sparse_normalized/AMR_Gene_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_norm, '50amr_matrices/normalized/AMR_Gene_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_raw, '50amr_matrices/raw/AMR_Gene_Raw.csv', sep=',', row.names = F, col.names = T)


