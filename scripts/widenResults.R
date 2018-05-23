
# Load packages -----------------------------------------------------------

library(tidyverse)
library(pavian)
library(here)



# Create directory for aggregated data ------------------------------------

aggregated_dir <- here("aggregated_data_for_analysis")

ifelse(
  !dir.exists(aggregated_dir), 
  dir.create((aggregated_dir), mode='777'), 
FALSE
)

# Read Kraken reports with Pavian ---------------------------------------------

# Obtain the filenames of all Kraken reports

krakenReportPaths <- Sys.glob(here("post_PhiXFilter_Kraken",'*.tabular'))

krakenReportNames <- list.files(path = here("post_PhiXFilter_Kraken"), pattern = '*.tabular')

krakenReportNames <- 
  krakenReportNames %>%
  map(function(x) str_replace(x, "\\.tabular$", ""))

krakenReportsPavian <-
  krakenReportPaths %>%
  map(function(x) pavian::read_report(x)) %>%
  set_names(nm = krakenReportNames)
  
krakenReportsPavian <-
  krakenReportsPavian %>%
  map(safely(function(x){
    filt <- pavian::filter_taxon(
    report = x, 
    filter_taxon = "Eukaryota", 
    rm_clade = TRUE, 
    do_message = TRUE
    )
    filt
}))

taxa_to_remove <- c("u_unclassified", "-_root", "-_cellular organisms")

krakenReportsPavianMerged <- 
  krakenReportsPavian %>%
  map_dfr(function(x){
    x <- x$result %>%
      filter(!name %in% taxa_to_remove) %>%
      filter(taxRank != "-") %>%
      mutate(taxLineage=str_replace(taxLineage, "-_root\\|-_cellular organisms\\|", "")) %>%
      mutate(taxLineage=str_replace(taxLineage, "-_root\\|", ""))
    x
    }, .id = "Sample")

make_kraken_analytical <- function(x, column){
  x <- x %>%
  select(Sample, column, taxLineage) %>%
  spread(key = Sample, value = column, fill = 0) %>%
  rename(Lineage = taxLineage)
  x
}

tax_columns <- c("cladeReads", "taxonReads")
  
krakenAnalytical <- map(tax_columns, ~ make_kraken_analytical(krakenReportsPavianMerged, .x)) %>%
  set_names(nm=tax_columns)
  
iwalk(krakenAnalytical,
      ~ write.csv(.x, here("aggregated_data_for_analysis", paste0("krakenAnalytical", "_",.y,".csv")), row.names = FALSE))


# Read AMR and MegaBio Coverage Sampler Results -------------------------------

# Parse the results with Python script first

# Then collect the names of the parsed files

# Make sure the fecal composite data is present as well.

amrCovSamplerPaths <- Sys.glob(here("AMR_CovSampler_parsed", '*CovSampler_parsed.tab'))

amrCovSamplerNames <- list.files(
  path = here("./AMR_CovSampler_parsed"),
  pattern = "*CovSampler_parsed.tab"
  ) %>%
  map(function(x) str_replace(x, "_CovSampler_parsed\\.tab$", ""))

amrCovSampler <- amrCovSamplerPaths %>%
  map(function(x) read_tsv(x)) %>% 
      set_names(nm=amrCovSamplerNames)

amrReportsMerged <- amrCovSampler %>%
  map_dfr(function(x) x, .id="Sample")

amrAnalytical <- amrReportsMerged %>%
  select(Sample, Header, Hits) %>%
  spread(key = Sample, value = Hits, fill = 0)

amrClassification <- amrAnalytical$Header

amrAnalytical <- amrAnalytical %>%
  select(-Header) %>%
  as.matrix(.)

row.names(amrAnalytical) <- amrClassification


# Reading MEGABio data ----------------------------------------------------

megaBioPaths <- Sys.glob(here("MegaBio_results", '*','*.tabular'))

megaBioNames <- list.files(
  path = Sys.glob("./MegaBio_results/*"),
  pattern = "*CovSampler_parsed\\.tab$"
  ) %>%
  map(function(x) str_replace(x, "_MBio_CovSampler_parsed\\.tab$", ""))

megaBioReports <- 
  megaBioPaths %>%
  map(~ read_tsv(.x)) %>% 
      set_names(nm=megaBioNames)

megaBioReportsMerged <- megaBioReports %>%
  map_dfr(function(x) x, .id="Sample") %>%
  rename(Header = `Gene Id`)

# Change sample names to reflect names in metadata files

megaBioReportsMerged <- megaBioReportsMerged %>%
  mutate(Sample = str_replace(Sample, "FC_A062_H_006", "FC_006_A062")) %>%
  mutate(Sample = str_replace(Sample, "FC_A062_H_007", "FC_007_A062")) %>%
  mutate(Sample = str_replace(Sample, "FC_A070_H_006", "FC_006_A070")) %>%
  mutate(Sample = str_replace(Sample, "FC_A070_H_007", "FC_007_A070")) %>%
  mutate(Sample = str_replace(Sample, "FC_N003_H_007", "FC_007_N003")) %>%
  mutate(Sample = str_replace(Sample, "FC_N003_H_008", "FC_008_N003")) %>%
  mutate(Sample = str_replace(Sample, "FC_N013_H_007", "FC_007_N013")) %>%
  mutate(Sample = str_replace(Sample, "FC_N013_H_008", "FC_008_N013")) %>%
  mutate(Sample = str_replace(Sample, "FC_S034_H_007", "FC_007_S034")) %>%
  mutate(Sample = str_replace(Sample, "FC_S034_H_008", "FC_008_S034")) %>%
  mutate(Sample = str_replace(Sample, "FC_S040_H_007", "FC_007_S040")) %>%
  mutate(Sample = str_replace(Sample, "FC_S040_H_008", "FC_008_S040")) %>%
  mutate(Sample = str_replace(Sample, "FC_S040_H_007", "FC_007_S040")) %>%
  mutate(Sample = str_replace(Sample, "FC_S040_H_008", "FC_008_S040")) %>%
  mutate(Sample = str_replace(Sample, "FC_V042_Nat", "FC_Nat_V042")) %>%
  mutate(Sample = str_replace(Sample, "FC_V045_Nat", "FC_Nat_V045")) %>%
  mutate(Sample = str_replace(Sample, "FC_V052_Nat", "FC_Nat_V052")) %>%
  mutate(Sample = str_replace(Sample, "FC_V053_H_006", "FC_006_V053")) %>%
  mutate(Sample = str_replace(Sample, "FC_V053_H_007", "FC_007_V053")) %>%
  mutate(Sample = str_replace(Sample, "FC_V055", "FC_Con_V055")) %>%
  mutate(Sample = str_replace(Sample, "FC_V046_H_006", "FC_006_V046")) %>%
  mutate(Sample = str_replace(Sample, "FC_V046_H_007", "FC_007_V046"))
  

# Concatenate MEGARes and MEGABio data ------------------------------------

amrReportsMerged <- 
  amrReportsMerged %>%
  select(-c(Class, Mechanism, Group))

amrBioConcat <- rbind(amrReportsMerged, megaBioReportsMerged)

amrBioAnalytical <- amrBioConcat %>%
  select(Sample, Header, Hits) %>%
  spread(key = Sample, value = Hits, fill = 0)

amrBioClassification <- amrBioAnalytical$Header

amrBioAnalytical <- amrBioAnalytical %>%
  select(-Header) %>%
  as.matrix(.)

row.names(amrBioAnalytical) <- amrBioClassification

write.csv(amrBioAnalytical, here('aggregated_data_for_analysis', 'amrBioAnalytical.csv'))

# Update annotations file with new MEGABio annotations ---------------------

megaresMegabioCSU <- read.csv('megares_annotations_v1.01.csv')

megaBioAAFC <- read.csv('megabio_AAFC_v0.2_annotation.csv')

megaresAMR <- megaresMegabioCSU %>%
  slice(1:3824) # row at which AMR ends. Biometal/biociode starts in next row

megaresMegabioUpdated <- rbind(megaresAMR, megaBioAAFC)

write.csv(megaresMegabioUpdated, here('aggregated_data_for_analysis', 'megaresMegabioUpdated.csv'), row.names = FALSE)
