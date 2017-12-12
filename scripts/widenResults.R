library(readr)
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)

# Read Kraken reports with Pavian ---------------------------------------------

# Obtain the filenames of all Kraken reports

krakenReportPaths <- Sys.glob(file.path(".",
                                    "post_PhiXFilter_Kraken",
                                    '*.tabular'))

krakenReportNames <- list.files(path = Sys.glob("./post_PhiXFilter_Kraken/"),
                            pattern = "*.tabular")

krakenReportNames <- krakenReportNames %>%
  map(function(x) str_replace(x, "\\.tabular$", ""))


krakenReports <- krakenReportPaths %>%
  map(read_tsv)
  # set_names(nm=krakenReportNames)

# krakenReportsPavian <- lapply(krakenReportPaths, function(x) pavian::read_report(x))
  
krakenReportsPavian <- krakenReportPaths %>%
  map(function(x) pavian::read_report(x)) %>%
  set_names(nm=krakenReportNames)

# krakenReportsPavianFilt <- krakenReportsPavian %>%
#   map(function(x) pavian::filter_cladeReads(cladeReads = "cladeReads",
#                                             tax_data=x, 
#                                             rm_taxa = c("u_unclassifed", 
#                                                         "p_Chordata")))

krakenReportsPavianMerged <- krakenReportsPavian %>%
  map_dfr(function(x){
    x <- x %>%
      filter(name != "u_unclassified") # Filter names of tax ranks (e.g. D,P)
    x
    }, .id = "Sample")

krakenAnalytical <- krakenReportsPavianMerged %>%
  select(Sample, cladeReads, taxLineage) %>%
  spread(key = Sample, value = cladeReads, fill = 0)

lineages <- krakenAnalytical$taxLineage

krakenAnalytical <- krakenAnalytical %>%
  select(-taxLineage) %>%
  as.matrix(.)

row.names(krakenAnalytical) <- lineages

# Read AMR and MegaBio Coverage Sampler Results -------------------------------

# Parse the results with Python script first

# Then collect the names of the parsed files

# Make sure the fecal composite data is present as well.

amrCovSamplerPaths <- Sys.glob(file.path(".",
                                    "AMR_CovSampler_parsed",
                                    '*CovSampler_parsed.tab'))

amrCovSamplerNames <- list.files(path = Sys.glob("./AMR_CovSampler_parsed/"),
                            pattern = "*CovSampler_parsed.tab")

amrCovSamplerNames <- amrCovSamplerNames %>%
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

megaBioPaths <- Sys.glob(file.path(".",
                                    "MegaBio_results",
                                    "*",
                                    '*CovSampler_parsed.tab'))

megaBioNames <- list.files(path = Sys.glob("./MegaBio_results/*/"),
                            pattern = "*CovSampler_parsed.tab")

megaBioNames <- megaBioNames %>%
  map(function(x) str_replace(x, "_MBio_CovSampler_parsed\\.tab$", ""))

megaBioReports <- megaBioPaths %>%
  map(function(x) read_tsv(x)) %>% 
      set_names(nm=megaBioNames)

megaBioReportsMerged <- megaBioReports %>%
  map_dfr(function(x) x, .id="Sample")


# Concatenate MEGARes and MEGABio data ------------------------------------

amrBioConcat <- rbind(amrReportsMerged, megaBioReportsMerged)

amrBioAnalytical <- amrBioConcat %>%
  select(Sample, Header, Hits) %>%
  spread(key = Sample, value = Hits, fill = 0)

amrBioClassification <- amrBioAnalytical$Header

amrBioAnalytical <- amrBioAnalytical %>%
  select(-Header) %>%
  as.matrix(.)

row.names(amrBioAnalytical) <- amrBioClassification
