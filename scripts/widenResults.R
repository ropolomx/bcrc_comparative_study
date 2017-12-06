library(readr)
library(stringr)
library(purrr)

# Read Kraken reports -----------------------------------------------------

# Obtain the filenames of all Kraken reports

krakenReportPaths <- Sys.glob(file.path(".",
                                    "post_PhiXFilter_Kraken",
                                    '*.tabular'))

krakenReportNames <- list.files(path = Sys.glob("./post_PhiXFilter_Kraken/"),
                            pattern = "*.tabular")

krakenReportNames <- krakenReportNames %>%
  map(function(x) str_replace(x, "_Kraken\\.tabular$", ""))

krakenReports <- krakenReportPaths %>%
  map(read_tsv)
  # set_names(nm=krakenReportNames)

# Read AMR and MegaBio Coverage Sampler Results -------------------------------

# Parse the results with Python script first

# Then collect the names of the parsed files

amrCovSamplerPaths <- Sys.glob(file.path(".",
                                    "*",
                                    '*CovSampler*.tabular'))

amrCovSamplerNames <- list.files(path = Sys.glob("./*"),
                            pattern = "*CovSampler*.tabular")

amrCovSamplerNames <- amrCovSamplerNames %>%
  map(function(x) str_replace(x, "_CovSampler\\.tabular$", ""))

amrCovSampler <- amrCovSamplerPaths %>%
  map(read_tsv) %>%
  set_names(nm=amrCovSamplerNames)

# Let's now read all the Coverage Sampler tabular files
# We are using the readr package (read_tsv)
# We are also using the list of sample names extracted in the previous function
# to set the names of the list elements
# This will make life so much easier!

# Join all the datasets into one dataframe that will be analyzed

amrResults <- do.call("rbind", amrResults)

# Need to add a column with sample name to the Coverage Sampler output
# Also need to add a column with sample depth

amrResults$SampleName <- row.names(amrResults)
amrResults$SampleID <- str_extract(amrResults$SampleName, "^.*_rarefied")
amrResults$SampleID <- str_replace(amrResults$SampleID, "_rarefied", "")
amrResults$Depth <- str_extract(amrResults$SampleName, "rarefied_.*_")
amrResults$Depth <- str_replace(amrResults$Depth, "rarefied_", "")
amrResults$Depth <- str_replace(amrResults$Depth, "_","")
amrResults$amrLevel <- str_extract(amrResults$SampleName, "(class|mechanism|group|gene)")

# The other (easier) option is to read csv file containing all rarefied datasets which were concatenated with Python-Pandas

amrRarefiedConcat <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefiedConcat.csv')

# The dataframes amrResults or amrRarefiedConcat that were generated above can be used for 
#plotting rarefaction curves with ggplot2