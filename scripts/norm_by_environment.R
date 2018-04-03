# Split and normalize by environment type ----------------------------------

# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore

# Split analytic data frames

transp_df <- function(x){
  analytic <- as.data.frame(t(x))
  analytic$ID <- row.names(analytic)
  row.names(analytic) <- NULL
  analytic
}

split_by_environment <- function(x){
  merge_meta <- left_join(amr_trans, metadata, by = "ID")
  split_meta <- split(merge_meta$Type)
  return(split_meta)
}

df_retrans <- function(x){
  row.names(x) <- x$ID
  retrans <- x %>%
    select_if(is.numeric)
  retrans <- as.data.frame(t(amr_retrans))
  retrans
}

normalize_split <- function(split_df){
  normalized <- split_df %>%
    map(~ df_retrans(.x)) %>%
    map(~ newMRexperiment(.x)) %>%
    map(~ cumNorm(.x))
}
   
# Drake plan for normalization of AMR data by environment -----------------

# AMR files

amr_filepath <- here('aggregated_data_for_analysis', 'amrBioAnalytical.csv')
amr_df <- read.csv(file = amr_filepath, header = TRUE, row.names = 1)

kraken_df <- temp_kraken_list

by_environment_plan <- drake_plan(
  transpose_analytic = {
    transp_df(df)
    },
  analytic_by_environment = {
    split_by_environment(transposed)
    },
  normalize_by_environment = {
    normalize_split(analytic_by_environment)
    },
  extract_norm = {
    map(
      normalize_by_environment,
      ~ data.table(MRcounts(.x, norm = TRUE)))
  },
  extract_raw = {
    map(
      normalize_by_environment,
      ~ data.table(MRcounts(.x, norm = FALSE)))
  } 
)

by_env_config <- drake_config(by_environment_plan)
vis_drake_graph(by_env_config, targets_only = TRUE, font_size = 12)
make(by_environment_plan)
