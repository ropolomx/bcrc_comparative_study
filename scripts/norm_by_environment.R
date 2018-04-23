
# Functions for splitting and normalizing microbiome and resistome --------

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
  merge_meta <- left_join(x, metadata, by = "ID")
  split_meta <- split(merge_meta, merge_meta$Type)
  split_meta
  }

df_retrans <- function(x){
  row.names(x) <- x$ID
  retrans <- x %>%
    select_if(is.numeric)
  retrans <- as.data.frame(t(retrans))
  retrans
}

normalize_split <- function(df){
  normalized <- 
    df %>%
    df_retrans(.) %>% 
    newMRexperiment(.) %>%
    cumNorm(.)
  normalized
}
   
# Drake plan for normalization by environment -----------------------------

# AMR files

amr_filepath <- here('aggregated_data_for_analysis', 'amrBioAnalytical.csv')
amr_df <- read.csv(file = amr_filepath, header = TRUE, row.names = 1)

# List of Kraken files 

kraken_df <- temp_kraken_list

kraken_taxon_df <- kraken_df$taxonReads
kraken_clade_df <- kraken_df$cladeReads

by_env_plan <- drake_plan(
  transpose_analytic_amr = {
    transp_df(amr_df)
    },
  transpose_analytic_kraken = {
    map(
      kraken_df, 
      ~ transp_df(.x)
    )
    },
  analytic_by_environment_amr = {
    split_by_environment(transpose_analytic_amr)
    },
  analytic_by_environment_kraken = {
    map(
      transpose_analytic_kraken,
      ~ split_by_environment(.x)
    )
    },
  normalize_by_environment_amr = {
    map(
      analytic_by_environment_amr,
      ~ normalize_split(.x)
    )
    },
  normalize_by_environment_kraken = {
    modify_depth(
      analytic_by_environment_kraken,
      .depth = 2,
      ~ normalize_split(.x)
    )
    },
  extract_norm_amr = {
    map(
      normalize_by_environment_amr,
      ~ data.table(MRcounts(.x, norm = TRUE)))
  },
  extract_norm_kraken = {
    modify_depth(
      normalize_by_environment_kraken,
      .depth = 2,
      ~ data.table(MRcounts(.x, norm = TRUE))
    )
  },
  extract_raw_amr = {
    map(
      normalize_by_environment_amr,
      ~ data.table(MRcounts(.x, norm = FALSE)))
  },
  extract_raw_kraken = {
    modify_depth(
      normalize_by_environment_kraken,
      .depth = 2,
      ~ data.table(MRcounts(.x, norm = FALSE)))
    },
  generate_analytic = {
    map(
      extract_norm_amr,
      ~ .x[, lapply(.SD, sum), by=.x$class, .SDcols=!c(.x$header, .x$mechanism, .x$group)]
    )
  },
  untidy_norm_amr = {
    map(
      generate_analytic,
      ~ widen_amr(.x)
    )
  },
  diversity_norm_amr = {
    map(
      untidy_norm_amr,
      ~ calc_diversity_df(.x)
    )
  },
  strings_in_dots = "literals"
)

# by_environment_eval <- evaluate_plan(
#   by_environment_plan, 
#   rules = list(df = c(amr_df, kraken_taxon_df, kraken_clade_df)), 
#   expand = TRUE
#   )

check_plan(by_env_plan)

by_env_config <- drake_config(by_env_plan)

vis_drake_graph(
  by_env_config, 
  from=c("kraken_df", "amr_df"), 
  to=c("diversity_norm_amr"), 
  font_size = 12
  )

make(by_env_plan, jobs = 2, verbose = 2)
