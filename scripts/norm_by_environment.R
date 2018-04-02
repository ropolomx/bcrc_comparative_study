# Split and normalize by environment type ----------------------------------

# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore

# Split AMR

transp_amr <- function(x){
  amr_analytic <- as.data.frame(t(x))
  amr_analytic$ID <- row.names(amr_analytic)
  row.names(amr_analytic) <- NULL
  amr_analytic
}

amr_df <- read.table(here('aggregated_data_for_analysis', 'amrBioAnalytical.csv'), 
                                  header=T, row.names=1, sep=',')

amr_trans <- transp_amr(amr_df)

amr_by_environment <- left_join(amr_trans, metadata, by = "ID") %>%
  split(.$Type)

df_retrans <- function(x){
  row.names(x) <- x$ID
  amr_retrans <- x %>%
    select_if(is.numeric)
  amr_retrans <- as.data.frame(t(amr_retrans))
  amr_retrans
}

amr_by_environment_norm <-
  amr_by_environment %>%
  map() %>%
  map(~ newMRexperiment(.x)) %>%
  map(~ cumNorm(.x))


# Drake plan for normalization of AMR data by environment -----------------

by_environment_plan <- drake::drake_plan(
  amr_by_environment_norm = {
    map(amr_by_environment, ~ df_retrans(.x)) %>%
      map(~ newMRexperiment(.x)) %>%
      map(~ cumNorm(.x))
    }
)
