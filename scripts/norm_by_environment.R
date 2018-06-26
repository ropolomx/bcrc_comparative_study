
# Load packages -----------------------------------------------------------

library(drake)

# Custom functions --------------------------------------------------------

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
    newMRexperiment(.) %>%
    cumNorm(.)
  normalized
}

merge_amr <- function(amr_norm){
  amr_norm$header <- row.names(amr_df)
  amr_norm <- left_join(amr_norm, annotations, by="header") # left outer join with dplyr
  amr_norm
}

group_by_amr_level <- function(amr_analytic, amr_level) {
  all_levels <- c("header", "mechanism", "group", "class")
  to_ignore <- all_levels[all_levels != amr_level]
  amr_dt <- as.data.table(amr_analytic)
  if(amr_level == "gene") {
    amr_by_level_analytic <- newMRexperiment(counts = amr_dt[!(group %in% snp_regex), .SD, .SDcols = !all_levels])
    rownames(amr_by_level_analytic) <- amr_dt$header
  } else{
  amr_by_level <- amr_dt[, lapply(.SD, sum), by = amr_level, .SDcols = !to_ignore]
  amr_by_level_analytic <- newMRexperiment(counts = amr_by_level[, .SD, .SDcols = !amr_level])
  rownames(amr_by_level_analytic) <- amr_by_level[[amr_level]]
  }
  amr_by_level_analytic
}

dt_metadata <- function(meta){
  meta_dt <- data.table(meta)
  setkeyv(meta_dt, sample_column_id)
  meta_dt
}

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

count_barplot <-function(melted, level){
  ggplot(data = subset(melted, Level_ID == "Class"), 
         aes_string(x='Type', y='Normalized_Count', fill = 'Name')) +
    geom_bar(stat="identity") + 
    scale_fill_brewer(palette = "Set3")
}


# Drake plan for normalization by environment -----------------------------

# AMR files

amr_filepath <-
  here('aggregated_data_for_analysis', 'amrBioAnalytical.csv')

amr_df <-
  read.csv(file = amr_filepath,
    header = TRUE,
    row.names = 1)

# List of Kraken files 

kraken_df <- temp_kraken_list

kraken_taxon_df <- kraken_df$taxonReads
kraken_clade_df <- kraken_df$cladeReads

# Metadata

meta_dt <- dt_metadata(metadata)

# Drake plan

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
  retranspose_amr = {
    map(
      analytic_by_environment_amr,
      ~ df_retrans(.x)
    )
  },
  retranspose_kraken = {
    modify_depth(
      analytic_by_environment_kraken,
      .depth = 2,
      ~ df_retrans(.x)
    )
  },
  normalize_by_environment_amr = {
    map(
      retranspose_amr,
      ~ normalize_split(.x)
    )
    },
  normalize_by_environment_kraken = {
    modify_depth(
      retranspose_kraken,
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
  generate_analytic_norm_amr = {
    map(
      extract_norm_amr,
      ~ merge_amr(.x)
    )
  },
  generate_analytic_norm_kraken = {
    modify_depth(
      extract_norm_kraken,
      .depth = 2,
      ~ merge_amr(.x) # Change function to tax here
    )
  },
  generate_analytic_raw_amr = {
    map(
      extract_raw_amr,
      ~ merge_amr(.x)
    )
  },
  remove_wild_type_norm = {
    map(
      generate_analytic_norm_amr,
      ~ .x[!(.x$group %in% snp_regex), ]
    )
  },
  remove_wild_type_raw = {
    map(
      generate_analytic_raw_amr,
      ~ .x[!(.x$group %in% snp_regex),] # Hack: needed to work with different notation than data table's
    )
  },
  group_by_class_norm = {
    map(
      remove_wild_type_norm,
      ~ group_by_amr_level(.x, "class")
    )
  },
  group_by_mech_norm = {
    map(
      remove_wild_type_norm,
      ~ group_by_amr_level(.x, "mechanism")
    )
  },
  group_by_group_norm = {
    map(
      remove_wild_type_norm,
      ~ group_by_amr_level(.x, "group")
    )
  },
  group_by_gene_norm = {
    map(
      remove_wild_type_norm,
      ~ group_by_amr_level(.x, "gene")
    )
  },
  group_by_class_raw = {
    map(
      remove_wild_type_raw,
      ~ group_by_amr_level(.x, "class")
    )
  },
  group_by_mech_raw = {
    map(
      remove_wild_type_raw,
      ~ group_by_amr_level(.x, "mechanism")
    )
  },
  group_by_group_raw = {
    map(
      remove_wild_type_raw,
      ~ group_by_amr_level(.x, "group")
    )
  },
  group_by_gene_raw = {
    map(
      remove_wild_type_raw,
      ~ group_by_amr_level(.x, "gene")
    )
  },
  amr_norm_list = {
    list(
      "Class" = group_by_class_norm,
      "Mechanism" = group_by_mech_norm,
      "Group" = group_by_group_norm,
      "Gene" = group_by_gene_norm
    )
  },
  amr_raw_list = {
    list(
      "Class" = group_by_class_raw,
      "Mechanism" = group_by_mech_raw,
      "Group" = group_by_group_raw,
      "Gene" = group_by_gene_raw
    )
  },
  amr_trans_norm_list = {
    purrr::transpose(amr_norm_list)
  },
  amr_trans_raw_list = {
    purrr::transpose(amr_raw_list)
  },
  melt_norm_amr = {
    map(
      amr_trans_norm_list,
      ~ rbind(
        melt_dt(MRcounts(.x$Class), "Class"),
        melt_dt(MRcounts(.x$Mechanism), "Mechanism"),
        melt_dt(MRcounts(.x$Group), "Group"),
        melt_dt(MRcounts(.x$Gene), "Gene")
      )
    )
  },
  melt_raw_amr = {
    map(
      amr_trans_raw_list,
      ~ rbind(
        melt_dt(MRcounts(.x$Class), "Class"),
        melt_dt(MRcounts(.x$Mechanism), "Mechanism"),
        melt_dt(MRcounts(.x$Group), "Group"),
        melt_dt(MRcounts(.x$Gene), "Gene")
      )
    )
  },
  match_meta_norm_amr = {
    modify_depth(
      as.vector(amr_norm_list),
        .depth = 2,
        ~ match_metadata(.x, meta_dt)
      )
  },
  match_meta_raw_amr = {
    modify_depth(
      as.vector(amr_raw_list),
      .depth = 2,
      ~ match_metadata(.x, meta_dt)
    )
  },
  plot_norm_class_amr ={
    map_dfr(
      melt_norm_amr,
      ~ .x,
      .id = "Type"
      ) %>%
  count_barplot(., "Class")
  },
  #   map2(
  #     meg_barplot(melted_data = melt_norm_amr,
  #     metadata = meta_dt,
  #     sample_var = sample_column_id,
  #     group_var = ,
  #     level_var = )
  # },
  
  # group_by_amr_levels_raw = {
  #   cross2(
  #     remove_wild_type_raw,
  #     all_levels) %>%
  #     map(~ group_by_amr_level(.x[[1]], .x[[2]]))
  # },
  strings_in_dots = "literals"
)
    
# amr_raw <- amr_raw[!(group %in% snp_regex), ]

# by_environment_eval <- evaluate_plan( 
#   by_environment_plan, 
#   rules = list(df = c(amr_df, kraken_taxon_df, kraken_clade_df)), 
#   expand = TRUE
#   )

check_plan(by_env_plan)

by_env_config <- drake_config(by_env_plan)

vis_drake_graph(by_env_config,
  from = c("kraken_df", "amr_df"),
  # to=c("remove_wild_type_raw"),
  font_size = 12)

make(by_env_plan, jobs = 2, verbose = 1)
