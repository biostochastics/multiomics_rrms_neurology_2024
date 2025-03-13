###############################################################################
# COMPLETE MULTI-OMIC ANALYSIS OF OCRELIZUMAB EFFECTS IN RRMS
#  
# This script implements the analysis pipeline described in:
# “Multi-Omic characterization of the effects of Ocrelizumab in patients with 
#  relapsing-remitting multiple sclerosis”
###############################################################################


###############################################################################
# 1. SETUP AND CONFIGURATION
###############################################################################

# List of required packages
required_packages <- c(
  "data.table",      # Fast I/O
  "tidyverse",       # Data manipulation / piping
  "lme4",            # Linear mixed effects models
  "lmerTest",        # p-values for mixed models
  "broom.mixed",     # Tidy model results
  "ggplot2",         # Main plotting framework
  "ggpubr",          # Publication-ready plots
  "igraph",          # Graph analysis
  "ggraph",          # Graph visualization
  "WGCNA",           # Weighted gene co-expression network analysis
  "MEGENA",          # Multi-scale Embedded Gene Co-expression Network Analysis
  "Rfast",           # Fast row/col operations
  "corrplot",        # Correlation visualization
  "pheatmap",        # Heatmaps
  "tidymodels",      # Tidy machine learning pipeline
  "missForest",      # RF-based missing data imputation
  "parallel",        # Parallel computing
  "doParallel",      # Parallel backend
  "Hmisc",           # rcorr and other utility functions
  "gridExtra",       # Combine multiple plots
  "OlinkAnalyze",    # Olink data analysis (if preprocessing)
  "BiocManager"      # Bioconductor package manager
)

# Install any missing CRAN packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only=TRUE)
  }
}

# Some specialized packages may need BiocManager
bio_packages <- c("TissueEnrich", "BayesPrism")
for(bp in bio_packages) {
  if(!require(bp, character.only = TRUE)) {
    BiocManager::install(bp, ask = FALSE)
    library(bp, character.only=TRUE)
  }
}

# Load libraries
sapply(required_packages, require, character.only = TRUE)

# Define “not in” helper
`%ni%` <- Negate(`%in%`)

# Set up parallel cluster
num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)


###############################################################################
# 2. DATA IMPORT AND PREPROCESSING
###############################################################################

#--- 2.1 Read and preprocess proteomics data (Olink or otherwise) ---#
# file_path: path to proteomics data (long format)
# qc_threshold: min proportion of above-LOD or non-NA needed to keep a feature
read_proteomics <- function(file_path, qc_threshold = 0.5) {
  prot_dt <- data.table::fread(file_path, sep = "\t", header = TRUE)
  
  # Basic filtering: keep only measurements that are not NA
  # Then keep biomarkers that pass the qc_threshold
  prot_filtered <- prot_dt %>%
    dplyr::group_by(Biomarker) %>%
    dplyr::filter(mean(!is.na(value)) >= qc_threshold) %>%
    dplyr::ungroup()
  
  return(prot_filtered)
}


#--- 2.2 Read and preprocess metabolomics data (untargeted) ---#

read_metabolomics <- function(file_path, qc_threshold = 0.5) {
  met_dt <- data.table::fread(file_path, sep = "\t", header = TRUE)
  
  # Subset to keep only features that have at least qc_threshold 
  # proportion of non-missing.
  keep_metabolites <- met_dt %>%
    dplyr::group_by(Biomarker) %>%
    dplyr::summarise(nonNA = mean(!is.na(value))) %>%
    dplyr::filter(nonNA >= qc_threshold) %>%
    dplyr::pull(Biomarker)
  
  met_dt <- met_dt %>% dplyr::filter(Biomarker %in% keep_metabolites)
  
  # Because these are in long format, for imputation we pivot wider temporarily
  met_wide <- met_dt %>%
    tidyr::pivot_wider(
      id_cols = c("public_client_id","Time"),
      names_from = "Biomarker",
      values_from = "value"
    )
  
  # Use missForest if missing data remain
  mat_for_impute <- as.data.frame(met_wide %>% dplyr::select(-public_client_id, -Time))
  imputed <- missForest(mat_for_impute)$ximp
  
  # Put back into long format
  imputed_long <- imputed %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="rowid") %>%
    tidyr::pivot_longer(
      cols = -rowid,
      names_to = "Biomarker",
      values_to = "value"
    )
  
  # match rowid to (public_client_id, Time)
  rowid_map <- data.frame(
    rowid = as.character(seq_len(nrow(met_wide))),
    public_client_id = met_wide$public_client_id,
    Time = met_wide$Time
  )
  
  merged_long <- dplyr::left_join(imputed_long, rowid_map, by = "rowid") %>%
    dplyr::select(public_client_id, Time, Biomarker, value)
  
  # log2 transform final values if needed
  merged_long <- merged_long %>%
    dplyr::mutate(value = log2(value + 1))
  
  return(merged_long)
}


#--- 2.3 Read and preprocess lipidomics data (targeted or semi-targeted) ---#
read_lipidomics <- function(file_path, qc_threshold = 0.5) {
  lip_dt <- data.table::fread(file_path, sep = "\t", header = TRUE)
  
  # Filter out lipid-class totals if needed, or keep only actual species
  lip_filtered <- lip_dt %>%
    dplyr::group_by(Biomarker) %>%
    dplyr::filter(mean(!is.na(value)) >= qc_threshold) %>%
    dplyr::ungroup()
  
  # Log2 transform
  lip_filtered <- lip_filtered %>%
    dplyr::mutate(value = log2(value + 1))
  
  return(lip_filtered)
}


###############################################################################
# 3. UTILITY FUNCTIONS FOR WIDE FORMAT, MERGING, ETC.
###############################################################################

# Pivots a long-format table of (ID, Time, Biomarker, value) into wide format
create_wide_data <- function(data, id_cols, value_col) {
  wide_data <- data %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(id_cols),
      names_from = Biomarker,
      values_from = dplyr::all_of(value_col)
    )
  return(wide_data)
}

# Merges proteomics, metabolomics, lipidomics into a single wide table
merge_multiomics <- function(prot_data, met_data, lip_data) {
  prot_wide <- create_wide_data(prot_data, c("public_client_id", "Time"), "value")
  met_wide  <- create_wide_data(met_data,  c("public_client_id", "Time"), "value")
  lip_wide  <- create_wide_data(lip_data,  c("public_client_id", "Time"), "value")
  
  # Merge by subject ID and Time
  merged <- prot_wide %>%
    dplyr::full_join(met_wide, by = c("public_client_id", "Time")) %>%
    dplyr::full_join(lip_wide, by = c("public_client_id", "Time"))
  
  return(merged)
}


###############################################################################
# 4. LINEAR MIXED MODELS FOR EACH BIOMARKER
###############################################################################

# Fits both a continuous-time and categorical-time LMM for each biomarker.
# It returns a nested list of results including p-values, etc.
run_lmm_for_biomarker <- function(biomarker_data, time_cont_col = "Time", 
                                  time_cat_col = "Time", subject_col = "public_client_id") {
  
  # Fit continuous-time model
  mod_cont <- lmerTest::lmer(
    value ~ as.numeric(get(time_cont_col)) + (1|get(subject_col)),
    data = biomarker_data, REML = TRUE
  )
  
  # Fit categorical-time model
  biomarker_data[[time_cat_col]] <- factor(biomarker_data[[time_cat_col]])  # ensure factor
  mod_cat <- lmerTest::lmer(
    value ~ get(time_cat_col) + (1|get(subject_col)),
    data = biomarker_data, REML = TRUE
  )
  
  # Summaries
  cont_summary <- broom.mixed::tidy(mod_cont)
  cat_summary <- broom.mixed::tidy(mod_cat)
  
  # Compare models by anova
  model_comparison <- anova(mod_cont, mod_cat)
  
  return(list(
    model_continuous  = mod_cont,
    model_categorical = mod_cat,
    summary_continuous = cont_summary,
    summary_categorical = cat_summary,
    anova_comparison = model_comparison
  ))
}


# Wrapper that applies run_lmm_for_biomarker to each biomarker in a data.frame
# that must be in long format with columns: public_client_id, Time, Biomarker, value
run_lmm_all_biomarkers <- function(long_data) {
  
  results_list <- list()
  biomarkers <- unique(long_data$Biomarker)
  
  # Parallel loop over biomarkers
  results_list <- foreach::foreach(bm = biomarkers, .combine = c, .packages = c("dplyr","lme4","lmerTest","broom.mixed")) %dopar% {
    sub_data <- long_data %>% dplyr::filter(Biomarker == bm)
    if(nrow(sub_data) < 3) {
      # Not enough data to do anything
      return(list(bm = NULL))
    }
    fit <- NULL
    try({
      fit <- run_lmm_for_biomarker(sub_data)
    }, silent = TRUE)
    
    if(is.null(fit)) {
      return(list(bm = NULL))
    } else {
      return(list(bm = fit))
    }
  }
  names(results_list) <- biomarkers
  return(results_list)
}


###############################################################################
# 5. MEGENA ANALYSIS (MULTI-SCALE EMBEDDED GENE CO-EXPRESSION NETWORK ANALYSIS)
###############################################################################

run_megena_analysis <- function(expression_df, corr_method = "spearman", 
                                n_perm = 1000, pval_cutoff = 0.05, 
                                hub_pval = 0.05, min_size = 10) {
  
  # expression_df: numeric dataframe of dimension samples x features
  # Convert to a matrix
  dat <- as.matrix(expression_df)
  
  # Step 1: correlation
  # doPerm=TRUE uses a permutation-based approach to estimate FDR
  # For large feature sets, this can be slow; see MEGENA docs for scaling hints.
  ijw <- MEGENA::calculate.correlation(
    dat,
    doPerm = n_perm,
    output.corTable = TRUE,
    output.permFDR = TRUE,
    method = corr_method
  )
  
  # Step 2: Construct Planar Filtered Network (PFN)
  pfn <- MEGENA::calculate.PFN(
    ijw[, 1:3],
    doPar = TRUE,
    num.cores = num_cores
  )
  
  # Step 3: Build igraph object
  g <- igraph::graph.data.frame(pfn, directed=FALSE)
  
  # Step 4: MEGENA
  megena_result <- MEGENA::do.MEGENA(
    g,
    mod.pval = pval_cutoff,
    hub.pval = hub_pval,
    remove.unsig = TRUE,
    min.size = min_size,
    max.size = igraph::vcount(g)/2,
    doPar = TRUE,
    num.cores = num_cores,
    n.perm = 100
  )
  
  return(megena_result)
}


###############################################################################
# 6. WGCNA ANALYSIS
###############################################################################

run_wgcna_analysis <- function(expression_df) {
  # expression_df: numeric data frame or matrix (samples x features)
  # Typically, you transpose to features x samples for WGCNA.
  
  dat <- t(expression_df)  # WGCNA wants genes in rows, samples in columns
  dat <- as.data.frame(dat)
  
  # Choose soft-threshold
  powers <- 2:20
  sft <- WGCNA::pickSoftThreshold(dat, powerVector = powers, verbose = 0)
  # pick the power around which scale-free topology ~ 0.9 (or your chosen threshold)
  chosen_power <- sft$powerEstimate
  if(is.na(chosen_power)) {
    # fallback
    chosen_power <- 6
  }
  
  # blockwiseModules
  netwk <- WGCNA::blockwiseModules(
    dat,
    power = chosen_power,
    networkType = "signed",
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    saveTOMs = FALSE
  )
  
  # Convert numeric labels -> color labels
  mergedColors <- WGCNA::labels2colors(netwk$colors)
  
  # Eigengenes
  MEs <- WGCNA::moduleEigengenes(dat, colors = netwk$colors)$eigengenes
  
  return(list(
    sft = sft,
    chosen_power = chosen_power,
    network = netwk,
    module_colors = mergedColors,
    MEs = MEs
  ))
}


###############################################################################
# 7. DIGITAL CYTOMETRY (OPTIONAL) WITH BAYESPRISM
###############################################################################

# Example digital cytometry with BayesPrism (requires scRNA-seq reference)
run_digital_cytometry <- function(bulk_expression_df, sc_reference, cell_types) {
  # sc_reference is a matrix or list-based reference with known cell_types
  # cell_types is the list of cell types to deconvolve
  
  # Convert the bulk expression to matrix
  bulk_mat <- as.matrix(bulk_expression_df)
  
  # Construct BayesPrism object 
  # (see BayesPrism docs for constructing the reference argument properly)
  bp_obj <- new.prism(
    reference = sc_reference,
    mixture   = bulk_mat,
    cell.types = cell_types
  )
  
  # Run inference
  bp_obj <- run.prism(
    bp = bp_obj, 
    n.cores = num_cores
  )
  
  # Extract proportions
  props <- get.proportion(bp_obj, which.theta = "final")
  
  return(list(
    bp_obj = bp_obj,
    proportions = props
  ))
}


###############################################################################
# 8. VISUALIZATION HELPERS: Volcano, Trajectory, Heatmap, Network
###############################################################################

#--- 8.1 Volcano plot for a generic LMM result table (estimate vs p-value) ---#
create_volcano_plot <- function(df, fc_col = "estimate", pval_col = "p.value",
                                label_col = "biomarker", title = "Volcano Plot") {
  # df is a tidy model result with columns = estimate, p.value, biomarker
  df <- df %>%
    dplyr::mutate(logp = -log10(.data[[pval_col]]),
                  direction = ifelse(.data[[fc_col]] > 0, "Up", "Down"))
  
  ggplot(df, aes(x = .data[[fc_col]], y = logp)) +
    geom_point(aes(color = direction), alpha=0.6) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey40") +
    geom_vline(xintercept = 0, linetype="dashed", color="grey40") +
    scale_color_manual(values = c("Up"="red","Down"="blue")) +
    theme_bw(base_size = 14) +
    labs(x="Coefficient (Time)", y="-log10(p-value)", title=title)
}


#--- 8.2 Simple trajectory plot for a single biomarker ---#
plot_biomarker_trajectory <- function(long_data, biomarker_name,
                                      time_col="Time", id_col="public_client_id",
                                      title=NULL) {
  df <- long_data %>% dplyr::filter(Biomarker == biomarker_name)
  if(is.null(title)){
    title <- paste("Trajectory:", biomarker_name)
  }
  p <- ggplot(df, aes_string(x = time_col, y="value", group = id_col)) +
    geom_line(alpha=0.3) +
    geom_point() +
    theme_bw(base_size = 14) +
    labs(title=title)
  return(p)
}


###############################################################################
# 9. MAIN ANALYSIS PIPELINE
###############################################################################

run_analysis_pipeline <- function(
  proteomics_file,
  metabolomics_file,
  lipidomics_file,
  metadata_file,
  output_dir = "analysis_output"
) {
  
  message("==== Starting Multi-Omic Analysis Pipeline ====")
  # Create output dirs
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
  fig_dir  <- file.path(output_dir, "figures")
  res_dir  <- file.path(output_dir, "results")
  if(!dir.exists(fig_dir)) dir.create(fig_dir)
  if(!dir.exists(res_dir)) dir.create(res_dir)
  
  # 1. Load metadata (must include columns: public_client_id, Time, etc.)
  message("1) Loading metadata...")
  meta <- data.table::fread(metadata_file)
  
  # 2. Read each omic dataset
  message("2) Reading omics data...")
  prot_data <- read_proteomics(proteomics_file, qc_threshold=0.5)
  met_data  <- read_metabolomics(metabolomics_file, qc_threshold=0.5)
  lip_data  <- read_lipidomics(lipidomics_file, qc_threshold=0.5)
  
  # Optionally filter Time points only to BL, 6M, 12M if needed:
  # prot_data  <- prot_data  %>% dplyr::filter(Time %in% c(0,6,12))
  # met_data   <- met_data   %>% dplyr::filter(Time %in% c(0,6,12))
  # lip_data   <- lip_data   %>% dplyr::filter(Time %in% c(0,6,12))
  
  # 3. Run LMM for each biomarker in each dataset
  message("3) Running LMM across biomarkers...")
  all_omics_long <- dplyr::bind_rows(prot_data, met_data, lip_data)
  
  # Make sure Time is numeric (if you have 0,6,12) or 0,1,2, etc.
  # Then run
  lmm_results <- run_lmm_all_biomarkers(all_omics_long)
  
  # Summaries for continuous-time effect
  continuous_summary <- purrr::map_dfr(
    names(lmm_results),
    function(x) {
      fit <- lmm_results[[x]]
      if(is.null(fit)) return(NULL)
      cres <- fit$summary_continuous
      # keep only row with term "as.numeric(get(time_cont_col))"
      row_c <- cres %>% dplyr::filter(stringr::str_detect(term, "as.numeric"))
      if(nrow(row_c)==0) return(NULL)
      row_c$biomarker <- x
      return(row_c)
    }
  )
  
  # 4. Merge multi-omics wide for potential correlation / network
  message("4) Merging multi-omics data for network analysis...")
  # We pivot to wide so that each row = sample, columns = biomarkers
  # We must handle that each subject may have multiple timepoints
  # For an example: we’ll do it just for BL/6M/12M and keep them separate or combine. 
  wide_merged <- merge_multiomics(prot_data, met_data, lip_data)
  
  # 5. Run MEGENA
  message("5) Running MEGENA (co-expression) - be mindful of time/size complexity...")
  # For MEGENA, we typically want each row = sample, each col = biomarker
  # Then we do correlation across biomarkers. 
  # Possibly restrict to one time point or to changes. 
  # Here, for demonstration, we just use all as is, after removing ID columns.
  
  # Filter out the ID columns
  megena_input <- wide_merged %>%
    dplyr::select(-public_client_id, -Time) %>%
    # Some samples might have missing values across some biomarkers
    # You can do a quick rowwise removal or imputation:
    mutate_all(~ ifelse(is.na(.), mean(., na.rm=TRUE), .))  # simple mean impute
  
  # Scale
  megena_input_scaled <- scale(megena_input)
  
  megena_res <- run_megena_analysis(
    expression_df = as.data.frame(megena_input_scaled),
    corr_method   = "spearman"
  )
  
  # 6. Run WGCNA (optional)
  message("6) Running WGCNA (optional, can be large for many features)...")
  wgcna_res <- run_wgcna_analysis(as.data.frame(megena_input_scaled))
  
  # 7. Basic visuals: Volcano plot of time coefficient, for instance
  message("7) Generating example volcano plot for continuous-time analysis...")
  vol_plot <- create_volcano_plot(
    df = continuous_summary,
    fc_col = "estimate",
    pval_col = "p.value",
    label_col = "biomarker",
    title = "Volcano Plot - Continuous Time Coefficient"
  )
  
  # Save the figure
  ggsave(filename = file.path(fig_dir, "volcano_plot_time_continuous.pdf"),
         plot = vol_plot, width=8, height=6)
  
  # 8. Save results to disk
  message("8) Saving results...")
  # LMM results as an RDS
  saveRDS(lmm_results, file = file.path(res_dir, "lmm_results.rds"))
  # Summaries
  data.table::fwrite(
    continuous_summary,
    file = file.path(res_dir,"continuous_time_summary.tsv"),
    sep="\t"
  )
  # MEGENA
  saveRDS(megena_res, file = file.path(res_dir,"megena_results.rds"))
  # WGCNA
  saveRDS(wgcna_res, file = file.path(res_dir,"wgcna_results.rds"))
  
  message("==== Pipeline completed successfully! ====")
  
