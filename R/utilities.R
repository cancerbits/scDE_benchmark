# define helper functions

# pretty print a time difference between two proc.time() calls
time_diff <- function(start_time, end_time = NULL) {
  if (is.null(end_time)) {
    end_time <- proc.time()
  }
  dt_cpu <- lubridate::make_difftime(num = sum(end_time[c('user.self', 'sys.self')] - start_time[c('user.self', 'sys.self')]))
  dt_elapsed <- lubridate::make_difftime(num = end_time['elapsed'] - start_time['elapsed'])

  sprintf('Elapsed time: %1.2f %s; CPU time: %1.2f %s',
          dt_elapsed, attr(x = dt_elapsed, which = 'units'),
          dt_cpu, attr(x = dt_cpu, which = 'units'))
}

# pretty print a time length
pprint_time <- function(dtime) {
  dtime <- lubridate::make_difftime(num = dtime)
  sprintf('%s time: %1.2f %s', names(dtime),
          dtime, attr(x = dtime, which = 'units'))
}


# load all the pipeline definitions from their csv files
load_pipelines <- function(config) {
  # read the pipeline definition files
  transformations <- readr::read_csv(file = file.path(config$project_root, 'pipelines', 'transformations.csv'))
  de_tests <- readr::read_csv(file = file.path(config$project_root, 'pipelines', 'de_tests.csv'))
  
  # expand all the transformations
  transformations <- tidyr::separate_rows(transformations, `size factor methods`) %>% 
    mutate(`size factor methods` = case_when(`size factor methods` == 'none' ~ NA_character_,
                                             TRUE ~ `size factor methods`),
           name = paste(`transformation method`, `size factor methods`, sep = '_'),
           name = stringr::str_remove(name, '_NA$'))
  # transformations <- tidyr::separate_rows(transformations, `size factor methods`) %>% 
  #   mutate(idx = 1:n()) %>%
  #   group_by(idx) %>%
  #   mutate(name = stringr::str_interp(string = name, env = list(sf = `size factor methods`))) %>%
  #   ungroup()
  
  # expand all the test methods
  de_tests <- tidyr::separate_rows(de_tests, `test type`) %>%
    tidyr::separate_rows(`size factor methods`) %>%
    tidyr::separate_rows(`pseudo replicates`) %>%
    tidyr::separate_rows(`aggregation strategy`) %>%
    mutate(`size factor methods` = case_when(`size factor methods` == 'none' ~ NA_character_,
                                             TRUE ~ `size factor methods`),
           name = paste(`de method`, `size factor methods`, `test type`, `pseudo replicates`, `aggregation strategy`, sep = '_'),
           name = stringr::str_remove_all(name, '_NA'))
    
  # 
  # de_tests <- tidyr::separate_rows(de_tests, `size factor methods`) %>% 
  #   tidyr::separate_rows(`pseudo replicates`) %>%
  #   mutate(idx = 1:n()) %>%
  #   group_by(idx) %>%
  #   mutate(name = stringr::str_interp(string = name, env = list(sf = `size factor methods`, pr = `pseudo replicates`))) %>%
  #   ungroup()
  
  # create data frame of pipelines
  pipelines <- cross_join(transformations, de_tests, suffix = c(".trans", ".de")) %>%
    filter(`use with transformation` == 'ANY' | 
          (`use with transformation` == 'COUNTS' & grepl(pattern = 'counts', x = name.trans)) | 
          (`use with transformation` == 'RAW_COUNTS' & name.trans == 'counts') |
          (`use with transformation` == 'SPARSE' & `result is sparse`) |
          (`use with transformation` == 'POS' & `result is pos`))
  
  # make sure the types are as expected
  pipelines <- mutate(pipelines,
                      `pseudo replicates` = case_when(`pseudo replicates` == 'NA' ~ NA_integer_,
                                                      TRUE ~ as.integer(`pseudo replicates`)))
  return(pipelines)
}


# helper for downsampling counts
downsample_counts <- function(counts, frac = 0.70) {
  if (frac >= 1) {
    return(counts)
  }
  target_umi <- round(colSums(counts) * frac)
  new_counts <- Seurat:::RunUMISamplingPerCell(data = counts, sample_val = target_umi, display_progress = FALSE)
  # Seurat:::RunUMISamplingPerCell creates non-standard sparse matrix with explicit zeros - fix here
  new_counts <- Matrix::drop0(new_counts)
  dimnames(new_counts) <- dimnames(counts)
  return(new_counts)
}

# helper to create pseudobulk from single-cell data
pseudobulk <- function(counts, grouping) {
  group_levels <- levels(grouping)
  has_group <- !is.na(grouping)
  mat <- sapply(group_levels, function(gr) {
    sparseMatrixStats::rowSums2(x = counts, cols = has_group & grouping == gr)
  })
  colnames(mat) <- group_levels
  rownames(mat) <- rownames(counts)
  return(mat)
}

# helper for adding an effect direction
# (used if a DE method does not return one)
effect_direction_via_mean <- function(mat, grouping) {
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  mat <- as.matrix(mat)
  m1 <- matrixStats::rowMeans2(x = mat, useNames = FALSE, cols = sel1)
  m2 <- matrixStats::rowMeans2(x = mat, useNames = FALSE, cols = sel2)
  return(sign(m1 / m2))
}

run_pipeline <- function(input_data, transformation, de_method,
                         test_type, size_factor_method, G, aggregation_strategy,
                         seed_trans = NULL, seed_de = NULL) {
  trans_out <- run_tr(mat = input_data$counts, transformation = transformation, seed = seed_trans)
  
  de_out <- run_de(mat = trans_out$res, 
                   grouping = input_data$cell_metadata$group_id, 
                   test_method = de_method, 
                   test_type = test_type, 
                   size_factor_method = size_factor_method, 
                   G = G, 
                   aggregation_strategy = aggregation_strategy, 
                   seed = seed_de)
  
  # join the de results with the simulation parameters if they are available
  if ('sim_metadata' %in% names(input_data)) {
    res <- dplyr::left_join(de_out$res, input_data$sim_metadata$gene_info, by = c('feature' = 'gene'))
  } else {
    res <- de_out$res
  }
  timing <- rbind(trans_out$time, de_out$time)
  rownames(timing) <- c('transformation', 'de_method')
  return(list(res = res, timing = timing))
}

run_pipeline_old <- function(input_data, transformation, de_method,
                         seed_trans = NULL, seed_de = NULL) {
  trans_out <- run_tr(mat = input_data$counts, transformation = transformation, seed = seed_trans)
  de_out <- run_de(mat = trans_out$res, grouping = input_data$cell_metadata$group_id,
                   methods = de_method, seed = seed_de)
  # join the de results with the simulation parameters if they are available
  if ('sim_metadata' %in% names(input_data)) {
    res <- dplyr::left_join(de_out$res, input_data$sim_metadata$gene_info, by = c('feature' = 'gene'))
  } else {
    res <- de_out$res
  }

  timing <- data.frame(step = c('transformation', 'de_method'),
                       name = c(transformation, de_method),
                       time = c(trans_out$timing$time, de_out$timing$time))
  return(list(res = res, timing = timing))
}

# for classification performance metrics

gilbert_skill_score <- function(TP, FN, FP, TN) {
  (TP * TN - FN * FP) / ((FN + FP) * (TP + FN + FP + TN) + (TP * TN - FN * FP))
}

jacc_sim <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

class_stats <- function(prediction, reference, FDR = NULL, JS = TRUE) {
  cm <- caret::confusionMatrix(data = prediction, reference = reference, positive = 'TRUE')
  res <- data.frame(t(c(cm$byClass, cm$overall)))
  if (!is.null(FDR)) {
    pr_curve <- PRROC::pr.curve(scores.class0 = FDR, weights.class0 = reference == 'FALSE')
    res$AUPR <- pr_curve$auc.integral
  }
  tab <- data.frame(t(as.numeric(cm$table)))
  colnames(tab) <- c('TN', 'FP', 'FN', 'TP')
  if (JS) {
    js = sum((prediction == 'TRUE') & (reference == 'TRUE')) / sum((prediction == 'TRUE') | (reference == 'TRUE'))
    tab <- cbind(tab, data.frame(JS = js))
  }
  tab <- cbind(tab, data.frame(GSS = gilbert_skill_score(tab$TP, tab$FN, tab$FP, tab$TN),
                               Bias = (tab$TP + tab$FP) / (tab$TP + tab$FN)))
  res <- cbind(tab, res)
  return(res)
}

# assumes three states: -1, 0, 1
three_way_confusion <- function(prediction, reference) {
  prediction <- factor(prediction, levels = c(-1, 0, 1))
  reference <- factor(reference, levels = c(-1, 0, 1))
  tab <- table(reference, prediction)
  OD <- tab[1,3] + tab[3,1]  # opposite effect direction
  TP <- tab[1,1] + tab[3,3]
  FN <- tab[1,2] + tab[3,2] + OD/2
  FP <- tab[2,1] + tab[2,3] + OD/2
  TN <- length(reference) - TP - FN - FP
  return(data.frame(TP=TP, FN=FN, FP=FP, TN=TN, OD=OD))
}

perf_metrics <- function(prediction, reference) {
  data.frame(n = length(reference)) %>%
    bind_cols(three_way_confusion(prediction, reference)) %>%
    mutate(Precision = TP / (TP + FP),
           Sensitivity = TP / (TP + FN),
           Specificity = TN / (TN + FP),
           Bias = (TP + FP) / (TP + FN),
           FDR = case_when(TP + FP == 0 ~ 0, TRUE ~ 1 - Precision),
           BA = (Sensitivity + Specificity) / 2,
           F1 = (2 * TP) / (2 * TP + FP + FN),
           MCC = case_when((TP + FP) == 0 | (TP + FN) == 0 | (TN + FP) == 0 | (TN + FN) == 0 ~ 0,
                           TRUE ~ ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))),
           FM = case_when((TP + FP) == 0 ~ 0, TRUE ~ sqrt(Precision * Sensitivity)),
           GSS = gilbert_skill_score(TP, FN, FP, TN),
           CSI = TP / (TP + FN + FP),
           PSS = (TP / (TP + FN)) - (FP / (FP + TN)),
           OR = (TP * TN) / (FP * FN),
           ORSS = (OR - 1) / (OR + 1))
}


## some basic scRNA-seq processing

# normalization, dimensionality reduction, hierarchical clustering
default_processing <- function(counts, sct_res_clip_range = c(-10, 10), pca_dim = 7, run_umap = FALSE) {
  vst_out <- sctransform::vst(umi = counts, method = 'glmGamPoi', res_clip_range = sct_res_clip_range, verbosity = 0)
  pca_out <- irlba::prcomp_irlba(x = t(vst_out$y), n = pca_dim, center = TRUE, scale. = FALSE)
  dmat <- dist(pca_out$x, method = 'euclidean')
  hc <- fastcluster::hclust(d = dmat, method = 'ward.D2')
  umap <- NA
  if (run_umap) {
    umap <- uwot::umap(X = dmat, n_neighbors = nrow(pca_out$x)/2, fast_sgd = TRUE, min_dist = 1)
  }
  return(list(vst_out = vst_out, pca_out = pca_out, dmat = dmat, hc = hc, umap = umap))
}

# hard clustering
cluster_counts <- function(counts, k, ...) {
  dp_out <- default_processing(counts, ...)
  membership <- as.numeric(droplevels(factor(cutree(tree = dp_out$hc, k = k))))-1
  tab <- table(membership)
  if (length(tab) != k) {
    stop('Fewer/more clusters than asked for')
  }
  if (sum(tab > 0) < k) {
    stop('Empty clusters')
  }
  if (min(tab) == 1) {
    warning('At least one cluster has only one member')
  }
  return(membership)
}

