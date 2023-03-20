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

run_pipeline <- function(input_data, transformation, de_method,
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
