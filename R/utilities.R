# example of an R file with a helper function

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
  mat <- sapply(levels(grouping), function(gr) {
    sparseMatrixStats::rowSums2(x = counts, cols = grouping == gr)
  })
  colnames(mat) <- levels(grouping)
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
