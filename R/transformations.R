# this script defines all the count data transformations

# the methods should accept the raw count matrix (columns are cells)
#
# return value should the transformed expression matrix

## main function to delegate by transformation name

run_tr <- function(mat, transformation, seed = NULL) {
  fn <- paste0('tr_', transformation)
  if (!exists(fn, mode='function')) {
    stop(fn, ' is not a known function')
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cat('run_tr:\t', fn, '\t')
  stime = system.time({
    res <- do.call(what = fn, args = list(counts = mat))
  })
  timing <- data.frame(transformation = transformation, time = stime[3])
  cat('elapsed time:\t', stime[3], '\n')
  return(list(transformation = transformation, timing = timing, res = res))
}

## helper functions for shifted log transformations

size_factors <- function(counts, method) {
  if (is.null(method)) {
    return(NULL)
  }
  if (method == FALSE) {
    return(FALSE)
  }
  if (method == 'poscounts') {
    sf <- transformGamPoi:::estimate_size_factors(Y = as.matrix(counts), method = method)
  } else if (method %in% c('deconvolution', 'normed_sum')) {
    sf <- transformGamPoi:::estimate_size_factors(Y = counts, method = method)
  } else if (method == '10k') {
    if (inherits(x = counts, what = 'dgCMatrix')) {
      sf <- sparseMatrixStats::colSums2(counts) / 1e4
    } else {
      sf <- matrixStats::colSums2(counts) / 1e4
    }
  } else if (method == 'pf') {
    if (inherits(x = counts, what = 'dgCMatrix')) {
      cell_sum <- sparseMatrixStats::colSums2(counts)
    } else {
      cell_sum <- matrixStats::colSums2(counts)
    }
    sf <- cell_sum / mean(cell_sum)
  } else if (method == 'gmean') {
    if (inherits(x = counts, what = 'dgCMatrix')) {
      feats <- diff(counts@p)
      cell_det_rate <- feats / nrow(counts)
      cell_gmean <- sctransform:::row_gmean_dgcmatrix(matrix = Matrix::t(counts), eps = median(cell_det_rate))
      sf <- cell_gmean / exp(mean(log(cell_gmean)))
    } else {
      cell_det_rate <- matrixStats::colMeans2(counts > 0)
      eps <- median(cell_det_rate)
      cell_gmean <- exp(colMeans(log(counts + eps))) - eps
      sf <- cell_gmean / exp(mean(log(cell_gmean)))
    }
  } else {
    stop('unknown method')
  }
  return(sf)
}

shifted_log_trans  <- function(counts, pseudocount = 1, size_factor_method = NULL) {
  stopifnot(inherits(x = counts, what = 'dgCMatrix'))
  stopifnot(pseudocount > 0)
  
  # set up the return value
  tr_counts <- counts
  
  # no size factors needed, just add pseudocount 
  if (is.null(size_factor_method)) {
    tr_counts@x <- log(tr_counts@x + pseudocount)
  } else {
    if (size_factor_method == 'scran') {
      size_factor_method <- 'deconvolution'
    }
    sf <- size_factors(counts = counts, method = size_factor_method)
    tr_counts@x <- log(tr_counts@x / rep.int(sf, diff(counts@p)) + pseudocount)
  }
  if (pseudocount != 1) {
    # we no longer have a sparse matrix, as all zeros are converted to log(pseudocount)
    zeros <- tr_counts == 0
    tr_counts <- as.matrix(tr_counts)
    tr_counts[as.vector(zeros)] <- log(pseudocount)
  }
  
  return(tr_counts)
}

## individual transformations below

# counts

# no transformation
tr_counts <- function(counts) {
  return(counts)
}

# sctransform corrected counts
tr_counts_sct <- function(counts) {
  vst_out <- sctransform::vst(umi = counts, method = 'glmGamPoi', min_cells = 3, 
                              residual_type = 'none', return_cell_attr = TRUE, verbosity = 0)
  sctransform::correct_counts(x = vst_out, umi = counts, verbosity = 0)
}

tr_counts_sct2 <- function(counts) {
  vst_out <- sctransform::vst(umi = counts, vst.flavor = 'v2', min_cells = 3, 
                              residual_type = 'none', return_cell_attr = TRUE, verbosity = 0)
  sctransform::correct_counts(x = vst_out, umi = counts, verbosity = 0)
}

# fractions

tr_sum1 <- function(counts) {
  if (inherits(x = counts, what = 'dgCMatrix')) {
    tr_counts <- counts
    tr_counts@x <- tr_counts@x / rep.int(sparseMatrixStats::colSums2(counts), diff(counts@p))
    return(tr_counts)
  } else {
    stop('need dgCMatrix')
  }
}

# shifted log

tr_sl_pc1 <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = NULL)
}
tr_sl_pc1_ns <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = 'normed_sum')
}
tr_sl_pc1_pc <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = 'poscounts')
}
tr_sl_pc1_dc <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = 'scran')
}
tr_sl_pc1_pf <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = 'pf')
}
tr_sl_pc1_10k <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = '10k')
}
tr_sl_pc1_gm <- function(counts) {
  shifted_log_trans(counts, pseudocount = 1, size_factor_method = 'gmean')
}

tr_sl_pc5 <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = NULL)
}
tr_sl_pc5_ns <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = 'normed_sum')
}
tr_sl_pc5_pc <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = 'poscounts')
}
tr_sl_pc5_dc <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = 'scran')
}
tr_sl_pc5_pf <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = 'pf')
}
tr_sl_pc5_10k <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = '10k')
}
tr_sl_pc5_gm <- function(counts) {
  shifted_log_trans(counts, pseudocount = 5, size_factor_method = 'gmean')
}


# scaled shifted log
# proportional fitting followed by log1p followed by proportional fitting
tr_pf_sl_pf <- function(counts) {
  if (inherits(x = counts, what = 'dgCMatrix')) {
    tr_counts <- shifted_log_trans(counts, pseudocount = 1, size_factor_method = 'pf')
    sf <- size_factors(counts = tr_counts, method = 'pf')
    tr_counts@x <- tr_counts@x / rep.int(sf, diff(tr_counts@p))
    return(tr_counts)
  } else {
    stop('need dgCMatrix')
  }
}

# acosh via transformGamPoi

tr_acosh <- function(counts) {
  transformGamPoi::acosh_transform(data = counts, size_factors = FALSE)
}

tr_acosh_ns <- function(counts) {
  transformGamPoi::acosh_transform(data = counts, size_factors = 'normed_sum')
}

tr_acosh_pc <- function(counts) {
  res <- transformGamPoi::acosh_transform(data = as.matrix(counts), size_factors = 'poscounts')
  Matrix::Matrix(res)
}

tr_acosh_dc <- function(counts) {
  transformGamPoi::acosh_transform(data = counts, size_factors = 'deconvolution')
}

tr_acosh_10k <- function(counts) {
  sf <- size_factors(counts = counts, method = '10k')
  transformGamPoi::acosh_transform(data = counts, size_factors = sf)
}

tr_acosh_gm <- function(counts) {
  sf <- size_factors(counts = counts, method = 'gmean')
  transformGamPoi::acosh_transform(data = counts, size_factors = sf)
}

tr_acosh_pf <- function(counts) {
  sf <- size_factors(counts = counts, method = 'pf')
  transformGamPoi::acosh_transform(data = counts, size_factors = sf)
}

# residuals

tr_rt_rq <- function(counts) {
  transformGamPoi::residual_transform(data = as.matrix(counts), residual_type = "randomized_quantile")
}

tr_rt_pr_tgp <- function(counts) {
  transformGamPoi::residual_transform(data = as.matrix(counts), residual_type = "pearson")
}

tr_rt_pr_sct <- function(counts) {
  vst_out <- sctransform::vst(umi = counts, method = 'glmGamPoi', min_cells = 3, verbosity = 0)
  vst_out$y
}

tr_rt_pr_sct2 <- function(counts) {
  vst_out <- sctransform::vst(umi = counts, vst.flavor = 'v2', min_cells = 3, verbosity = 0)
  vst_out$y
}

# shifted log of corrected counts

tr_sl_pc1_sct <- function(counts) {
  vst_out <- sctransform::vst(umi = counts, method = 'glmGamPoi', min_cells = 3, 
                              residual_type = 'none', return_cell_attr = TRUE, verbosity = 0)
  cc <- sctransform::correct_counts(x = vst_out, umi = counts, verbosity = 0)
  tr_sl_pc1(cc)
}

tr_sl_pc1_sct2 <- function(counts) {
  vst_out <- sctransform::vst(umi = counts, vst.flavor = 'v2', min_cells = 3, 
                              residual_type = 'none', return_cell_attr = TRUE, verbosity = 0)
  cc <- sctransform::correct_counts(x = vst_out, umi = counts, verbosity = 0)
  tr_sl_pc1(cc)
}

