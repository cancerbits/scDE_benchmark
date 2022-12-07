# this script defines all the DE methods

# single-cell methods should accept an expression matrix (columns are cells)
# and a grouping variable (factor)

## main function to delegate by method name

# methods may be a vector
run_de <- function(mat, grouping, methods, seed = NULL) {
  
  timing <- data.frame(method = methods, time = NA_real_)
  res_lst <- list()
  i <- 0
  for (method in methods) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    i <- i + 1
    fn <- paste0('de_', method)
    cat('run_de:\t', fn, '\t')
    stime = system.time({
      res_lst[[method]] <- do.call(what = fn, args = list(mat = mat, grouping = grouping))
    })
    res_lst[[method]]$method <- factor(method, levels = methods)
    timing$time[i] <- stime[3]
    cat('elapsed time:\t', stime[3], '\n')
  }
  res <- do.call(rbind, res_lst)
  rownames(res) <- NULL
  return(list(timing = timing, res = res))
}

# simple effect size estimation: difference of medians
effect_size <- function(mat, grouping) {
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  mat1 <- mat[, sel1, drop = FALSE]
  mat2 <- mat[, sel2, drop = FALSE]
  if (inherits(x = mat, what = 'dgCMatrix')) {
    median_diff <- sparseMatrixStats::rowMedians(mat1) - sparseMatrixStats::rowMedians(mat2)
  } else {
    median_diff <- matrixStats::rowMedians(mat1) - matrixStats::rowMedians(mat2)
  }
  return(median_diff)
}

## single-cell methods operating on expression data of any kind (any transformation)

# riffle with empirical p-value
de_random_emp <- function(mat, grouping) {
  ro <- riffle::diff_mean_test(y = mat, group_labels = grouping, 
                               R = 499, log2FC_th = -Inf, mean_th = -Inf, cells_th = 0, verbosity = 0)
  res <- select(ro, gene, emp_pval, emp_pval_adj, zscore) %>%
    mutate(zscore = sign(-zscore)) %>%
    dplyr::rename(feature = gene, pval = emp_pval, FDR = emp_pval_adj, effect_direction = zscore)
  return(res)
}

# riffle with approximated p-value
de_random_app <- function(mat, grouping) {
  ro <- riffle::diff_mean_test(y = mat, group_labels = grouping, 
                               R = 49, log2FC_th = -Inf, mean_th = -Inf, cells_th = 0, verbosity = 0)
  res <- select(ro, gene, pval, pval_adj, zscore) %>%
    mutate(zscore = sign(-zscore)) %>%
    dplyr::rename(feature = gene, FDR = pval_adj, effect_direction = zscore)
  return(res)
}

# wilcoxon rank-sum test as implemented in presto
de_wilcox <- function(mat, grouping) {
  ref_grp <- levels(grouping)[2]
  res <- presto::wilcoxauc(X = mat, y = grouping) %>%
    filter(group == ref_grp) %>%
    mutate(effect_direction = sign(auc - 0.5)) %>%
    dplyr::rename(FDR = padj) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

# wilcoxon rank-sum test as implemented in limma
de_wilcox_limma <- function(mat, grouping) {
  mat <- as.matrix(mat)
  j <- which(grouping == levels(grouping)[1])
  pval <- apply(mat, 1, function(x) {
    return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = x)), 1))
  })
  res <- data.frame(feature = rownames(mat), pval = pval) %>%
    mutate(FDR = p.adjust(p = pval, method = 'fdr'))
  return(res)
}

# wilcoxon rank-sum test as implemented in stats
de_wilcox_stats <- function(mat, grouping) {
  mat <- as.matrix(mat)
  pval <- apply(mat, 1, function(x) {
    return(wilcox.test(x ~ grouping)$p.value)
  })
  res <- data.frame(feature = rownames(mat), pval = pval) %>%
    mutate(FDR = p.adjust(p = pval, method = 'fdr'))
  return(res)
}

# many ttests; main idea and most code is taken from Rfast ttests
de_ttest <- function(mat, grouping) {
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  if (inherits(x = mat, what = 'dgCMatrix')) {
    mat1 <- mat[, sel1, drop = FALSE]
    mat2 <- mat[, sel2, drop = FALSE]
    n1 <- ncol(x = mat1)
    n2 <- ncol(x = mat2)
    m1 <- sparseMatrixStats::rowMeans2(x = mat1, useNames = FALSE)
    m2 <- sparseMatrixStats::rowMeans2(x = mat2, useNames = FALSE)
    f1 <- sparseMatrixStats::rowVars(x = mat1, center = m1, useNames = FALSE) / n1
    f2 <- sparseMatrixStats::rowVars(x = mat2, center = m2, useNames = FALSE) / n2
  } else {
    n1 <- sum(sel1)
    n2 <- sum(sel2)
    m1 <- matrixStats::rowMeans2(x = mat, useNames = FALSE, cols = sel1)
    m2 <- matrixStats::rowMeans2(x = mat, useNames = FALSE, cols = sel2)
    f1 <- matrixStats::rowVars(x = mat, center = m1, cols = sel1, useNames = FALSE) / n1
    f2 <- matrixStats::rowVars(x = mat, center = m2, cols = sel2, useNames = FALSE) / n2
  }
  fac <- f1 + f2
  dof <- fac^2 / ( f1^2 / (n1 - 1) + f2^2 / (n2 - 1) )
  stat <- ( m1 - m2 ) / sqrt(fac)
  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )
  res <- data.frame(feature = rownames(mat), pval = pvalue) %>%
    mutate(FDR = p.adjust(p = pval, method = 'fdr'), effect_direction = sign(-stat))
  return(res)
}

# logistic regression likelihood ratio test
de_logreg_slow <- function(mat, grouping) {
  mat <- as.matrix(mat)
  pval <- apply(mat, 1, function(x) {
    suppressWarnings(mod1 <- glm(grouping ~ x, family = 'binomial'))
    suppressWarnings(mod2 <- glm(grouping ~ 1, family = 'binomial'))
    suppressWarnings(this_p <- lmtest::lrtest(mod1, mod2)$Pr[2])
    return(this_p)
  })
  res <- data.frame(feature = rownames(mat), pval = pval) %>%
    mutate(FDR = p.adjust(p = pval, method = 'fdr'))
  return(res)
}

# Rfast-style logistic regression
de_logreg <- function(mat, grouping) {
  y <- as.numeric(grouping) - 1
  x <- t(as.matrix(mat))
  dm <- dim(x)
  n <- dm[1]
  d <- dm[2]
  p <- sum(y)/n
  ini <-  - 2 * ( n * p * log(p) + (n - n * p) * log(1 - p) )
  mod <- Rfast::logistic_only(x, y, b_values = TRUE)
  stat <- ini - mod['deviance', ]
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = FALSE)
  res <- data.frame(feature = rownames(mat), pval = pval) %>%
    mutate(FDR = p.adjust(p = pval, method = 'fdr'), effect_direction = sign(mod['slope', ]))
  return(res)
}

# LRT model proposed in McDavid et al, Bioinformatics, 2013; uses Seurat
de_lrt <- function(mat, grouping) {
  mat <- as.matrix(mat)
  xmat <- t(mat[, grouping == levels(grouping)[1]])
  ymat <- t(mat[, grouping == levels(grouping)[2]])
  pval <- sapply(1:nrow(mat), function(i) {
    Seurat:::DifferentialLRT(x = xmat[, i], y = ymat[, i])
  })
  res <- data.frame(feature = rownames(mat), pval = pval) %>%
    mutate(FDR = p.adjust(p = pval, method = 'fdr'))
  return(res)
}

# MAST
de_mast <- function(mat, grouping) {
  ref_grp <- levels(grouping)[2]
  contr <- paste0('group', ref_grp)
  # MAST is noisy (several messages) - shut it up
  suppressMessages({
    sca <- MAST::FromMatrix(
      exprsArray = as.matrix(mat),
      check_sanity = FALSE,
      cData = data.frame(group = grouping),
      fData = data.frame(name = rownames(mat))
    )
    zlmCond <- MAST::zlm(formula = ~ group, sca = sca)
    summaryCond <- MAST::summary(object = zlmCond, doLRT = contr)}
  )
  summaryDt <- filter(summaryCond$datatable, contrast == contr, component %in% c('D', 'H', 'logFC')) %>% 
    tidyr::pivot_wider(id_cols = primerid, names_from = component, values_from = c(`Pr(>Chisq)`, coef))
  res <- data.frame(feature = rownames(mat)) %>% 
    left_join(summaryDt, by = c('feature' = 'primerid')) %>%
    mutate(pval = `Pr(>Chisq)_H`, FDR = p.adjust(p = pval, method = 'fdr'), 
           effect_direction = case_when(is.nan(coef_logFC) ~ sign(coef_D), TRUE ~ sign(coef_logFC))) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

## single-cell methods operating on counts

# quasi-likelihood ratio testing as implemented in glmGamPoi
de_qlrt <- function(mat, grouping, size_factors = FALSE) {
  # I have not figured out how to specify the contrast if there are spaces in the group levels
  levels(grouping) <- gsub(pattern = ' ', replacement = '_', x = levels(grouping))
  ref_grp <- levels(grouping)[2]
  fit <- glmGamPoi::glm_gp(data = as.matrix(mat), design = ~ grouping, size_factors = size_factors)
  res <- glmGamPoi::test_de(fit = fit, contrast = paste0('grouping', ref_grp))
  res <- dplyr::rename(res, feature = name, FDR = adj_pval) %>%
    mutate(effect_direction = sign(lfc)) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

de_qlrt_10k <- function(mat, grouping) {
  sf <- size_factors(counts = mat, method = '10k')
  de_qlrt(mat, grouping, size_factors = sf)
}

de_qlrt_gm <- function(mat, grouping) {
  sf <- size_factors(counts = mat, method = 'gmean')
  de_qlrt(mat, grouping, size_factors = sf)
}

de_qlrt_pf <- function(mat, grouping) {
  sf <- size_factors(counts = mat, method = 'pf')
  de_qlrt(mat, grouping, size_factors = sf)
}

de_qlrt_ns <- function(mat, grouping) {
  de_qlrt(mat, grouping, size_factors = 'normed_sum')
}

de_qlrt_dc <- function(mat, grouping) {
  de_qlrt(mat, grouping, size_factors = 'deconvolution')
}

de_qlrt_pc <- function(mat, grouping) {
  de_qlrt(mat, grouping, size_factors = 'poscounts')
}

# based on http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
de_scdeseq <- function(mat, grouping) {
  mat <- as.matrix(mat)
  dds <- DESeq2::DESeqDataSetFromMatrix(mat, data.frame(grouping = grouping), ~ grouping)
  DESeq2::sizeFactors(dds) <- scran::calculateSumFactors(mat)
  dds <- DESeq2::estimateDispersions(dds, fitType = 'glmGamPoi')
  dds <- DESeq2::nbinomLRT(object = dds, reduced= ~ 1, type = 'glmGamPoi')
  
  res <- as.data.frame(DESeq2::results(dds)) %>% tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename(pval = pvalue) %>%
    mutate(FDR = p.adjust(pval, method = 'fdr'), effect_direction = sign(log2FoldChange)) %>%
    select(feature, pval, FDR, effect_direction)
}

de_scedger_lrt <- function(mat, grouping) {
  de_edger(mat, grouping, test_type = 'LRT')
}

de_scedger_qlf <- function(mat, grouping) {
  de_edger(mat, grouping, test_type = 'QLF')
}

de_scedger_exact <- function(mat, grouping) {
  de_edger(mat, grouping, test_type = 'exact')
}

# taken from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
de_scedger_qlfdetrate <- function(mat, grouping) {
  dge <- edgeR::DGEList(mat, group = grouping)
  dge <- edgeR::calcNormFactors(dge)
  cdr <- scale(sparseMatrixStats::colMeans2(mat > 0))
  design <- model.matrix(~ cdr + grouping)
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design)
  test_res <- edgeR::glmQLFTest(fit)
  res <- test_res$table %>% tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename(pval = PValue) %>%
    mutate(FDR = p.adjust(pval, method = 'fdr'), effect_direction = sign(logFC)) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}


## bulk methods applied to pseudo-bulk data

pseudobulk_de <- function(mat, grouping, G, test_method, ...) {
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  n1 <- sum(sel1)
  n2 <- sum(sel2)
  if (n1 < G | n2 < G) {
    stop('too few cells')
  }
  # create replicate labels per group
  replicate <- numeric(length = n1 + n2)
  replicate[sel1] <- sample(1:n1 %% G)
  replicate[sel2] <- sample(1:n2 %% G) + G
  replicate <- factor(replicate, levels = 0:(2*G-1))
  # create pseudobulk
  pb <- pseudobulk(counts = mat, grouping = replicate)
  # create pseudobulk grouping factor
  pb_group <- factor(1:(2*G) > G)
  
  if (test_method == 'edgeR') {
    res <- de_edger(mat = pb, grouping = pb_group, ...)
  } else if (test_method == 'DESeq2') {
    res <- de_deseq(mat = pb, grouping = pb_group, ...)
  } else if (test_method == 'glmGamPoi') {
    ellipsis_args <- list(...)
    if ('size_factors' %in% names(ellipsis_args)) {
      sf <- size_factors(counts = pb, method = list(...)$size_factors)
      res <- de_qlrt(mat = pb, grouping = pb_group, size_factors = sf)
    } else {
      stop('need size_factors argument')
    }
  } else if (test_method == 'limma') {
    res <- de_limma(mat = pb, grouping = pb_group, ...)
  } else {
    stop('test_method unknown')
  }
  return(res)
}

# edgeR 

# run edgeR on bulk or pseudobulk data
de_edger <- function(mat, grouping, test_type) {
  design <- model.matrix(~ grouping)
  
  y <- edgeR::DGEList(counts = mat, group = grouping)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  
  if (test_type == 'exact') {
    test_res <- edgeR::exactTest(y)
  }
  else if (test_type == 'LRT') {
    fit <- edgeR::glmFit(y, design)
    test_res <- edgeR::glmLRT(fit)
  }
  else if (test_type == 'QLF') {
    fit <- edgeR::glmQLFit(y, design)
    test_res <- edgeR::glmQLFTest(fit)
  }
  else {
    stop('test_type must be "exact", "LRT" or "QLF"')
  }
  
  res <- test_res$table %>% tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename(pval = PValue) %>%
    mutate(FDR = p.adjust(pval, method = 'fdr'), effect_direction = sign(logFC)) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

de_limma <- function(mat, grouping, test_type) {
  design <- model.matrix(~ grouping)
  
  allow_trend <- FALSE
  if (test_type == 'voom') {
    y <- limma::voom(counts = mat, design = design)
  }
  else if (test_type == 'trend') {
    allow_trend <- TRUE
    y <- edgeR::DGEList(counts = mat, group = grouping)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::cpm(y, log = TRUE, prior.count = 3)
  } 
  else {
    stop('test_type must be "voom" or "trend"')
  }
  
  fit <- limma::lmFit(y, design)
  fit <- limma::eBayes(fit, trend = allow_trend, robust = allow_trend)
  res <- data.frame(feature = rownames(fit$p.value),
                    pval = fit$p.value[, ncol(design)], 
                    row.names = NULL) %>% 
    mutate(FDR = p.adjust(pval, method = 'fdr'), effect_direction = sign(fit$t[, ncol(design)])) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

de_edger_exact_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'edgeR', test_type = 'exact')
}

de_edger_exact_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'edgeR', test_type = 'exact')
}

de_edger_exact_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'edgeR', test_type = 'exact')
}

de_edger_lrt_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'edgeR', test_type = 'LRT')
}

de_edger_lrt_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'edgeR', test_type = 'LRT')
}

de_edger_lrt_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'edgeR', test_type = 'LRT')
}

de_edger_qlf_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'edgeR', test_type = 'QLF')
}

de_edger_qlf_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'edgeR', test_type = 'QLF')
}

de_edger_qlf_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'edgeR', test_type = 'QLF')
}

# DESeq2

de_deseq <- function(mat, grouping, test_type) {
  dds <- DESeq2::DESeqDataSetFromMatrix(mat, data.frame(grouping = grouping), ~ grouping)
  if (test_type == 'wald') {
    dds <- DESeq2::DESeq(dds, test = 'Wald')
  } else if (test_type == 'lrt') {
    dds <- DESeq2::DESeq(dds, test = 'LRT', reduced= ~ 1)
  } else {
    stop('test_type needs to be either wald or lrt')
  }
  
  res <- as.data.frame(DESeq2::results(dds)) %>% tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename(pval = pvalue) %>%
    mutate(FDR = p.adjust(pval, method = 'fdr'), effect_direction = sign(log2FoldChange)) %>%
    select(feature, pval, FDR, effect_direction)
}

de_deseq_wald_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'DESeq2', test_type = 'wald')
}

de_deseq_wald_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'DESeq2', test_type = 'wald')
}

de_deseq_wald_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'DESeq2', test_type = 'wald')
}

de_deseq_lrt_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'DESeq2', test_type = 'lrt')
}

de_deseq_lrt_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'DESeq2', test_type = 'lrt')
}

de_deseq_lrt_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'DESeq2', test_type = 'lrt')
}

# limma
de_limma_voom_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'limma', test_type = 'voom')
}

de_limma_voom_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'limma', test_type = 'voom')
}

de_limma_voom_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'limma', test_type = 'voom')
}

de_limma_trend_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'limma', test_type = 'trend')
}

de_limma_trend_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'limma', test_type = 'trend')
}

de_limma_trend_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'limma', test_type = 'trend')
}

# use glmGamPoi qlrt on pseudobulk data (test since this could be an alternative to edger/deseq)

de_qlrt_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = FALSE)
}

de_qlrt_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = FALSE)
}

de_qlrt_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = FALSE)
}

de_qlrt_10k_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = '10k')
}

de_qlrt_10k_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = '10k')
}

de_qlrt_10k_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = '10k')
}

de_qlrt_pf_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = 'pf')
}

de_qlrt_pf_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = 'pf')
}

de_qlrt_pf_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = 'pf')
}

de_qlrt_ns_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = 'normed_sum')
}

de_qlrt_ns_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = 'normed_sum')
}

de_qlrt_ns_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = 'normed_sum')
}

de_qlrt_pc_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = 'poscounts')
}

de_qlrt_pc_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = 'poscounts')
}

de_qlrt_pc_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = 'poscounts')
}

de_qlrt_dc_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = 'deconvolution')
}

de_qlrt_dc_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = 'deconvolution')
}

de_qlrt_dc_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = 'deconvolution')
}

de_qlrt_gm_3pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 3, test_method = 'glmGamPoi', size_factors = 'gmean')
}

de_qlrt_gm_5pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 5, test_method = 'glmGamPoi', size_factors = 'gmean')
}

de_qlrt_gm_7pr <- function(mat, grouping) {
  pseudobulk_de(mat, grouping, G = 7, test_method = 'glmGamPoi', size_factors = 'gmean')
}

