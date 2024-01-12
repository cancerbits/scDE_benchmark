

# riffle with empirical p-value
de_random_emp <- function(mat, grouping, ...) {
  ro <- riffle::diff_mean_test(y = mat, group_labels = grouping, 
                               R = 499, log2FC_th = -Inf, mean_th = -Inf, cells_th = 0, verbosity = 0)
  res <- select(ro, gene, emp_pval, emp_pval_adj, zscore) %>%
    mutate(zscore = sign(-zscore)) %>%
    dplyr::rename(feature = gene, pval = emp_pval, FDR = emp_pval_adj, effect_direction = zscore)
  return(res)
}

# riffle with approximated p-value
de_random_app <- function(mat, grouping, ...) {
  ro <- riffle::diff_mean_test(y = mat, group_labels = grouping, 
                               R = 49, log2FC_th = -Inf, mean_th = -Inf, cells_th = 0, verbosity = 0)
  res <- select(ro, gene, pval, pval_adj, zscore) %>%
    mutate(zscore = sign(-zscore)) %>%
    dplyr::rename(feature = gene, FDR = pval_adj, effect_direction = zscore)
  return(res)
}

# wilcoxon rank-sum test as implemented in presto
de_wilcox <- function(mat, grouping, ...) {
  ref_grp <- levels(grouping)[2]
  res <- presto::wilcoxauc(X = mat, y = grouping) %>%
    filter(group == ref_grp) %>%
    mutate(effect_direction = sign(auc - 0.5)) %>%
    dplyr::rename(FDR = padj) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

# many ttests; main idea and most code is taken from Rfast ttests
de_ttest <- function(mat, grouping, ...) {
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

# Rfast-style logistic regression
de_logreg <- function(mat, grouping, ...) {
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

# MAST
de_mast <- function(mat, grouping, ...) {
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

# scDD
de_scdd <- function(mat, grouping, ...) {
  mat <- as.matrix(mat)
  # create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays=list(normcounts = mat), 
    colData = data.frame(group_id = grouping))
  
  prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
  scDD_res <- scDD::scDD(SCdat = sce, 
                         condition = 'group_id', 
                         prior_param = prior_param, 
                         testZeroes = TRUE, 
                         categorize = FALSE, 
                         param = BiocParallel::SerialParam()) %>%
    scDD::results()
  
  res <- dplyr::select(scDD_res, gene, combined.pvalue) %>%
    dplyr::rename(feature = gene, pval = combined.pvalue) %>%
    dplyr::mutate(FDR = p.adjust(p = pval, method = 'fdr'))
  
  # add simple mean-based effect direction estimate
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  m1 <- matrixStats::rowMeans2(x = mat, useNames = FALSE, cols = sel1)
  m2 <- matrixStats::rowMeans2(x = mat, useNames = FALSE, cols = sel2)
  res$effect_direction = sign(m2 - m1)
  
  return(res)
}

# for single-cell count data

# quasi-likelihood ratio testing as implemented in glmGamPoi
# valid size factor methods: 10k, gmean, pf, normed_sum, deconvolution, poscounts
de_qlrt <- function(mat, grouping, size_factor_method = 'none', ...) {
  if (is.na(size_factor_method) | is.null(size_factor_method) | size_factor_method == 'none') {
    size_factors <- FALSE
  } else {
    size_factors <- size_factors(counts = mat, method = size_factor_method)
  }
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

# single-cell parameterization of deseq
# based on http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
de_scdeseq <- function(mat, grouping, ...) {
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

de_scedger <- function(mat, grouping, test_type, ...) {
  de_edger(mat, grouping, test_type = test_type)
}

# taken from https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
de_scedger_qlfdetrate <- function(mat, grouping, ...) {
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


## methods designed for bulk data

de_edger <- function(mat, grouping, test_type, ...) {
  design <- model.matrix(~ grouping)
  
  y <- edgeR::DGEList(counts = mat, group = grouping)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  
  if (test_type == 'exact') {
    test_res <- edgeR::exactTest(y)
  }
  else if (test_type %in% c('LRT', 'lrt')) {
    fit <- edgeR::glmFit(y, design)
    test_res <- edgeR::glmLRT(fit)
  }
  else if (test_type %in% c('QLF', 'qlf')) {
    fit <- edgeR::glmQLFit(y, design)
    test_res <- edgeR::glmQLFTest(fit)
  }
  else {
    stop('test_type must be "exact", "LRT", "lrt", "QLF", or "qlf"')
  }
  
  res <- test_res$table %>% tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename(pval = PValue) %>%
    mutate(FDR = p.adjust(pval, method = 'fdr'), effect_direction = sign(logFC)) %>%
    select(feature, pval, FDR, effect_direction)
  return(res)
}

de_limma <- function(mat, grouping, test_type, ...) {
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

de_deseq <- function(mat, grouping, test_type, ...) {
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

pseudobulk_de <- function(mat, grouping, G, aggregation_strategy, test_method, ...) {
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  n1 <- sum(sel1)
  n2 <- sum(sel2)
  if (n1 < G | n2 < G) {
    stop('too few cells')
  }
  
  # create replicate labels per group
  replicate <- numeric(length = n1 + n2)
  if (aggregation_strategy == 'cluster') {
    replicate[sel1] <- cluster_counts(mat[, sel1, drop = FALSE], G)
    replicate[sel2] <- cluster_counts(mat[, sel2, drop = FALSE], G) + G
  } else {
    replicate[sel1] <- sample(1:n1 %% G)
    replicate[sel2] <- sample(1:n2 %% G) + G
  }
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
    res <- de_qlrt(mat = pb, grouping = pb_group, ...)
  } else if (test_method == 'limma') {
    res <- de_limma(mat = pb, grouping = pb_group, ...)
  } else {
    stop('test_method unknown')
  }
  return(res)
}



## refactored code to run one de method 

run_de <- function(mat, grouping, test_method, test_type = NULL, 
                       size_factor_method = NULL, G = NULL, 
                       aggregation_strategy = NULL, seed = NULL) {
  
  if (is.numeric(G) && G > 0) {
    fn <- 'pseudobulk_de'
  } else {
    fn <- paste0('de_', test_method)
  }
  #args = c(list(mat = mat, grouping = grouping), as.list(...))
  
  if (!is.null(seed)) {
    set.seed(seed = seed, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')
  }
  stime = system.time({
    res <- do.call(what = fn, args = list(mat = mat, grouping = grouping, test_method = test_method, 
                                          test_type = test_type, size_factor_method = size_factor_method,
                                          G = G, aggregation_strategy = aggregation_strategy))
  })
  message(pprint_time(stime[3]))
  res <- list(res = res, time = stime)
  return(res)
}
