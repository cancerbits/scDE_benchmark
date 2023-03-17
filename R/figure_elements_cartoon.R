# only run once to create toy data figures for introductory figure

config <- yaml::read_yaml("config.yaml")

FIG_DIR <- file.path(config$out_root, 'figures')

set.seed(9)
# number of genes and cells
G <- 12
C <- 9

# Create raw counts for two groups
# use poisson model
mu_range <- c(0.1, 3)
mu_vec1 <- runif(n = G, min = mu_range[1], max = mu_range[2])
mu_vec2 <- mu_vec1
mu_vec2[9] <- 0.2
mu_vec2[2] <- 2

mat1 <- t(sapply(1:G, function(i) rpois(n = C, lambda = mu_vec1[i])))
mat2 <- t(sapply(1:G, function(i) rpois(n = C, lambda = mu_vec2[i])))

suppressPackageStartupMessages(library('ComplexHeatmap'))

cf <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, 
            gp = gpar(col = "grey33", fill = NA))
}
reds <- hcl.colors(n = 9, palette = 'Reds')
blues <- hcl.colors(n = 9, palette = 'Blues')

hm1 <- Heatmap(mat1, cluster_rows = FALSE, cluster_columns = FALSE, 
               col = rev(reds), border = TRUE, 
               cell_fun = cf, row_title = 'Genes', show_heatmap_legend = FALSE,
               column_title = 'Cells', column_title_side = 'bottom',
               top_annotation = HeatmapAnnotation(
                 foo = anno_block(labels = 'Group1', gp = gpar(fill = reds[5]))))

hm2 <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE, 
               col = rev(blues), border = TRUE, 
               cell_fun = cf, row_title = 'Genes', show_heatmap_legend = FALSE,
               column_title = 'Cells', column_title_side = 'bottom',
               top_annotation = HeatmapAnnotation(
                 foo = anno_block(labels = 'Group2', gp = gpar(fill = blues[5]))))

draw(hm1 + hm2,
     column_title = "Count matrices", column_title_gp = gpar(fontsize = 14))


cairo_pdf(filename = file.path(FIG_DIR, 'cartoon_counts.pdf'), width = 2, height = 1.8, pointsize = 13)
draw(hm1 + hm2)
dev.off()

svg(filename = file.path(FIG_DIR, 'cartoon_counts.svg'), width = 2, height = 1.8, pointsize = 13)
draw(hm1 + hm2)
dev.off()

# same for the transformed values

suppressPackageStartupMessages(library('Matrix'))

# proportional fitting followed by log1p
tr_log_pf <- function(counts) {
  counts <- as(counts, 'dgCMatrix')
  if (inherits(x = counts, what = 'dgCMatrix')) {
    pf <- sparseMatrixStats::colSums2(counts)
    sf <- mean(pf)
    tr_counts <- counts
    tr_counts@x <- log(tr_counts@x * rep.int(sf / pf, diff(tr_counts@p)) + 1)
    return(as.matrix(tr_counts))
  } else {
    stop('need dgCMatrix')
  }
}

mat1tr <- tr_log_pf(mat1)
mat2tr <- tr_log_pf(mat2)

hm1 <- Heatmap(mat1tr, cluster_rows = FALSE, cluster_columns = FALSE, 
               col = rev(reds), border = TRUE, 
               cell_fun = cf, row_title = 'Genes', show_heatmap_legend = FALSE,
               column_title = 'Cells', column_title_side = 'bottom',
               top_annotation = HeatmapAnnotation(
                 foo = anno_block(labels = 'Group1', gp = gpar(fill = reds[5]))))

hm2 <- Heatmap(mat2tr, cluster_rows = FALSE, cluster_columns = FALSE, 
               col = rev(blues), border = TRUE, 
               cell_fun = cf, row_title = 'Genes', show_heatmap_legend = FALSE,
               column_title = 'Cells', column_title_side = 'bottom',
               top_annotation = HeatmapAnnotation(
                 foo = anno_block(labels = 'Group2', gp = gpar(fill = blues[5]))))

draw(hm1 + hm2)


cairo_pdf(filename = file.path(FIG_DIR, 'cartoon_norm.pdf'), width = 2, height = 1.8, pointsize = 13)
draw(hm1 + hm2)
dev.off()

svg(filename = file.path(FIG_DIR, 'cartoon_norm.svg'), width = 2, height = 1.8, pointsize = 13)
draw(hm1 + hm2)
dev.off()

# pseudo-bulk of pseudo-replicates
pseudobulk <- function(counts, grouping) {
  counts <- as(counts, 'dgCMatrix')
  mat <- sapply(levels(grouping), function(gr) {
    sparseMatrixStats::rowSums2(x = counts, cols = grouping == gr)
  })
  colnames(mat) <- NULL
  return(mat)
}

mat1gr <- pseudobulk(counts = mat1, grouping = as.factor(ceiling((1:C)/3)))
mat2gr <- pseudobulk(counts = mat2, grouping = as.factor(ceiling((1:C)/3)))

hm1 <- Heatmap(mat1gr, cluster_rows = FALSE, cluster_columns = FALSE, 
               col = rev(reds), border = TRUE, 
               cell_fun = cf, row_title = 'Genes', show_heatmap_legend = FALSE,
               column_title = 'Pseudo-replicates', column_title_side = 'bottom',
               top_annotation = HeatmapAnnotation(
                 foo = anno_block(labels = 'Group1', gp = gpar(fill = reds[5]))))

hm2 <- Heatmap(mat2gr, cluster_rows = FALSE, cluster_columns = FALSE, 
               col = rev(blues), border = TRUE, 
               cell_fun = cf, row_title = 'Genes', show_heatmap_legend = FALSE,
               column_title = 'Pseudo-replicates', column_title_side = 'bottom',
               top_annotation = HeatmapAnnotation(
                 foo = anno_block(labels = 'Group2', gp = gpar(fill = blues[5]))))

draw(hm1 + hm2)


cairo_pdf(filename = file.path(FIG_DIR, 'cartoon_bulk.pdf'), width = 1.8, height = 1.8, pointsize = 13)
draw(hm1 + hm2)
dev.off()

svg(filename = file.path(FIG_DIR, 'cartoon_bulk.svg'), width = 1.8, height = 1.8, pointsize = 13)
draw(hm1 + hm2)
dev.off()

# Evaluation cartoon

M <- 4
P <- 22

evalmat <- t(sapply(1:P, function(i) rpois(n = M, lambda = i)))
evalmat <- log(evalmat + 5)
hm3 <- Heatmap(evalmat, col = hcl.colors(n = 9),
               row_title = '343 DE pipelines', 
               column_title = 'Metrics', column_title_side = 'bottom',
               show_row_dend = FALSE, show_column_dend = FALSE,
               border = TRUE, cell_fun = cf,
               heatmap_legend_param = list(
                 title = '', at = range(evalmat), 
                 labels = c("Worst", "Best"), border = FALSE
               ))
draw(hm3)

fname <- 'cartoon_eval'
cairo_pdf(filename = file.path(FIG_DIR, sprintf('%s.pdf', fname)), width = 1.8, height = 4, pointsize = 13)
draw(hm3)
dev.off()

svg(filename = file.path(FIG_DIR, sprintf('%s.svg', fname)), width = 1.8, height = 4, pointsize = 13)
draw(hm3)
dev.off()

# Results cartoon

pvals <- sapply(1:nrow(mat1tr), function(i) {
  t.test(mat1tr[i, ], mat2tr[i, ])$p.value
})
resmat <- matrix(c(pvals, abs((pvals < 0.1) - 1)), nrow = length(pvals))
colnames(resmat) <- c('p-val', 'is DE')
hm4 <- Heatmap(resmat, col = hcl.colors(n = 9, palette = 'grays'),
               cluster_rows = FALSE, cluster_columns = FALSE,
               column_names_side = "top", 
               show_heatmap_legend = FALSE,
               border = TRUE, cell_fun = cf)
draw(hm4)

fname <- 'cartoon_res'
cairo_pdf(filename = file.path(FIG_DIR, sprintf('%s.pdf', fname)), width = 0.7, height = 1.8, pointsize = 13)
draw(hm4)
dev.off()

svg(filename = file.path(FIG_DIR, sprintf('%s.svg', fname)), width = 0.7, height = 1.8, pointsize = 13)
draw(hm4)
dev.off()
