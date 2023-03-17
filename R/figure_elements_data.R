# only run once to create cartoon-like figures of the data

source('R/setup.R')
source('R/utilities.R')
source('R/transformations.R')
source('R/de_methods.R')

library('dplyr')
library('ggplot2')
library('scales')
library('patchwork')

## muscat data

set.seed(42)

DATA_DIR <- file.path(config$out_root, 'simulated_data')
FIG_DIR <- file.path(config$out_root, 'figures')

data_files <- list.files(path = DATA_DIR, full.names = TRUE)

daf_lst <- list()
for (f in data_files) {
  data_name <- gsub(pattern = 'sim_data_', replacement = '', x = basename(f))
  data_name <- gsub(pattern = '\\.Rds$', replacement = '', x = data_name)
  message('simulation experiment ', data_name)
  
  muscat_sample <- readRDS(f)[[1]]
  
  expr_mat <- tr_sl_pc1_dc(muscat_sample$counts)
  # remove constant features
  is_const <- sparseMatrixStats::rowSds(expr_mat) <= .Machine$double.eps
  expr_mat <- as.matrix(expr_mat[!is_const, ])
  
  pca_res <- irlba::prcomp_irlba(x = t(expr_mat), n = 10, center = TRUE, scale. = TRUE)
  nn <- floor(min(table(muscat_sample$cell_metadata$group_id)) * 0.2)
  umap <- uwot::umap(X = pca_res$x, n_neighbors = nn, metric = 'cosine', min_dist = 1)
  umap <- prcomp(umap)$x
  
  daf <- data.frame(experiment = rep(data_name, nrow(umap)),
                    grp = muscat_sample$cell_metadata$group_id,
                    umap1 = umap[, 1],
                    umap2 = umap[, 2])
  
  # make sure group A is left
  umean <- group_by(daf, grp) %>% summarise(umean = mean(umap1)) %>% pull(umean)
  if (umean[2] < umean[1]) {
    daf$umap1 <- -daf$umap1
    daf$umap2 <- -daf$umap2
  }
  
  daf_lst[[data_name]] <- daf
}

daf <- as_tibble(data.table::rbindlist(daf_lst))
daf$experiment <- as.factor(daf$experiment)

# Recode the experiment names and reorder
exp_names_new <- list("20 vs 20 cells" = "20_vs_20", "200 vs 200 cells" = "200_vs_200", 
                      "2k vs 2k cells" = "2k_vs_2k", "200 vs 2k cells" = "200_vs_2k")
exp_names_new_order <- names(exp_names_new)
levels(daf$experiment) <- exp_names_new
daf$experiment <- factor(daf$experiment, levels = exp_names_new_order)
daf <- arrange(daf, experiment)

reds <- hcl.colors(n = 9, palette = 'Reds')
blues <- hcl.colors(n = 9, palette = 'Blues')

g <- ggplot(daf, aes(umap1, umap2, color = grp)) +
  geom_point(color = 'black', size = 1.6) +
  geom_point(size = 1.4) +
  facet_wrap(~ experiment, scales = 'free', ncol = 2) +
  scale_color_manual(values = c(reds[5], blues[5])) +
  theme_classic(base_size = 12) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.line = element_blank(), axis.title = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 14))

fname <- 'simulated_data_umap'
size_factor <- 1.5
w <- 2.6 * size_factor
h <- 1.3 * size_factor
cairo_pdf(filename = file.path(FIG_DIR, sprintf('%s.pdf', fname)), width = w, height = h, pointsize = 12)
show(g)
dev.off()

svg(filename = file.path(FIG_DIR, sprintf('%s.svg', fname)), width = w, height = h, pointsize = 12)
show(g)
dev.off()

## immune data

library('ggrepel')
set.seed(42)

meta_data <- readRDS(file = file.path(config$out_root, 'immune_data/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix_seurat_metadata.Rds'))
#counts <- readRDS(file = file.path(config$project_root, 'data/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix_seurat_counts.Rds'))
md <- mutate(meta_data, Celltype = broad_cell_type)

# labels
lab <- group_by(md, Celltype) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

g <- ggplot(md, aes(UMAP_1, UMAP_2, color = Celltype)) +
  geom_point(color = 'black', size = 1.6) +
  geom_point(size = 1.4) +
  geom_label_repel(data = lab, aes(label = Celltype)) +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.line = element_blank(), axis.title = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

fname <- 'immune_data_umap'
w <- 2.4
h <- 2
cairo_pdf(filename = file.path(FIG_DIR, sprintf('%s.pdf', fname)), width = w, height = h, pointsize = 13)
show(g)
dev.off()

svg(filename = file.path(FIG_DIR, sprintf('%s.svg', fname)), width = w, height = h, pointsize = 13)
show(g)
dev.off()

## some detailed results for simulation experiments

# define the pipelines of interest
poi <- c('counts-deseq_wald_3pr', 'counts-edger_qlf_3pr', 
         'sl_pc1_10k-wilcox', 'counts-mast')

experiment <- '2k_vs_2k'
RES_DIR <- file.path(config$out_root, 'results', 'simulated_data')

rds_files <- list.files(path = RES_DIR, pattern = sprintf('^%s-.*\\.Rds$', experiment), full.names = TRUE)
# only read the files of the pipelines of interest
in_poi <- sapply(poi, function(x) grepl(pattern = x, basename(rds_files)))
rds_files <- rds_files[rowSums(in_poi) > 0]
gene_mean_breaks <- c(-Inf, 0.1, 1, Inf)

res_list <- lapply(X = rds_files, FUN = function(f) {
  f_parts <- stringr::str_split(string = stringr::str_remove(basename(f), pattern = '\\.Rds$'), pattern = '-')[[1]]
  replicate <- as.numeric(f_parts[2])
  transformation <- f_parts[3]
  de_method <- f_parts[4]
  pipeline <- paste(transformation, de_method, sep = '-')
  
  res <- readRDS(f) #$res %>% mutate(file = basename(f))
  transformation_time = dplyr::filter(res$timing, step == 'transformation')$time
  method_time = dplyr::filter(res$timing, step == 'de_method')$time
  
  # add geometric mean of simulation, and pipeline label, and fix missing p-values
  de_res <- res$res %>% 
    mutate(sim_mean = expm1((log1p(sim_mean.A) + log1p(sim_mean.B)) / 2),
           expr_grp = cut(sim_mean, breaks = gene_mean_breaks, labels = c('low', 'medium', 'high')),
           pval = case_when(is.na(pval) ~ 1, TRUE ~ pval),
           FDR = p.adjust(pval, method = 'fdr'),
           pred_de = factor(FDR < 0.05, levels = c(FALSE, TRUE)),
           pred_de3 = factor(effect_direction * (FDR < 0.05), levels = c(-1, 0, 1)),
           true_de = factor(!is.na(logFC), levels = c(FALSE, TRUE)),
           true_de3 = factor(case_when(is.na(logFC) ~ 0, logFC > 0 ~ 1, logFC < 0 ~ -1), levels = c(-1, 0, 1)),
           replicate = replicate,
           transformation = transformation,
           de_method = de_method,
           pipeline = pipeline)
  
  return(de_res)
})
res <- as_tibble(data.table::rbindlist(res_list)) %>%
  mutate(across(where(is.character), as.factor))

# prec as function of mean
plot_data <- filter(res, pred_de3 != 0) %>%
  mutate(TP = pred_de3 == true_de3)

p1 <- ggplot(plot_data, aes(log10(sim_mean), as.numeric(TP), color = pipeline)) +
  geom_hline(yintercept = 0.95, alpha = 1, linetype = 2, color = 'black') +
  geom_smooth(se = TRUE, method = "gam", method.args = list(family = "binomial"), formula = y ~ s(x, bs = "cs")) +
  coord_cartesian(expand = FALSE) +
  ylab('Precision') +
  xlab('Gene mean') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(0.05, 0.35, 0.65, 0.95)) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), 
                     labels = label_math(expr = 10^.x, format = force)) +
  annotation_logticks(sides = 'b')

# sens as function of mean
plot_data <- filter(res, true_de3 != 0) %>%
  mutate(TP = pred_de3 == true_de3)

p2 <- ggplot(plot_data, aes(log10(sim_mean), as.numeric(TP), color = pipeline)) +
  geom_hline(yintercept = 0.95, alpha = 1, linetype = 2, color = 'black') +
  geom_smooth(se = FALSE, method = "gam", method.args = list(family = "binomial"), formula = y ~ s(x, bs = "cs")) +
  coord_cartesian(expand = FALSE) +
  ylab('Sensitivity') +
  xlab('Gene mean') +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(0.05, 0.35, 0.65, 0.95)) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), 
                     labels = label_math(expr = 10^.x, format = force)) +
  annotation_logticks(sides = 'b')

sf <- 1.2
cairo_pdf(filename = file.path(FIG_DIR, 'simulated_data_prec_sens.pdf'), 
          width = 3.9 * sf, height = 1.33 * sf)
show(p1 + p2)
dev.off()

