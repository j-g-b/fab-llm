#
source("helpers.R")
library(tidyverse)
library(magrittr)
set.seed(139875)

#
ras_mut <- readr::read_csv("data/diff_dep/ras_mut_status.csv") %>%
  reshape2::melt(id.vars = "status", variable.name = "cl") %>%
  dplyr::group_by(cl) %>%
  dplyr::filter(value[status == "is AML"] == 1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(status != "is AML") %>%
  dplyr::group_by(cl) %>%
  dplyr::summarise(is_ras_mut = any(value == 1)) %>%
  dplyr::ungroup()
table <- readr::read_csv("data/diff_dep/focused_screen.csv") %>%
  reshape2::melt(id.vars = "Gene", variable.name = "cl") %>%
  dplyr::right_join(ras_mut) %>%
  dplyr::group_by(Gene, is_ras_mut) %>%
  dplyr::summarise(r = sum(!is.na(value)),
                   s = sd(value, na.rm = T),
                   m = mean(value, na.rm = T)) %>%
  dplyr::ungroup() %>%
  plyr::ddply(plyr::.(Gene), function(d){
    data.frame(r = d$r[d$is_ras_mut],
               s = d$s[d$is_ras_mut],
               m = d$m[d$is_ras_mut],
               r1 = d$r[!d$is_ras_mut],
               s1 = d$s[!d$is_ras_mut],
               m1 = d$m[!d$is_ras_mut])
  })
#
Ybar <- table$m %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
Ybar1 <- table$m1 %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
S <- table$s %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
S1 <- table$s1 %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
R <- table$r %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
R1 <- table$r1 %>%
  matrix(nrow = 1) %>%
  magrittr::set_colnames(table$Gene)
#
embeds <- reticulate::py_load_object("data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle")
#
d <- 100
E <- sapply(1:length(table[["Gene"]]), function(j){
  embeds[[table[["Gene"]][j]]]
}) |>
  t() |>
  magrittr::set_rownames(table[["Gene"]])
E <- scale(E)
E <- svd(E)[["u"]][, 1:d]%*%diag(sqrt(svd(E)[["d"]][1:d]))
gns_to_use <- table$Gene

#
Ybar <- Ybar[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
S <- S[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
R <- R[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
#
Ybar1 <- Ybar1[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
S1 <- S1[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
R1 <- R1[, gns_to_use, drop = F] %>% magrittr::set_rownames("diff_score")
#
p_vals <- fabp_lin_reg(Y = Ybar, S = S, R = R, 
                       U = matrix(1, nrow = 1, ncol = 1), V = E, 
                       pool_sampling_var = F,
                       Y1 = Ybar1,
                       S1 = S1,
                       R1 = R1)

#
Ybarnull <- Ybar
Ybarnull1 <- Ybar
Ybarnull[1, ] <- rnorm(ncol(Ybar), sd = S[1, ]/sqrt(R[1, ]))
Ybarnull1[1, ] <- rnorm(ncol(Ybar), sd = S1[1, ]/sqrt(R1[1, ]))
Snull <- S
Snull1 <- S
Snull[1, ] <- S[1, ]*sqrt(rchisq(ncol(Ybar), df = R[1, ]-1))/sqrt(R[1, ])
Snull1[1, ] <- S1[1, ]*sqrt(rchisq(ncol(Ybar), df = R1[1, ]-1))/sqrt(R1[1, ])
null_p_vals <- fabp_lin_reg(Y = Ybarnull, S = Snull, R = R, 
                            U = matrix(1, nrow = 1, ncol = 1), V = E, 
                            pool_sampling_var = F,
                            Y1 = Ybarnull1,
                            S1 = Snull1,
                            R1 = R1)
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
rank_p <- data.frame(p = p_vals$p, fdr = p_vals$fdr_p, rank = rank(p_vals$p), type = "UMPU", row = p_vals$row)
rank_fabp <- data.frame(p = p_vals$fabp, fdr = p_vals$fdr_fabp, rank = rank(p_vals$fabp), type = "FAB", row = p_vals$row)
adapt_p <- adaptMT::adapt_glmnet(x = E, pvals = p_vals$p, alphas = seq(0.01, 0.25, by = 0.01))

#
comparison_res <- sapply(1:length(adapt_p$alphas), function(i){
  c(sum(p_vals$fdr_p < adapt_p$alphas[i]),
    sum(p_vals$fdr_fabp < adapt_p$alphas[i]),
    adapt_p$nrejs[i])
}) |>
  t()
pdf(file = "figures/adapt_compare_t.pdf", family = "Times", height = 5, width = 5)
plot(comparison_res[, 1], adapt_p$alphas, type = 'l', lwd = 1.5, col = 'black', 
     xlim = c(0, 50), xlab = "Rank", ylab = "FDR")
lines(comparison_res[, 2], adapt_p$alphas, lwd = 1.5, col = "#08519c")
lines(comparison_res[, 3], adapt_p$alphas, lwd = 1.5, col = "#756bb1")
abline(h = 0.1, col = '#a50f15')
legend(0, y = 0.25, legend = c("t test", "FAB test", "AdaPT"),
       lty = 1, lwd = 1.5, col = c('black', "#08519c", "#756bb1"))
dev.off()

comparison_tab <- sapply(c(), function(i){
  c(sum(p_vals$fdr_p < adapt_p$alphas[i]),
    sum(p_vals$fdr_fabp < adapt_p$alphas[i]),
    adapt_p$nrejs[i])
}) |>
  t()


#
pos_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_fabp < 0.1, observed > 0) %>%
  dplyr::arrange(observed)
neg_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_fabp < 0.1, observed < 0) %>%
  dplyr::arrange(observed)
pl1 <- p_vals %>%
  dplyr::mutate(x = rank(observed, ties.method = "first")) %>%
  ggplot2::ggplot() +
  geom_bar(aes(x = x, y = observed, fill = fdr_fabp < 0.1), stat = "identity", colour = "white") +
  geom_text(data = pos_df, x = 100, y = seq(0.5, 2, length.out = nrow(pos_df)), aes(label = column), size = 2, hjust = 1.1, fontface = "bold") +
  geom_text(data = neg_df, x = 25, y = seq(-2, -0.5, length.out = nrow(neg_df)), aes(label = column), size = 2, hjust = -0.1, fontface = "bold") +
  geom_segment(data = pos_df, aes(x = x, y = observed), xend = 100, yend = seq(0.5, 2, length.out = nrow(pos_df)), colour = wes_pal[5], alpha = 0.5) +
  geom_segment(data = neg_df, aes(x = x, y = observed), xend = 25, yend = seq(-2, -0.5, length.out = nrow(neg_df)), colour = wes_pal[5], alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("#969696", "#08519c"), labels = c("FAB t FDR < 0.1", "FAB t FDR > 0.1")) +
  ylim(c(-2, 2)) +
  guides(fill=guide_legend(title="")) +
  labs(x = "Genes ranked by differential dependency", y = "Differential dependency")
pos_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_p < 0.1, observed > 0) %>%
  dplyr::arrange(observed)
neg_df <- dplyr::filter(p_vals %>%
                          dplyr::mutate(x = rank(observed, ties.method = "first")), fdr_p < 0.1, observed < 0) %>%
  dplyr::arrange(observed)
pl2 <- p_vals %>%
  dplyr::mutate(x = rank(observed, ties.method = "first")) %>%
  ggplot2::ggplot() +
  geom_bar(aes(x = x, y = observed, fill = fdr_p < 0.1), stat = "identity", colour = "white") +
  geom_text(data = pos_df, x = 100, y = seq(0.5, 2, length.out = nrow(pos_df)), aes(label = column), size = 2, hjust = 1.1, fontface = "bold") +
  geom_text(data = neg_df, x = 25, y = seq(-2, -0.5, length.out = nrow(neg_df)), aes(label = column), size = 2, hjust = -0.1, fontface = "bold") +
  geom_segment(data = pos_df, aes(x = x, y = observed), xend = 100, yend = seq(0.5, 2, length.out = nrow(pos_df)), colour = wes_pal[5], alpha = 0.5) +
  geom_segment(data = neg_df, aes(x = x, y = observed), xend = 25, yend = seq(-2, -0.5, length.out = nrow(neg_df)), colour = wes_pal[5], alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("#969696", 'black'), labels = c("t FDR < 0.1", "t FDR > 0.1")) +
  ylim(c(-2, 2)) +
  guides(fill = guide_legend(title="")) +
  labs(x = "Genes ranked by differential dependency", y = "Differential dependency")
#
pp <- cowplot::plot_grid(pl1, pl2, nrow = 2, align = "hv")
pp
ggplot2::ggsave("figures/diff_dep_example.pdf", height = 6, width = 7, family = "Times")

#
pdf(file = "figures/diff_dep_null.pdf", family = "Times", height = 5, width = 5)
p_f <- p_vals$p
p_fab <- p_vals$fabp
np_f <- null_p_vals$p
np_fab <- null_p_vals$fabp
plot(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(p_f < q)}),
     type = 'l', col = 'black', lwd = 2, ylab = "Fraction p-values below threshold", xlab = "Threshold",
     xlim = c(0, 1), ylim = c(0, 1), main = "t-test")
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(p_fab < q)}),
      type = 'l', col = '#08519c', lwd = 2)
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(np_f < q)}),
      type = 'l', col = 'black', lty = 2, lwd = 1.5)
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(np_fab < q)}),
      type = 'l', col = '#08519c', lty = 2, lwd = 1.5)
abline(0, 1, col = '#a50f15')
legend(x = 0, y = 0.9, legend = c("Null data", "Real data"), lty = c(2, 1), lwd = c(1.5, 2))
dev.off()

# Sensitivity analysis
fdrs <- c(0.01, 0.025, 0.05, 0.1)
Ds <- c(10, 50, 100)
sens_res <- array(NA, dim = c(length(fdrs), 2, length(Ds)))
E <- sapply(1:length(table[["Gene"]]), function(j){
  embeds[[table[["Gene"]][j]]]
}) |>
  t() |>
  magrittr::set_rownames(table[["Gene"]])
E <- scale(E)
svE <- svd(E)
for(j in 1:length(Ds)){
  #
  d <- Ds[j]
  
  #
  Ed <- svE[["u"]][, 1:d]%*%diag(sqrt(svE[["d"]][1:d]))
  row.names(Ed) <- gns_to_use
  
  #
  p_vals <- fabp_lin_reg(Y = Ybar, S = S, R = R, 
                         U = matrix(1, nrow = 1, ncol = 1), V = Ed, 
                         pool_sampling_var = F,
                         Y1 = Ybar1,
                         S1 = S1,
                         R1 = R1)
  
  #
  for(k in 1:length(fdrs)){
    sens_res[k, 1, j] <- sum(p_vals$fdr_p < fdrs[k])
    sens_res[k, 2, j] <- sum(p_vals$fdr_fabp < fdrs[k])
  }
}
