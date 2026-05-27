#
set.seed(2501)
source("helpers.R")

#
cl_info <- readr::read_csv("data/Model.csv")
gd <- readr::read_csv("data/CRISPRGeneDependency.csv")
cols_to_use <- gd[,-1] |> 
               as.matrix() |> 
               apply(2, function(x){
                 sum(x > 0.5, na.rm = T) > 50 & sum(x > 0.5, na.rm = T) < nrow(gd) - 50
               })
colnames(gd)[1] <- "ModelID"

#
gene_cols <- which(cols_to_use) + 1
res_list <- list()
res_null_list <- list()
for(j in 1:length(gene_cols)){
  #
  df <- dplyr::left_join(cl_info, gd[, c(1, gene_cols[j])])
  df <- df[, c(1, 6, 14, 16, 18, ncol(df))] |>
        dplyr::distinct() |>
        dplyr::filter(!is.na(Age), !is.na(Sex),
                      !is.na(PrimaryOrMetastasis),
                      PrimaryOrMetastasis != "Unknown",
                      Sex != "Unknown",
                      )
  df <- df[!is.na(c(df[, ncol(df)])[[1]]), ]
  gene_name <- colnames(df)[ncol(df)]
  colnames(df)[ncol(df)] <- "gene_name"
  df$gene_name <- ifelse(df$gene_name > 0.5, 1, 0)
  
  #
  res <- glm(gene_name~Age+Sex+PrimaryOrMetastasis, data = df, family = "binomial")
  if(res$converged & !res$boundary){
    res_list[[gene_name]] <- summary(res)$coefficients[4, ]
  }
  
  #
  df_null <- df
  pri_bin <- ifelse(df$PrimaryOrMetastasis == "Primary", 1, 0)
  null_logit <- log(res$fitted.values/(1-res$fitted.values)) - pri_bin*res$coefficients[4]
  null_probs <- exp(null_logit)/(1+exp(null_logit))
  df_null$gene_name <- sapply(1:nrow(df), function(i){
    sample(c(0, 1), 1, prob = c(1-null_probs[i], null_probs[i]))
  })
  res_null <- glm(gene_name~Age+Sex+PrimaryOrMetastasis, data = df_null, family = "binomial")
  if(res_null$converged & !res_null$boundary){
    res_null_list[[gene_name]] <- summary(res_null)$coefficients[4, ]
  }
}

#
res_table <- do.call(rbind, res_list)
res_table <- res_table[res_table[, 2] < 5, ]

#
res_null_table <- do.call(rbind, res_null_list)
res_null_table <- res_null_table[res_null_table[, 2] < 5, ]

#
embeds <- reticulate::py_load_object("data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle")

#
gns_to_use <- intersect(gsub(" [(].*", "", row.names(res_table)), names(embeds)) |>
              intersect(gsub(" [(].*", "", row.names(res_null_table)))

#
d <- 100
E <- sapply(1:length(gns_to_use), function(j){
  embeds[[gns_to_use[j]]]
}) |>
  t() |>
  magrittr::set_rownames(gns_to_use)
E <- scale(E)
svE <- svd(E)
E <- svE[["u"]][, 1:d]%*%diag(sqrt(svE[["d"]][1:d]))
row.names(E) <- gns_to_use
row.names(res_table) <- gsub(" [(].*", "", row.names(res_table))
row.names(res_null_table) <- gsub(" [(].*", "", row.names(res_null_table))

#
p_vals <- fabp_lin_reg_z(Y = res_table[gns_to_use, 1, drop=F], 
                              R = matrix(1, nrow=length(gns_to_use)) |>
                                  magrittr::set_rownames(gns_to_use),
                              S = res_table[gns_to_use, 2, drop=F],
                              U = E, V = matrix(1, 1, 1),
                              pool_sampling_var = F)

#
null_p_vals <- fabp_lin_reg_z(Y = res_null_table[gns_to_use, 1, drop=F], 
                              R = matrix(1, nrow=length(gns_to_use)) |>
                                magrittr::set_rownames(gns_to_use),
                              S = res_null_table[gns_to_use, 2, drop=F],
                              U = E, V = matrix(1, 1, 1),
                              pool_sampling_var = F)

#
adapt_p <- adaptMT::adapt_glmnet(x = E, pvals = p_vals$p, alphas = seq(0.01, 0.25, by = 0.01))
comparison_res <- sapply(1:length(adapt_p$alphas), function(i){
  c(sum(p_vals$fdr_p < adapt_p$alphas[i]),
    sum(p_vals$fdr_fabp < adapt_p$alphas[i]),
    adapt_p$nrejs[i])
}) |>
  t()
pdf(file = "figures/adapt_compare_z.pdf", family = "Times", height = 5, width = 5)
plot(comparison_res[, 1], adapt_p$alphas, type = 'l', lwd = 1.5, col = 'black', 
     xlim = c(0, 118), xlab = "Rank", ylab = "FDR")
lines(comparison_res[, 2], adapt_p$alphas, lwd = 1.5, col = "#08519c")
lines(comparison_res[, 3], adapt_p$alphas, lwd = 1.5, col = "#756bb1")
abline(h = 0.1, col = '#a50f15')
legend(0, y = 0.25, legend = c("z test", "FAB test", "AdaPT"),
       lty = 1, lwd = 1.5, col = c('black', "#08519c", "#756bb1"))
dev.off()

#
pdf(file = "figures/logistic_reg.pdf", family = "Times", height = 5, width = 9)
layout(matrix(c(1, 2), nrow = 1))
plot(p_vals$observed, p_vals$predicted, main = "",
     xlab = "Observed", ylab = "Predicted")
abline(0, 1, col = "#a50f15", lty = 2)
plot(-log10(p_vals$fdr_fabp), -log10(p_vals$fdr_p),
     col = ifelse(-log10(p_vals$fdr_fabp) > 1 & -log10(p_vals$fdr_p) < 1,
                  "#08519c", ifelse(-log10(p_vals$fdr_fabp) < 1 & -log10(p_vals$fdr_p) > 1,
                                    "#a50f15", ifelse(-log10(p_vals$fdr_fabp) > 1 & -log10(p_vals$fdr_p) > 1, "black", "gray"))),
     xlab = "FDR FAB z (-log10 scale)", ylab = "FDR Wald z (-log10 scale)")
abline(h = 1, lty = 3, col = "darkgray")
abline(v = 1, lty = 3, col = "darkgray")
abline(0, 1, col = "darkgray", lty = 2)
points(x=rep(2.25, 3), y=1.1*c(0.75, 0.5, 0.25)-0.1, col = c("#08519c", 'black', 'gray'), lty = c(1,1,1))
text(x=rep(2.5, 3), y=1.1*c(0.75, 0.5, 0.25)-0.1, c("FAB only", "Both", "Neither"),
     cex = 0.9, adj = 0)
dev.off()

#
pdf(file = "figures/logistic_reg_null.pdf", family = "Times", height = 5, width = 5)
plot(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(p_vals$p < q)}),
     type = 'l', col = 'black', lwd = 2, ylab = "Fraction p-values below threshold", xlab = "Threshold",
     xlim = c(0, 1), ylim = c(0, 1), main = "z-test")
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(p_vals$fabp < q)}),
      type = 'l', col = '#08519c', lwd = 2)
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(null_p_vals$p < q)}),
      type = 'l', col = 'black', lty = 2, lwd = 1.5)
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(null_p_vals$fabp < q)}),
      type = 'l', col = '#08519c', lty = 2, lwd = 1.5)
abline(0, 1, col = '#a50f15')
legend(x = 0, y = 0.9, legend = c("Null data", "Real data"), lty = c(2, 1), lwd = c(1.5, 2))
dev.off()

# Sensitivity analysis
E <- sapply(1:length(gns_to_use), function(j){
  embeds[[gns_to_use[j]]]
}) |>
  t() |>
  magrittr::set_rownames(gns_to_use)
E <- scale(E)
svE <- svd(E)
fdrs <- c(0.01, 0.025, 0.05, 0.1)
Ds <- c(10, 50, 100, 200)
sens_res <- array(NA, dim = c(length(fdrs), 2, length(Ds)))
for(j in 1:length(Ds)){
  #
  d <- Ds[j]
  
  #
  E <- svE[["u"]][, 1:d]%*%diag(sqrt(svE[["d"]][1:d]))
  row.names(E) <- gns_to_use
  row.names(res_table) <- gsub(" [(].*", "", row.names(res_table))
  row.names(res_null_table) <- gsub(" [(].*", "", row.names(res_null_table))
  
  #
  p_vals <- fabp_lin_reg_z(Y = res_table[gns_to_use, 1, drop=F], 
                           R = matrix(1, nrow=length(gns_to_use)) |>
                             magrittr::set_rownames(gns_to_use),
                           S = res_table[gns_to_use, 2, drop=F],
                           U = E, V = matrix(1, 1, 1),
                           pool_sampling_var = F)
  
  #
  for(k in 1:length(fdrs)){
    sens_res[k, 1, j] <- sum(p_vals$fdr_p < fdrs[k])
    sens_res[k, 2, j] <- sum(p_vals$fdr_fabp < fdrs[k])
  }
}
