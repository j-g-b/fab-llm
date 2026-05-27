#
set.seed(137)

#
pSeq <- SeuratDisk::LoadH5Seurat("data/K562_gwps_normalized_bulk_01.h5seurat")
embeds <- reticulate::py_load_object("data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle")

#
gene_map <- sapply(row.names(pSeq@meta.data), function(x){
  splitl <- strsplit(x, "_")
  c(splitl[[1]][2], splitl[[1]][4])
}) |>
  t() |>
  unique()
gene_map <- gene_map[grepl("ENS", row.names(gene_map)), ]
row.names(gene_map) <- gene_map[, 2]
gene_map <- gene_map[, 1, drop = F]
gns_to_use <- intersect(row.names(pSeq@assays$RNA$scale.data), row.names(gene_map[gene_map[, 1] %in% names(embeds),, drop=F]))
gene_map <- gene_map[gns_to_use, , drop=F]

#
E <- sapply(1:length(gene_map[,1]), function(j){
  embeds[[gene_map[j,1]]]
}) |>
  t() |>
  magrittr::set_rownames(gene_map[,1])
E <- scale(E)

#
d <- 10
UE <- svd(E, nv = 0, nu = d)[["u"]][, 1:d]

#
cntrl_profiles <- pSeq@assays$RNA$scale.data[row.names(gene_map), grepl("non-targeting", colnames(pSeq@assays$RNA$scale.data))]
pert_profiles <- pSeq@assays$RNA$scale.data[row.names(gene_map), !grepl("non-targeting", colnames(pSeq@assays$RNA$scale.data))]
cntrl_profiles[is.infinite(cntrl_profiles)] <- 0
pert_profiles[is.infinite(pert_profiles)] <- 0

#
p <- nrow(pert_profiles)
n1 <- ncol(cntrl_profiles)
n2 <- 1
Sigmai <- solve(cov(t(crossprod(UE, cntrl_profiles))))
cntrl_mean <- rowMeans(cntrl_profiles)

rsamp <- sample(1:ncol(cntrl_profiles), floor(2*ncol(cntrl_profiles)/3))
svU <- svd(cntrl_profiles[, setdiff(1:ncol(cntrl_profiles), rsamp)] - rowMeans(cntrl_profiles[, setdiff(1:ncol(cntrl_profiles), rsamp)])%*%t(rep(1, ncol(cntrl_profiles) - length(rsamp))))
Sigma_tildU <- svU[["u"]][, 1:(ncol(cntrl_profiles) - length(rsamp) - 1)]%*%diag(1/svU[["d"]][1:(ncol(cntrl_profiles) - length(rsamp) - 1)])
UEtild <- UE - Sigma_tildU%*%crossprod(Sigma_tildU, UE)
Sigma_hati <- solve(cov(t(crossprod(UEtild, cntrl_profiles[, rsamp]))))
nsplt <- length(rsamp)
spltmean <- rowMeans(cntrl_profiles[, rsamp])

p_fab <- rep(NA, ncol(pert_profiles))
p_fab_splt <- rep(NA, ncol(pert_profiles))
p_rp <- rep(NA, ncol(pert_profiles))
for(j in 1:ncol(pert_profiles)){
  #
  py <- crossprod(UE, pert_profiles[, j] - cntrl_mean)
  T2FAB <- c(((n1+n2-d-1)/(d*(n1+n2-2)))*(n1*n2/(n1+n2))*crossprod(py, crossprod(Sigmai, py)))
  p_fab[j] <- 1 - pf(T2FAB, d, n1 + n2 - 1 - d)
  
  #
  py <- crossprod(UEtild, pert_profiles[, j] - spltmean)
  T2FABsplt <- c(((nsplt+n2-d-1)/(d*(nsplt+n2-2)))*(nsplt*n2/(nsplt+n2))*crossprod(py, crossprod(Sigma_hati, py)))
  p_fab_splt[j] <- 1 - pf(T2FABsplt, d, nsplt + n2 - 1 - d)
  
  #
  Z <- matrix(rnorm(p*d), nrow = p)
  py <- crossprod(Z, pert_profiles[, j] - cntrl_mean)
  #Si <- solve(cov(t(crossprod(Z, cntrl_profiles))))
  T2 <- c(((n1+n2-d-1)/(d*(n1+n2-2)))*(n1*n2/(n1+n2))*crossprod(py, solve(cov(t(crossprod(Z, cntrl_profiles))), py)))
  p_rp[j] <- 1 - pf(T2, d, n1 + n2 - 1 - d)
  print(j)
}

#
fdr_rp <- p.adjust(p_rp, method = "BH")
fdr_fab <- p.adjust(p_fab, method = "BH")
fdr_fab_splt <- p.adjust(p_fab_splt, method = "BH")

#
sens_res <- matrix(c(sum(fdr_rp < 0.01), sum(fdr_fab < 0.01),
                     sum(fdr_rp < 0.025), sum(fdr_fab < 0.025),
                     sum(fdr_rp < 0.05), sum(fdr_fab < 0.05),
                     sum(fdr_rp < 0.1), sum(fdr_fab < 0.1)), ncol = 2, byrow = T)
colnames(sens_res) <- c("RP", "FAB")
row.names(sens_res) <- paste0("FDR = ", c(0.01, 0.025, 0.05, 0.1))
saveRDS(sens_res, file = "results/perturb_seq10_sens.rds")


