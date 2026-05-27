#
library(reticulate)
source("helpers.R")
set.seed(5234891)

#
gene_names <- readr::read_csv("data/Census_allThu Sep 26 16_48_50 2024.csv")[["Gene Symbol"]]
# Read in gene embeddings
embeds <- reticulate::py_load_object("data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle")
# Read in gene expression data
CCLE <- readr::read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
colnames(CCLE) <- gsub(" [(].*", "", colnames(CCLE))
gene_names <- intersect(gene_names, intersect(names(embeds), colnames(CCLE)))

#
E <- sapply(1:length(gene_names), function(j){
  embeds[[gene_names[j]]]
}) |>
  t() |>
  magrittr::set_rownames(gene_names)

#
cl_names <- CCLE[[1]]
CCLE <- CCLE[, colnames(CCLE) %in% gene_names]
CCLE <- as.matrix(CCLE)
row.names(CCLE) <- cl_names

# Read CCLE metadata
cl_metadata <- readr::read_csv("data/Model.csv")
cl_metadata <- dplyr::left_join(data.frame(ModelID = row.names(CCLE)), cl_metadata)
trt_metadata <- readr::read_csv("data/Repurposing_Public_24Q2_Treatment_Meta_Data.csv")

# Read Repurposing data
REP <- readr::read_csv("data/Repurposing_Public_24Q2_LFC_COLLAPSED.csv")
REP <- dplyr::filter(REP, !grepl("Failure", broad_id))
REP <- dplyr::mutate(REP, cell_line_id = gsub("::.*", "", row_id))
REP <- reshape2::acast(REP, cell_line_id ~ broad_id, value.var = "LFC", fun.aggregate = function(x){mean(x, na.rm=T)})

# Get common dims
common_cls <- intersect(row.names(REP), row.names(CCLE))
REP <- REP[common_cls, ]
CCLE <- CCLE[common_cls, ]

# Study
complete_data_cls <- apply(REP, 1, function(x){all(!is.na(x))})
complete_data_gns <- apply(CCLE, 2, function(x){all(!is.na(x))})
CCLE <- CCLE[complete_data_cls, complete_data_gns]
REP <- REP[complete_data_cls, ]
E <- E[complete_data_gns, ]
E <- scale(E)

#
lineage_matrix <- model.matrix(~OncotreeLineage, cl_metadata[!is.na(cl_metadata[["OncotreeLineage"]]), ])
row.names(lineage_matrix) <- cl_metadata$ModelID[!is.na(cl_metadata[["OncotreeLineage"]])]
lineage_matrix <- lineage_matrix[row.names(REP), ]
lmU <- svd(lineage_matrix)[["u"]][, 1:(ncol(lineage_matrix)-1)]
lmQ <- svd(diag(nrow(lmU)) - tcrossprod(lmU))[["u"]][,1:(nrow(lmU) - ncol(lmU))]

#
CCLE <- crossprod(lmQ, CCLE)
REP <- crossprod(lmQ, REP)

#
n <- nrow(CCLE)
p <- ncol(CCLE)
sst <- 1e-6
nullS <- 10000

#
d <- 100
svE <- svd(E, nu = d, nv = 0)
svV <- svd(CCLE, nu = 0, nv = d)
GE <- CCLE%*%svE[["u"]][, 1:d]%*%diag(svE[["d"]][1:d])
iGE <- solve((1/sst)*crossprod(GE) + diag(ncol(GE)))
null_fab_dist <- sapply(1:nullS, function(s){
  x <- rnorm(n)
  y <- x / sqrt(sum(x^2))
  -(n/2)*log((1/sst) - ((1/sst)^2)*crossprod(y, GE)%*%iGE%*%crossprod(GE, y))
})

#
df <- ncol(CCLE)
svC <- svd(CCLE, nu = df, nv = 0)
CE <- svC[["u"]]
iCE <- solve(crossprod(CE))
null_f_dist <- qf(seq(1/(nullS), 1, length.out = nullS), df, n-df-1)

#
p_f <- rep(NA, ncol(REP))
p_fab <- rep(NA, ncol(REP))
np_f <- rep(NA, ncol(REP))
np_fab <- rep(NA, ncol(REP))
for(j in 1:ncol(REP)){
  #
  x <- REP[, j]
  y <- x / sqrt(sum(x^2))
  u <- rnorm(n)
  u <- u / sqrt(sum(u^2))
  
  #
  F_FAB <- c(-(n/2)*log((1/sst) - ((1/sst)^2)*crossprod(y, GE)%*%iGE%*%crossprod(GE, y)))
  yhat <- y - CE%*%(iCE%*%crossprod(CE, y))
  F_STAT <- c(((n-df)/df)*(1 - t(y)%*%yhat)/(t(y)%*%yhat))
  p_fab[j] <- mean(F_FAB <= null_fab_dist)
  p_f[j] <- mean(F_STAT <= null_f_dist)
  
  #
  F_FAB <- c(-(n/2)*log((1/sst) - ((1/sst)^2)*crossprod(u, GE)%*%iGE%*%crossprod(GE, u)))
  uhat <- u - CE%*%(iCE%*%crossprod(CE, u))
  F_STAT <- c(((n-df)/df)*(1 - t(u)%*%uhat)/(t(u)%*%uhat))
  np_fab[j] <- mean(F_FAB <= null_fab_dist)
  np_f[j] <- mean(F_STAT <= null_f_dist)
  
  print(j)
}

#
nrp <- 100
null_frp_dist <- qf(seq(1/(nullS), 1, length.out = nullS), d, n-d-1)
p_frp <- matrix(NA, ncol(REP), nrp)
for(s in 1:nrp){
  #
  CE <- svC[["u"]]%*%svd(matrix(rnorm(df*d), nrow = df))[["u"]]
  iCE <- solve(crossprod(CE))
  
  for(j in 1:ncol(REP)){
    #
    x <- REP[, j]
    y <- x / sqrt(sum(x^2))
    
    #
    yhat <- y - CE%*%(iCE%*%crossprod(CE, y))
    F_STAT <- c(((n-d)/d)*(1 - t(y)%*%yhat)/(t(y)%*%yhat))
    
    #
    p_frp[j, s] <- mean(F_STAT <= null_frp_dist)
  
  }
  #
  print(s)
}

#
fdr_frp <- apply(p_frp, 2, function(x){p.adjust(x, method = "BH")})
fdr_f <- p.adjust(p_f, method = "BH")
fdr_fab <- p.adjust(p_fab, method = "BH")

#
pdf(file = "figures/repurposing.pdf", family = "Times", height = 5, width = 5)
plot(1:length(fdr_f), fdr_f[order(fdr_f)], type = 'l', lwd = 2,
     ylim = c(0, 1), xlab = "Rank", ylab = "FDR")
lines(1:length(fdr_fab), fdr_fab[order(fdr_fab)], type = 'l', lwd = 2, col = '#08519c')
for(s in 1:nrp){
  lines(1:length(fdr_frp[, s]), fdr_frp[, s][order(fdr_frp[, s])], type = 'l', lwd = 0.1, col = scales::alpha('lightblue', 0.25))
}
abline(h = 0.1, col = '#a50f15')
legend(x = 900, y = 0.4, legend = c("F test", "FAB test", "Rnd. proj."), col = c("black", '#08519c', scales::alpha('lightblue', 0.25)), lwd = c(2, 2, 1))
dev.off()

pdf(file = "figures/repurposing_null.pdf", family = "Times", height = 5, width = 5)
plot(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(p_f < q)}),
     type = 'l', col = 'black', lwd = 2, ylab = "Fraction p-values below threshold", xlab = "Threshold",
     xlim = c(0, 1), ylim = c(0, 1), main = "F-test")
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(p_fab < q)}),
      type = 'l', col = '#08519c', lwd = 2)
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(np_f < q)}),
     type = 'l', col = 'black', lty = 2, lwd = 1.5)
lines(seq(1/100, 1, length.out = 100), sapply(seq(1/100, 1, length.out = 100), function(q){mean(np_fab < q)}),
      type = 'l', col = '#08519c', lty = 2, lwd = 1.5)
abline(0, 1, col = '#a50f15')
legend(x = 0, y = 0.9, legend = c("Null data", "Real data"), lty = c(2, 1), lwd = c(1.5, 2))
dev.off()

#
fdrs <- c(0.01, 0.05, 0.1, 0.2)
Discov <- matrix(NA, nrow = length(fdrs), ncol = 4)
colnames(Discov) <- c("F only", "FAB only", "Both", "Neither")
for(j in 1:length(fdrs)){
  ft <- fdrs[j]
  F_discov <- which(fdr_f < ft)
  FAB_discov <- which(fdr_fab < ft)
  all_discov <- union(F_discov, FAB_discov)
  both_discov <- intersect(FAB_discov, F_discov)
  Fonly <- setdiff(all_discov, FAB_discov)
  FABonly <- setdiff(all_discov, F_discov)
  Discov[j, ] <- c(length(Fonly), length(FABonly), length(both_discov), length(fdr_f) - length(all_discov))
}

# Sensitivity analysis
sens_res <- matrix(c(sum(fdr_f < 0.01), sum(fdr_fab < 0.01),
                     sum(fdr_f < 0.025), sum(fdr_fab < 0.025),
                     sum(fdr_f < 0.05), sum(fdr_fab < 0.05),
                     sum(fdr_f < 0.1), sum(fdr_fab < 0.1)), ncol = 2, byrow = T)
colnames(sens_res) <- c("RP", "FAB")
row.names(sens_res) <- paste0("FDR = ", c(0.01, 0.025, 0.05, 0.1))
saveRDS(sens_res, file = "results/repurposing_exper100_sens.rds")