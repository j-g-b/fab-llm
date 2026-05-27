#
library(reticulate)

embeds <- reticulate::py_load_object("data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle")
cgc_genes <- c(readr::read_tsv("data/archive/gene_names.txt")$row.names.E.)
cgc_genes <- intersect(names(embeds), cgc_genes)

#
cgc_embeds <- sapply(1:length(cgc_genes), function(j){
    embeds[[cgc_genes[j]]]
}) |>
  t() |>
  magrittr::set_rownames(cgc_genes)
cgc_embeds <- scale(cgc_embeds)

#
CCLE <- readr::read_csv("data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
cl_names <- CCLE[[1]]
CCLE <- CCLE[, gsub(" [(].*", "", colnames(CCLE)) %in% cgc_genes]
colnames(CCLE) <- gsub(" [(].*", "", colnames(CCLE))
CCLE <- as.matrix(CCLE)
row.names(CCLE) <- cl_names

#
cl_info <- readr::read_csv("data/Model.csv")
aml_cls <- cl_info[["ModelID"]][cl_info[["DepmapModelType"]] == "AML"]
bll_cls <- cl_info[["ModelID"]][cl_info[["DepmapModelType"]] == "BLL"]

#
delta <- colMeans(CCLE[intersect(aml_cls, row.names(CCLE)), ]) -
         colMeans(CCLE[intersect(bll_cls, row.names(CCLE)), ])
U <- svd(cgc_embeds[colnames(CCLE), ])[["u"]]
proj_delta <- crossprod(U, delta)

S <- 100
rproj_delta <- matrix(NA, S, length(delta))
for(s in 1:S){
  Z <- matrix(rnorm(nrow(U)^2), nrow = nrow(U))
  Uz <- svd(Z)[["u"]]
  rproj_delta[s, ] <- crossprod(Uz, delta)
  print(s)
}

#
pdf(file = "figures/case_study.pdf", width = 10, height = 5, family = "Times")
layout(matrix(c(1, 2), nrow = 1))
plot(cumsum(proj_delta^2)/sum(delta^2), type = 'l', lwd = 1.5, col = '#08519c',
     ylab = "Fraction of total signal norm", xlab = "Projection dimension")
abline(0, 1/length(delta), col ='#a50f15')
rproj_delta_frac <- apply(rproj_delta, 1, function(x){
  cumsum(x^2)/sum(delta^2)
}) |> t()
for(s in 1:S){
  lines(1:length(delta), rproj_delta_frac[s, ],
        lwd = 0.1, col = scales::alpha('black', 0.1))
}
lines(1:length(delta), apply(rproj_delta_frac, 2, function(x){quantile(x, 0.025)}), lty = 2)
lines(1:length(delta), apply(rproj_delta_frac, 2, function(x){quantile(x, 0.975)}), lty = 2)
legend(x = 5, y = 0.9, c("LLM embeddings", "Random proj."), col = c('#08519c', 'black'), lwd = c(1.5, 1))

plot((cumsum(proj_delta^2)/sum(delta^2))[1:200], type = 'l', lwd = 1.5, col = '#08519c',
     ylab = "Fraction of total signal norm", xlab = "Projection dimension",
     main = "Zoom-in on first 200 dimensions")
abline(0, 1/length(delta), col ='#a50f15')
rproj_delta_frac <- apply(rproj_delta, 1, function(x){
  cumsum(x^2)/sum(delta^2)
}) |> t()
for(s in 1:S){
  lines(1:200, rproj_delta_frac[s, 1:200],
        lwd = 0.1, col = scales::alpha('black', 0.1))
}
lines(1:200, apply(rproj_delta_frac[,1:200], 2, function(x){quantile(x, 0.025)}), lty = 2)
lines(1:200, apply(rproj_delta_frac[,1:200], 2, function(x){quantile(x, 0.975)}), lty = 2)
dev.off()