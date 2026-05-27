#
library(limma)
library(affy)
library(tidyverse)
library(magrittr)

#
d <- 10

#
targets <- readTargets("data/diff_exp/Targets.txt")
x <- read.ilmn(files="data/diff_exp/probe profile.txt", ctrlfiles="data/diff_exp/control probe profile.txt", other.columns="Detection")
y <- neqc(x)
expressed <- rowSums(y$other$Detection < 0.05) >= 3
y <- y[expressed,]
ct <- factor(targets$CellType)

#
embeds <- reticulate::py_load_object("data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle")

#
p1 <- list()
p2 <- list()
p3 <- list()
#
wes_pal <- wesanderson::wes_palette("Zissou1", 5)
#
contrasts <- matrix(c("ML", "LP", "ML", "MS", "MS", "LP"), nrow = 3)
#
for(k in 1:3){
  Y <- rowMeans(y$E[, ct == contrasts[k, 1]]) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  Y1 <- rowMeans(y$E[, ct == contrasts[k, 2]]) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  S <- apply(y$E[, ct == contrasts[k, 1]], 1, sd) %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  S1 <- apply(y$E[, ct == contrasts[k, 2]], 1, sd)  %>%
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  R <- rep(3, nrow(Y)) %>% 
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  R1 <- rep(3, nrow(Y)) %>% 
    matrix(ncol = 1) %>%
    magrittr::set_rownames(y$genes$SYMBOL) %>%
    magrittr::set_colnames("diff_score")
  #
  
  #
  gns <- readr::read_csv("data/diff_exp/gene_map.csv")$gene
  gns_to_use <- intersect(intersect(gns, names(embeds)), row.names(S))
  E <- sapply(1:length(gns_to_use), function(j){
    embeds[[gns_to_use[j]]]
  }) |>
    t() |>
    magrittr::set_rownames(gns_to_use)
  E <- scale(E)
  svE <- svd(E, nu = d,nv = d)
  E <- svE[["u"]]%*%diag(sqrt(svE[["d"]][1:d]))
  #
  p_vals <- fabp_lin_reg(Y = t(Y[gns_to_use, , drop=F]), 
                         S = t(S[gns_to_use, , drop=F]), 
                         R = t(R[gns_to_use, , drop=F]), V = E, U = matrix(1, nrow = 1, ncol = 1),
                         pool_sampling_var = T,
                         Y1 = t(Y1[gns_to_use, , drop=F]), 
                         S1 = t(S1[gns_to_use, , drop=F]), 
                         R1 = t(R1[gns_to_use, , drop=F]))
  #
  rank_p <- data.frame(p = p_vals$p, fdr = p_vals$fdr_p, rank = rank(p_vals$p), type = "UMPU", observed = p_vals$statistic)
  rank_fabp <- data.frame(p = p_vals$fabp, fdr = p_vals$fdr_fabp, rank = rank(p_vals$fabp), type = "FAB", observed = p_vals$statistic)
  if(k==1){
    p1[[paste0(k, "")]] <- data.frame(FDR = seq(0, 0.2, length.out = 50)) %>%
      dplyr::group_by(FDR) %>%
      dplyr::summarise(UMPU = sum(p_vals$fdr_p < FDR),
                       FAB = sum(p_vals$fdr_fabp < FDR)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(FDR < 0.2) %>%
      reshape2::melt(id.vars = "FDR", variable.name = "Test", value.name = "Discoveries") %>%
      ggplot2::ggplot() +
      geom_line(aes(x = Discoveries, y = FDR, group = Test, colour = Test)) +
      geom_hline(yintercept = 0.1, size = 0.5, colour = "#a50f15", alpha = 0.5) +
      theme_light() +
      scale_colour_manual(values = c(wes_pal[1], wes_pal[4]), name = "Test") +
      labs(x = "Rank", y = "FDR") +
      theme(legend.position = c(0.2, 0.5))
  } else {
    p1[[paste0(k, "")]] <- data.frame(FDR = seq(0, 0.2, length.out = 100)) %>%
      dplyr::group_by(FDR) %>%
      dplyr::summarise(UMPU = sum(p_vals$fdr_p < FDR),
                       FAB = sum(p_vals$fdr_fabp < FDR)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(FDR < 0.2) %>%
      reshape2::melt(id.vars = "FDR", variable.name = "Test", value.name = "Discoveries") %>%
      ggplot2::ggplot() +
      geom_line(aes(x = Discoveries, y = FDR, group = Test, colour = Test)) +
      geom_hline(yintercept = 0.1, size = 0.5, colour = "#a50f15", alpha = 0.5) +
      theme_light() +
      scale_colour_manual(values = c(wes_pal[1], wes_pal[4]), name = "Test") +
      labs(x = "Rank", y = "FDR") +
      guides(colour = F)
  }
  #
  p2[[paste0(k, "")]] <- p_vals %>%
    ggplot2::ggplot() +
    geom_hex(aes(y = predicted, x = observed), bins = 100) +
    theme_light() + 
    scale_fill_viridis_c(option = "magma") +
    guides(fill = F) +
    labs(x = "Observed", y = "Predicted")
  #
}
#
cowplot::plot_grid(cowplot::plot_grid(p1$`1`, p2$`1`, ncol = 1),
                   cowplot::plot_grid(p1$`2`, p2$`2`, ncol = 1),
                   cowplot::plot_grid(p1$`3`, p2$`3`, ncol = 1),
                   ncol = 3)
ggplot2::ggsave("figures/limma_example.pdf", width = 12, height = 8)      

