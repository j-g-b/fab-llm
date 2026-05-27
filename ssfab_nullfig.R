#
pwrd5n1 <- readRDS("results/ssfab_decay5_noise1.rds")
pwrd5n06 <- readRDS("results/ssfab_decay5_noise0.6.rds")
pwrd5n03 <- readRDS("results/ssfab_decay5_noise0.3.rds")
pwrd20n1 <- readRDS("results/ssfab_decay20_noise1.rds")
pwrd20n06 <- readRDS("results/ssfab_decay20_noise0.6.rds")
pwrd20n03 <- readRDS("results/ssfab_decay20_noise0.3.rds")

#
nsigs <- 15
ntot <- 20

#
apwrd5n1 <- sapply(1:3, function(k){
  pwrd5n1[[2]][[k]][,(nsigs+1):ntot] |> rowMeans()
}) |>
  magrittr::set_colnames(c("r = 0.8", "r = 0.66", "r = 0.5")) |>
  magrittr::set_rownames(c("RP", "AFAB", "SS AFAB"))
apwrd5n06 <- sapply(1:3, function(k){
  pwrd5n06[[2]][[k]][,(nsigs+1):ntot] |> rowMeans()
}) |>
  magrittr::set_colnames(c("r = 0.8", "r = 0.66", "r = 0.5")) |>
  magrittr::set_rownames(c("RP", "AFAB", "SS AFAB"))
apwrd5n03 <- sapply(1:3, function(k){
  pwrd5n03[[2]][[k]][,(nsigs+1):ntot] |> rowMeans()
}) |>
  magrittr::set_colnames(c("r = 0.8", "r = 0.66", "r = 0.5")) |>
  magrittr::set_rownames(c("RP", "AFAB", "SS AFAB"))
apwrd20n1 <- sapply(1:3, function(k){
  pwrd20n1[[2]][[k]][,(nsigs+1):ntot] |> rowMeans()
}) |>
  magrittr::set_colnames(c("r = 0.8", "r = 0.66", "r = 0.5")) |>
  magrittr::set_rownames(c("RP", "AFAB", "SS AFAB"))
apwrd20n06 <- sapply(1:3, function(k){
  pwrd20n06[[2]][[k]][,(nsigs+1):ntot] |> rowMeans()
}) |>
  magrittr::set_colnames(c("r = 0.8", "r = 0.66", "r = 0.5")) |>
  magrittr::set_rownames(c("RP", "AFAB", "SS AFAB"))
apwrd20n03 <- sapply(1:3, function(k){
  pwrd20n03[[2]][[k]][,(nsigs+1):ntot] |> rowMeans()
}) |>
  magrittr::set_colnames(c("r = 0.8", "r = 0.66", "r = 0.5")) |>
  magrittr::set_rownames(c("RP", "AFAB", "SS AFAB"))

#
pdf("figures/sim_study_null.pdf", height = 7, width = 12, family = "Times")
layout(matrix(1:6, nrow = 2, byrow = T))
barplot(c(mean(apwrd5n1[1, ]), mean(apwrd5n1[2, ]), apwrd5n1[3, ]),
        names.arg = c("RP", "AFAB", "SS AFAB \n r = 0.8", "SS AFAB \n r = 0.66", "SS AFAB \n r = 0.5"),
        main = expression(paste("Slow decay, ", gamma == 1)),
        ylim = c(0, 0.1), ylab = "Avg. Power", col = c('black', "#08519c", "#cbc9e2", "#9e9ac8", "#756bb1"), border = NA)
abline(h = 0.05, col = '#a50f15', lty = 2)
barplot(c(mean(apwrd5n06[1, ]), mean(apwrd5n06[2, ]), apwrd5n06[3, ]),
        names.arg = c("RP", "AFAB", "SS AFAB \n r = 0.8", "SS AFAB \n r = 0.66", "SS AFAB \n r = 0.5"),
        main = expression(paste("Slow decay, ", gamma == 2/3)),
        ylim = c(0, 0.1), ylab = "Avg. Power", col = c('black', "#08519c", "#cbc9e2", "#9e9ac8", "#756bb1"), border = NA)
abline(h = 0.05, col = '#a50f15', lty = 2)
barplot(c(mean(apwrd5n03[1, ]), mean(apwrd5n03[2, ]), apwrd5n03[3, ]),
        names.arg = c("RP", "AFAB", "SS AFAB \n r = 0.8", "SS AFAB \n r = 0.66", "SS AFAB \n r = 0.5"),
        main = expression(paste("Slow decay, ", gamma == 1/3)),
        ylim = c(0, 0.1), ylab = "Avg. Power", col = c('black', "#08519c", "#cbc9e2", "#9e9ac8", "#756bb1"), border = NA)
abline(h = 0.05, col = '#a50f15', lty = 2)
barplot(c(mean(apwrd20n1[1, ]), mean(apwrd20n1[2, ]), apwrd20n1[3, ]),
        names.arg = c("RP", "AFAB", "SS AFAB \n r = 0.8", "SS AFAB \n r = 0.66", "SS AFAB \n r = 0.5"),
        main = expression(paste("Fast decay, ", gamma == 1)),
        ylim = c(0, 0.1), ylab = "Avg. Power", col = c('black', "#08519c", "#cbc9e2", "#9e9ac8", "#756bb1"), border = NA)
abline(h = 0.05, col = '#a50f15', lty = 2)
barplot(c(mean(apwrd20n06[1, ]), mean(apwrd20n06[2, ]), apwrd20n06[3, ]),
        names.arg = c("RP", "AFAB", "SS AFAB \n r = 0.8", "SS AFAB \n r = 0.66", "SS AFAB \n r = 0.5"),
        main = expression(paste("Fast decay, ", gamma == 2/3)),
        ylim = c(0, 0.1), ylab = "Avg. Power", col = c('black', "#08519c", "#cbc9e2", "#9e9ac8", "#756bb1"), border = NA)
abline(h = 0.05, col = '#a50f15', lty = 2)
barplot(c(mean(apwrd20n03[1, ]), mean(apwrd20n03[2, ]), apwrd20n03[3, ]),
        names.arg = c("RP", "AFAB", "SS AFAB \n r = 0.8", "SS AFAB \n r = 0.66", "SS AFAB \n r = 0.5"),
        main = expression(paste("Fast decay, ", gamma == 1/3)),
        ylim = c(0, 0.1), ylab = "Avg. Power", col = c('black', "#08519c", "#cbc9e2", "#9e9ac8", "#756bb1"), border = NA)
abline(h = 0.05, col = '#a50f15', lty = 2)
dev.off()

