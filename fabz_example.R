#
ss <- 1
n <- 3
alpha <- 0.05
tau <- 1

#
thetas <- seq(-2, 2, length.out = 1000)
zpwr <- pnorm(qnorm(alpha/2), mean = sqrt(n)*thetas/sqrt(ss), sd = 1) + pnorm(qnorm(1-alpha/2), mean = sqrt(n)*thetas/sqrt(ss), sd = 1, lower.tail = F)
ospwr1 <- pnorm(qnorm(alpha), mean = sqrt(n)*thetas/sqrt(ss), sd = 1)
ospwr2 <- pnorm(qnorm(1-alpha), mean = sqrt(n)*thetas/sqrt(ss), sd = 1, lower.tail = F)

#
pdf(file = "figures/pwr_illust.pdf", family = "Times", height = 6, width = 9)
par(mar=c(4, 4, 2, 2))
layout(matrix(c(1, 2, 3, 4), nrow = 2), widths = c(1, 1), heights = c(1, 0.5))
Mus <- c(0.15, 1/3, 0.75)
Bs <- c("#bdd7e7", "#6baed6", "#3182bd", "#08519c")
Rs <- c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
plot(thetas, zpwr, type = 'l', lwd = 1, lty = 3, ylim = c(0, 1),
     xlab = expression(theta), ylab = "Power", col = 'darkgray')
abline(h = alpha, lwd = 0.75)
text(x = 1.5, y = 0.08, expression(paste(alpha,  "= 0.05")), cex=0.9)
lines(thetas, ospwr2, lwd = 1, lty = 1, col = Bs[4])
for(j in 1:length(Mus)){
  mu <- Mus[j]
  ncp <- (sqrt(n)*thetas/sqrt(ss) + mu*sqrt(ss)/(sqrt(n)*tau^2))^2
  fabzpwr <- pchisq(qchisq(1-alpha, df = 1, ncp = (mu*sqrt(ss)/(sqrt(n)*tau^2))^2), df = 1, ncp = ncp, lower.tail = F)
  lines(thetas, fabzpwr, col = Bs[j], lty = 2)
}
segments(x0=rep(-1, 5)+0.3, y0=1.1*c(0.95, 0.9, 0.85, 0.8, 0.75)-0.1,
         x1=rep(-0.5, 5)+0.3, y1=1.1*c(0.95, 0.9, 0.85, 0.8, 0.75)-0.1, 
       col = c("darkgray", Bs), lty = c(3, 2, 2, 2, 1))
text(x=rep(0, 5)+0.3, y=1.1*c(0.95, 0.9, 0.85, 0.8, 0.75)-0.1, c("2-sided", expression(paste(mu,  "= 0.15")), expression(paste(mu,  "= 0.33")), expression(paste(mu,  "= 0.75")), "1-sided"),
     cex = 0.9)
plot(thetas, dnorm(thetas), type = 'l', col = 'darkgray', lty = 3,
     ylab = "Density", xlab = expression(theta))
for(j in 1:length(Mus)){
  mu <- Mus[j]
  lines(thetas, dnorm(thetas, mean = mu, sd = tau), col = Bs[j], lty = 2)
}

#
plot(thetas, zpwr, type = 'l', lwd = 1, lty = 3, ylim = c(0, 1),
     xlab = expression(theta), ylab = "Power", col = 'darkgray')
abline(h = alpha, lwd = 0.75)
text(x = 1.5, y = 0.08, expression(paste(alpha,  "= 0.05")), cex=0.9)
lines(thetas, ospwr1, lwd = 1, lty = 1, col = Rs[4])
for(j in 1:length(Mus)){
  mu <- -Mus[j]
  ncp <- (sqrt(n)*thetas/sqrt(ss) + mu*sqrt(ss)/(sqrt(n)*tau^2))^2
  fabzpwr <- pchisq(qchisq(1-alpha, df = 1, ncp = (mu*sqrt(ss)/(sqrt(n)*tau^2))^2), df = 1, ncp = ncp, lower.tail = F)
  lines(thetas, fabzpwr, col = Rs[j], lty = 2)
}
segments(x0=rep(-1, 5)+0.3, y0=1.1*c(0.95, 0.9, 0.85, 0.8, 0.75)-0.1,
         x1=rep(-0.5, 5)+0.3, y1=1.1*c(0.95, 0.9, 0.85, 0.8, 0.75)-0.1, 
         col = c("darkgray", Rs), lty = c(3, 2, 2, 2, 1))
text(x=rep(0, 5)+0.3, y=1.1*c(0.95, 0.9, 0.85, 0.8, 0.75)-0.1, c("2-sided", expression(paste(mu,  "= -0.15")), expression(paste(mu,  "= -0.33")), expression(paste(mu,  "= -0.75")), "1-sided"),
     cex = 0.9)
plot(thetas, dnorm(thetas), type = 'l', col = 'darkgray', lty = 3,
     ylab = "Density", xlab = expression(theta))
for(j in 1:length(Mus)){
  mu <- -Mus[j]
  lines(thetas, dnorm(thetas, mean = mu, sd = tau), col = Rs[j], lty = 2)
}
dev.off()
