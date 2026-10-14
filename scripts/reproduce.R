# Execute the coursework, then compare with R's independent GLM implementation.
if (!file.exists("analysis/coursework.R")) {
  stop("Run from the repository root: Rscript scripts/reproduce.R")
}
dir.create("results", showWarnings = FALSE)
pdf("results/coursework-plots.pdf", width = 8, height = 6)
course <- new.env()
tryCatch(sys.source("analysis/coursework.R", envir = course), finally = dev.off())
links <- c("logit", "probit", "cloglog")
coefficients <- list(course$b_logit, course$b_probit, course$b_log)
probabilities <- list(course$p_logit, course$p_probit, course$p_log)
dat <- data.frame(x = course$x, y = course$y, n = course$n)
models <- lapply(links, function(link) {
  glm(cbind(y, n - y) ~ x, data = dat, family = binomial(link = link),
      control = glm.control(epsilon = 1e-12, maxit = 100))
})
for (i in seq_along(links)) {
  stopifnot(models[[i]]$converged)
  stopifnot(max(abs(coef(models[[i]]) - coefficients[[i]])) < 1e-4)
  stopifnot(max(abs(fitted(models[[i]]) - probabilities[[i]])) < 1e-6)
  stopifnot(all(is.finite(probabilities[[i]])),
            all(probabilities[[i]] > 0 & probabilities[[i]] < 1))
}
comparison <- data.frame(
  link = links,
  intercept = vapply(coefficients, function(b) b[1], numeric(1)),
  slope = vapply(coefficients, function(b) b[2], numeric(1)),
  deviance = vapply(probabilities, function(p) course$deviances(dat$y, dat$n, p), numeric(1)),
  pearson = vapply(probabilities, function(p) course$chi2_pearson(dat$y, dat$n, p), numeric(1))
)
stopifnot(max(abs(comparison$deviance - vapply(models, deviance, numeric(1)))) < 1e-6)
write.csv(comparison, "results/comparison.csv", row.names = FALSE)
predictions <- data.frame(dat, observed = dat$y / dat$n)
for (i in seq_along(links)) predictions[[links[i]]] <- as.numeric(probabilities[[i]])
write.csv(predictions, "results/predictions.csv", row.names = FALSE)

# A common scale makes differences between the three links visible.
png("results/dose-response.png", width = 1440, height = 900, res = 150)
par(mar = c(4.5, 4.5, 3, 1), bg = "#faf8f2", fg = "#102d2b")
plot(dat$x, dat$y / dat$n, pch = 19, cex = 1.2, ylim = c(0, 1),
     xlab = "Dose", ylab = "Observed / fitted proportion",
     main = "Binomial GLMs: three links, one dose-response study")
grid(col = "#e3e3dc")
colors <- c("#176076", "#df5b2d", "#7c426b")
dense <- data.frame(x = seq(min(dat$x), max(dat$x), length.out = 300))
for (i in seq_along(links)) {
  lines(dense$x, predict(models[[i]], dense, type = "response"),
        col = colors[i], lty = i, lwd = 2.5)
}
points(dat$x, dat$y / dat$n, pch = 19, cex = 1.2)
legend("topleft", c("Observed", links), pch = c(19, NA, NA, NA),
       lty = c(NA, 1:3), lwd = 2, col = c("#102d2b", colors), bty = "n")
dev.off()
print(comparison, row.names = FALSE)
cat("PASS: all three manual fits agree with stats::glm.\n")
