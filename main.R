# Libraries --------------------------------------------------------------------

## TeX+markdown+knitr
library(rmarkdown)
library(tinytex)
library(knitr)
library(kableExtra)
library(bookdown)

## Essential
library(tidyverse)

## Plotting extra
library(Matrix)
library(ggpubr)
library(gridExtra)
library(latex2exp)

## Reference implementations
source("qrde.R")
source("jittering.R")

# Helpers ----------------------------------------------------------------------

## A color palette adopted for color-blind people based on https://jfly.uni-koeln.de/color/
cbp <- list(
  red = "#D55E00", blue = "#56B4E9", green = "#009E73", orange = "#E69F00",
  navy = "#0072B2", pink = "#CC79A7", yellow = "#F0E442", grey = "#999999"
)
cbp$values <- unname(unlist(cbp))

# Setup ------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, fig.align = "center",
                      fig.pos = "H", fig.height = 6)
options(knitr.kable.NA = '-', scipen = 999)
theme_set(theme_bw())

# Functions --------------------------------------------------------------------
getBetaHdi <- function(a, b, width) {
  eps <- 1e-9
  if (a < 1 + eps & b < 1 + eps) # Degenerate case
    return(c(NA, NA))
  if (a < 1 + eps & b > 1) # Left border case
    return(c(0, width))
  if (a > 1 & b < 1 + eps) # Right border case
    return(c(1 - width, 1))
  if (width > 1 - eps)
    return(c(0, 1))
  
  # Middle case
  mode <- (a - 1) / (a + b - 2)
  pdf <- function(x) dbeta(x, a, b)
  
  l <- uniroot(
    f = function(x) pdf(x) - pdf(x + width),
    lower = max(0, mode - width),
    upper = min(mode, 1 - width),
    tol = 1e-9
  )$root
  r <- l + width
  return(c(l, r))
}

thdquantile <- function(x, probs, width = 1 / sqrt(length(x))) sapply(probs, function(p) {
  n <- length(x)
  if (n == 0) return(NA)
  if (n == 1) return(x)
  x <- sort(x)
  a <- (n + 1) * p
  b <- (n + 1) * (1 - p)
  hdi <- getBetaHdi(a, b, width)
  hdiCdf <- pbeta(hdi, a, b)
  cdf <- function(xs) {
    xs[xs <= hdi[1]] <- hdi[1]
    xs[xs >= hdi[2]] <- hdi[2]
    (pbeta(xs, a, b) - hdiCdf[1]) / (hdiCdf[2] - hdiCdf[1])
  }
  iL <- floor(hdi[1] * n)
  iR <- ceiling(hdi[2] * n)
  cdfs <- cdf(iL:iR/n)
  W <- tail(cdfs, -1) - head(cdfs, -1)
  sum(x[(iL+1):iR] * W)
})

kdequantile <- function(x, bw, prob) {
  uniroot(function(u) sum(pnorm(u, x, bw)) / length(x) - prob, c(-100, 100))$root
}
bind <- function(fig1, fig2) {
  plot1 <- fig1()
  plot2 <- fig2()
  grid.arrange(plot1, plot2, nrow = 2)
}
build_df <- function(type1, type2, n1, n2, x1, x2) {
  format <- function(let, type, n) {
    paste0("(", let, ") n=", n, "; ", type)
  }
  rbind(
    data.frame(x = x1, type = format("a", type1, n1)),
    data.frame(x = x2, type = format("b", type2, n2))
  )
}

# Figures ----------------------------------------------------------------------
figure_intro <- function() {
  set.seed(7353)
  x <- rnorm(30)
  
  ggplot(data.frame(x), aes(x, 1)) + 
    geom_violin(bw = 0.9, trim = FALSE, draw_quantiles = 0.5,
                col = cbp$blue, fill = "transparent", linewidth = 1.1) +
    geom_boxplot(width = 0.3, col = cbp$red, fill = "transparent", linewidth = 1.1) +
    geom_rug(sides = "b") +
    scale_x_continuous(limits = c(-3.5, 3.5), breaks = -3:3) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    labs(x = "x")
}

figure_densities <- function(df, names = c()) {
  colors <- c(
    "QRDE-HD" = cbp$green,
    "QRDE-THD" = cbp$orange,
    "QRDE-HF7" = cbp$pink,
    "KDE" = cbp$blue,
    "Histogram" = cbp$red
  )
  if (is.numeric(df))
    df <- data.frame(x = df)
  p <- ggplot(df, aes(x))
  for (name in names) {
    if (name == "QRDE-HD")
      p <- p + geom_qrdensity(aes(color = "QRDE-HD"), linewidth = 1.1, bincount = 5000)
    if (name == "QRDE-THD")
      p <- p + geom_qrdensity(aes(color = "QRDE-THD"), linewidth = 1.1,
                              Q = thdquantile, bincount = 5000)
    if (name == "QRDE-HF7")
      p <- p + geom_qrdensity(aes(color = "QRDE-HF7"), Q = quantile, bincount = 5000)
    if (name == "KDE")
      p <- p + geom_density(aes(color = "KDE"))
    if (name == "Histogram")
      p <- p + geom_histogram(aes(color = "Histogram", y = after_stat(density)),
                              fill = "transparent", bins = 30)
  }
  p <- p +
    geom_rug(sides = "b", linewidth = 1.1) +
    scale_color_manual(values = colors) +
    labs(color = "") +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(hjust = 0)
    )
  if (any(names(df) == "type"))
    p <- p + facet_wrap(vars(type), ncol = 1, scales = "free")
  p
}

figure_hf7 <- function() {
  set.seed(1729)
  n1 <- 10
  n2 <- 13
  df <- build_df(
    "Norm(0, 1)",
    "Norm(0, 1)",
    n1, n2, rnorm(n1), rnorm(n2)
  )
  figure_densities(df, c("QRDE-HF7", "QRDE-HD"))
}

figure_hist <- function() {
  set.seed(42)
  n1 <- 150
  x1 <- rnorm(n1)
  
  n2 <- 100
  x2 <- rep(NA, n2)
  shift <- 1000
  mix <- sample(0:2, n2, TRUE)
  x2[mix == 0] <- rnorm(n2, 0)[mix == 0]
  x2[mix == 1] <- rnorm(n2, shift, 0.1)[mix == 1]
  x2[mix == 2] <- rnorm(n2, shift + 10, 0.1)[mix == 2]

  df <- build_df(
    "Norm(0, 1)",
    paste0("Mix of Norm(0, 1), ",
           "Norm(", shift, ", 0.1^2), ",
           "Norm(", shift + 10, ", 0.1^2)"),
    n1, n2, x1, x2)
  figure_densities(df, c("Histogram", "QRDE-HD"))
}

figure_tail <- function() {
  set.seed(1729)
  x1 <- 1:3
  n1 <- length(x1)
  n2 <- 30
  x2 <- runif(30)
  df <- build_df(
    "{1, 2, 3}",
    "Unif(0, 1)",
    n1, n2, x1, x2)
  figure_densities(df, c("KDE", "QRDE-HD"))
}

figure_hist_kde <- function() {
  set.seed(42)
  n1 <- 20
  x1 <- rep(NA, n1)
  mix <- sample(0:1, n1, TRUE)
  x1[mix == 0] <- rnorm(n1, 0)[mix == 0]
  x1[mix == 1] <- rnorm(n1, 10)[mix == 1]
  
  set.seed(5)
  n2 <- 30
  x2 <- rlnorm(n2)

  df <- build_df(
    "Mix of Norm(0, 1) and Norm(10, 1)",
    "LogNorm(0, 1)",
    n1, n2, x1, x2
  )
  figure_densities(df, c("Histogram", "KDE", "QRDE-HD"))
}


figure_bimodality <- function() {
  set.seed(42)
  n1 <- 200
  x1 <- rep(NA, n1)
  shift <- 1000
  mix <- sample(c(rep(0, 8), 1, 2), n1, TRUE)
  x1[mix == 0] <- runif(n1, 0, 100)[mix == 0]
  x1[mix == 1] <- rnorm(n1, 45)[mix == 1]
  x1[mix == 2] <- rnorm(n1, 55)[mix == 2]

  set.seed(42)
  n2 <- 200
  x2 <- rep(NA, n2)
  shift <- 1000
  mix <- sample(c(rep(0, 4), 1, 2), n2, TRUE)
  x2[mix == 0] <- runif(n2, 0, 1000)[mix == 0]
  x2[mix == 1] <- rnorm(n2, 495)[mix == 1]
  x2[mix == 2] <- rnorm(n2, 505)[mix == 2]

  df <- build_df(
    "Mix of Unif(0, 100), Norm(45, 1), Norm(55, 1)",
    "Mix of Unif(0, 1000), Norm(495, 1), Norm(505, 1)",
    n1, n2, x1, x2
  )
  figure_densities(df, c("Histogram", "KDE", "QRDE-HD"))
}

figure_multimodality <- function() {
  modeCount1 <- 10
  n1 <- 1000
  modeCount2 <- 30
  n2 <- 10000
  gen <- function(modeCount, n) {
    set.seed(1729)
    x <- rep(NA, n)
    shift <- 1000
    modes <- 0:(modeCount - 1)
    mix <- sample(modes, n, TRUE)
    for (i in modes)
      x[mix == i] <- rnorm(n, i * 4, 0.3)[mix == i]
    x
  }
  title <- function(modeCount) {
    paste0("Mix of Norm(i * 4, 1) for i=0..", modeCount - 1)
  }
  df <- build_df(
    title(modeCount1), title(modeCount2),
    n1, n2,
    gen(modeCount1, n1), gen(modeCount2, n2)
  )
  figure_densities(df, c("Histogram", "KDE", "QRDE-HD"))
}

figure_overfitting <- function() {
  n <- 500

  set.seed(1729)
  x1 <- rnorm(n)
  set.seed(1729)
  x2 <- runif(500)

  df <- build_df(
    "Norm(0, 1)", "Unif(0, 1)",
    n, n, x1, x2
  )
  figure_densities(df, c("KDE", "QRDE-HD"))
}

figure_robustness <- function() {
  n <- 50
  set.seed(23)
  x1 <- c(runif(n), 10^3)
  set.seed(23)
  x2 <- c(runif(n), 10^9)
  
  df <- build_df(
    "Unif(0, 1) with added {10^3}",
    "Unif(0, 1) with added {10^9}",
    n, n, x1, x2
  )
  figure_densities(df, c("QRDE-HD", "QRDE-THD")) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    geom_rect(aes(xmin = 0.82, xmax = 0.95, ymin = -Inf, ymax = Inf),
              fill = cbp$yellow, alpha = 0.005)
}

figure_discretization <- function() {
  set.seed(1729)
  n <- 2000
  x0 <- round(rnorm(n), 1)

  x <- round(x0, 1)
  df <- data.frame(type = "n = 2000; Norm(0, 1)", n, x)
  p1 <- figure_densities(df, c("KDE", "QRDE-HD")) +
    coord_cartesian(xlim = c(-2.5, 2.5)) +
    scale_x_continuous(breaks = seq(-2.5, 2.5, by = 0.5))
  
  x <- round(x0 / 5, 1)
  df <- data.frame(type = "n = 2000; Norm(0, 0.2^2)", n, x)
  p2 <- figure_densities(df, c("KDE")) +
    coord_cartesian(xlim = c(-0.5, 0.5)) +
    scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1))
  grid.arrange(p1, p2, nrow = 2)
}

figure_jittering <- function() {
  set.seed(1729)
  x <- rnorm(2000)
  xr <- round(x, 1)
  xj <- jitter(xr, 0.1)
  df <- rbind(
    data.frame(type = "(a) Original", x = x, type_d = "(d) Original vs. Jittered"),
    data.frame(type = "(b) Rounded", x = xr, type_d = ""),
    data.frame(type = "(c) Jittered", x = xj, type_d = "(d) Original vs. Jittered")
  )
  p1 <- ggplot(df, aes(x, col = type)) +
    facet_wrap(vars(type), nrow = 1) +
    geom_qrdensity() +
    geom_rug(sides = "b") +
    scale_color_manual(values = c(cbp$orange, cbp$red, cbp$green)) +
    scale_x_continuous(breaks = -3:3) +
    labs(x = "") +
    theme(legend.position = "none")
  p2 <- ggplot(df[df$type != "(b) Rounded",], aes(x, col = type)) +
    facet_wrap(vars(type_d)) +
    geom_qrdensity() +
    scale_color_manual(values = c(cbp$orange, cbp$green)) +
    scale_x_continuous(breaks = -3:3) +
    theme(legend.position = "none")
  grid.arrange(p1, p2, nrow = 2)
}