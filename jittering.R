#' @param x sample
#' @param s resolution of the measurements
jitter <- function(x, s) {
  x <- sort(x)
  n <- length(x)
  # Searching for intervals [i;j] of tied values
  i <- 1
  while (i <= n) {
    j <- i
    while (j < n && x[j + 1] - x[i] < s / 2) {
      j <- j + 1
    }
    if (i < j && j - i + 1 < n) {
      k <- j - i + 1
      u <- 0:(k - 1) / (k - 1)
      xi <- u - 0.5
      if (i == 1)
        xi <- u / 2
      if (j == n)
        xi <- (u - 1) / 2
      if (i == 1 && j == n)
        xi <- u - 0.5
      x[i:j] <- x[i:j] + xi * s
    
    }
    i <- j + 1
  }
  return(x)
}

# Demo
set.seed(1729)
x <- rnorm(2000)
xr <- round(x, 1)
xj <- jitter(xr, 0.1)
df <- rbind(
  data.frame(type = "(a) Original", x = x),
  data.frame(type = "(b) Rounded", x = xr),
  data.frame(type = "(c) Jittered", x = xj)
)
ggplot(df, aes(x)) +
  facet_wrap(vars(type), nrow = 1) +
  geom_qrdensity() +
  geom_rug(sides = "b")
ggplot(df[df$type != "(b) Rounded",], aes(x, col = type)) +
  geom_qrdensity()
