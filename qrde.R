library(ggplot2) # Plotting
library(Hmisc)   # Harrell-Davis quantile estimator

StatQDensity <- ggproto("StatQDensity", Stat,
  compute_group = function(data, scales, bincount, Q) {
    # Transforming the input data$x to QRDE-HD
    p <- seq(0, 1, length.out = bincount + 1)
    q <- Q(data$x, p)
    h <- pmax(1 / bincount / (tail(q, -1) - head(q, -1)), 0)
    den_x <- rep(q, each = 2)
    den_y <- c(0, rep(h, each = 2), 0)
    data.frame(x = den_x, y = den_y)
  }
)

#' @param bincount the number of bins in the pseudo-histogram
#' @param Q the target quantile estimator (default: Harrell-Davis)
geom_qrdensity <- function(mapping = NULL, data = NULL,
                           stat = "qdensity", position = "identity",
                           bincount = 1000, Q = hdquantile, ...) {
  layer(
    stat = StatQDensity,
    data = data,
    mapping = mapping,
    geom = GeomLine,
    position = position,
    params = list(bincount = bincount, Q = Q, ...),
  )
}

# Demo
set.seed(42)
x <- numeric(200)
mix <- sample(c(rep(0, 8), 1, 2), 200, TRUE)
x[mix == 0] <- runif(sum(mix == 0), 0, 100)
x[mix == 1] <- rnorm(sum(mix == 1), 45)
x[mix == 2] <- rnorm(sum(mix == 2), 55)
ggplot(data.frame(x), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "transparent", col = "black") +
  geom_density(col = "red") +
  geom_qrdensity(col = "#00AA00", linewidth = 1.2) +
  geom_rug(sides = "b")
