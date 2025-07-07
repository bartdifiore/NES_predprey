# Scaling/unscaling function to facilitate model convergence
set.seed(13)
x <- rnorm(100)

scaled_x <- scale(x)
unscaled_x <- as.numeric((scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center"))
all.equal(x, unscaled_x)

scaled<- function(x, center, scale){
  (x - center) / scale
}

unscale <- function(scaled_x, center, scale) {
  if (is.null(attr(scaled_x, "scaled:scale")) == F) {
    # (scaled_x * sd) + m
    (scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center")
  }
  if (is.null(attr(scaled_x, "scaled:scale")) == T) {
    (scaled_x * scale) + center
  }
}

unscale_aja <- function(scaled_x, orig_mean, orig_sd) {
  (scaled_x * orig_sd) + orig_mean
}