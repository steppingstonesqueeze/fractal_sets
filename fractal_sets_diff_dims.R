# fractal_sets.R — Mandelbrot & Julia with smooth coloring and complex arithmetic
# Usage examples (from R console):
#   source("fractal_sets.R")
#   img <- mandelbrot(width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800)
#   plot_fractal(img, file="mandelbrot.png")
#
#   img2 <- julia(c=0.355+0.355i, width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800)
#   plot_fractal(img2, file="julia.png")

# Key difference here : you can try other integral dimensions such that z <- z^d + c 
# (default is 2 as always)

plot_dir <- "/Users/gn/work/projects/github_blitz_aug13_18/day5/mandelbrot_julia/"

grid_complex <- function(width=1600, height=1200, center=0+0i, scale=1.5) {
  cx <- Re(center); cy <- Im(center)
  aspect <- width / height
  x <- seq(cx - scale*aspect, cx + scale*aspect, length.out=width)
  y <- seq(cy - scale,        cy + scale,        length.out=height)
  # matrix with row-major y, col-major x; we want z[row, col] = x[col] + i*y[row]
  Z <- outer(y*1i, x, FUN="+")
  return(Z)
}

escape_time_smooth <- function(z_abs, n, bailout=2.0) {
  # n + 1 - log2(log|z|)
  mu <- n + 1 - log(log(pmax(z_abs, bailout) + 1e-16), base=2)
  return(mu)
}

mandelbrot <- function(power_dim=3, width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800, bailout=2.0) {
  C <- grid_complex(width, height, center, scale)
  Z <- matrix(0+0i, nrow=nrow(C), ncol=ncol(C))
  escaped <- matrix(FALSE, nrow=nrow(C), ncol=ncol(C))
  escape_iter <- matrix(max_iter, nrow=nrow(C), ncol=ncol(C))
  z_abs_when_escape <- matrix(0.0, nrow=nrow(C), ncol=ncol(C))
  
  for (n in 0:(max_iter-1)) {
    Z <- Z^power_dim + C
    absZ <- Mod(Z)
    mask <- (!escaped) & (absZ > bailout)
    if (any(mask)) {
      escaped[mask] <- TRUE
      escape_iter[mask] <- n
      z_abs_when_escape[mask] <- absZ[mask]
    }
  }
  mu <- escape_time_smooth(z_abs_when_escape, escape_iter, bailout)
  mu[!escaped] <- max_iter
  return(mu)
}

julia <- function(power_dim=3, c=0.355+0.355i, width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800, bailout=2.0) {
  Z <- grid_complex(width, height, center, scale)
  escaped <- matrix(FALSE, nrow=nrow(Z), ncol=ncol(Z))
  escape_iter <- matrix(max_iter, nrow=nrow(Z), ncol=ncol(Z))
  z_abs_when_escape <- matrix(0.0, nrow=nrow(Z), ncol=ncol(Z))
  
  for (n in 0:(max_iter-1)) {
    Z <- Z^power_dim + c
    absZ <- Mod(Z)
    mask <- (!escaped) & (absZ > bailout)
    if (any(mask)) {
      escaped[mask] <- TRUE
      escape_iter[mask] <- n
      z_abs_when_escape[mask] <- absZ[mask]
    }
  }
  mu <- escape_time_smooth(z_abs_when_escape, escape_iter, bailout)
  mu[!escaped] <- max_iter
  return(mu)
}

# Simple palette builder (no external packages required)
palette_fractal <- function(n=256) {
  # HSV sweep with gentle value curve
  h <- seq(0, 1, length.out=n)
  s <- rep(1, n)
  v <- (seq(0, 1, length.out=n))^0.7
  cols <- hsv(h=h, s=s, v=v)
  return(cols)
}


plot_fractal <- function(mu, file=NULL, palette=palette_fractal(512)) {
  if (!is.matrix(mu)) stop("`mu` must be a matrix (height x width).")
  nr <- nrow(mu); nc <- ncol(mu)
  if (nr <= 0 || nc <= 0) stop("`mu` has zero rows or columns.")
  
  # No transpose/flip gymnastics; let image() handle it with ylim.
  mu_plot <- mu
  
  # Percentile clip for contrast
  finite_vals <- mu_plot[is.finite(mu_plot)]
  if (length(finite_vals) > 0) {
    lo <- as.numeric(quantile(finite_vals, probs=0.02, names=FALSE))
    hi <- as.numeric(quantile(finite_vals, probs=0.98, names=FALSE))
  } else {
    lo <- suppressWarnings(min(mu_plot, na.rm=TRUE))
    hi <- suppressWarnings(max(mu_plot, na.rm=TRUE))
  }
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    lo <- suppressWarnings(min(mu_plot, na.rm=TRUE))
    hi <- suppressWarnings(max(mu_plot, na.rm=TRUE)) + 1e-9
  }
  
  mu_norm <- pmax(pmin(mu_plot, hi), lo)
  idx_vec <- 1L + as.integer((length(palette)-1) * (mu_norm - lo) / (hi - lo + 1e-12))
  idx_vec[!is.finite(idx_vec)] <- 1L
  idx <- matrix(idx_vec, nrow=nr, ncol=nc)
  # clamp to valid palette range
  idx[idx < 1L] <- 1L
  idx[idx > length(palette)] <- length(palette)
  
  
  draw <- function() {
    op <- par(mar=c(0,0,0,0)); on.exit(par(op), add=TRUE)
    # image(x,y,z): nrow(z) must equal length(x); ncol(z) must equal length(y).
    # Our idx is [nr x nc], so pass t(idx) with x=1:nc, y=1:nr.
    image(x=1:nc, y=1:nr, z=t(idx), col=palette, axes=FALSE, useRaster=TRUE,
          xlab="", ylab="", xaxs="i", yaxs="i", ylim=c(nr, 1))
  }
  if (!is.null(file)) {
    png(filename=file, width=nc, height=nr)
    draw(); dev.off()
    message(sprintf("Saved image to %s", file))
  } else {
    draw()
  }
}

# plot both and save to local files - use plot_dir on top to choose where to save 
img <- mandelbrot(power_dim=3, width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800)
fname <- paste0(plot_dir, "mandelbrot_800_dim3.png")
plot_fractal(img, file=fname)

img2 <- julia(power_dim=3, c=0.355+0.355i, width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800)
fname <- paste0(plot_dir, "julia_800_dim3.png")
plot_fractal(img2, file=fname)


