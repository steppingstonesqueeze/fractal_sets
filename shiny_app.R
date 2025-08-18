# Shiny app for exploring Mandelbrot and Julia sets 

# app.R — Hardened Mandelbrot & Julia (no NULL/length-0 issues)

suppressPackageStartupMessages(library(shiny))

`%||%` <- function(a, b) if (length(a)) a else b

# ----- Core fractal functions -----
grid_complex <- function(width=800, height=600, center=0+0i, scale=1.4) {
  cx <- Re(center); cy <- Im(center)
  aspect <- width/height
  x <- seq(cx - scale*aspect, cx + scale*aspect, length.out=width)
  y <- seq(cy - scale,        cy + scale,        length.out=height)
  outer(y*1i, x, FUN = "+")
}

escape_time_smooth <- function(z_abs, n, bailout=2.0) {
  n + 1 - log(log(pmax(z_abs, bailout) + 1e-16), base = 2)
}

mandelbrot <- function(width, height, center=0+0i, scale=1.4, max_iter=800, bailout=2.0) {
  C <- grid_complex(width, height, center, scale)
  Z <- matrix(0+0i, nrow = nrow(C), ncol = ncol(C))
  escaped <- matrix(FALSE, nrow = nrow(C), ncol = ncol(C))
  escape_iter <- matrix(max_iter, nrow = nrow(C), ncol = ncol(C))
  z_abs_when_escape <- matrix(0.0, nrow = nrow(C), ncol = ncol(C))
  for (n in 0:(max_iter-1)) {
    Z <- Z^2 + C
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
  mu
}

julia <- function(c=0.355+0.355i, width, height, center=0+0i, scale=1.4, max_iter=800, bailout=2.0) {
  Z <- grid_complex(width, height, center, scale)
  escaped <- matrix(FALSE, nrow = nrow(Z), ncol = ncol(Z))
  escape_iter <- matrix(max_iter, nrow = nrow(Z), ncol = ncol(Z))
  z_abs_when_escape <- matrix(0.0, nrow = nrow(Z), ncol = ncol(Z))
  for (n in 0:(max_iter-1)) {
    Z <- Z^2 + c
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
  mu
}

palette_fractal <- function(n=1024) {
  h <- seq(0, 1, length.out=n)
  s <- rep(1, n)
  v <- (seq(0, 1, length.out=n))^0.7
  hsv(h, s, v)
}

# ----- UI -----
ui <- fluidPage(
  titlePanel("Fractal Explorer (Hardened)"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode", "Mode", c("Mandelbrot", "Julia"), inline = TRUE),
      numericInput("width",  "Width",  900, min = 50, max = 4000, step = 50),
      numericInput("height", "Height", 700, min = 50, max = 3000, step = 50),
      sliderInput("max_iter", "Max Iterations", min = 50, max = 5000, value = 800, step = 50),
      sliderInput("bailout",  "Bailout radius", min = 2, max = 8, value = 2, step = 0.5),
      
      conditionalPanel("input.mode == 'Mandelbrot'",
                       numericInput("cx", "Center Re(c)", 0, step = 0.01),
                       numericInput("cy", "Center Im(c)", 0, step = 0.01),
                       sliderInput("scale_m", "Scale (half-height Im)", min = 0.0005, max = 2, value = 1.4, step = 0.0005)
      ),
      conditionalPanel("input.mode == 'Julia'",
                       numericInput("cr", "c real", 0.355, step = 0.01),
                       numericInput("ci", "c imag", 0.355, step = 0.01),
                       numericInput("jcx", "Viewport center Re(z)", 0, step = 0.01),
                       numericInput("jcy", "Viewport center Im(z)", 0, step = 0.01),
                       sliderInput("scale_j", "Scale (half-height Im)", min = 0.0005, max = 2, value = 1.4, step = 0.0005)
      )
    ),
    mainPanel(
      plotOutput("plot", height = "80vh")
    )
  )
)

# ----- Server -----
server <- function(input, output, session) {
  
  output$plot <- renderPlot({
    # Require core inputs to exist (prevents length-0/NULL)
    req(input$mode, input$width, input$height, input$max_iter, input$bailout)
    
    # Safe pulls with defaults
    w <- as.integer(input$width %||% 900);  if (!is.finite(w) || w < 2) w <- 900
    h <- as.integer(input$height %||% 700); if (!is.finite(h) || h < 2) h <- 700
    it <- as.integer(input$max_iter %||% 800); if (!is.finite(it) || it < 1) it <- 800
    R <- as.numeric(input$bailout %||% 2);    if (!is.finite(R) || R <= 0) R <- 2
    
    pal <- palette_fractal(1024)
    
    if (identical(input$mode, "Mandelbrot")) {
      cx <- as.numeric(input$cx %||% 0); if (!is.finite(cx)) cx <- 0
      cy <- as.numeric(input$cy %||% 0); if (!is.finite(cy)) cy <- 0
      sc <- as.numeric(input$scale_m %||% 1.4); if (!is.finite(sc) || sc <= 0) sc <- 1.4
      
      mu <- mandelbrot(width = w, height = h,
                       center = complex(real = cx, imaginary = cy),
                       scale = sc, max_iter = it, bailout = R)
      
    } else {
      cr <- as.numeric(input$cr %||% 0.355); if (!is.finite(cr)) cr <- 0.355
      ci <- as.numeric(input$ci %||% 0.355); if (!is.finite(ci)) ci <- 0.355
      jcx <- as.numeric(input$jcx %||% 0);   if (!is.finite(jcx)) jcx <- 0
      jcy <- as.numeric(input$jcy %||% 0);   if (!is.finite(jcy)) jcy <- 0
      sc  <- as.numeric(input$scale_j %||% 1.4); if (!is.finite(sc) || sc <= 0) sc <- 1.4
      
      mu <- julia(c = complex(real = cr, imaginary = ci),
                  width = w, height = h,
                  center = complex(real = jcx, imaginary = jcy),
                  scale = sc, max_iter = it, bailout = R)
    }
    
    # Contrast and palette mapping (robust against NAs/Inf)
    finite_vals <- mu[is.finite(mu)]
    if (length(finite_vals)) {
      lo <- as.numeric(stats::quantile(finite_vals, probs = 0.02, names = FALSE))
      hi <- as.numeric(stats::quantile(finite_vals, probs = 0.98, names = FALSE))
      if (!is.finite(lo) || !is.finite(hi) || hi <= lo) { lo <- min(finite_vals); hi <- max(finite_vals) + 1e-9 }
    } else {
      lo <- min(mu, na.rm = TRUE); hi <- max(mu, na.rm = TRUE) + 1e-9
    }
    
    mu_norm <- pmax(pmin(mu, hi), lo)
    idx_vec <- 1L + as.integer((length(pal)-1) * (mu_norm - lo) / (hi - lo + 1e-12))
    idx_vec[!is.finite(idx_vec)] <- 1L
    idx <- matrix(idx_vec, nrow = nrow(mu_norm), ncol = ncol(mu_norm))
    idx[idx < 1L] <- 1L; idx[idx > length(pal)] <- length(pal)
    
    op <- par(mar = c(0,0,0,0)); on.exit(par(op), add = TRUE)
    image(x = 1:ncol(idx), y = 1:nrow(idx), z = t(idx), col = pal,
          axes = FALSE, useRaster = TRUE,
          xlab = "", ylab = "", xaxs = "i", yaxs = "i", ylim = c(nrow(idx), 1))
  })
}

shinyApp(ui, server)
