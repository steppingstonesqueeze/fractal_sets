Mandelbrot & Julia — Python and R
=================================

Files:
- `fractal_sets.py` — Python CLI & library for Mandelbrot/Julia with smooth coloring.
- `fractal_sets.R`  — R functions `mandelbrot()`, `julia()`, and `plot_fractal()`.

Python quickstart:
------------------
pip install numpy matplotlib
python fractal_sets.py --mode mandelbrot --outfile mandel.png
python fractal_sets.py --mode julia --c 0.355+0.355j --outfile julia.png

R quickstart:
-------------
source("fractal_sets.R")
img <- mandelbrot(width=1600, height=1200, center=0+0i, scale=1.5, max_iter=800)
plot_fractal(img, file="mandel_r.png")

img2 <- julia(c=0.355+0.355i, width=1600, height=1200, center=0+0i, scale=1.2, max_iter=800)
plot_fractal(img2, file="julia_r.png")

Notes:
- Bailout radius is set to 2.0 (standard); adjust with `--bailout` or function arg.
- Smooth coloring uses: n + 1 - log2(log |z|).
- Viewport: specify `--center` (complex) and `--scale` (half-height). Aspect ratio follows image size.
- For speed, both implementations use vectorized complex arrays; increase `max_iter` for deep zooms.

Several other additions : can now choose other dimensions to compute Julia + Mandelbrot sets for.
We also hasve a shiny app that can be run using shiny::runApp("shiny_app.R") # use full app path #
