# Common theme and helpers for EIA Nature Methods figures
# Author: generated analysis pipeline
# Date: 2026-09-14
#
# Font note: this script uses "sans" as the font-family token. On macOS,
# Quartz is remapped so that "sans" -> Arial (regular/bold/italic/bold-italic).
# This keeps PNG (ragg) and PDF (Quartz) rendering consistent and avoids the
# "font category error" that occurs when asking the R Quartz device for an
# unregistered family name such as "Arial" directly. On non-Quartz systems
# the default sans-serif family is used.

# Ensure user library is on the search path (where install_pkgs.R places packages)
rv <- paste0(R.version$major, ".", strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1])
user_lib <- file.path(Sys.getenv("HOME"), "Library/R", rv, "library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

# PKG points at the v27 package root (two levels up from code/R) so that
# figures are written to <package>/figures. Set by each figure script before
# sourcing common.R (base_dir <- .../code/R).
if (!exists("PKG") && exists("base_dir")) PKG <- dirname(dirname(base_dir))

require(ggplot2, quietly = TRUE)
require(ragg, quietly = TRUE)

# Map Quartz sans family to Arial on macOS (best-effort, silently ignored elsewhere)
tryCatch({
  if (capabilities("aqua")) {
    quartzFonts(sans = quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic")))
  }
}, error = function(e) invisible(NULL))

FONT <- "sans"

# Theme tuned to Nature-style: Arial via sans token, clean axes, no grid, black axis lines
theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.title = element_text(colour = "black", size = base_size, face = "bold"),
      plot.title = element_text(colour = "black", size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(colour = "black", size = base_size - 1, hjust = 0.5),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

# Panel label helper (a/b/c/d) placed at top-left inside panel
panel_label <- function(label, x = 0.02, y = 0.98, size = 13, face = "bold") {
  annotation_custom(
    grob = grid::textGrob(label, x = x, y = y, just = c("left", "top"),
                          gp = grid::gpar(fontsize = size, fontface = face, fontfamily = FONT, col = "black")),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
}

# Save helper: PNG + PDF at 300 dpi
save_fig <- function(p, name, w, h, out_dir = file.path(PKG, "figures")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  png_path <- file.path(out_dir, paste0(name, ".png"))
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = w, height = h, units = "in", dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(pdf_path, p, width = w, height = h, units = "in", bg = "white")
  message("Saved ", name, ": ", png_path, " and ", pdf_path)
}
