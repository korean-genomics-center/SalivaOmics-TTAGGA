# Shared palette and theme for SalivaOmics manuscript figures.
# Colour-blind-safe palette: Wong (2011) Nature Methods 8:441.
# Values are those used for the published figures; do not change.

COL_BLOOD  <- "#D55E00"   # blood
COL_SALIVA <- "#0072B2"   # saliva
COL_SHARED <- "#009E73"   # shared
COL_HG38   <- "#7F7F7F"   # hg38 reference
COL_KOREF  <- "#2C7BB6"   # KOREF-TTAGGA reference
COL_MATCH  <- "#009E73"
COL_LOST   <- "#D55E00"
COL_RECALL <- "#2C7BB6"

theme_paper <- function() {
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text        = ggplot2::element_text(color = "black"),
    plot.tag         = ggplot2::element_text(face = "bold", size = 12)
  )
}
