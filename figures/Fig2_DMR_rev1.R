library(ggplot2)
library(patchwork)
library(data.table)

# ---------------------------------------------------------------
# Path configuration. Set SALIVA_ROOT (and WORK_DIR) in your environment,
# or edit the two lines below to point at your local data directory.
#   export SALIVA_ROOT=/path/to/your/dedup_rerun
# Figure scripts read inputs under ${SALIVA_ROOT}; no data are bundled here.
# ---------------------------------------------------------------
SALIVA_ROOT <- Sys.getenv("SALIVA_ROOT", unset = "/path/to/dedup_rerun")
WORK_DIR    <- Sys.getenv("WORK_DIR",    unset = "/path/to/workdir")

OUTDIR <- paste0(SALIVA_ROOT,"/figures_rev")
dir.create(OUTDIR, showWarnings=FALSE)

COL_BLOOD  <- "#D55E00"
COL_SALIVA <- "#0072B2"
COL_SHARED <- "#009E73"

theme_paper <- function() {
  theme_bw(base_size=10) +
  theme(
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(color="black"),
    plot.tag          = element_text(face="bold", size=12)
  )
}

# pos: c(x, y) anchor; just: c(hjust, vjust); bg_alpha: legend background opacity
legend_inset <- function(pos=c(0.98,0.98), just=c(1,1), bg_alpha=0.6) {
  theme(
    legend.position      = pos,
    legend.justification = just,
    legend.direction     = "horizontal",
    legend.key.size      = unit(0.25,"cm"),
    legend.text          = element_text(size=7.5),
    legend.background    = element_rect(fill=alpha("white",bg_alpha),
                                        color="grey80", linewidth=0.3),
    legend.margin        = margin(1,3,1,3)
  )
}

HG38 <- paste0(SALIVA_ROOT,"")
DMR  <- file.path(HG38, "dss/hg38/dmr_results.txt")

cat("Loading DMR results...\n")
dmr <- fread(DMR)
cat("DMR loaded:", nrow(dmr), "rows\n")
print(colnames(dmr))

# DSS callDMR output: chr start end length nCG meanMethy1 meanMethy2 diff.Methy areaStat
# Verified: meanMethy1 = Blood, meanMethy2 = Saliva
# diff.Methy = meanMethy1 - meanMethy2 ; diff > 0 => Blood Hyper
dmr[, diff_methy := diff.Methy]
dmr[, area_stat  := areaStat]
dmr[, dmr_len    := length]
dmr[, direction  := ifelse(diff_methy > 0, "Blood Hyper", "Saliva Hyper")]
dmr[, direction  := factor(direction, levels=c("Blood Hyper","Saliva Hyper"))]

cat("Direction counts:\n")
print(table(dmr$direction))

COLDIR <- c("Blood Hyper"=COL_BLOOD, "Saliva Hyper"=COL_SALIVA)

# Panel A: Volcano plot -> legend top-left, raised above point cloud, more transparent bg
# expand y-axis upward slightly so legend sits in empty headroom
pA <- ggplot(dmr, aes(x=diff_methy,
                      y=log10(abs(area_stat)+1),
                      color=direction)) +
  geom_point(alpha=0.3, size=0.3) +
  scale_color_manual(values=COLDIR) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=expansion(mult=c(0.02, 0.12))) +
  labs(tag="A",
       x="Mean methylation difference (Blood - Saliva)",
       y="log10(|areaStat| + 1)",
       color="") +
  theme_paper() + legend_inset(pos=c(0.02,1.0), just=c(0,1), bg_alpha=0.55) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1)))

# Panel B: DMR length distribution -> legend top-right (empty at high bp)
pB <- ggplot(dmr, aes(x=dmr_len, fill=direction)) +
  geom_histogram(bins=60, alpha=0.7, position="identity") +
  scale_fill_manual(values=COLDIR) +
  scale_x_log10() +
  labs(tag="B",
       x="DMR length (bp, log scale)",
       y="Count", fill="") +
  theme_paper() + legend_inset(pos=c(0.98,0.98), just=c(1,1), bg_alpha=0.6)

# Panel C: DMR summary bar (counts derived from data)
cnt <- as.data.frame(table(dmr$direction))
colnames(cnt) <- c("Category","Count")
cnt$Pct <- round(100 * cnt$Count / sum(cnt$Count), 1)
cnt$Category <- factor(cnt$Category, levels=c("Blood Hyper","Saliva Hyper"))

pC <- ggplot(cnt, aes(x=Category, y=Count/1000, fill=Category)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual(values=COLDIR) +
  geom_text(aes(label=paste0(format(Count, big.mark=","),
                             "\n(", Pct, "%)")),
            vjust=-0.3, size=3) +
  scale_y_continuous(limits=c(0,145)) +
  labs(tag="C", x="", y="DMR count (thousands)") +
  theme_paper() +
  theme(legend.position="none")

# Panel D: Genomic annotation
anno_df <- data.frame(
  Region   = factor(c("Gene body","Intergenic",
                      "Promoter\n(±2kb)","CpG Island"),
                    levels=c("Gene body","Intergenic",
                             "Promoter\n(±2kb)","CpG Island")),
  Count    = c(74695, 53247, 14806, 3627),
  Pct      = c(58.4, 41.6, 11.6, 2.8)
)

pD <- ggplot(anno_df, aes(x=Region, y=Pct)) +
  geom_bar(stat="identity", fill="#555555", width=0.6) +
  geom_text(aes(label=paste0(Pct,"%\n(n=",
                format(Count, big.mark=","),")")),
            vjust=-0.3, size=2.6) +
  scale_y_continuous(limits=c(0,70)) +
  labs(tag="D", x="", y="DMRs (%)") +
  theme_paper()

# Combine
fig2 <- (pA | pB) / (pC | pD) +
  plot_layout(heights=c(1,1))

ggsave(file.path(OUTDIR, "Fig2_DMR_rev1.pdf"),
       fig2, width=170, height=160, units="mm", dpi=300)
ggsave(file.path(OUTDIR, "Fig2_DMR_rev1.png"),
       fig2, width=170, height=160, units="mm", dpi=300)

cat("Fig2_rev1 saved!\n")
