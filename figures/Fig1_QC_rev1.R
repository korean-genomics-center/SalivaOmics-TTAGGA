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

cat("Loading beta matrix...\n")
HG38 <- paste0(SALIVA_ROOT,"/matrices/hg38")
mat <- fread(file.path(HG38, "final_hg38_union_CpG_matrix_with_NA_and_coverage.txt"),
             select=c("Blood_hg38_beta","Saliva_hg38_beta"))
set.seed(42)
idx <- sample(nrow(mat), 500000)
beta_df <- data.frame(
  beta   = c(mat$Blood_hg38_beta[idx], mat$Saliva_hg38_beta[idx]),
  tissue = c(rep("Blood",500000), rep("Saliva",500000))
)
beta_df$beta <- suppressWarnings(as.numeric(beta_df$beta)); beta_df <- beta_df[!is.na(beta_df$beta),]
cat("Done!\n")

# ── Panel A ───────────────────────────────────────
pA <- ggplot() +
  annotate("rect", xmin=0.15, xmax=0.85, ymin=0.82, ymax=0.99,
           fill="#F5F7FA", color="#1C2B3A", linewidth=0.7) +
  annotate("text", x=0.5, y=0.905,
           label="Single Individual (KOREF1)",
           size=2.8, fontface="bold", color="#1C2B3A") +
  annotate("rect", xmin=0.02, xmax=0.44, ymin=0.53, ymax=0.75,
           fill="#FFF0EB", color=COL_BLOOD, linewidth=0.7) +
  annotate("text", x=0.23, y=0.70,
           label="Blood", size=2.8, fontface="bold", color=COL_BLOOD) +
  annotate("text", x=0.23, y=0.62,
           label="276 GB BAM\n59.5× cov",
           size=2.3, color=COL_BLOOD, lineheight=0.9) +
  annotate("rect", xmin=0.56, xmax=0.98, ymin=0.53, ymax=0.75,
           fill="#EBF4FF", color=COL_SALIVA, linewidth=0.7) +
  annotate("text", x=0.77, y=0.70,
           label="Saliva", size=2.8, fontface="bold", color=COL_SALIVA) +
  annotate("text", x=0.77, y=0.62,
           label="191 GB BAM\n32.6× cov",
           size=2.3, color=COL_SALIVA, lineheight=0.9) +
  annotate("segment", x=0.38, xend=0.23, y=0.82, yend=0.75,
           arrow=arrow(length=unit(0.1,"cm")), color="grey50", linewidth=0.5) +
  annotate("segment", x=0.62, xend=0.77, y=0.82, yend=0.75,
           arrow=arrow(length=unit(0.1,"cm")), color="grey50", linewidth=0.5) +
  annotate("rect", xmin=0.05, xmax=0.95, ymin=0.05, ymax=0.42,
           fill="#EBFAF4", color=COL_SHARED, linewidth=0.7) +
  annotate("text", x=0.5, y=0.235,
           label="EM-seq (Blood 59.5×, Saliva 32.6× mean CpG)\nhg38 + KOREF-TTAGGA reference\nDSS Beta-binomial DMR calling",
           size=2.5, color="#1C2B3A", lineheight=1.1) +
  annotate("segment", x=0.23, xend=0.33, y=0.53, yend=0.42,
           arrow=arrow(length=unit(0.1,"cm")), color="grey50", linewidth=0.5) +
  annotate("segment", x=0.77, xend=0.67, y=0.53, yend=0.42,
           arrow=arrow(length=unit(0.1,"cm")), color="grey50", linewidth=0.5) +
  xlim(0,1) + ylim(0,1) +
  labs(tag="A") +
  theme_void() +
  theme(plot.tag=element_text(face="bold", size=12),
        plot.margin=margin(4,4,4,4))

# ── Panel B — 왼쪽 상단, 왼쪽 앵커 ──────────────
pB <- ggplot(beta_df, aes(x=beta, color=tissue, fill=tissue)) +
  geom_density(alpha=0.25, linewidth=0.8) +
  scale_color_manual(values=c("Blood"=COL_BLOOD,"Saliva"=COL_SALIVA)) +
  scale_fill_manual(values=c("Blood"=COL_BLOOD,"Saliva"=COL_SALIVA)) +
  scale_x_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
  labs(tag="B", x="DNA methylation (β)", y="Density",
       color="", fill="") +
  theme_paper() +
  theme(
    legend.position      = c(0.02, 0.98),  # 왼쪽 상단 패널 내부
    legend.justification = c(0, 1),         # 왼쪽 위 앵커 → 안으로 들어옴
    legend.direction     = "horizontal",
    legend.key.size      = unit(0.25,"cm"),
    legend.key.width     = unit(0.35,"cm"),
    legend.text          = element_text(size=7.5),
    legend.background    = element_rect(fill=alpha("white",0.85),
                                        color="grey80", linewidth=0.3),
    legend.margin        = margin(1,3,1,3)
  )

# ── Panel C — 오른쪽 상단, 오른쪽 앵커 ──────────
cov_df <- data.frame(
  Tissue    = rep(c("Blood","Saliva"), each=3),
  Threshold = rep(c("cov ≥ 1×","cov ≥ 10×","cov ≥ 20×"), 2),
  Pct       = c(97.3, 78.0, 47.5, 97.0, 81.8, 47.7)
)
cov_df$Threshold <- factor(cov_df$Threshold,
  levels=c("cov ≥ 1×","cov ≥ 10×","cov ≥ 20×"))

pC <- ggplot(cov_df, aes(x=Threshold, y=Pct, fill=Tissue)) +
  geom_bar(stat="identity", position="dodge", width=0.65) +
  scale_fill_manual(values=c("Blood"=COL_BLOOD,"Saliva"=COL_SALIVA)) +
  scale_y_continuous(limits=c(0,118), breaks=seq(0,100,25)) +
  geom_text(aes(label=sprintf("%.1f%%",Pct)),
            position=position_dodge(0.65),
            vjust=-0.25, size=2.2) +
  labs(tag="C", x="", y="CpG sites (%)", fill="") +
  theme_paper() +
  theme(
    legend.position      = c(0.98, 0.98),  # 오른쪽 상단 패널 내부
    legend.justification = c(1, 1),         # 오른쪽 위 앵커 → 안으로 들어옴
    legend.direction     = "horizontal",
    legend.key.size      = unit(0.25,"cm"),
    legend.key.width     = unit(0.35,"cm"),
    legend.text          = element_text(size=7.5),
    legend.background    = element_rect(fill=alpha("white",0.85),
                                        color="grey80", linewidth=0.3),
    legend.margin        = margin(1,3,1,3)
  )

# ── Panel D ───────────────────────────────────────
venn_df <- data.frame(
  Category = factor(c("Blood-only","Shared","Saliva-only"),
                    levels=c("Blood-only","Shared","Saliva-only")),
  Count    = c(1930332, 55874385, 1676564),
  Pct      = c(3.2, 93.9, 2.8)
)

pD <- ggplot(venn_df, aes(x=Category, y=Count/1e6, fill=Category)) +
  geom_bar(stat="identity", width=0.6) +
  scale_fill_manual(values=c(
    "Blood-only"  = COL_BLOOD,
    "Shared"      = COL_SHARED,
    "Saliva-only" = COL_SALIVA)) +
  geom_text(aes(label=paste0(round(Count/1e6,1),"M\n(",Pct,"%)")),
            vjust=-0.3, size=2.8) +
  scale_y_continuous(limits=c(0,65)) +
  labs(tag="D", x="", y="CpG sites (millions)") +
  theme_paper() +
  theme(legend.position="none")

# ── Combine ───────────────────────────────────────
fig1 <- (pA | pB) / (pC | pD) +
  plot_layout(heights=c(1,1))

ggsave(file.path(OUTDIR, "Fig1_QC_rev1.pdf"), fig1, width=170, height=160, units="mm", dpi=300, device=cairo_pdf)
ggsave(file.path(OUTDIR, "Fig1_QC_rev1.png"),
       fig1, width=170, height=160, units="mm", dpi=300)

cat("Fig1_rev1 saved!\n")
