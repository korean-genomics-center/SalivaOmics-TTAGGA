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

COL_REGLINE <- "#E41A1C"; COL_IDENT <- "#000000"
COL_NON <- "#2C7BB6"; COL_DMR <- "#D7191C"

theme_paper <- function() {
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        axis.text=element_text(color="black"),
        plot.tag=element_text(face="bold", size=12))
}

SAMP <- paste0(SALIVA_ROOT,"/null_baseline_v2/cpg_sample_dmrid.bed")
d <- fread(SAMP, header=FALSE,
           col.names=c("chr","start","end","Blood","Saliva","dmr_id"))
d <- d[is.finite(Blood) & is.finite(Saliva)]
d[, Region := ifelse(dmr_id!=".", "DMR", "non-DMR")]
d[, Region := factor(Region, levels=c("non-DMR","DMR"))]
rmse <- function(a,b) sqrt(mean((a-b)^2))

# Panel A
fitA <- lm(Blood ~ Saliva, data=d); r2A <- summary(fitA)$r.squared; rA <- cor(d$Saliva,d$Blood)
set.seed(1); dp <- d[sample(.N, min(80000,.N))]
pA <- ggplot(dp, aes(Saliva, Blood)) +
  geom_point(alpha=0.07, size=0.3, color="grey45") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color=COL_IDENT, linewidth=0.6) +
  geom_smooth(method="lm", se=FALSE, color=COL_REGLINE, linewidth=0.9) +
  annotate("text", x=0.04, y=0.97, hjust=0, label=sprintf("r = %.3f\nR² = %.3f", rA, r2A), size=3) +
  annotate("text", x=0.62, y=0.18, hjust=0, color=COL_REGLINE, label="regression", size=2.6) +
  annotate("text", x=0.62, y=0.10, hjust=0, color=COL_IDENT, label="y = x (identity)", size=2.6) +
  coord_equal(xlim=c(0,1), ylim=c(0,1)) +
  labs(tag="A", x="Saliva methylation (β)", y="Blood methylation (β)") + theme_paper()

# Panel B
corr_by <- d[, .(R2=cor(Saliva,Blood)^2, n=.N), by=Region]
pB <- ggplot(corr_by, aes(Region, R2, fill=Region)) +
  geom_bar(stat="identity", width=0.55) +
  geom_text(aes(label=sprintf("R²=%.3f\n(n=%s)", R2, format(n,big.mark=","))), vjust=-0.3, size=2.8) +
  scale_fill_manual(values=c("non-DMR"=COL_NON, "DMR"=COL_DMR)) +
  scale_y_continuous(limits=c(0,1.08), expand=expansion(mult=c(0,0.02))) +
  labs(tag="B", x=NULL, y="Saliva–Blood R²") +
  theme_paper() + theme(legend.position="none")

# Panel C: boxplot, high-contrast traffic-light colors
dmr_ids <- unique(d[Region=="DMR", dmr_id])
N_REP <- 50; test_frac <- 0.3
res <- data.table(); set.seed(2025)
for (rep in 1:N_REP) {
  test_dmr <- sample(dmr_ids, size=round(length(dmr_ids)*test_frac))
  is_test <- d$dmr_id %in% test_dmr
  tr <- d[!is_test]; te <- d[Region=="DMR" & is_test]
  r_u <- rmse(te$Saliva, te$Blood)
  fit <- lm(Blood ~ Saliva, data=tr)
  r_g <- rmse(pmin(pmax(predict(fit, newdata=te),0),1), te$Blood)
  off <- tr[Region=="DMR", mean(Blood - Saliva)]
  r_o <- rmse(pmin(pmax(te$Saliva+off,0),1), te$Blood)
  # null baseline: predict held-out blood by training-set mean blood (ignores saliva)
  mblood <- tr[Region=="DMR", mean(Blood)]
  r_n <- rmse(rep(mblood, nrow(te)), te$Blood)
  res <- rbind(res, data.table(
        Method=c("Uncorrected","Blood-mean\n(null, no saliva)","Global\nlinear","Global\nmean-shift"),
        RMSE=c(r_u,r_n,r_g,r_o)))
}
res[, Method := factor(Method, levels=c("Uncorrected","Blood-mean\n(null, no saliva)","Global\nlinear","Global\nmean-shift"))]
summ <- res[, .(mean=mean(RMSE), sd=sd(RMSE)), by=Method]
cat("Summary:\n"); print(summ)

# traffic-light: worse=red, mid=orange, best=green (high contrast)
cols <- c("Uncorrected"="#D73027",
          "Blood-mean\n(null, no saliva)"="#4575B4",
          "Global\nlinear"="#FC8D59",
          "Global\nmean-shift"="#1A9850")

ymin <- floor((min(res$RMSE)-0.005)*100)/100
ymax <- ceiling((max(res$RMSE)+0.005)*100)/100

pC <- ggplot(res, aes(Method, RMSE, fill=Method)) +
  geom_boxplot(width=0.5, outlier.shape=NA, linewidth=0.45) +
  geom_jitter(width=0.1, size=0.6, alpha=0.4, color="grey20") +
  scale_fill_manual(values=cols) +
  coord_cartesian(ylim=c(ymin, ymax)) +
  scale_y_continuous(breaks=seq(ymin, ymax, by=0.01)) +
  labs(tag="C", x=NULL, y="RMSE within held-out DMRs") +
  theme_paper() +
  theme(legend.position="none",
        axis.text.x=element_text(size=6.5, lineheight=0.85))

fig5 <- (pA | pB | pC) + plot_layout(widths=c(1.0,0.9,1.4))
ggsave(file.path(OUTDIR, "Fig5_Correction_rev1.pdf"), fig5, width=220, height=85, units="mm", dpi=300, device=cairo_pdf)
ggsave(file.path(OUTDIR, "Fig5_Correction_rev1.png"), fig5, width=220, height=85, units="mm", dpi=300)
cat("Fig5_rev1 saved!\n")
