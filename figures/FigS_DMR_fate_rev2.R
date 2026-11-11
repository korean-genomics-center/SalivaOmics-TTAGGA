library(ggplot2); library(patchwork); library(data.table)
# ---------------------------------------------------------------
# Path configuration. Set SALIVA_ROOT (and WORK_DIR) in your environment,
# or edit the two lines below to point at your local data directory.
#   export SALIVA_ROOT=/path/to/your/dedup_rerun
# Figure scripts read inputs under ${SALIVA_ROOT}; no data are bundled here.
# ---------------------------------------------------------------
SALIVA_ROOT <- Sys.getenv("SALIVA_ROOT", unset = "/path/to/dedup_rerun")
WORK_DIR    <- Sys.getenv("WORK_DIR",    unset = "/path/to/workdir")

OUT <- paste0(SALIVA_ROOT,"/figures_rev")
COL_MATCH <- "#009E73"; COL_LOST <- "#D55E00"; COL_RECALL <- "#2C7BB6"
tp <- function() theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(), axis.text=element_text(color="black"),
        plot.tag=element_text(face="bold", size=12))

# 확정값 (서버 재집계, 표준염색체 기준)
HG38_STD <- 127587   # 표준염색체 hg38 DMR (= matched 106,678 + lost 20,930 의 모집단)
MATCHED  <- 106678
LOST     <- 20930
KRECALL  <- 28508

# ── Panel A: Fate of hg38 DMRs ──────────────────────────
fate <- data.table(
  cat = factor(c("Matched in KOREF","Lost in KOREF"),
               levels=c("Lost in KOREF","Matched in KOREF")),
  n   = c(MATCHED, LOST),
  pct = c(MATCHED/HG38_STD*100, LOST/HG38_STD*100))
pA <- ggplot(fate, aes(x="", y=n, fill=cat)) +
  geom_bar(stat="identity", width=0.5) +
  geom_text(aes(label=sprintf("%s\n%.1f%%", format(n,big.mark=","), pct)),
            position=position_stack(vjust=0.5), size=3, color="white", fontface="bold") +
  scale_fill_manual(values=c("Matched in KOREF"=COL_MATCH,"Lost in KOREF"=COL_LOST),
                    breaks=c("Matched in KOREF","Lost in KOREF")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(tag="A", x=sprintf("hg38 DMR\n(standard chromosomes; n=%s)", format(HG38_STD,big.mark=",")),
       y="DMR count", fill="", title="Fate of hg38 DMRs in KOREF") +
  tp() + theme(legend.position="top", plot.title=element_text(size=10, hjust=0.5))

# ── Panel B: DMR composition by reference ───────────────
comp <- data.table(
  ref  = factor(rep(c("hg38","KOREF-TTAGGA"), c(2,2)), levels=c("hg38","KOREF-TTAGGA")),
  part = c("Shared (matched)","hg38-only (lost)","Shared (matched)","KOREF-only (recalled)"),
  n    = c(MATCHED, LOST, MATCHED, KRECALL))
comp[, part := factor(part, levels=c("hg38-only (lost)","Shared (matched)","KOREF-only (recalled)"))]
pB <- ggplot(comp, aes(x=ref, y=n, fill=part)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=format(n,big.mark=",")),
            position=position_stack(vjust=0.5), size=2.9, color="white", fontface="bold") +
  scale_fill_manual(values=c("Shared (matched)"=COL_MATCH,
                             "hg38-only (lost)"=COL_LOST,
                             "KOREF-only (recalled)"=COL_RECALL)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(tag="B", x=NULL, y="DMR count", fill="", title="DMR composition by reference") +
  guides(fill=guide_legend(ncol=1)) +
  tp() + theme(legend.position="right",
               legend.text=element_text(size=7),
               legend.key.size=unit(0.35,"cm"),
               plot.title=element_text(size=10, hjust=0.5))

figS <- (pA | pB) + plot_layout(widths=c(1,1.15))
ggsave(file.path(OUT,"FigS_DMR_fate_rev2.pdf"), figS, width=190, height=95, units="mm", dpi=300, device=cairo_pdf)
ggsave(file.path(OUT,"FigS_DMR_fate_rev2.png"), figS, width=190, height=95, units="mm", dpi=300)
cat("FigS_DMR_fate_rev2 saved (n=127,587)!\n")
