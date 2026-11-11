library(ggplot2)
# ---------------------------------------------------------------
# Path configuration. Set SALIVA_ROOT (and WORK_DIR) in your environment,
# or edit the two lines below to point at your local data directory.
#   export SALIVA_ROOT=/path/to/your/dedup_rerun
# Figure scripts read inputs under ${SALIVA_ROOT}; no data are bundled here.
# ---------------------------------------------------------------
SALIVA_ROOT <- Sys.getenv("SALIVA_ROOT", unset = "/path/to/dedup_rerun")
WORK_DIR    <- Sys.getenv("WORK_DIR",    unset = "/path/to/workdir")

COL_HG38 <- "#7F7F7F"; COL_KOREF <- "#2C7BB6"
OUT <- paste0(SALIVA_ROOT,"/figures_rev")

# 통일 dedup 값 (표준염색체, n=135,098, CAT _unique 주석)
#   hg38 RefGene: gene body 58.6 / promoter 11.6 / intergenic 41.4
#   CAT 개인주석: gene body 77.7 / promoter 15.3 / intergenic 22.3
anno_df <- data.frame(
  Region     = factor(rep(c("Gene body","Promoter","Intergenic"), each=2),
                       levels=c("Gene body","Promoter","Intergenic")),
  Annotation = factor(rep(c("hg38 annotation","Personal (KOREF-TTAGGA) annotation"), 3),
                      levels=c("hg38 annotation","Personal (KOREF-TTAGGA) annotation")),
  Pct        = c(58.6, 77.7,   11.6, 15.3,   41.4, 22.3)
)

pB <- ggplot(anno_df, aes(x=Region, y=Pct, fill=Annotation)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.6) +
  geom_text(aes(label=paste0(Pct,"%")), position=position_dodge(0.7),
            vjust=-0.3, size=3) +
  scale_fill_manual(values=c("hg38 annotation"=COL_HG38,
                             "Personal (KOREF-TTAGGA) annotation"=COL_KOREF)) +
  scale_y_continuous(limits=c(0,90), expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL, y="DMRs (%)", fill=NULL,
       caption="n = 135,098 dedup KOREF-TTAGGA DMRs (standard chromosomes)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        axis.text=element_text(color="black"),
        legend.position="top",
        legend.text=element_text(size=8),
        plot.caption=element_text(size=7, color="grey40"))

ggsave(file.path(OUT,"FigB_genomic_classification_rev.pdf"), pB, width=120, height=100, units="mm", dpi=300, device=cairo_pdf)
ggsave(file.path(OUT,"FigB_genomic_classification_rev.png"), pB, width=120, height=100, units="mm", dpi=300)
cat("Panel B (3분류, dedup 통일값) saved!\n")
