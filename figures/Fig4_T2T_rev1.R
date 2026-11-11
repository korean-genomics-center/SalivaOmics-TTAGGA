library(ggplot2); library(patchwork); library(data.table); library(scales)
# ---------------------------------------------------------------
# Path configuration. Set SALIVA_ROOT (and WORK_DIR) in your environment,
# or edit the two lines below to point at your local data directory.
#   export SALIVA_ROOT=/path/to/your/dedup_rerun
# Figure scripts read inputs under ${SALIVA_ROOT}; no data are bundled here.
# ---------------------------------------------------------------
SALIVA_ROOT <- Sys.getenv("SALIVA_ROOT", unset = "/path/to/dedup_rerun")
WORK_DIR    <- Sys.getenv("WORK_DIR",    unset = "/path/to/workdir")

OUT<-paste0(SALIVA_ROOT,"/figures_rev")
V  <-paste0(SALIVA_ROOT,"/fig4_rev")
COL_HG38<-"#7F7F7F"; COL_KOREF<-"#2C7BB6"; COL_BLOOD<-"#D55E00"; COL_SALIVA<-"#0072B2"
tp<-function() theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),axis.text=element_text(color="black"),plot.tag=element_text(face="bold",size=12))
lh<-function() theme(legend.position=c(0.98,0.98),legend.justification=c(1,1),legend.direction="horizontal",legend.key.size=unit(0.35,"cm"),legend.text=element_text(size=7.5),legend.title=element_text(size=8),legend.background=element_rect(fill=alpha("white",0.6),color="grey80",linewidth=0.3))

# A: dedup CpG 수
cpg<-data.frame(Reference=factor(rep(c("hg38","KOREF-TTAGGA"),each=2),levels=c("hg38","KOREF-TTAGGA")),
  Metric=rep(c("Shared CpG (M)","Union CpG (M)"),2), Value=c(55.9,59.5,61.1,65.8))
pA<-ggplot(cpg,aes(Metric,Value,fill=Reference))+geom_bar(stat="identity",position=position_dodge(0.7),width=0.6)+
  geom_text(aes(label=Value),position=position_dodge(0.7),vjust=-0.3,size=2.8)+
  scale_fill_manual(values=c("hg38"=COL_HG38,"KOREF-TTAGGA"=COL_KOREF))+
  scale_y_continuous(limits=c(0,80),expand=expansion(mult=c(0,0.02)))+
  labs(tag="A",x=NULL,y="CpG count (millions)",fill="Reference")+tp()+lh()

# B: Genomic classification (gene body/promoter/intergenic, hg38 vs CAT) — dedup 통일값 n=135,098
anno_df<-data.frame(
  Region=factor(rep(c("Gene body","Promoter","Intergenic"),each=2),levels=c("Gene body","Promoter","Intergenic")),
  Annotation=factor(rep(c("hg38 anno","KOREF-TTAGGA anno"),3),levels=c("hg38 anno","KOREF-TTAGGA anno")),
  Pct=c(58.6,77.7, 11.6,15.3, 41.4,22.3))
pB<-ggplot(anno_df,aes(Region,Pct,fill=Annotation))+
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.6)+
  geom_text(aes(label=paste0(Pct,"%")),position=position_dodge(0.7),vjust=-0.3,size=2.6)+
  scale_fill_manual(values=c("hg38 anno"=COL_HG38,"KOREF-TTAGGA anno"=COL_KOREF))+
  scale_y_continuous(limits=c(0,90),expand=expansion(mult=c(0,0.02)))+
  labs(tag="B",x=NULL,y="DMRs (%)",fill=NULL)+tp()+
  theme(legend.position="top",legend.text=element_text(size=7),legend.title=element_blank())

# C(기존 D): EWAS top genes
rt<-function(path,topn=12){d<-fread(path);d<-d[1:min(topn,.N)];d[,gene:=factor(gene,levels=rev(gene))];d}
ewb<-rt(file.path(V,"blood_ewas_top_genes.txt")); ews<-rt(file.path(V,"saliva_ewas_top_genes.txt"))
pCb<-ggplot(ewb,aes(n_cpg,gene))+geom_bar(stat="identity",width=0.7,fill=COL_BLOOD)+
  labs(x="EWAS associations",y=NULL,subtitle="Blood (8.8% DMR overlap)")+tp()+
  theme(axis.text.y=element_text(size=7.5),plot.subtitle=element_text(face="bold",size=9,hjust=0.5))
pCs<-ggplot(ews,aes(n_cpg,gene))+geom_bar(stat="identity",width=0.7,fill=COL_SALIVA)+
  labs(x="EWAS associations",y=NULL,subtitle="Saliva (8.0% DMR overlap)")+tp()+
  theme(axis.text.y=element_text(size=7.5),plot.subtitle=element_text(face="bold",size=9,hjust=0.5))
pC<-wrap_elements((pCb|pCs))

fig4<-((pA|pB)/pC)+plot_layout(heights=c(0.9,1.1))+
  plot_annotation(tag_levels=list(c("A","B","C")))&theme(plot.tag=element_text(face="bold",size=12))
ggsave(file.path(OUT,"Fig4_T2T_rev1.pdf"),fig4,width=190,height=170,units="mm",dpi=300)
ggsave(file.path(OUT,"Fig4_T2T_rev1.png"),fig4,width=190,height=170,units="mm",dpi=300)
cat("Fig4_rev1 saved (A+B+C, Panel C[GO] 제거)!\n")
