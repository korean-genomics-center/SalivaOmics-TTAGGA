suppressMessages({library(ggplot2);library(patchwork);library(data.table);library(scales)})
# ---------------------------------------------------------------
# Path configuration. Set SALIVA_ROOT (and WORK_DIR) in your environment,
# or edit the two lines below to point at your local data directory.
#   export SALIVA_ROOT=/path/to/your/dedup_rerun
# Figure scripts read inputs under ${SALIVA_ROOT}; no data are bundled here.
# ---------------------------------------------------------------
SALIVA_ROOT <- Sys.getenv("SALIVA_ROOT", unset = "/path/to/dedup_rerun")
WORK_DIR    <- Sys.getenv("WORK_DIR",    unset = "/path/to/workdir")

OUT<-paste0(SALIVA_ROOT,"/figures_rev")
V<-paste0(SALIVA_ROOT,"/annotation_rev")
M<-paste0(SALIVA_ROOT,"/marker_enrich")
tp<-function() theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),axis.text=element_text(color="black"),plot.tag=element_text(face="bold",size=12))
rn<-function(x)sapply(strsplit(x,"/"),function(v)as.numeric(v[1])/as.numeric(v[2]))
wl<-function(x,w=26)sapply(x,function(s)paste(strwrap(s,width=w),collapse="\n"))
le<-function(path,topn=8,ww=26){ if(!file.exists(path)||length(readLines(path))<2) return(NULL)
  d<-fread(path); if(nrow(d)<1) return(NULL); d<-d[order(p.adjust)][1:min(topn,.N)]
  d[,GeneRatioNum:=rn(GeneRatio)]; d[,Desc:=wl(Description,ww)]; d[,Desc:=factor(Desc,levels=Desc[order(GeneRatioNum)])]
  d[,.(Desc,GeneRatioNum,Count,p.adjust)]}
dp<-function(df,sub,leg=TRUE){ if(is.null(df)) return(ggplot()+annotate("text",0,0,label=paste0(sub,":\nno enriched term\n(after background control)"),size=2.8)+theme_void())
  g<-ggplot(df,aes(GeneRatioNum,Desc,size=Count,color=p.adjust))+geom_point()+
   scale_color_gradient(low="#D7191C",high="#2C7BB6",trans="log10",labels=scientific)+
   scale_size_continuous(range=c(2,6))+scale_x_continuous(n.breaks=4)+
   labs(x="Gene ratio",y=NULL,subtitle=sub,color="p.adjust",size="Count")+tp()+
   theme(axis.text.y=element_text(size=7,lineheight=0.85),plot.subtitle=element_text(face="bold",size=9,hjust=0.5),
     legend.key.size=unit(0.3,"cm"),legend.text=element_text(size=6.5),legend.title=element_text(size=7.5))
  if(!leg) g<-g+theme(legend.position="none"); g}

# A: GO BP (매칭배경)
go_b<-dp(le(file.path(V,"GO_BP_blood_hyper_matchedbg.txt")),"Blood Hyper",FALSE)
go_s<-dp(le(file.path(V,"GO_BP_saliva_hyper_matchedbg.txt")),"Saliva Hyper",TRUE)
pA<-wrap_elements((go_b|go_s))

# B: KEGG (매칭배경)
kg_b<-dp(le(file.path(V,"KEGG_blood_hyper_matchedbg.txt"),6),"Blood Hyper",FALSE)
kg_s<-dp(le(file.path(V,"KEGG_saliva_hyper_matchedbg.txt"),6),"Saliva Hyper",TRUE)
pB<-wrap_elements((kg_b|kg_s))

# C: 마커 중첩 (deconv 대체)  ★tag="C" 제거 — plot_annotation이 C를 붙이므로 중복 방지
r<-fread(file.path(M,"marker_overlap_summary.tsv"))
r[,`:=`(blood_pct=round(blood_hyper/n_marker*100,1),saliva_pct=round(saliva_hyper/n_marker*100,1))]
r[,celltype:=factor(celltype,levels=c("Blood-T","Blood-NK","Blood-B","Blood-Mono+Macro","Blood-Granul","Head-Neck-Ep","Epid-Kerat","Lung-Ep-Bron","Breast-Luminal-Ep","Colon-Ep"))]
mm<-melt(r[,.(celltype,`Blood-hyper`=blood_pct,`Saliva-hyper`=saliva_pct)],id.vars="celltype",variable.name="DMR",value.name="Pct")
pC<-ggplot(mm,aes(celltype,Pct,fill=DMR))+geom_bar(stat="identity",position="dodge",width=0.72)+
  scale_fill_manual(values=c("Blood-hyper"="#D55E00","Saliva-hyper"="#0072B2"))+
  scale_y_continuous(expand=expansion(mult=c(0,0.06)))+
  labs(x=NULL,y="% markers\noverlapping DMRs",fill="DMR",subtitle="DMR overlap with cell-type markers (Loyfer atlas)")+
  tp()+theme(axis.text.x=element_text(angle=35,hjust=1,size=7.5),plot.subtitle=element_text(size=9,hjust=0.5),legend.position="top")
pC<-wrap_elements(pC)

fig3<-(pA/pB/pC)+plot_layout(heights=c(1.1,1.0,1.1))+
  plot_annotation(tag_levels=list(c("A","B","C")))&theme(plot.tag=element_text(face="bold",size=12))
ggsave(file.path(OUT,"Fig3_Biology_rev2.pdf"),fig3,width=210,height=300,units="mm",dpi=300)
ggsave(file.path(OUT,"Fig3_Biology_rev2.png"),fig3,width=210,height=300,units="mm",dpi=300)
cat("Fig3 rev2 saved (tag=C 제거, 라벨 중복 해소)\n")
