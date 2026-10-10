source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))

options(scipen=100)
suppressPackageStartupMessages(library(infercnv))
for(p in commandArgs(TRUE)){
 d=file.path('work/cnv',p);if(file.exists(file.path(d,'cell_labels.tsv')) && file.exists(file.path(d,'primary_bins.tsv.gz')) && file.exists(file.path(d,'hepatocyte_reference_bins.tsv.gz')))next;pd=file.path(d,'primary');sd=file.path(d,'hepatocyte_reference')
 if(!file.exists(file.path(pd,'COMPLETED'))||!file.exists(file.path(sd,'COMPLETED')))next
 primary=readRDS(file.path(pd,'final_infercnv.rds'));second=readRDS(file.path(sd,'final_infercnv.rds'))
 regfile=list.files(pd,pattern='17_HMM_pred.*pred_cnv_regions.dat$',full.names=TRUE);mapfile=list.files(pd,pattern='17_HMM_pred.*cell_groupings$',full.names=TRUE)
 stopifnot(length(regfile)==1,length(mapfile)==1)
 reg=read.delim(regfile);mapping=read.delim(mapfile);ann=read.delim(file.path(d,'annotations.tsv'),header=FALSE,col.names=c('cell','annotation'))
 heps=ann$cell[grepl('^hep_',ann$annotation)];stopifnot(all(heps %in% colnames(primary@expr.data)),all(heps %in% colnames(second@expr.data)));refindices=unlist(second@reference_grouped_cell_indices);go=second@gene_order;E=second@expr.data
 group=setNames(mapping$cell_group_name,mapping$cell);lab=data.frame(cell=heps,patient=p,tissue=ifelse(grepl(paste0(p,'T_'),heps),'Tumor','Adjacent'),cnv_group=unname(group[heps]),broad_primary_chromosomes=0L,concordant_chromosomes=0L,label='CNV_low_or_undetected',stringsAsFactors=FALSE)
 audit=list()
 for(gr in unique(lab$cnv_group)){
  ix=which(lab$cnv_group==gr);ci=match(lab$cell[ix],colnames(E));rr=reg[reg$cell_group_name==gr & reg$end-reg$start>=1e7 & reg$state!=3,,drop=FALSE]
  lab$broad_primary_chromosomes[ix]=length(unique(rr$chr));if(nrow(rr)==0)next
  concord=lapply(seq_along(ci),function(j)character());valid=0
  for(i in seq_len(nrow(rr))){
   a=rr[i,];gi=which(go$chr==a$chr & go$start>=a$start & go$stop<=a$end)
   if(length(gi)<10)next
   valid=valid+1;refmean=colMeans(E[gi,refindices,drop=FALSE])-1;threshold=max(.05,as.numeric(quantile(abs(refmean),.99)));delta=colMeans(E[gi,ci,drop=FALSE])-1;sign_expected=ifelse(a$state>3,1,-1);pass=sign_expected*delta>threshold
   for(j in which(pass))concord[[j]]=unique(c(concord[[j]],as.character(a$chr)))
   audit[[length(audit)+1]]=data.frame(patient=p,cnv_group=gr,chr=a$chr,start=a$start,end=a$end,state=a$state,genes=length(gi),threshold=threshold,group_cells=length(ix),concordant_cells=sum(pass))
  }
  nchr=lengths(concord);lab$concordant_chromosomes[ix]=nchr;lab$label[ix]='CNV_signal_indeterminate'
  ok=nchr>=2 & lab$broad_primary_chromosomes[ix]>=2
  if(sum(ok)>=20)lab$label[ix[ok]]='RNA_CNV_supported_candidate'
 }
 write.table(lab,file.path(d,'cell_labels.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(length(audit))write.table(do.call(rbind,audit),file.path(d,'segment_concordance.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 # Ten-Mb summaries for visual audit; does not define the HMM states or labels.
 for(mode in c('primary','hepatocyte_reference')){
  ob=if(mode=='primary')primary else second;go=ob@gene_order;bin=paste(go$chr,floor(go$start/1e7),sep=':');groups=split(seq_len(nrow(go)),bin);groups=groups[lengths(groups)>=10];B=do.call(rbind,lapply(groups,function(ix)colMeans(ob@expr.data[ix,,drop=FALSE])));rownames(B)=names(groups)
  con=gzfile(file.path(d,paste0(mode,'_bins.tsv.gz')),'wt');write.table(B,con,sep='\t',quote=FALSE,col.names=NA);close(con)
 }
 cat(p,'labels\n');print(table(lab$tissue,lab$label));rm(primary,second,E);gc()
}
