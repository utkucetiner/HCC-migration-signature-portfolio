source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))

options(scipen=100)
suppressPackageStartupMessages(library(infercnv))
excluded=readLines('work/cnv_no_signature/excluded_signature_genes.txt');stopifnot(length(excluded)==47)
rows=list()
for(p in sprintf('HCC%02d',3:10)){
 expected=readLines(file.path('work/cnv_no_signature',p,'cells.txt'))
 for(mode in c('primary','hepatocyte_reference')){
  ob=readRDS(file.path('work/cnv_no_signature',p,mode,'final_infercnv.rds'))
  stopifnot(!any(rownames(ob@expr.data)%in%excluded),setequal(colnames(ob@expr.data),expected),all(is.finite(ob@expr.data)))
  rows[[length(rows)+1]]=data.frame(patient=p,mode=mode,genes=nrow(ob@expr.data),cells=ncol(ob@expr.data),excluded_found=0)
  rm(ob);gc()
 }
}
write.table(do.call(rbind,rows),'outputs/CNV_47_gen_disarida_nesne_dogrulamasi.tsv',sep='\t',row.names=FALSE,quote=FALSE)
print(do.call(rbind,rows))
