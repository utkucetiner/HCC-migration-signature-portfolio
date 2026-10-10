source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))
options(scipen=100)
suppressPackageStartupMessages(library(infercnv));source(file.path(source_root, 'scripts/cnv_consensus_fastpath.R'))
d='work/cnv/HCC06/primary';f=list.files(d,pattern='^17_HMM_pred.*infercnv_obj$',full.names=TRUE);stopifnot(length(f)==1);ob=readRDS(f)
groups=unlist(ob@tumor_subclusters$subclusters,recursive=FALSE);errors=0
for(g in groups){m=ob@expr.data[,g,drop=FALSE];stopifnot(identical(original_state_consensus(m),fast_state_consensus(m)))}
cat('All',length(groups),'native HCC06 HMM subgroup consensus vectors match exactly\n')
writeLines(paste(length(groups),'native subgroup vectors identical'),'work/cnv_fastpath_verification.txt')
