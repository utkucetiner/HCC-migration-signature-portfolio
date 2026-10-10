source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))

options(scipen=100)
suppressPackageStartupMessages(library(infercnv))
for(p in commandArgs(TRUE)){
 set.seed(20261010);d=file.path('work/cnv',p);counts=as(Matrix::readMM(gzfile(file.path(d,'counts.mtx.gz'))),'CsparseMatrix');rownames(counts)=readLines(file.path(d,'genes.txt'));colnames(counts)=readLines(file.path(d,'cells.txt'));ann=read.delim(file.path(d,'annotations.tsv'),header=FALSE,row.names=1)
 print(table(ann$V2));ref=grep('^hep_',unique(ann$V2),value=TRUE);ref=ref[!grepl('Tumor',ref)]
 if(length(ref)!=1)stop('cannot identify adjacent hepatocyte reference')
 obj=CreateInfercnvObject(raw_counts_matrix=counts,gene_order_file='work/cnv/gene_order.tsv',annotations_file=ann,ref_group_names=ref,min_max_counts_per_cell=c(100,Inf))
 out=file.path(d,'hepatocyte_reference');obj=infercnv::run(obj,cutoff=0.1,out_dir=out,cluster_by_groups=TRUE,denoise=TRUE,HMM=FALSE,num_threads=2,no_plot=TRUE,no_prelim_plot=TRUE,save_rds=TRUE)
 saveRDS(obj,file.path(out,'final_infercnv.rds'));writeLines(p,file.path(out,'COMPLETED'))
}
