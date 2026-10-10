source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))

options(scipen=100)
suppressPackageStartupMessages(library(infercnv))
source(file.path(source_root, 'scripts/cnv_consensus_fastpath.R'))
args=commandArgs(TRUE);patients=if(length(args))args else sprintf('HCC%02d',3:10)
for(p in patients){
 if(file.exists(file.path('work/cnv_no_signature',p,'primary','COMPLETED')))next
 set.seed(20261010)
 d=file.path('work/cnv_no_signature',p);con=gzfile(file.path(d,'counts.mtx.gz'));counts=as(Matrix::readMM(con),'CsparseMatrix')
 rownames(counts)=readLines(file.path(d,'genes.txt'));colnames(counts)=readLines(file.path(d,'cells.txt'))
 ann=read.delim(file.path(d,'annotations.tsv'),header=FALSE,row.names=1)
 refs=grep('^ref_',unique(ann$V2),value=TRUE)
 obj=CreateInfercnvObject(raw_counts_matrix=counts,gene_order_file='work/cnv_no_signature/gene_order.tsv',annotations_file=ann,ref_group_names=refs,min_max_counts_per_cell=c(100,Inf))
 out=file.path(d,'primary');dir.create(out,showWarnings=FALSE)
 obj=infercnv::run(obj,cutoff=0.1,out_dir=out,cluster_by_groups=TRUE,denoise=TRUE,HMM=TRUE,HMM_type='i6',HMM_report_by='subcluster',analysis_mode='subclusters',tumor_subcluster_partition_method='leiden',num_threads=2,no_plot=TRUE,write_expr_matrix=FALSE,BayesMaxPNormal=0,save_rds=TRUE,no_prelim_plot=TRUE,inspect_subclusters=FALSE)
 saveRDS(obj,file.path(out,'final_infercnv.rds'))
 writeLines(p,file.path(out,'COMPLETED'))
}
