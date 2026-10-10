source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))
# Exact reporting optimization: constant HMM rows have the first-column consensus.
# Nonconstant rows/matrices are delegated to the unmodified inferCNV function.
original_state_consensus=get('.get_state_consensus',asNamespace('infercnv'))
fast_state_consensus=function(cell_group_matrix){
 if(ncol(cell_group_matrix)>0 && all(is.finite(cell_group_matrix)) && all(cell_group_matrix==cell_group_matrix[,1]))return(as.numeric(cell_group_matrix[,1]))
 original_state_consensus(cell_group_matrix)
}
set.seed(20261010)
for(i in 1:30){
 m=matrix(sample(1:6,120,replace=TRUE),nrow=20)
 stopifnot(identical(original_state_consensus(m),fast_state_consensus(m)))
 m=matrix(rep(sample(1:6,20,replace=TRUE),6),nrow=20)
 stopifnot(identical(original_state_consensus(m),fast_state_consensus(m)))
}
assignInNamespace('.get_state_consensus',fast_state_consensus,ns='infercnv')
