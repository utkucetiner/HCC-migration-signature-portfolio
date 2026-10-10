source(file.path(Sys.getenv('HCC_SOURCE_ROOT', '.'), 'scripts/runtime.R'))
rows=list()
for(k in c(sprintf('HCC%02d',3:10),paste0('sn',1:3))){
 d=read.delim(paste0('work/overlap_',k,'.tsv'))
 for(mode in c('linear','quadratic')){
  formula=if(mode=='linear')T~x1+x2+x3+x4+x5 else T~x1+x2+x3+x4+x5+I(x1^2)+I(x2^2)+I(x3^2)+I(x4^2)+I(x5^2)
  f=glm(formula,data=d,family=binomial(),control=glm.control(epsilon=1e-12,maxit=200));p=fitted(f);w=ifelse(d$T==1,1-p,p)
  diff=weighted.mean(d$y[d$T==1],w[d$T==1])-weighted.mean(d$y[d$T==0],w[d$T==0])
  rows[[length(rows)+1]]=data.frame(patient=k,model=mode,difference=diff,converged=f$converged)
 }
}
write.table(do.call(rbind,rows),'work/overlap_R_verification.tsv',sep='\t',row.names=FALSE,quote=FALSE)
