"""Within-patient overlap and entropy weighting; no data access at import."""
import numpy as np,json
from scipy.optimize import minimize,linprog
from scipy.special import expit,logsumexp
names=['log_depth','S','G2M','identity','portal_minus_central']
def analyze(C,y,t):
 m=C.mean(0);sd=C.std(0);X=(C-m)/np.where(sd>1e-12,sd,1);D=np.c_[np.ones(len(t)),X]
 def fg(b):
  eta=D@b
  return np.logaddexp(0,eta).sum()-t@eta,D.T@(expit(eta)-t)
 fit=minimize(fg,np.zeros(6),jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':2000})
 p=expit(D@fit.x);w=np.where(t==1,1-p,p)
 def diag(w):
  a=t==0;b=t==1;wa=w[a]/w[a].sum();wb=w[b]/w[b].sum()
  means=wb@C[b]-wa@C[a]; pooled=np.sqrt((C[a].var(0)+C[b].var(0))/2)
  return {'difference':float(wb@y[b]-wa@y[a]),'SMD':dict(zip(names,np.divide(means,pooled,out=np.zeros_like(means),where=pooled>1e-12).tolist())),'ESS_N':float(1/(wa@wa)),'ESS_T':float(1/(wb@wb)),'max_weight_N':float(wa.max()),'max_weight_T':float(wb.max())}
 out={'n_N':int(sum(t==0)),'n_T':int(sum(t==1)),'raw':diag(np.ones(len(t))),'overlap':diag(w),'logistic_gradient_max':float(abs(fg(fit.x)[1]).max()),'logistic_success':bool(fit.success),'PS_N_quantiles':np.quantile(p[t==0],[0,.05,.5,.95,1]).tolist(),'PS_T_quantiles':np.quantile(p[t==1],[0,.05,.5,.95,1]).tolist()}
 # Symmetric entropy balancing to midpoint of unweighted tissue means; constraints on five first moments.
 target=(X[t==0].mean(0)+X[t==1].mean(0))/2;ew=np.zeros(len(t));eb=[]
 for g in [0,1]:
  mask=t==g;A=X[mask]-target
  def dual(l):
   v=A@l;z=np.exp(v-logsumexp(v));return logsumexp(v)-np.log(len(v)),z@A
  feasible=linprog(np.zeros(len(A)),A_eq=np.vstack([np.ones(len(A)),A.T]),b_eq=np.r_[1.,np.zeros(5)],bounds=(0,None),method='highs')
  if not feasible.success:
   eb.append({'group':g,'converged':False,'constraint_error':None,'infeasible':True});continue
  f=minimize(dual,np.zeros(5),jac=True,method='L-BFGS-B',bounds=[(-100,100)]*5,options={'gtol':1e-10,'ftol':1e-14,'maxiter':2000})
  v=A@f.x;ww=np.exp(v-logsumexp(v));ew[mask]=ww
  eb.append({'group':g,'converged':bool(f.success),'constraint_error':float(abs(ww@A).max()),'lambda_norm':float(np.linalg.norm(f.x))})
 out['entropy_diagnostics']=eb
 if all(r['constraint_error'] is not None and r['constraint_error']<1e-6 for r in eb):out['entropy']=diag(ew)
 else:out['entropy_error']='Target moments infeasible or optimizer failed; no outcome estimate reported.'
 return out,w,p
