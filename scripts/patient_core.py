"""Exact resampling of fixed patient contrasts; no data access at import."""
import math
import numpy as np
from scipy.stats import t

def compositions(total, length):
    if length==1:
        yield (total,); return
    for k in range(total+1):
        for rest in compositions(total-k,length-1): yield (k,)+rest

def summarize(cohort, method, patients, values, excluded=()):
    x=np.asarray(values,dtype=float); n=len(x)
    assert n>=3 and len(set(patients))==n and np.isfinite(x).all()
    counts=np.array(list(compositions(n,n)))
    weights=np.array([math.factorial(n)/math.prod(math.factorial(int(k)) for k in c)/n**n for c in counts])
    means=counts@x/n; order=np.argsort(means); means=means[order];weights=weights[order]
    assert np.isclose(weights.sum(),1,atol=1e-12)
    average=float(x.mean()); variance=float(np.sum(weights*(means-average)**2))
    assert np.isclose(variance, float(x.var(ddof=0)/n),rtol=1e-10,atol=1e-15)
    cdf=np.cumsum(weights); bounds=[float(means[min(np.searchsorted(cdf,q),len(means)-1)]) for q in [.025,.975]]
    se=float(x.std(ddof=1)/np.sqrt(n)); delta=float(t.ppf(.975,n-1)*se)
    loo=[{'excluded_patient':p,'mean':float(np.delete(x,i).mean()),'change_from_full':float(np.delete(x,i).mean()-average)} for i,p in enumerate(patients)]
    return {'cohort':cohort,'method':method,'n_patients':n,'patients':patients,'patient_differences':x.tolist(),'excluded_unavailable_patients':list(excluded),'mean':average,'exact_patient_bootstrap_percentile_95':bounds,'student_t_95':[average-delta,average+delta],'patient_sd':float(x.std(ddof=1)),'positive_patients':int((x>0).sum()),'negative_patients':int((x<0).sum()),'bootstrap_unique_count_vectors':len(counts),'bootstrap_probability_sum':float(weights.sum()),'leave_one_patient_out':loo,'loo_min':min(r['mean'] for r in loo),'loo_max':max(r['mean'] for r in loo),'loo_sign_flip':any(np.sign(r['mean'])!=np.sign(average) for r in loo)}
