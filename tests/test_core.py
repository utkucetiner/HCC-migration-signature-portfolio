"""Small independent checks of the public numerical core; no research data."""
import itertools, math, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from patient_core import compositions, summarize
from balance_core import analyze

class PatientResamplingTests(unittest.TestCase):
    def test_eight_patient_count_vectors(self):
        vectors=list(compositions(8,8))
        self.assertEqual(len(vectors),6435)
        self.assertEqual(len(set(vectors)),6435)
        self.assertTrue(all(sum(v)==8 for v in vectors))

    def test_three_patient_distribution_against_ordered_draws(self):
        x=np.array([-0.4,0.1,0.8])
        r=summarize('toy','fixed',['a','b','c'],x)
        draws=np.array([x[list(i)].mean() for i in itertools.product(range(3),repeat=3)])
        expected=np.quantile(draws,[.025,.975],method='inverted_cdf')
        np.testing.assert_allclose(r['exact_patient_bootstrap_percentile_95'],expected,atol=1e-14)
        self.assertAlmostEqual(draws.var(),x.var()/3)
        self.assertAlmostEqual(r['bootstrap_probability_sum'],1)

    def test_patient_influence_is_equal_weighted(self):
        x=np.array([-0.8,0.2,0.3,0.4])
        r=summarize('toy','fixed',list('abcd'),x)
        for i,row in enumerate(r['leave_one_patient_out']):
            self.assertAlmostEqual(row['mean'],sum(v for j,v in enumerate(x) if j!=i)/3)
        self.assertTrue(r['loo_sign_flip'])

    def test_invalid_patient_inputs_are_rejected(self):
        for labels,values in [(['a','a','b'],[1,2,3]),(['a','b','c'],[1,2,float('nan')])]:
            with self.assertRaises(AssertionError): summarize('toy','fixed',labels,values)

class WeightingTests(unittest.TestCase):
    def test_balanced_covariates_recover_known_contrast(self):
        rng=np.random.default_rng(18);C=rng.normal(size=(60,5));C=np.vstack([C,C])
        t=np.repeat([0.,1.],60); y=C@np.array([.1,-.2,.3,.4,-.1])+.4*t
        r,w,p=analyze(C,y,t)
        self.assertAlmostEqual(r['overlap']['difference'],.4,places=10)
        self.assertAlmostEqual(r['entropy']['difference'],.4,places=10)
        np.testing.assert_allclose(p,.5,atol=1e-12)
        self.assertAlmostEqual(r['overlap']['ESS_N'],60)
        self.assertLess(max(abs(v) for v in r['overlap']['SMD'].values()),1e-12)

    def test_outcome_shift_preserves_contrast(self):
        rng=np.random.default_rng(19);C=rng.normal(size=(160,5));t=np.repeat([0.,1.],80)
        y=C[:,0]+.2*t
        a,_,_=analyze(C,y,t);b,_,_=analyze(C,y+7.5,t)
        self.assertAlmostEqual(a['overlap']['difference'],b['overlap']['difference'],places=12)
        self.assertAlmostEqual(a['entropy']['difference'],b['entropy']['difference'],places=12)

if __name__=='__main__': unittest.main()
