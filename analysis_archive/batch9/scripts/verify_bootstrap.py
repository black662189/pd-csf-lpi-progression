"""Verify the completed paired bootstrap using independent NumPy calculations."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

B = Path(__file__).resolve().parents[1]
O = B / 'results/classes_v1'
checks = []
def check(name, value):
    assert value, name
    checks.append(name)
def read(name):
    return pd.read_csv(O / (name + '.tsv'), sep='\t')

lock = B / 'protocol/analysis_lock_v1.json'
check('analysis lock SHA256 unchanged', hashlib.sha256(lock.read_bytes()).hexdigest() ==
      (B / 'protocol/analysis_lock_v1.sha256').read_text().strip())
boot, con, pilot, baseline = [read(x) for x in
    ['paired_bootstrap', 'paired_class_contrasts', 'paired_pilot', 'paired_score_baseline']]
check('8000 unique draw score pairs', len(boot) == 8000 and not boot.duplicated(['iter','score']).any())
check('four contrasts', set(con.score) == {'LPI','LPC','LPE','LYSO_ALL'} and len(con) == 4)
check('all bootstrap deltas equal paired coefficient subtraction',
      np.allclose(boot.delta, boot.CSF_beta - boot.plasma_beta, atol=1e-12))
for row in con.itertuples():
    z = boot[boot.score.eq(row.score)]
    check(row.score + ' all 2000 iterations present', set(z.iter) == set(range(1,2001)))
    check(row.score + ' model qualification matches status',
          z.qualified.eq(z.CSF_status.eq('tested') & z.plasma_status.eq('tested')).all())
    check(row.score + ' reported failure fraction and count',
          z.qualified.sum() == row.bootstrap_qualified and
          abs((~z.qualified).mean() - row.failure_fraction) < 1e-12)
    check(row.score + ' <= 5 percent failed', row.failure_fraction <= .05)
    ci = np.quantile(z.loc[z.qualified,'delta'], [.025,.975,.00625,.99375])
    check(row.score + ' percentile intervals reproduced',
          np.allclose(ci, [row.lo95,row.hi95,row.lo_family98_75,row.hi_family98_75], atol=1e-12))
    b = baseline[baseline.score.eq(row.score)].set_index('assay')
    check(row.score + ' paired baseline reproduced',
          np.allclose([row.beta_CSF,row.beta_plasma,row.delta],
          [b.loc['CSF','beta'],b.loc['plasma','beta'],b.loc['CSF','beta']-b.loc['plasma','beta']],atol=1e-12))
    check(row.score + ' simultaneous interval classification',
          row.family_supported == (row.lo_family98_75 > 0 or row.hi_family98_75 < 0))
    check(row.score + ' same 137 people and 923 observations', row.people == 137 and row.rows == 923)
cols = [c for c in pilot.columns if c != 'seconds']
pd.testing.assert_frame_equal(boot.loc[boot.iter.le(20),cols].sort_values(['score','iter']).reset_index(drop=True),
                              pilot[cols].sort_values(['score','iter']).reset_index(drop=True))
checks.append('first twenty draws match pilot exactly')
for name in ['batch9_prepare_pilot_v1','batch9_bootstrap_v1']:
    status=json.loads((B/'jobs'/(name+'.status.json')).read_text())
    check(name+' completed successfully', status['state']=='complete' and status['last_returncode']==0)
out={'state':'passed','checks_passed':len(checks),'checks':checks,
     'replicates_per_contrast':2000,'all_contrast_draws':8000,
     'failed_contrast_draws':int((~boot.qualified).sum())}
(O/'python_bootstrap_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='checks'}))
