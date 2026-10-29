# Changelog

## 1.0.0

- First public code deposit for the CSF LPI longitudinal study.
- Added portable anchor-to-score and model-ready-to-contrast workflows, explicit formulas, a synthetic end-to-end test and synthetic-only GitHub Actions.
- Included 28 study source scripts with configurable machine paths and source/adaptation hashes.
- Corrected exposure-grid argument shadowing in the public diagnostic source and independently checked the exposure-contrast vectors.
- Refit 10 core/ML models from authorized PPMI model-ready inputs. All 36 contrasts reproduce the earlier corrected results at the manuscript's three-decimal precision; effect/CI, p-value and degrees-of-freedom numerical tolerances are reported separately.
- Distributed code and aggregate verification records only. No participant-level PPMI inputs or manuscript documents are included.

A public code release does not change the study's exploratory status or establish independent cohort replication.
