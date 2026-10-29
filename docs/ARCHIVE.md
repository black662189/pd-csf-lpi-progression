# Executed-study source archive

The archive contains 28 original analysis/verification scripts. Private machine paths have been replaced with configuration variables; `archive_provenance.json` records hashes and adaptations. Original private source files and results were not overwritten.

Set `PD_PROJECT_ROOT` to a working copy of the original project layout and `PD_EVIDENCE_ROOT` to the upstream clinical/lipid evidence directory. Optional `PD_R_LIB` and `PD_COLOC_R_LIB` designate private R package libraries. Copy the relevant `batch*/scripts` directories under the configured project root before invoking scripts that resolve their root relative to `__file__`. All outputs should be directed to new versioned directories. These scripts retain project-specific input filenames and are not a second standalone demo.

| Batch | Role and upstream dependency |
| --- | --- |
| 5 | Endpoint-specific 547-slot preparation and screen; requires qualified clinical and molecular long tables and original slot/reference tables |
| 6 | Original LPI feature models, paired contrast and influence analysis; consumes prepared frozen input cache |
| 8 | Clinical subparts and source-field qualification |
| 9 | Fixed lipid scores and paired bootstrap; includes score/bootstrap numerical checks |
| 10 | Measurement, transcriptomic and clinical sensitivity context; corrected background preparation retained |
| 11 | Fixed parent-score construction and conditional models |
| 14–15 | MBOAT7 and LPI-led regional genetic analyses; requires source-specific GWAS/QTL summaries |
| 16 | Model checks, whole-participant bootstrap, external postmortem CSF analysis and diagnostic figures |

Some scripts require separately frozen protocol files and source-derived inputs described in the manuscript supplement. They are not replaced with synthetic values. The portable root-level pipeline is the supported entry point for rerunning the core longitudinal estimators from qualified model-ready data.

## Contrast correction

The public copy of batch16 `diagnostics.R` corrects the historical `grid(t,x)` argument/column collision using `exposure_value` and `data.table::set`. The earlier `export_contrasts_v2.R` is retained as the executed correction source. The portable implementation independently constructs the same contrast and tests its analytic vector; it can refit the models directly rather than requiring undistributed fitted R objects. Semantic model names replace expressions such as `fits[[ids[1]]]` in incremental comparison rows. These adaptations do not change the model formulas, exposure definitions or family thresholds.

## Reproduction status

See `validation/` for actual executed checks. These records distinguish synthetic software testing from fresh refitting with authorized clinical inputs. They do not certify unrestricted raw-data reproduction, biological replication, or independent external validation.
