# Input contracts

## `model_ready.csv`

One row is one included participant's subsequent observation. Required columns:

| Column | Definition |
| --- | --- |
| PATNO | Local participant key; character |
| time | Strictly positive years after original molecular index |
| y | Valid Part I total, 0–52 |
| x | Frozen baseline LPI score, constant within person |
| parent | Frozen baseline PI score; NA permitted for score-only participants |
| y0, up0 | Original baseline Part I and OFF/untreated Part III |
| age0 | Original baseline age divided by 10 |
| sex | Original categorical sex coding |
| duration0 | Disease duration in years divided by 5 |
| ledd0 | log(1 + original baseline LEDD) |
| genetic_group | Original genetic/sporadic category |

Participant/time keys must be unique and each person must have at least two distinct observations. All baseline fields must remain constant within person. No missing value imputation or outcome rescaling is performed. The score analysis uses all supplied eligible rows; parent-conditioned and ML comparisons use the common subset with finite parent scores.

## `anchors.csv` and `followup.csv`

The anchor table has one row per **already clinically qualified original molecular index**, including participants who later lack sufficient outcome observations. Required fields are PATNO, LPI16, LPI18, PI1–PI6, y0, up0, age, sex, duration, ledd and genetic_group. Lipid columns contain positive released relative signals on the pre-log2 ratio scale. PI1–PI6 must correspond to the manuscript's fixed species, not newly selected species. `age`, `duration` and `ledd` are untransformed baseline values.

The follow-up table contains PATNO, time and y. The preparation script freezes exposure scales using the anchor table first, then joins strictly subsequent valid outcomes and requires two observations. It does not choose a replacement index, infer clinical diagnoses, resolve conflicting source visit IDs, or reconstruct PPMI raw-field definitions. The archive contains the original project-specific input qualification code.

For external assays, molecular identity, normalization and calibration need separate qualification. The PPMI score calibration cannot be transferred numerically by merely renaming external lipid columns.
