# CSF LPI and longitudinal nonmotor change in Parkinson disease

Analysis code accompanying **Cerebrospinal fluid lysophosphatidylinositol signals and subsequent nonmotor change in Parkinson disease**.

The study examines baseline CSF LPI(16:0)/LPI(18:0), subsequent MDS-UPDRS Part I change, matched CSF–plasma association contrasts, and the association conditional on measured phosphatidylinositol (PI). This is an exploratory observational analysis. Internal sensitivity analyses are not independent validation or evidence of causation.

## Start with the synthetic example

The example creates entirely synthetic anchors and longitudinal observations, freezes lipid standardization before follow-up exclusions, fits the core models, and exports clinical contrasts. No PPMI login or participant data are required.

```bash
conda env create -f environment.yml
conda activate pd-csf-lpi
Rscript tests/smoke.R outputs/synthetic-01
```

Each invocation requires a **new output directory**. The test checks the original-anchor scaling denominator, duplicate participant/time rejection, model qualification and 18 exposure-contrast vectors against independent analytic expressions. It does not test for a significant biological finding. Package versions used for the executed check are recorded in `validation/`.

## Analyze authorized model-ready data

```bash
Rscript scripts/run_models.R /path/to/model_ready.csv outputs/authorized-refit-01
```

This runner performs six REML fits (LPI score and PI-conditional score, each with diagonal random effects, correlated random effects and quadratic time) and four same-observation ML fits for incremental model comparisons. Outputs include fixed effects, 36 clinical contrasts, model-fit summaries, semantic comparison labels, explicit formulas and R session information.

For already-qualified molecular anchors and follow-up observations:

```bash
Rscript scripts/prepare_scores.R /path/to/anchors.csv /path/to/followup.csv outputs/prepared-01
Rscript scripts/run_models.R outputs/prepared-01/model_ready.csv outputs/refit-01
```

The input contract and clinical eligibility responsibilities are described in [DATA_SCHEMA.md](docs/DATA_SCHEMA.md). The preparation example starts **after** clinical anchor qualification; it is not a substitute for raw PPMI field reconciliation.

## Repository contents

| Path | Purpose |
| --- | --- |
| `R/models.R` | Explicit model formulas, qualified fits, metrics and corrected contrasts |
| `scripts/` | Portable score preparation and model entry points |
| `tests/smoke.R` | Synthetic end-to-end test |
| `analysis_archive/` | 28 study analysis and verification scripts, with configurable paths |
| `docs/METHODS.md` | Estimands, scaling, formulas and selection context |
| `docs/ARCHIVE.md` | Archive execution order and required upstream inputs |
| `docs/archive_provenance.json` | Source hashes, public hashes and exact adaptation categories |
| `validation/` | Executed synthetic and controlled-data numerical reproduction summaries |
| `.github/workflows/` | Synthetic-only continuous integration |

The portable runner reproduces the core model and contrast calculations. The archive also documents the original screen, paired bootstrap, clinical and measurement sensitivities, transcriptomic context, genetic analyses and external CSF tests. Those modules require the versioned inputs described in the manuscript and archive guide. The repository does not claim that all historical analyses can be recreated from synthetic data or from this download alone.

## Data access

PPMI participant-level data must be obtained through the [PPMI data-access process](https://www.ppmi-info.org/access-data-specimens/download-data). No participant-level PPMI data, saved subject-level model objects, credentials or unpublished manuscript files are distributed here. The independent postmortem CSF source is [Fernández-Irigoyen et al., 2021](https://doi.org/10.3390/biomedicines9050491). Its cross-sectional evaluation did not supply positive external association evidence and does not test future Part I change.

## Versioning and citation

Use a tagged GitHub release and its commit identifier when citing a specific code version. `CITATION.cff` identifies the code deposit; it is not a substitute for the manuscript's final author list. A DOI will be added only after an archive service has actually registered it. A GitHub URL is not a DOI.

## License

Original code is distributed under the [MIT License](LICENSE). This license does not grant rights to PPMI data, third-party datasets or dependency packages.
