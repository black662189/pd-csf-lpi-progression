# Statistical definitions

## Study context

Part I was selected during project exploration before the retained endpoint-specific 547-feature/biofluid screen. Subsequent fixed-score, parent-conditioned and diagnostic analyses are candidate-selected internal follow-ups. Their family-specific corrections do not remove the earlier endpoint selection. The portable runner therefore does not assign a discovery q value to its focused tests.

## Frozen exposure

For each annotated species, positive released relative signal is log2-transformed and standardized using the mean and sample SD among clinically eligible original index samples, before exclusion for future outcome completeness. The LPI score is the equal-weight mean of LPI(16:0) and LPI(18:0) species z scores, standardized again among complete eligible anchors. Both species are required. The original CSF scaling population contained 227 anchors; the longitudinal score model included 200 participants and 1,389 observations. Do not restandardize the 200-person subset.

The PI score uses the six fixed annotated PI members in the manuscript, following the same species-z/mean/restandardization scheme, with all members required. It is a relative-signal score, not a molar ratio or enzyme-activity measurement. Its conditional model includes 199 participants and 1,379 observations. Matched CSF/plasma bootstrap analyses use their separately documented within-pair scales; these are not interchangeable with the original index scales.

## Core formula

```r
y ~ time * (x + y0 + up0 + age0 + sex + duration0 + ledd0 + genetic_group) +
    (1 + time || PATNO)
```

`time` is years since the molecular index; only strictly subsequent observations enter. `y` is MDS-UPDRS Part I total. `x` is the frozen baseline LPI score. `y0` and `up0` are baseline Part I and OFF/untreated Part III. `age0=age/10`, `duration0=duration/5`, and `ledd0=log1p(LEDD)`. The PI-conditional formula adds `parent` inside the time-interacted covariate expression.

Models use REML, the bobyqa optimizer, a 100,000 evaluation limit, and lmerTest Satterthwaite inference. Rank deficiency, optimizer convergence messages, singular random effects or nonfinite inference stop the portable run. No fallback formula is selected automatically.

The correlated specification replaces `||` with `|`. The quadratic specification replaces `time` outside the parentheses with `(time + time_sq)`, where `time_sq=(time/5)^2`, interacting both time terms with **every** covariate. The random slope remains in linear time. ML comparisons of different fixed effects use identical observations. Incremental R-squared is in-sample fit, not validated prediction.

## Contrasts

For a horizon `h`, the difference in change per one-SD exposure contrast is evaluated by the fixed-effects design:

```text
[X(h,+0.5)-X(h,-0.5)] - [X(0,+0.5)-X(0,-0.5)]
```

The implementation uses a distinct `exposure_value` argument and explicit assignment to avoid data.table column/argument shadowing. An independent check requires the contrast vector to contain `h` at `time:x`, `(h/5)^2` at `time_sq:x` when present, and zeros elsewhere. CIs use the complete fixed-effect covariance and Satterthwaite degrees of freedom. Cohort-standardized expected change gives each included participant one equal-weight baseline covariate row. These quantities are between-person model contrasts, not treatment effects or MCID estimates.

## Bootstrap and other modules

The archived scripts retain whole-participant resampling, distinct cluster labels for repeated selections, fixed exposure scales, seeds and original family definitions. The paired bootstrap resamples both fluids jointly. The synthetic smoke test is a software test; it is neither an external cohort nor evidence that the biological model is correct.
