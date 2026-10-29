# Method handoff, 2026-10-29. Not executed on project data by the method reviewer.
# Input rows must already pass the main analysis' explicit identity and finite-data QC.

estimate_sdy_variable_n <- function(se, maf, n) {
    stopifnot(length(se) == length(maf), length(se) == length(n))
    stopifnot(all(is.finite(se) & se > 0),
              all(is.finite(maf) & maf > 0 & maf <= 0.5),
              all(is.finite(n) & n > 0))
    v <- se^2
    inv_v <- 1 / v
    n_var_g <- 2 * n * maf * (1 - maf)
    # Algebraically identical to coloc::sdY.est's no-intercept regression,
    # with the actual per-variant n supplied outside the coloc dataset object.
    variance_estimate <- sum(inv_v * n_var_g) / sum(inv_v^2)
    stopifnot(is.finite(variance_estimate), variance_estimate > 0)
    per_variant_scale <- sqrt(n_var_g * v)
    list(sdY = sqrt(variance_estimate),
         n_range = range(n),
         scale_diagnostic_quantiles = quantile(
             per_variant_scale, c(0.05, 0.25, 0.5, 0.75, 0.95)),
         note = paste("Approximate effective scale from summary statistics;",
                      "not a measured population SD or confidence interval."))
}

make_quant_dataset <- function(d, sd_y) {
    stopifnot(is.finite(sd_y), sd_y > 0, !anyDuplicated(d$snp))
    list(beta = d$beta, varbeta = d$se^2,
         snp = d$snp, position = d$position,
         type = "quant", sdY = sd_y)
}

make_pd_dataset <- function(d) {
    stopifnot(!anyDuplicated(d$snp))
    # beta is log-odds; direct beta/SE ABF does not use N or case fraction.
    list(beta = d$beta, varbeta = d$se^2,
         snp = d$snp, position = d$position, type = "cc")
}

prior_grid <- data.frame(
    p1 = rep(1e-4, 3), p2 = rep(1e-4, 3),
    p12 = c(1e-6, 1e-5, 5e-5))

# Main expression sdY=1 is an INT-scale approximation, not a universal default.
# Lipid sdY must be estimated/known separately for each lipid before any pairwise merge.
# Recommended lipid scale sensitivity: multiply its fixed estimate by c(0.5, 1, 2).
# Call coloc::check_dataset() then coloc::coloc.abf() with beta/varbeta inputs.
# Log sessionInfo(), packageVersion('coloc'), warnings and every H0-H4 posterior.
# Do not interpret pairwise posteriors as a joint three-trait posterior.
