# Gamma Model Simulation: MLE vs. Method of Moments

> **Key result:** maximum likelihood is more efficient than the method of moments (MM) at every sample size tested. A Monte Carlo study (10,000 replications per sample size) also exposed an error in the MM covariance matrix given in the course assignment. With the corrected formula, MM coverage at n = 500 rises from 0.896 to 0.951 for the shape and from 0.860 to 0.947 for the scale.

![Coverage by sample size](results/coverage.png)

## Problem

How well do asymptotic results for the MLE and MM hold in finite samples of a Gamma distribution? The study compares bias, MSE, interval width, and 95% coverage.

## Data

The data are simulated: Gamma(α₀ = 2, β₀ = 3), with n ∈ {20, 50, 100, 200, 500} and seed 42.

## Method

- **MM:** closed-form estimates, with delta-method standard errors.
- **MLE:** Newton-Raphson on the log scale for the shape, started at the MM estimate, with standard errors from expected and observed information.
- **Intervals:** 95% Wald intervals. Tail errors are recorded separately.

## The covariance correction

The assignment PDF, kept in [`assignment/`](assignment/) in its original Portuguese, states the MM covariance with Var(α̂) = α(2+α) and Var(β̂) = β²(2+α)/α. The delta method gives a different result:

| Term | Assignment | Delta method (used here) |
|---|---|---|
| Var(α̂) | α(2+α) | 2α(α+1) |
| Var(β̂) | β²(2+α)/α | β²(2α+3)/α |
| Cov(α̂, β̂) | β(2+α) | −2β(α+1) |

The derivation is in the [report](report/gamma_simulation_report.pdf). The correction changes only MM standard errors, widths, and coverage. MM estimates, bias, and MSE stay the same because the seed is unchanged.

| n = 500 | Empirical SD | Mean SE, assignment | Mean SE, corrected | Coverage, assignment | Coverage, corrected |
|---|---:|---:|---:|---:|---:|
| α̂ (MM) | 0.155 | 0.127 | 0.156 | 0.896 | 0.951 |
| β̂ (MM) | 0.248 | 0.189 | 0.251 | 0.860 | 0.947 |

## Results

- **Efficiency:** MLE has lower MSE than MM at every n, for both parameters.
- **Small samples:** coverage for β is low at n = 20 (0.852 for MM, 0.867 for MLE), with most misses on one side.
- **Large samples:** at n = 500, all methods are close to the nominal 95%.
- **Fisher vs. observed information:** the two coincide, which follows algebraically from the Gamma likelihood at the MLE.

Full results are in [`results/full_results.csv`](results/full_results.csv).

## Limitations

- Only one parameter pair (α₀ = 2, β₀ = 3) was studied.
- Monte Carlo error applies: with 10,000 replications, coverage near 0.95 has a standard error of about 0.002.

## How to reproduce

Requires R with `ggplot2`. Tested with R 4.5.3.

```bash
Rscript src/gamma_simulation.R
```

Run it from the repository root. It rewrites `results/`.

## Repository structure

```text
assignment/   Original course assignment (Portuguese)
report/       Report in English (LaTeX source and PDF)
results/      Table and figures
src/          R simulation
```
