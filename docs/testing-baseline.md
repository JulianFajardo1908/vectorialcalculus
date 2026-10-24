# PR-01: baseline tests and known limitations

**Scope:** initial testthat suite and CI only. No changes to numerical implementations.
**Base:** vectorialcalculus 1.0.5, master at 3a826de810bc555636afd3c1796ff711d249cab9.
**Execution:** Tests have not been executed locally because R is unavailable in the current environment. Consult GitHub Actions results before claiming the suite passes.

## Test infrastructure
- tests/testthat.R uses the standard testthat R CMD check entry point.
- testthat (>=3.0.0) and edition 3 are already specified in DESCRIPTION.
- The new workflow runs release R package checks on Ubuntu and Windows for pull requests and pushes to master.
- No DESCRIPTION, NAMESPACE, R implementation, or man files are modified.
- Initial tests use plot=FALSE whenever the current API permits.

## Analytic oracles

| Function | Test | Analytic value |
| --- | --- | --- |
| gradient_scalar | x^2 + 3y^2 at (2,-1) | (4,-6) |
| divergence_field | (x,2y,3z) | 6 |
| curl3d | (-y,x,0) | (0,0,2) |
| curl3d with time | (-y,x+t,z) | (0,0,2) |
| curve_sample3d | (t,t^2,0), five nodes | Exact tibble data |
| arc_length3d | Unit circle on [0,2*pi] | 2*pi |
| integrate_double_xy | x+y over triangle x,y>=0, x+y<=1 | 1/3, both orders |

## Confirmed source-level mathematical issue

The third-derivative finite-difference formula in R/curvature_torsion3d.R uses:

```r
(f(t - 2 * h) - 2 * f(t - h) + 2 * f(t + h) - f(t + 2 * h)) /
  (2 * h^3)
```

For f(t)=t^3 this returns -6 instead of +6, reversing the sign of the third derivative. It can reverse torsion. Analytically the helix r(t)=(cos(t), sin(t), t) has curvature 1/2 and torsion +1/2.

A test records the desired behavior in test-known-issues.R but is deliberately skipped until a separate mathematical correction PR changes the implementation. Do not mistake the skipped test for passing coverage.

## Further issues for later investigation
- arc_length3d defaults to pracma::romberg; baseline test explicitly selects stats::integrate, so that suggested dependency is not required for this test.
- integrate_double_xy requests purrr for numerical computation and checks ggplot2 even when plot_domain=FALSE. Tests skip when these suggested dependencies are absent.
- integrate_double_xy returns zero for reversed inner limits, potentially hiding an invalid region; do not correct until mathematical semantics have been approved.
- This suite samples only a small fraction of current public exports. It is not evidence of complete coverage.

## Instructions for reproducible execution

Install R >= 4.1, the package dependencies and testthat. Then, from the package directory, run:

```r
testthat::test_local(".")
```

Run R CMD check (e.g. through the GitHub Actions workflow) and inspect logs. Record failure/skip counts and platform-specific results before approving further refactoring.

## Next implementation tasks
1. Triage CI failures; distinguish existing failures from mistakes in new tests.
2. Expand tests for each public function and edge cases.
3. Fix the torsion sign in an independently reviewed mathematical correction PR.
4. Begin new S3 classes only after the baseline is trustworthy.
