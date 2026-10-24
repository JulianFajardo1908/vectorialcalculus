# AGENTS.md — vectorialcalculus engineering rules

Scope: whole repository. This file records design constraints for Codex and contributors.

## Mathematical correctness
- Treat mathematical meaning, dimensions, hypotheses, orientation, and domains as public API contracts.
- Distinguish exact mathematics from numerical approximation. Do not claim a property is proved merely because sampled points pass.
- For every numerical method specify its assumptions, stopping conditions, error estimate (if any), and limitations.
- Never change integration regions or ignore parts of a domain because a visual focus is elsewhere.
- Add analytic test oracles (where possible), edge cases, and regression tests before substantial refactors.

## S3 architecture (design in docs/architecture-s3.md)
- Four working classes: vc_univariate, vc_scalar_field, vc_trajectory, vc_vector_field; vc_univariate inherits vc_scalar_field.
- Keep domain, study, focus, view, numerics and discretization conceptually distinct.
- Implement classes using explicit constructors, validators, and documented S3 dispatch.
- Keep the mathematical engine independent from plotting and rendering; visualization cannot be an implicit side effect of a numerical operation in the new API.
- The new API must be compatible with base R |> pipes, using the mathematical object as the first argument.
- Use natural return types; do not wrap every numeric result into a tibble.
- Treat mesh policies and mesh instances as distinct. Separate visualization sampling from numerical integration/error control.
- Preserve open/closed endpoints, valid curve orientation, and dimensions.
- A planar type III region in this project means a region with valid type I and type II descriptions; verify equivalence or record that it is unverified.
- For unions with interior-disjoint pieces, avoid double-counting. Surface and singular measures require separate rules.

## Compatibility and dependencies
- Existing CRAN 1.0.5 public functions and documented returns remain stable unless a change is explicitly approved.
- Never silently rename exports, reorder parameters or change return classes/attributes.
- Prefer base R and existing dependencies. plotly/ggplot2 must remain optional for calculation; do not add dependencies casually.
- Do not edit generated NAMESPACE or man/*.Rd by hand; use roxygen2 in implementation phases.
- R minimum is currently 4.1.0, so base |> is available.
- Do not modify master directly; use small focused branches and pull requests.
- PR-00 is documentation only; no R/, man/, DESCRIPTION, NAMESPACE, tests or CI changes in PR-00.

## Performance and numerical verification
- Measure before and after: evaluation counts, elapsed time, bytes, result error and graphics object size.
- A refined mesh is not automatically better; include uniform baselines and cases where adaptive loses.
- Respect max_points and related budgets and report when tolerance cannot be reached.
- Cache evaluations only when correct for the function, arguments and context; be cautious of random or stateful functions.
- Avoid unnecessary plotly traces and large eager grid allocations.

## Documentation and collaboration
- Public R function names and source-code identifiers are in English; architecture discussion/report may be in Spanish.
- Explain every public argument, its defaults, meaning, mathematical conditions, limitations and output.
- Include basic, intermediate and advanced examples where appropriate. Use base R |> where it benefits readability.
- Keep code and filenames portable; avoid non-ASCII identifiers/filenames; accents in natural-language prose are acceptable.
- Produce a summary of changes, mathematical justification, tests executed, failures, and compatibility risks in each PR.
- Do not claim tests were run unless actually executed.
- Escalate uncertain mathematical design or any breaking change for the package author to decide.
