# PR-05: internal multivariate scalar and planar-region prototype

This PR builds on PR-04 and, transitively, PR-02. All new functions are **internal** and provisional. No NAMESPACE exports or legacy implementations are modified.

## Objects
- `.vc_new_box()`: finite n-dimensional Cartesian product of nondegenerate intervals, with independently open/closed endpoints.
- `.vc_new_scalar()`: S3 scalar field with known box domain and box study region; accepts `function(v)` or `function(x,y,...)` by explicit `input` setting.
- `.vc_new_region_type1()`: closed-boundary Type I planar region with vertical bounds.
- `.vc_new_region_type2()`: closed-boundary Type II planar region with horizontal bounds.
- `.vc_new_region_type3()`: paired Type I/II descriptions, marked `unverified`.

## Mathematical contract
- Domain and study must be dimensionally compatible and satisfy exact box inclusion.
- Evaluating inside the mathematical domain but outside study remains permitted.
- A Type III region inherits Type I for classification and membership. Its Type II description is retained as supplementary information; no global equivalence is asserted.
- Membership in a Type III region is determined by its Type I representation. The supplementary Type II representation is not required to agree during evaluation. Agreement at finitely many points would not prove global equivalence.
- Planar boundaries are checked at query points, not proven ordered for every real argument. Integration is **not** implemented yet.
- Nonfinite scalar outputs and points outside known domains are rejected.

## Tests
- Box endpoint behavior and dimension.
- Scalar fields with vector/coordinate function signatures.
- Inclusion of study in a partially open box.
- Membership in a triangle represented by Type I, II and III.
- Type III inheritance from Type I and independence of supplementary Type II membership.
- Rejection of inverted boundaries at a queried point.

## Explicit limitations
- Boxes and intervals currently require finite nondegenerate bounds.
- Type I/II regions are closed-boundary prototypes; open boundary and orientation semantics need a separate design.
- No union of regions, implicit domains, meshes, adaptive discretization or rendering yet.
- No claim of proof of equivalence or global boundary ordering.
- Tests need to pass on GitHub Actions before the PR can be considered verified.
- The internal S3 generic is a proof of concept and will need a deliberate NAMESPACE registration plan before public export.

## Merge order
PR-02 (testing) -> PR-04 (univariate) -> PR-05 (scalar regions).
PR-01 (architecture document) remains a separate design review.
PR-03 (torsion fix) is independent of these S3 prototypes.
