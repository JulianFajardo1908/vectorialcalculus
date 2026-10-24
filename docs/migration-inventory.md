# Legacy migration inventory (PR-00)

Source: repository `master` at commit `3a826de810bc555636afd3c1796ff711d249cab9`, 2026-10-24. This is a **filename-and-export inventory**, not a complete behavior audit. Check each implementation, formal parameters, return type, side effects and numerical assumptions before changing code. Existing files should not be rewritten in PR-00.

Current repository has 47 source files in `R/`; the 46 names exported in `NAMESPACE` differ because `R/partial_derivative_surface.R` exports `partial_derivatives_surface`. Source files are listed below exactly as found.

| R source | Candidate target module |
| --- | --- |
| `R/arc_length3d.R` | Trajectory / differential geometry |
| `R/binormal3d.R` | Trajectory / differential geometry |
| `R/critical_points_2d.R` | Scalar fields / derivatives / optimization |
| `R/critical_points_nd.R` | Scalar fields / derivatives / optimization |
| `R/curl3d.R` | Vector fields / line integrals |
| `R/curvature_torsion3d.R` | Trajectory / differential geometry |
| `R/curve_sample3d.R` | Trajectory / differential geometry |
| `R/cylindrical_surface3d.R` | Regions / surfaces / visualization |
| `R/directional_derivative3d.R` | Scalar fields / derivatives / optimization |
| `R/divergence_field.R` | Vector fields / line integrals |
| `R/frenet_frame3d.R` | Trajectory / differential geometry |
| `R/gradient_direction2d.R` | Scalar fields / derivatives / optimization |
| `R/gradient_scalar.R` | Scalar fields / derivatives / optimization |
| `R/integrate_double_polar.R` | Integration / pedagogical quadrature |
| `R/integrate_double_xy.R` | Integration / pedagogical quadrature |
| `R/integrate_triple_general.R` | Integration / pedagogical quadrature |
| `R/lagrange_check.R` | Scalar fields / derivatives / optimization |
| `R/line_integral2d.R` | Vector fields / line integrals |
| `R/line_integral3d_work.R` | Vector fields / line integrals |
| `R/line_integral_vector2d.R` | Vector fields / line integrals |
| `R/newton_raphson2d.R` | Scalar fields / derivatives / optimization |
| `R/newton_raphson_anim.R` | Scalar fields / derivatives / optimization |
| `R/normal3d.R` | Trajectory / differential geometry |
| `R/osculating_circle3d.R` | Trajectory / differential geometry |
| `R/osculating_ribbon3d.R` | Trajectory / differential geometry |
| `R/partial_derivatives_surface.R` | Scalar fields / derivatives / optimization |
| `R/plot_curve3d.R` | Trajectory / differential geometry |
| `R/plot_surface_with_tangents.R` | Regions / surfaces / visualization |
| `R/region_xyz0.R` | Regions / surfaces / visualization |
| `R/related_rates_grad.R` | Scalar fields / derivatives / optimization |
| `R/riemann_prisms3d.R` | Integration / pedagogical quadrature |
| `R/riemann_rectangles2d.R` | Integration / pedagogical quadrature |
| `R/riemann_sum_1d_plot.R` | Integration / pedagogical quadrature |
| `R/riemann_sum_2d_plot.R` | Integration / pedagogical quadrature |
| `R/secant_tangent.R` | Scalar fields / derivatives / optimization |
| `R/solid_cylindrical3d.R` | Regions / surfaces / visualization |
| `R/solid_of_revolution_y.R` | Regions / surfaces / visualization |
| `R/solid_spherical3d.R` | Regions / surfaces / visualization |
| `R/solid_xyz3d.R` | Regions / surfaces / visualization |
| `R/streamline_and_field3d.R` | Vector fields / line integrals |
| `R/surface_integral_z.R` | Regions / surfaces / visualization |
| `R/surface_parametric_area.R` | Regions / surfaces / visualization |
| `R/tangent3d.R` | Trajectory / differential geometry |
| `R/tangent_plane3d.R` | Scalar fields / derivatives / optimization |
| `R/total_differential_nd.R` | Scalar fields / derivatives / optimization |
| `R/vector_field3d.R` | Vector fields / line integrals |
| `R/xy_region.R` | Regions / surfaces / visualization |

## Migration review checklist per function

- Capture original formals/default values and accepted input kinds.
- Capture exact scalar/vector/list/tibble/plotly return schema and print/plot side effects.
- Record numeric integration and derivative methods, accuracy constraints and optional dependencies.
- Add analytic or metamorphic regression tests, including orientation and boundary cases.
- Decide whether to reuse a `vc_*` evaluation/region/sampling engine without altering old behavior.
- Document any incompatible bug correction separately and request author approval.

## Initial crosswalk examples

| Legacy function | Planned new operation |
| --- | --- |
| `curve_sample3d()` | `vc_sample()` for trajectories |
| `arc_length3d()` | `vc_arc_length()` |
| `gradient_scalar()` | `vc_gradient()` |
| `curl3d()` | `vc_curl()` in three dimensions |
| `divergence_field()` | `vc_divergence()` |
| `integrate_double_xy()` | `vc_integrate()` with planar regions |
| `integrate_triple_general()` | Three-dimensional region integration |
| `xy_region()` | Type I region display |
| `line_integral2d()` | Scalar line integral with trajectory |
| `line_integral3d_work()` | Vector line integral/work |
| `vector_field3d()` | `vc_sample()` and `vc_plot()` for vector fields |
| `solid_xyz3d()` | Region geometry, volume and rendering |

**Open task for PR-01:** automated export-vs-source inventory, signatures/returns/side effects, contract tests and reproducible baseline. No claim of API-level equivalence before these tests.
