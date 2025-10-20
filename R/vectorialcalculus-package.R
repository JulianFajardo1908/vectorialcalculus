#' vectorialcalculus: 3D Vector Calculus Tools
#'
#' Functions for visualization and analysis in vector calculus:
#' parametric curves, scalar fields, and vector fields. The package
#' provides utilities to compute and illustrate objects such as
#' gradients, divergences, curls, line integrals, and surface integrals,
#' alongside 2D/3D pedagogical graphics.
#'
#' @details
#' \itemize{
#'   \item \strong{Parametric curves}: tangent/normal/binormal frames,
#'         curvature and torsion, arc-length and work along a path.
#'   \item \strong{Scalar fields}: gradient and directional derivatives,
#'         level sets and isosurfaces.
#'   \item \strong{Vector fields}: divergence, curl, and flux/ circulation
#'         via surface and line integrals.
#'   \item \strong{Visualization}: 2D/3D plots to support teaching/learning
#'         (arrows, streamlines, surfaces, and overlaid trajectories).
#' }
#'
#' Typical workflow:
#' \enumerate{
#'   \item Define the curve/field.
#'   \item Compute geometric/analytic quantities (e.g., gradient, curl).
#'   \item Visualize and interpret the results in 2D/3D.
#' }
#'
#' @section Authors:
#' \itemize{
#'   \item Julio Lizarazo Osorio (\emph{Author}) - Universidad Pedagógica y Tecnológica de Colombia (UPTC).
#'   \item Julian Mauricio Fajardo (\emph{Author, Maintainer}) - Universidad El Bosque.
#' }
#'
#' @author
#' Julio Lizarazo Osorio \email{julio.lizarazo@uptc.edu.co}
#'
#' Julian Mauricio Fajardo \email{julian.fajardo@unbosque.edu.co}
#'
#' @seealso
#' Project page and issue tracker (if applicable):
#' \itemize{
#'   \item \url{https://github.com/tuusuario/vectorialcalculus}
#'   \item \url{https://github.com/tuusuario/vectorialcalculus/issues}
#' }
#'
#' @references
#' Standard vector calculus texts and teaching references.
#'
#' @keywords package
#' @docType package
#' @name vectorialcalculus
#' @aliases vectorialcalculus-package
#'
#' @examples
#' \dontrun{
#' # Example (sketch) - gradient of a simple scalar field:
#' f <- function(x, y, z) x^2 + y^2 + z
#' # gradient(f, c("x","y","z"))  # assuming exported by the package
#'
#' # Visual aid (sketch) - plot a parametric curve and its tangent:
#' # r <- function(t) cbind(sin(t), cos(t), t/5)
#' # plot_curve3d(r, t = seq(0, 6*pi, length.out = 200))
#' }
"_PACKAGE"
