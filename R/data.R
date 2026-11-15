#' Wisconsin Prognostic Breast Cancer data
#'
#' Follow-up data for 198 patients with invasive breast cancer and no evidence of distant
#' metastases at diagnosis, seen consecutively by Dr. William H. Wolberg (University of
#' Wisconsin) from 1984. The outcome is the time to recurrence. Thirty predictors describe the
#' cell nuclei on a digitised image of a fine needle aspirate of the breast mass: ten
#' characteristics are computed for each nucleus, and each characteristic is summarised over
#' the nuclei of the image by its mean, its standard error and its "worst" (largest) value.
#' The ten characteristics form ten natural groups of three predictors, which makes the data
#' a convenient example for group selection in the Cox model.
#'
#' @format A data frame with 198 rows and 34 columns:
#' \describe{
#'   \item{status}{factor: \code{"R"} recurrence, \code{"N"} no recurrence (censored).}
#'   \item{time}{recurrence time (if \code{status == "R"}) or disease-free time (if
#'     \code{status == "N"}), in months.}
#'   \item{mean_radius, SE_radius, worst_radius}{radius: mean of the distances from the centre
#'     to points on the perimeter.}
#'   \item{mean_texture, SE_texture, worst_texture}{texture: standard deviation of the
#'     gray-scale values.}
#'   \item{mean_perimeter, SE_perimeter, worst_perimeter}{perimeter.}
#'   \item{mean_area, SE_area, worst_area}{area.}
#'   \item{mean_smoothness, SE_smoothness, worst_smoothness}{smoothness: local variation in
#'     radius lengths.}
#'   \item{mean_compactness, SE_compactness, worst_compactness}{compactness.}
#'   \item{mean_concavity, SE_concavity, worst_concavity}{concavity: severity of the concave
#'     portions of the contour.}
#'   \item{mean_concavepoints, SE_concavepoints, worst_concavepoints}{concave points: number of
#'     concave portions of the contour.}
#'   \item{mean_symmetry, SE_symmetry, worst_symmetry}{symmetry.}
#'   \item{mean_fractaldim, SE_fractaldim, worst_fractaldim}{fractal dimension.}
#'   \item{tsize}{diameter of the excised tumor in centimetres.}
#'   \item{pnodes}{number of positive axillary lymph nodes observed at the time of surgery
#'     (missing for 4 patients).}
#' }
#'
#' @details The columns are ordered as in the original file: the ten means, the ten standard
#'   errors, the ten worst values, then \code{tsize} and \code{pnodes}. \code{\link{wpbc_groups}}
#'   returns the predictors arranged by characteristic together with the group vector.
#'
#' @source Wolberg, W., Street, W. and Mangasarian, O. (1995). Breast Cancer Wisconsin
#'   (Prognostic) [Dataset]. UCI Machine Learning Repository.
#'   \doi{10.24432/C5GK50}. Distributed under the Creative Commons Attribution 4.0
#'   International (CC BY 4.0) licence. The copy in this package was taken from the data set
#'   \code{wpbc} of the R package \pkg{TH.data} (version 1.1.5), which uses the variable names
#'   shown above.
#'
#' @references Street, W. N., Mangasarian, O. L. and Wolberg, W. H. (1995). An inductive
#'   learning approach to prognostic prediction. In A. Prieditis and S. Russell (eds.),
#'   \emph{Proceedings of the Twelfth International Conference on Machine Learning},
#'   pages 522-530. San Francisco: Morgan Kaufmann.
#'
#' @examples
#' data(wpbc)
#' table(wpbc$status)
#' w <- wpbc_groups()                     # complete cases, predictors arranged by group
#' fit <- rism_bma(w$time, w$status, scale(w$X), w$g_id, v = 0.25, group_names = w$group_names)
#' round(fit$pip, 2)
"wpbc"

#' Predictors and groups of the Wisconsin Prognostic Breast Cancer data
#'
#' Arranges \code{\link{wpbc}} for group selection: complete cases, the survival outcome, the
#' 32 predictors ordered by group, and the group vector. There are 12 groups: one for each of
#' the ten nuclear characteristics (its mean, standard error and worst value), and tumor size
#' and the number of positive lymph nodes as groups of one predictor each.
#'
#' @param complete if \code{TRUE} (default), the 4 patients with a missing lymph node count
#'   are removed.
#' @return A list with \code{time} (months), \code{status} (1 = recurrence, 0 = censored),
#'   \code{X} (numeric matrix, 32 columns), \code{g_id} (group of each column, 1 to 12) and
#'   \code{group_names}.
#' @examples
#' w <- wpbc_groups()
#' dim(w$X); table(w$g_id); sum(w$status)
#' @export
wpbc_groups = function(complete = TRUE) {
  e = new.env()
  utils::data("wpbc", package = "rism", envir = e)
  d = e$wpbc
  if (complete) d = stats::na.omit(d)
  feat = c("radius", "texture", "perimeter", "area", "smoothness", "compactness",
           "concavity", "concavepoints", "symmetry", "fractaldim")
  cols = c(unlist(lapply(feat, function(f) paste(c("mean", "SE", "worst"), f, sep = "_"))), "tsize", "pnodes")
  list(time = d$time, status = as.numeric(d$status == "R"), X = as.matrix(d[, cols]),
       g_id = c(rep(1:10, each = 3), 11, 12),
       group_names = c("Radius", "Texture", "Perimeter", "Area", "Smoothness", "Compactness", "Concavity",
                       "Concave points", "Symmetry", "Fractal dim.", "Tumor size", "Lymph nodes"))
}
