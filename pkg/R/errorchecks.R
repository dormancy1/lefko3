#' Summary of Class "lefkoCondMat"
#' 
#' This function provides basic information summarizing the characteristics of
#' conditional matrices derived from a \code{lefkoCondMat} object.
#' 
#' @name summary.lefkoCondMat
#' 
#' @param object An object of class \code{lefkoCondMat}.
#' @param ... Other parameters.
#' 
#' @return A text summary of the object shown on the console, showing the number
#' of historical matrices, as well as the number of conditional matrices nested
#' within each historical matrix.
#' 
#' @examples
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' lathcondmats <- cond_hmpm(ehrlen3)
#' summary(lathcondmats)
#' 
#' # Cypripedium  example
#' data(cypdata)
#'  
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4, 
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04", 
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
#'     "D", "XSm", "Sm", "D", "XSm", "Sm", "mat", "mat", "mat", "SD", "P1"),
#'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
#'     "SL", "SL", "D", "XSm", "Sm", "rep", "rep"),
#'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
#'     "SL", "SL", "SL", "SL", "SL", "SL", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
#'     "mat", "mat", "mat", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
#'     "XSm", "D", "XSm", "Sm", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
#'     "XSm", "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
#'     NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'     NA, 0.5, 0.5),
#'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#'   stageframe = cypframe_raw, historical = TRUE)
#' 
#' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added", "size1added"), 
#'   supplement = cypsupp3r, yearcol = "year2", patchcol = "patchid",
#'   indivcol = "individ")
#' 
#' cypcondmats <- cond_hmpm(cypmatrix3r)
#' 
#' summary(cypcondmats)
#' 
#' @export
summary.lefkoCondMat <- function(object, ...) {
  
  histmatrices <- object$Mcond
  condmatrices <- histmatrices[[1]]
  firstcondmat <- condmatrices[[1]]
  
  numhistmats <- length(histmatrices)
  prevstages <- length(condmatrices)
  matdim <- dim(firstcondmat)
  
  writeLines(paste0("\nThis lefkoCondMat object contains ", prevstages,
      " conditional matrices per historical matrix."))
  writeLines(paste0("It covers ", numhistmats, " main historical matrices."))
  writeLines(paste0("Each conditional matrix is a square matrix with ", matdim[1],
      " rows and columns, and a total of ", matdim[1]*matdim[1], " elements."))
  writeLines(paste0("\nThe order of conditional matrices corresponding to stage in occasion t-1 is:\n",
      paste(object$ahstages$stage, collapse = " ")))
  writeLines("\nThe order of historical matrices is: \n")
  print.data.frame(object$labels)
  
  writeLines("\nThe order of conditional matrices matches the stage column in object $ahstages.")
  writeLines("The order of historical matrices follows that shown in object $labels.")
}

#' Create Matrix Image
#' 
#' Function \code{image3()} is a generic function that creates matrix plots.
#' 
#' @name image3
#' 
#' @param mats A lefkoMat object, or a single projection matrix, for which the
#' dominant eigenvalue is desired.
#' @param ... Other parameters
#' 
#' @return Produces a single matrix image, or a series of images, depending on
#' the input. Non-zero elements appear as red space, while zero elements appear
#' as white space.
#' 
#' @seealso \code{\link{image3.lefkoMat}()}
#' @seealso \code{\link{image3.matrix}()}
#' 
#' @examples
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' image3(ehrlen3, used = 1, type = "U")
#' 
#' @export
image3 <- function(mats, ...) UseMethod("image3")

#' Create Matrix Image(s) for lefkoMat Object
#' 
#' Function \code{image3.lefkoMat} plots matrix images for matrices supplied
#' within \code{lefkoMat} objects.
#' 
#' @name image3.lefkoMat
#' 
#' @param mats A \code{lefkoMat} object.
#' @param used A numeric value or vector designating the matrices to plot. Can
#' also take the value \code{"all"}, which plots all matrices. Defaults to
#' \code{"all"}.
#' @param type Character value indicating whether to plot \code{A}, \code{U}, or
#' \code{F} matrices. Defaults to \code{"A"}.
#' @param ... Other parameters.
#' 
#' @return Plots a matrix image, or series of matrix images, denoting non-zero
#' elements as red space and zero elements as white space.
#' 
#' @examples 
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' image3(ehrlen3, used = 1, type = "U")
#' 
#' @export
image3.lefkoMat <- function(mats, used = "all", type = "A", ...) {
  
  allmats <- c(1:length(mats$A))
  
  if (!is.character(type)) {
    stop("Please enter A, F, or U for type option.", call. = FALSE)
  }
  
  type <- tolower(type)
  if (!is.element(type, c("a", "u", "f"))) {
    stop("Please enter A, F, or U for type option.", call. = FALSE)
  }
  
  if (all(is.character(used))) {
    if (all(tolower(used) != "all")) {
      stop("Value entered for matrix option not recognized.", call. = FALSE)
    } else {
      chosen_mat <- allmats
    }
  } else if (is.numeric(used) & is.element(used, allmats)) {
    chosen_mat <-  used
  } else {
    stop("Value entered for matrix option not recognized.", call. = FALSE)
  }
  
  if (type == "u") {
    chosen_list <- mats$U[chosen_mat]
  } else if (type == "f") {
    chosen_list <- mats$F[chosen_mat]
  } else {
    chosen_list <- mats$A[chosen_mat]
  }
  
  if (is(chosen_list[[1]], "dgCMatrix")) {
    lapply(chosen_list, function(X) {Matrix::image(X, xlab = "", ylab = "", 
      sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  } else {
    lapply(chosen_list, function(X) {Matrix::image(Matrix::Matrix(X, sparse = TRUE),
      xlab = "", ylab = "", sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  }
}

#' Create a Matrix Image for a Single Matrix
#' 
#' Function \code{image3.matrix} plots a matrix image for a single matrix.
#' 
#' @name image3.matrix
#' 
#' @param mats A \code{matrix} class object.
#' @param ... Other parameters.
#' 
#' @return Plots a matrix image, or series of matrix images, denoting non-zero
#' elements as red space and zero elements as white space.
#' 
#' @examples 
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ", sparse_output = FALSE)
#' 
#' image3(ehrlen3$U[[1]])
#' 
#' @export
image3.matrix <- function(mats, ...) {
  
  Matrix::image(Matrix::Matrix(mats, sparse = TRUE), xlab = "", ylab = "",
    sub = "", col.regions = c("white", "red"), lwd = 0,
    at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)
}

#' Create a Matrix Image for a Single Sparse Matrix
#' 
#' Function \code{image3.dgCMatrix} plots a matrix image for a single sparse
#' matrix.
#' 
#' @name image3.dgCMatrix
#' 
#' @param mats A \code{matrix} class object.
#' @param ... Other parameters.
#' 
#' @return Plots a matrix image, or series of matrix images, denoting non-zero
#' elements as red space and zero elements as white space.
#' 
#' @examples 
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ", sparse_output = TRUE)
#' 
#' image3(ehrlen3$U[[1]])
#' 
#' @export
image3.dgCMatrix <- function(mats, ...) {
  
  Matrix::image(mats, xlab = "", ylab = "", sub = "",
    col.regions = c("white", "red"), lwd = 0, at = c(0, 0.0000000000001, Inf),
    drop.unused.levels = FALSE)
}

#' Create Matrix Images for Matrices in a List
#' 
#' Function \code{image3.list} plots matrix images for matrices contained in a
#' list of matrices.
#' 
#' @name image3.list
#' 
#' @param mats A \code{list} class object.
#' @param used A numeric vector of projection matrices within \code{mats} to
#' represent as matrix images. Can also take the text value \code{"all"}, which
#' will produce images of all matrices. Defaults to \code{"all"}.
#' @param ... Other parameters.
#' 
#' @return Plots a matrix image, or series of matrix images, denoting non-zero
#' elements as red space and zero elements as white space.
#' 
#' @examples 
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' image3(ehrlen3$A, used = 1)
#' 
#' @export
image3.list <- function(mats, used = "all", ...) {
  
  allmats <- c(1:length(mats))
  
  if (all(is.character(used))) {
    if (all(tolower(used) != "all")) {
      stop("Value entered for matrix option not recognized.", call. = FALSE)
    } else {
      chosen_mat <- allmats
    }
  } else if (is.numeric(used) & is.element(used, allmats)) {
    chosen_mat <-  used
  } else {
    stop("Value entered for matrix option not recognized.", call. = FALSE)
  }
  
  chosen_list <- mats[chosen_mat]
  
  if (is(chosen_list[[1]], "dgCMatrix")) {
    lapply(chosen_list, function(X) {Matrix::image(X, xlab = "", ylab = "", 
      sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  } else if (is.matrix(chosen_list[[1]])) {
    lapply(chosen_list, function(X) {Matrix::image(Matrix::Matrix(X, sparse = TRUE),
      xlab = "", ylab = "", sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  } else {
    stop("Chosen elements include non-matrix objects. Please choose only list
      elements containing matrix objects.", call. = FALSE)
  }
}

#' Create Matrix Image(s) for lefkoSens Object
#' 
#' Function \code{image3.lefkoSens} plots matrix images for sensitivity matrices
#' supplied within \code{lefkoSens} objects.
#' 
#' @name image3.lefkoSens
#' 
#' @param mats A \code{lefkoSens} object.
#' @param used A numeric value or vector designating the matrices to plot. Can
#' also take the value \code{"all"}, which plots all matrices. Defaults to
#' \code{"all"}.
#' @param type Character value indicating whether to plot \code{"a"}historical or
#' \code{"h"}istorical sensitivity matrices. Defaults to \code{"a"}historical,
#' but will plot a historical sensitivity matrix image if no ahistorical
#' sensitivity matrix exists.
#' @param ... Other parameters.
#' 
#' @return Plots a matrix image, or series of matrix images, denoting non-zero
#' elements as red space and zero elements as white space.
#' 
#' @examples 
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
#'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
#'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen_sens <- sensitivity3(ehrlen3)
#' 
#' image3(ehrlen_sens, used = 1, type = "h")
#' 
#' @export
image3.lefkoSens <- function(mats, used = "all", type = "a", ...) {
  
  allahmats <- c(1:length(mats$ah_sensmats))
  allhmats <- c(1:length(mats$h_sensmats))
  
  allmats <- c(1:max(c(allahmats, allhmats)))
  
  if (!is.character(type)) {
    stop("Please enter a or h for type option.", call. = FALSE)
  }
  
  type <- tolower(type)
  if (!is.element(type, c("a", "h"))) {
    stop("Please enter a or h for type option.", call. = FALSE)
  }
  
  if (all(is.character(used))) {
    if (all(tolower(used) != "all")) {
      stop("Value entered for matrix option not recognized.", call. = FALSE)
    } else {
      chosen_mat <- allmats
    }
  } else if (is.numeric(used) & is.element(used, allmats)) {
    chosen_mat <-  used
  } else {
    stop("Value entered for matrix option not recognized.", call. = FALSE)
  }
  
  if (type == "h") {
    if (any(is.null(mats$h_sensmats))) {
      stop("This object does not appear to have historical sensitivity matrices. Please try ahistorical option.",
        call. = FALSE)
    }
    chosen_list <- mats$h_sensmats[chosen_mat]
  } else {
    if (any(is.null(mats$ah_sensmats))) {
      warning("This object does not appear to have ahistorical sensitivity matrices. Will use historical sensitivity matrices instead.",
        call. = FALSE)
      
      chosen_list <- mats$h_sensmats[chosen_mat]
    } else {
      chosen_list <- mats$ah_sensmats[chosen_mat]
    }
  }
  
  if (is(chosen_list[[1]], "dgCMatrix")) {
    lapply(chosen_list, function(X) {Matrix::image(X, xlab = "", ylab = "", 
      sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  } else {
    lapply(chosen_list, function(X) {Matrix::image(Matrix::Matrix(X, sparse = TRUE),
      xlab = "", ylab = "", sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  }
}

#' Create Matrix Image(s) for lefkoElas Object
#' 
#' Function \code{image3.lefkoElas} plots matrix images for elasticity matrices
#' supplied within \code{lefkoElas} objects.
#' 
#' @name image3.lefkoElas
#' 
#' @param mats A \code{lefkoElas} object.
#' @param used A numeric value or vector designating the matrices to plot. Can
#' also take the value \code{"all"}, which plots all matrices. Defaults to
#' \code{"all"}.
#' @param type Character value indicating whether to plot \code{"a"}historical or
#' \code{"h"}istorical elasticity matrices. Defaults to \code{"a"}historical,
#' but will plot a historical elasticity matrix image if no ahistorical
#' elasticity matrix exists.
#' @param ... Other parameters.
#' 
#' @return Plots a matrix image, or series of matrix images, denoting non-zero
#' elements as red space and zero elements as white space.
#' 
#' @examples 
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
#'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
#'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen_elas <- elasticity3(ehrlen3)
#' 
#' image3(ehrlen_elas, used = 1, type = "h")
#' 
#' @export
image3.lefkoElas <- function(mats, used = "all", type = "a", ...) {
  
  allahmats <- c(1:length(mats$ah_elasmats))
  allhmats <- c(1:length(mats$h_elasmats))
  
  allmats <- c(1:max(c(allahmats, allhmats)))
  
  if (!is.character(type)) {
    stop("Please enter a or h for type option.", call. = FALSE)
  }
  
  type <- tolower(type)
  if (!is.element(type, c("a", "h"))) {
    stop("Please enter a or h for type option.", call. = FALSE)
  }
  
  if (all(is.character(used))) {
    if (all(tolower(used) != "all")) {
      stop("Value entered for matrix option not recognized.", call. = FALSE)
    } else {
      chosen_mat <- allmats
    }
  } else if (is.numeric(used) & is.element(used, allmats)) {
    chosen_mat <-  used
  } else {
    stop("Value entered for matrix option not recognized.", call. = FALSE)
  }
  
  if (type == "h") {
    if (any(is.null(mats$h_elasmats))) {
      stop("This object does not appear to have historical sensitivity matrices. Please try ahistorical option.",
        call. = FALSE)
    }
    chosen_list <- mats$h_elasmats[chosen_mat]
  } else {
    if (any(is.null(mats$ah_elasmats))) {
      warning("This object does not appear to have ahistorical sensitivity matrices. Will use historical sensitivity matrices instead.",
        call. = FALSE)
      
      chosen_list <- mats$h_elasmats[chosen_mat]
    } else {
      chosen_list <- mats$ah_elasmats[chosen_mat]
    }
  }
  
  if (is(chosen_list[[1]], "dgCMatrix")) {
    lapply(chosen_list, function(X) {Matrix::image(X, xlab = "", ylab = "", 
      sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  } else {
    lapply(chosen_list, function(X) {Matrix::image(Matrix::Matrix(X, sparse = TRUE),
      xlab = "", ylab = "", sub = "", col.regions = c("white", "red"), lwd = 0,
      at = c(0, 0.0000000000001, Inf), drop.unused.levels = FALSE)})
  }
}

#' Calculate Difference Matrices Between lefkoMat Objects of Equal Dimensions
#' 
#' Function \code{diff_lM()} takes two \code{lefkoMat} objects with completely
#' equal dimensions, including both the size and number of matrices, and
#' gives the matrix differences between each corresponding set.
#' 
#' @name diff_lM
#' 
#' @param mpm1 The first \code{lefkoMat} object.
#' @param mpm2 The second \code{lefkoMat} object.
#' 
#' @return An object of class \code{lefkoDiff}, which is a set of \code{A},
#' \code{U}, and \code{F} matrices corresponding to the differences between each
#' set of matrices, followed by the \code{hstages}, \code{ahstages}, and
#' \code{labels} elements from each input \code{lefkoMat} object. Elements
#' labelled with a \code{1} at the end refer to \code{mpm1}, while those
#' labelled \code{2} at the end refer to \code{mpm2}.
#' 
#' @section Notes:
#' The exact difference is calculated as the respective matrix in \code{mpm1}
#' minus the corresponding matrix in \code{mpm2}.
#' 
#' This function first checks to see if the number of matrices is the same, and
#' then whether the matrix dimensions are the same. If the two sets differ in at
#' least one of these characteristics, then the function will yield a fatal
#' error.
#' 
#' If the lengths and dimensions of the input \code{lefkoMat} objects are the
#' same, then this will check if the \code{labels} element is essentially the
#' same. If not, then the function will yield a warning, but will still operate.
#' 
#' @examples
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 3, 6, 11, 19.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1.5, 1.5, 3.5, 5)
#' comments <- c("Dormant seed", "1st yr protocorm", "2nd yr protocorm",
#'   "3rd yr protocorm", "Seedling", "Dormant adult",
#'   "Extra small adult (1 shoot)", "Small adult (2-4 shoots)",
#'   "Medium adult (5-7 shoots)", "Large adult (8-14 shoots)",
#'   "Extra large adult (>14 shoots)")
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset, 
#'   binhalfwidth = binvec, comments = comments)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' seeds_per_pod <- 5000
#' 
#' cypsupp2_raw <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
#'   givenrate = c(0.03, 0.15, 0.1, 0.1, 0.1, 0.05, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, (0.5 * seeds_per_pod),
#'     (0.5 * seeds_per_pod)),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' cypsupp3_raw <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3",
#'     "SL", "SL", "SL", "D", "D", "SD", "P1"),
#'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
#'     "rep", "rep"),
#'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "SL", "P3",
#'     "SL", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "D", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
#'   givenrate = c(0.01, 0.05, 0.10, 0.20, 0.1, 0.1, 0.05, 0.05, 0.05, NA, NA,
#'     NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'     (0.5 * seeds_per_pod), (0.5 * seeds_per_pod)),
#'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#'   stageframe = cypframe_raw, historical = TRUE)
#' 
#' cypmatrix2rp <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw,
#'   year = "all", patch = "all", stages = c("stage3", "stage2"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2_raw, 
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw,
#'   year = "all", stages = c("stage3", "stage2"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2_raw, 
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypmatrix3rp <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"), 
#'   size = c("size3added", "size2added", "size1added"), supplement = cypsupp3_raw, 
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
#'   year = "all", stages = c("stage3", "stage2", "stage1"), 
#'   size = c("size3added", "size2added", "size1added"), supplement = cypsupp3_raw, 
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypmatrix2r_3 <- hist_null(cypmatrix2r)
#' cypmatrix2r_3 <- delete_lM(cypmatrix2r_3, year = 2004)
#' diff_r <- diff_lM(cypmatrix3r, cypmatrix2r_3)
#' 
#' cypmatrix2rp_3 <- hist_null(cypmatrix2rp)
#' cypmatrix2rp_3 <- delete_lM(cypmatrix2rp_3, year = 2004)
#' diff_rp <- diff_lM(cypmatrix3rp, cypmatrix2rp_3)
#' 
#' @export
diff_lM <- function(mpm1, mpm2) {
  if (is.null(mpm1) | is.null(mpm2)) {
    stop("Function diff_lM() requires two lefkoMat objects as input.",
      call. = FALSE)
  } else if (all(is.na(mpm1)) | all(is.na(mpm2))) {
    stop("Function diff_lM() requires two lefkoMat objects as input.",
      call. = FALSE)
  }
  if (!is(mpm1, "lefkoMat") | !is(mpm2, "lefkoMat")) {
    stop("Function diff_lM() requires two lefkoMat objects as input.",
      call. = FALSE)
  }
  
  if (length(mpm1$A) != length(mpm2$A)) {
    stop("Objects mpm1 and mpm2 must have the same number of matrices.",
      call. = FALSE)
  }
  if (dim(mpm1$A[[1]])[1] != dim(mpm2$A[[1]])[1]) {
    stop("Objects mpm1 and mpm2 must include matrices of the same dimensions.",
      call. = FALSE)
  }
  
  new_diffs_A <- lapply(c(1:length(mpm1$A)), function(X) {
    newmat <- mpm1$A[[X]] - mpm2$A[[X]]
    
    return(newmat)
  })
  
  new_diffs_U <- lapply(c(1:length(mpm1$A)), function(X) {
    newmat <- mpm1$U[[X]] - mpm2$U[[X]]
    
    return(newmat)
  })
  
  new_diffs_F <- lapply(c(1:length(mpm1$A)), function(X) {
    newmat <- mpm1$F[[X]] - mpm2$F[[X]]
    
    return(newmat)
  })
  
  if (any((mpm1$labels$year2 != mpm2$labels$year2))) {
    warning("Input lefkoMat objects have seemingly different labels objects.",
      call. = FALSE)
  }
  
  output <- list(A = new_diffs_A, U = new_diffs_U, F = new_diffs_F,
    hstages1 = mpm1$hstages, hstages2 = mpm2$hstages, ahstages1 = mpm1$ahstages,
    ahstages2 = mpm2$ahstages, labels1 = mpm1$labels, labels2 = mpm2$labels)
  
  class(output) <- "lefkoDiff"
  
  return(output)
}

#' Summary of Class "hfvdata"
#'
#' A function to simplify the viewing of basic information describing
#' demographic data in historical vertical format (data frames of class
#' \code{hfvdata}).
#' 
#' @name summary_hfv
#' 
#' @param object An object of class \code{hfvdata}.
#' @param popid A string denoting the name of the variable denoting population
#' identity.
#' @param patchid A string denoting the name of the variable denoting patch
#' identity.
#' @param individ A string denoting the name of the variable denoting individual
#' identity.
#' @param year2id A string denoting the name of the variable denoting the year
#' in time \emph{t}.
#' @param full A logical value indicating whether to include basic data frame
#' summary information in addition to hfvdata-specific summary information.
#' Defaults to \code{TRUE}.
#' @param err_check A logical value indicating whether to check for errors in
#' stage assignment.
#' @param ... Other parameters.
#' 
#' @return A summary of the object. The first line shows the numbers of
#' populations, patches, individuals, and time steps. If \code{full = TRUE}, 
#' then this is followed by a standard data frame summary of the hfv dataset.
#' If \code{err_check = TRUE}, then a subset of the original data frame input
#' as \code{object} is exported with only rows showing stage assignment issues.
#' 
#' @section Notes:
#' Stage assignment issue identified by option \code{err_check} fall under two
#' categories. First, all rows showing \code{NoMatch} as the identified stage
#' for \code{stage1}, \code{stage2}, or \code{stage3} are identified. Second,
#' all rows showing \code{stage1 = "NotAlive"} and \code{alive1 = 1},
#' \code{stage2 = "NotAlive"} and \code{alive2 = 1}, or
#' \code{stage3 = "NotAlive"} and \code{alive3 = 1} are identified.
#' 
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' summary_hfv(cypraw_v1)
#' 
#' @export
summary_hfv <- function(object, popid = "popid", patchid = "patchid",
  individ = "individ", year2id = "year2", full = TRUE, err_check = TRUE, ...) {
  
  identified_problems <- NULL
  need_return <- FALSE
  
  demodata <- object
  
  if (!is(demodata, "hfvdata")) {
    stop("This summary function requires an object of class hfvdata as input.",
      call. = FALSE)
  }
  
  matdim <- dim(demodata)
  
  if (!is.element(popid, names(demodata))) {
    stop("Population ID variable not found.", call. = FALSE)
  }
  if (!is.element(patchid, names(demodata))) {
    stop("Patch ID variable not found.", call. = FALSE)
  }
  if (!is.element(individ, names(demodata))) {
    stop("Individual ID variable not found.", call. = FALSE)
  }
  if (!is.element(year2id, names(demodata))) {
    stop("Year at time t ID variable not found.", call. = FALSE)
  }
  
  totalpops <- length(unique(demodata[, popid]))
  totalpatches <- length(unique(demodata[, patchid]))
  totalindivs <- length(unique(demodata[, individ]))
  totalyears <- length(unique(demodata[, year2id]))
  
  grammar_rows <- " rows, "
  grammar_vars <- " variables, "
  grammar_pops <- " populations, "
  grammar_patches <- " patches, "
  grammar_indivs <- " individuals, and "
  grammar_years <- " time steps."
  if (matdim[1] == 1) grammar_rows <- " row, "
  if (matdim[2] == 1) grammar_vars <- " variable, "
  if (totalpops == 1) grammar_pops <- " population, "
  if (totalpatches == 1) grammar_patches <- " patch, "
  if (totalindivs == 1) grammar_indivs <- " individual, and "
  if (totalyears == 1) grammar_years <- " time step."
  
  writeLines(paste0("\nThis hfv dataset contains ", matdim[1], grammar_rows,
      matdim[2], grammar_vars, totalpops, grammar_pops))
  writeLines(paste0(totalpatches, grammar_patches, totalindivs, grammar_indivs,
      totalyears, grammar_years))
  
  if (err_check) {
    stage1_NoMatches <- which(demodata[, "stage1"] == "NoMatch")
    stage2_NoMatches <- which(demodata[, "stage2"] == "NoMatch")
    stage3_NoMatches <- which(demodata[, "stage3"] == "NoMatch")
    
    stage1_NotAlive <- which(demodata[, "stage1"] == "NotAlive")
    stage2_NotAlive <- which(demodata[, "stage2"] == "NotAlive")
    stage3_NotAlive <- which(demodata[, "stage3"] == "NotAlive")
    
    alive1 <- which(demodata[, "alive1"] == 1)
    alive2 <- which(demodata[, "alive2"] == 1)
    alive3 <- which(demodata[, "alive3"] == 1)
    
    s1_NotA_al <- intersect(stage1_NotAlive, alive1)
    s2_NotA_al <- intersect(stage2_NotAlive, alive2)
    s3_NotA_al <- intersect(stage3_NotAlive, alive3)
    
    problem_rows <- sort(unique(c(stage1_NoMatches, stage2_NoMatches,
      stage3_NoMatches, s1_NotA_al, s2_NotA_al, s3_NotA_al)))
    
    if (length(problem_rows) > 0) {
      need_return <- TRUE
      writeLines(paste0("Problems in stage assignment identified in rows:\n"))
      print(problem_rows)
      
      identified_problems <- demodata[problem_rows,]
    }
  }
  
  if (full) {
    dethonthetoilet <- summary.data.frame(demodata)
    print(dethonthetoilet, digits = 3)
  }
  
  if (err_check & need_return) return(identified_problems)
}

