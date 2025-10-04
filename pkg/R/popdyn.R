#' Estimate Stable Stage Distribution
#' 
#' \code{stablestage3()} is a generic function that returns the stable stage 
#' distribution for a population projection matrix or set of matrices. This
#' function is made to handle very large and sparse matrices supplied as 
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as ahistorical
#' matrices.
#' 
#' @name stablestage3
#' 
#' @param mats A lefkoMat object, a population projection matrix, or a list of
#' population projection matrices for which the stable stage distribution is
#' desired.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for details.
#' 
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.list}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' @seealso \code{\link{stablestage3.dgCMatrix}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' # Here we use supplemental() to provide overwrite and reproductive info
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' stablestage3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
stablestage3 <- function(mats, ...) UseMethod("stablestage3")

#' Estimate Stable Stage Distribution of Matrices in lefkoMat Object
#' 
#' \code{stablestage3.lefkoMat()} returns the deterministic stable stage
#' distributions of all \code{A} matrices in an object of class \code{lefkoMat},
#' as well as the long-run projected mean stage distribution in stochastic
#' analysis. This function can handle large and sparse matrices, and so can be
#' used with large historical matrices, IPMs, age x stage matrices, as well as
#' ahistorical matrices.
#' 
#' @name stablestage3.lefkoMat
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of occasions to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param force_sparse A text string indicating whether to use sparse matrix
#' encoding (\code{"yes"}) if standard matrices are provided. Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distributions (and long-run
#' mean stage distributions in stochastic analysis) corresponding to the
#' matrices in a \code{lefkoMat} object.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical, and whether the analysis is deterministic or
#' stochastic. If deterministic and ahistorical, then a single data frame is
#' output, which includes the number of the matrix within the \code{A} element
#' of the input \code{lefkoMat} object, followed by the stage id (numeric and
#' assigned through \code{\link{sf_create}()}), the stage name, and the
#' estimated proportion of the stable stage distribution (\code{ss_prop}). If
#' stochastic and ahistorical, then a single data frame is output starting with
#' the number of the population-patch (\code{matrix_set}), a string
#' concatenating the names of the population and the patch (\code{poppatch}),
#' the assigned stage id number (\code{stage_id}), and the stage name
#' (\code{stage}), and the long-run average stage distribution (\code{ss_prop}).
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{hist} element describes the historical
#' stage-pair distribution, while the \code{ahist} element describes the stage
#' distribution. If deterministic, then \code{hist} contains a data frame
#' including the matrix number (\code{matrix}), the numeric stage designations for
#' stages in occasions \emph{t} and \emph{t}-1, (\code{stage_id_2} and
#' \code{stage_id_1}, respectively), followed by the respective stage names (
#' \code{stage_2} and \code{stage_1}), and ending with the estimated stable
#' stage-pair distribution. The associated \code{ahist} element is as before. If
#' stochastic, then the \code{hist} element contains a single data frame with
#' the number of the population-patch (\code{matrix_set}), a string
#' concatenating the names of the population and the patch (\code{poppatch}),
#' the assigned stage id numbers in times \emph{t} and \emph{t}-1 (
#' \code{stage_id_2} and \code{stage_id_2}, respectively), and the associated
#' stage names (\code{stage_2} and \code{stage_1}, respectively), and the
#' long-run average stage distribution (\code{ss_prop}). The associated
#' \code{ahist} element is as before in the ahistorical, stochastic case.
#'
#' In addition to the data frames noted above, stochastic analysis will result
#' in the additional output of a list of matrices containing the actual
#' projected stage distributions across all projected occasions, in the order of
#' population-patch combinations in the \code{lefkoMat} input.
#'
#' @section Notes:
#' In stochastic analysis, the projected mean distribution is the arithmetic
#' mean across the final 1000 projected occasions if the simulation is at least
#' 2000 projected occasions long. If between 500 and 2000 projected occasions
#' long, then only the final 200 are used, and if fewer than 500 occasions are
#' used, then all are used. Note that because stage distributions in stochastic
#' simulations can change greatly in the initial portion of the run, we
#' encourage a minimum of 2000 projected occasions per simulation, with 10000
#' preferred.
#' 
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.list}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' @seealso \code{\link{stablestage3.dgCMatrix}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' # Here we use supplemental() to provide overwrite and reproductive info
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' stablestage3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
stablestage3.lefkoMat <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, force_sparse = "auto", ...) {
  
  matrix_set <- theprophecy <- NULL
  sparsemethod <- sparse_auto <- sparse_input <- FALSE
  
  if (is(mats$A[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (!sparse_input) {
    if (is.logical(force_sparse)) {
      if (force_sparse) {
        sparsemethod <- TRUE
      } else sparsemethod <- FALSE
    } else if (is.element(tolower(force_sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
      sparsemethod <- TRUE
    } else if (is.element(tolower(force_sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
      sparsemethod <- FALSE
    } else {
      if (is.element(tolower(force_sparse), c("au", "aut", "auto"))) sparse_auto <- TRUE
      
      elements_total <- length(mats$A[[1]])
      dense_elements <- length(which(mats$A[[1]] != 0))
      
      if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
        sparsemethod <- 1
      } else sparsemethod <- 0
    }
  }
  
  if (!stochastic) {
    baldrick <- if (is.matrix(mats$A)) {
      if (!sparse_input) {
        .ss3matrix(mats$A, sparsemethod)
      } else .ss3matrix_sp(mats$A)
      
    } else if (is.list(mats$A)) {
      if (!sparse_input) {
        unlist(lapply(mats$A, .ss3matrix, sparsemethod))
      } else unlist(lapply(mats$A, .ss3matrix_sp))
      
    } else {
      stop("Input not recognized.")
    }
    
  } else {
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    mats$labels$poppatch <- paste(mats$labels$pop, mats$labels$patch)
    used_poppatches <- as.list(unique(mats$labels$poppatch))
    
    #Here we get the full stage distribution series for all occasions, as a list
    princegeorge <- lapply(used_poppatches, function(X) {
      used_slots <- which(mats$labels$poppatch == X)
      
      if (length(used_slots) < 2) {
        warning("Only 1 annual matrix found for some population-patch combinations.
          Stochastic analysis requires multiple annual matrices per population-patch combination.",
          call. = FALSE)
      }
      
      if (!any(is.na(tweights))) {
        tw_error_1 <- "Option tweights must be NA, a numeric vector equal to the"
        tw_error_2 <- "number of years or matrices supplied, or a square matrix"
        tw_error_3 <- "with dimensions equal to the number of matrices per patch."
        tw_error <- paste(tw_error_1, tw_error_2, tw_error_3)
        
        if (is.matrix(tweights)) {
          theprophecy <- rep(0, times)
          used_weights <- tweights[, 1]
        
          for (i in c(1:times)) {
            used_weights <- used_weights / sum(used_weights)
            chosen_timepath <- sample(used_slots, 1, replace = TRUE, prob = used_weights)
            theprophecy[i] <- chosen_timepath[1] - 1
            used_weights <- tweights[, chosen_timepath[1]]
          }
          
          
        } else if (is.numeric(tweights)) {
          if (length(tweights) != length(used_slots)) {
            if (length(tweights) == length(mats$A)) {
              used_weights <- tweights[used_slots] / sum(tweights[used_slots])
            } else {
              stop(tw_error, call. = FALSE)
            }
          } else {
            used_weights <- tweights / sum(tweights)
          }
          
          theprophecy <- sample(used_slots, times, replace = TRUE,
            prob = used_weights) - 1
        } else {
          stop(tw_error, call. = FALSE)
        }
      } else {
        used_weights <- rep(1, length(used_slots))
        used_weights <- used_weights / sum(used_weights)
        
        theprophecy <- sample(used_slots, times, replace = TRUE,
          prob = used_weights) - 1
      }
      
      
      
      
      #theprophecy <- sample(used_slots, times, replace = TRUE, prob = used_weights) - 1
      starter <- if (!sparse_input) {
        .ss3matrix(mats$A[[used_slots[1]]], sparsemethod)
      } else .ss3matrix_sp(mats$A[[used_slots[1]]])
      
      theseventhmatrix <- if (!sparse_input) {
        .proj3(starter, mats$A, theprophecy, 1, 0, 0, sparse_auto, sparsemethod)
      } else .proj3sp(starter, mats$A, theprophecy, 1, 0, 0)
      ssonly <- theseventhmatrix[((dim(mats$A[[1]])[1]) + 1):(2 *(dim(mats$A[[1]])[1])),]
      
      return(ssonly)
    })
    
    # Create mean distributions
    baldrick <- unlist(
      lapply(princegeorge, function(X) {
        if (times > 2000) {
          usedX <- X[,(times - 999):(times)]
        } else if (times > 500) {
          usedX <- X[,(times-199):(times)]
        } else {
          usedX <- X
        }
        apply(usedX, 1, mean)
      })
    )
  }
  
  # The final bits sort everything and clean it up, and create the ahistorical
  # version if a historical input was used
  if (is.list(mats$A)) {
    if (!stochastic) {
      multiplier <- length(mats$A)
    } else {
      multiplier <- length(used_poppatches)
    }
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels_orig <- mats$ahstages[,1:2]
    mat_dims <- dim(mats$A[[1]])[1]
    if (dim(labels_orig)[1] == mat_dims) {
      labels <- labels_orig
    } else {
      newmult <- mat_dims / dim(labels_orig)[1]
      if (mat_dims %% dim(labels_orig)[1] != 0) {
        stop("Matrices do not appear to be ahistorical, historical, or age x stage.
          Cannot proceed. Please make sure that matrix dimensions match stage descriptions.",
          call. = FALSE)
      }
      
      if (!all(is.na(mats$agestages))) {
        check_min <- min(mats$agestages$age)
      }
      age_bit <- c(apply(as.matrix(c((0 + check_min):(newmult + check_min - 1))), 1, rep, dim(labels_orig)[1]))
      core_labels <- cbind.data.frame(age_bit, do.call("rbind.data.frame", replicate(newmult, labels_orig, simplify = FALSE)))
      core_labels$agestage_id <- c(1:length(core_labels$stage_id))
      core_labels$agestage <- apply(as.matrix(core_labels$agestage_id), 1, function(X) {
        paste(core_labels$age_bit[X], core_labels$stage[X])
      })
      labels <- core_labels
      names(labels) <- c("age", "stage_id", "stage", "agestage_id", "agestage")
    }
    
    if (!stochastic) {
      modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      output <- cbind.data.frame(modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix", "age", "stage_id", "stage", "agestage_id", "agestage", "ss_prop")
      } else {
        names(output) <- c("matrix", "stage_id", "stage", "ss_prop")
      }
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      unipoppatches <- unique(modlabels[,1])
      
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      output <- cbind.data.frame(modnums, modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix_set", "poppatch", "age", "stage_id", "stage",
          "agestage_id", "agestage", "ss_prop")
      } else {
        names(output) <- c("matrix_set", "poppatch", "stage_id", "stage",
          "ss_prop")
      }
    }
    rownames(output) <- c(1:dim(output)[1])
  } else {
    labels <- mats$hstages
    
    if (!stochastic) {
      modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      
      outputh <- cbind.data.frame(modlabels, baldrick)
      
      names(outputh) <- c("matrix", "stage_id_2", "stage_id_1", "stage_2",
        "stage_1", "ss_prop")
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      unipoppatches <- unique(modlabels[,1])
      
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      outputh <- cbind.data.frame(modnums, modlabels, baldrick)
      names(outputh) <- c("matrix_set", "poppatch", "stage_id_2", "stage_id_1",
        "stage_2", "stage_1", "ss_prop")
    }
    
    rownames(outputh) <- c(1:dim(outputh)[1])
    
    ahlabels <- mats$ahstages[,c("stage_id", "stage")]
    ss2 <- c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- if (!stochastic) {
        subset(outputh, matrix == X)
      } else {
        subset(outputh, matrix_set == X)
      }
      apply(as.matrix(ahlabels[,1]), 1, function(Y) {
        sum(rightset$ss_prop[which(rightset$stage_id_2 == Y)])
      })
    }))
    
    if (!stochastic) {
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])), 
        rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), ss2)
      names(outputah) <- c("matrix", "stage_id", "stage", "ss_prop")
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(ahlabels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(ahlabels))))
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      outputah <- cbind.data.frame(modnums, modlabels, ss2)
      names(outputah) <- c("matrix_set", "poppatch", "stage_id", "stage", "ss_prop")
    }
    
    rownames(outputah) <- c(1:dim(outputah)[1])
    
    if (!stochastic) {
      output <- list(hist = outputh, ahist = outputah)
    } else {
      output <- list(hist = outputh, ahist = outputah, projections = princegeorge)
    }
  }
  
  return(output)
}

#' Estimate Stable Stage Distribution of a Single Population Projection Matrix
#' 
#' \code{stablestage3.matrix()} returns the stable stage distribution for a 
#' population projection matrix. This function can handle large and sparse
#' matrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @name stablestage3.matrix
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#' @param force_sparse A text string indicating whether to use sparse matrix
#' encoding (\code{"yes"}) when supplied with standard matrices. Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distribution corresponding to
#' the input matrix.
#' 
#' @section Notes:
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.list}()}
#' @seealso \code{\link{stablestage3.dgCMatrix}()}
#' 
#' @examples
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
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean$A[[1]])
#' 
#' @export
stablestage3.matrix <- function(mats, force_sparse = "auto", ...)
{
  sparsemethod <- 0
  sparse_input <- FALSE
  
  if (is(mats, "dgCMatrix")) sparse_input <- TRUE
  
  if (is.logical(force_sparse) & !sparse_input) {
    if (force_sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(force_sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(force_sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!sparse_input) {
    wcorr <- .ss3matrix(mats, sparsemethod)
  } else wcorr <- .ss3matrix_sp(mats)
  
  return(wcorr)
}

#' Estimate Stable Stage Distribution of a Single Population Projection Matrix
#' 
#' \code{stablestage3.dgCMatrix()} returns the stable stage distribution for a 
#' sparse population projection matrix.
#' 
#' @name stablestage3.dgCMatrix
#' 
#' @param mats A population projection matrix of class \code{dgCMatrix}.
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distribution corresponding to
#' the input matrix.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.list}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
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
#' stablestage3(ehrlen3$A[[1]])
#' 
#' @export
stablestage3.dgCMatrix <- function(mats, ...)
{
  
  wcorr <- .ss3matrix_sp(mats)
  
  return(wcorr)
}

#' Estimate Stable Stage Distribution of a List of Projection Matrices
#' 
#' \code{stablestage3.list()} returns the stable stage distributions for stages
#' in population projection matrices arranged in a general list. The function
#' makes no assumptions about whether the matrix is ahistorical and simply
#' provides stable stage distribution values corresponding to each row, meaning
#' that the overall stable stage distribution of basic life history stages in a
#' historical matrix are not provided (the \code{\link{stablestage3.lefkoMat}()}
#' historical estimates these on the basis of stage description information
#' provided in the \code{lefkoMat} object used as input in that function). This
#' provided in the handle large and sparse matrices, and so can be used with
#' large historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @name stablestage3.list
#' 
#' @param mats A list of population projection matrices, all in either class
#' \code{matrix} or class \code{dgCMatrix}.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of occasions to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param force_sparse A text string indicating whether to use sparse matrix
#' encoding (\code{"yes"}) when supplied with standard matrices. Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a list of vector data frames characterizing the 
#' stable stage distributions for stages of each population projection matrix.
#' 
#' @section Notes:
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' @seealso \code{\link{stablestage3.dgCMatrix}()}
#' 
#' @examples
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
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean$A)
#' 
#' @export
stablestage3.list <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, force_sparse = "auto", ...) {
  
  sparsemethod <- 0
  sparse_initial <- FALSE
  w_list <- Xlist <- NULL
  
  if (is(mats[[1]], "dgCMatrix")) sparse_initial <- TRUE
  
  if (!sparse_initial) {
    if (is.logical(force_sparse)) {
      if (force_sparse) {
        sparsemethod <- 1
      } else sparsemethod <- 0
    } else if (is.element(tolower(force_sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
      sparsemethod <- 1
    } else if (is.element(tolower(force_sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
      sparsemethod <- 0
    } else {
      elements_total <- length(mats[[1]])
      dense_elements <- length(which(mats[[1]] != 0))
      
      if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
        sparsemethod <- 1
      } else sparsemethod <- 0
    }
  }
  
  list_length <- length(mats)
  
  massive_panic <- apply(as.matrix(c(1:list_length)), 1, function(X) {
    if (is.matrix(mats[[X]])) {
      outchunk <- 1
    } else if (is(mats[[X]], "dgCMatrix")) {
      outchunk <- 2
    } else outchunk <- 3
  })
  
  if (!all(massive_panic == 1) & !all(massive_panic == 2)) {
    stop("Matrix list must be entirely in standard matrix format, or entirely in dgCMatrix format.",
      call. = FALSE)
  }
  
  if (!stochastic) {
    w_list <- lapply(mats, function(X) {
      if (is.matrix(X)) {
        output <- .ss3matrix(X, sparsemethod)
      } else if (is(X, "dgCMatrix")) {
        output <- .ss3matrix_sp(X)
      } else {
        stop("Unrecognized input in object mats.", FALSE)
      }
    })
  } else {
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    used_slots <- c(1:list_length)
    if (!any(is.na(tweights))) {
      tw_error_1 <- "Option tweights must be NA, a numeric vector equal to the"
      tw_error_2 <- "number of matrices supplied."
      tw_error <- paste(tw_error_1, tw_error_2)
      
      if (is.matrix(tweights)) {
        theprophecy <- rep(0, times)
        used_weights <- tweights[, 1]
      
        for (i in c(1:times)) {
          used_weights <- used_weights / sum(used_weights)
          chosen_timepath <- sample(used_slots, 1, replace = TRUE, prob = used_weights)
          theprophecy[i] <- chosen_timepath[1] - 1
          used_weights <- tweights[, chosen_timepath[1]]
        }
      } else if (is.numeric(tweights)) {
        if (length(tweights) != length(used_slots)) {
          if (length(tweights) == list_length) {
            used_weights <- tweights[used_slots] / sum(tweights[used_slots])
          } else {
            stop(tw_error, call. = FALSE)
          }
        } else {
          used_weights <- tweights / sum(tweights)
        }
        
        theprophecy <- sample(used_slots, times, replace = TRUE,
          prob = used_weights) - 1
      } else {
        stop(tw_error, call. = FALSE)
      }
    } else {
      used_weights <- rep(1, length(used_slots))
      used_weights <- used_weights / sum(used_weights)
      
      theprophecy <- sample(used_slots, times, replace = TRUE,
        prob = used_weights) - 1
    }
    
    starter <- if (!sparse_initial) {
      .ss3matrix(mats[[1]], sparsemethod)
    } else .ss3matrix_sp(mats[[1]])
    
    theseventhmatrix <- if (!sparse_initial) {
      .proj3(starter, mats, theprophecy, 1, 0, 0, TRUE, sparsemethod)
    } else .proj3sp(starter, mats, theprophecy, 1, 0, 0)
    ssonly <- theseventhmatrix[((dim(mats[[1]])[1]) + 1):(2 *(dim(mats[[1]])[1])),]
    
    if (times > 2000) {
      Xlist <- ssonly[,(times - 999):(times)]
    } else if (times > 500) {
      Xlist <- ssonly[,(times-199):(times)]
    } else {
      Xlist <- ssonly
    }
    w_list <- apply(Xlist, 1, mean)
  }
  
  return(w_list)
}

#' Estimate Reproductive Value
#' 
#' \code{repvalue3()} is a generic function that estimates returns the
#' reproductive values of stages in a population projection matrix or a set of
#' matrices. The specifics of estimation vary with the class of input object.
#' This function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as ahistorical
#' matrices.
#' 
#' @name repvalue3
#' 
#' @param mats A lefkoMat object, a population projection matrix, or a list of
#' population projection matrices for which the reproductive value vector is
#' desired.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for details.
#' 
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' @seealso \code{\link{repvalue3.dgCMatrix}()}
#' @seealso \code{\link{repvalue3.list}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' # Here we use supplemental() to provide overwrite and reproductive info
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' repvalue3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
repvalue3 <- function(mats, ...) UseMethod("repvalue3")

#' Estimate Reproductive Value Vectors of Matrices in a lefkoMat Object
#' 
#' \code{repvalue3.lefkoMat()} returns the reproductive values for stages in a
#' set of population projection matrices provided as a \code{lefkoMat} object.
#' This function can handle large and sparse matrices, and so can be used with
#' large historical matrices, IPMs, age x stage matrices, as well as ahistorical
#' matrices.
#' 
#' @name repvalue3.lefkoMat
#' 
#' @param mats An object of class \code{lefkoMat} object.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of occasions to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param seed A number to use as a random number seed.
#' @param force_sparse A text string indicating whether to use sparse matrix
#' encoding (\code{"yes"}) when supplied with standard matrices. Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns the asymptotic reproductive value vectors if
#' deterministic analysis is chosen, and long-run mean reproductive value
#' vectors if stochastic analysis is chosen.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical, and whether the analysis is deterministic or
#' stochastic. If deterministic and ahistorical, then a single data frame is
#' output, which includes the number of the matrix within the \code{A} element
#' of the input \code{lefkoMat} object, followed by the stage id (numeric and
#' assigned through \code{\link{sf_create}()}), the stage name, and the
#' estimated proportion of the reproductive value vector (\code{rep_value}). If
#' stochastic and ahistorical, then a single data frame is output starting with
#' the number of the population-patch (\code{matrix_set}), a string
#' concatenating the names of the population and the patch (\code{poppatch}),
#' the assigned stage id number (\code{stage_id}), and the stage name
#' (\code{stage}), and the long-run mean reproductive value vector
#' (\code{rep_value}).
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{hist} element describes the historical
#' stage-pair reproductive values, while the \code{ahist} element describes the
#' stage reproductive values. If deterministic, then \code{hist} contains a data
#' frame including the matrix number (\code{matrix}), the numeric stage
#' designations for stages in occasions \emph{t} and \emph{t}-1,
#' (\code{stage_id_2} and \code{stage_id_1}, respectively), followed by the
#' respective stage names (\code{stage_2} and \code{stage_1}), and ending with
#' the estimated reproductive values (\code{rep_value}). The associated
#' \code{ahist} element is as before. If stochastic, then the \code{hist}
#' element contains a single data frame with the number of the population-patch
#' (\code{matrix_set}), a string concatenating the names of the population and
#' the patch (\code{poppatch}), the assigned stage id numbers in times \emph{t}
#' and \emph{t}-1 (\code{stage_id_2} and \code{stage_id_2}, respectively), and
#' the associated stage names (\code{stage_2} and \code{stage_1}, respectively),
#' and the long-run mean reproductive values (\code{rep_value}). The associated
#' \code{ahist} element is as before in the ahistorical, stochastic case.
#'
#' In addition to the data frames noted above, stochastic analysis will result
#' in the additional output of a list of matrices containing the actual
#' projected reproductive value vectors across all projected occasions, in the
#' order of population-patch combinations in the \code{lefkoMat} input.
#'
#' @section Notes:
#' In stochastic analysis, the projected mean reproductive value vector is the
#' arithmetic mean across the final projected 1000 occasions if the simulation
#' is at least 2000 projected occasions long. If between 500 and 2000 projected
#' occasions long, then only the final 200 are used, and if fewer than 500
#' occasions are used, then all are used. Note that because reproductive values
#' in stochastic simulations can change greatly in the initial portion of the
#' run, we encourage a minimum 2000 projected occasions per simulation, with
#' 10000 preferred.
#' 
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have several hundred rows and columns. Defaults
#' work best when matrices are very small and dense, or very large and sparse.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' @seealso \code{\link{repvalue3.dgCMatrix}()}
#' @seealso \code{\link{repvalue3.list}()}
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
#' # Here we use supplemental() to provide overwrite and reproductive info
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' repvalue3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
repvalue3.lefkoMat <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, force_sparse = "auto", ...) {
  
  matrix_set <- poppatch <- NULL
  sparsemethod <- sparse_auto <- sparse_input <- FALSE
  
  if (is(mats$A[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (!sparse_input) {
    if (is.logical(force_sparse)) {
      if (force_sparse) {
        sparsemethod <- TRUE
      } else sparsemethod <- FALSE
    } else if (is.element(tolower(force_sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
      sparsemethod <- TRUE
    } else if (is.element(tolower(force_sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
      sparsemethod <- FALSE
    } else {
      if (is.element(tolower(force_sparse), c("au", "aut", "auto"))) sparse_auto <- TRUE
      
      elements_total <- length(mats$A[[1]])
      dense_elements <- length(which(mats$A[[1]] != 0))
      
      if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
        sparsemethod <- 1
      } else sparsemethod <- 0
    }
  }
  
  if (!stochastic) {
    baldrick <- if (is.matrix(mats$A)) {
      
      almost_final <- if (!sparse_input) {
        .rv3matrix(mats$A, sparsemethod)
      } else .rv3matrix_sp(mats$A)
      
    } else if (is.list(mats$A)) {
      final <- unlist(lapply(mats$A, function(X) {
        almost_final <- if (!sparse_input) {
          .rv3matrix(X, sparsemethod)
        } else .rv3matrix_sp(X)
        
        return(almost_final/almost_final[which(almost_final == min(almost_final[which(almost_final > 0)])[1])[1]])
      }))
    } else {
      stop("Input not recognized.")
    }
  } else {
    
    if (!is.na(seed)) {
      set.seed(seed)
    }
    
    mats$labels$poppatch <- paste(mats$labels$pop, mats$labels$patch)
    used_poppatches <- as.list(unique(mats$labels$poppatch))
    
    # Full stage distribution and reproductive value vector series for all occasions, as list
    # Stage distributions in top half of matrix; reproductive value is at bottom
    princegeorge <- lapply(used_poppatches, function(X) {
      used_slots <- which(mats$labels$poppatch == X)
      
      if (length(used_slots) < 2) {
        warning("Only 1 annual matrix found for some pop-patch combinations.
          Stochastic analysis requires multiple annual matrices per pop-patch combination.",
          call. = FALSE)
      }
      
      if (!any(is.na(tweights))) {
        tw_error_1 <- "Option tweights must be NA, a numeric vector equal to the"
        tw_error_2 <- "number of years or matrices supplied, or a square matrix"
        tw_error_3 <- "with dimensions equal to the number of matrices per patch."
        tw_error <- paste(tw_error_1, tw_error_2, tw_error_3)
        
        if (is.matrix(tweights)) {
          theprophecy <- rep(0, times)
          used_weights <- tweights[, 1]
        
          for (i in c(1:times)) {
            used_weights <- used_weights / sum(used_weights)
            chosen_timepath <- sample(used_slots, 1, replace = TRUE, prob = used_weights)
            theprophecy[i] <- chosen_timepath[1] - 1
            used_weights <- tweights[, chosen_timepath[1]]
          }
          
        } else if (is.numeric(tweights)) {
          if (length(tweights) != length(used_slots)) {
            if (length(tweights) == length(mats$A)) {
              used_weights <- tweights[used_slots] / sum(tweights[used_slots])
            } else {
              stop(tw_error, call. = FALSE)
            }
          } else {
            used_weights <- tweights / sum(tweights)
          }
          
          theprophecy <- sample(used_slots, times, replace = TRUE,
            prob = used_weights) - 1
        } else {
          stop(tw_error, call. = FALSE)
        }
      } else {
        used_weights <- rep(1, length(used_slots))
        used_weights <- used_weights / sum(used_weights)
        
        theprophecy <- sample(used_slots, times, replace = TRUE,
          prob = used_weights) - 1
      }
      
      if (!sparse_input) {
        starter <- .ss3matrix(mats$A[[used_slots[1]]], sparsemethod)
      } else starter <- .ss3matrix_sp(mats$A[[used_slots[1]]])
      
      theseventhmatrix <- if (!sparse_input) {
        .proj3(starter, mats$A, theprophecy, 1, 0, 0, sparse_auto, sparsemethod)
      } else .proj3sp(starter, mats$A, theprophecy, 1, 0, 0)
      
      almostall <- theseventhmatrix[((dim(mats$A[[1]])[1]) + 1):(3 * (dim(mats$A[[1]])[1])),]
      
      return(almostall)
    })
    
    # Create mean distributions
    baldrick <- unlist(
      lapply(princegeorge, function(X) {
        if (times > 2000) {
          usedX1 <- X[1:(dim(mats$A[[1]])[1]), (times - 998):(times+1)]
          usedX2 <- X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])),1:1000]
        } else if (times > 500) {
          usedX1 <- X[1:(dim(mats$A[[1]])[1]), (times - 198):(times+1)]
          usedX2 <- X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])), 1:200]
        } else {
          usedX1 <- X[1:(dim(mats$A[[1]])[1]),]
          usedX2 <- X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])),]
        }
        meanX1 <- apply(usedX1, 1, mean)
        meanX2 <- apply(usedX2, 1, mean)
        meanX2 <- zapsmall(meanX2)
        meanX2 <- meanX2 / meanX2[(which(meanX2 > 0)[1])]
        meanX <- c(meanX1, meanX2)
        
        return(meanX)
      })
    )
    
    princegeorge <- lapply(princegeorge, function(X) {
      return(X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])),])
    })
  }
  
  if (is.list(mats$A)) {
    if (!stochastic) {
      multiplier <- length(mats$A)
    } else {
      multiplier <- length(used_poppatches)
    }
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels_orig <- mats$ahstages[,1:2]
    mat_dims <- dim(mats$A[[1]])[1]
    if (dim(labels_orig)[1] == mat_dims) {
      labels <- labels_orig
    } else {
      newmult <- mat_dims / dim(labels_orig)[1]
      if (mat_dims %% dim(labels_orig)[1] != 0) {
        stop("Matrices do not appear to be ahistorical, historical, or age x stage.
          Cannot proceed. Please make sure that matrix dimensions match stage descriptions.",
          call. = FALSE)
      }
      
      if (!all(is.na(mats$agestages))) {
        check_min <- min(mats$agestages$age)
      }
      age_bit <- c(apply(as.matrix(c((0 + check_min):(newmult + check_min - 1))), 1, rep, dim(labels_orig)[1]))
      core_labels <- cbind.data.frame(age_bit, do.call("rbind.data.frame",
        replicate(newmult, labels_orig, simplify = FALSE)))
      core_labels$agestage_id <- c(1:length(core_labels$stage_id))
      core_labels$agestage <- apply(as.matrix(core_labels$agestage_id), 1, function(X) {
        paste(core_labels$age_bit[X], core_labels$stage[X])
      })
      labels <- core_labels
      names(labels) <- c("age", "stage_id", "stage", "agestage_id", "agestage")
    }
    
    if (!stochastic) {
      modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      output <- cbind.data.frame(modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix", "age", "stage_id", "stage", "agestage_id", "agestage", "rep_value")
      } else {
        names(output) <- c("matrix", "stage_id", "stage", "rep_value")
      }
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      unipoppatches <- unique(modlabels[,1])
      
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      output <- cbind.data.frame(modnums, modlabels, baldrick[(mat_dims+1):(2*mat_dims)])
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix_set", "poppatch", "age", "stage_id", "stage", "agestage_id", "agestage", "rep_value")
      } else {
        names(output) <- c("matrix_set", "poppatch", "stage_id", "stage", "rep_value")
      }
    }
    
    rownames(output) <- c(1:dim(output)[1])
    
  } else {
    
    # This section translates historical results to ahistorical and then cleans everything up
    if (!stochastic) {
      ss3 <- stablestage3.lefkoMat(mats, force_sparse = force_sparse) #stablestage3.lefkoMat
      rahist <- ss3$ahist
      rhist <-ss3$hist
      rhist$ss3sum <- apply(as.matrix(c(1:dim(rhist)[1])), 1, function(X) {
        rahist$ss_prop[intersect(which(rahist$stage_id == rhist$stage_id_2[X]), 
          which(rahist$matrix == rhist$matrix[X]))]
      })
      rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
      rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
      rhist$rep_value <- baldrick
      
      rhist$rv3raw <- rhist$sscorr * rhist$rep_value
      
      outputh <- rhist[,c("matrix", "stage_id_2", "stage_id_1", "stage_2", "stage_1", "rep_value")]
      
      ahlabels <- mats$ahstages[,c("stage_id", "stage")]
      rv2 <- Re(c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
        rightset <- subset(rhist, matrix == X)
        apply(as.matrix(ahlabels[,1]), 1, function(Y) {
          sum(rightset$rv3raw[which(rightset$stage_id_2 == Y)])
        })
      })))
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])),
        rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), rv2)
      names(outputah) <- c("matrix", "stage_id", "stage", "rep_value_unc")
      
      outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
        matsub <- subset(outputah, matrix == outputah$matrix[X])
        entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
        return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
      })
      
      outputah <- outputah[,c("matrix", "stage_id", "stage", "rep_value")]
      rownames(outputah) <- c(1:dim(outputah)[1])
      
    } else {
      
      ss_unlisted <- apply(as.matrix(c(1:multiplier)), 1, function(X) {
        return(baldrick[((2 * (X - 1) * (dim(mats$A[[1]])[1])) + 1): (((2 * (X - 1)) + 1) * (dim(mats$A[[1]])[1]))])
      })
      rv_unlisted <- apply(as.matrix(c(1:multiplier)), 1, function(X) {
        return(baldrick[((((2 * (X - 1)) + 1) * (dim(mats$A[[1]])[1])) + 1): (X * 2 * (dim(mats$A[[1]])[1]))])
      })
      
      ss_sums <- apply(as.matrix(c(1:multiplier)), 1, function(X) {
        sum_vecs <- apply(as.matrix(mats$hstages$stage_id_2), 1, function(Y) {
          sum(ss_unlisted[(which(mats$hstages$stage_id_2 == Y)) + ((X - 1) * dim(mats$hstages)[1])])
        })
        return(sum_vecs)
      })
      
      ahistsize <- length(mats$ahstages$stage_id)
      histsize <- length(mats$hstages$stage_id_2)
      
      poppatchcol <- rep(unlist(used_poppatches), each = ahistsize)
      unipoppatches <- unique(poppatchcol)
      poppatchnum <- apply(as.matrix(poppatchcol), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      rahist <- cbind.data.frame(matrix_set = poppatchnum, poppatch = poppatchcol, 
        stage_id = rep(mats$ahstages$stage_id, times = multiplier),
        stage = rep(mats$ahstages$stage, times = multiplier))
      
      poppatchcol <- rep(unlist(used_poppatches), each = histsize)
      poppatchnum <- apply(as.matrix(poppatchcol), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      rhist <- cbind.data.frame(matrix_set = poppatchnum, poppatch = poppatchcol,
        stage_id_2 = rep(mats$hstages$stage_id_2, times = multiplier),
        stage_id_1 = rep(mats$hstages$stage_id_1, times = multiplier),
        stage_2 = rep(mats$hstages$stage_2, times = multiplier),
        stage_1 = rep(mats$hstages$stage_1, times = multiplier))
      
      rhist$ss_prop <- as.vector(ss_unlisted)
      rhist$ss3sum <- as.vector(ss_sums)
      rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
      rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
      rhist$rep_value <- as.vector(rv_unlisted)
      
      rhist$rv3raw <- rhist$sscorr * rhist$rep_value
      
      outputh <- rhist[,c("matrix_set", "poppatch", "stage_id_2", "stage_id_1",
        "stage_2", "stage_1", "rep_value")]
      
      ahlabels <- mats$ahstages[,c("stage_id", "stage")]
      rv2 <- Re(c(apply(as.matrix(unlist(used_poppatches)), 1, function(X) {
        rightset <- subset(rhist, poppatch == X)
        apply(as.matrix(ahlabels[,1]), 1, function(Y) {
          sum(rightset$rv3raw[which(rightset$stage_id_2 == Y)])
        })
      })))
      
      poppatchcol <- rep(unlist(used_poppatches), each = ahistsize)
      
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])),
        poppatchcol, rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), rv2)
      names(outputah) <- c("matrix_set", "poppatch", "stage_id", "stage",
        "rep_value_unc")
      
      outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
        matsub <- subset(outputah, poppatch == outputah$poppatch[X])
        
        if (!all(matsub$rep_value_unc == 0)) {
          entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
          return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
        } else return(0)
      })
      
      outputah <- outputah[,c("matrix_set", "poppatch", "stage_id", "stage", "rep_value")]
      rownames(outputah) <- c(1:dim(outputah)[1])
    }
    
    if (!stochastic) {
      output <-list(hist = outputh, ahist = outputah)
    } else {
      output <-list(hist = outputh, ahist = outputah, projections = princegeorge)
    }
  }
  return(output)
}

#' Estimate Reproductive Value Vector for a Single Population Projection Matrix
#' 
#' \code{repvalue3.matrix()} returns the reproductive values for stages in a 
#' population projection matrix. The function makes no assumptions about whether
#' the matrix is ahistorical and simply provides standard reproductive values
#' corresponding to each row, meaning that the overall reproductive values of
#' basic life history stages in a historical matrix are not provided (the 
#' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis of
#' stage description information provided in the \code{lefkoMat} object used as
#' input in that function). This function can handle large and sparse matrices,
#' and so can be used with large historical matrices, IPMs, age x stage
#' matrices, as well as smaller ahistorical matrices.
#' 
#' @name repvalue3.matrix
#' 
#' @param mats A population projection matrix.
#' @param force_sparse A text string indicating whether to use sparse matrix
#' encoding (\code{"yes"}) when supplied with standard matrices. Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a vector data frame characterizing the 
#' reproductive values for stages of a population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue, divided by the first non-zero element of the left 
#' eigenvector.
#' 
#' @section Notes:
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.dgCMatrix}()}
#' @seealso \code{\link{repvalue3.list}()}
#' 
#' @examples
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
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean$A[[1]])
#' 
#' @export
repvalue3.matrix <- function(mats, force_sparse = "auto", ...)
{
  sparsemethod <- 0
  
  if (is.logical(force_sparse)) {
    if (force_sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(force_sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(force_sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  v <- .rv3matrix(mats, sparsemethod)

  return(v)
}

#' Estimate Reproductive Value Vector for a Single Population Projection Matrix
#' 
#' \code{repvalue3.dgCMatrix()} returns the reproductive values for stages in a 
#' sparse population projection matrix. The function makes no assumptions about
#' whether the matrix is ahistorical and simply provides standard reproductive
#' values corresponding to each row, meaning that the overall reproductive
#' values of basic life history stages in a historical matrix are not provided
#' (the \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis
#' of stage description information provided in the \code{lefkoMat} object used
#' as input in that function).
#' 
#' @name repvalue3.dgCMatrix
#' 
#' @param mats A population projection matrix.
#' @param ... Other parameters.
#' 
#' @return This function returns a vector data frame characterizing the 
#' reproductive values for stages of a population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue, divided by the first non-zero element of the left 
#' eigenvector.
#' 
#' @section Notes:
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have several hundred rows and columns. Defaults
#' work best when matrices are very small and dense, or very large and sparse.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' @seealso \code{\link{repvalue3.list}()}
#' 
#' @examples
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
#' repvalue3(ehrlen3$A[[1]])
#' 
#' @export
repvalue3.dgCMatrix <- function(mats, ...)
{
  
  v <- .rv3matrix_sp(mats)
  
  return(v)
}

#' Estimate Reproductive Value Vector for a List of Projection Matrices
#' 
#' \code{repvalue3.list()} returns the reproductive values for stages in
#' population projection matrices arranged in a general list. The function makes
#' no assumptions about whether the matrix is ahistorical and simply provides
#' standard reproductive values corresponding to each row, meaning that the
#' overall reproductive values of basic life history stages in a historical
#' matrix are not provided (the \code{\link{repvalue3.lefkoMat}()} function
#' estimates these on the basis of stage description information provided in the
#' \code{lefkoMat} object used as input in that function). This function can
#' handle large and sparse matrices, and so can be used with large historical
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical
#' matrices.
#' 
#' @name repvalue3.list
#' 
#' @param mats A list of population projection matrices, all in either class
#' \code{matrix} or class \code{dgCMatrix}.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of occasions to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param force_sparse A text string indicating whether to use sparse matrix
#' encoding (\code{"yes"}) when supplied with standard matrices. Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a list of vector data frames characterizing the 
#' reproductive values for stages of each population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue, divided by the first non-zero element of the left 
#' eigenvector.
#' 
#' @section Notes:
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have several hundred rows and columns. Defaults
#' work best when matrices are very small and dense, or very large and sparse.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.dgCMatrix}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
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
#' repvalue3(ehrlen3$A)
#' 
#' @export
repvalue3.list <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, force_sparse = "auto", ...) {
  
  sparsemethod <- 0
  sparse_initial <- FALSE
  v_list <- Xlist <- NULL
  
  if (is(mats[[1]], "dgCMatrix")) sparse_initial <- TRUE
  
  if (!sparse_initial) {
    if (is.logical(force_sparse)) {
      if (force_sparse) {
        sparsemethod <- 1
      } else sparsemethod <- 0
    } else if (is.element(tolower(force_sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
      sparsemethod <- 1
    } else if (is.element(tolower(force_sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
      sparsemethod <- 0
    } else {
      elements_total <- length(mats[[1]])
      dense_elements <- length(which(mats[[1]] != 0))
      
      if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
        sparsemethod <- 1
      } else sparsemethod <- 0
    }
  }
  
  list_length <- length(mats)
  
  massive_panic <- apply(as.matrix(c(1:list_length)), 1, function(X) {
    if (is.matrix(mats[[X]])) {
      outchunk <- 1
    } else if (is(mats[[X]], "dgCMatrix")) {
      outchunk <- 2
    } else outchunk <- 3
  })
  
  if (!all(massive_panic == 1) & !all(massive_panic == 2)) {
    stop("Matrix list must be entirely in standard matrix format, or entirely in dgCMatrix format.",
      call. = FALSE)
  }
  
  if (!stochastic) {
    v_list <- lapply(mats, function(X) {
      if (is.matrix(X)) {
        output <- .rv3matrix(X, sparsemethod)
      } else if (is(X, "dgCMatrix")) {
        output <- .rv3matrix_sp(X)
      } else {
        stop("Unrecognized input in object mats.", FALSE)
      }
    })
  } else {
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    used_slots <- c(1:list_length)
    if (!any(is.na(tweights))) {
      tw_error_1 <- "Option tweights must be NA, a numeric vector equal to the"
      tw_error_2 <- "number of matrices supplied."
      tw_error <- paste(tw_error_1, tw_error_2)
      
      if (is.matrix(tweights)) {
        theprophecy <- rep(0, times)
        used_weights <- tweights[, 1]
      
        for (i in c(1:times)) {
          used_weights <- used_weights / sum(used_weights)
          chosen_timepath <- sample(used_slots, 1, replace = TRUE, prob = used_weights)
          theprophecy[i] <- chosen_timepath[1] - 1
          used_weights <- tweights[, chosen_timepath[1]]
        }
      } else if (is.numeric(tweights)) {
        if (length(tweights) != length(used_slots)) {
          if (length(tweights) == list_length) {
            used_weights <- tweights[used_slots] / sum(tweights[used_slots])
          } else {
            stop(tw_error, call. = FALSE)
          }
        } else {
          used_weights <- tweights / sum(tweights)
        }
        
        theprophecy <- sample(used_slots, times, replace = TRUE,
          prob = used_weights) - 1
      } else {
        stop(tw_error, call. = FALSE)
      }
    } else {
      used_weights <- rep(1, length(used_slots))
      used_weights <- used_weights / sum(used_weights)
      
      theprophecy <- sample(used_slots, times, replace = TRUE,
        prob = used_weights) - 1
    }
    
    starter <- if (!sparse_initial) {
      .rv3matrix(mats[[1]], sparsemethod)
    } else .rv3matrix_sp(mats[[1]])
    
    theseventhmatrix <- if (!sparse_initial) {
      .proj3(starter, mats, theprophecy, 1, 0, 0, TRUE, sparsemethod)
    } else .proj3sp(starter, mats, theprophecy, 1, 0, 0)
    almostall <- theseventhmatrix[((dim(mats[[1]])[1]) + 1):(3 * (dim(mats[[1]])[1])),]
    
    if (times > 2000) {
      usedX1 <- almostall[1:(dim(mats[[1]])[1]), (times - 998):(times+1)]
      usedX2 <- almostall[((dim(mats[[1]])[1]) + 1):(2*(dim(mats[[1]])[1])),1:1000]
    } else if (times > 500) {
      usedX1 <- almostall[1:(dim(mats[[1]])[1]), (times - 198):(times+1)]
      usedX2 <- almostall[((dim(mats[[1]])[1]) + 1):(2*(dim(mats[[1]])[1])), 1:200]
    } else {
      usedX1 <- almostall[1:(dim(mats[[1]])[1]),]
      usedX2 <- almostall[((dim(mats[[1]])[1]) + 1):(2*(dim(mats[[1]])[1])),]
    }
    meanX1 <- apply(usedX1, 1, mean)
    meanX2 <- apply(usedX2, 1, mean)
    meanX2 <- zapsmall(meanX2)
    meanX2 <- meanX2 / meanX2[(which(meanX2 > 0)[1])]
    v_list <- c(meanX1, meanX2)
  }
  
  return(v_list)
}

#' Estimate Sensitivity of Population Growth Rate to Matrix Elements
#' 
#' \code{sensitivity3()} is a generic function that returns the sensitivity of
#' the population growth rate to the elements of the matrices in a matrix
#' population model. Currently, this function estimates both deterministic and
#' stochastic sensitivities, where the growth rate is \eqn{\lambda} in the
#' former case and the log of the stochastic \eqn{\lambda} in the latter case.
#' This function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects, as lists of matrices, and as individual matrices.
#' 
#' @name sensitivity3
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which
#' the stable stage distribution is desired.
#' @param ... Other parameters
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' @seealso \code{\link{sensitivity3.dgCMatrix}()}
#' @seealso \code{\link{sensitivity3.list}()}
#' 
#' @examples
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
#' sensitivity3(ehrlen3)
#' 
#' @export
sensitivity3 <- function(mats, ...) UseMethod("sensitivity3")

#' Estimate Sensitivity of Population Growth Rate of a lefkoMat Object
#' 
#' \code{sensitivity3.lefkoMat()} returns the sensitivities of population growth
#' rate to elements of all \code{$A} matrices in an object of class
#' \code{lefkoMat}. If deterministic, then \eqn{\lambda} is taken as the
#' population growth rate. If stochastic, then the log of stochastic
#' \eqn{\lambda}, or the log stochastic growth rate, is taken as the population
#' growth rate. This function can handle large and sparse matrices, and so can
#' be used with large historical matrices, IPMs, age x stage matrices, as well
#' as smaller ahistorical matrices.
#' 
#' @name sensitivity3.lefkoMat
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) sensitivity analysis. Defaults to
#' FALSE.
#' @param times The number of occasions to project forward in stochastic
#' simulation. Defaults to \code{10000}.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param append_mats A logical value indicating whether to include the original
#' A, U, and F matrices in the output \code{lefkoSens} object.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoSens}, which is a
#' list of 8 elements. The first, \code{h_sensmats}, is a list of historical
#' sensitivity matrices (\code{NULL} if an ahMPM is used as input). The second,
#' \code{ah_elasmats}, is a list of either ahistorical sensitivity matrices if
#' an ahMPM is used as input, or, if an hMPM is used as input, then the result
#' is a list of ahistorical matrices based on the equivalent historical
#' dependencies assumed in the input historical matrices. The third element,
#' \code{hstages}, is a data frame showing historical stage pairs (\code{NULL}
#' if an ahMPM used as input). The fourth element, \code{agestages}, show the
#' order of age-stage combinations, if age-by-stage MPMs have been supplied. The
#' fifth element, \code{ahstages}, is a data frame showing the order of
#' ahistorical stages. The last 3 elements are the A, U, and F portions of the
#' input.
#' 
#' @section Notes:
#' All sensitivity matrix outputs from this function are in standard matrix
#' format.
#' 
#' Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
#' Matrix Population Models). Stochastic sensitivities are estimated as eqn.
#' 14.97 in Caswell (2001). Note that stochastic sensitivities are of the log of
#' the stochastic \eqn{\lambda}.
#'
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' The \code{time_weights} and \code{steps} arguments are now deprecated.
#' Instead, please use the \code{tweights} and \code{times} arguments.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' @seealso \code{\link{sensitivity3.dgCMatrix}()}
#' @seealso \code{\link{sensitivity3.list}()}
#' 
#' @examples
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
#' sensitivity3(ehrlen3, stochastic = TRUE)
#' 
#' @export
sensitivity3.lefkoMat <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, sparse = "auto", append_mats = FALSE, ...) {
  
  sparsemethod <- 0
  sparse_input <- FALSE
  
  possible_args <- c("mats", "stochastic", "times", "tweights", "historical",
    "seed", "force_sparse", "append_mats", "time_weights", "steps", "sparse")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  if ("time_weights" %in% further_args_names) {
    stop("Argument time_weights is deprecated. Please use argument tweights instead",
      call. = FALSE)
  }
  if ("steps" %in% further_args_names) {
    stop("Argument steps is deprecated. Please use argument times instead",
      call. = FALSE)
  }
  
  if (is(mats$A[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- if (!sparse_input ) {
      length(which(mats$A[[1]] != 0))
    } else {
      sum(diff(mats$A[[1]]@p))
    }
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (stochastic & length(mats$A) < 2) {
    stop("Stochastic sensitivity estimation cannot be completed with fewer than 2 annual matrices.",
      call. = FALSE)
  }
  
  if (!is.element("agestages", names(mats))) {
    mats$agestages <- NA
  }
  
  if (!stochastic) {
    # Deterministic sensitivity analysis
    message("Running deterministic analysis...")
    
    baldrick <- if (all(is.na(mats$hstages))) {
      if (!sparse_input) {
        lapply(mats$A, .sens3matrix, sparsemethod)
      } else {
        lapply(mats$A, .sens3matrix_spinp)
      }
    } else {
      if (!sparse_input) {
        lapply(mats$A, .sens3hlefko, mats$ahstages, mats$hstages)
      } else {
        lapply(mats$A, .sens3hlefko_sp, mats$ahstages, mats$hstages)
      }
    }
    
    hlabels <- mats$hstages
    ahlabels <- mats$ahstages
    agelabels <- mats$agestages
    new_labels <- mats$labels
    
    if (all(is.na(mats$hstages))) {
      
      output <- list(h_sensmats = NULL, ah_sensmats = baldrick, hstages = NULL, 
        ahstages = ahlabels, agestages = agelabels, labels = new_labels)
      
    } else {
      he_list <- lapply(baldrick, function(X) {X$h_smat})
      ahe_list <- lapply(baldrick, function(X) {X$ah_smat})
      
      output <- list(h_sensmats = he_list, ah_sensmats = ahe_list,
        hstages = hlabels, ahstages = ahlabels, agestages = agelabels,
        labels = new_labels)
    }
  } else {
    # Stochastic sensitivity analysis
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    if(!any(is.na(tweights))) {
      message("Running stochastic analysis with tweights option...")
      
      returned_cubes <- .stoch_senselas(mats, times = times, historical = FALSE,
        style = 1, sparsemethod, lefkoProj = TRUE, tweights = tweights) 
    } else {
      message("Running stochastic analysis...")
      
      returned_cubes <- .stoch_senselas(mats, times = times, historical = FALSE,
        style = 1, sparsemethod, lefkoProj = TRUE)
    }
    
    old_labels <- mats$labels
    old_labels$concat <- paste(old_labels$pop, old_labels$patch)
    labs_unique <- unique(old_labels$concat)
    lab_indices <- apply(as.matrix(labs_unique), 1, function(X) {
      found_guys <- which(old_labels$concat == X)
      min_guy <- min(found_guys)
      
      return(min_guy)
    })
    new_labels <- mats$labels[lab_indices, c(1, 2)]
    
    if (!all(is.na(mats$hstages))) {
      output <- list(h_sensmats = returned_cubes[[1]], ah_sensmats = returned_cubes[[2]],
        hstages = mats$hstages, agestages = mats$agestages,
        ahstages = mats$ahstages, labels = new_labels)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = returned_cubes[[1]],
        hstages = mats$hstages, agestages = mats$agestages,
        ahstages = mats$ahstages, labels = new_labels)
    }
  }
  
  if (append_mats) {
    output$A <- mats$A
    output$U <- mats$U
    output$F <- mats$F
  }
  
  class(output) <- "lefkoSens"
  
  return(output)
}

#' Estimate Sensitivity of Population Growth Rate of a Single Matrix
#' 
#' \code{sensitivity3.matrix()} returns the sensitivities of \eqn{\lambda} to
#' elements of a single matrix. Because this handles only one matrix, the
#' sensitivities are inherently deterministic and based on the dominant eigen
#' value as the best metric of the population growth rate. This function can
#' handle large and sparse matrices, and so can be used with large historical
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical
#' matrices.
#' 
#' @name sensitivity3.matrix
#' 
#' @param mats An object of class \code{matrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a single deterministic sensitivity matrix.
#' 
#' @section Notes:
#' All sensitivity matrix outputs from this function are in standard matrix
#' format.
#' 
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.dgCMatrix}()}
#' @seealso \code{\link{sensitivity3.list}()}
#' 
#' @examples
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
#' ehrlen3mean <- lmean(ehrlen3)
#' sensitivity3(ehrlen3mean$A[[1]])
#' 
#' @export
sensitivity3.matrix <- function(mats, sparse = "auto", ...)
{
  sparsemethod <- 0
  
  possible_args <- c("mats", "sparse")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  wcorr <- .sens3matrix(mats, FALSE)
  
  if (sparsemethod == 1) wcorr <- as(wcorr, "dgCMatrix")  
  
  return(wcorr)
}

#' Estimate Sensitivity of Population Growth Rate of a Single Matrix
#' 
#' \code{sensitivity3.dgCMatrix()} returns the sensitivities of \eqn{\lambda} to
#' elements of a single, sparse matrix. Because this handles only one matrix,
#' sensitivities are inherently deterministic and based on the dominant eigen
#' value as the best metric of the population growth rate.
#' 
#' @name sensitivity3.dgCMatrix
#' 
#' @param mats An object of class \code{dgCMatrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a single deterministic sensitivity matrix.
#' 
#' @section Notes:
#' All sensitivity matrix outputs from this function are in standard matrix
#' format.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.list}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' 
#' @examples
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
#' sensitivity3(ehrlen3$A[[1]])
#' 
#' @export
sensitivity3.dgCMatrix <- function(mats, sparse = "auto", ...)
{
  
  possible_args <- c("mats", "sparse")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- sum(diff(mats@p))
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  wcorr <- .sens3matrix_spinp(mats)
  if (sparsemethod == 1) {
    wcorr <- as(wcorr, "dgCMatrix")
  }
  
  return(wcorr)
}

#' Estimate Sensitivity of Population Growth Rate of a List of Matrices
#' 
#' \code{sensitivity3.list()} returns the sensitivities of population growth
#' rate to elements of matrices supplied in a list. The sensitivity analysis can
#' be deterministic or stochastic, but if the latter then at least two A
#' matrices must be included in the list. This function can handle large and
#' sparse matrices, and so can be used with large historical matrices, IPMs,
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @name sensitivity3.list
#' 
#' @param mats An object of class \code{matrix}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) sensitivity analysis. Defaults to
#' FALSE.
#' @param times The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param historical A logical value indicating whether matrices are historical.
#' Defaults to \code{FALSE}.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param append_mats A logical value indicating whether to include the original
#' matrices input as object \code{mats} in the output \code{lefkoSense} object.
#' Defaults to FALSE.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoSens}, which is a
#' list of 8 elements. The first, \code{h_sensmats}, is a list of historical
#' sensitivity matrices (\code{NULL} if an ahMPM is used as input). The second,
#' \code{ah_elasmats}, is a list of ahistorical sensitivity matrices if an ahMPM
#' is used as input (\code{NULL} if an hMPM is used as input). The third
#' element, \code{hstages}, the fourth element, \code{agestages}, and the fifth
#' element, \code{ahstages}, are \code{NULL}. The last 3 elements include the
#' original A matrices supplied (as the \code{A} element), followed by
#' \code{NULL}s for the U and F elements.
#' 
#' @section Notes:
#' All sensitivity matrix outputs from this function are in standard matrix
#' format.
#' 
#' Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
#' Matrix Population Models). Stochastic sensitivities are estimated as eqn.
#' 14.97 in Caswell (2001). Note that stochastic sensitivities are with regard
#' to the log of the stochastic \eqn{\lambda}.
#'
#' Currently, this function does not estimate equivalent ahistorical stochastic
#' sensitivities for input historical matrices, due to the lack of guidance
#' input on the order of stages (guidance is provided within \code{lefkoMat}
#' objects).
#'
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' The \code{time_weights} and \code{steps} arguments are now deprecated.
#' Instead, please use the \code{tweights} and \code{times} arguments.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' @seealso \code{\link{sensitivity3.dgCMatrix}()}
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
#' sensitivity3(ehrlen3$A)
#' 
#' # Cypripedium example
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' sensitivity3(cypmatrix2r$A)
#' 
#' @export
sensitivity3.list <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, historical = FALSE, seed = NA, sparse = "auto",
  append_mats = FALSE, ...) {
  
  sparsemethod <- 0
  sparse_input <- FALSE
  
  possible_args <- c("mats", "stochastic", "times", "tweights", "historical",
    "seed", "force_sparse", "append_mats", "time_weights", "steps", "sparse")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  if ("time_weights" %in% further_args_names) {
    stop("Argument time_weights is deprecated. Please use argument tweights instead",
      call. = FALSE)
  }
  if ("steps" %in% further_args_names) {
    stop("Argument steps is deprecated. Please use argument times instead",
      call. = FALSE)
  }
  
  if (is(mats[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats[[1]])
    dense_elements <- if (!sparse_input ) {
      length(which(mats[[1]] != 0))
    } else {
      sum(diff(mats[[1]]@p))
    }
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if(length(setdiff(unlist(lapply(mats, class)), c("matrix", "array", "dgCMatrix"))) > 0) {
    stop("Input list must be composed only of numeric matrices.", call. = FALSE)
  }
  if (stochastic & length(mats) < 2) {
    stop("Stochastic sensitivity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic sensitivity analysis
    message("Running deterministic analysis...")
    
    baldrick <- if (!sparse_input) {
      lapply(mats, .sens3matrix, sparsemethod)
    } else {
      lapply(mats, .sens3matrix_spinp)
    }
    
    if (historical) {
      output <- list(h_sensmats = baldrick, ah_sensmats = NULL, hstages = NULL,
        agestages = NULL, ahstages = NULL, U = NULL, F = NULL)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = baldrick, hstages = NULL,
        agestages = NULL, ahstages = NULL, U = NULL, F = NULL)
    }
    
  } else {
    # Stochastic sensitivity analysis
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    if(!any(is.na(tweights))) {
      message("Running stochastic analysis with tweights option...")
      
      returned_cube <- .stoch_senselas(mats, times = times, historical = historical,
        style = 1, sparsemethod, lefkoProj = FALSE, tweights = tweights)[[1]]
    } else {
      message("Running stochastic analysis...")
      
      returned_cube <- .stoch_senselas(mats, times = times, historical = historical,
        style = 1, sparsemethod, lefkoProj = FALSE)[[1]]
    }
    
    if (historical) {
      output <- list(h_sensmats = list(returned_cube[[1]]), ah_sensmats = NULL,
        hstages = NULL, agestages = NULL, ahstages = NULL, U = NULL, F = NULL)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = list(returned_cube[[1]]),
        hstages = NULL, agestages = NULL, ahstages = NULL, U = NULL, F = NULL)
    }
  }
  
  if (append_mats) {
    output$A <- mats
  }
  
  class(output) <- "lefkoSens"
  
  return(output)
}

#' Estimate Elasticity of Population Growth Rate to Matrix Elements
#' 
#' \code{elasticity3()} is a generic function that returns the elasticity of
#' the population growth rate to the elements of the matrices in a matrix
#' population model. Currently, this function estimates both deterministic and
#' stochastic elasticities, where the growth rate is \eqn{\lambda} in the former
#' case and the log of the stochastic \eqn{\lambda} in the latter case. This
#' function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects, as lists of matrices, and as individual matrices.
#' 
#' @name elasticity3
#' 
#' @param mats A lefkoMat object, a population projection matrix, or a list of
#' population projection matrices for which the stable stage distribution is
#' desired.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' @seealso \code{\link{elasticity3.dgCMatrix}()}
#' @seealso \code{\link{elasticity3.list}()}
#' @seealso \code{\link{summary.lefkoElas}()}
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
#' ehrlen3mean <- lmean(ehrlen3)
#' elasticity3(ehrlen3mean)
#' 
#' # Cypripedium example
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' elasticity3(cypmatrix2r)
#' 
#' @export
elasticity3 <- function(mats, ...) UseMethod("elasticity3")

#' Estimate Elasticity of Population Growth Rate of a lefkoMat Object
#' 
#' \code{elasticity3.lefkoMat()} returns the elasticities of population growth
#' rate to elements of all \code{$A} matrices in an object of class
#' \code{lefkoMat}. If deterministic, then \eqn{\lambda} is taken as the
#' population growth rate. If stochastic, then stochastic \eqn{\lambda}, or
#' the stochastic growth rate, is taken as the population growth rate. This
#' function can handle large and sparse matrices, and so can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @name elasticity3.lefkoMat
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) elasticity analysis. Defaults to
#' FALSE.
#' @param times The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param append_mats A logical value indicating whether to include the original
#' A, U, and F matrices in the output \code{lefkoElas} object.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoElas}, which is a
#' list with 8 elements. The first, \code{h_elasmats}, is a list of historical
#' elasticity matrices (\code{NULL} if an ahMPM is used as input). The second,
#' \code{ah_elasmats}, is a list of either ahistorical elasticity matrices if an
#' ahMPM is used as input, or, if an hMPM is used as input, then the result is a
#' list of elasticity matrices in which historical elasticities have been summed
#' by the stage in occasions \emph{t} and \emph{t}+1 to produce
#' historically-corrected elasticity matrices, which are equivalent in dimension
#' to ahistorical elasticity matrices but reflect the effects of stage in
#' occasion \emph{t}-1. The third element, \code{hstages}, is a data frame
#' showing historical stage pairs (NULL if ahMPM used as input). The fourth
#' element, \code{agestages}, shows age-stage combinations in the order used in
#' age-by-stage MPMs, if suppled. The fifth element, \code{ahstages}, is a data
#' frame showing the order of ahistorical stages. The last 3 elements are the A,
#' U, and F portions of the input.
#' 
#' @section Notes:
#' Deterministic elasticities are estimated as eqn. 9.72 in Caswell (2001,
#' Matrix Population Models). Stochastic elasticities are estimated as eqn.
#' 14.99 in Caswell (2001). Note that stochastic elasticities are of the
#' stochastic \eqn{\lambda}, while stochastic sensitivities are with regard to
#' the log of the stochastic \eqn{\lambda}.
#' 
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' The \code{time_weights}, \code{steps}, and \code{force_sparse} arguments are
#' now deprecated. Instead, please use the \code{tweights}, \code{times}, and
#' \code{sparse} arguments.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.dgCMatrix}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' @seealso \code{\link{elasticity3.list}()}
#' @seealso \code{\link{summary.lefkoElas}()}
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
#' elasticity3(ehrlen3, stochastic = TRUE)
#' 
#' # Cypripedium example
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' elasticity3(cypmatrix2r)
#' 
#' @export
elasticity3.lefkoMat <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, sparse = "auto", append_mats = FALSE, ...) {
  
  sparsemethod <- 0
  sparse_input <- FALSE
  
  possible_args <- c("mats", "stochastic", "times", "tweights", "historical",
    "seed", "force_sparse", "sparse", "append_mats", "time_weights", "steps")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  if ("time_weights" %in% further_args_names) {
    stop("Argument time_weights is deprecated. Please use argument tweights instead",
      call. = FALSE)
  }
  if ("steps" %in% further_args_names) {
    stop("Argument steps is deprecated. Please use argument times instead",
      call. = FALSE)
  }
  if ("force_sparse" %in% further_args_names) {
    stop("Argument force_sparse is deprecated. Please use argument sparse instead",
      call. = FALSE)
  }
  
  if (is(mats$A[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- if (!sparse_input ) {
      length(which(mats$A[[1]] != 0))
    } else {
      sum(diff(mats$A[[1]]@p))
    }
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!is.element("agestages", names(mats))) {
    mats$agestages <- NA
  }
  
  if (!stochastic) {
    # Deterministic elasticity analysis
    message("Running deterministic analysis...")
    
    baldrick <- if (is.matrix(mats$A)) {
      if (!sparse_input) {
        .elas3matrix(mats$A, sparsemethod)
      } else {
        .elas3sp_matrix(mats$A)
      }
    } else if (is.list(mats$A)) {
      
      if (all(is.na(mats$hstages))) {
        if(!sparse_input) {
          lapply(mats$A, .elas3matrix, sparsemethod)
        } else {
          lapply(mats$A, .elas3sp_matrix)
        }
      } else {
        if (!sparse_input) {
          lapply(mats$A, .elas3hlefko, mats$ahstages, mats$hstages)
        } else {
          lapply(mats$A, .elas3sp_hlefko, mats$ahstages, mats$hstages)
        }
      }
    } else {
      stop("Input not recognized.")
    }
    
    if (is.list(mats$A)) {
      multiplier <- length(mats$A)
    } else multiplier <- 1
    
    hlabels <- mats$hstages
    ahlabels <- mats$ahstages
    agelabels <- mats$agestages
    new_labels <- mats$labels
    
    if (all(is.na(mats$hstages))) {
      
      output <- list(h_elasmats = NULL, ah_elasmats = baldrick, hstages = NULL,
        agestages = mats$agestages, ahstages = ahlabels, agestages = agelabels,
        labels = new_labels)
    } else {
      
      he_list <- lapply(baldrick, function(X) {X$h_emat})
      ahe_list <- lapply(baldrick, function(X) {X$ah_emat})
      
      hlabels <- mats$hstages
      ahlabels <- mats$ahstages #Originally only the first two columns
      
      output <- list(h_elasmats = he_list, ah_elasmats = ahe_list, 
        hstages = hlabels, agestages = mats$agestages, ahstages = ahlabels,
        labels = new_labels)
    }
  } else {
    # Stochastic elasticity analysis
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    if(!any(is.na(tweights))) {
      message("Running stochastic analysis with tweights option...")
      
      returned_cubes <- .stoch_senselas(mats, times = times, historical = FALSE,
        style = 2, sparsemethod, lefkoProj = TRUE, tweights = tweights) 
    } else {
      message("Running stochastic analysis...")
      
      returned_cubes <- .stoch_senselas(mats, times = times, historical = FALSE,
        style = 2, sparsemethod, lefkoProj = TRUE) 
    }
    
    old_labels <- mats$labels
    old_labels$concat <- paste(old_labels$pop, old_labels$patch)
    labs_unique <- unique(old_labels$concat)
    lab_indices <- apply(as.matrix(labs_unique), 1, function(X) {
      found_guys <- which(old_labels$concat == X)
      min_guy <- min(found_guys)
      
      return(min_guy)
    })
    new_labels <- mats$labels[lab_indices, c(1, 2)]
    
    if (!all(is.na(mats$hstages))) {
      output <- list(h_elasmats = returned_cubes[[1]], ah_elasmats = returned_cubes[[2]],
        hstages = mats$hstages, agestages = mats$agestages, 
        ahstages = mats$ahstages, labels = new_labels)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = returned_cubes[[1]],
        hstages = mats$hstages, agestages = mats$agestages, 
        ahstages = mats$ahstages, labels = new_labels)
    }
  }
  
  if (append_mats) {
    output$A <- mats$A
    output$U <- mats$U
    output$F <- mats$F
  }
  
  class(output) <- "lefkoElas"
  
  return(output)
}

#' Estimate Elasticity of Population Growth Rate of a Single Matrix
#' 
#' \code{elasticity3.matrix()} returns the elasticities of lambda to elements
#' of a single matrix.  Because this handles only one matrix, the elasticities
#' are inherently deterministic and based on the dominant eigen value as the
#' best metric of the population growth rate. This function can handle large and
#' sparse matrices, and so can be used with large historical matrices, IPMs,
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @name elasticity3.matrix
#' 
#' @param mats An object of class \code{matrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a single elasticity matrix.
#' 
#' @section Notes:
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' The \code{force_sparse} argument is now deprecated. Please use \code{sparse}
#' instead.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.list}()}
#' @seealso \code{\link{elasticity3.dgCMatrix}()}
#' @seealso \code{\link{summary.lefkoElas}()}
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
#' ehrlen3mean <- lmean(ehrlen3)
#' elasticity3(ehrlen3mean$A[[1]])
#' 
#' @export
elasticity3.matrix <- function(mats, sparse = "auto", ...)
{
  sparsemethod <- 0
  
  possible_args <- c("mats", "force_sparse", "sparse")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  if ("force_sparse" %in% further_args_names) {
    stop("Argument force_sparse is deprecated. Please use argument sparse instead",
      call. = FALSE)
  }
  
  if (is.logical(sparse[1])) {
    if (sparse[1]) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  wcorr <- .elas3matrix(mats, FALSE)
  
  if (sparsemethod == 1) wcorr <- as(wcorr, "dgCMatrix")  
  
  return(wcorr)
}

#' Estimate Elasticity of Population Growth Rate of a Single Sparse Matrix
#' 
#' \code{elasticity3.dgCMatrix()} returns the elasticities of lambda to elements
#' of a single matrix.  Because this handles only one matrix, the elasticities
#' are inherently deterministic and based on the dominant eigen value as the
#' best metric of the population growth rate.
#' 
#' @name elasticity3.dgCMatrix
#' 
#' @param mats An object of class \code{dgCMatrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param ... Other parameters.
#' 
#' @return This function returns a single elasticity matrix in \code{dgCMatrix}
#' format.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.list}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' @seealso \code{\link{summary.lefkoElas}()}
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
#' ehrlen3mean <- lmean(ehrlen3)
#' elasticity3(ehrlen3mean$A[[1]])
#' 
#' @export
elasticity3.dgCMatrix <- function(mats, sparse = "auto", ...)
{
  
  possible_args <- c("mats", "sparse")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- sum(diff(mats@p))
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  wcorr <- .elas3sp_matrix(mats)
  
  if (sparsemethod == 0) wcorr <- as.matrix(wcorr)
  
  return(wcorr)
}

#' Estimate Elasticity of Population Growth Rate of a List of Matrices
#' 
#' \code{elasticity3.list()} returns the elasticities of lambda to elements
#' of a single matrix. This function can handle large and sparse matrices, and 
#' so can be used with large historical matrices, IPMs, age x stage matrices,
#' as well as smaller ahistorical matrices.
#' 
#' @name elasticity3.list
#' 
#' @param mats A list of objects of class \code{matrix} or \code{dgCMatrix}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) elasticity analysis. Defaults to
#' FALSE.
#' @param times The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices.
#' @param historical A logical value denoting whether the input matrices are
#' historical. Defaults to FALSE.
#' @param seed A number to use as a random number seed in stochastic projection.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param append_mats A logical value indicating whether to include the original
#' matrices input as object \code{mats} in the output \code{lefkoElas} object.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoElas}, which is a
#' list with 8 elements. The first, \code{h_elasmats}, is a list of historical
#' elasticity matrices, though in the standard list case it returns a NULL
#' value. The second, \code{ah_elasmats}, is a list of ahistorical elasticity
#' matrices. The third element, \code{hstages}, the fourth element,
#' \code{agestages}, and the fifth element, \code{ahstages}, are set to NULL.
#' The last 3 elements are the original A matrices in element A, followed by
#' NULL values for the U and F elements.
#' 
#' @section Notes:
#' Deterministic elasticities are estimated as eqn. 9.72 in Caswell (2001,
#' Matrix Population Models). Stochastic elasticities are estimated as eqn.
#' 14.99 in Caswell (2001). Note that stochastic elasticities are of stochastic
#' \eqn{\lambda}, while stochastic sensitivities are with regard to the log of
#' the stochastic \eqn{\lambda}.
#' 
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 30 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' The \code{time_weights}, \code{steps}, and \code{force_sparse} arguments are
#' now deprecated. Instead, please use the \code{tweights}, \code{times}, and
#' \code{sparse} arguments.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' @seealso \code{\link{elasticity3.dgCMatrix}()}
#' @seealso \code{\link{summary.lefkoElas}()}
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
#' elasticity3(ehrlen3$A, stochastic = TRUE)
#' 
#' # Cypripedium example
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' elasticity3(cypmatrix2r$A)
#' 
#' @export
elasticity3.list <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, historical = FALSE, seed = NA, sparse = "auto",
  append_mats = FALSE, ...) {
  
  sparsemethod <- 0
  sparse_input <- FALSE
  
  possible_args <- c("mats", "stochastic", "times", "tweights", "historical",
    "seed", "force_sparse", "append_mats", "time_weights", "steps")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  if ("time_weights" %in% further_args_names) {
    stop("Argument time_weights is deprecated. Please use argument tweights instead",
      call. = FALSE)
  }
  if ("steps" %in% further_args_names) {
    stop("Argument steps is deprecated. Please use argument times instead",
      call. = FALSE)
  }
  if ("force_sparse" %in% further_args_names) {
    stop("Argument force_sparse is deprecated. Please use argument sparse instead",
      call. = FALSE)
  }
  
  if (is(mats[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats[[1]])
    dense_elements <- if (!sparse_input ) {
      length(which(mats[[1]] != 0))
    } else {
      sum(diff(mats[[1]]@p))
    }
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if(length(setdiff(unlist(lapply(mats, class)), c("matrix", "array", "dgCMatrix"))) > 0) {
    stop("Input list must be composed only of numeric matrices.", call. = FALSE)
  }
  if (stochastic & length(mats) < 2) {
    stop("Stochastic elasticity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic elasticity analysis
    message("Running deterministic analysis...")
    
    baldrick <- if (!sparse_input) {
      lapply(mats, .elas3matrix, sparsemethod)
    } else {
      lapply(mats, .elas3sp_matrix)
    }
    
    if (historical) {
      output <- list(h_elasmats = baldrick, ah_elasmats = NULL, hstages = NULL,
        ahstages = NULL)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = baldrick, hstages = NULL,
        ahstages = NULL)
    }
  } else {
    # Stochastic elasticity analysis
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    if(!any(is.na(tweights))) {
      message("Running stochastic analysis with tweights option...")
      
      returned_cube <- .stoch_senselas(mats, times = times, historical = historical,
        style = 2, sparsemethod, lefkoProj = FALSE, tweights = tweights)[[1]]
    } else {
      message("Running stochastic analysis...")
      
      returned_cube <- .stoch_senselas(mats, times = times, historical = historical,
        style = 2, sparsemethod, lefkoProj = FALSE)[[1]]
    }
    
    if (historical) {
      output <- list(h_elasmats = list(returned_cube[[1]]), ah_elasmats = NULL,
        hstages = NULL, ahstages = NULL)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = list(returned_cube[[1]]),
        hstages = NULL, ahstages = NULL)
    }
  }
  
  if (append_mats) {
    output$A <- mats
  }
  
  class(output) <- "lefkoElas"
  
  return(output)
}

#' Conduct a Life Table Response Experiment
#' 
#' \code{ltre3()} returns a set of matrices of one-way LTRE (life table response
#' experiment), stochastic LTRE (sLTRE) matrices, or small noise approximation
#' LTRE (sna-LTRE) contributions.
#' 
#' @name ltre3
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param refmats A reference lefkoMat object, or matrix, for use as the
#' control. Default is \code{NA}, which sets to the same object as \code{mats}.
#' @param ref A numeric value indicating which matrix or matrices in
#' \code{refmats} to use as the control. The numbers used must correspond to the
#' number of the matrices in the \code{labels} element of the associated
#' \code{lefkoMat} object. The default setting, \code{NA}, uses all entries in
#' \code{refmats}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (\code{FALSE}) or stochastic (\code{TRUE}) elasticity analysis.
#' Defaults to \code{FALSE}.
#' @param times The number of occasions to project forward in stochastic
#' simulation. Defaults to \code{10000}.
#' @param burnin The number of initial steps to ignore in stochastic projection
#' when calculating stochastic elasticities. Must be smaller than \code{steps}.
#' Defaults to \code{3000}.
#' @param tweights An optional numeric vector or matrix denoting the
#' probabilities of choosing each matrix in a stochastic projection. If a matrix
#' is input, then a first-order Markovian environment is assumed, in which the
#' probability of choosing a specific annual matrix depends on which annual
#' matrix is currently chosen. If a vector is input, then the choice of annual
#' matrix is assumed to be independent of the current matrix. Defaults to equal
#' weighting among matrices. Note that SNA-LTRE analysis cannot take matrix
#' input.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}, in which case sparse matrix encoding is used with square
#' matrices with at least 50 rows and no more than 50\% of elements with values
#' greater than zero.
#' @param seed Optional numeric value corresponding to the random seed for
#' stochastic simulation.
#' @param append_mats A logical value denoting whether to include the original
#' \code{A}, \code{U}, and \code{F} matrices in the returned \code{lefkoLTRE}
#' object. Defaults to \code{FALSE}.
#' @param sna_ltre A logical value indicating whether to treat stochastic LTRE
#' via the sna-LTRE approach from Davison et al. (2019) (\code{TRUE}), or the
#' stochastic LTRE approximation from Davison et al. (2010) (\code{FALSE}).
#' Defaults to \code{FALSE}.
#' @param tol A numeric value indicating a lower positive limit to matrix
#' element values when applied to stochastic and small noise approximation LTRE
#' estimation protocols. Matrix element values lower than this will be treated
#' as \code{0.0} values. Defaults to \code{1e-30}.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoLTRE}. This
#' includes a list of LTRE matrices as object \code{cont_mean} if a
#' deterministic LTRE is called for, or a list of mean-value LTRE matrices as
#' object \code{cont_mean} and a list of SD-value LTRE matrices as object
#' \code{cont_sd} if a stochastic LTRE is called for. If a small-noise
#' approximation LTRE (SNA-LTRE) is performed, then the output includes six
#' objects: \code{cont_mean}, which provides the contributions of shifts in mean
#' matrix elements; \code{cont_elas}, which provides the contributions of shifts
#' in the elasticities of matrix elements; \code{cont_cv}, which provides the
#' contributions of temporal variation in matrix elements; \code{cont_corr},
#' which provides the contributions of temporal correlations in matrix elements;
#' \code{r_values_m}, which provides a vector of log deterministic lambda values
#' for treatment populations; and \code{r_values_ref}, which provides the log
#' deterministic lambda of the mean reference matrix.This is followed by the
#' stageframe as object \code{ahstages}, the order of historical stages as
#' object \code{hstages}, the age-by-stage order as object \code{agestages}, the
#' order of matrices as object \code{labels}, and, if requested, the original A,
#' U, and F matrices.
#' 
#' @section Notes:
#' Deterministic LTRE is one-way, fixed, and based on the sensitivities of the
#' matrix midway between each input matrix and the reference matrix, per Caswell
#' (2001, Matrix Population Models, Sinauer Associates, MA, USA). Stochastic
#' LTRE is performed via two methods. The stochastic LTRE approximation is
#' simulated per Davison et al. (2010) Journal of Ecology 98:255-267
#' (doi: 10.1111/j.1365-2745.2009.01611.x). The small noise approximation
#' (sna-LTRE) is analyzed per Davison et al. (2019) Ecological Modelling 408:
#' 108760 (doi: 10.1016/j.ecolmodel.2019.108760).
#' 
#' All stochastic and small noise approximation LTREs conducted without
#' reference matrices are performed as spatial tests of the population dynamics
#' among patches.
#' 
#' Default behavior for stochastic LTRE uses the full population provided in
#' \code{mats} as the reference if no \code{refmats} and \code{ref} is provided.
#' If no \code{refmats} is provided but \code{ref} is, then the matrices noted
#' in \code{ref} are used as the reference matrix set. Year and patch order is
#' utilized from object \code{mats}, but not from object \code{refmats}, in
#' which each matrix is assumed to represent a different year from one
#' population. This function cannot currently handle multiple populations within
#' the same \code{mats} object (although such analysis is possible if these
#' populations are designated as patches instead).
#' 
#' If \code{sparse = "auto"}, the default, then sparse matrix encoding
#' will be used if the size of the input matrices is at least 50 columns by 50
#' rows for deterministic and stochastic LTREs and 10 columns by 10 rows for
#' small noise approximation LTREs, in all cases as long as 50\% of the elements
#' in the first matrix are non-zero.
#' 
#' Stochastic LTREs do not test for the impact of temporal change in vital
#' rates. An MPM with a single population, a single patch, and only annual
#' matrices will produce contributions of 0 to stochastic \eqn{\lambda}.
#' 
#' Speed can sometimes be increased by shifting from automatic sparse matrix
#' determination to forced dense or sparse matrix projection. This will most
#' likely occur when matrices have between 10 and 300 rows and columns.
#' Defaults work best when matrices are very small and dense, or very large and
#' sparse.
#' 
#' SNA-LTRE analysis cannot test the impact of first-order Markovian
#' environments. However, different random weightings of annual matrices are
#' allowed if given in vector format.
#' 
#' The \code{time_weights}, \code{steps}, \code{force_sparse}, and \code{rseed}
#' arguments are now deprecated. Instead, please use the \code{tweights},
#' \code{times}, \code{sparse}, and \code{seed} arguments.
#' 
#' @seealso \code{\link{summary.lefkoLTRE}()}
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' ltre3(cypmatrix2r, sna_ltre = TRUE)
#' 
#' @export
ltre3 <- function(mats, refmats = NA, ref = NA, stochastic = FALSE,
  times = 10000, burnin = 3000, tweights = NA, sparse = "auto", seed = NA,
  append_mats = FALSE, sna_ltre = FALSE, tol = 1e-30, ...) {
  
  sparsemethod <- 0
  sparse_input <- FALSE
  
  possible_args <- c("mats", "stochastic", "times", "tweights", "rseed",
    "seed", "force_sparse", "sparse", "append_mats", "time_weights", "steps",
    "refmats", "ref", "sna_ltre", "tol")
  further_args <- list(...)
  further_args_names <- names(further_args)
  if (length(setdiff(further_args_names, possible_args)) > 0) {
    stop("Some arguments not recognized.", call. = FALSE)
  }
  if ("time_weights" %in% further_args_names) {
    stop("Argument time_weights is deprecated. Please use argument tweights instead",
      call. = FALSE)
  }
  if ("steps" %in% further_args_names) {
    stop("Argument steps is deprecated. Please use argument times instead",
      call. = FALSE)
  }
  if ("force_sparse" %in% further_args_names) {
    stop("Argument force_sparse is deprecated. Please use argument sparse instead",
      call. = FALSE)
  }
  if ("rseed" %in% further_args_names) {
    stop("Argument rseed is deprecated. Please use argument seed instead",
      call. = FALSE)
  }
  
  if (!is(mats, "lefkoMat")) stop("Function ltre3() requires a lefkoMat object as input.",
    call. = FALSE)
  if (is(mats$A[[1]], "dgCMatrix")) sparse_input = TRUE
  
  if (is.logical(sparse)) {
    if (sparse) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  } else if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true", "ja", "tak"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false", "nein", "nie"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- if (!sparse_input ) {
      length(which(mats$A[[1]] != 0))
    } else {
      sum(diff(mats$A[[1]]@p))
    }
    
    if ((dense_elements / elements_total) < 0.5 & elements_total > 2500) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!is.element("agestages", names(mats))) {
    mats$agestages <- NULL
  }
  
  if (all(is.na(refmats))) {
    warning("Matrices input as mats will also be used in reference matrix calculation.",
      call. = FALSE)
  } else if (is(refmats, "lefkoMat")) {
    refmats <- refmats$A
  } else if (is.list(refmats)) {
    if (!is.matrix(refmats[[1]]) & !is(refmats[[1]], "dgCMatrix")) {
      stop("Object refmats not recognized. Use only objects of
        class lefkoMat, list, matrix, or dgCMatrix.", call. = FALSE)
    }
  } else if (is.matrix(refmats) | is(refmats, "dgCMatrix")) {
    refmats <- list(refmats)
  } else {
    stop("Object refmats not recognized. Use only objects of class lefkoMat, list, or matrix.",
      call. = FALSE)
  }
  
  if (all(is.numeric(ref))) {
    if (length(ref) > 1) {
      if (!stochastic) {
        message("This function currently conducts deterministic LTREs against only single
          matrices, or mean matrices. Will use the mean.")
      }
    }
  } else if (all(is.na(ref))) {
    message("Using all refmats matrices in reference matrix calculation.")
    
    if (!all(is.na(refmats))) {
      ref <- c(1:length(refmats))
    } else {
      ref <- c(1:length(mats$A))
    }
  }
  
  if (sna_ltre & !stochastic) {
    stochastic <- TRUE
  }
  
  if (!stochastic & length(ref) > 1) {
    meanout <- 1
  } else if (stochastic & length(ref) == 1) {
    stop("Stochastic LTRE requires multiple matrices in object refmats.", call. = FALSE)
  } else if (stochastic & length(mats$A) == 1) {
    stop("Stochastic LTRE requires multiple matrices in object mats.", call. = FALSE)
  } else {
    meanout <- 0;
  }
  
  if (all(is.na(refmats))) {
    if (any(ref > length(mats$A)) | any(ref < 1)) {
      stop("Numbers of matrices provided in object ref must be integers between 1 and
        the total number of matrices in object mats (or refmats, if given).",
        call. = FALSE)
    }
    
  } else {
    if (dim(mats$A[[1]])[1] != dim(refmats[[1]])[1] | dim(mats$A[[1]])[2] != dim(refmats[[1]])[2]) {
      stop("Matrices used in objects mats and refmats must be of equal dimension.",
        call. = FALSE)
    }
  }
  
  if (!stochastic) {
    # Deterministic LTRE analysis
    
    if (all(is.na(refmats))) {
      baldrick <- .ltre3matrix(mats$A, refnum = ref, mean = meanout,
        sparse = sparsemethod)
    } else {
      baldrick <- .ltre3matrix(mats$A, refnum = ref, refmats_ = refmats,
        mean = meanout, sparse = sparsemethod)
    }
    
    baldrick$ahstages <- mats$ahstages
    baldrick$agestages <- mats$agestages
    baldrick$hstages <- mats$hstages
    baldrick$labels <- mats$labels
    
  } else {
    # Stochastic and SNA LTRE analysis
    if (burnin >= times) {
      stop("Option burnin must be smaller than option steps.", call. = FALSE)
    }
    
    if (!is.na(seed[1])) {
      set.seed(seed[1])
    }
    
    if (!sna_ltre) {
      if (all(is.na(refmats))) {
        if (all(is.na(tweights))) {
          baldrick <- .sltre3matrix(mats$A, labels = mats$labels, refnum = ref,
            steps = times, burnin = burnin, sparse = sparsemethod,
            tol_used = tol)
        } else {
          baldrick <- .sltre3matrix(mats$A, labels = mats$labels, refnum = ref,
            tweights_ = tweights, steps = times, burnin = burnin,
            sparse = sparsemethod, tol_used = tol)
        }
      } else {
        if (all(is.na(tweights))) {
          baldrick <- .sltre3matrix(mats$A, labels = mats$labels, refnum = ref,
            refmats_ = refmats, steps = times, burnin = burnin,
            sparse = sparsemethod, tol_used = tol)
        } else {
          baldrick <- .sltre3matrix(mats$A, labels = mats$labels, refnum = ref,
            refmats_ = refmats, tweights_ = tweights, steps = times,
            burnin = burnin, sparse = sparsemethod, tol_used = tol)
        }
      }
    } else {
      if (is.matrix(tweights)) {
        stop("SNA-LTRE analysis cannot proceed with matrix input in the tweights argument. Use vectors only.",
          call. = FALSE)
      }
      
      if (all(is.na(refmats))) {
        if (all(is.na(tweights))) {
          baldrick <- .snaltre3matrix(mats$A, labels = mats$labels,
            refnum = ref, sparse = sparsemethod, tol_used = tol)
        } else {
          baldrick <- .snaltre3matrix(mats$A, labels = mats$labels,
            refnum = ref, tweights_ = tweights, sparse = sparsemethod,
            tol_used = tol)
        }
      } else {
        if (all(is.na(tweights))) {
          baldrick <- .snaltre3matrix(mats$A, labels = mats$labels,
            refnum = ref, refmats_ = refmats, sparse = sparsemethod,
            tol_used = tol)
        } else {
          baldrick <- .snaltre3matrix(mats$A, labels = mats$labels,
            refnum = ref, refmats_ = refmats, tweights_ = tweights,
            sparse = sparsemethod, tol_used = tol)
        }
      }
    }
    
    alabels <- mats$labels[,c("pop", "patch")]
    alabels <- unique(alabels)
    
    baldrick$ahstages <- mats$ahstages
    baldrick$agestages <- mats$agestages
    baldrick$hstages <- mats$hstages
    baldrick$labels <- alabels
  }
  
  if (append_mats) {
    baldrick$A <- mats$A
    baldrick$U <- mats$U
    baldrick$F <- mats$F
  }
  
  return(baldrick)
}

#' Summarize lefkoElas Objects
#' 
#' Function \code{summary.lefkoElas()} summarizes \code{lefkoElas} objects.
#' Particularly, it breaks down elasticity values by the kind of ahistorical
#' and, if applicable, historical transition.
#' 
#' @name summary.lefkoElas
#' 
#' @param object A \code{lefkoElas} object.
#' @param ... Other parameters currently not utilized.
#' 
#' @return A list composed of 2 data frames. The first, \code{hist}, is a data
#' frame showing the summed elasticities for all 16 kinds of historical
#' transition per matrix, with each column corresponding to each elasticity
#' matrix in order. The second, \code{ahist}, is a data frame showing the
#' summed elasticities for all 4 kinds of ahistorical transition per matrix,
#' with each column corresponding to each elasticity matrix in order.
#' 
#' @examples
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
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframe, historical = FALSE)
#'   
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen2 <- rlefko2(data = lathvert, stageframe = lathframe, year = "all",
#'   stages = c("stage3", "stage2"), supplement = lathsupp2,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3elas <- elasticity3(ehrlen3)
#' ehrlen2elas <- elasticity3(ehrlen2)
#' 
#' summary(ehrlen3elas)
#' summary(ehrlen2elas)
#' 
#' @export
summary.lefkoElas <- function(object, ...) {
  
  sparse_input <- FALSE
  num_h_mats <- length(object$h_elasmats)
  num_ah_mats <- length(object$ah_elasmats)
  
  used_emats <- if(num_h_mats == 0) object$ah_elasmats else object$h_elasmats
  
  used_iterations <- if(num_h_mats > 0) num_h_mats else num_ah_mats
  if (is(used_emats[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if (num_h_mats == 0) {
    if (!all(is.null(object$agestages))) {
      if (!all(is.na(object$agestages))) {
        if (is.element("stage_id", names(object$agestages))) {
          new_ahstages_list <- apply(as.matrix(c(1:length(object$agestages$stage_id))), 
            1, function(X) {
              return(object$ahstages[which(object$ahstages$stage_id == object$agestages$stage_id[X]),])
            })
          new_ahstages <- do.call("rbind.data.frame", new_ahstages_list)
          indices <- .bambi2(new_ahstages)
        } else {
          indices <- .bambi2(object$ahstages)
        }
      } else {
        indices <- .bambi2(object$ahstages)
      }
    } else {
      if (all(is.null(object$ahstages))) {
        mat_size <- dim(object$ah_elasmats[[1]])[1]
        new_sf <- sf_skeleton(mat_size, standard = FALSE)
        
        indices <- .bambi2(new_sf)
      } else {
        indices <- .bambi2(object$ahstages)
      }
    }
    
  } else {
    indices <- .bambi3(object$ahstages, object$hstages)
  }
  
  for (i in c(1:used_iterations)) {
    if (!sparse_input) {
      trialguy <- .demolition3(used_emats[[i]], indices)
    } else {
      trialguy <- .demolition3sp(used_emats[[i]], indices)
    }
    
    if (i == 1) {
      if (num_h_mats > 0) hist <- trialguy$hist[,c(1,2)]
      ahist <- trialguy$ahist[,c(1,2)]
      
      if (num_h_mats > 0) names(hist)[2] <- "matrix1"
      names(ahist)[2] <- "matrix1"
    } else {
      if (num_h_mats > 0) hist <- cbind.data.frame(hist, trialguy$hist[,2])
      ahist <- cbind.data.frame(ahist, trialguy$ahist[,2])
      if (num_h_mats > 0) names(hist)[(i+1)] <- paste0("matrix", i)
      names(ahist)[(i+1)] <- paste0("matrix", i)
    }
  }
  
  if (num_h_mats > 0) {
    output <- list(hist = hist, ahist = ahist)
  } else {
    output <- list(hist = NULL, ahist = ahist)
  }
  
  return (output)
}

#' Summarize lefkoLTRE Objects
#' 
#' Function \code{summary.lefkoLTRE()} summarizes \code{lefkoLTRE} objects.
#' Particularly, it breaks down LTRE contributions by the kind of ahistorical
#' and, if applicable, historical transition.
#' 
#' @name summary.lefkoLTRE
#' 
#' @param object A \code{lefkoLTRE} object.
#' @param ... Other parameters currently not utilized.
#' 
#' @return A list of data frames. In all cases, the first data frame is one
#' showing the positive, negative, and total contributions of elements in
#' each LTRE contribution matrix. If not a SNA-LTRE, then there are an
#' additional two (if deterministic) or four (if stochastic) data frames. If
#' deterministic, then \code{hist_det} is a data frame showing the summed LTRE
#' contributions for all 16 kinds of historical transition per matrix, with each
#' column corresponding to each A matrix in order, followed by all summed
#' positive and all summed negative contributions. Object \code{ahist_det} is a
#' data frame showing the summed LTRE contributions for all four kinds of
#' ahistorical transition per matrix, with order as before, followed by summed
#' positive and summed negative contributions. If stochastic, then
#' \code{hist_mean} and \code{hist_sd} are the summed LTRE contributions for the
#' mean vital rates and variability in vital rates, respectively, according to
#' all 16 historical transition types, followed by summed positive and negative
#' contributions, and \code{ahist_mean} and \code{ahist_sd} are the equivalent
#' ahistorical versions. The output for the SNA-LTRE also includes the
#' logs of the deterministic lambda estimated through function \code{ltre3()}.
#' 
#' @examples
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
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframe, historical = FALSE)
#'   
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen2 <- rlefko2(data = lathvert, stageframe = lathframe, year = "all",
#'   stages = c("stage3", "stage2"), supplement = lathsupp2,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3ltre <- ltre3(ehrlen3)
#' summary(ehrlen3ltre)
#' 
#' @export
summary.lefkoLTRE <- function(object, ...) {
  
  trialguy1 <- trialguy2 <- NULL
  hist1 <- hist2 <- ahist1 <- ahist2 <- NULL
  sparse_input <- FALSE
  
  if (is.element("cont_cv", names(object))) {
    ltretype <- 3 # sna-LTRE
    
  } else if (is.element("cont_sd", names(object))) {
    ltretype <- 2 # Stochastic LTRE
    
  } else if (is.element("cont_mean", names(object))) {
    ltretype <- 1 # Deterministic LTRE
  } else {
    stop("Input object is of unrecognized type. Use only lefkoLTRE obects with this function.",
      call. = FALSE)
  }
  
  numstages <- dim(object$cont_mean[[1]])[1]
  if (is(object$cont_mean[[1]], "dgCMatrix")) sparse_input <- TRUE
  
  if(!all(is.na(object$hstages))) {
    if(numstages == dim(object$hstages)[1]) {
      historical <- TRUE
    } else {
      historical <- FALSE
    }
  } else {
    historical <- FALSE
  }
  
  if (!historical) {
    if (!all(is.null(object$agestages))) {
      if (!all(is.na(object$agestages))) {
        if (is.element("stage_id", names(object$agestages))) {
          new_ahstages_list <- apply(as.matrix(c(1:length(object$agestages$stage_id))), 
            1, function(X) {
              return(object$ahstages[which(object$ahstages$stage_id == object$agestages$stage_id[X]),])
            })
          new_ahstages <- do.call("rbind.data.frame", new_ahstages_list)
          indices <- .bambi2(new_ahstages)
        } else {
          indices <- .bambi2(object$ahstages)
        }
      } else {
        indices <- .bambi2(object$ahstages)
      }
    } else {
      if (all(is.null(object$ahstages))) {
        mat_size <- dim(object$ah_elasmats[[1]])[1]
        new_sf <- sf_skeleton(mat_size, standard = FALSE)
        
        indices <- .bambi2(new_sf)
      } else {
        indices <- .bambi2(object$ahstages)
      }
    }
  } else {
    indices <- .bambi3(object$ahstages, object$hstages)
  }
  
  used_iterations <- length(object$cont_mean)
  
  # General summary for all types of LTRE
  general_df <- .demolition4(object)
  
  # Additional summaries for LTRE and sLTRE
  for (i in c(1:used_iterations)) {
    if (!sparse_input) {
      trialguy1 <- .demolition3(object$cont_mean[[i]], indices)
    } else {
      trialguy1 <- .demolition3sp(object$cont_mean[[i]], indices)
    }
    
    if (ltretype == 2) {
      if (!sparse_input) {
        trialguy2 <- .demolition3(object$cont_sd[[i]], indices)
      } else {
        trialguy2 <- .demolition3sp(object$cont_sd[[i]], indices)
      }
    }
    
    if (i == 1) {
      hist1 <- trialguy1$hist
      ahist1 <- trialguy1$ahist
      if (historical) names(hist1)[2] <- "matrix1"
      if (historical) names(hist1)[3] <- "matrix1_pos"
      if (historical) names(hist1)[4] <- "matrix1_neg"
      names(ahist1)[2] <- "matrix1"
      names(ahist1)[3] <- "matrix1_pos"
      names(ahist1)[4] <- "matrix1_neg"
      
      if (ltretype == 2) {
        hist2 <- trialguy2$hist
        ahist2 <- trialguy2$ahist
        if (historical) names(hist2)[2] <- "matrix1"
        if (historical) names(hist2)[3] <- "matrix1_pos"
        if (historical) names(hist2)[4] <- "matrix1_neg"
        names(ahist2)[2] <- "matrix1"
        names(ahist2)[3] <- "matrix1_pos"
        names(ahist2)[4] <- "matrix1_neg"
      }
      
    } else {
      if (historical) hist1 <- cbind.data.frame(hist1, trialguy1$hist[,c(2:4)])
      ahist1 <- cbind.data.frame(ahist1, trialguy1$ahist[,c(2:4)])
      if (historical) names(hist1)[(((i-1)*3)+2)] <- paste0("matrix", i)
      if (historical) names(hist1)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
      if (historical) names(hist1)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
      names(ahist1)[(((i-1)*3)+2)] <- paste0("matrix", i)
      names(ahist1)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
      names(ahist1)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
      
      if (ltretype == 2) {
        if (historical) hist2 <- cbind.data.frame(hist2, trialguy2$hist[,c(2:4)])
        ahist2 <- cbind.data.frame(ahist2, trialguy2$ahist[,c(2:4)])
        if (historical) names(hist2)[(((i-1)*3)+2)] <- paste0("matrix", i)
        if (historical) names(hist2)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
        if (historical) names(hist2)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
        names(ahist2)[(((i-1)*3)+2)] <- paste0("matrix", i)
        names(ahist2)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
        names(ahist2)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
      }
    }
  }
  
  output <- if (ltretype == 1) {
    list(overall = general_df, hist_mean = hist1, ahist_mean = ahist1)
  } else if (ltretype == 2) {
    list(overall = general_df, hist_mean = hist1, hist_sd = hist2,
      ahist_mean = ahist1, ahist_sd = ahist2)
  } else {
    list(overall = general_df, hist_mean = hist1, ahist_mean = ahist1,
      r_values_m = object$r_values_m, r_value_ref = object$r_value_ref)
  }
  
  return (output)
}

#' Summarize lefkoProj Objects
#' 
#' Function \code{summary.lefkoProj()} summarizes \code{lefkoProj} objects.
#' Particularly, it breaks down the data frames provided in the 
#' \code{projection} element in ways meaningful for those running simulations.
#' 
#' @name summary.lefkoProj
#' 
#' @param object A \code{lefkoProj} object.
#' @param threshold A threshold population size to be searched for in
#' projections. Defaults to 1.
#' @param inf_alive A logical value indicating whether to treat infinitely
#' large population size as indicating that the population is still extant.
#' If \code{FALSE}, then the population is considered extinct. Defaults to
#' \code{TRUE}.
#' @param milepost A numeric vector indicating at which points in the projection
#' to assess detailed results. Can be input as integer values, in which case
#' each number must be between 1 and the total number of occasions projected in
#' each projection, or decimals between 0 and 1, which would then be translated
#' into the corresponding projection steps of the total. Defaults to
#' \code{c(0, 0.25, 0.50, 0.75, 1.00)}.
#' @param ext_time A logical value indicating whether to output extinction times
#' per population-patch. Defaults to \code{FALSE}.
#' @param ... Other parameters currently not utilized.
#' 
#' @return Apart from a statement of the results, this function outputs a list
#' with the following elements:
#' \item{milepost_sums}{A data frame showing the number of replicates at each
#' of the milepost times that is above the threshold population/patch size.}
#' \item{extinction_times}{A dataframe showing the numbers of replicates going
#' extinct (\code{ext_reps}) and mean extinction time (\code{ext_time}) per
#' population-patch. If \code{ext_time = FALSE}, then only outputs \code{NA}.}
#' 
#' @section Notes:
#' The \code{inf_alive} and \code{ext_time} options both assess whether
#' replicates have reached a value of \code{NaN} or \code{Inf}. If
#' \code{inf_alive = TRUE} or \code{ext_time = TRUE} and one of these values is
#' found, then the replicate is counted in the \code{milepost_sums} object if
#' the last numeric value in the replicate is above the \code{threshold} value,
#' and is counted as extant and not extinct if the last numeric value in the
#' replicate is above the extinction threshold of a single individual.
#' 
#' Extinction time is calculated on the basis of whether the replicate ever
#' falls below a single individual. A replicate with a positive population size
#' below 0.0 that manages to rise above 1.0 individual is still considered to
#' have gone extinct the first time it crossed below 1.0.
#' 
#' If the input \code{lefkoProj} object is a mixture of two or more other
#' \code{lefkoProj} objects, then mileposts will be given relative to the
#' maximum number of time steps noted.
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
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
#'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
#'   repmatrix = lathrepm, supplement = lathsupp3, yearcol = "year2",
#'   indivcol = "individ")
#' 
#' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
#' summary(lathproj)
#' 
#' # Cypripedium example
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
#'   supplement = cypsupp3r, yearcol = "year2", 
#'   patchcol = "patchid", indivcol = "individ")
#' 
#' cypstoch <- projection3(cypmatrix3r, nreps = 5, stochastic = TRUE)
#' summary(cypstoch, ext_time = TRUE)
#' 
#' @export
summary.lefkoProj <- function(object, threshold = 1, inf_alive = TRUE,
  milepost = c(0, 0.25, 0.50, 0.75, 1.00), ext_time = FALSE, ...) {
  
  num_reps_vec <- num_times_vec <- NULL
  appended <- FALSE
  max_times <- max_reps <- 1L
  ave_times <- ave_reps <- 1.0
  
  poppatches <- length(object$labels[,1])
  
  if (is.matrix(object$control)) {
    unique_indices <- unique(object$control[,1])
    dup_counts <- apply(as.matrix(unique_indices), 1, function(X) {
      found <- length(which(object$control[,1] == X))
      return(found)
    })
    
    appended <- TRUE
    max_reps <- max(dup_counts, na.rm = TRUE)
    ave_reps <- mean(dup_counts, na.rm = TRUE)
    
    max_times <- max(object$control[,3], na.rm = TRUE)
    total_reps_across <- sum(object$control[,2], na.rm = TRUE)
    
    times_contribs <- apply(as.matrix(c(1:(dim(object$control)[1]))), 1, function(X) {
      mean_div <- object$control[X, 2] / total_reps_across
      current_contrib <- object$control[X, 3] * mean_div
      
      return(current_contrib)
    })
    ave_times <- sum(times_contribs)
    
    num_reps_vec <- apply(as.matrix(c(1:(dim(object$labels)[1]))), 1, function(X) {
      found_reps <- sum(object$control[which(object$control[,1] == X), 2], na.rm = TRUE)
    })
    
    num_times_vec <- apply(as.matrix(c(1:(dim(object$labels)[1]))), 1, function(X) {
      found_times <- max(object$control[which(object$control[,1] == X), 3], na.rm = TRUE)
    })
  } else {
    max_reps <- object$control[1]
    max_times <- ave_times <- object$control[2]
    
    num_reps_vec <- rep(max_reps, poppatches)
    num_times_vec <- rep(max_times, poppatches)
  }
  
  if (any(milepost < 0)) {
    stop("Option milepost may not take negative values.", call. = FALSE)
  }
  if (any(milepost > max_times)) {
    stop("Option milepost may not take values higher than the number of actual
      number of projected occasions.", call. = FALSE)
  }
  
  if (inf_alive | ext_time) {
    for (i in c(1:poppatches)) {
      for (j in c(1:num_reps_vec[i])) {
        for (k in c(1:(num_times_vec[i] + 1))) {
          if ((is.nan(object$pop_size[[i]][j, k]) | is.infinite(object$pop_size[[i]][j, k])) & k > 1) {
            object$pop_size[[i]][j, k] <- object$pop_size[[i]][j, (k - 1)]   #. max_found
          }
        }
      }
    }
  }
  
  if (ext_time) {
    the_numbers <- apply(as.matrix(c(1:poppatches)), 1, function(X) {
      freemasonry <- apply(as.matrix(c(1:num_reps_vec[X])), 1, function(Y) {
        ext_points <- which(object$pop_size[[X]][Y,] < 1)
        if (length(ext_points) > 0) return (min(ext_points)) else return (NA)
      })
      ext_varmints <- length(which(!is.na(freemasonry) & !is.nan(freemasonry)))
      if (ext_varmints > 0) {
        ext_time <- mean(freemasonry, na.rm = TRUE)
      } else {
        ext_time <- NA
      }
      return (c(ext_varmints, ext_time))
    })
    
    the_numbers <- t(the_numbers)
    the_numbers <- as.data.frame(the_numbers)
    colnames(the_numbers) <- c("ext_reps", "ext_time")
    
    if (dim(object$labels)[1] > 1) {
      row_labels <- apply(object$labels, 1, function(X) {
        paste(X[1], X[2])
      })
      rownames(the_numbers) <- row_labels
    }
  } else {
    the_numbers <- NA
  }
  
  for (i in c(1:poppatches)) {
    if (any(milepost > num_times_vec[i])) {
      stop("Entered milepost values are outside the allowable range.", call. = FALSE)
    }
  }
  
  milepost_sums <- apply(as.matrix(c(1:poppatches)), 1, function (X) {
    
    used_milepost <- milepost
    
    if (all(milepost >= 0) & all(milepost <= 1)) {
      used_milepost <- floor(used_milepost * num_times_vec[X]) + 1
    } else if (any(milepost == 0)) {
      used_milepost <- used_milepost + 1
    }
    
    if (num_reps_vec[X] > 1) {
      phew <- apply(as.matrix(object$pop_size[[X]][,used_milepost]), 2, function(Y) {
        above_vector <- which(as.vector(Y) >= threshold)
        
        return(length(above_vector))
      })
      return(phew)
    } else {
      phew <- apply(as.matrix(object$pop_size[[X]][,used_milepost]), 1, function(Y) {
        above_vector <- which(as.vector(Y) >= threshold)
        
        return(length(above_vector))
      })
      return(phew)
    }
  })
  
  if (is.matrix(milepost_sums)) {
    rownames(milepost_sums) <- milepost
    
    col_labels <- apply(object$labels, 1, function(X) {
      paste(X[1], X[2])
    })
    colnames(milepost_sums) <- col_labels
  } else {
    rownames(milepost_sums) <- milepost
  }
  
  writeLines(paste0("\nThe input lefkoProj object covers ", poppatches,
    " population-patches."), con = stdout())
  if (appended) {
    writeLines(paste0("It is an appended projection, including an average and maximum of ",
      format(ave_times, digits = 5), " and ", max_times, " steps per "))
    writeLines(paste0("replicate, and an average and maximum of ",
      format(ave_reps, digits = 5), " and ", max_reps, " replicates."))
  } else {
  writeLines(paste0("It is a single projection including ", max_times,
    " projected steps per replicate, and ", max_reps, " replicates, respectively."),
    con = stdout())
  }
  writeLines(paste0("The number of replicates with population size above the threshold size of ",
    threshold, " is as in"), con = stdout())
  writeLines(paste0("the following matrix, with pop-patches given by column and milepost times given by row: \n"),
    con = stdout())
  
  output <- list(milepost_sums = milepost_sums, extinction_times = the_numbers)
  
  return (output)
}

#' Plot Projection Simulations
#' 
#' Function \code{plot.lefkoProj()} produces plots of \code{lefkoProj} objects.
#' Acts as a convenient wrapper for the \code{plot.default()} function.
#' 
#' @name plot.lefkoProj
#' 
#' @param x A \code{lefkoProj} object.
#' @param variable The focus variable of the plot to produce. Defaults to
#' \code{"popsize"}, which produces line plots of the \code{popsize} element in
#' object \code{x}.
#' @param style A string denoting ther kind of plot to produce. Currently
#' limited to \code{"timeseries"}, which shows \code{variable} against time on
#' the x axis. Other choices include \code{"statespace"}, which plots
#' \code{variable} at one time on the x axis against the same variable in the
#' next time on the y axis.
#' @param repl The replicate to plot. Defaults to \code{"all"}, in which case
#' all replicates are plotted.
#' @param patch The patch to plot, as labeled in the \code{labels} element in
#' object \code{x}. Defaults to \code{"pop"}, in which case only the final
#' population-level projection is plotted. Can also be set to \code{"all"}, in
#' which case projections for all patches and population in the \code{labels}
#' element are plotted.
#' @param auto_ylim A logical value indicating whether the maximum of the y axis
#' should be determined automatically. Defaults to \code{TRUE}, but reverts to
#' \code{FALSE} if any setting for \code{ylim} is given.
#' @param auto_col A logical value indicating whether to shift the color of
#' lines associated with each patch automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{col} is given.
#' @param auto_lty A logical value indicating whether to shift the line type
#' associated with each replicate automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{lty} is given.
#' @param auto_title A logical value indicating whether to add a title to each
#' plot. The plot is composed of the concatenated population and patch names.
#' Defaults to \code{FALSE}.
#' @param ... Other parameters used by functions \code{plot.default()} and
#' \code{lines()}.
#' 
#' @return A plot of the results of a \code{\link{projection3}()} run.
#' 
#' @section Notes:
#' Output plots are currently limited to time series and state space plots of
#' population size.
#' 
#' The default settings will preferentially plot any projections marked as
#' \code{0} in the \code{patch} portion of the \code{labels} element of the
#' input MPM. This can produce confusing results if a mean MPM resulting from
#' the \code{lmean()} function is used as input and the \code{add_mean} setting
#' is set to the default, which is \code{TRUE}.
#' 
#' @examples
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
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
#'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
#'   repmatrix = lathrepm, supplement = lathsupp3, yearcol = "year2",
#'   indivcol = "individ")
#' 
#' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
#' plot(lathproj)
#' 
#' @export
plot.lefkoProj <- function(x, variable = "popsize", style = "time",
  repl = "all", patch = "pop", auto_ylim = TRUE, auto_col = TRUE,
  auto_lty = TRUE, auto_title = FALSE, ...) {
  
  appended <- FALSE
  
  further_args <- list(...)
  
  if (length(further_args) == 0) further_args <- list()
  
  if (!is.element("type", names(further_args))) {
    further_args$type <- "l"
  }
  if (is.element("col", names(further_args))) {
    auto_col <- FALSE
  }
  if (is.element("lty", names(further_args))) {
    auto_lty <- FALSE
  }
  if (is.element("ylim", names(further_args))) {
    auto_ylim <- FALSE
  }
  if (is.element("main", names(further_args))) {
    auto_title <- FALSE
  }
  basal_args <- further_args
  
  actual_patches <- c(1:length(x$labels$patch))
  if (is.matrix(x$control)) {
    appended <- TRUE
    actual_replicates <- apply(as.matrix(actual_patches), 1, function(X) {
      num_reps <- sum(x$control[which(x$control[, 1] == X), 2])
    })
  } else {
    actual_replicates <- rep(x$control[1], length(x$labels$patch))
  }
  
  if (length(grep("pop", variable)) > 0 | length(grep("size", variable)) > 0) {
    variable <- "popsize"
  }
  
  if (length(grep("stat", style)) > 0) {
    style <- "statespace"
    
    if (!is.element("ylab", names(further_args))) {
      further_args$ylab <- "State in time t+1"
    }
    if (!is.element("xlab", names(further_args))) {
      further_args$xlab <- "State in time t"
    }
  } else if (length(grep("tim", style)) > 0) {
    style <- "timeseries"
    
    if (!is.element("ylab", names(further_args))) {
      further_args$ylab <- "Population size"
    }
    if (!is.element("xlab", names(further_args))) {
      further_args$xlab <- "Time"
    }
  }
  
  if (all(is.na(patch))) {
    patch <- actual_patches
  } else if (is.character(patch)) {
    if (any(grep("al", patch))) {
      patch <- actual_patches
    } else if (length(patch) == 1) {
      if (grep("po", patch)) {
        if (length(actual_patches) == 1) {
          patch <- actual_patches
        } else {
          patch <- which(x$labels$patch == 0)
        }
      }
    } else {
      patch <- as.numeric(patch)
      
      if (any(is.na(patch))) {
        stop("Setting patch not understood. Please do not combine text with numbers.",
          call. = FALSE)
      }
      if (!all(is.element(patch, actual_patches))) {
        stop("Setting patch not understood. Please check the labels element of
          the input lefkoProj object.", call. = FALSE)
      }
    }
  } else if (is.numeric(patch)) {
    if (any(!is.element(patch, actual_patches))) {
      stop("Setting patch not understood. Please check the labels element of
          your input lefkoProj object.", call. = FALSE)
    }
  }
  
  if (is.character(repl)) {
    if (any(grep("al", repl))) {
      repl <- actual_replicates
    } else {
      stop("Setting repl not understood.", call. = FALSE)
    }
  }
  if (variable != "popsize") {
    stop("Function plot.lefkoProj() does not currently handle plots of elements
      other than `popsize`.", call. = FALSE)
  }
  
  used_col <- 1
  
  for (i in patch) {
    if (auto_title) {
      used_string <- paste("pop", x$labels[i, 1], "patch", x$labels[i, 2])
      further_args$main <- used_string
    }
    
    used_lty <- 1
    
    if (auto_ylim) {
      further_args$ylim <- c(0, max(x$pop_size[[i]], na.rm = TRUE))
    }
    if (auto_col) {
      further_args$col <- used_col
      basal_args$col <- used_col
    }
    if (auto_lty) {
      further_args$lty <- used_lty
      basal_args$lty <- used_lty
    }
      
    if (style == "timeseries") {
      c_xy <- xy.coords(x = c(1:length(x$pop_size[[i]][1,])), y = x$pop_size[[i]][1,])
      further_args$x <- c_xy
      
    } else if (style == "statespace") {
      c_xy <- xy.coords(x = x$pop_size[[i]][1,1:(dim(x$pop_size[[i]])[2] - 1)],
        y = x$pop_size[[i]][1,2:dim(x$pop_size[[i]])[2]])
      further_args$x <- c_xy
    }
    
    do.call("plot.default", further_args)
    
    if (repl[i] > 1) {
      for (j in c(2:actual_replicates[i])) {
        if (style == "timeseries") {
          basal_args$x <- c(1:length(x$pop_size[[i]][j,]))
          basal_args$y <- x$pop_size[[i]][j,]
          
        } else if (style == "statespace") {
          basal_args$y <- x$pop_size[[i]][j,2:dim(x$pop_size[[i]])[2]]
          basal_args$x <- x$pop_size[[i]][j,1:(dim(x$pop_size[[i]])[2] - 1)]
        
        }
        do.call("lines", basal_args)
      }
      used_lty <- used_lty + 1;
    }
    used_col <- used_col + 1;
    if (used_col > length(palette())) used_col <- 1;
  }
}

