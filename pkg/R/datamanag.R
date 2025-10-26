#' Create Historical Vertical Data Frame from Horizontal Data Frame
#' 
#' Function \code{verticalize3()} returns a vertically formatted demographic
#' data frame organized to create historical projection matrices, given a
#' horizontally formatted input data frame. It also handles stage assignments
#' if given an appropriate stageframe.
#' 
#' @name verticalize3
#' 
#' @param data The horizontal data file. A valid data frame is required as
#' input.
#' @param noyears The number of years or observation occasions in the dataset. A
#' valid integer is required as input.
#' @param firstyear The first year or occasion of observation. Defaults to
#' \code{1}.
#' @param popidcol A variable name or column number corresponding to the 
#' identity of the population for each individual.
#' @param patchidcol A variable name or column number corresponding to the 
#' identity of the patch or subpopulation for each individual, if patches have
#' been designated within populations.
#' @param individcol A variable name or column number corresponding to the 
#' identity of each individual.
#' @param blocksize The number of variables corresponding to each occasion in
#' the input dataset designated in \code{data}, if a set pattern of variables is
#' used for each observation occasion in the data frame used as input. If such a
#' pattern is not used, and all variable names are properly noted as character
#' vectors in the other input variables, then this may be set to \code{NA}.
#' Defaults to \code{NA}.
#' @param xcol A variable name(s) or column number(s) corresponding to the X 
#' coordinate of each individual, or of each individual at each occasion, in
#' Cartesian space. Can refer to the only instance, the first instance, or all
#' instances of X variables. In the last case, the values should be entered as a
#' vector.
#' @param ycol A variable name(s) or column number(s) corresponding to the Y
#' coordinate of each individual, or of each individual at each occasion, in
#' Cartesian space. Can refer to the only instance, the first instance, or all
#' instances of Y variables. In the last case, the values should be entered as a
#' vector.
#' @param juvcol A variable name(s) or column number(s) that marks individuals
#' in immature stages within the dataset. This function assumes that immature
#' individuals are identified in this variable marked with a number equal to or
#' greater than \code{1}, and that mature individuals are marked as \code{0} or
#' \code{NA}. Can refer to the first instance, or all instances of these
#' variables. In the latter case, the values should be entered as a vector.
#' @param sizeacol A variable name(s) or column number(s) corresponding to the
#' size entry associated with the first year or observation occasion in the
#' dataset. Can refer to the first instance, or all instances of these
#' variables. In the latter case, the values should be entered as a vector.
#' This variable should refer to the first size variable in the stageframe,
#' unless \code{stagesize = "sizeadded"}.
#' @param sizebcol A second variable name(s) or column number(s) corresponding
#' to the size entry associated with the first year or observation occasion in
#' the dataset. Can refer to the first instance, or all instances of these
#' variables. In the latter case, the values should be entered as a vector.
#' This variable should refer to the second size variable in the stageframe,
#' unless \code{stagesize = "sizeadded"}.
#' @param sizeccol A third variable name(s) or column number(s) corresponding to
#' the size entry associated with the first year or observation occasion in the
#' dataset. Can refer to the first instance, or all instances of these variables.
#' In the latter case, the values should be entered as a vector. This variable
#' should refer to the third size variable in the stageframe, unless
#' \code{stagesize = "sizeadded"}.
#' @param repstracol A variable name(s) or column number(s) corresponding to the
#' production of reproductive structures, such as flowers, associated with the 
#' first year or observation period in the input dataset. This can be binomial 
#' or count data, and is used to analyze the probability of reproduction. Can
#' refer to the first instance, or all instances of these variables. In the
#' latter case, the values should be entered as a vector.
#' @param repstrbcol A second variable name(s) or column number(s) corresponding
#' to the production of reproductive structures, such as flowers, associated
#' with the first year or observation period in the input dataset. This can be 
#' binomial or count data, and is used to analyze the probability of
#' reproduction. Can refer to the first instance, or all instances of these
#' variables. In the latter case, the values should be entered as a vector.
#' @param fecacol A variable name(s) or column number(s) denoting fecundity
#' associated with the first year or observation occasion in the input dataset.
#' This may represent egg counts, fruit counts, seed production, etc. Can refer
#' to the first instance, or all instances of these variables. In the latter
#' case, the values should be entered as a vector.
#' @param fecbcol A second variable name(s) or column number(s) denoting
#' fecundity associated with the first year or observation occasion in the input
#' dataset. This may represent egg counts, fruit counts, seed production, etc.
#' Can refer to the first instance, or all instances of these variables. In the
#' latter case, the values should be entered as a vector.
#' @param indcovacol A variable name(s) or column number(s) corresponding to an
#' individual covariate to be used in analysis. Can refer to the only instance,
#' the first instance, or all instances of these variables. In the last case,
#' the values should be entered as a vector.
#' @param indcovbcol A variable name(s) or column number(s) corresponding to an
#' individual covariate to be used in analysis. Can refer to the only instance,
#' the first instance, or all instances of these variables. In the last case,
#' the values should be entered as a vector.
#' @param indcovccol A second variable name(s) or column number(s) corresponding
#' to an individual covariate to be used in analysis. Can refer to the only
#' instance, the first instance, or all instances of these variables. In the
#' last case, the values should be entered as a vector.
#' @param aliveacol Variable name(s) or column number(s) providing information
#' on whether an individual is alive at a given occasion. If used, living status
#' must be designated as binomial (living = \code{1}, dead = \code{0}). Can
#' refer to the first instance of a living status variable in the dataset, or
#' a full vector of all living status variables in temporal order.
#' @param deadacol Variable name(s) or column number(s) providing information on
#' whether an individual is alive at a given occasion. If used, dead status must
#' be designated as binomial (dead = \code{1}, living = \code{0}).  Can refer to
#' the first instance of a dead status variable in the dataset, or a full vector
#' of all dead status variables in temporal order.
#' @param obsacol A variable name(s) or column number(s) providing information
#' on whether an individual is in an observable stage at a given occasion. If
#' used, observation status must be designated as binomial (observed = \code{1}, 
#' not observed = \code{0}). Can refer to the first instance of an observation
#' status variable in the dataset, or a full vector of all observation status
#' variables in temporal order.
#' @param nonobsacol A variable name(s) or column number(s) providing
#' information on whether an individual is in an unobservable stage at a given
#' occasion. If used, observation status must be designated as binomial (not
#' observed = \code{1}, observed = \code{0}). Can refer to the first instance of
#' a non-observation status variable in the dataset, or a full vector of all
#' non-observation status variables in temporal order.
#' @param censorcol A variable name(s) or column number(s) corresponding to the
#' first entry of a censor variable, used to distinguish between entries to use
#' and entries not to use, or to designate entries with special issues that
#' require further attention. Can refer to the first instance of a censor status
#' variable in the dataset, or a full vector of all censor status variables in
#' temporal order. Can also refer to a single censor status variable used for
#' the entire individual, if \code{singlecensor = TRUE}.
#' @param repstrrel This is a scalar multiplier on variable \code{repstrbcol} to
#' make it equivalent to \code{repstracol}. This can be useful if two 
#' reproductive status variables have related but unequal units, for example if
#' \code{repstracol} refers to one-flowered stems while \code{repstrbcol} refers
#' to two-flowered stems. Defaults to \code{1}.
#' @param fecrel This is a scalar multiplier on variable \code{fecbcol} to make
#' it equivalent to \code{fecacol}. This can be useful if two fecundity 
#' variables have related but unequal units. Defaults to \code{1}.
#' @param stagecol Optional variable name(s) or column number(s) corresponding
#' to life history stage at a given occasion. Can refer to the first instance of
#' a stage identity variable in the dataset, or a full vector of all stage
#' identity variables in temporal order.
#' @param stageassign The stageframe object identifying the life history model
#' being operationalized. Note that if \code{stagecol} is provided, then this
#' stageframe is not used for stage designation.
#' @param stagesize A variable name or column number describing which size 
#' variable to use in stage estimation. Defaults to NA, and can also take 
#' \code{sizea}, \code{sizeb}, \code{sizec}, \code{sizeab}, \code{sizebc},
#' \code{sizeac}, \code{sizeabc}, or \code{sizeadded}, depending on
#' which size variable within the input dataset is chosen. Note that the
#' variable(s) chosen should be presented in the order of the primary,
#' secondary, and tertiary variables in the stageframe input with
#' \code{stageassign}. For example, choosing \code{sizeb} assumes that this size
#' is the primary variable in the stageframe.
#' @param censorkeep The value of the censor variable identifying data to be
#' included in analysis. Defaults to \code{0}, but may take any value including
#' \code{NA}. Note that if \code{NA} is the value to keep, then this function
#' will alter all \code{NA}s to \code{0} values, and all other values to
#' \code{1}, treating \code{0} as the new value to keep.
#' @param censorRepeat A logical value indicating whether the censor variable
#' is a single column, or whether it repeats across occasion blocks. Defaults to
#' \code{FALSE}.
#' @param censor A logical variable determining whether the output data should 
#' be censored using the variable defined in \code{censorcol}. Defaults to 
#' \code{FALSE}.
#' @param coordsRepeat A logical value indicating whether X and Y coordinates
#' correspond to single X and Y columns. If \code{TRUE}, then each observation
#' occasion has its own X and Y variables. Defaults to \code{FALSE}.
#' @param spacing The spacing at which density should be estimated, if density
#' estimation is desired and X and Y coordinates are supplied. Given in the same
#' units as those used in the X and Y coordinates given in \code{xcol} and 
#' \code{ycol}. Defaults to \code{NA}.
#' @param NAas0 If \code{TRUE}, then all \code{NA} entries for size and
#' fecundity variables will be set to 0. This can help increase the sample size
#' analyzed by \code{\link{modelsearch}()}, but should only be used when it is
#' clear that this substitution is biologically realistic. Defaults to
#' \code{FALSE}.
#' @param NRasRep If \code{TRUE}, then will treat non-reproductive but mature 
#' individuals as reproductive during stage assignment. This can be useful when
#' a MPM is desired without separation of reproductive and non-reproductive but
#' mature stages of the same size. Only used if \code{stageassign} is set to a
#' stageframe. Defaults to \code{FALSE}.
#' @param NOasObs If \code{TRUE}, then will treat individuals that are
#' interpreted as not observed in the dataset as though they were observed
#' during stage assignment. This can be useful when a MPM is desired without
#' separation of observable and unobservable stages. Only used if
#' \code{stageassign} is set to a stageframe. Defaults to \code{FALSE}.
#' @param prebreeding A logical term indicating whether the life history model
#' is pre-breeding. If so, then \code{1} is added to all ages. Defaults to
#' \code{TRUE}.
#' @param age_offset A number to add automatically to all values of age at time
#' \emph{t}. Defaults to \code{0}.
#' @param reduce A logical variable determining whether unused variables and 
#' some invariant state variables should be removed from the output dataset.
#' Defaults to \code{TRUE}.
#' @param a2check A logical variable indicating whether to retain all data with
#' living status at occasion \emph{t}. Defaults to \code{FALSE}, in which case
#' data for occasions in which the individual is not alive in time \emph{t} is
#' not retained. This option should be kept \code{FALSE}, except to inspect
#' potential errors in the dataset.
#' @param quiet A logical variable indicating whether to silence warnings.
#' Defaults to \code{FALSE}.
#' 
#' @return If all inputs are properly formatted, then this function will output
#' a historical vertical data frame (class \code{hfvdata}), meaning that the
#' output data frame will have three consecutive occasions of size and
#' reproductive data per individual per row. This data frame is in standard
#' format for all functions used in \code{lefko3}, and so can be used without
#' further modification.
#' 
#' Variables in this data frame include the following:
#' \item{rowid}{Unique identifier for the row of the data frame.}
#' \item{popid}{Unique identifier for the population, if given.}
#' \item{patchid}{Unique identifier for patch within population, if given.}
#' \item{individ}{Unique identifier for the individual.}
#' \item{year2}{Year or time at occasion \emph{t}.}
#' \item{firstseen}{Occasion of first observation.}
#' \item{lastseen}{Occasion of last observation.}
#' \item{obsage}{Observed age in occasion \emph{t}, assuming first observation
#' corresponds to age = 0.}
#' \item{obslifespan}{Observed lifespan, given as \code{lastseen - firstseen + 1}.}
#' \item{xpos1,xpos2,xpos3}{X position in Cartesian space in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively, if provided.}
#' \item{ypos1,ypos2,ypos3}{Y position in Cartesian space in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively, if provided.}
#' \item{sizea1,sizea2,sizea3}{Main size measurement in occasions \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizeb1,sizeb2,sizeb3}{Secondary size measurement in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizec1,sizec2,sizec3}{Tertiary measurement in occasions \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{size1added,size2added,size3added}{Sum of primary, secondary, and 
#' tertiary size measurements in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, 
#' respectively.}
#' \item{repstra1,repstra2,repstra3}{Main numbers of reproductive structures in
#' occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstrb1,repstrb2,repstrb3}{Secondary numbers of reproductive 
#' structures in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstr1added,repstr2added,repstr3added}{Sum of primary and secondary
#' reproductive structures in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, 
#' respectively.}
#' \item{feca1,feca2,feca3}{Main numbers of offspring in occasions \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecb1,fecb2, fecb3}{Secondary numbers of offspring in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fec1added,fec2added,fec3added}{Sum of primary and secondary fecundity
#' in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{censor1,censor2,censor3}{Censor state values in occasions \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{juvgiven1,juvgiven2,juvgiven3}{Binomial variable indicating whether
#' individual is juvenile in occasions \emph{t}-1, \emph{t}, and \emph{t}+1.
#' Only given if \code{juvcol} is provided.}
#' \item{obsstatus1,obsstatus2,obsstatus3}{Binomial observation state in
#' occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstatus1,repstatus2,repstatus3}{Binomial reproductive state in
#' occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecstatus1,fecstatus2,fecstatus3}{Binomial offspring production state
#' in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{matstatus1,matstatus2,matstatus3}{Binomial maturity state in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{alive1,alive2,alive3}{Binomial state as alive in occasions \emph{t}-1,
#'  \emph{t}, and \emph{t}+1, respectively.}
#' \item{density}{Radial density of individuals per unit designated in
#' \code{spacing}. Only given if \code{spacing} is not NA.}
#' 
#' @section Notes:
#' In some datasets on species with unobservable stages, observation status
#' (\code{obsstatus}) might not be inferred properly if a single size variable
#' is used that does not yield sizes greater than 0 in all cases in which
#' individuals were observed. Such situations may arise, for example, in plants
#' when leaf number is the dominant size variable used, but individuals
#' occasionally occur with inflorescences but no leaves. In this instances,
#' it helps to mark related variables as \code{sizeb} and \code{sizec}, because
#' observation status will be interpreted in relation to all 3 size variables.
#' Further analysis can then utilize only a single size variable, of the user's
#' choosing. Similar issues can arise in reproductive status (\code{repstatus}).
#' 
#' Juvenile designation should only be used when juveniles fall outside of the
#' size classification scheme used in determining stages. If juveniles are to be
#' size classified along the size spectrum that adults also fall on, then
#' it is best to treat juveniles as mature but not reproductive.
#' 
#' Warnings that some individuals occur in state combinations that do not match
#' any stages in the stageframe used to assign stages are common when first
#' working with a dataset. Typically, these situations can be identified as
#' \code{NoMatch} entries in \code{stage3}, although such entries may crop up in
#' \code{stage1} and \code{stage2}, as well. In rare cases, these warnings will
#' arise with no concurrent \code{NoMatch} entries, which indicates that the
#' input dataset contained conflicting state data at once suggesting that the
#' individual is in some stage but is also dead. The latter is removed if the
#' conflict occurs in occasion \emph{t} or \emph{t}-1, as only living entries
#' are allowed in time \emph{t} and time \emph{t}-1 may involve living entries
#' as well as non-living entries immediately prior to birth.
#' 
#' Care should be taken to avoid variables with negative values indicating size,
#' fecundity, or reproductive or observation status. Negative values can be
#' interpreted in different ways, typically reflecting estimation through other
#' algorithms rather than actual measured data. Variables holding negative
#' values can conflict with data management algorithms in ways that are
#' difficult to predict.
#' 
#' Unusual errors (e.g. \code{"Error in .pfj..."}) may occur in cases where the
#' variables are improperly passed, where seemingly numeric variables include
#' text, or where the \code{blocksize} is improperly set.
#' 
#' Density estimation is performed as a count of individuals alive and within
#' the radius specified in \code{spacing} of the respective individual at some
#' point in time.
#' 
#' If a censor variable is included for each monitoring occasion, and the
#' \code{blocksize} option is set, then the user must set
#' \code{censorRepeat = TRUE} in order to censor the correct transitions.
#' Failing this step will likely lead to the loss of a large portion of the data
#' as all data for entire individuals will be excluded.
#' 
#' @examples
#' # Lathyrus example using blocksize - when repeated patterns exist in variable
#' # order
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
#' # Cypripedium example using partial repeat patterns with blocksize and part
#' # explicit variable name cast
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
#'   repstracol = c("Inf.04", "Inf.05", "Inf.06", "Inf.07", "Inf.08", "Inf.09"),
#'   repstrbcol = c("Inf2.04", "Inf2.05", "Inf2.06", "Inf2.07", "Inf2.08", "Inf2.09"), 
#'   fecacol = "Pod.04", stageassign = cypframe_raw, stagesize = "sizeadded",
#'   NAas0 = TRUE, NRasRep = TRUE)
#' 
#' @export
verticalize3 <- function(data, noyears, firstyear = 1, popidcol = 0,
  patchidcol = 0, individcol= 0, blocksize = NA, xcol = 0, ycol = 0, juvcol = 0,
  sizeacol, sizebcol = 0, sizeccol = 0, repstracol = 0, repstrbcol = 0,
  fecacol = 0, fecbcol = 0, indcovacol = 0, indcovbcol = 0, indcovccol = 0,
  aliveacol = 0, deadacol = 0, obsacol = 0, nonobsacol = 0, censorcol = 0,
  repstrrel = 1, fecrel = 1, stagecol = 0, stageassign = NA, stagesize = NA,
  censorkeep = 0, censorRepeat = FALSE, censor = FALSE,
  coordsRepeat = FALSE, spacing = NA, NAas0 = FALSE, NRasRep = FALSE,
  NOasObs = FALSE, prebreeding = TRUE, age_offset = 0, reduce = TRUE,
  a2check = FALSE, quiet = FALSE) {
  
  stassign <- rowid <- alive2 <- indataset <- censor1 <- censor2 <- NULL
  censor3 <- censbool <- NULL
  
  popid <- patchid <- individ <- NA
  RepasObs <- FALSE
  
  #This first section tests the input for valid entries
  data.limits <- dim(data)
  
  if (length(blocksize) != 1) {
    stop("The blocksize option must equal a single number, or NA.", call. = FALSE)
  }
  
  if (is.na(blocksize)) {
    blocksize <- 0
  }
  
  if (!all(is.logical(c(censorRepeat, censor, coordsRepeat, NAas0, NRasRep, 
      NOasObs, reduce, a2check, quiet, prebreeding)))) {
    stop("Some logical variables have been set to non-logical values.", call. = FALSE)
  }
  
  if (prebreeding) age_offset <- age_offset + 1
  
  if (is.character(popidcol)) {
    if (is.element(popidcol, names(data))) {
      true.popidcol <- which(names(data) == popidcol)
      popidcol <- true.popidcol
    } else {
      stop("Please enter popidcol exactly as it appears in the dataset.", 
        call. = FALSE)}
  }
  
  if (is.character(patchidcol)) {
    if (is.element(patchidcol, names(data))) {
      true.patchidcol <- which(names(data) == patchidcol)
      patchidcol <- true.patchidcol
    } else {
      stop("Please enter patchidcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(individcol)) {
    if (is.element(individcol, names(data))) {
      true.individcol <- which(names(data) == individcol)
      individcol <- true.individcol
    } else {
      stop("Please enter individcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (all(is.character(xcol))) {
    if (all(is.element(xcol, names(data)))) {
      true.xcol <- match(xcol, names(data))
      xcol <- true.xcol
    } else {
      stop("Please enter xcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(xcol < 0) | any(xcol > data.limits[2])) {
    stop("Variable xcol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(ycol))) {
    if (all(is.element(ycol, names(data)))) {
      true.ycol <- match(ycol, names(data))
      ycol <- true.ycol
    } else {
      stop("Please enter ycol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(ycol < 0) | any(ycol > data.limits[2])) {
    stop("Variable ycol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(juvcol))) {
    if (all(is.element(juvcol, names(data)))) {
      true.juvcol <- match(juvcol, names(data))
      juvcol <- true.juvcol
    } else {
      stop("Please enter juvcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(juvcol < 0) | any(juvcol > data.limits[2])) {
    stop("Variable juvcol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(sizeacol))) {
    if (all(is.element(sizeacol, names(data)))) {
      true.sizeacol <- match(sizeacol, names(data))
      sizeacol <- true.sizeacol
    } else {
      stop("Please enter sizeacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(sizeacol < 0) | any(sizeacol > data.limits[2])) {
    stop("Variable sizeacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(sizebcol))) {
    if (all(is.element(sizebcol, names(data)))) {
      true.sizebcol <- match(sizebcol, names(data))
      sizebcol <- true.sizebcol
    } else {
      stop("Please enter sizebcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(sizebcol < 0) | any(sizebcol > data.limits[2])) {
    stop("Variable sizebcol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(sizeccol))) {
    if (all(is.element(sizeccol, names(data)))) {
      true.sizeccol <- match(sizeccol, names(data))
      sizeccol <- true.sizeccol
    } else {
      stop("Please enter sizeccol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(sizeccol < 0) | any(sizeccol > data.limits[2])) {
    stop("Variable sizeccol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(repstracol))) {
    if (all(is.element(repstracol, names(data)))) {
      true.repstracol <- match(repstracol, names(data))
      repstracol <- true.repstracol
    } else {
      stop("Please enter repstracol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(repstracol < 0) | any(repstracol > data.limits[2])) {
    stop("Variable repstracol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(repstrbcol))) {
    if (all(is.element(repstrbcol, names(data)))) {
      true.repstrbcol <- match(repstrbcol, names(data))
      repstrbcol <- true.repstrbcol
    } else {
      stop("Please enter repstrbcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(repstrbcol < 0) | any(repstrbcol > data.limits[2])) {
    stop("Variable repstrbcol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(fecacol))) {
    if (all(is.element(fecacol, names(data)))) {
      true.fecacol <- match(fecacol, names(data))
      fecacol <- true.fecacol
    } else {
      stop("Please enter fecacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(fecacol < 0) | any(fecacol > data.limits[2])) {
    stop("Variable fecacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(fecbcol))) {
    if (all(is.element(fecbcol, names(data)))) {
      true.fecbcol <- match(fecbcol, names(data))
      fecbcol <- true.fecbcol
    } else {
      stop("Please enter fecbcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(fecbcol < 0) | any(fecbcol > data.limits[2])) {
    stop("Variable fecbcol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(indcovacol))) {
    if (all(is.element(indcovacol, names(data)))) {
      true.indcovacol <- match(indcovacol, names(data))
      indcovacol <- true.indcovacol
    } else {
      stop("Please enter indcovacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(indcovacol < 0) | any(indcovacol > data.limits[2])) {
    stop("Variable indcovacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(indcovbcol))) {
    if (all(is.element(indcovbcol, names(data)))) {
      true.indcovbcol <- match(indcovbcol, names(data))
      indcovbcol <- true.indcovbcol
    } else {
      stop("Please enter indcovbcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(indcovbcol < 0) | any(indcovbcol > data.limits[2])) {
    stop("Variable indcovbcol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(indcovccol))) {
    if (all(is.element(indcovccol, names(data)))) {
      true.indcovccol <- match(indcovccol, names(data))
      indcovccol <- true.indcovccol
    } else {
      stop("Please enter indcovccol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(indcovccol < 0) | any(indcovccol > data.limits[2])) {
    stop("Variable indcovccol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(aliveacol))) {
    if (all(is.element(aliveacol, names(data)))) {
      true.aliveacol <- match(aliveacol, names(data))
      aliveacol <- true.aliveacol
    } else {
      stop("Please enter aliveacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(aliveacol < 0) | any(aliveacol > data.limits[2])) {
    stop("Variable aliveacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(deadacol))) {
    if (all(is.element(deadacol, names(data)))) {
      true.deadacol <- match(deadacol, names(data))
      deadacol <- true.deadacol
    } else {
      stop("Please enter deadacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(deadacol < 0) | any(deadacol > data.limits[2])) {
    stop("Variable deadacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(obsacol))) {
    if (all(is.element(obsacol, names(data)))) {
      true.obsacol <- match(obsacol, names(data))
      obsacol <- true.obsacol
    } else {
      stop("Please enter obsacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(obsacol < 0) | any(obsacol > data.limits[2])) {
    stop("Variable obsacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(nonobsacol))) {
    if (all(is.element(nonobsacol, names(data)))) {
      true.nonobsacol <- match(nonobsacol, names(data))
      nonobsacol <- true.nonobsacol
    } else {
      stop("Please enter nonobsacol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(nonobsacol < 0) | any(nonobsacol > data.limits[2])) {
    stop("Variable nonobsacol designation is out of bounds.", call. = FALSE)
  }
  
  if (all(is.character(censorcol))) {
    if (all(is.element(censorcol, names(data)))) {
      true.censorcol <- match(censorcol, names(data))
      censorcol <- true.censorcol
    } else {
      stop("Please enter censorcol exactly as it appears in the dataset.",
        call. = FALSE)}
  } else if (any(censorcol < 0) | any(censorcol > data.limits[2])) {
    stop("Variable censorcol designation is out of bounds.", call. = FALSE)
  }
  
  stagesize <- tolower(stagesize)
  if(!all(is.na(stageassign))) {
    if(!is(stageassign, "stageframe")) {
      stop("The stageassign option can only take NA or a stageframe object as input.",
        call. = FALSE)
    }
    stagesize_options <- c("sizea", "sizeb", "sizec", "sizeadded", "sizeab",
      "sizeac", "sizebc", "sizeabc")
    
    if(!is.element(stagesize, stagesize_options)) {
      stop("The stagesize option must equal 'NA', 'sizea', 'sizeb', 'sizec', 'sizeab',
        'sizeac', 'sizebc', 'sizeabc', or 'sizeadded'.", call. = FALSE)
    }
  }
  
  if (all(is.character(stagecol))) {
    if (all(is.element(stagecol, names(data)))) {
      true.stagecol <- match(stagecol, names(data))
      stagecol <- true.stagecol
    } else {stop("Please enter stagecol exactly as it appears in the dataset.", call. = FALSE)}
  } else if (any(stagecol < 0) | any(stagecol > data.limits[2])) {
    stop("Variable stagecol designation is out of bounds.", call. = FALSE)
  }
  
  input.lengths <- c(length(xcol), length(ycol), length(juvcol), length(sizeacol),
    length(sizebcol), length(sizeccol), length(repstracol), length(repstrbcol),
    length(fecacol), length(fecbcol), length(indcovacol), length(indcovbcol),
    length(indcovccol), length(aliveacol), length(deadacol), length(obsacol),
    length(nonobsacol), length(censorcol), length(stagecol))
  
  if (any(!is.element(input.lengths, c(1, noyears)))) {
    stop("All input variables must be either single constants, or vectors of length noyears.",
      call. = FALSE)
  }
  
  # Here we will modify our approach to verticalization based on the input stageframe
  if (!all(is.na(stageassign)) & all(stagecol == 0)) {
    
    stassign <- TRUE 
    
    if (stagesize == "sizeabc") {
      stagesizecol <- 8
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size == 0),
        intersect(which(stageassign$size_b == 0), which(stageassign$size_c == 0))))
    } else if (stagesize == "sizebc") {
      stagesizecol <- 7
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size_b == 0),
        which(stageassign$size_c == 0)))
    } else if (stagesize == "sizeac") {
      stagesizecol <- 6
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size == 0),
        which(stageassign$size_c == 0)))
    } else if (stagesize == "sizeab") {
      stagesizecol <- 5
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size == 0),
        which(stageassign$size_b == 0)))
    } else if (stagesize == "sizeadded") {
      stagesizecol <- 4
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size == 0))
    } else if (stagesize == "sizec") {
      stagesizecol <- 3
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size_c == 0))
    } else if (stagesize == "sizeb") {
      stagesizecol <- 2
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size_b == 0))
    } else {
      stagesizecol <- 1
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size == 0))
    }
    
    repcheck <- intersect(repcheck1, which(stageassign$obsstatus == 1))
    
    if (length(repcheck) > 0) {
      if (!quiet) {
        message("Stageframe indicates the presence of an observeable, reproductive stage with a size of 0.")
      }
      RepasObs <- TRUE
    } else if (length(repcheck1) > 0) {
      if (!quiet) {
        message("Stageframe indicates the presence of an unobserveable, reproductive stage with a size of 0.")
      }
    }
    
  } else {
    stassign <- FALSE
    
    stageassign <- as.data.frame(matrix(c(NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA), ncol = 29),
      stringsAsFactors = FALSE)
    names(stageassign) <- c("stage", "size", "size_b", "size_c", "min_age", "max_age",
      "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
      "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center", "sizebin_width",
      "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
      "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
      "group", "comments")
    
    stagesizecol <- 0
  }
  
  if (!is.na(spacing)) {
    if (any(xcol == 0) | any(ycol == 0)) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  if (censor) {
    if (length(censorcol) == 1 & blocksize != 0) {
      fullcenvec <- as.vector(apply(as.matrix(c(1:noyears)), 1, function(X) {
        c(data[,(censorcol + (X - 1) * blocksize)])
      }))
    } else {
      fullcenvec <- unlist(data[,censorcol])
    }
    
    if (!is.element(censorkeep, fullcenvec)) {
      stop("Please enter a valid value for censorkeep. This value should occur in the censor variable within the dataset.", 
        call. = FALSE)
    }
  }
  
  if (is.na(censorkeep) & censor) { # This section checks to see if NA is the censor value to keep
    censbool <- TRUE
    censorkeep <- 0
  } else {
    censbool <- FALSE
  }
  
  popdatareal <- .pfj(data, stageassign, noyears, firstyear, (popidcol - 1),
    (patchidcol - 1), (individcol - 1), blocksize, (xcol - 1), (ycol - 1),
    (juvcol - 1), (sizeacol - 1), (sizebcol - 1), (sizeccol - 1),
    (repstracol - 1), (repstrbcol - 1), (fecacol - 1), (fecbcol - 1),
    (indcovacol - 1), (indcovbcol - 1), (indcovccol - 1), (aliveacol - 1),
    (deadacol - 1), (obsacol - 1), (nonobsacol - 1), (censorcol - 1),
    (stagecol - 1), repstrrel, fecrel, NAas0, NRasRep, RepasObs, NOasObs,
    stassign, stagesizecol, censorkeep, censbool, censorRepeat, coordsRepeat,
    a2check, reduce, quiet)
  
  if ("censor2" %in% colnames(popdatareal) & censor) {
    popdatareal <- subset(popdatareal, censor1 == censorkeep)
    popdatareal <- subset(popdatareal, censor2 == censorkeep)
    popdatareal <- subset(popdatareal, censor3 == censorkeep)
  }
  
  if (!is.na(spacing)) {
    popdatareal$density <- .density3(popdatareal, which(names(popdatareal) == "xpos2"),
      which(names(popdatareal) == "ypos2"), which(names(popdatareal) == "year2"), spacing)
  }
  
  popdatareal$obsage <- popdatareal$obsage + age_offset
  
  if (censor & !censorRepeat & blocksize > 0) {
    message("hfv dataset has been censored using a non-repeating censor variable.
      Please check the resulting dataset with the summary_hfv() function to make
      sure that not too much has been removed from the original dataset.")
  }
  
  return(popdatareal)
}

#' Create Historical Vertical Data Frame from Ahistorical Vertical Data Frame
#' 
#' Function \code{historicalize3()} returns a vertically formatted demographic
#' data frame organized to create historical projection matrices, given a
#' vertically but ahistorically formatted data frame. This data frame is in
#' standard \code{hfvdata} format and can be used in all functions in the
#' package.
#' 
#' @name historicalize3
#' 
#' @param data The horizontal data file.
#' @param popidcol A variable name or column number corresponding to the
#' identity of the population for each individual.
#' @param patchidcol A variable name or column number corresponding to the
#' identity of the patch or subpopulation for each individual, if patches have
#' been designated within populations.
#' @param individcol A variable name or column number corresponding to the
#' unique identity of each individual.
#' @param year2col A variable name or column number corresponding to occasion
#' \emph{t} (year or time).
#' @param year3col A variable name or column number corresponding to occasion
#' \emph{t}+1 (year or time).
#' @param xcol A variable name or column number corresponding to the X
#' coordinate of each individual in Cartesian space.
#' @param ycol A variable name or column number corresponding to the Y
#' coordinate of each individual in Cartesian space.
#' @param sizea2col A variable name or column number corresponding to the
#' primary size entry in occasion \emph{t}.
#' @param sizea3col A variable name or column number corresponding to the
#' primary size entry in occasion \emph{t}+1.
#' @param sizeb2col A variable name or column number corresponding to the
#' secondary size entry in occasion \emph{t}.
#' @param sizeb3col A variable name or column number corresponding to the
#' secondary size entry in occasion \emph{t}+1.
#' @param sizec2col A variable name or column number corresponding to the
#' tertiary size entry in occasion \emph{t}.
#' @param sizec3col A variable name or column number corresponding to the
#' tertiary size entry in occasion \emph{t}+1.
#' @param repstra2col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in occasion \emph{t}.
#' This can be binomial or count data, and is used to in analysis of the
#' probability of reproduction.
#' @param repstra3col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in occasion
#' \emph{t}+1. This can be binomial or count data, and is used to in analysis
#' of the probability of reproduction.
#' @param repstrb2col A second variable name or column number corresponding to
#' the production of reproductive structures, such as flowers, in occasion
#' \emph{t}. This can be binomial or count data.
#' @param repstrb3col A second variable name or column number corresponding to
#' the production of reproductive structures, such as flowers, in occasion 
#' \emph{t}+1. This can be binomial or count data.
#' @param feca2col A variable name or column number corresponding to fecundity
#' in occasion \emph{t}. This may represent egg counts, fruit counts, seed 
#' production, etc.
#' @param feca3col A variable name or column number corresponding to fecundity
#' in occasion \emph{t}+1. This may represent egg counts, fruit counts, seed
#' production, etc.
#' @param fecb2col A second variable name or column number corresponding to 
#' fecundity in occasion \emph{t}. This may represent egg counts, fruit counts,
#' seed production, etc.
#' @param fecb3col A second variable name or column number corresponding to 
#' fecundity in occasion \emph{t}+1. This may represent egg counts, fruit
#' counts, seed production, etc.
#' @param indcova2col A variable name or column number corresponding to an
#' individual covariate to be used in analysis, in occasion \emph{t}.
#' @param indcova3col A variable name or column number corresponding to an
#' individual covariate to be used in analysis, in occasion \emph{t}+1.
#' @param indcovb2col A variable name or column number corresponding to a second
#' individual covariate to be used in analysis, in occasion \emph{t}.
#' @param indcovb3col A variable name or column number corresponding to a second
#' individual covariate to be used in analysis, in occasion \emph{t}+1.
#' @param indcovc2col A variable name or column number corresponding to a third
#' individual covariate to be used in analysis, in occasion \emph{t}.
#' @param indcovc3col A variable name or column number corresponding to a third
#' individual covariate to be used in analysis, in occasion \emph{t}+1.
#' @param alive2col A variable name or column number that provides information
#' on whether an individual is alive in occasion \emph{t}. If used, living
#' status must be designated as binomial (living = \code{1}, dead = \code{0}).
#' @param alive3col A variable name or column number that provides information
#' on whether an individual is alive in occasion \emph{t}+1. If used, living
#' status must be designated as binomial (living = \code{1}, dead = \code{0}).
#' @param dead2col A variable name or column number that provides information on
#' whether an individual is dead in occasion \emph{t}. If used, dead status
#' must be designated as binomial (living = \code{0}, dead = \code{1}).
#' @param dead3col A variable name or column number that provides information on
#' whether an individual is dead in occasion \emph{t}+1. If used, dead status
#' must be designated as binomial (living = \code{0}, dead = \code{1}).
#' @param obs2col A variable name or column number providing information on
#' whether an individual is in an observable stage in occasion \emph{t}. If
#' used, observation status must be designated as binomial (observed = \code{1},
#' not observed = \code{0}).
#' @param obs3col A variable name or column number providing information on
#' whether an individual is in an observable stage in occasion \emph{t}+1. If
#' used, observation status must be designated as binomial (observed = \code{1},
#' not observed = \code{0}).
#' @param nonobs2col A variable name or column number providing information on
#' whether an individual is in an unobservable stage in occasion \emph{t}. If
#' used, observation status must be designated as binomial (observed = \code{0},
#' not observed = \code{1}).
#' @param nonobs3col A variable name or column number providing information on
#' whether an individual is in an unobservable stage in occasion \emph{t}+1. If
#' used, observation status must be designated as binomial (observed = \code{0},
#' not observed = \code{1}).
#' @param repstrrel This is a scalar multiplier making the variable represented
#' by \code{repstrb2col} equivalent to the variable represented by 
#' \code{repstra2col}. This can be useful if two reproductive status variables
#' have related but unequal units, for example if \code{repstrb2col} refers to
#' one-flowered stems while \code{repstra2col} refers to two-flowered stems.
#' @param fecrel This is a scalar multiplier making the variable represented by
#' \code{fecb2col} equivalent to the variable represented by \code{feca2col}.
#' This can be useful if two fecundity variables have related but unequal units.
#' @param stage2col Optional variable name or column number corresponding to
#' life history stage in occasion \emph{t}.
#' @param stage3col Optional variable name or column number corresponding to
#' life history stage in occasion \emph{t}+1.
#' @param juv2col A variable name or column number that marks individuals in
#' immature stages in occasion \emph{t}. Function \code{historicalize3()}
#' assumes that immature individuals are identified in this variable marked with
#' a number equal to or greater than \code{1}, and that mature individuals are
#' marked as \code{0} or \code{NA}.
#' @param juv3col A variable name or column number that marks individuals in
#' immature stages in occasion \emph{t}+1. Function \code{historicalize3()}
#' assumes that immature individuals are identified in this variable marked with
#' a number equal to or greater than \code{1}, and that mature individuals are
#' marked as \code{0} or \code{NA}.
#' @param stageassign The stageframe object identifying the life history model
#' being operationalized. Note that if \code{stage2col} is provided, then this
#' stageframe is not utilized in stage designation.
#' @param stagesize A variable name or column number describing which size 
#' variable to use in stage estimation. Defaults to \code{NA}, and can also take 
#' \code{sizea}, \code{sizeb}, \code{sizec}, \code{sizeab}, \code{sizebc},
#' \code{sizeac}, \code{sizeabc}, or \code{sizeadded}, depending on
#' which size variable within the input dataset is chosen. Note that the
#' variable(s) chosen should be presented in the order of the primary,
#' secondary, and tertiary variables in the stageframe input with
#' \code{stageassign}. For example, choosing \code{sizeb} assumes that this size
#' variable in the dataset is the primary variable in the stageframe.
#' @param censor A logical variable determining whether the output data should
#' be censored using the variable defined in \code{censorcol}. Defaults to 
#' \code{FALSE}.
#' @param censorcol A variable name or column number corresponding to a censor
#' variable within the dataset, used to distinguish between entries to use and
#' those to discard from analysis, or to designate entries with special issues 
#' that require further attention.
#' @param censorkeep The value of the censoring variable identifying data that
#' should be included in analysis. Defaults to \code{0}, but may take any value
#' including \code{NA}.
#' @param spacing The spacing at which density should be estimated, if density
#' estimation is desired and X and Y coordinates are supplied. Given in the same
#' units as those used in the X and Y coordinates given in \code{xcol} and
#' \code{ycol}. Defaults to \code{NA}.
#' @param NAas0 If TRUE, then all \code{NA} entries for size and fecundity
#' variables will be set to \code{0}. This can help increase the sample size
#' analyzed by  \code{\link{modelsearch}()}, but should only be used when it is
#' clear that this substitution is biologically realistic. Defaults to
#' \code{FALSE}.
#' @param NRasRep If set to \code{TRUE}, then this function will treat
#' non-reproductive but mature individuals as reproductive during stage
#' assignment. This can be useful when a matrix is desired without separation of
#' reproductive and non-reproductive but mature stages of the same size. Only
#' used if \code{stageassign} is set to a valid stageframe. Defaults to
#' \code{FALSE}.
#' @param NOasObs If \code{TRUE}, then will treat individuals that are
#' interpreted as not observed in the dataset as though they were observed
#' during stage assignment. This can be useful when a MPM is desired without
#' separation of observable and unobservable stages. Only used if
#' \code{stageassign} is set to a stageframe. Defaults to \code{FALSE}.
#' @param prebreeding A logical term indicating whether the life history model
#' is pre-breeding. If so, then \code{1} is added to all ages. Defaults to
#' \code{TRUE}.
#' @param age_offset A number to add automatically to all values of age at time
#' \emph{t}. Defaults to \code{0}.
#' @param reduce A logical variable determining whether unused variables and
#' some invariant state variables should be removed from the output dataset.
#' Defaults to \code{TRUE}.
#' @param a2check A logical variable indicating whether to retain all data with
#' living status at occasion \emph{t}. Defaults to \code{FALSE}, in which case
#' data for occasions in which the individual is not alive in time \emph{t} is
#' not retained. This option should be kept \code{FALSE}, except to inspect
#' potential errors in the dataset.
#' @param quiet A logical variable indicating whether to silence warnings.
#' Defaults to \code{FALSE}.
#' 
#' @return If all inputs are properly formatted, then this function will output
#' a historical vertical data frame (class \code{hfvdata}), meaning that the
#' output data frame will have three consecutive years of size and reproductive
#' data per individual per row. This data frame is in standard format for all
#' functions used in \code{lefko3}, and so can be used without further 
#' modification. Note that determination of state in occasions \emph{t}-1 and
#' \emph{t}+1 gives preference to condition in occasion \emph{t} within the
#' input dataset. Conflicts in condition in input datasets that have both
#' occasions \emph{t} and \emph{t}+1 listed per row are resolved by using
#' condition in occasion \emph{t}.
#' 
#' Variables in this data frame include the following:
#' \item{rowid}{Unique identifier for the row of the data frame.}
#' \item{popid}{Unique identifier for the population, if given.}
#' \item{patchid}{Unique identifier for patch within population, if given.}
#' \item{individ}{Unique identifier for the individual.}
#' \item{year2}{Year or time in occasion \emph{t}.}
#' \item{firstseen}{Occasion of first observation.}
#' \item{lastseen}{Occasion of last observation.}
#' \item{obsage}{Observed age in occasion \emph{t}, assuming first observation
#' corresponds to age = 0.}
#' \item{obslifespan}{Observed lifespan, given as \code{lastseen - firstseen + 1}.}
#' \item{xpos1,xpos2,xpos3}{X position in Cartesian space in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively, if provided.}
#' \item{ypos1,ypos2,ypos3}{Y position in Cartesian space in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively, if provided.}
#' \item{sizea1,sizea2,sizea3}{Main size measurement in occasions \emph{t}-1,
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizeb1,sizeb2,sizeb3}{Secondary size measurement in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizec1,sizec2,sizec3}{Tertiary size measurement in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{size1added,size2added,size3added}{Sum of primary, secondary, and
#' tertiary size measurements in occasions \emph{t}-1, \emph{t}, and \emph{t}+1,
#' respectively.}
#' \item{repstra1,repstra2,repstra3}{Main numbers of reproductive structures in
#' occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstrb1,repstrb2,repstrb3}{Secondary numbers of reproductive
#' structures in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstr1added,repstr2added,repstr3added}{Sum of primary and secondary
#' reproductive structures in occasions \emph{t}-1, \emph{t}, and \emph{t}+1,
#' respectively.}
#' \item{feca1,feca2,feca3}{Main numbers of offspring in occasions \emph{t}-1,
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecb1,fecb2, fecb3}{Secondary numbers of offspring in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fec1added,fec2added,fec3added}{Sum of primary and secondary fecundity
#' in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{censor1,censor2,censor3}{Censor status values in occasions \emph{t}-1,
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{juvgiven1,juvgiven2,juvgiven3}{Binomial variable indicating whether
#' individual is juvenile in occasions \emph{t}-1, \emph{t}, and \emph{t}+1.
#' Only given if \code{juvcol} is provided.}
#' \item{obsstatus1,obsstatus2,obsstatus3}{Binomial observation status in
#' occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstatus1,repstatus2,repstatus3}{Binomial reproductive status in
#' occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecstatus1,fecstatus2,fecstatus3}{Binomial offspring production status
#' in occasions \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{matstatus1,matstatus2,matstatus3}{Binomial maturity status in occasions
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{alive1,alive2,alive3}{Binomial status as alive in occasions \emph{t}-1,
#'  \emph{t}, and \emph{t}+1, respectively.}
#' \item{density}{Density of individuals per unit designated in \code{spacing}.
#' Only given if spacing is not \code{NA}.}
#' 
#' @section Notes:
#' Warnings that some individuals occur in state combinations that do not match
#' any stages in the stageframe used to assign stages, and that some individuals
#' match characteristics of several stages in the stageframe, are common when
#' first working with a dataset. Typically, these situations can be identified as
#' \code{NoMatch} entries in \code{stage3}, although such entries may crop up in
#' \code{stage1} and \code{stage2}, as well. In some cases, these warnings will
#' arise with no concurrent \code{NoMatch} entries. These are important warnings
#' and suggest that there is likely a problem with the stageframe. The most
#' common such problems are: 1) stages have significant overlap in
#' characteristics, with the most common being overlapping size bins caused by
#' erroneous definitions of size bin halfwidths; and 2) some individuals exist
#' in states not defined within the stageframe.
#' 
#' In some datasets with unobservable stages, observation status
#' (\code{obsstatus}) might not be inferred properly if a single size variable
#' is used that does not yield sizes greater than 0 in all cases in which
#' individuals were observed. Such situations may arise, for example, in plants
#' when leaf number is the dominant size variable used, but individuals
#' occasionally occur with inflorescences but no leaves. In this instances,
#' it helps to mark related variables as \code{sizeb} and \code{sizec}, because
#' observation status will be interpreted in relation to all 3 size variables.
#' Alternatively, observation status may be input via \code{obs2col} and
#' \code{obs3col} to force computation with given values (although this requires
#' all instances of observation and non-observation to be known and coded ahead
#' of time). Further analysis can then utilize only a single size variable, of
#' the user's choosing. Similar issues can arise in reproductive status
#' (\code{repstatus}).
#' 
#' Juvenile designation should only be used when juveniles fall outside of the
#' size classification scheme used in determining stages. If juveniles are to be
#' size classified along the size spectrum that adults also fall on, then
#' it is best to treat juveniles as mature but not reproductive.
#' 
#' Care should be taken to avoid variables with negative values indicating size,
#' fecundity, or reproductive or observation status. Negative values can be
#' interpreted in different ways, typically reflecting estimation through other
#' algorithms rather than actual measured data. Variables holding negative
#' values can conflict with data management algorithms in ways that are
#' difficult to predict.
#' 
#' Unusual errors (e.g. \code{"Error in pjf..."}) may occur in cases where the
#' variables are improperly passed, or where seemingly numeric variables include
#' text and so get automatically converted to string variables.
#' 
#' Density estimation is performed as a count of individuals alive and within
#' the radius specified in \code{spacing} of the respective individual at some
#' point in time.
#' 
#' @examples
#' data(cypvert)
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
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", 
#'   individcol = "plantid", year2col = "year2", sizea2col = "Inf2.2", 
#'   sizea3col = "Inf2.3", sizeb2col = "Inf.2", sizeb3col = "Inf.3", 
#'   sizec2col = "Veg.2", sizec3col = "Veg.3", repstra2col = "Inf2.2", 
#'   repstra3col = "Inf2.3", repstrb2col = "Inf.2", repstrb3col = "Inf.3", 
#'   feca2col = "Pod.2", feca3col = "Pod.3", repstrrel = 2, 
#'   stageassign = cypframe_raw, stagesize = "sizeadded", censorcol = "censor",
#'   censor = FALSE, NAas0 = TRUE, NRasRep = TRUE, reduce = TRUE)
#'
#' @export
historicalize3 <- function(data, popidcol = 0, patchidcol = 0, individcol, 
  year2col = 0, year3col = 0, xcol = 0, ycol = 0, sizea2col = 0, sizea3col = 0,
  sizeb2col = 0, sizeb3col = 0, sizec2col = 0, sizec3col = 0, repstra2col = 0, 
  repstra3col = 0, repstrb2col = 0, repstrb3col = 0, feca2col = 0, feca3col = 0,
  fecb2col = 0, fecb3col = 0, indcova2col = 0, indcova3col = 0, indcovb2col = 0,
  indcovb3col = 0, indcovc2col = 0, indcovc3col = 0, alive2col = 0,
  alive3col = 0, dead2col = 0, dead3col = 0, obs2col = 0, obs3col = 0,
  nonobs2col = 0, nonobs3col = 0, repstrrel = 1, fecrel = 1, stage2col = 0,
  stage3col = 0, juv2col = 0, juv3col = 0, stageassign = NA, stagesize = NA,
  censor = FALSE, censorcol = 0, censorkeep = 0, spacing = NA, NAas0 = FALSE,
  NRasRep = FALSE, NOasObs = FALSE, prebreeding = TRUE, age_offset = 0,
  reduce = TRUE, a2check = FALSE, quiet = FALSE) {
  
  alive2 <- indataset <- censor1 <- censor2 <- censor3 <- censbool <- NULL
  
  if (!all(is.logical(c(NAas0, NRasRep, NOasObs, reduce, quiet, prebreeding)))) {
    stop("Some logical variables have been assigned non-logical values.",
      call. = FALSE)
  }
  
  if (prebreeding) age_offset <- age_offset + 1
  
  if (is.na(individcol)) {
    stop("Individual ID variable is required.", call. = FALSE)
  }
  
  if (is.na(year2col) & is.na(year3col)) {
    stop("Variable identifying either year2 (occasion t) or year3 (occasion t+1) is required.",
      call. = FALSE)
  }
  
  if (is.character(popidcol)) {
    if (is.element(popidcol, names(data))) {
      true.popidcol <- which(names(data) == popidcol)
      popidcol <- true.popidcol
    } else {
      stop("Please enter popidcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(patchidcol)) {
    if (is.element(patchidcol, names(data))) {
      true.patchidcol <- which(names(data) == patchidcol)
      patchidcol <- true.patchidcol
    } else {
      stop("Please enter patchidcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(individcol)) {
    if (is.element(individcol, names(data))) {
      true.individcol <- which(names(data) == individcol)
      individcol <- true.individcol
    } else {
      stop("Please enter individcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (year2col != 0) {
    if (is.character(year2col)) {
      if (is.element(year2col, names(data))) {
        true.year2col <- which(names(data) == year2col)
        year2col <- true.year2col
      } else {
        stop("Please enter year2col exactly as it appears in the dataset.",
          call. = FALSE)
      }
    }
  } 
  
  if (year3col != 0) {
    if (is.character(year3col)) {
      if (is.element(year3col, names(data))) {
        true.year3col <- which(names(data) == year3col)
        year3col <- true.year3col
      } else {
        stop("Please enter year3col exactly as it appears in the dataset.",
          call. = FALSE)
      }
    }
  }
  
  if (year2col != 0 & year3col == 0) {
    data$year3 <- data[,year2col] + 1
    year3col <- which(names(data) == "year3")
  } else if (year2col == 0 & year3col != 0) {
    data$year2 <- data[,year3col] - 1
    year2col <- which(names(data) == "year2")
  }
  
  if (is.character(xcol)) {
    if (is.element(xcol, names(data))) {
      true.xcol <- which(names(data) == xcol)
      xcol <- true.xcol
    } else {
      stop("Please enter xcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(ycol)) {
    if (is.element(ycol, names(data))) {
      true.ycol <- which(names(data) == ycol)
      ycol <- true.ycol
    } else {
      stop("Please enter ycol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(sizea2col)) {
    if (is.element(sizea2col, names(data))) {
      true.sizea2col <- which(names(data) == sizea2col)
      sizea2col <- true.sizea2col
    } else {
      stop("Please enter sizea2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(sizea3col)) {
    if (is.element(sizea3col, names(data))) {
      true.sizea3col <- which(names(data) == sizea3col)
      sizea3col <- true.sizea3col
    } else {
      stop("Please enter sizea3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(sizeb2col)) {
    if (is.element(sizeb2col, names(data))) {
      true.sizeb2col <- which(names(data) == sizeb2col)
      sizeb2col <- true.sizeb2col
    } else {
      stop("Please enter sizeb2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(sizeb3col)) {
    if (is.element(sizeb3col, names(data))) {
      true.sizeb3col <- which(names(data) == sizeb3col)
      sizeb3col <- true.sizeb3col
    } else {
      stop("Please enter sizeb3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(sizec2col)) {
    if (is.element(sizec2col, names(data))) {
      true.sizec2col <- which(names(data) == sizec2col)
      sizec2col <- true.sizec2col
    } else {
      stop("Please enter sizec2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(sizec3col)) {
    if (is.element(sizec3col, names(data))) {
      true.sizec3col <- which(names(data) == sizec3col)
      sizec3col <- true.sizec3col
    } else {
      stop("Please enter sizec3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(repstra2col)) {
    if (is.element(repstra2col, names(data))) {
      true.repstra2col <- which(names(data) == repstra2col)
      repstra2col <- true.repstra2col
    } else {
      stop("Please enter repstra2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(repstra3col)) {
    if (is.element(repstra3col, names(data))) {
      true.repstra3col <- which(names(data) == repstra3col)
      repstra3col <- true.repstra3col
    } else {
      stop("Please enter repstra3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(repstrb2col)) {
    if (is.element(repstrb2col, names(data))) {
      true.repstrb2col <- which(names(data) == repstrb2col)
      repstrb2col <- true.repstrb2col
    } else {
      stop("Please enter repstrb2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(repstrb3col)) {
    if (is.element(repstrb3col, names(data))) {
      true.repstrb3col <- which(names(data) == repstrb3col)
      repstrb3col <- true.repstrb3col
    } else {
      stop("Please enter repstrb3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(feca2col)) {
    if (is.element(feca2col, names(data))) {
      true.feca2col <- which(names(data) == feca2col)
      feca2col <- true.feca2col
    } else {
      stop("Please enter feca2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(feca3col)) {
    if (is.element(feca3col, names(data))) {
      true.feca3col <- which(names(data) == feca3col)
      feca3col <- true.feca3col
    } else {
      stop("Please enter feca3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(fecb2col)) {
    if (is.element(fecb2col, names(data))) {
      true.fecb2col <- which(names(data) == fecb2col)
      fecb2col <- true.fecb2col
    } else {
      stop("Please enter fecb2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(fecb3col)) {
    if (is.element(fecb3col, names(data))) {
      true.fecb3col <- which(names(data) == fecb3col)
      fecb3col <- true.fecb3col
    } else {
      stop("Please enter fecb3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(indcova2col)) {
    if (is.element(indcova2col, names(data))) {
      true.indcova2col <- which(names(data) == indcova2col)
      indcova2col <- true.indcova2col
    } else {
      stop("Please enter indcova2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(indcova3col)) {
    if (is.element(indcova3col, names(data))) {
      true.indcova3col <- which(names(data) == indcova3col)
      indcova3col <- true.indcova3col
    } else {
      stop("Please enter indcova3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(indcovb2col)) {
    if (is.element(indcovb2col, names(data))) {
      true.indcovb2col <- which(names(data) == indcovb2col)
      indcovb2col <- true.indcovb2col
    } else {
      stop("Please enter indcovb2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(indcovb3col)) {
    if (is.element(indcovb3col, names(data))) {
      true.indcovb3col <- which(names(data) == indcovb3col)
      indcovb3col <- true.indcovb3col
    } else {
      stop("Please enter indcovb3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(indcovc2col)) {
    if (is.element(indcovc2col, names(data))) {
      true.indcovc2col <- which(names(data) == indcovc2col)
      indcovc2col <- true.indcovc2col
    } else {
      stop("Please enter indcovc2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(indcovc3col)) {
    if (is.element(indcovc3col, names(data))) {
      true.indcovc3col <- which(names(data) == indcovc3col)
      indcovc3col <- true.indcovc3col
    } else {
      stop("Please enter indcovc3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(alive2col)) {
    if (is.element(alive2col, names(data))) {
      true.alive2col <- which(names(data) == alive2col)
      alive2col <- true.alive2col
    } else {
      stop("Please enter alive2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(alive3col)) {
    if (is.element(alive3col, names(data))) {
      true.alive3col <- which(names(data) == alive3col)
      alive3col <- true.alive3col
    } else {
      stop("Please enter alive3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(dead2col)) {
    if (is.element(dead2col, names(data))) {
      true.dead2col <- which(names(data) == dead2col)
      dead2col <- true.dead2col
    } else {
      stop("Please enter dead2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(dead3col)) {
    if (is.element(dead3col, names(data))) {
      true.dead3col <- which(names(data) == dead3col)
      dead3col <- true.dead3col
    } else {
      stop("Please enter dead3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(obs2col)) {
    if (is.element(obs2col, names(data))) {
      true.obs2col <- which(names(data) == obs2col)
      obs2col <- true.obs2col
    } else {
      stop("Please enter obs2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(obs3col)) {
    if (is.element(obs3col, names(data))) {
      true.obs3col <- which(names(data) == obs3col)
      obs3col <- true.obs3col
    } else {
      stop("Please enter obs3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(nonobs2col)) {
    if (is.element(nonobs2col, names(data))) {
      true.nonobs2col <- which(names(data) == nonobs2col)
      nonobs2col <- true.nonobs2col
    } else {
      stop("Please enter nonobs2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(nonobs3col)) {
    if (is.element(nonobs3col, names(data))) {
      true.nonobs3col <- which(names(data) == nonobs3col)
      nonobs3col <- true.nonobs3col
    } else {
      stop("Please enter nonobs3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(stage2col)) {
    if (is.element(stage2col, names(data))) {
      true.stage2col <- which(names(data) == stage2col)
      stage2col <- true.stage2col
    } else {
      stop("Please enter stage2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(stage3col)) {
    if (is.element(stage3col, names(data))) {
      true.stage3col <- which(names(data) == stage3col)
      stage3col <- true.stage3col
    } else {
      stop("Please enter stage3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(juv2col)) {
    if (is.element(juv2col, names(data))) {
      true.juv2col <- which(names(data) == juv2col)
      juv2col <- true.juv2col
    } else {
      stop("Please enter juv2col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(juv3col)) {
    if (is.element(juv3col, names(data))) {
      true.juv3col <- which(names(data) == juv3col)
      juv3col <- true.juv3col
    } else {
      stop("Please enter juv3col exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (is.character(censorcol)) {
    if (is.element(censorcol, names(data))) {
      true.censorcol <- which(names(data) == censorcol)
      censorcol <- true.censorcol
    } else {
      stop("Please enter censorcol exactly as it appears in the dataset.",
        call. = FALSE)
    }
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.",
        call. = FALSE)
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.",
        call. = FALSE)
    }
  }
  
  if (!all(is.na(stageassign))) {
    
    stagesize <- tolower(stagesize)
    stagesize_options <- c("sizea", "sizeb", "sizec", "sizeadded", "sizeab",
      "sizeac", "sizebc", "sizeabc")
    
    if(!is.element(stagesize, stagesize_options)) {
      stop("The stagesize option must equal 'NA', 'sizea', 'sizeb', 'sizec', 'sizeab',
        'sizeac', 'sizebc', 'sizeabc', or 'sizeadded'.", call. = FALSE)
    }
    
    stassign <- TRUE 
    
    if (stagesize == "sizeabc") {
      stagesizecol <- 8
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size == 0),
        intersect(which(stageassign$size_b == 0), which(stageassign$size_c == 0))))
    } else if (stagesize == "sizebc") {
      stagesizecol <- 7
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size_b == 0),
        which(stageassign$size_c == 0)))
    } else if (stagesize == "sizeac") {
      stagesizecol <- 6
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size == 0),
        which(stageassign$size_c == 0)))
    } else if (stagesize == "sizeab") {
      stagesizecol <- 5
      repcheck1 <- intersect(which(stageassign$repstatus == 1), intersect(which(stageassign$size == 0),
        which(stageassign$size_b == 0)))
    } else if (stagesize == "sizeadded") {
      stagesizecol <- 4
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size == 0))
    } else if (stagesize == "sizec") {
      stagesizecol <- 3
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size_c == 0))
    } else if (stagesize == "sizeb") {
      stagesizecol <- 2
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size_b == 0))
    } else {
      stagesizecol <- 1
      repcheck1 <- intersect(which(stageassign$repstatus == 1), which(stageassign$size == 0))
    }
    
    repcheck <- intersect(repcheck1, which(stageassign$obsstatus == 1))
    
    if (length(repcheck) > 0) {
      message("Stageframe indicates the presence of an observeable, reproductive stage with a size of 0.")
      RepasObs <- TRUE
    } else if (length(repcheck1) > 0) {
      message("Stageframe indicates the presence of an unobserveable, reproductive stage with a size of 0.")
    }
    
  } else {
    stassign <- FALSE
    
    stageassign <- sf_skeleton(1)
    stageassign$stage[1] <- NA
    stageassign$size[1] <- 0
    stageassign$obsstatus[1] <- 0
    stageassign$indataset[1] <- 0
    stageassign$binhalfwidth_raw[1] <- 0
    stageassign$sizebin_center[1] <- 0
    stageassign$sizebin_width[1] <- 0
    stageassign$comments[1] <- NA
    
    stagesizecol <- 0
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  if (censor) {
    if (!is.element(censorkeep, data[,censorcol])) {
      stop("Please enter a valid value for censorkeep. This value should occur
        in the censor variable within the dataset.", call. = FALSE)
    }
  }
  
  if (is.na(censorkeep) & censor) {
    censbool <- TRUE
    censorkeep <- 0
  } else {
    censbool <- FALSE
  }
  
  popdata <- .jpf(data, stageassign, (popidcol - 1), (patchidcol - 1),
    (individcol - 1), (year2col - 1), (year3col - 1), (xcol - 1), (ycol - 1),
    (juv2col - 1), (juv3col - 1), (sizea2col - 1), (sizea3col - 1),
    (sizeb2col - 1), (sizeb3col - 1), (sizec2col - 1), (sizec3col - 1),
    (repstra2col - 1), (repstra3col - 1), (repstrb2col - 1), (repstrb3col - 1),
    (feca2col - 1), (feca3col - 1), (fecb2col - 1), (fecb3col - 1),
    (indcova2col - 1), (indcova3col - 1), (indcovb2col - 1), (indcovb3col - 1),
    (indcovc2col - 1), (indcovc3col - 1), (alive2col - 1), (alive3col - 1),
    (dead2col - 1), (dead3col - 1), (obs2col - 1), (obs3col - 1),
    (nonobs2col - 1), (nonobs3col - 1), repstrrel, fecrel, (stage2col - 1),
    (stage3col - 1), (censorcol - 1), NAas0, NRasRep, NOasObs, stassign,
    stagesizecol, censorkeep, censbool, a2check, reduce, quiet)
  
  if ("censor2" %in% colnames(popdata) & censor) {
    popdata <- subset(popdata, censor1 == censorkeep & censor2 == censorkeep)
    popdata <- subset(popdata, censor3 == censorkeep)
  }
  
  if (!is.na(spacing) & "xpos2" %in% colnames(popdata)) {
    popdata$density <- .density3(popdata, which(names(popdata) == "xpos2"), 
      which(names(popdata) == "ypos2"), which(names(popdata) == "year2"), spacing)
  }
  
  popdata$obsage <- popdata$obsage + age_offset
  
  return(popdata)
}

#' Import MPM from COMPADRE or COMADRE Database
#' 
#' Function \code{.import_Com()} imports matrices from the COMPADRE and COMADRE
#' databases, given knowledge of the corresponding \code{MatrixID}, and builds
#' them into a \code{lefkoMat} object. Users are encourage to explore package
#' \code{Rcompadre} for further details of the use of these databases.
#' 
#' @name .import_Com
#' 
#' @param matrix_id The values of \code{MatrixID} from the used database
#' corresponding to the matrices to import.
#' @param database The holding the COMPADRE or COMADRE database, as in the
#' global environment.
#' @param add_FC A logical value indicating whether to sum the \code{matF} and
#' \code{matC} matrices to produce the \code{F} matrix. If \code{FALSE}, then
#' only uses the \code{matF} matrix. Defaults to \code{TRUE}.
#' 
#' @return A \code{lefkoMat} object with the following elements:
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and occasions. All matrices output in the matrix class.}
#' \item{U}{A list of survival transition matrices sorted as in A. All matrices
#' output in the matrix class.}
#' \item{F}{A list of fecundity matrices sorted as in A. All matrices output in
#' the matrix class.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes status
#' as an entry stage through reproduction.}
#' \item{hstages}{A single value of \code{NA}.}
#' \item{agestages}{A single value of \code{NA}.}
#' \item{labels}{A data frame giving the population, patch, and year of each
#' matrix in order. Taken from the \code{MatrixPopulation},
#' \code{MatrixTreatment}, and \code{MatrixStartYear} variables in the database
#' metadata.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{When output from COMPADRE or COMADRE, set a a vector of two NA
#' values.}
#' 
#' @keywords internal
#' @noRd
.import_Com <- function(matrix_id, database, add_FC = TRUE) {
  output <- core_indices <- meta_data <- stage_lists <- matrices <- NULL
  A_mats <- U_mats <- F_mats <- NULL
  db_format <- num_mats <- stages_num <- 0
  
  if (!is.null(database)) {
    if (is(database, "CompadreDB")) {
      core_indices <- which(is.element(database@data$MatrixID, matrix_id))
      meta_data <- as.data.frame(database@data[core_indices,-1])
      format <- 1
    } else if (is.list(database)) {
      core_indices <- which(is.element(database$metadata$MatrixID, matrix_id))
      meta_data <- database$metadata[core_indices,]
      format <- 2
    } else {
      stop("Object database not recognized.", call. = FALSE)
    }
    
    if (length(core_indices) == 0) stop("Entered matrix_id cannot be found.", call. = FALSE)
    num_mats <- length(core_indices)
    num_mats_vec <- as.matrix(c(1:num_mats))
    
    if (format == 1) {
      stage_lists <- apply(num_mats_vec, 1, function(X) {
        return(database@data$mat[[core_indices[X]]]@matrixClass)
      })
      
      matrices <- apply(num_mats_vec, 1, function(X) {
        return(database@data$mat[[core_indices[X]]])
      })
      
      A_mats <- lapply(matrices, function(X) {
        return(X@matA)
      })
      
      U_mats <- lapply(matrices, function(X) {
        return(X@matU)
      })
      
      F_mats <- lapply(matrices, function(X) {
        if (add_FC) {
          return(X@matF + X@matC)
        } else {
          return(X@matF)
        }
      })
    } else if (format == 2) {
      stage_lists <- apply(num_mats_vec, 1, function(X) {
        return(database$matrixClass[[core_indices[X]]])
      })
      
      matrices <- apply(num_mats_vec, 1, function(X) {
        return(database$mat[[core_indices[X]]])
      })
      
      A_mats <- lapply(matrices, function(X) {
        return(X$matA)
      })
      
      U_mats <- lapply(matrices, function(X) {
        return(X$matU)
      })
      
      F_mats <- lapply(matrices, function(X) {
        if (add_FC) {
          return(X$matF + X$matC)
        } else {
          return(X$matF)
        }
      })
    }
    
    check_rows <- apply(num_mats_vec, 1, function(X) {
      return(dim(A_mats[[X]])[1])
    })
    if (length(unique(check_rows)) > 1) {
      stop("Input matrices appear to be of unequal dimension.", call. = FALSE)
    }
    stages_num <- check_rows[1]
    
    found_elems_U <- apply(num_mats_vec, 1, function(X) {
      return(length(which(U_mats[[X]] > 0.0)))
    })
    found_elems_F <- apply(num_mats_vec, 1, function(X) {
      return(length(which(F_mats[[X]] > 0.0)))
    })
    matrix_qc <- c(sum(found_elems_U), sum(found_elems_F), num_mats)
    
    reprod_stages_sum <- t(apply(num_mats_vec, 1, function(X) {
      workup1 <- colSums(F_mats[[X]])
    }))
    reprod_stages <- colSums(reprod_stages_sum)
    reprod_stages[which(reprod_stages > 0)] <- 1
    
    entry_stages_sum <- t(apply(num_mats_vec, 1, function(X) {
      workup1 <- rowSums(F_mats[[X]])
    }))
    entry_stages <- colSums(entry_stages_sum)
    entry_stages[which(entry_stages > 0)] <- 1
    
    mature_stages <- rep(0, stages_num)
    mature_stages[which(entry_stages == 0)] <- 1
    mature_stages[which(reprod_stages > 0)] <- 1
    
    pops <- meta_data$MatrixPopulation
    patches <- meta_data$MatrixTreatment
    years <- meta_data$MatrixStartYear
    labels <- data.frame(pop = pops, patch = patches, year2 = years)
    
    new_sf <- sf_skeleton(stages_num)
    if (!is.element("stage_id", names(new_sf$ahstages))) {
      no_cols <- length(new_sf)
      new_sf$stage_id <- c(1:length(new_sf$stage))
      new_sf <- new_sf[, c((no_cols+1), 1:no_cols)]
    }
    new_sf$repstatus <- reprod_stages
    new_sf$immstatus <- entry_stages
    new_sf$matstatus <- mature_stages
    new_sf$comments <- stage_lists[[1]]$MatrixClassAuthor
    
    dataqc <- c(NA, NA)
    
    output <- list(A = A_mats, U = U_mats, F = F_mats, ahstages = new_sf,
      hstages = as.data.frame(NA), agestages = as.data.frame(NA),
      labels = labels, matrixqc = matrix_qc, dataqc = dataqc)
    class(output) <- "lefkoMat"
  } else {
    message("Please load either the COMPADRE or COMADRE database.")
  }
  
  return(output)
}

#' Create lefkoMat Object from Given Input Matrices
#' 
#' Function \code{.import_mats()} creates lefkoMat objects from supplied
#' matrices and extra information.
#' 
#' @name .import_mats
#' 
#' @param mats A list of A matrices.
#' @param stageframe A stageframe describing all stages utilized.
#' @param hstages A data frame outlining the order of historical stages, if
#' matrices provided in \code{mats} are historical. Defaults to NA.
#' @param agestages A data frame outlining the order of ahistorical age-stages,
#' if age-by-stage matrices are provided.
#' @param historical A logical value indicating whether input matrices are
#' historical or not. Defaults to FALSE.
#' @param agebystage A logical value indicating whether input matrices are
#' ahistorical age-by-stage matrices. If TRUE, then object \code{agestages} is
#' required. Defaults to FALSE.
#' @param UFdecomp A logical value indicating whether U and F matrices should be
#' inferred. Defaults to TRUE.
#' @param entrystage The stage or stages produced by reproductive individuals.
#' Used to determine which transitions are reproductive for U-F decomposition.
#' Defaults to \code{1}, which corresponds to the first stage in the stageframe.
#' @param poporder The order of populations in the list supplied in object
#' \code{mats}. Defaults to 1.
#' @param patchorder The order of patches in the list supplied in object
#' \code{mats}. Defaults to 1.
#' @param yearorder The order of monitoring occasions in the list supplied in
#' object \code{mats}. Defaults to NA, which leads to each matrix within each
#' population-patch combination being a different monitoring occasion.
#' @param sparse_output A logical value indicating whether to output matrices in
#' sparse format. Defaults to \code{FALSE}, in which case all matrices are
#' output in standard matrix format.
#' 
#' @return A \code{lefkoMat} object incorporating the matrices input in object
#' \code{mats} as object \code{A}, their U and F decompositions in objects
#' \code{U} and \code{F} (if requested), the provided stageframe as object
#' \code{ahstages}, the order of historical stages as object \code{hstages} (if
#' \code{historical = TRUE}), the order of matrices as object \code{labels},
#' a short quality control section used by the \code{\link{summary.lefkoMat}()}
#' function (\code{matrixqc}), and an empty \code{dataqc} object.
#' 
#' @keywords internal
#' @noRd
.import_mats <- function(mats, stageframe, hstages = NA, agestages = NA,
  historical = FALSE, agebystage = FALSE, UFdecomp = TRUE, entrystage = 1,
  poporder = 1, patchorder = 1, yearorder = NA, sparse_output = FALSE) {
  
  F_indices <- NULL
  
  if (!is.list(mats)) {
    stop("Object mats must be an object of class list.", call. = FALSE)
  }
  if (!(is.matrix(mats[[1]]) | is(mats[[1]], "dgCMatrix"))) {
    stop("Object mats must be a list composed of objects of class matrix or dgCMatrix.",
      call. = FALSE)
  }
  mat_length <- length(mats)
  
  if (!is(stageframe, "stageframe")) {
    warning("Object stageframe is not of class stageframe.", call. = FALSE)
  }
  if(all(is.na(agestages)) & agebystage) {
    stop("Age-by-stage matrix inputs require a valid agestages input.", 
      call. = FALSE)
  }
  
  numstages <- dim(stageframe)[1]
  
  dimtest <- unique(as.vector(apply(as.matrix(c(1:mat_length)), 1, function(X) {
    dim(as.matrix(mats[[X]]))
  })))
  if (length(dimtest) != 1) {
    stop("Supplied matrices must be square matrices of equal size.",
      call. = FALSE)
  }
  if (!historical & !agebystage) {
    if (dimtest != numstages) {
      stop("Ahistorical matrices must have the same number of rows and columns as the number of stages in the stageframe",
        call. = FALSE)
    }
  } else if (agebystage) {
    if (dimtest != dim(agestages)[1]) {
      stop("Age-by-stage matrices must have the same number of rows and columns as the number of rows in the agestages data frame.",
        call. = FALSE)
    }
  }
  
  if (historical & all(is.na(hstages))) { # Automatic hstages creation
    stage_vec <- stageframe$stage
    stageid_vec <- c(1:numstages)
    hstages <- cbind.data.frame(stage_id_2 = rep(stageid_vec, numstages),
      stage_id_1 = rep(stageid_vec, each = numstages),
      stage_2 = rep(stage_vec, numstages),
      stage_1 = rep(stage_vec, each = numstages))
    
    if (length(hstages$stage_2) != unique(as.vector(dimtest))) {
      stop("If hstages is not provided, then input historical matrices must have rows and columns equal to the number of stages in the stageframe squared.",
        call. = FALSE)
    }
  }
  
  new_frame <- .sf_reassess(stageframe, supplement = NULL, overwrite = NULL,
    repmatrix = NULL, agemat = agebystage, historical = historical, format = 1)
  stageframe <- new_frame$stageframe
  stageframe <- stageframe[-(dim(stageframe)[1]),]
  # Entry stage issue
  
  true_poporder <- poporder
  true_patchorder <- patchorder
  true_yearorder <- yearorder
  
  if (length(poporder) == 1) {
    if (!is.numeric(poporder)) poporder <- 1
    true_poporder <- rep(poporder, mat_length)
  }
  if (length(patchorder) == 1) {
    if (!is.numeric(patchorder)) patchorder <- 1
    true_patchorder <- rep(patchorder, mat_length)
  }
  if (all(is.na(yearorder))) {
    poppatch_combos <- apply(as.matrix(c(1:mat_length)), 1, function(X) {
      paste(true_poporder[X], true_patchorder[X])
      }
    )
    unique_poppatches <- unique(poppatch_combos)
    total_poppatches <- length(unique_poppatches)
    
    true_yearorder <- apply(as.matrix(c(1:total_poppatches)), 1, function(X) {
      total_years <- length(which(poppatch_combos == unique_poppatches[X]))
      return(c(1:total_years))
    })
  }
  
  if (length(true_poporder) != length(true_patchorder)) {
    stop("Objects poporder and patchorder are not of equal length.", call. = FALSE)
  }
  if (length(true_poporder) != length(true_yearorder)) {
    stop("Objects poporder and yearorder are not of equal length.", call. = FALSE)
  }
  
  labels <- cbind.data.frame(pop = true_poporder, patch = true_patchorder,
    year2 = true_yearorder)
  
  if (!any(is.numeric(entrystage)) & UFdecomp) {
    if (all(is.na(entrystage))) {
      warning("No entry stage provided, so assuming that first stage is entry stage.", call. = FALSE)
      entrystage <- 1
    } else if (all(is.element(entrystage, stageframe$stage))) {
      total_entries <- length(entrystage)
      
      entrystage_proxy <- apply(as.matrix(c(1:total_entries)), 1, function(X) {
        which(stageframe$stage == entrystage[X])}
      )
      entrystage <- entrystage_proxy
    } else {
      stop("Unable to interpret entry stage designations.", call. = FALSE)
    }
  }
  
  if (entrystage > length(stageframe$entrystage)) {
    stop("Chosen entry stage does not exist.", call. = FALSE)
  }
  new_entries <- rep(0, length(stageframe$entrystage))
  new_entries[entrystage] <- 1
  stageframe$entrystage <- new_entries
  
  matrixqc <- c(NA, NA, mat_length)
  
  if (UFdecomp) {
    rep_from <- which(stageframe$repstatus == 1)
    
    # Indexing
    if (!historical & !agebystage) { #Ahistorical
      for (i in c(1:length(rep_from))) {
        for (j in c(1:length(entrystage))) {
          if (i == 1 & j == 1) {
            F_indices <- ((numstages * (rep_from[i] - 1)) + entrystage[j])
          } else {
            F_indices <- append(F_indices, ((numstages * (rep_from[i] - 1)) + entrystage[j]))
          }
        }
      }
    } else if (historical & !agebystage) { #Historical
      h_numstages <- length(hstages$stage_2)
      
      for (time2o1 in c(1:h_numstages)) {
        for (time32n in c(1:h_numstages)) {
          if (hstages$stage_id_2[time2o1] == hstages$stage_id_1[time32n]) {
            if (is.element(hstages$stage_id_2[time2o1], rep_from) & 
              is.element(hstages$stage_id_2[time32n], entrystage)) {
              
              if (length(F_indices) == 0) {
                F_indices <- ((time2o1 - 1) * h_numstages) + time32n
              } else {
                F_indices <- append(F_indices, (((time2o1 - 1) * h_numstages) + time32n))
              }
            }
          }
        }
      }
    } else if (agebystage & !historical) { #Age-by-stage ahistorical
      age_numstages <- length(agestages$stage_id)
      
      for (fromas in c(1:age_numstages)) {
        for (toas in c(1:age_numstages)) {
          if (is.element(agestages$stage_id[fromas], rep_from) &
            is.element(agestages$stage_id[toas], entrystage)) {
            
            if (length(F_indices) == 0) {
              F_indices <- ((fromas - 1) * age_numstages) + toas
            } else {
              F_indices <- append(F_indices, (((fromas - 1) * age_numstages) + toas))
            }
          }
        }
      }
    } else {
      stop("Requested format unrecognized. Cannot proceed.", call. = FALSE)
    }
    
    U <- lapply(mats, function(X) {
      newmat <- X
      newmat[F_indices] <- 0
      
      if (!sparse_output) {
        newmat <- as.matrix(newmat)
      } else {
        newmat <- as(newmat, "dgCMatrix")
      }
      return(newmat)
    })
    
    total_mat_elems <- length(as.matrix(mats[[1]]))
    nonFindices <- c(1:total_mat_elems)[!is.element(c(1:total_mat_elems), F_indices)]
    F <- lapply(mats, function(X) {
      newmat <- X
      newmat[nonFindices] <- 0
      
      if (!sparse_output) {
        newmat <- as.matrix(newmat)
      } else {
        newmat <- as(newmat, "dgCMatrix")
      }
      return(newmat)
    })
    
    for (i in c(1:mat_length)) {
      if (!sparse_output) {
        mats[[i]] <- as.matrix(mats[[i]])
      } else {
        mats[[i]] <- as(mats[[i]], "dgCMatrix")
      }
    }
    
    Utrans <- sum(apply(as.matrix(c(1:mat_length)), 1, function(X) {
      length(which(as.matrix(U[[X]]) > 0))
    }))
    
    Ftrans <- sum(apply(as.matrix(c(1:mat_length)), 1, function(X) {
      length(which(as.matrix(F[[X]]) > 0))
    }))
    
    matrixqc[1] <- Utrans
    matrixqc[2] <- Ftrans
    
  } else {
    U <- NA
    F <- NA
    
    Atrans <- sum(apply(as.matrix(c(1:mat_length)), 1, function(X) {
      length(which(as.matrix(mats[[X]]) > 0))
    }))
    
    matrixqc[1] <- Atrans
  }
  
  for (i in c(1:mat_length)) {
    if (!sparse_output) {
      mats[[i]] <- as.matrix(mats[[i]])
    } else {
      mats[[i]] <- as(mats[[i]], "dgCMatrix")
    }
  }
  
  if (is.element("size", names(stageframe))) {
    names(stageframe)[which(names(stageframe) == "size")] <- "original_size"
  }
  
  if (!is.element("entrystage", names(stageframe))) {
    stageframe$entrystage <- 0
    stageframe$entrystage[entrystage] <- 1
  }
  
  dataqc <- c(NA, NA)
  
  output <- list(A = mats, U = U, F = F, hstages = hstages, agestages = NA,
    ahstages = cbind.data.frame(stage_id = c(1:numstages), stageframe),
    labels = labels, matrixqc = matrixqc, dataqc = dataqc)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create lefkoMat Object from Given Input Matrices or an MPM Database
#' 
#' Function \code{create_lM()} creates lefkoMat objects from supplied matrices
#' and extra information, or from a supplied MPM database such as COMPADRE or
#' COMADRE.
#' 
#' @name create_lM
#' 
#' @param mats A list of A matrices, or, if importing from a matrix database
#' such as COMPADRE or COMADRE, then the object holding the database.
#' @param stageframe A stageframe describing all stages utilized.
#' @param hstages A data frame outlining the order of historical stages, if
#' matrices provided in \code{mats} are historical. Defaults to NA.
#' @param agestages A data frame outlining the order of ahistorical age-stages,
#' if age-by-stage matrices are provided.
#' @param historical A logical value indicating whether input matrices are
#' historical or not. Defaults to FALSE.
#' @param agebystage A logical value indicating whether input matrices are
#' ahistorical age-by-stage matrices. If TRUE, then object \code{agestages} is
#' required. Defaults to FALSE.
#' @param UFdecomp A logical value indicating whether U and F matrices should be
#' inferred. Defaults to TRUE.
#' @param entrystage The stage or stages produced by reproductive individuals.
#' Used to determine which transitions are reproductive for U-F decomposition.
#' Defaults to \code{1}, which corresponds to the first stage in the stageframe.
#' @param poporder The order of populations in the list supplied in object
#' \code{mats}. Defaults to 1.
#' @param patchorder The order of patches in the list supplied in object
#' \code{mats}. Defaults to 1.
#' @param yearorder The order of monitoring occasions in the list supplied in
#' object \code{mats}. Defaults to NA, which leads to each matrix within each
#' population-patch combination being a different monitoring occasion.
#' @param matrix_id The values of \code{MatrixID} from the used database
#' corresponding to the matrices to import, if importing from a database. Not
#' used if importing a list of matrices.
#' @param add_FC A logical value indicating whether to sum the \code{matF} and
#' \code{matC} matrices to produce the \code{F} matrix. If \code{FALSE}, then
#' only uses the \code{matF} matrix. Only used if importing from the COMPADRE or
#' COMADRE database. Defaults to \code{TRUE}.
#' @param sparse_output A logical value indicating whether to output matrices in
#' sparse format. Defaults to \code{FALSE}, in which case all matrices are
#' output in standard matrix format. Does not apply to matrices imported from
#' COMPADRE or COMADRE, which are always in standard format.
#' 
#' @return A \code{lefkoMat} object incorporating the matrices input in object
#' \code{mats} as object \code{A}, their U and F decompositions in objects
#' \code{U} and \code{F} (if requested), the provided stageframe as object
#' \code{ahstages}, the order of historical stages as object \code{hstages} (if
#' \code{historical = TRUE}), the order of matrices as object \code{labels}, and
#' a short quality control section used by the \code{\link{summary.lefkoMat}()}
#' function.
#' 
#' @section Notes for importing lists of matrices:
#' Lists may be composed of a mix of matrices in standard R format (i.e. created
#' via the \code{matrix()} function), and in \code{dgCMatrix} sparse format
#' (i.e. created via the \code{Matrix::Matrix()} function with
#' \code{sparse = TRUE}.) All matrices will be forced to one or the other,
#' depending on the value given for the \code{sparse_output} argument.
#' 
#' U and F decomposition assumes that elements holding fecundity values are
#' to be interpreted solely as fecundity rates. Users wishing to split these
#' elements between fecundity and survival should do so manually after running
#' this function.
#' 
#' Age-by-stage MPMs require an \code{agestages} data frame outlining the order
#' of age-stages. This data frame has 3 variables: \code{stage_id}, which is the
#' number of the stage as labelled by the equivalently named variable in the
#' \code{stageframe}; \code{stage}, which is the official name of the stage as
#' given in the equivalently named variable in the \code{stageframe}; and
#' \code{age}, which of course gives the age associated with the stage at that
#' time. The number of rows must be equal to the number of rows and columns of
#' each entered matrix.
#' 
#' Users may edit the \code{dataqc} object, setting the first \code{NA} to the
#' number of individuals sampled, and the second \code{NA} to the number of
#' rows in a vertical version of the demographic dataset. This is not required,
#' however.
#' 
#' @section Notes for importing from COMPADRE or COMADRE:
#' For this function to operate, users must have either the COMPADRE database
#' or the COMADRE database loaded into the global environment. Note that the
#' sample databases supplied within package \code{Rcompadre} will not work with
#' this function.
#' 
#' This function does not and cannot replace the wonderful tools offered to
#' explore the COMPADRE and COMADRE packages. Please see package
#' \code{Rcompadre} to use those tools. Note that function \code{import_Com()}
#' has no relationship to the \code{Rcompadre} development team.
#' 
#' Function \code{import_Com()} requires that the dimensions of all matrices
#' imported into a single \code{lefkoMat} object be equal.
#' 
#' The reproductive and maturity status of each stage is determined by patterns
#' assessed within the \code{F} matrices. Users should check that these values
#' make sense.
#' 
#' Stage names may be edited manually afterward.
#' 
#' Users may edit the \code{dataqc} object, setting the first \code{NA} to the
#' number of individuals sampled, and the second \code{NA} to the number of
#' rows in a vertical version of the demographic dataset. This is not required,
#' however.
#' 
#' @seealso \code{\link{add_lM}()}
#' @seealso \code{\link{delete_lM}()}
#' @seealso \code{\link{subset_lM}()}
#' 
#' @examples
#' # These matrices are of 9 populations of the plant species Anthyllis
#' # vulneraria, and were originally published in Davison et al. (2010) Journal
#' # of Ecology 98:255-267 (doi: 10.1111/j.1365-2745.2009.01611.x).
#' 
#' sizevector <- c(1, 1, 2, 3) # These sizes are not from the original paper
#' stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
#' repvector <- c(0, 0, 1, 1)
#' obsvector <- c(1, 1, 1, 1)
#' matvector <- c(0, 1, 1, 1)
#' immvector <- c(1, 0, 0, 0)
#' propvector <- c(0, 0, 0, 0)
#' indataset <- c(1, 1, 1, 1)
#' binvec <- c(0.5, 0.5, 0.5, 0.5)
#' 
#' anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' # POPN C 2003-2004
#' XC3 <- matrix(c(0, 0, 1.74, 1.74,
#' 0.208333333, 0, 0, 0.057142857,
#' 0.041666667, 0.076923077, 0, 0,
#' 0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
#' 
#' # 2004-2005
#' XC4 <- matrix(c(0, 0, 0.3, 0.6,
#' 0.32183908, 0.142857143, 0, 0,
#' 0.16091954, 0.285714286, 0, 0,
#' 0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
#' 
#' # 2005-2006
#' XC5 <- matrix(c(0, 0, 0.50625, 0.675,
#' 0, 0, 0, 0.035714286,
#' 0.1, 0.068965517, 0.0625, 0.107142857,
#' 0.3, 0.137931034, 0, 0.071428571), 4, 4, byrow = TRUE)
#' 
#' # POPN E 2003-2004
#' XE3 <- matrix(c(0, 0, 2.44, 6.569230769,
#' 0.196428571, 0, 0, 0,
#' 0.125, 0.5, 0, 0,
#' 0.160714286, 0.5, 0.133333333, 0.076923077), 4, 4, byrow = TRUE)
#' 
#' XE4 <- matrix(c(0, 0, 0.45, 0.646153846,
#' 0.06557377, 0.090909091, 0.125, 0,
#' 0.032786885, 0, 0.125, 0.076923077,
#' 0.049180328, 0, 0.125, 0.230769231), 4, 4, byrow = TRUE)
#' 
#' XE5 <- matrix(c(0, 0, 2.85, 3.99,
#' 0.083333333, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.416666667, 0.1, 0, 0.1), 4, 4, byrow = TRUE)
#' 
#' mats_list <- list(XC3, XC4, XC5, XE3, XE4, XE5)
#' yr_ord <- c(1, 2, 3, 1, 2, 3)
#' pch_ord <- c(1, 1, 1, 2, 2, 2)
#' 
#' anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA,
#'   historical = FALSE, poporder = 1, patchorder = pch_ord, yearorder = yr_ord)
#' 
#' # A theoretical example showcasing historical matrices
#' sizevector <- c(1, 2, 3) # These sizes are not from the original paper
#' stagevector <- c("Sdl", "Veg", "Flo")
#' repvector <- c(0, 0, 1)
#' obsvector <- c(1, 1, 1)
#' matvector <- c(0, 1, 1)
#' immvector <- c(1, 0, 0)
#' propvector <- c(1, 0, 0)
#' indataset <- c(1, 1, 1)
#' binvec <- c(0.5, 0.5, 0.5)
#' 
#' exframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' A1 <- matrix(c(0.10, 0, 0, 0.12, 0, 0, 0.15, 0, 0,
#'   0.15, 0, 0, 0.17, 0, 0, 0.20, 0, 0,
#'   0.20, 0, 0, 0.22, 0, 0, 0.25, 0, 0,
#'   0, 0.20, 0, 0, 0.22, 0, 0, 0.25, 0,
#'   0, 0.25, 0, 0, 0.27, 0, 0, 0.30, 0,
#'   0, 0.30, 0, 0, 0.32, 0, 0, 0.35, 0,
#'   0, 0, 2.00, 0, 0, 3.00, 0, 0, 4.00,
#'   0, 0, 0.35, 0, 0, 0.37, 0, 0, 0.40,
#'   0, 0, 0.40, 0, 0, 0.42, 0, 0, 0.45), 9, 9, byrow = TRUE)
#' 
#' A2 <- matrix(c(0.10, 0, 0, 0.12, 0, 0, 0.15, 0, 0,
#'   0.15, 0, 0, 0.17, 0, 0, 0.20, 0, 0,
#'   0.20, 0, 0, 0.22, 0, 0, 0.25, 0, 0,
#'   0, 0.20, 0, 0, 0.22, 0, 0, 0.25, 0,
#'   0, 0.25, 0, 0, 0.27, 0, 0, 0.30, 0,
#'   0, 0.30, 0, 0, 0.32, 0, 0, 0.35, 0,
#'   0, 0, 5.00, 0, 0, 6.00, 0, 0, 7.00,
#'   0, 0, 0.35, 0, 0, 0.37, 0, 0, 0.40,
#'   0, 0, 0.40, 0, 0, 0.42, 0, 0, 0.45), 9, 9, byrow = TRUE)
#' 
#' A3 <- matrix(c(0.10, 0, 0, 0.12, 0, 0, 0.15, 0, 0,
#'   0.15, 0, 0, 0.17, 0, 0, 0.20, 0, 0,
#'   0.20, 0, 0, 0.22, 0, 0, 0.25, 0, 0,
#'   0, 0.20, 0, 0, 0.22, 0, 0, 0.25, 0,
#'   0, 0.25, 0, 0, 0.27, 0, 0, 0.30, 0,
#'   0, 0.30, 0, 0, 0.32, 0, 0, 0.35, 0,
#'   0, 0, 8.00, 0, 0, 9.00, 0, 0, 10.00,
#'   0, 0, 0.35, 0, 0, 0.37, 0, 0, 0.40,
#'   0, 0, 0.40, 0, 0, 0.42, 0, 0, 0.45), 9, 9, byrow = TRUE)
#' 
#' B1 <- matrix(c(0.10, 0, 0, 0.12, 0, 0, 0.15, 0, 0,
#'   0.15, 0, 0, 0.17, 0, 0, 0.20, 0, 0,
#'   0.20, 0, 0, 0.22, 0, 0, 0.25, 0, 0,
#'   0, 0.20, 0, 0, 0.22, 0, 0, 0.25, 0,
#'   0, 0.25, 0, 0, 0.27, 0, 0, 0.30, 0,
#'   0, 0.30, 0, 0, 0.32, 0, 0, 0.35, 0,
#'   0, 0, 11.00, 0, 0, 12.00, 0, 0, 13.00,
#'   0, 0, 0.35, 0, 0, 0.37, 0, 0, 0.40,
#'   0, 0, 0.40, 0, 0, 0.42, 0, 0, 0.45), 9, 9, byrow = TRUE)
#' 
#' B2 <- matrix(c(0.10, 0, 0, 0.12, 0, 0, 0.15, 0, 0,
#'   0.15, 0, 0, 0.17, 0, 0, 0.20, 0, 0,
#'   0.20, 0, 0, 0.22, 0, 0, 0.25, 0, 0,
#'   0, 0.20, 0, 0, 0.22, 0, 0, 0.25, 0,
#'   0, 0.25, 0, 0, 0.27, 0, 0, 0.30, 0,
#'   0, 0.30, 0, 0, 0.32, 0, 0, 0.35, 0,
#'   0, 0, 14.00, 0, 0, 15.00, 0, 0, 16.00,
#'   0, 0, 0.35, 0, 0, 0.37, 0, 0, 0.40,
#'   0, 0, 0.40, 0, 0, 0.42, 0, 0, 0.45), 9, 9, byrow = TRUE)
#' 
#' B3 <- matrix(c(0.10, 0, 0, 0.12, 0, 0, 0.15, 0, 0,
#'   0.15, 0, 0, 0.17, 0, 0, 0.20, 0, 0,
#'   0.20, 0, 0, 0.22, 0, 0, 0.25, 0, 0,
#'   0, 0.20, 0, 0, 0.22, 0, 0, 0.25, 0,
#'   0, 0.25, 0, 0, 0.27, 0, 0, 0.30, 0,
#'   0, 0.30, 0, 0, 0.32, 0, 0, 0.35, 0,
#'   0, 0, 17.00, 0, 0, 18.00, 0, 0, 19.00,
#'   0, 0, 0.35, 0, 0, 0.37, 0, 0, 0.40,
#'   0, 0, 0.40, 0, 0, 0.42, 0, 0, 0.45), 9, 9, byrow = TRUE)
#' 
#' histmats <- list(A1, A2, A3, B1, B2, B3)
#' stageframe <- exframe
#' pch_ord <- c("A", "A", "A", "B", "B", "B")
#' yr_ord <- c(1, 2, 3, 1, 2, 3)
#' 
#' hist_trial <- create_lM(histmats, exframe, historical = TRUE,
#'   UFdecomp = TRUE, entrystage = 1, patchorder = pch_ord, yearorder = yr_ord)
#'   
#' @export
create_lM <- function(mats, stageframe = NULL, hstages = NA, agestages = NA,
  historical = FALSE, agebystage = FALSE, UFdecomp = TRUE, entrystage = 1,
  poporder = 1, patchorder = 1, yearorder = NA, matrix_id = NULL, add_FC = TRUE,
  sparse_output = FALSE) {
  
  output <- NULL
  
  if (is(mats, "CompadreDB")) {
    output <- .import_Com(matrix_id = matrix_id, database = mats,
      add_FC = add_FC)
  } else if (is.list(mats)) {
    if (is.matrix(mats[[1]]) | is(mats[[1]], "dgCMatrix")) {
      output <- .import_mats(mats = mats, stageframe = stageframe,
        hstages = hstages, agestages = agestages, historical = historical,
        agebystage = agebystage, UFdecomp = UFdecomp, entrystage = entrystage,
        poporder = poporder, patchorder = patchorder, yearorder = yearorder,
        sparse_output = sparse_output)
    } else if ("metadata" %in% names(mats)) {
      output <- .import_Com(matrix_id = matrix_id, database = mats,
        add_FC = add_FC)
    }
  } else {
    stop("Object mats must be an object of class list, or a database holding matrices.", call. = FALSE)
  }
  
  return(output)
}

#' Add Matrices to lefkoMat Object
#' 
#' Function \code{add_lM()} adds matrices to lefkoMat objects.
#' 
#' @param lM The lefkoMat object to add matrices to.
#' @param Amats Either a single \code{A} matrix, or a list of \code{A} matrices.
#' Not necessary if \code{Umats} and \code{Fmats} are both provided.
#' @param Umats Either a single \code{U} matrix, or a list of \code{U} matrices.
#' Not necessary if \code{Amats} and \code{Fmats} are both provided, or if
#' \code{UFdecomp = TRUE} and \code{entrystage} is provided.
#' @param Fmats Either a single \code{F} matrix, or a list of \code{U} matrices.
#' Not necessary if \code{Amats} and \code{Umats} are both provided, or if
#' \code{UFdecomp = TRUE} and \code{entrystage} is provided.
#' @param UFdecomp A logical value indicating whether U and F matrices should be
#' inferred from A matrices and the given \code{entrystage}. Defaults to TRUE.
#' @param entrystage The stage or stages produced by reproductive individuals.
#' Used to determine which transitions are reproductive for U-F decomposition.
#' Defaults to \code{1}, which corresponds to the first stage in the stageframe.
#' @param pop The population designation for each matrix. If object \code{lM}
#' includes only a single population, then defaults to that designation.
#' Otherwise requires a designation as input.
#' @param patch The patch designation for each matrix. If object \code{lM}
#' includes only a single patch, then defaults to that designation. Otherwise
#' requires a designation as input.
#' @param year The designation for occasion at time \emph{t} corresponding to
#' each matrix. Cannot be left empty.
#' 
#' @return A \code{lefkoMat} object incorporating the new matrices within the
#' object input in \code{lM}. 
#' 
#' @section Notes:
#' This function will not allow matrices of different dimension from those input
#' in object \code{lM} to be added to that object.
#' 
#' Two of \code{Amats}, \code{Umats}, and \code{Fmats} must be provided for this
#' function to proceed. Also, if \code{Amats}, \code{Umats}, and \code{Fmats}
#' are all provided, then this function will default to replacing \code{Amats}
#' with the sum of the respective \code{Umats} and \code{Fmats}.
#' 
#' @seealso \code{\link{create_lM}()}
#' @seealso \code{\link{delete_lM}()}
#' @seealso \code{\link{subset_lM}()}
#' 
#' @examples
#' # These matrices are of 9 populations of the plant species Anthyllis
#' # vulneraria, and were originally published in Davison et al. (2010) Journal
#' # of Ecology 98:255-267 (doi: 10.1111/j.1365-2745.2009.01611.x).
#' 
#' sizevector <- c(1, 1, 2, 3) # These sizes are not from the original paper
#' stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
#' repvector <- c(0, 0, 1, 1)
#' obsvector <- c(1, 1, 1, 1)
#' matvector <- c(0, 1, 1, 1)
#' immvector <- c(1, 0, 0, 0)
#' propvector <- c(0, 0, 0, 0)
#' indataset <- c(1, 1, 1, 1)
#' binvec <- c(0.5, 0.5, 0.5, 0.5)
#' 
#' anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' # POPN C 2003-2004
#' XC3 <- matrix(c(0, 0, 1.74, 1.74,
#' 0.208333333, 0, 0, 0.057142857,
#' 0.041666667, 0.076923077, 0, 0,
#' 0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
#' 
#' # 2004-2005
#' XC4 <- matrix(c(0, 0, 0.3, 0.6,
#' 0.32183908, 0.142857143, 0, 0,
#' 0.16091954, 0.285714286, 0, 0,
#' 0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
#' 
#' # 2005-2006
#' XC5 <- matrix(c(0, 0, 0.50625, 0.675,
#' 0, 0, 0, 0.035714286,
#' 0.1, 0.068965517, 0.0625, 0.107142857,
#' 0.3, 0.137931034, 0, 0.071428571), 4, 4, byrow = TRUE)
#' 
#' # POPN E 2003-2004
#' XE3 <- matrix(c(0, 0, 2.44, 6.569230769,
#' 0.196428571, 0, 0, 0,
#' 0.125, 0.5, 0, 0,
#' 0.160714286, 0.5, 0.133333333, 0.076923077), 4, 4, byrow = TRUE)
#' 
#' XE4 <- matrix(c(0, 0, 0.45, 0.646153846,
#' 0.06557377, 0.090909091, 0.125, 0,
#' 0.032786885, 0, 0.125, 0.076923077,
#' 0.049180328, 0, 0.125, 0.230769231), 4, 4, byrow = TRUE)
#' 
#' XE5 <- matrix(c(0, 0, 2.85, 3.99,
#' 0.083333333, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.416666667, 0.1, 0, 0.1), 4, 4, byrow = TRUE)
#' 
#' mats_list <- list(XC3, XC4, XC5, XE3, XE4, XE5)
#' yr_ord <- c(1, 2, 3, 1, 2, 3)
#' pch_ord <- c(1, 1, 1, 2, 2, 2)
#' 
#' anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA,
#'   historical = FALSE, poporder = 1, patchorder = pch_ord, yearorder = yr_ord)
#'   
#' XH3 <- matrix(c(0, 0, 0.1125, 1.05,
#' 0.2, 0, 0, 0,
#' 0, 0.5, 0, 0,
#' 0.2, 0.5, 0, 0), 4, 4, byrow = TRUE)
#' 
#' XH3u <- matrix(c(0, 0, 0, 0,
#' 0.2, 0, 0, 0,
#' 0, 0.5, 0, 0,
#' 0.2, 0.5, 0, 0), 4, 4, byrow = TRUE)
#' 
#' XH4 <- matrix(c(0, 0, 0, 0,
#' 0, 0, 0.5, 0,
#' 0.8, 0.5, 0.25, 0.25,
#' 0.2, 0, 0, 0.75), 4, 4, byrow = TRUE)
#' 
#' XH4u <- matrix(c(0, 0, 0, 0,
#' 0, 0, 0.5, 0,
#' 0.8, 0.5, 0.25, 0.25,
#' 0.2, 0, 0, 0.75), 4, 4, byrow = TRUE)
#' 
#' XH5 <- matrix(c(0, 0, 0.2, 1.05,
#' 0, 0, 0, 0,
#' 0.001, 0.001, 0.333333333, 0,
#' 0.001, 0, 0, 0), 4, 4, byrow = TRUE)
#' 
#' XH5u <- matrix(c(0, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.001, 0.001, 0.333333333, 0,
#' 0.001, 0, 0, 0), 4, 4, byrow = TRUE)
#' 
#' anth_lefkoMat <- add_lM(anth_lefkoMat, Amats = list(XH3, XH4, XH5),
#'   Umats = list(XH3u, XH4u, XH5u), patch = c(3, 3, 3), year = c(1, 2, 3))
#'   
#' @export
add_lM <- function(lM, Amats = NA, Umats = NA, Fmats = NA, UFdecomp = FALSE,
  entrystage = 1, pop = NA, patch = NA, year = NA) {
  
  F_indices <- numstages <- NULL
  
  if (!is(lM, "lefkoMat")) {
    stop("This function requires a lefkoMat object as input.", call. = FALSE)
  }
  
  if (all(is.na(Amats)) & all(is.na(Umats))) {
    stop("Please add either Amats or Umats.", call. = FALSE)
  }
  if (all(is.na(Amats)) & all(is.na(Fmats))) {
    stop("Please add either Amats or Fmats.", call. = FALSE)
  }
  if (all(is.na(Umats)) & all(is.na(Fmats)) & !UFdecomp) {
    stop("Please add either Umats or Fmats, or add U-F decomposition information.",
      call. = FALSE)
  }
  
  mat_dim <- dim(lM$A[[1]])[1]
  numstages <- length(lM$ahstages$stage_id)
  
  if (!all(is.na(lM$hstages))) {
    historical <- TRUE
  } else {
    historical <- FALSE
  }
  
  if (!all(is.na(lM$agestages))) {
    agebystage <- TRUE
  } else {
    agebystage <- FALSE
  }
  
  if (historical & agebystage) {
    stop("This function cannot proceed with historical age-by-stage MPMs.",
      call. = FALSE)
  }
  
  if (!all(is.na(Amats))) {
    if (is.matrix(Amats)) {
      Amats <- list(Amats)
    }
    
    if (!all(unique(unlist(lapply(Amats, dim))) == mat_dim)) {
      stop("Input A matrices must be of the same dimensions as the matrices in object lM", call. = FALSE)
    }
  }
  
  if (UFdecomp & all(is.na(Umats))) {
    if (all(is.na(entrystage))) {
      warning("No entry stage provided, so assuming that first stage is entry stage.", call. = FALSE)
      entrystage <- 1
    } else if (all(is.element(entrystage, lM$ahstages$stage)) | 
      all(is.element(entrystage, lM$ahstages$stage_id))) {
      
      total_entries <- length(entrystage)
      
      if (!all(is.element(entrystage, lM$ahstages$stage_id))) {
        entrystage <- apply(as.matrix(entrystage), 1, function(X) {
          return(which(lM$ahstages$stage == X))
        })
      }
    } else {
      stop("Unable to interpret entry stage designations.", call. = FALSE)
    }
    
    rep_from <- which(lM$ahstages$repstatus == 1)
    
    # Indexing
    if (!historical & !agebystage) { #Ahistorical
      for (i in c(1:length(rep_from))) {
        for (j in c(1:length(entrystage))) {
          if (i == 1 & j == 1) {
            F_indices <- ((numstages * (rep_from[i] - 1)) + entrystage[j])
          } else {
            F_indices <- append(F_indices, ((numstages * (rep_from[i] - 1)) + entrystage[j]))
          }
        }
      }
    } else if (historical & !agebystage) { #Historical
      h_numstages <- length(lM$hstages$stage_2)
      
      for (time2o1 in c(1:h_numstages)) {
        for (time32n in c(1:h_numstages)) {
          if (lM$hstages$stage_id_2[time2o1] == lM$hstages$stage_id_1[time32n]) {
            if (is.element(lM$hstages$stage_id_2[time2o1], rep_from) & 
              is.element(lM$hstages$stage_id_2[time32n], entrystage)) {
              
              if (length(F_indices) == 0) {
                F_indices <- ((time2o1 - 1) * h_numstages) + time32n
              } else {
                F_indices <- append(F_indices, (((time2o1 - 1) * h_numstages) + time32n))
              }
            }
          }
        }
      }
    } else if (agebystage & !historical) { #Age-by-stage ahistorical
      age_numstages <- length(lM$agestages$stage_id)
      
      for (fromas in c(1:age_numstages)) {
        for (toas in c(1:age_numstages)) {
          if (is.element(lM$agestages$stage_id[fromas], rep_from) &
            is.element(lM$agestages$stage_id[toas], entrystage)) {
            
            if (length(F_indices) == 0) {
              F_indices <- ((fromas - 1) * age_numstages) + toas
            } else {
              F_indices <- append(F_indices, (((fromas - 1) * age_numstages) + toas))
            }
          }
        }
      }
    } else {
      stop("Requested format unrecognized. Cannot proceed.", call. = FALSE)
    }
        
    Umats <- lapply(Amats, function(X) {
      newmat <- X
      newmat[F_indices] <- 0
      
      return(newmat)
    })
    
    nonFindices <- c(1:length(Amats[[1]]))[!is.element(c(1:length(Amats[[1]])), F_indices)]
    Fmats <- lapply(Amats, function(X) {
      newmat <- X
      newmat[nonFindices] <- 0
      
      return(newmat)
    })
  }
  
  if (!all(is.na(Umats))) {
    if (is.matrix(Umats)) {
      Umats <- list(Umats)
    }
    
    if (!all(unique(unlist(lapply(Umats, dim))) == mat_dim)) {
      stop("Input U matrices must be of the same dimensions as the matrices in object lM", call. = FALSE)
    }
  }
  if (!all(is.na(Fmats))) {
    if (is.matrix(Fmats)) {
      Fmats <- list(Fmats)
    }
    
    if (!all(unique(unlist(lapply(Fmats, dim))) == mat_dim)) {
      stop("Input F matrices must be of the same dimensions as the matrices in object lM", call. = FALSE)
    }
  }
  
  if (!UFdecomp) {
    if (!all(is.na(Amats)) & !all(is.na(Umats))) {
      list_lengthA <- length(Amats)
      list_lengthU <- length(Umats)
      
      if (list_lengthA != list_lengthU) {
        stop("Input Amats and Umats objects must include the same number of matrices.", call. = FALSE)
      }
      
      Fmats <- lapply(as.list(1:list_lengthA), function(X) {
        return(Amats[[X]] - Umats[[X]])
      })
      
    } else if (!all(is.na(Amats)) & !all(is.na(Fmats))) {
      list_lengthA <- length(Amats)
      list_lengthF <- length(Fmats)
      
      if (list_lengthA != list_lengthF) {
        stop("Input Amats and Fmats objects must include the same number of matrices.", call. = FALSE)
      }
      
      Umats <- lapply(as.list(1:list_lengthA), function(X) {
        return(Amats[[X]] - Fmats[[X]])
      })
      
    } else if (!all(is.na(Umats)) & !all(is.na(Fmats))) {
      list_lengthU <- length(Umats)
      list_lengthF <- length(Fmats)
      
      if (list_lengthU != list_lengthF) {
        stop("Input Umats and Fmats objects must include the same number of matrices.", call. = FALSE)
      }
      
      Amats <- lapply(as.list(1:list_lengthU), function(X) {
        return(Umats[[X]] + Fmats[[X]])
      })
      
      list_lengthA <- list_lengthU
    }
  } else {
    list_lengthA <- length(Amats)
  }
  
  if (all(is.na(year)) | length(year) != list_lengthA) {
    stop("This function requires occasion in time t for each matrix to be added.", call. = FALSE)
  }
  if(all(is.na(pop))) {
    popsinlM <- unique(lM$labels$pop)
    
    if (length(popsinlM) == 1) {
      pop <- rep(popsinlM, list_lengthA)
    } else {
      stop("Please input the population designation for each matrix to be added.", call. = FALSE)
    }
  }
  if(all(is.na(patch))) {
    patchesinlM <- unique(lM$labels$patch)
    
    if (length(patchesinlM) == 1) {
      patch <- rep(patchesinlM, list_lengthA)
    } else {
      stop("Please input the patch designation for each matrix to be added.", call. = FALSE)
    }
  }
  
  #Now the appending
  lM$A <- append(lM$A, Amats)
  lM$U <- append(lM$U, Umats)
  lM$F <- append(lM$F, Fmats)
  
  surv_additions <- sum(unlist(lapply(Umats, function(X) {
    length(which(as.matrix(X) > 0))
  })))
  fec_additions <- sum(unlist(lapply(Fmats, function(X) {
    length(which(as.matrix(X) > 0))
  })))
  
  lM$matrixqc[1] <- lM$matrixqc[1] + surv_additions
  lM$matrixqc[2] <- lM$matrixqc[2] + fec_additions
  lM$matrixqc[3] <- lM$matrixqc[3] + list_lengthA
  
  newlabels <- cbind.data.frame(pop = pop, patch = patch, year2 = year)
  lM$labels <- rbind.data.frame(lM$labels, newlabels)
  
  return(lM)
}

#' Delete Matrices from lefkoMat Object
#' 
#' Function \code{delete_lM()} deletes matrices from \code{lefkoMat} objects.
#' 
#' @param lM The \code{lefkoMat} object to delete matrices from.
#' @param mat_num Either a single integer corresponding to the matrix to remove
#' within the \code{labels} element of \code{lM}, or a vector of such integers.
#' @param pop The population designation for matrices to remove. Only used if
#' \code{mat_num} is not given.
#' @param patch The patch designation for matrices to remove. Only used if
#' \code{mat_num} is not given.
#' @param year The time \emph{t} designation for matrices to remove. Only used
#' if \code{mat_num} is not given.
#' 
#' @return A \code{lefkoMat} object in which the matrices specified in \code{lM}
#' have been removed. 
#' 
#' @section Notes:
#' If \code{mat_num} is not provided, then at least one of \code{pop},
#' \code{patch}, or \code{year} must be provided. If at least two of \code{pop},
#' \code{patch}, and \code{year} are provided, then function \code{detele_lM()}
#' will identify matrices to remove as the intersection of provided inputs.
#' 
#' @seealso \code{\link{create_lM}()}
#' @seealso \code{\link{add_lM}()}
#' @seealso \code{\link{subset_lM}()}
#' 
#' @examples
#' # These matrices are of 9 populations of the plant species Anthyllis
#' # vulneraria, and were originally published in Davison et al. (2010) Journal
#' # of Ecology 98:255-267 (doi: 10.1111/j.1365-2745.2009.01611.x).
#' 
#' sizevector <- c(1, 1, 2, 3) # These sizes are not from the original paper
#' stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
#' repvector <- c(0, 0, 1, 1)
#' obsvector <- c(1, 1, 1, 1)
#' matvector <- c(0, 1, 1, 1)
#' immvector <- c(1, 0, 0, 0)
#' propvector <- c(0, 0, 0, 0)
#' indataset <- c(1, 1, 1, 1)
#' binvec <- c(0.5, 0.5, 0.5, 0.5)
#' 
#' anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' # POPN C 2003-2004
#' XC3 <- matrix(c(0, 0, 1.74, 1.74,
#' 0.208333333, 0, 0, 0.057142857,
#' 0.041666667, 0.076923077, 0, 0,
#' 0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
#' 
#' # 2004-2005
#' XC4 <- matrix(c(0, 0, 0.3, 0.6,
#' 0.32183908, 0.142857143, 0, 0,
#' 0.16091954, 0.285714286, 0, 0,
#' 0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
#' 
#' # 2005-2006
#' XC5 <- matrix(c(0, 0, 0.50625, 0.675,
#' 0, 0, 0, 0.035714286,
#' 0.1, 0.068965517, 0.0625, 0.107142857,
#' 0.3, 0.137931034, 0, 0.071428571), 4, 4, byrow = TRUE)
#' 
#' # POPN E 2003-2004
#' XE3 <- matrix(c(0, 0, 2.44, 6.569230769,
#' 0.196428571, 0, 0, 0,
#' 0.125, 0.5, 0, 0,
#' 0.160714286, 0.5, 0.133333333, 0.076923077), 4, 4, byrow = TRUE)
#' 
#' XE4 <- matrix(c(0, 0, 0.45, 0.646153846,
#' 0.06557377, 0.090909091, 0.125, 0,
#' 0.032786885, 0, 0.125, 0.076923077,
#' 0.049180328, 0, 0.125, 0.230769231), 4, 4, byrow = TRUE)
#' 
#' XE5 <- matrix(c(0, 0, 2.85, 3.99,
#' 0.083333333, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.416666667, 0.1, 0, 0.1), 4, 4, byrow = TRUE)
#' 
#' mats_list <- list(XC3, XC4, XC5, XE3, XE4, XE5)
#' yr_ord <- c(1, 2, 3, 1, 2, 3)
#' pch_ord <- c(1, 1, 1, 2, 2, 2)
#' 
#' anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA,
#'   historical = FALSE, poporder = 1, patchorder = pch_ord, yearorder = yr_ord)
#'   
#' smaller_anth_lM <- delete_lM(anth_lefkoMat, patch = 2)
#' 
#' @export
delete_lM <- function(lM, mat_num = NA, pop = NA, patch = NA, year = NA) {
  
  if (!is(lM, "lefkoMat")) {
    stop("This function requires a lefkoMat object as input.", call. = FALSE)
  }
  if (all(is.na(mat_num)) & all(is.na(pop)) & all(is.na(patch)) & all(is.na(year))) {
    stop("Please add some input to determine which matrices to remove.", call. = FALSE)
  }
  
  mat_length <- length(lM$A)
  mat_dim <- dim(lM$A[[1]])[1]
  pop_set <- unique(lM$labels$pop)
  patch_set <- unique(lM$labels$patch)
  year_set <- unique(lM$labels$year2)
  
  if (!all(is.na(mat_num))) {
    if (max(mat_num) > mat_length | min(mat_num) < 1) {
      stop("Matrices chosen for deletion are outside the range of matrices input in object lM.",
        call. = FALSE)
    }
  } else {
    
    pop_subset <- patch_subset <- year_subset <- seq(from = 1, to = mat_length)
    
    if (!all(is.na(pop))) {
      if (!all(is.element(pop, pop_set))) {
        stop("Option pop includes terms not representing known populations in object lM.",
          call. = FALSE)
      }
      
      pop_subset <- which(is.element(lM$labels$pop, pop))
    }
    if (!all(is.na(patch))) {
      if (!all(is.element(patch, patch_set))) {
        stop("Option patch includes terms not representing known patches in object lM.",
          call. = FALSE)
      }
      
      patch_subset <- which(is.element(lM$labels$patch, patch))
    }
    if (!all(is.na(year))) {
      if (!all(is.element(year, year_set))) {
        stop("Option year includes terms not representing known times in object lM.",
          call. = FALSE)
      }
      
      year_subset <- which(is.element(lM$labels$year2, year))
    }
    
    mat_num <- intersect(intersect(pop_subset, patch_subset), year_subset)
  }
  
  if (length(unique(mat_num)) >= mat_length) {
    stop("Options suggest that all matrices are to be deleted. Cannot proceed.",
      call. = FALSE)
  }
  
  #Now the deletion
  lM$A <- lM$A[-mat_num]
  lM$U <- lM$U[-mat_num]
  lM$F <- lM$F[-mat_num]
  
  surv_portions <- sum(unlist(lapply(lM$U, function(X) {
    length(which(as.matrix(X) > 0))
  })))
  fec_portions <- sum(unlist(lapply(lM$F, function(X) {
    length(which(as.matrix(X) > 0))
  })))
  
  lM$matrixqc[1] <- surv_portions
  lM$matrixqc[2] <- fec_portions
  lM$matrixqc[3] <- lM$matrixqc[3] - length(mat_num)
  
  lM$labels <- lM$labels[-mat_num,]
  rownames(lM$labels) <- seq(from = 1, to = length(lM$A))
  
  return(lM)
}

#' Create New lefkoMat Object as Subset of Another lefkoMat Object
#' 
#' Function \code{subset_lM()} creates a new \code{lefkoMat} object from a
#' subset of matrices in another \code{lefkoMat} object.
#' 
#' @param lM The \code{lefkoMat} object to select matrices from.
#' @param mat_num Either a single integer corresponding to the matrix to select
#' within the \code{labels} element of \code{lM}, or a vector of such integers.
#' @param pop The population designation for matrices to select. Only used if
#' \code{mat_num} is not given.
#' @param patch The patch designation for matrices to select. Only used if
#' \code{mat_num} is not given.
#' @param year The time \emph{t} designation for matrices to select. Only used
#' if \code{mat_num} is not given.
#' 
#' @return A \code{lefkoMat} object composed of the matrices specified in the
#' options. 
#' 
#' @section Notes:
#' If \code{mat_num} is not provided, then at least one of \code{pop},
#' \code{patch}, or \code{year} must be provided. If at least two of \code{pop},
#' \code{patch}, and \code{year} are provided, then function \code{subset_lM()}
#' will identify matrices as the intersection of provided inputs.
#' 
#' @seealso \code{\link{create_lM}()}
#' @seealso \code{\link{add_lM}()}
#' @seealso \code{\link{delete_lM}()}
#' 
#' @examples
#' # These matrices are of 9 populations of the plant species Anthyllis
#' # vulneraria, and were originally published in Davison et al. (2010) Journal
#' # of Ecology 98:255-267 (doi: 10.1111/j.1365-2745.2009.01611.x).
#' 
#' sizevector <- c(1, 1, 2, 3) # These sizes are not from the original paper
#' stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
#' repvector <- c(0, 0, 1, 1)
#' obsvector <- c(1, 1, 1, 1)
#' matvector <- c(0, 1, 1, 1)
#' immvector <- c(1, 0, 0, 0)
#' propvector <- c(0, 0, 0, 0)
#' indataset <- c(1, 1, 1, 1)
#' binvec <- c(0.5, 0.5, 0.5, 0.5)
#' 
#' anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' # POPN C 2003-2004
#' XC3 <- matrix(c(0, 0, 1.74, 1.74,
#' 0.208333333, 0, 0, 0.057142857,
#' 0.041666667, 0.076923077, 0, 0,
#' 0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
#' 
#' # 2004-2005
#' XC4 <- matrix(c(0, 0, 0.3, 0.6,
#' 0.32183908, 0.142857143, 0, 0,
#' 0.16091954, 0.285714286, 0, 0,
#' 0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
#' 
#' # 2005-2006
#' XC5 <- matrix(c(0, 0, 0.50625, 0.675,
#' 0, 0, 0, 0.035714286,
#' 0.1, 0.068965517, 0.0625, 0.107142857,
#' 0.3, 0.137931034, 0, 0.071428571), 4, 4, byrow = TRUE)
#' 
#' # POPN E 2003-2004
#' XE3 <- matrix(c(0, 0, 2.44, 6.569230769,
#' 0.196428571, 0, 0, 0,
#' 0.125, 0.5, 0, 0,
#' 0.160714286, 0.5, 0.133333333, 0.076923077), 4, 4, byrow = TRUE)
#' 
#' XE4 <- matrix(c(0, 0, 0.45, 0.646153846,
#' 0.06557377, 0.090909091, 0.125, 0,
#' 0.032786885, 0, 0.125, 0.076923077,
#' 0.049180328, 0, 0.125, 0.230769231), 4, 4, byrow = TRUE)
#' 
#' XE5 <- matrix(c(0, 0, 2.85, 3.99,
#' 0.083333333, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.416666667, 0.1, 0, 0.1), 4, 4, byrow = TRUE)
#' 
#' # POPN F 2003-2004
#' XF3 <- matrix(c(0, 0, 1.815, 7.058333333,
#' 0.075949367, 0, 0.05, 0.083333333,
#' 0.139240506, 0, 0, 0.25,
#' 0.075949367, 0, 0, 0.083333333), 4, 4, byrow = TRUE)
#' 
#' XF4 <- matrix(c(0, 0, 1.233333333, 7.4,
#' 0.223880597, 0, 0.111111111, 0.142857143,
#' 0.134328358, 0.272727273, 0.166666667, 0.142857143,
#' 0.119402985, 0.363636364, 0.055555556, 0.142857143), 4, 4, byrow = TRUE)
#' 
#' XF5 <- matrix(c(0, 0, 1.06, 3.372727273,
#' 0.073170732, 0.025, 0.033333333, 0,
#' 0.036585366, 0.15, 0.1, 0.136363636,
#' 0.06097561, 0.225, 0.166666667, 0.272727273), 4, 4, byrow = TRUE)
#' 
#' # POPN G 2003-2004
#' XG3 <- matrix(c(0, 0, 0.245454545, 2.1,
#' 0, 0, 0.045454545, 0,
#' 0.125, 0, 0.090909091, 0,
#' 0.125, 0, 0.090909091, 0.333333333), 4, 4, byrow = TRUE)
#' 
#' XG4 <- matrix(c(0, 0, 1.1, 1.54,
#' 0.111111111, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.111111111, 0, 0, 0), 4, 4, byrow = TRUE)
#' 
#' XG5 <- matrix(c(0, 0, 0, 1.5,
#' 0, 0, 0, 0,
#' 0.090909091, 0, 0, 0,
#' 0.545454545, 0.5, 0, 0.5), 4, 4, byrow = TRUE)
#' 
#' # POPN L 2003-2004
#' XL3 <- matrix(c(0, 0, 1.785365854, 1.856521739,
#' 0.128571429, 0, 0, 0.010869565,
#' 0.028571429, 0, 0, 0,
#' 0.014285714, 0, 0, 0.02173913), 4, 4, byrow = TRUE)
#' 
#' XL4 <- matrix(c(0, 0, 14.25, 16.625,
#' 0.131443299, 0.057142857, 0, 0.25,
#' 0.144329897, 0, 0, 0,
#' 0.092783505, 0.2, 0, 0.25), 4, 4, byrow = TRUE)
#' 
#' XL5 <- matrix(c(0, 0, 0.594642857, 1.765909091,
#' 0, 0, 0.017857143, 0,
#' 0.021052632, 0.018518519, 0.035714286, 0.045454545,
#' 0.021052632, 0.018518519, 0.035714286, 0.068181818), 4, 4, byrow = TRUE)
#' 
#' # POPN O 2003-2004
#' XO3 <- matrix(c(0, 0, 11.5, 2.775862069,
#' 0.6, 0.285714286, 0.333333333, 0.24137931,
#' 0.04, 0.142857143, 0, 0,
#' 0.16, 0.285714286, 0, 0.172413793), 4, 4, byrow = TRUE)
#' 
#' XO4 <- matrix(c(0, 0, 3.78, 1.225,
#' 0.28358209, 0.171052632, 0, 0.166666667,
#' 0.084577114, 0.026315789, 0, 0.055555556,
#' 0.139303483, 0.447368421, 0, 0.305555556), 4, 4, byrow = TRUE)
#' 
#' XO5 <- matrix(c(0, 0, 1.542857143, 1.035616438,
#' 0.126984127, 0.105263158, 0.047619048, 0.054794521,
#' 0.095238095, 0.157894737, 0.19047619, 0.082191781,
#' 0.111111111, 0.223684211, 0, 0.356164384), 4, 4, byrow = TRUE)
#' 
#' # POPN Q 2003-2004
#' XQ3 <- matrix(c(0, 0, 0.15, 0.175,
#' 0, 0, 0, 0,
#' 0, 0, 0, 0,
#' 1, 0, 0, 0), 4, 4, byrow = TRUE)
#' 
#' XQ4 <- matrix(c(0, 0, 0, 0.25,
#' 0, 0, 0, 0,
#' 0, 0, 0, 0,
#' 1, 0.666666667, 0, 1), 4, 4, byrow = TRUE)
#' 
#' XQ5 <- matrix(c(0, 0, 0, 1.428571429,
#' 0, 0, 0, 0.142857143,
#' 0.25, 0, 0, 0,
#' 0.25, 0, 0, 0.571428571), 4, 4, byrow = TRUE)
#' 
#' # POPN R 2003-2004
#' XR3 <- matrix(c(0, 0, 0.7, 0.6125,
#' 0.25, 0, 0, 0.125,
#' 0, 0, 0, 0,
#' 0.25, 0.166666667, 0, 0.25), 4, 4, byrow = TRUE)
#' 
#' XR4 <- matrix(c(0, 0, 0, 0.6,
#' 0.285714286, 0, 0, 0,
#' 0.285714286, 0.333333333, 0, 0,
#' 0.285714286, 0.333333333, 0, 1), 4, 4, byrow = TRUE)
#' 
#' XR5 <- matrix(c(0, 0, 0.7, 0.6125,
#' 0, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0.333333333, 0, 0.333333333, 0.625), 4, 4, byrow = TRUE)
#' 
#' # POPN S 2003-2004
#' XS3 <- matrix(c(0, 0, 2.1, 0.816666667,
#' 0.166666667, 0, 0, 0,
#' 0, 0, 0, 0,
#' 0, 0, 0, 0.166666667), 4, 4, byrow = TRUE)
#' 
#' XS4 <- matrix(c(0, 0, 0, 7,
#' 0.333333333, 0.5, 0, 0,
#' 0, 0, 0, 0,
#' 0.333333333, 0, 0, 1), 4, 4, byrow = TRUE)
#' 
#' XS5 <- matrix(c(0, 0, 0, 1.4,
#' 0, 0, 0, 0,
#' 0, 0, 0, 0.2,
#' 0.111111111, 0.75, 0, 0.2), 4, 4, byrow = TRUE)
#' 
#' mats_list <- list(XC3, XC4, XC5, XE3, XE4, XE5, XF3, XF4, XF5, XG3, XG4, XG5,
#'   XL3, XL4, XL5, XO3, XO4, XO5, XQ3, XQ4, XQ5, XR3, XR4, XR5, XS3, XS4, XS5)
#' 
#' yr_ord <- c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1,
#'   2, 3, 1, 2, 3)
#' 
#' pch_ord <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7,
#'   8, 8, 8, 9, 9, 9)
#' 
#' anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA, historical = FALSE,
#'   poporder = 1, patchorder = pch_ord, yearorder = yr_ord)
#'   
#' smaller_anth_lM <- subset_lM(anth_lefkoMat, patch = c(1, 2, 3), 
#'   year = c(1, 2))
#' 
#' @export
subset_lM <- function(lM, mat_num = NA, pop = NA, patch = NA, year = NA) {
  
  if (!is(lM, "lefkoMat")) {
    stop("This function requires a lefkoMat object as input.", call. = FALSE)
  }
  if (all(is.na(mat_num)) & all(is.na(pop)) & all(is.na(patch)) & all(is.na(year))) {
    stop("Please add some input to determine which matrices to select.", call. = FALSE)
  }
  
  mat_length <- length(lM$A)
  mat_dim <- dim(lM$A[[1]])[1]
  pop_set <- unique(lM$labels$pop)
  patch_set <- unique(lM$labels$patch)
  year_set <- unique(lM$labels$year2)
  
  if (!all(is.na(mat_num))) {
    if (max(mat_num) > mat_length | min(mat_num) < 1) {
      stop("Matrices chosen for selection are outside the range of matrices input in object lM.",
        call. = FALSE)
    }
  } else {
    
    pop_subset <- patch_subset <- year_subset <- seq(from = 1, to = mat_length)
    
    if (!all(is.na(pop))) {
      if (!all(is.element(pop, pop_set))) {
        stop("Option pop includes terms not representing known populations in object lM.",
          call. = FALSE)
      }
      
      pop_subset <- which(is.element(lM$labels$pop, pop))
    }
    if (!all(is.na(patch))) {
      if (!all(is.element(patch, patch_set))) {
        stop("Option patch includes terms not representing known patches in object lM.",
          call. = FALSE)
      }
      
      patch_subset <- which(is.element(lM$labels$patch, patch))
    }
    if (!all(is.na(year))) {
      if (!all(is.element(year, year_set))) {
        stop("Option year includes terms not representing known times in object lM.",
          call. = FALSE)
      }
      
      year_subset <- which(is.element(lM$labels$year2, year))
    }
    
    mat_num <- intersect(intersect(pop_subset, patch_subset), year_subset)
  }
  
  if (length(unique(mat_num)) >= mat_length) {
    stop("Options suggest that all matrices are to be selected. Cannot proceed.",
      call. = FALSE)
  }
  
  #Now the selection
  lM$A <- lM$A[mat_num]
  lM$U <- lM$U[mat_num]
  lM$F <- lM$F[mat_num]
  
  surv_portions <- sum(unlist(lapply(lM$U, function(X) {
    length(which(as.matrix(X) > 0))
  })))
  fec_portions <- sum(unlist(lapply(lM$F, function(X) {
    length(which(as.matrix(X) > 0))
  })))
  
  lM$matrixqc[1] <- surv_portions
  lM$matrixqc[2] <- fec_portions
  lM$matrixqc[3] <- length(mat_num)
    
  lM$labels <- lM$labels[mat_num,]
  rownames(lM$labels) <- seq(from = 1, to = length(lM$A))
  
  return(lM)
}

