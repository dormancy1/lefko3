#' Create Overwrite Table for MPM Development
#' 
#' \code{overwrite()} returns a data frame describing which particular
#' transitions within an ahistorical or historical projection matrix to
#' overwrite with either given rates and probabilities, or other estimated
#' transitions. This function is now deprecated in favor of function
#' \code{\link{supplemental}()}.
#' 
#' @name overwrite
#' 
#' @param stage3 The name of the stage in occasion \emph{t}+1 in the transition
#' to be replaced. Abbreviations for groups of stages are also allowed
#' (see Notes).
#' @param stage2 The name of the stage in occasion \emph{t} in the transition to
#' be replaced. Abbreviations for groups of stages are also allowed (see Notes).
#' @param stage1 The name of the stage in occasion \emph{t}-1 in the transition
#' to be replaced. Only needed if a historical matrix is to be produced.
#' Abbreviations for groups of stages are also allowed (see Notes).
#' @param eststage3 The name of the stage to replace \code{stage3}. Only needed
#' if a transition will be replaced by another estimated transition.
#' @param eststage2 The name of the stage to replace \code{stage2}. Only needed
#' if a transition will be replaced by another estimated transition.
#' @param eststage1 The name of the stage to replace \code{stage1}. Only needed
#' if a transition will be replaced by another estimated transition, and the
#' matrix to be estimated is historical.
#' @param givenrate A fixed rate or probability to replace for the transition
#' described by \code{stage3}, \code{stage2}, and \code{stage1}.
#' @param type A vector denoting the kind of transition between occasions
#' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
#' \code{S}, or \code{s} for the replacement of a survival transition; or
#' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
#' transition. If empty or not provided, then defaults to \code{1} for survival
#' transition.
#' @param type_t12 An optional vector denoting the kind of transition between
#' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
#' deVries format is desired. This should be entered as \code{1}, \code{S}, or
#' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
#' fecundity transitions. Defaults to \code{1} for survival transition, with
#' impacts only on the construction of deVries-format hMPMs.
#'
#' @return A data frame that puts the above vectors together and can be used as
#' input in \code{\link{flefko3}()}, \code{\link{flefko2}()},
#' \code{\link{rlefko3}()},\code{\link{rlefko2}()}, and
#' \code{\link{aflefko2}()}.
#' 
#' Variables in this data frame include the following:
#' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
#' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
#' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced.}
#' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace the
#' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{givenrate}{A constant to be used as the value of the transition.}
#' \item{convtype}{Designates whether the transition from occasion \emph{t} to
#' occasion \emph{t}+1 is a survival-transition probability (1) or a fecundity
#' rate (2).}
#' \item{convtype_t12}{Designates whether the transition from occasion
#' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), a
#' fecundity rate (2).}
#' 
#' @section Notes:
#' This function is deprecated. Please use \code{\link{supplemental}()}.
#' 
#' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
#' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
#' are to be used, \code{nrep} if all mature but non-reproductive stages are to
#' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
#' immature stages are to be used, \code{prop} if all propagule stages are to be
#' used, \code{npr} if all non-propagule stages are to be used, and leave empty
#' or use \code{all} if all stages in stageframe are to be used.
#' 
#' @examples
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm"),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm"),
#'   givenrate = c(0.1, 0.2, 0.2, 0.2, 0.25, NA, NA, NA),
#'   type = c("S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypover2r
#' 
#' cypover3r <- overwrite(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL", 
#'     "D", "XSm", "Sm", "D", "XSm", "Sm"),
#'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL",
#'     "SL", "SL", "SL"),
#'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
#'     "SL", "SL", "SL"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm",
#'     "Sm"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm",
#'     "XSm", "XSm"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm",
#'     "XSm", "XSm"),
#'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA),
#'   type = c("S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypover3r
#' 
#' @export
overwrite <- function(stage3, stage2, stage1 = NA, eststage3 = NA, 
  eststage2 = NA, eststage1 = NA, givenrate = NA, type = NA, type_t12 = NA) {
  
  message("This function is deprecated. Please use supplemental() instead")
  
  if (length(stage3) != length(stage2)) {
    stop("All transitions to overwrite require information at least for stage2
      and stage3. These inputs must also be of equal length.", call. = FALSE)
  }
  
  fulllength <- max(length(stage3), length(stage2), length(stage1), length(eststage3), 
    length(eststage2), length(eststage1), length(givenrate), length(type))
  
  if (length(stage3) != fulllength) {
    stop("Please provide all input vectors in the same order.", call. = FALSE)
  }
  
  if (length(stage1) < fulllength) {
    missinglength <- fulllength - length(stage1)
    stage1 <- as.character(append(stage1, rep(NA, missinglength)))
  }
  if (length(eststage3) < fulllength) {
    missinglength <- fulllength - length(eststage3)
    eststage3 <- as.character(append(eststage3, rep(NA, missinglength)))
  }
  if (length(eststage2) < fulllength) {    
    missinglength <- fulllength - length(eststage2)
    eststage2 <- as.character(append(eststage2, rep(NA, missinglength)))
  }
  if (length(eststage1) < fulllength) {
    missinglength <- fulllength - length(eststage1)
    eststage1 <- as.character(append(eststage1, rep(NA, missinglength)))
  }
  if (length(givenrate) < fulllength) {
    missinglength <- fulllength - length(givenrate)
    givenrate <- as.numeric(append(givenrate, rep(NA, missinglength)))
  }
  if (length(type) < fulllength) {
    missinglength <- fulllength - length(type)
    type <- as.character(append(type, rep(NA, missinglength)))
  }
  if (length(type_t12) < fulllength) {
    missinglength <- fulllength - length(type_t12)
    type_t12 <- as.character(append(type_t12, rep(NA, missinglength)))
  }
  
  if(!all(is.na(type))) {
    convtype <- rep(1, length(type))
    convtype[which(type == "F")] <- 2
    convtype[which(type == "f")] <- 2
  } else {
    convtype <- rep(1, length(stage3))
  }
  if(!all(is.na(type_t12))) {
    convtype_t12 <- rep(1, length(type_t12))
    convtype_t12[which(type_t12 == "F")] <- 2
    convtype_t12[which(type_t12 == "f")] <- 2
    convtype_t12[which(type_t12 == "2")] <- 2
  } else {
    convtype_t12 <- rep(1, length(stage3))
  }
  
  fullpack <- cbind.data.frame(stage3, stage2, stage1, eststage3, eststage2,
    eststage1, givenrate, convtype, convtype_t12, stringsAsFactors = FALSE)
  
  return(fullpack)
}

#' Test Overdispersion and Zero Inflation in Size and Fecundity Distributions
#' 
#' Function \code{sf_distrib} takes a historically formatted vertical data as
#' input and tests whether size and fecundity data are dispersed according to a
#' Poisson distribution (where mean = variance), and whether the number of 0s
#' exceeds expectations. This function is now deprecated in favor of function
#' \code{\link{hfv_qc}()}.
#' 
#' @name sf_distrib
#' 
#' @param data A historical vertical data file, which is a data frame of class
#' \code{hfvdata}.
#' @param sizea A vector holding the name or column number of the variables
#' corresponding to primary size in occasions \emph{t}+1 and \emph{t}. Input 
#' only if \code{sizea} is to be tested.
#' @param sizeb A vector holding the name or column number of the variables
#' corresponding to secondary size in occasions \emph{t}+1 and \emph{t}. Input 
#' only if \code{sizeb} is to be tested.
#' @param sizec A vector holding the name or column number of the variables
#' corresponding to tertiary size in occasions \emph{t}+1 and \emph{t}. Input 
#' only if \code{sizec} is to be tested.
#' @param obs3 The name or column number of the variable corresponding to
#' observation status in occasion \emph{t}+1. This should be used if observation
#' status will be used as a vital rate to absorb states of size = 0.
#' @param fec A vector holding the names or column numbers of the variables
#' corresponding to in occasions \emph{t}+1 and \emph{t}. Input only if 
#' \code{fec} is to be tested.
#' @param repst A vector holding the names or column numbers of the variables
#' corresponding to reproductive status in occasions \emph{t}+1 and \emph{t}.
#' If not provided, then fecundity will be tested without subsetting to only
#' reproductive individuals.
#' @param zisizea A logical value indicating whether to conduct a test of zero
#' inflation in primary size. Defaults to \code{TRUE}.
#' @param zisizeb A logical value indicating whether to conduct a test of zero
#' inflation in secondary size. Defaults to \code{TRUE}.
#' @param zisizec A logical value indicating whether to conduct a test of zero
#' inflation in tertiary size. Defaults to \code{TRUE}.
#' @param zifec A logical value indicating whether to conduct a test of zero
#' inflation in fecundity. Defaults to TRUE.
#' @param fectime An integer indicating whether to treat fecundity as occurring
#' in time \emph{t} (\code{2}) or time \emph{t}+1 (\code{3}). Defaults to
#' \code{2}.
#' @param show.size A logical value indicating whether to show the output for
#' tests of size. Defaults to \code{TRUE}.
#' @param show.fec A logical value indicating whether to show the output for
#' tests of fecundity. Defaults to \code{TRUE}.
#'
#' @return Produces text describing the degree and significance of difference
#' from expected dispersion, and the degree and significance of zero inflation.
#' The tests are chi-squared score tests based on the expectations of 
#' mean = variance, and 0s as abundant as predicted by the value of lambda
#' estimated from the dataset. See van der Broek (1995) for more details.
#' 
#' @section Notes:
#' This function subsets the data in the same way as \code{\link{modelsearch}()}
#' before testing underlying distributions, making the output much more
#' appropriate than a simple analysis of size and fecundity variables in
#' \code{data}.
#' 
#' The specific test used for overdispersion is a chi-squared test of the
#' dispersion parameter estimated using a generalized linear model predicting
#' the response given size in occasion \emph{t}, under a quasi-Poisson
#' distribution.
#' 
#' The specific test used for zero-inflation is the chi-squared test presented
#' in van der Broek (1995).
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
#'   9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
#'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r",
#'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'   0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' sf_distrib(lathvertln, sizea = c("sizea3", "sizea2"), fec = c("feca3", "feca2"),
#'   repst = c("repstatus3", "repstatus2"), zifec = FALSE)
#' 
#' @export
sf_distrib <- function(data, sizea = NA, sizeb = NA, sizec = NA, obs3 = NA,
  fec = NA, repst = NA, zisizea = TRUE, zisizeb = TRUE, zisizec = TRUE,
  zifec = TRUE, fectime = 2, show.size = TRUE, show.fec = TRUE) {
  
  alive3 <- NULL
  
  if (!is(data, "hfvdata")) {
    stop("Function sf_distrib requires an object of class hfvdata as input.",
      call. = FALSE)
  }
  
  if (all(is.na(c(sizea, sizeb, sizec)))) {
    stop("Function sf_distrib requires a size variable, even if only fecundity is
      to be tested. Please designate at least one such variable", call. = FALSE)
  }
  
  if (!all(is.na(fec)) & all(is.na(repst))) {
    warning("Fecundity will be assessed without subsetting to only reproductive
      individuals. If incorrect, then please provide reproductive status
      variables (repst) in time t+1 and t.", call. = FALSE)
  }
  
  sdata <- subset(data, alive3 == 1)
  
  if (length(obs3) > 1) {
    stop("Please input only one column number or name in obs3.", call. = FALSE)
  }
  if (!all(is.na(obs3))) {
    if (is.numeric(obs3)) {
      if (!is.element(obs3, c(1:length(names(data))))) {
        stop("Column number given in obs is not recognized.", call. = FALSE)
      } else {
        sdata <- sdata[which(sdata[,obs3] == 1),]
      }
    } else if (is.character(obs3)) {
      obs3low <- tolower(obs3)
      datanames <- tolower(names(sdata))
      
      obs3proxy <- grep(obs3low, datanames, fixed = TRUE)
      
      if (length(obs3proxy) == 0) {
        stop("Name of obs3 variable does not match any variable in dataset.",
          call. = FALSE)
      } else if (length(obs3proxy) > 1) {
        stop("Variable obs3 appears to match several variable names in dataset.",
          call. = FALSE)
      }
      
      sdata <- sdata[which(sdata[,obs3proxy] == 1),]
    }
  }
  
  sizeatest <- .knightswhosaynee(sdata, sizea, 1, zisizea, sizea, sizeb, sizec,
    repst, fectime, show.size)
  sizebtest <- .knightswhosaynee(sdata, sizeb, 2, zisizeb, sizea, sizeb, sizec,
    repst, fectime, show.size)
  sizectest <- .knightswhosaynee(sdata, sizec, 3, zisizec, sizea, sizeb, sizec,
    repst, fectime, show.size)
  
  fectest <- .knightswhosaynee(sdata, fec, 4, zifec, sizea, sizeb, sizec, repst,
    fectime, show.fec)
}

#' Internal Test of Size and Fecundity For Overdispersion and Zero-Inflation
#' 
#' This internal function automates the testing of size and fecundity variables
#' for overdispersion and zero-inflation. It is used within
#' \code{\link{sf_distrib}()}.
#' 
#' @name .knightswhosaynee
#' 
#' @param data_used The modified data frame to be used, typically \code{sdata}.
#' @param variable_vec A vector giving the names or numbers of the size or
#' fecundity variables to be tested in occasions \emph{t}+1 and \emph{t},
#' respectively.
#' @param term_used An integer denoting which variable is being worked with. If
#' \code{1}, then \code{sizea}; if \code{2}, then \code{sizeb}, if \code{3},
#' then \code{sizec}, and if \code{4}, then \code{fec}.
#' @param zi_state A logical value indicating whether to test for
#' zero-inflation.
#' @param size_a The vector of variable names or column numbers denoting
#' primary size (\code{sizea}). Used in fecundity assessment.
#' @param size_b The vector of variable names or column numbers denoting
#' secondary size (\code{sizeb}). Used in fecundity assessment.
#' @param size_c The vector of variable names or column numbers denoting
#' tertiary size (\code{sizec}). Used in fecundity assessment.
#' @param repst The vector of variable names or column numbers denoting
#' reproductive status (\code{repst}). Used in fecundity assessment.
#' @param fectime An integer denoting whether fecundity is assessed in time
#' \emph{t}+1 (\code{3}) or time \emph{t} (\code{2}). Used in fecundity
#' assessment.
#' @param show_var A logical variable indicating whether to show the results of
#' tests for the particular variable in question.
#' 
#' @return This function produces text in the console giving the results of the
#' tests of overdispersion and zero inflation. No specific object is returned.
#' 
#' @section Notes:
#' This function will not test for overdispersion and zero inflation in
#' non-integer variables.
#' 
#' @keywords internal
#' @noRd
.knightswhosaynee <- function(data_used, variable_vec, term_used, zi_state,
  size_a, size_b, size_c, repst, fectime, show_var) {
  
  var_used <- FALSE
  
  full_term <- c("sizea", "sizeb", "sizec", "fec")
  full_inenglish <- c("Primary size", "Secondary size", "Tertiary size", "Fecundity")
  full_inenglish_small <- tolower(full_inenglish)
  
  if (!all(is.na(variable_vec))) {
    if (length(variable_vec) < 2) {
      
      stop(paste0("Function sf_distrib requires ", full_term[term_used], " in occasions t and t+1
        to function. Please designate these variables within a 2-element vector as input for",
        full_term[term_used]), call. = FALSE)
    } else if (length(variable_vec) > 2) {
      if (any(is.na(variable_vec[c(1,2)]))) {
        stop(paste0("NAs cannot be included within the first two elements of ", full_term[term_used], "."),
          call. = FALSE)
      }
      warning(paste0("Only the first two entries will be used in option ", full_term[term_used],
        " corresponding to occasions t+1 and t, respectively."), call. = FALSE)
    }
    if (all(is.numeric(variable_vec))) {
      if (any(!is.element(variable_vec, c(1:length(names(data_used)))))) {
        stop(paste0("Column numbers given in ", full_term[term_used]," are not recognized."), call. = FALSE)
      }
    } else if (all(is.character(variable_vec))) {
      if (any(!is.element(variable_vec, names(data_used)))) {
        stop(paste0("Column names given in ", full_term[term_used]," are not recognized."), call. = FALSE)
      }
    }
    
    var_used <- TRUE
    if (term_used < 4) {
      var3data <- data_used[, variable_vec[1]]
      var2data <- data_used[, variable_vec[2]]
    } else {
      if (!all(is.na(size_a)) & length(size_a) == 2) {
        var2 <- size_a
      } else if (!all(is.na(size_b)) & length(size_b) == 2) {
        var2 <- size_b
      } else if (!all(is.na(size_c)) & length(size_c) == 2) {
        var2 <- size_c
      }
      
      if (!all(is.na(repst))) {
        if (length(repst) < 2) {
          stop("Fecundity cannot be tested without repst holding the names of the reproductive status
            variables for times t+1 and t.", call. = FALSE)
        }
        if (fectime == 3) {
          repst_data <- data_used[which(data_used[,repst[1]] == 1),]
          var3data <- repst_data[, variable_vec[1]]
          var2data <- repst_data[, var2[1]]
        } else if (fectime == 2) {
          repst_data <- data_used[which(data_used[,repst[2]] == 1),]
          var3data <- repst_data[, variable_vec[2]]
          var2data <- repst_data[, var2[2]]
        } else {
          stop("Timing of fecundity not recognized.", call. = FALSE)
        }
        
      } else {
        if (fectime == 3) {
          var3data <- data_used[, variable_vec[1]]
          var2data <- data_used[, var2[1]]
        } else {
          var3data <- data_used[, variable_vec[2]]
          var2data <- data_used[, var2[2]]
        }
      }
    }
  
    #Here is the test of overdispersion
    jvmean <- mean(var3data)
    jvvar <- stats::var(var3data)
    
    v_pmodel <- stats::glm(var3data ~ var2data)
    v_qpmodel <- stats::glm(var3data ~ var2data)
    v_disp <- summary(v_qpmodel)$dispersion
    v_df <- summary(v_pmodel)$df.residual
    
    jvodchip <- stats::pchisq(v_disp * v_df, v_df, lower = FALSE)
    
    jvintcheck <- var3data%%1
    if (length(which(jvintcheck != 0)) > 0) {
      writeLines(paste0("Non-integer values detected, so will not test for overdispersion and zero-inflation in ",
          full_term[term_used]))
      show_var <- FALSE
      zi_state <- FALSE
    }
    
    if (show_var) {
      writeLines(paste0("Mean ", full_term[term_used]," is ", signif(jvmean, digits = 4)))
      writeLines(paste0("\nThe variance in ", full_term[term_used]," is ", signif(jvvar, digits = 4)))
      writeLines("\nThe probability of this dispersion level by chance assuming that")
      writeLines(paste0("the true mean ", full_term[term_used]," = variance in ", full_term[term_used], ","))
      writeLines(paste0("and an alternative hypothesis of overdispersion, is ", signif(jvodchip, digits = 4)))
      
      if (jvodchip <= 0.05 & jvvar > jvmean) {
        writeLines(paste0("\n", full_inenglish[term_used]," is significantly overdispersed.\n"))
      } else if (jvodchip <= 0.05 & jvvar < jvmean) {
        writeLines(paste0("\n", full_inenglish[term_used]," is significantly underdispersed.\n"))
      } else {
        writeLines(paste0("\nDispersion level in ", full_inenglish_small[term_used]," matches expectation.\n"))
      }
    }
    
    #Here is the test of zero inflation
    if (zi_state) {
      
      v0est <- exp(-jvmean) #Estimated lambda
      v0n0 <- sum(var3data == 0) #Actual no of zeroes
      
      v0exp <- length(var3data) * v0est #Expected no of zeroes
      
      jvdbs <- (v0n0 - v0exp)^2 / (v0exp * (1 - v0est) - length(var3data) * jvmean * (v0est^2))
      jvzichip <- stats::pchisq(jvdbs, df = 1, lower.tail = FALSE)
      
      if (v0n0 < v0exp & jvzichip < 0.50) { #Correction for lower than expected numbers of 0s
        jvzichip <- 1 - jvzichip
      }
      
      if (show_var) {
        writeLines(paste0("\nMean lambda in ", full_term[term_used]," is ", signif(v0est, digits = 4)))
        writeLines(paste0("The actual number of 0s in ", full_term[term_used]," is ", v0n0))
        writeLines(paste0("The expected number of 0s in ", full_term[term_used]," under the null hypothesis is ", signif(v0exp, digits = 4)))
        writeLines(paste0("The probability of this deviation in 0s from expectation by chance is ", signif(jvzichip, digits = 4)))
        
        if (jvzichip <= 0.05 & v0n0 > v0exp) {
          writeLines(paste0("\n", full_inenglish[term_used]," is significantly zero-inflated.\n"))
        } else {
          writeLines(paste0("\n", full_inenglish[term_used]," is not significantly zero-inflated.\n"))
          
          if (v0n0 == 0) {
            writeLines(paste0(full_inenglish[term_used]," does not appear to include 0s, suggesting
              that a zero-truncated distribution may be warranted."))
          }
          writeLines("\n")
        }
        
        writeLines("\n\n")
      }
    }
  }
}

#' Set Density Dependence Relationships in Vital Rates
#' 
#' Function \code{density_vr()} provides all necessary data to incorporate
#' density dependence into the vital rate functions used to create matrices in
#' function-based projections using function \code{f_projection3()}. Four forms
#' of density dependence are allowed, including the Ricker function, the
#' Beverton-Holt function, the Usher function, and the logistic function. In
#' each case, density must have an effect with at least a one time-step delay
#' (see Notes).
#'
#' @name density_vr
#' 
#' @param density_yn A 14 element logical vector denoting whether each vital
#' rate is subject to density dependence. The order of vital rates is: survival
#' probability, observation probability, primary size transition, secondary size
#' transition, tertiary size transition, reproductive status probability,
#' fecundity rate, juvenile survival probability, juvenile observation
#' probability, juvenile primary size transition, juvenile secondary size
#' transition, juvenile tertiary size transition, juvenile reproductive status
#' probability, and juvenile maturity status probability. Defaults to a vector
#' of 14 \code{FALSE} values.
#' @param style A 14 element vector coding for the style of density dependence
#' on each vital rate. Options include \code{0}: no density dependence,
#' \code{1}, \code{ricker}, \code{ric}, or \code{r} for the Ricker function;
#' \code{2}, \code{beverton}, \code{bev}, and \code{b} for the Beverton-Holt
#' function; \code{3}, \code{usher}, \code{ush}, and \code{u} for the Usher
#' function; and \code{4}, \code{logistic}, \code{log}, and \code{l} for the
#' logistic function. Defaults to 14 values of \code{0}.
#' @param time_delay A 14 element vector indicating the number of occasions back
#' on which density dependence operates. Defaults to 14 values of \code{1}, and
#' may not include any number less than 1.
#' @param alpha A 14 element vector indicating the numeric values to use as the
#' alpha term in the two parameter Ricker, Beverton-Holt, or Usher function, or
#' the value of the carrying capacity \emph{K} to use in the logistic equation
#' (see \code{Notes} for more on this term). Defaults to 14 values of \code{0}.
#' @param beta A 14 element vector indicating the numeric values to use as the
#' beta term in the two parameter Ricker, Beverton-Holt, or Usher function. Used
#' to indicate whether to use \emph{K} as a hard limit in the logistic equation
#' (see \code{Notes} below). Defaults to 14 values of \code{0}.
#' 
#' @return A data frame of class \code{lefkoDensVR} with 14 rows, one for each
#' vital rate in the order of: survival probability, observation probability,
#' primary size transition, secondary size transition, tertiary size transition,
#' reproductive status probability, fecundity rate, juvenile survival
#' probability, juvenile observation probability, juvenile primary size
#' transition, juvenile secondary size transition, juvenile tertiary size
#' transition, juvenile reproductive status probability, and juvenile maturity
#' status probability. This object can be used as input in function
#' \code{\link{f_projection3}()}.
#' 
#' Variables in this object include the following:
#' \item{vital_rate}{The vital rate to be modified.}
#' \item{density_yn}{Logical value indicating whether vital rate will be subject
#' to density dependence.}
#' \item{style}{Style of density dependence, coded as \code{1}, \code{2},
#' \code{3}, \code{4}, or \code{0} for the Ricker, Beverton-Holt, Usher, or
#' logistic function, or no density dependence, respectively.}
#' \item{time_delay}{The time delay on density dependence, in time steps.}
#' \item{alpha}{The value of alpha in the Ricker, Beverton-Holt, or Usher
#' function, or the value of carrying capacity, \emph{K}, in the logistic
#' function.}
#' \item{beta}{The value of beta in the Ricker, Beverton-Holt, or Usher
#' function.}
#' 
#' @section Notes:
#' This function provides inputs when density dependence is operationalized
#' directly on vital rates. It can be used only in function
#' \code{f_projection3()}. Users wishing to modify matrix elements directly by
#' density dependence functions for use in function-based or raw projections
#' with functions \code{projection3()} and \code{f_projection3()} should use
#' function \code{density_input()} to provide the correct inputs.
#' 
#' The parameters \code{alpha} and \code{beta} are applied according to the
#' two-parameter Ricker function, the two-parameter Beverton-Holt function, the
#' two-parameter Usher function, or the one-parameter logistic function.
#' Although the default is that a 1 time step delay is assumed, greater time
#' delays can be set through the \code{time_delay} option.
#' 
#' When using the logistic function, it is possible that the time delay used in
#' density dependent simulations will cause matrix elements to become negative.
#' To prevent this behavior, set the associated \code{beta} term to \code{1.0}.
#' Doing so will set \code{K} as the hard limit in the logistic equation,
#' essentially setting a minimum limit at \code{0} for all matrix elements
#' modified.
#' 
#' @seealso \code{\link{density_input}()}
#' @seealso \code{\link{f_projection3}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
#'   9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
#'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
#'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#'   0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
#'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
#'   propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
#'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988", 
#'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathvertln_adults <- subset(lathvertln, stage2index > 2)
#' surv_model <- glm(alive3 ~ sizea2 + sizea1 + as.factor(patchid) +
#'   as.factor(year2), data = lathvertln_adults, family = "binomial")
#' 
#' obs_data <- subset(lathvertln_adults, alive3 == 1)
#' obs_model <- glm(obsstatus3 ~ as.factor(patchid), data = obs_data,
#'   family = "binomial")
#' 
#' size_data <- subset(obs_data, obsstatus3 == 1)
#' siz_model <- lm(sizea3 ~ sizea2 + sizea1 + repstatus1 + as.factor(patchid) +
#'   as.factor(year2), data = size_data)
#' 
#' reps_model <- glm(repstatus3 ~ sizea2 + sizea1 + as.factor(patchid) +
#'   as.factor(year2), data = size_data, family = "binomial")
#' 
#' fec_data <- subset(lathvertln_adults, repstatus2 == 1)
#' fec_model <- glm(feca2 ~ sizea2 + sizea1 + repstatus1 + as.factor(patchid),
#'   data = fec_data, family = "poisson")
#' 
#' lathvertln_juvs <- subset(lathvertln, stage2index < 3)
#' jsurv_model <- glm(alive3 ~ as.factor(patchid), data = lathvertln_juvs,
#'   family = "binomial")
#' 
#' jobs_data <- subset(lathvertln_juvs, alive3 == 1)
#' jobs_model <- glm(obsstatus3 ~ 1, family = "binomial", data = jobs_data)
#' 
#' jsize_data <- subset(jobs_data, obsstatus3 == 1)
#' jsiz_model <- lm(sizea3 ~ as.factor(year2), data = jsize_data)
#' 
#' jrepst_model <- 0
#' jmatst_model <- 1
#' 
#' mod_params <- create_pm(name_terms = TRUE)
#' mod_params$modelparams[3] <- "patchid"
#' mod_params$modelparams[4] <- "alive3"
#' mod_params$modelparams[5] <- "obsstatus3"
#' mod_params$modelparams[6] <- "sizea3"
#' mod_params$modelparams[9] <- "repstatus3"
#' mod_params$modelparams[11] <- "feca2"
#' mod_params$modelparams[12] <- "sizea2"
#' mod_params$modelparams[13] <- "sizea1"
#' mod_params$modelparams[18] <- "repstatus2"
#' mod_params$modelparams[19] <- "repstatus1"
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "mat", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "Sdl", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "Sd", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, "mat", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, "Sdl", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, "Sdl", NA, NA),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframeln, historical = TRUE)
#' 
#' # While we do not use MPMs to initialize f_projections3(), we do use MPMs to
#' # initialize functions start_input() and density_input().
#' lathmat3ln <- flefko3(year = "all", patch = "all", data = lathvertln,
#'   stageframe = lathframeln, supplement = lathsupp3, paramnames = mod_params,
#'   surv_model = surv_model, obs_model = obs_model, size_model = siz_model,
#'   repst_model = reps_model, fec_model = fec_model, jsurv_model = jsurv_model,
#'   jobs_model = jobs_model, jsize_model = jsiz_model,
#'   jrepst_model = jrepst_model, jmatst_model = jmatst_model, reduce = FALSE)
#' 
#' e3m_sv <- start_input(lathmat3ln, stage2 = "Sd", stage1 = "Sd", value = 1000)
#' 
#' dyn7 <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#'   FALSE, FALSE, FALSE, FALSE, FALSE)
#' dst7 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' dal7 <- c(0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' dbe7 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' 
#' e3d_vr <- density_vr(density_yn = dyn7, style = dst7, alpha = dal7,
#'   beta = dbe7)
#' 
#' trial7_dvr_1 <- f_projection3(format = 1, data = lathvertln, supplement = lathsupp3,
#'   paramnames = mod_params, stageframe = lathframeln, nreps = 2,
#'   surv_model = surv_model, obs_model = obs_model, size_model = siz_model,
#'   repst_model = reps_model, fec_model = fec_model, jsurv_model = jsurv_model,
#'   jobs_model = jobs_model, jsize_model = jsiz_model,
#'   jrepst_model = jrepst_model, jmatst_model = jmatst_model,
#'   times = 100, stochastic = TRUE, standardize = FALSE, growthonly = TRUE,
#'   integeronly = FALSE, substoch = 0, sp_density = 0, start_frame = e3m_sv,
#'   density_vr = e3d_vr)
#' }
#' 
#' @export
density_vr <- function(density_yn = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  style = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  time_delay = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  alpha = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  beta = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) {
  
  full_length <- NULL
  
  if (!all(is.logical(density_yn))) {
    if (all(is.numeric(density_yn))) {
      if (any(!is.element(density_yn, c(0, 1)))) {
        stop("Option density_yn must be a 14 element vector of TRUE and FALSE values.",
          call. = FALSE)
      } else {
        density_yn <- as.logical(density_yn)
      }
    } else {
      stop("Option density_yn must be a 14 element vector of TRUE and FALSE values.",
        call. = FALSE)
    }
  }
  
  if (all(is.character(style))) {
    style <- tolower(style)
    
    ricker_style <- c("1", "ricker", "ricke", "rick", "ric", "ri", "r")
    beverton_style <- c("2", "beverton", "beverto", "bevert", "bever", "beve",
      "bev", "be", "b", "holt", "hol", "ho", "h")
    usher_style <- c("3", "usher", "ushe", "ush", "us", "u")
    logistic_style <- c("4", "logistic", "logisti", "logist", "logis", "logi",
      "log", "lo", "l")
    no_style <- c("0")
    
    unknown_style <- which(!is.element(style, c(ricker_style, beverton_style,
        logistic_style, no_style)))
    if (length(unknown_style) > 0) {
      stop(paste0("Unknown style code used: ", style[unknown_style], ". Cannot
        process."), call. = FALSE)
    }
    
    lstyle <- rep(0, full_length)
    lstyle[which(is.element(style, ricker_style))] <- 1
    lstyle[which(is.element(style, beverton_style))] <- 2
    lstyle[which(is.element(style, usher_style))] <- 3
    lstyle[which(is.element(style, logistic_style))] <- 4
    
    style <- lstyle
  } else if (all(is.numeric(style))) {
    if (any(style > 4) | any(style < 0) | any(style %% 1 > 0)) {
      stop("Unknown style code used. Please only use numbers 0, 1, 2, 3, or 4.",
        call. = FALSE)
    }
  } else {
    stop("Input style codes do not conform to accepted inputs.", call. = FALSE)
  }
  
  if (!all(as.integer(time_delay) == time_delay)) {
    stop("Input for time_delay must be integer.", call. = FALSE)
  } else if (any(time_delay < 1)) {
    stop("Input for time_delay must be an integer no smaller than 1.",
      call. = FALSE)
  }
  
  if (any(!is.numeric(alpha)) & !all(is.na(alpha))) {
    stop("Option alpha must be either NA or a numeric value.", call. = FALSE)
  }
  if (any(!is.numeric(beta)) & !all(is.na(beta))) {
    stop("Option beta must be either NA or a numeric value.", call. = FALSE)
  }
  
  if (length(density_yn) != 14) {
    stop("Vector density_yn must be 14 elements long.",
      call. = FALSE)
  }
  if (length(alpha) != 14) {
    stop("Vector alpha must be 14 elements long.",
      call. = FALSE)
  }
  if (length(beta) != 14) {    
    stop("Vector beta must be 14 elements long.",
      call. = FALSE)
  }
  if (length(style) != 14) {
    stop("Vector style must be 14 elements long.",
      call. = FALSE)
  }
  if (length(time_delay) != 14) {
    stop("Vector time_delay must be 14 elements long.",
      call. = FALSE)
  }
  
  vital_rate <- c("survival", "observation", "size", "sizeb", "sizec",
    "repstatus", "fecundity", "juv_survival", "juv_observation", "juv_size",
    "juv_sizeb", "juv_sizec", "juv_reproduction", "juv_maturity")
  
  output <- cbind.data.frame(vital_rate, density_yn, style, time_delay, alpha,
    beta, stringsAsFactors = FALSE)
  
  class(output) <- append(class(output), "lefkoDensVR")
  
  return(output)
}

#' Create a Starting Vector for Population Projection
#' 
#' Function \code{start_input()} creates a data frame summarizing the non-zero
#' elements of the start vector for use in population projection analysis via
#' function \code{\link{projection3}()}.
#' 
#' @name start_input
#' 
#' @param mpm The lefkoMat object to be used in projection analysis.
#' @param stage2 A vector showing the name or number of a stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Abbreviations for groups of stages are also usable (see Notes).
#' This input is required for all stage-based and age-by-stage MPMs. Defaults to
#' \code{NA}.
#' @param stage1 A vector showing the name or number of a stage in occasion
#' \emph{t}-1 that should be set to a positive number of individuals in the
#' start vector. Abbreviations for groups of stages are also usable (see Notes).
#' This is only used for historical MPMs, since the rows of hMPMs correspond to
#' stage-pairs in times \emph{t} and \emph{t}-1 together. Only required for
#' historical MPMs, and will result in errors if otherwise used.
#' @param age2 A vector showing the age of each respective stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Only used for Leslie and age-by-stage MPMs. Defaults to \code{NA}.
#' @param value A vector showing the values, in order, of the number of
#' individuals set for the stage or stage-pair in question. Defaults to 1.
#' 
#' @return A list of class \code{lefkoSV}, with four objects, which can be used
#' as input in function \code{\link{projection3}()}. The last three include the
#' \code{ahstages}, \code{hstages}, and \code{agestages} objects from the
#' \code{lefkoMat} object supplied in \code{mpm}. The first element in the list
#' is a data frame with the following variables:
#' 
#' \item{stage2}{Stage at occasion \emph{t}.}
#' \item{stage_id_2}{The stage number associated with \code{stage2}.}
#' \item{stage1}{Stage at occasion \emph{t}-1, if historical. Otherwise NA.}
#' \item{stage_id_1}{The stage number associated with \code{stage1}.}
#' \item{age2}{The age of individuals in \code{stage2} and, if applicable,
#' \code{stage1}. Only used in age-by-stage MPMs.}
#' \item{row_num}{A number indicating the respective starting vector element.}
#' \item{value}{Number of individuals in corresponding stage or stage-pair.}
#' 
#' @section Notes:
#' Entries in \code{stage2}, and \code{stage1} can include abbreviations for
#' groups of stages. Use \code{rep} if all reproductive stages are to be used,
#' \code{nrep} if all mature but non-reproductive stages are to be used,
#' \code{mat} if all mature stages are to be used, \code{immat} if all immature
#' stages are to be used, \code{prop} if all propagule stages are to be used,
#' \code{npr} if all non-propagule stages are to be used, \code{obs} if all
#' observable stages are to be used, \code{nobs} if all unobservable stages are
#' to be used, and leave empty or use \code{all} if all stages in stageframe are
#' to be used.
#' 
#' @seealso \code{\link{density_input}()}
#' @seealso \code{\link{projection3}()}
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
#' 
#' e3m_sv <- start_input(ehrlen3mean, stage2 = "Sd", stage1 = "Sd", value = 1000)
#' 
#' lathproj <- projection3(ehrlen3, nreps = 5, times = 100, stochastic = TRUE,
#'   start_frame = e3m_sv)
#' 
#' @export
start_input <- function(mpm, stage2 = NA, stage1 = NA, age2 = NA, value = 1) {
  
  mpmrows <- stage2_id <- stage1_id <- start_vec <- full_length <- NULL
  
  if (all(!is(mpm, "lefkoMat"))) {
    stop("A regular lefkoMat object is required as input.", call. = FALSE)
  }
  
  if (all(is.na(stage2)) & all(is.na(age2))) {
    stop("Options stage2 and age2 cannot both be set to NA.", call. = FALSE)
  }
  if (all(is.null(stage2)) & all(is.null(age2))) {
    stop("Options stage2 and age2 cannot both be empty.", call. = FALSE)
  }
  
  if (!is.element("stage", names(mpm$ahstages))) {
    stop("Stageframe appears to be modified. Please make sure that a stage
      column exists holding stage names.", call. = FALSE)
  }
  
  if (all(is.na(mpm$hstages)) | all(is.null(mpm$hstages))) {
    historical <- FALSE
  } else {
    historical <- TRUE
  }
  if (all(is.na(mpm$agestages)) | all(is.null(mpm$agestages))) {
    agebystage <- FALSE
  } else {
    agebystage <- TRUE
  }
  
  if (historical & all(is.na(stage1))) {
    stop("Historical projection analysis requires that stage in time t-1 be
      designated for all stage pairs.", call. = FALSE)
  } else if (!historical & !all(is.na(stage1))) {
    stop("Ahistorical projection analysis cannot include designated stages in
      time t-1.", call. = FALSE)
  }
  
  full_length <- max(length(stage2), length(stage1), length(age2), length(value))
  
  if (length(value) == 1 & full_length > 1) {
    value <- rep(value, full_length)
  }
  
  if (length(stage2) == 1 & full_length > 1) {
    stage2 <- rep(stage2, full_length)
  }
  
  if (length(stage1) == 1 & full_length > 1) {
    stage1 <- rep(stage1, full_length)
  }
  
  if (length(age2) == 1 & full_length > 1) {
    age2 <- rep(age2, full_length)
  }
  
  if ((all(is.na(stage2)) | all(is.null(stage2))) & (all(is.na(age2)) | all(is.null(age2)))) {
    stop("Either stage2 or age2 must be provided.", call. = FALSE)
  }
  
  if (all(is.character(stage2))) {
    unknown_stage2 <- which(!is.element(tolower(stage2), c(tolower(mpm$ahstages$stage),
        c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "obs", "nobs"))))
    if (length(unknown_stage2) > 0) {
      stop(paste0("Unknown stage designations used in stage2: ",
        stage2[unknown_stage2]), call. = FALSE)
    }
    
    reassessed <- apply(as.matrix(c(1:length(stage2))), 1, function(X) {
      if (!is.na(stage2[X])) {
        if (is.element(stage2[X], mpm$ahstages$stage)) {
          shrubbery.small <- cbind.data.frame(stage2 = stage2[X], stage1 = stage1[X],
            age2 = age2[X], value = value[X], stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (is.element(stage2[X], as.character(mpm$ahstages$stage_id))) {
          shrubbery.small <- cbind.data.frame(stage2 = as.numeric(stage2[X]),
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "rep") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "all") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage,
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nrep") {
          shrubbery.small <- cbind.data.frame(
            stage2 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                which(mpm$ahstages$matstatus == 1))],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) =="mat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "immat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "prop") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "npr") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "obs") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nobs") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        }
      }
    })
    shrubbery <- do.call(rbind.data.frame, reassessed)
    
  } else if (all(is.numeric(stage2)) & !any(is.na(stage2))) {
    stage2_id <- stage2
    
    if (any(stage2_id > max(mpm$ahstages$stage_id)) | any(stage2_id < min(mpm$ahstages$stage_id))) {
      stop("Unknown stage2 codes used.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(stage2_id), 1, function(X) {
      return(mpm$ahstages$stage[X])
    })
    
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
    
  } else if (all(is.numeric(age2)) & !any(is.na(age2))) {
    if (!(all(is.na(stage2)))) {
      stop("Leslie MPMs should not be entered with the stage2 option set to
        values other than NA.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(age2), 1, function(X) {
      cross_ref <- which(mpm$ahstages$stage_id == X)
      return(mpm$ahstages$stage[cross_ref])
    })
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
    
  } else {
    stop("Input stage2 codes do not conform to accepted inputs.", call. = FALSE)
  }
  
  if (historical) {
    if (all(is.character(shrubbery$stage1)) & !all(is.na(shrubbery$stage1))) {
      
      unknown_stage1 <- which(!is.element(tolower(stage1), c(tolower(mpm$ahstages$stage),
          c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "obs", "nobs", "almostborn"))))
      if (length(unknown_stage1) > 0) {
      stop(paste0("Unknown stage designations used in stage1: ",
        stage1[unknown_stage1]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$stage1[X])) {
          if (is.element(shrubbery$stage1[X], mpm$ahstages$stage)) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (is.element(stage1[X], as.character(mpm$ahstages$stage_id))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = as.numeric(shrubbery$stage1[X]), age2 = age2[X],
              value = value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "rep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage, age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nrep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                  which(mpm$ahstages$matstatus == 1))],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) =="mat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "immat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "prop") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "npr") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "obs") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nobs") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      })
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
      
    } else if (all(is.numeric(shrubbery$stage1)) & !any(is.na(shrubbery$stage1))) {
      stage1_id <- shrubbery$stage1
      
      if (any(stage1_id > max(mpm$ahstages$stage_id)) | any(stage1_id < min(mpm$ahstages$stage_id))) {
        stop("Unknown stage1 codes used.", call. = FALSE)
      }
      
      stage1 <- apply(as.matrix(stage1_id), 1, function(X) {
        return(mpm$ahstages$stage[X])
      })
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (agebystage) {
    if (all(is.character(shrubbery$age2)) & !all(is.na(shrubbery$age2))) {
      
      unknown_age2 <- which(!is.element(tolower(age2), c(tolower(mpm$agestages$age),
          "all")))
      if (length(unknown_age2) > 0) {
        stop(paste0("Unknown age designations used in age2: ",
          stage1[unknown_age2]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$age2[X])) {
          common_ages <- unique(mpm$agestages$age[which(mpm$agestages$stage == shrubbery$stage2[X])])
          
          if (is.element(shrubbery$age2[X], as.character(common_ages))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = as.numeric(shrubbery$age2[X]),
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$age2[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = common_ages,
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      }) 
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
    } else if (all(is.numeric(shrubbery$age2)) & !any(is.na(shrubbery$age2))) {
      if (any(shrubbery$age2 > max(mpm$agestages$age)) | any(shrubbery$age2 < min(mpm$agestages$age))) {
        stop("Unknown age2 values used.", call. = FALSE)
      }
      
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = shrubbery$stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (!all(is.numeric(value))) {
    stop("Object value must be composed only of valid numbers.", call. = FALSE)
  }
  
  shrubbery$stage2_id <- apply(as.matrix(shrubbery$stage2), 1, function(X) {
    return(mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)])
  })
  shrubbery$stage1_id <- apply(as.matrix(shrubbery$stage1), 1, function(X) {
    possible_option <- mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)]
    if (length(possible_option) > 0) return(possible_option) else return(NA)
  })
  
  full_length <- dim(shrubbery)[1]
  
  if (!historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$ahstages)[1]) {
      stop("This ahistorical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- shrubbery$stage2_id
    
  } else if (agebystage & !historical) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$agestages)[1]) {
      stop("This age-by-stage mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    if (any(is.na(shrubbery$age2)) | any(!is.numeric(shrubbery$age2))) {
      stop("Option age2 must include only numbers for age-by-stage MPMs.",
        call. = FALSE)
    }
    
    if (any(shrubbery$age2 < min(mpm$agestages$age)) | any(shrubbery$age2 > max(mpm$agestages$age))) {
      stop("Option age2 can only take ages shown in element agestages within the input MPM.",
        call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$agestages$stage_id == shrubbery$stage2_id[X])
      vec1 <- which(mpm$agestages$age == shrubbery$age2[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else if (historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$hstages)[1]) {
      stop("This historical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$hstages$stage_id_2 == shrubbery$stage2_id[X])
      vec1 <- which(mpm$hstages$stage_id_1 == shrubbery$stage1_id[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else {
    stop("Format of mpm not recognized.", call. = FALSE)
  }
  
  output_tab <- cbind.data.frame(shrubbery$stage2, shrubbery$stage2_id,
    shrubbery$stage1, shrubbery$stage1_id, shrubbery$age2, start_vec,
    shrubbery$value, stringsAsFactors = FALSE)
  
  names(output_tab) <- c("stage2", "stage_id_2", "stage1", "stage_id_1", "age2",
    "row_num", "value")
  
  if (!historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2")])
  } else if (historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "stage1", "stage_id_1")])
  } else if (!historical & agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "age2")])
  }
  if (dim(out_check)[1] < dim(output_tab)[1]) {
    warning("Some stages, stage-pairs, or age-stages appear to be listed
      multiple times. This may cause errors in analysis.", call. = FALSE)
  }

  class(output_tab) <- append(class(output_tab), "lefkoSV")
  
  return(output_tab)
}

#' Function to Test Integer and Binomial Status of Vectors and Output Results
#' 
#' Function \code{.intbin_check()} tests for integer and binomial status of a
#' specific variable within a data frame.
#' 
#' @name .intbin_check
#' 
#' @param data The data frame to check.
#' @param term The variable to assess within object \code{data}.
#' 
#' @return This function outputs text assessing the integer and binomial status
#' of the variable tested. No values or objects are returned.
#' 
#' @keywords internal
#' @noRd
.intbin_check <- function(data, term) {
  term_col <- which(names(data) == term)
  if (length(term_col) == 0) stop("Unrecognized term in dataset.", call. = FALSE)
  
  check_term <- .integer_test(data[, term_col])
  
  writeLines(paste0("  Variable ", term, " has ", length(which(is.na(data[, term_col]))),
                    " missing values."))
  
  if (check_term == 1) {
    writeLines(paste0("  Variable ", term," may be a floating point variable rather than
      a binomial variable. One element is not an integer.\n"))
  } else if (check_term > 0) {
    writeLines(paste0("  Variable ", term," appears to be a floating point variable rather than
      a binomial variable. ", check_term, " elements are not integers.\n"))
  } else {
    term_int <- as.integer(data[, term_col])
    term_bin_check <- .binomial_test(term_int)
    
    if (term_bin_check == 1) {
      writeLines(paste0("  Variable ", term," does not appear to be a binomial variable.
        One element is not 0 or 1.\n"))
    } else if (term_bin_check > 0) {
      writeLines(paste0("  Variable ", term," does not appear to be a binomial variable. ",
                        term_bin_check, " elements are not 0 or 1.\n"))
    } else {
      writeLines(paste0("  Variable ", term," is a binomial variable.\n"))
    }
  }
}

#' Function to Check Size and Fecundity Variables for Distribution Assumptions
#' 
#' Function to assess the distribution of variables coding size or fecundity.
#' 
#' @name .sf_dist_check
#' 
#' @param data The data frame to check.
#' @param term The variable to assess within object \code{data}.
#' @param term2 A variable to use as an independent factor in overdispersion
#' tests.
#' 
#' @return This function outputs text to check whether the variable being tested
#' conforms to Gaussian, Gamma, Poisson, or negative binomial distribution
#' assumptions. No value or object is returned.
#' 
#' @keywords internal
#' @noRd
.sf_dist_check <- function(data, term, term2) {
  gaus.check <- NULL
  
  term_col <- which(names(data) == term)
  if (length(term_col) == 0) stop("Unrecognized term in dataset.", call. = FALSE)
  
  correct_var <- data[, term_col]
  check_term <- .integer_test(correct_var)
  
  writeLines(paste0("  Variable ", term, " has ", length(which(is.na(correct_var))),
                    " missing values."))
  
  if (check_term == 1) {
    writeLines(paste0("  Variable ", term,
                      " may be a floating point variable rather than a binomial variable."))
    writeLines("  One element is not an integer.\n")
  } else if (check_term > 0) {
    writeLines(paste0("  Variable ", term,
                      " appears to be a floating point variable."))
    writeLines(paste0("  ", check_term, " elements are not integers."))
    writeLines(paste0("  The minimum value of ", term, " is ",
                      signif(min(correct_var, na.rm = TRUE), digits = 4), " and the maximum is ",
                      signif(max(correct_var, na.rm = TRUE), digits = 4), "."))
    writeLines(paste0("  The mean value of ", term, " is ",
                      signif(mean(correct_var, na.rm = TRUE), digits = 4), " and the variance is ",
                      signif(var(correct_var, na.rm = TRUE), digits = 4), "."))
    
    if (length(correct_var) > 2 & length(correct_var) < 5001) {
      gaus.check <- stats::shapiro.test(correct_var)
    } else if (length(correct_var) > 2) {
      used_sample <- sample(correct_var, 5000)
      gaus.check <- stats::shapiro.test(used_sample)
    } else {
      stop("Variable to test Gaussian assumptions for appears to be too short.",
           call. = FALSE)
    }
    
    writeLines(paste0("  The value of the Shapiro-Wilk test of normality is ",
                      signif(gaus.check$statistic, digits = 4), " with P = ",
                      signif(gaus.check$p.value, digits = 4), "."))
    if (gaus.check$p.value <= 0.05) {
      writeLines(paste0("  Variable ", term,
                        " differs significantly from a Gaussian distribution.\n"))
    } else {
      writeLines(paste0("  Variable ", term,
                        " does not differ significantly from a Gaussian distribution.\n"))
    }
  } else {
    no_elems <- length(unique(correct_var))
    
    if (no_elems == 1) {
      writeLines(paste0("  Variable ", term," appears to be constant.\n"))
    } else if (no_elems == 2) {
      writeLines(paste0("  Variable ", term," has only two values, and so appears to be a binary variable.\n"))
    } else {
      writeLines(paste0("  Variable ", term," appears to be an integer variable.\n"))
    }
  }
  
  if (any(correct_var < 0)) {
    writeLines(paste0("  Some elements of ", term, " are negative.\n"))
  } else if (any(correct_var == 0)) {
    writeLines(paste0("  Variable ", term, " is fully non-negative.\n"))
  } else {
    writeLines(paste0("  Variable ", term, " is fully positive, lacking even 0s.\n"))
  }
  
  if (check_term == 0) {
    term2_col <- which(names(data) == term2)
    if (length(term2_col) == 0) stop("Unrecognized term in dataset.", call. = FALSE)
    
    correct_var2 <- data[, term2_col]
    
    # Overdispersion test
    jvmean <- mean(correct_var)
    jvvar <- stats::var(correct_var)
    
    v_pmodel <- stats::glm(correct_var ~ correct_var2)
    v_qpmodel <- stats::glm(correct_var ~ correct_var2)
    v_disp <- summary(v_qpmodel)$dispersion
    v_df <- summary(v_pmodel)$df.residual
    
    jvodchip <- stats::pchisq(v_disp * v_df, v_df, lower = FALSE)
    
    writeLines("  Overdispersion test:")
    
    writeLines(paste0("    Mean ", term," is ", signif(jvmean, digits = 4)))
    writeLines(paste0("    The variance in ", term," is ", signif(jvvar, digits = 4)))
    writeLines("    The probability of this dispersion level by chance assuming that")
    writeLines(paste0("    the true mean ", term," = variance in ", term, ","))
    writeLines(paste0("    and an alternative hypothesis of overdispersion, is ",
                      signif(jvodchip, digits = 4)))
    
    if (jvodchip <= 0.05 & jvvar > jvmean) {
      writeLines(paste0("    Variable ", term," is significantly overdispersed.\n"))
    } else if (jvodchip <= 0.05 & jvvar < jvmean) {
      writeLines(paste0("    Variable ", term," is significantly underdispersed.\n"))
    } else {
      writeLines(paste0("    Dispersion level in ", term," matches expectation.\n"))
    }
    
    # Test for zero-inflation and zero-truncation
    writeLines("  Zero-inflation and truncation tests:")
    
    v0est <- exp(-jvmean) #Estimated lambda
    v0n0 <- sum(correct_var == 0) #Actual no of zeroes
    
    v0exp <- length(correct_var) * v0est #Expected no of zeroes
    
    jvdbs <- (v0n0 - v0exp)^2 / (v0exp * (1 - v0est) - length(correct_var) * jvmean * (v0est^2))
    jvzichip <- stats::pchisq(jvdbs, df = 1, lower.tail = FALSE)
    
    if (v0n0 < v0exp & jvzichip < 0.50) { #Correction for lower than expected numbers of 0s
      jvzichip <- 1 - jvzichip
    }
    
    writeLines(paste0("    Mean lambda in ", term," is ", signif(v0est, digits = 4)))
    writeLines(paste0("    The actual number of 0s in ", term," is ", v0n0))
    writeLines(paste0("    The expected number of 0s in ", term,
                      " under the null hypothesis is ", signif(v0exp, digits = 4)))
    writeLines(paste0("    The probability of this deviation in 0s from expectation by chance is ",
                      signif(jvzichip, digits = 4)))
    
    if (jvzichip <= 0.05 & v0n0 > v0exp) {
      writeLines(paste0("    Variable ", term," is significantly zero-inflated.\n"))
    } else {
      writeLines(paste0("    Variable ", term," is not significantly zero-inflated.\n"))
      
      if (v0n0 == 0 & v0exp >= 1.0) {
        writeLines(paste0("    Variable ", term,
                          " does not include 0s, suggesting that a zero-truncated distribution may be warranted.\n"))
      }
    }
  }
}

#' Assess Quality of hfv Datasets
#' 
#' Function \code{hfv_qc()} tests the overall quality of hfv datasets, and also
#' runs a series of tests to assess which statistical distributions match the
#' variables within these datasets. The input format is equivalent to the input
#' format of function \code{\link{modelsearch}()}, allowing users to assess
#' vital rate variable distributions assuming the same internal dataset
#' subsetting used by the latter function and simply copy and pasting the
#' parameter options from one function to the other.
#' 
#' @name hfv_qc
#' 
#' @param data The vertical dataset to be used for analysis. This dataset should 
#' be of class \code{hfvdata}, but can also be a data frame formatted similarly
#' to the output format provided by functions \code{\link{verticalize3}()} or
#' \code{\link{historicalize3}()}, as long as all needed variables are properly
#' designated.
#' @param stageframe The stageframe characterizing the life history model used.
#' Optional unless \code{test.group = TRUE}, in which case it is required.
#' Defaults to \code{NULL}.
#' @param historical A logical variable denoting whether to assess the effects
#' of state in occasion \emph{t}-1, in addition to state in occasion \emph{t}.
#' Defaults to \code{TRUE}.
#' @param suite This describes the global model for each vital rate estimation,
#' and has the following possible values: \code{full}, includes main effects and
#' all two-way interactions of size and reproductive status; \code{main},
#' includes main effects only of size and reproductive status; \code{size},
#' includes only size (also interactions between size in historical model);
#' \code{rep}, includes only reproductive status (also interactions between
#' status in historical model); \code{age}, all vital rates estimated with age
#' and y-intercepts only; \code{cons}, all vital rates estimated only as
#' y-intercepts. Defaults to \code{size}.
#' @param vitalrates A vector describing which vital rates will be estimated via
#' linear modeling, with the following options: \code{surv}, survival
#' probability; \code{obs}, observation probability; \code{size}, overall size;
#' \code{repst}, probability of reproducing; and \code{fec}, amount of
#' reproduction (overall fecundity). May also be set to
#' \code{vitalrates = "leslie"}, which is equivalent to setting
#' \code{c("surv", "fec")} for a Leslie MPM. This choice also determines how
#' internal data subsetting for vital rate model estimation will work. Defaults
#' to \code{c("surv", "size", "fec")}.
#' @param surv A vector indicating the variable names coding for status as alive
#' or dead in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c("alive3", "alive2", "alive1")}.
#' @param obs A vector indicating the variable names coding for observation
#' status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c("obsstatus3", "obsstatus2", "obsstatus1")}.
#' @param size A vector indicating the variable names coding for the primary
#' size variable on occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
#' respectively. Defaults to \code{c("sizea3", "sizea2", "sizea1")}.
#' @param sizeb A vector indicating the variable names coding for the secondary
#' size variable on occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
#' respectively. Defaults to \code{c(NA, NA, NA)}, in which case \code{sizeb} is
#' not used.
#' @param sizec A vector indicating the variable names coding for the tertiary
#' size variable on occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
#' respectively. Defaults to \code{c(NA, NA, NA)}, in which case \code{sizec} is
#' not used.
#' @param repst A vector indicating the variable names coding for reproductive
#' status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c("repstatus3", "repstatus2", "repstatus1")}.
#' @param fec A vector indicating the variable names coding for fecundity in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to
#' \code{c("feca3", "feca2", "feca1")}.
#' @param stage A vector indicating the variable names coding for stage in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1. Defaults to
#' \code{c("stage3", "stage2", "stage1")}.
#' @param matstat A vector indicating the variable names coding for maturity
#' status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1. Defaults to
#' \code{c("matstatus3", "matstatus2", "matstatus1")}.
#' @param indiv A text value indicating the variable name coding individual
#' identity. Defaults to \code{"individ"}.
#' @param patch A text value indicating the variable name coding for patch,
#' where patches are defined as permanent subgroups within the study population.
#' Defaults to \code{NA}.
#' @param year A text value indicating the variable coding for observation
#' occasion \emph{t}. Defaults to \code{"year2"}.
#' @param density A text value indicating the name of the variable coding for
#' spatial density, should the user wish to test spatial density as a fixed
#' factor affecting vital rates. Defaults to \code{NA}.
#' @param patch.as.random If set to \code{TRUE} and \code{approach = "mixed"},
#' then \code{patch} is included as a random factor. If set to \code{FALSE} and
#' \code{approach = "glm"}, then \code{patch} is included as a fixed factor. All
#' other combinations of logical value and \code{approach} lead to \code{patch}
#' not being included in modeling. Defaults to \code{TRUE}.
#' @param year.as.random If set to \code{TRUE} and \code{approach = "mixed"},
#' then \code{year} is included as a random factor. If set to \code{FALSE}, then
#' \code{year} is included as a fixed factor. All other combinations of logical
#' value and \code{approach} lead to \code{year} not being included in modeling.
#' Defaults to \code{TRUE}.
#' @param juvestimate An optional variable denoting the stage name of the
#' juvenile stage in the vertical dataset. If not \code{NA}, and \code{stage} is
#' also given (see below), then vital rates listed in \code{vitalrates} other
#' than \code{fec} will also be estimated from the juvenile stage to all adult
#' stages. Defaults to \code{NA}, in which case juvenile vital rates are not
#' estimated.
#' @param juvsize A logical variable denoting whether size should be used as a
#' term in models involving transition from the juvenile stage. Defaults to
#' \code{FALSE}, and is only used if \code{juvestimate} does not equal
#' \code{NA}.
#' @param fectime A variable indicating which year of fecundity to use as the
#' response term in fecundity models. Options include \code{2}, which refers to
#' occasion \emph{t}, and \code{3}, which refers to occasion \emph{t}+1.
#' Defaults to \code{2}.
#' @param censor A vector denoting the names of censoring variables in the
#' dataset, in order from occasion \emph{t}+1, followed by occasion \emph{t},
#' and lastly followed by occasion \emph{t}-1. Defaults to \code{NA}.
#' @param age Designates the name of the variable corresponding to age in time
#' \emph{t} in the vertical dataset. Defaults to \code{NA}, in which case age
#' is not included in linear models. Should only be used if building Leslie or
#' age x stage matrices.
#' @param indcova Vector designating the names in occasions \emph{t}+1,
#' \emph{t}, and \emph{t}-1 of an individual covariate. Defaults to \code{NA}.
#' @param indcovb Vector designating the names in occasions \emph{t}+1,
#' \emph{t}, and \emph{t}-1 of a second individual covariate. Defaults to
#' \code{NA}.
#' @param indcovc Vector designating the names in occasions \emph{t}+1,
#' \emph{t}, and \emph{t}-1 of a third individual covariate. Defaults to
#' \code{NA}.
#' @param random.indcova A logical value indicating whether \code{indcova}
#' should be treated as a random categorical factor, rather than as a fixed
#' factor. Defaults to \code{FALSE}.
#' @param random.indcovb A logical value indicating whether \code{indcovb}
#' should be treated as a random categorical factor, rather than as a fixed
#' factor. Defaults to \code{FALSE}.
#' @param random.indcovc A logical value indicating whether \code{indcovc}
#' should be treated as a random categorical factor, rather than as a fixed
#' factor. Defaults to \code{FALSE}.
#' @param test.group A logical value indicating whether to include the
#' \code{group} variable from the input \code{stageframe} as a fixed categorical
#' variable in linear models. Defaults to \code{FALSE}.
#' @param ... Other parameters.
#' 
#' @return This function yields text output describing the subsets to be used in
#' linear vital rate modeling. No value or object is returned.
#' 
#' @section Notes:
#' This function is meant to handle input as would be supplied to function
#' \code{modelsearch()}. To use most easily, users may copy all input parameters
#' from a call to function \code{modelsearch()}, and paste directly within this
#' function. The exact subsets used in the \code{modelsearch()} run will also be
#' created here.
#' 
#' Tests of Gaussian normality are conducted as Shapiro-Wilk tests via base R's
#' \code{shapiro.test()} function. If datasets with more than 5000 rows are
#' supplied, function \code{hfv_qc()} will sample 5000 rows from the dataset and
#' conduct the Shapiro-Wilk test on the data sample.
#' 
#' Random factor variables are also tested for the presence of singleton
#' categories, which are factor values that occur only once in the used data
#' subset. Singleton categories may cause problems with estimation under mixed
#' modeling.
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
#'   9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
#'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
#'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#'   0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
#'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
#'   propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
#'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988", 
#'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' hfv_qc(lathvertln, historical = TRUE, suite = "main", 
#'   vitalrates = c("surv", "obs", "size", "repst", "fec"), juvestimate = "Sdl",
#'   indiv = "individ", patch = "patchid", year = "year2",year.as.random = TRUE,
#'   patch.as.random = TRUE)
#' 
#' @export
hfv_qc <- function(data, stageframe = NULL, historical = TRUE, suite = "size",
  vitalrates = c("surv", "size", "fec"), surv = c("alive3", "alive2", "alive1"),
  obs = c("obsstatus3", "obsstatus2", "obsstatus1"),
  size = c("sizea3", "sizea2", "sizea1"), sizeb = c(NA, NA, NA), 
  sizec = c(NA, NA, NA), repst = c("repstatus3", "repstatus2", "repstatus1"),
  fec = c("feca3", "feca2", "feca1"), stage = c("stage3", "stage2", "stage1"),
  matstat = c("matstatus3", "matstatus2", "matstatus1"),
  indiv = "individ", patch = NA, year = "year2", density = NA,
  patch.as.random = TRUE, year.as.random = TRUE, juvestimate = NA,
  juvsize = FALSE, fectime = 2, censor = NA, age = NA, indcova = NA,
  indcovb = NA, indcovc = NA, random.indcova = FALSE, random.indcovb = FALSE,
  random.indcovc = FALSE, test.group = FALSE, ...) {
  
  censor1 <- censor2 <- censor3 <- surv.data <- obs.data <- size.data <- NULL
  sizeb.data <- sizec.data <- juvsizeb.data <- juvsizec.data <- NULL
  repst.data <- fec.data <- juvsurv.data <- juvobs.data <- NULL
  juvsize.data <- juvrepst.data <- usedfec <- NULL
  patchcol <- yearcol <- extra_factors <- 0
  
  sizeb_used <- sizec_used <- density_used <- indcova_used <- FALSE
  indcovb_used <- indcovc_used <- indiv_used <- year_used <- patch_used <- FALSE
  
  # Some vectors for random variables
  # Names of random vars: indiv, year, patch, inda, indb, indc
  ran_vars <- c("none", "none", "none", "none", "none", "none")
  
  total_vars <- length(names(data))
  
  #Input testing, input standardization, and exception handling
  if (all(!is(data, "hfvdata"))) {
    warning("This function was made to work with standardized historically
      formatted vertical datasets, as provided by the verticalize() and
      historicalize() functions. Failure to format the input data properly and
      designate needed variables appropriately may result in nonsensical output.",
      call. = FALSE)
  }
  
  if (test.group) {
    if (is.null(stageframe)) {
      stop("Cannot test groups without inclusion of appropriate stageframe.",
        call. = FALSE)
    } else if (!is(stageframe, "stageframe")) {
      stop("Cannot test groups without inclusion of appropriate stageframe.",
        call. = FALSE)
    } else {
      all_groups <- unique(stageframe$group)
      
      if (length(all_groups) > 1) {
        extra_factors <- extra_factors + 1;
        
        data$group2 <- apply(as.matrix(data$stage2), 1, function(X){
          found_group <- which(stageframe$stage == X)
          if (length(found_group) != 1) {
            if (!is.element(tolower(X), c("notalive", "dead", "almostborn"))) {
              warning("Some group calls appear to lack a positive ID. Please check
                group identifications.", call. = FALSE)
            }
            group_call <- 0
          } else {
            group_call <- stageframe$group[found_group]
          }
          return(group_call)
        })
        data$group3 <- apply(as.matrix(data$stage3), 1, function(X){
          found_group <- which(stageframe$stage == X)
          if (length(found_group) != 1) {
            if (!is.element(tolower(X), c("notalive", "dead", "almostborn"))) {
              warning("Some group calls appear to lack a positive ID. Please check
                group identifications.", call. = FALSE)
            }
            group_call <- 0
          } else {
            group_call <- stageframe$group[found_group]
          }
          return(group_call)
        })
        data$group1 <- apply(as.matrix(data$stage1), 1, function(X){
          found_group <- which(stageframe$stage == X)
          if (length(found_group) != 1) {
            if (!is.element(tolower(X), c("notalive", "dead", "almostborn"))) {
              warning("Some group calls appear to lack a positive ID. Please check
                group identifications.", call. = FALSE)
            }
            group_call <- 0
          } else {
            group_call <- stageframe$group[found_group]
          }
          return(group_call)
        })
      } else {
        test.group <- FALSE
        warning("Only one stage group found, so will not test group.", call. = FALSE)
      }
    }
  }
  
  if (!requireNamespace("MuMIn", quietly = TRUE)) {
    stop("Package MuMIn is required. Please install it.",
      call. = FALSE)
  }
  
  # Here we will use text matching to identify variables
  vitalrates <- tolower(vitalrates)
  suite <- tolower(suite)
  
  if(length(grep("fu", suite)) > 0) {
    suite <- "full"
  } else if(length(grep("fl", suite)) > 0) {
    suite <- "full"
  } else if(length(grep("ma", suite)) > 0) {
    suite <- "main"
  } else if(length(grep("mn", suite)) > 0) {
    suite <- "main"
  } else if(length(grep("si", suite)) > 0) {
    suite <- "size"
  } else if(length(grep("sz", suite)) > 0) {
    suite <- "size"
  } else if(length(grep("re", suite)) > 0) {
    suite <- "rep"
  } else if(length(grep("rp", suite)) > 0) {
    suite <- "rep"
  } else if(length(grep("co", suite)) > 0) {
    suite <- "cons"
  } else if(length(grep("ag", suite)) > 0) {
    if (is.na(age)) {
      stop("Age variable required for age-based and age-by-stage MPMs.", call. = FALSE)
    }
    suite <- "cons"
  }
  
  models_to_estimate <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  if (length(grep("su", vitalrates)) > 0) {
    vitalrates[grep("su", vitalrates)] <- "surv"
    models_to_estimate[1] <- 1
  }
  if (length(grep("sr", vitalrates)) > 0) {
    vitalrates[grep("sr", vitalrates)] <- "surv"
    models_to_estimate[1] <- 1
  }
  if (length(grep("ob", vitalrates)) > 0) {
    vitalrates[grep("ob", vitalrates)] <- "obs"
    models_to_estimate[2] <- 1
  }
  if (length(grep("si", vitalrates)) > 0) {
    vitalrates[grep("si", vitalrates)] <- "size"
    models_to_estimate[3] <- 1
  }
  if (length(grep("sz", vitalrates)) > 0) {
    vitalrates[grep("sz", vitalrates)] <- "size"
    models_to_estimate[3] <- 1
  }
  if (length(grep("re", vitalrates)) > 0) {
    vitalrates[grep("re", vitalrates)] <- "repst"
    models_to_estimate[6] <- 1
  }
  if (length(grep("rp", vitalrates)) > 0) {
    vitalrates[grep("rp", vitalrates)] <- "repst"
    models_to_estimate[6] <- 1
  }
  if (length(grep("fe", vitalrates)) > 0) {
    vitalrates[grep("fe", vitalrates)] <- "fec"
    models_to_estimate[7] <- 1
  }
  if (length(grep("fc", vitalrates)) > 0) {
    vitalrates[grep("fc", vitalrates)] <- "fec"
    models_to_estimate[7] <- 1
  }
  if (length(vitalrates) == 1 & length(grep("lesl", vitalrates)) > 0) {
    vitalrates <- c("surv", "fec")
    models_to_estimate[1] <- 1
    models_to_estimate[7] <- 1
  }
  
  distoptions <- c("gaussian", "poisson", "negbin", "gamma")
  packoptions <- c("mixed", "glm") #The mixed option now handles all mixed models
  
  if (length(censor) > 3) {
    stop("Censor variables should be included either as 1 variable per row in the
      historical data frame (1 variable in the dataset), or as 1 variable for each of
      occasions t+1, t, and, if historical, t-1 (2 or 3 variables in the data frame).
      No more than 3 variables are allowed. If more than one are supplied, then they
      are assumed to be in order of occasion t+1, t, and t-1, respectively.",
      call. = FALSE)
  }
  if (length(indiv) > 1) {
    stop("Only one individual identity variable is allowed.", call. = FALSE)
  }
  if (length(year) > 1) {
    stop("Only one time variable is allowed. It must refer to time t.", call. = FALSE)
  }
  if (length(patch) > 1) {
    stop("Only one patch variable is allowed.", call. = FALSE)
  }
  if (length(density) > 1) {
    stop("Only one density variable is allowed.", call. = FALSE)
  }
  if (length(age) > 1) {
    stop("Only one age variable is allowed.", call. = FALSE)
  }
  
  if (is.element("surv", vitalrates)) {
    if (length(surv) > 3 | length(surv) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        survival variables as input parameters.", call. = FALSE)}
    if (all(is.numeric(surv))) {
      if (any(surv < 1) | any(surv > total_vars)) {
        stop("Survival variables do not match data frame.",
          call. = FALSE)
      } else {
        surv <- names(data)[surv]
      }
    }
    if (any(!is.element(surv, names(data)))) {
      stop("Survival variables must match data frame.", call. = FALSE)
    }
  }
  
  if (is.element("obs", vitalrates)) {
    if (length(obs) > 3 | length(obs) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        observation variables as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(obs))) {
      if (any(obs < 1) | any(obs > total_vars)) {
        stop("Observation variables do not match data frame.",
          call. = FALSE)
      } else {
       obs <- names(data)[obs]
      }
    }
    if (any(!is.element(obs, names(data)))) {
      stop("Observation status variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (is.element("size", vitalrates)) {
    if (length(size) > 3 | length(size) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical) size
        variables as input parameters.", call. = FALSE)}
    if (all(is.numeric(size))) {
      if (any(size < 1) | any(size > total_vars)) {
        stop("Size variables do not match data frame.",
          call. = FALSE)
      } else {
        size <- names(data)[size]
      }
    }
    if (any(!is.element(size, names(data)))) {
      stop("Size variables must match data frame.", call. = FALSE)
    }
    if (all(!is.na(sizeb))) {
      if (length(sizeb) > 3 | length(sizeb) == 1) {
        stop("This function requires 2 (if ahistorical) or 3 (if historical)
          secondary size variables as input parameters.",
          call. = FALSE)}
      if (all(is.numeric(sizeb))) {
        if (any(sizeb < 1) | any(sizeb > total_vars)) {
          stop("Secondary size variables do not match data frame.",
            call. = FALSE)
        } else {
          sizeb <- names(data)[sizeb]
        }
      }
      if (any(!is.element(sizeb, names(data)))) {
        stop("Secondary size variables must match data frame.",
          call. = FALSE)
      }
      sizeb_used <- 1
    }
    if (all(!is.na(sizec))) {
      if (length(sizec) > 3 | length(sizec) == 1) {
        stop("This function requires 2 (if ahistorical) or 3 (if historical)
          tertiary size variables as input parameters.",
          call. = FALSE)}
      if (all(is.numeric(sizec))) {
        if (any(sizec < 1) | any(sizec > total_vars)) {
          stop("Tertiary size variables do not match data frame.",
            call. = FALSE)
        } else {
          sizec <- names(data)[sizec]
        }
      }
      if (any(!is.element(sizec, names(data)))) {
        stop("Tertiary size variables must match data frame.",
          call. = FALSE)
      }
      sizec_used <- 1
    }
  }
  
  if (is.element("repst", vitalrates)) {
    if (length(repst) > 3 | length(repst) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        reproductive status variables as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(repst))) {
      if (any(repst < 1) | any(repst > total_vars)) {
        stop("Reproductive status variables do not match data frame.",
          call. = FALSE)
      } else {
        repst <- names(data)[repst]
      }
    }
    if (any(!is.element(repst, names(data)))) {
      stop("Reproductive status variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (is.element("fec", vitalrates)) {
    if (length(fec) > 3 | length(fec) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        fecundity variables as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(fec))) {
      if (any(fec < 1) | any(fec > total_vars)) {
        stop("Fecundity variables do not match data frame.",
          call. = FALSE)
      } else {
        fec <- names(data)[fec]
      }
    }
    if (any(!is.element(fec, names(data)))) {
      stop("Fecundity variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (fectime != 2 & fectime != 3) {
    stop("fectime must equal 2 or 3, depending on whether fecundity occurs in
      time t or t+1, respectively (the default is 2).", call. = FALSE)
  }
  
  if (all(!is.na(age))) {
    if (all(is.numeric(age))) {
      if (age > 0 & age <= total_vars) {
        agecol <- age
        age <- names(data)[agecol]
      }
    } else if (length(which(names(data) == age)) == 0) {
      stop("Variable age must either equal the name of the variable denoting age,
        or be set to NA.", call. = FALSE)
    } else {
      agecol <- which(names(data) == age)
    }
    extra_factors <- extra_factors + 1
  } else {age <- "none"}
  
  if (any(!is.na(indcova))) {
    if (length(indcova) > 3) {
      warning("Vector indcova holds the exact names of an individual covariate across
        times t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcova <- indcova[1:3]
      indcova_used <- TRUE
    } else if (length(indcova) == 1) {
      warning("Vector indcova requires the names of an individual covariate across
        times t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcova <- c("none", "none", "none")
    }
    
    if (all(is.numeric(indcova), na.rm = TRUE)) {
      if (all(indcova > 0, na.rm = TRUE) & all(indcova <= total_vars)) {
        indcova2col <- indcova[2]
        if (!is.na(indcova[3])) indcova1col <- indcova[3]
        indcova_used <- TRUE
      }
    } else {
      if (length(which(names(data) == indcova[2])) == 0 && indcova[2] != "none") {
        stop("Vector indcova must either equal either the exact names of an
          individual covariate across occasions t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcova2col <- which(names(data) == indcova[2])
        indcova_used <- TRUE
      }
      
      if (length(indcova) == 3) {
        if (length(which(names(data) == indcova[3])) == 0 && indcova[3] != "none") {
          stop("Vector indcova must either equal either the exact names of an
            individual covariate across occasions t+1, t, and t-1, or be set to NA.",
            call. = FALSE)
        } else {
          indcova1col <- which(names(data) == indcova[3])
          indcova_used <- TRUE
        }
      } else {
        indcova[3] <- "none"
      }
    }
    ran_vars[4] = names(data)[indcova2col]
  } else {indcova <- c("none", "none", "none")}
  
  if (indcova_used & !random.indcova) {
    extra_factors <- extra_factors + 1
  }
  
  if (any(!is.na(indcovb))) {
    if (length(indcovb) > 3) {
      warning("Vector indcovb holds the exact names of an individual covariate
        across times t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcovb <- indcovb[1:3]
      indcovb_used <- TRUE
    } else if (length(indcovb) == 1) {
      warning("Vector indcovb requires the names of an individual covariate across
        times t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcovb <- c("none", "none", "none")
    }
    
    if (all(is.numeric(indcovb), na.rm = TRUE)) {
      if (all(indcovb > 0, na.rm = TRUE) & all(indcovb <= total_vars)) {
        indcovb2col <- indcovb[2]
        if (!is.na(indcovb[3])) indcovb1col <- indcovb[3]
        indcovb_used <- TRUE
      }
    } else {
      if (length(which(names(data) == indcovb[2])) == 0 && indcovb[2] != "none") {
        stop("Vector indcovb must either equal either the exact names of an
          individual covariate across times t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcovb2col <- which(names(data) == indcovb[2])
        indcovb_used <- TRUE
      }
      
      if (length(indcovb) == 3) {
        if (length(which(names(data) == indcovb[3])) == 0 && indcovb[3] != "none") {
          stop("Vector indcovb must either equal either the exact names of an
            individual covariate across times t+1, t, and t-1, or be set to NA.",
            call. = FALSE)
        } else {
          indcovb1col <- which(names(data) == indcovb[3])
          indcovb_used <- TRUE
        }
      } else {
        indcovb[3] <- "none"
      }
    }
    ran_vars[5] = names(data)[indcovb2col]
  } else {indcovb <- c("none", "none", "none")}
  
  if (indcovb_used & !random.indcovb) {
    extra_factors <- extra_factors + 1
  }
  
  if (any(!is.na(indcovc))) {
    if (length(indcovc) > 3) {
      warning("Vector indcovc holds the exact names of an individual covariate
        across times t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcovc <- indcovc[1:3]
      indcovc_used <- TRUE
    } else if (length(indcovc) == 1) {
      warning("Vector indcovc requires the names of an individual covariate across
        times t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcovc <- c("none", "none", "none")
    } 
    
    if (all(is.numeric(indcovc), na.rm = TRUE)) {
      if (all(indcovc > 0, na.rm = TRUE) & all(indcovc <= total_vars)) {
        indcovc2col <- indcovc[2]
        if (!is.na(indcovc[3])) indcovc1col <- indcovc[3]
        indcovc_used <- TRUE
      }
    } else {
      if (length(which(names(data) == indcovc[2])) == 0 && indcovc[2] != "none") {
        stop("Vector indcovc must either equal either the exact names of an
          individual covariate across times t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcovc2col <- which(names(data) == indcovc[2])
        indcovc_used <- TRUE
      }
      
      if (length(indcovc) == 3) {
        if (length(which(names(data) == indcovc[3])) == 0 && indcovc[3] != "none") {
          stop("Vector indcovc must either equal either the exact names of an
            individual covariate across times t+1, t, and t-1, or be set to NA.",
            call. = FALSE)
        } else {
          indcovc1col <- which(names(data) == indcovc[3])
          indcovc_used <- TRUE
        }
      } else {
        indcovc[3] <- "none"
      }
    }
    ran_vars[6] = names(data)[indcovc2col]
  } else {indcovc <- c("none", "none", "none")}
  
  if (indcovc_used & !random.indcovc) {
    extra_factors <- extra_factors + 1
  }
  
  if (!is.na(indiv)) {
    if (!is.numeric(indiv) & length(which(names(data) == indiv)) == 0) {
      stop("Variable indiv must either equal the exact name of the variable denoting
        individual identity in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(indiv)){
      indivcol <- which(names(data) == indiv)
      
      if (any(is.na(data[,indivcol])) ) {
        warning("NAs in individual ID variable may cause unexpected behavior in mixed
          model building. Please rename all individuals with unique names, avoiding NAs.",
          call. = FALSE)
      }
      indiv_used <- TRUE
      
    } else if (is.numeric(indiv)) {
      if (any(indiv < 1) | any(indiv > total_vars)) {
        stop("Unable to interpret indiv variable.", call. = FALSE)
      } else {
        indivcol <- indiv
        indiv <- names(data)[indivcol]
        indiv_used <- TRUE
      }
    } else {
      stop("Unable to interpret indiv variable.", call. = FALSE)
    }
    ran_vars[1] = names(data)[indivcol]
  } else {indiv <- "none"}
  
  if (!is.na(patch)) {
    if (!is.numeric(patch) & length(which(names(data) == patch)) == 0) {
      stop("Variable patch must either equal the exact name of the variable denoting
        patch identity in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(patch)) {
      patchcol <- which(names(data) == patch)
    } else if (is.numeric(patch)) {
      if (any(patch < 1) | any(patch > total_vars)) {
        stop("Unable to interpret patch variable.", call. = FALSE)
      } else {
        patchcol <- patch  # Used to be names(data)[patch]
        patch <- names(data)[patchcol] 
      }
    } else {
      stop("Unable to interpret patch variable.", call. = FALSE)
    }
    ran_vars[3] = names(data)[patchcol]
    patch_used <- TRUE
  } else {patch <- "none"}
  
  if (patchcol > 0 & !patch.as.random) {
    extra_factors <- extra_factors + 1
  }
  
  if (!is.na(year)) {
    if (!is.numeric(year) & length(which(names(data) == year)) == 0) {
      stop("Variable year must either equal the exact name of the variable denoting
        occasion t in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(year)) {
      yearcol <- which(names(data) == year)
    } else if (is.numeric(year)) {
      if (any(year < 1) | any(year > total_vars)) {
        stop("Unable to interpret year variable.", call. = FALSE)
      } else {
        yearcol <- year  # Used to be names(data)[year]
        year <- names(data)[yearcol] 
      }
    } else {
      stop("Unable to interpret year variable.", call. = FALSE)
    }
    ran_vars[2] = names(data)[yearcol]
    year_used <- TRUE
  } else {year <- "none"}
  
  if (yearcol > 0 & !year.as.random) {
    extra_factors <- extra_factors + 1
  }
  
  if (!is.na(density)) {
    if (!is.numeric(density) & length(which(names(data) == density)) == 0) {
      stop("Variable density must either equal the exact name of the variable
        denoting spatial density in occasion t in the dataset, or be set to NA.",
        call. = FALSE)
    } else if (is.character(density)) {
      if (is.element(density, names(data))) {
        density_used <- TRUE
      } else {
        stop("Unable to interpret density variable.", call. = FALSE)
      }
    } else if (is.numeric(density)) {
      if (any(density < 1) | any(density > total_vars)) {
        stop("Unable to interpret density variable.", call. = FALSE)
      } else {
        density <- names(data)[density]
        density_used <- TRUE
      }
    } else {
      stop("Unable to interpret density variable.", call. = FALSE)
    }
  } else {density <- "none"}
  
  if (density_used) {
    extra_factors <- extra_factors + 1
  }
  
  # Here we test the dataset for appropriate stage names
  if (!is.na(juvestimate)) {
    if (!any(is.element(stage, names(data)))) {
      stop("Names of stage variables do not match the dataset.", call. = FALSE)}
    
    stage3col <- which(names(data) == stage[1])
    stage2col <- which(names(data) == stage[2])
    if (length(stage) == 3) {stage1col <- which(names(data) == stage[3])} else {stage1col <- 0}
    
    if (!is.element(juvestimate, unique(data[,stage2col]))) {
      stop("The stage declared as juvenile via juvestimate is not recognized within
        the stage2 variable in the dataset.", call. = FALSE)
    }
    if (length(matstat) > 3 | length(matstat) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical) maturity status
        variables as input parameters if juvenile parameters are to be estimated.",
        call. = FALSE)}
    if (all(is.numeric(matstat))) {
      if (any(matstat < 1) | any(matstat > total_vars)) {
        stop("Maturity status variables do not match data frame.",
          call. = FALSE)
      } else {
        matstat <- names(data)[matstat]
      }
    }
    if (any(!is.element(matstat, names(data)))) {
      stop("Maturity status variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (!is.logical(juvsize)) {
    stop("Option juvsize must be set to either TRUE or FALSE.", call. = FALSE)
  }
  
  {
    if (random.indcova | random.indcovb | random.indcovc) {
      warning("Random covariates can only be included in mixed models. Setting
        random.indcova, random.indcovb, and random.indcovc to FALSE.",
        call. = FALSE)
      random.indcova <- FALSE
      random.indcovb <- FALSE
      random.indcovc <- FALSE
    }
    
  }
  
  #Now we need to create the input datasets
  if (!all(is.na(censor))) {
    if (length(censor) == 1) {
      data$censor2 <- data[, which(names(data) == censor[1])]
    } else {
      data$censor3 <- data[, which(names(data) == censor[1])]
      data$censor2 <- data[, which(names(data) == censor[2])]
      if (length(censor) > 2) {
        data$censor1 <- data[, which(names(data) == censor[3])]
      }
    }
  } else {
    data$censor2 <- 1
  }
  data <- subset(data, censor2 == 1)
  
  if (!all(is.na(censor))) {
    if (length(censor) > 1) {
      data <- subset(data, censor3 == 1)
      if (length(censor) > 2) {
        data <- subset(data, censor1 == 1)
      }
    }
  }
  
  juvsurv.sole <- juvobs.sole <- juvsize.sole <- c(0, 0, 0, 0, 0, 0)
  juvsizeb.sole <- juvsizec.sole <- juvrepst.sole <- c(0, 0, 0, 0, 0, 0)
  surv.sole <- obs.sole <- size.sole <- sizeb.sole <- c(0, 0, 0, 0, 0, 0)
  sizec.sole <- repst.sole <- fec.sole <- c(0, 0, 0, 0, 0, 0)
  
  if (!is.na(juvestimate)) {
    juvindivs <- which(data[,stage2col] == juvestimate)
    adultindivs <- setdiff(c(1:length(data[,stage2col])), juvindivs)
    
    if (models_to_estimate[1] == 1) {models_to_estimate[8] = 1}
    juvsurv.data <- subset(data, data[,stage2col] == juvestimate & data[,which(names(data) == surv[2])] == 1)
    juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
    juvsurv.trans <- dim(juvsurv.data)[1]
    
    if (indiv_used) juvsurv.sole[1] <- length(which(table(juvsurv.data[, which(names(juvsurv.data) == indiv)]) == 1))
    if (year_used) juvsurv.sole[2] <- length(which(table(juvsurv.data[, which(names(juvsurv.data) == year)]) == 1))
    if (patch_used) juvsurv.sole[3] <- length(which(table(juvsurv.data[, which(names(juvsurv.data) == patch)]) == 1))
    if (indcova_used & random.indcova) juvsurv.sole[4] <- length(which(table(juvsurv.data[, which(names(juvsurv.data) == indcova[2])]) == 1))
    if (indcovb_used & random.indcovb) juvsurv.sole[5] <- length(which(table(juvsurv.data[, which(names(juvsurv.data) == indcovb[2])]) == 1))
    if (indcovc_used & random.indcovc) juvsurv.sole[6] <- length(which(table(juvsurv.data[, which(names(juvsurv.data) == indcovc[2])]) == 1))
    
    if (suite == "full" | suite == "main" | suite == "size") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[2])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[3])]))) {
          warning("NAs in size variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
    }
    
    if (suite == "full" | suite == "main" | suite == "rep") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[2])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.",
          call. = FALSE)
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[3])]))) {
          warning("NAs in reproductive status variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
    }
    
    if (is.element(0, juvsurv.data$matstatus3)) {
      warning("Function modelsearch() assumes that all juveniles either die or transition
        to maturity within 1 year. Some individuals in this dataset appear to live
        longer as juveniles than assumptions allow.", call. = FALSE)
    }
    
    juvobs.data <- subset(juvsurv.data, juvsurv.data[, which(names(juvsurv.data) == surv[1])] == 1)
    juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
    juvobs.trans <- dim(juvobs.data)[1]
    
    if (indiv_used) juvobs.sole[1] <- length(which(table(juvobs.data[, which(names(juvobs.data) == indiv)]) == 1))
    if (year_used) juvobs.sole[2] <- length(which(table(juvobs.data[, which(names(juvobs.data) == year)]) == 1))
    if (patch_used) juvobs.sole[3] <- length(which(table(juvobs.data[, which(names(juvobs.data) == patch)]) == 1))
    if (indcova_used & random.indcova) juvobs.sole[4] <- length(which(table(juvobs.data[, which(names(juvobs.data) == indcova[2])]) == 1))
    if (indcovb_used & random.indcovb) juvobs.sole[5] <- length(which(table(juvobs.data[, which(names(juvobs.data) == indcovb[2])]) == 1))
    if (indcovc_used & random.indcovc) juvobs.sole[6] <- length(which(table(juvobs.data[, which(names(juvobs.data) == indcovc[2])]) == 1))
    
    if (models_to_estimate[2] == 1) {
      models_to_estimate[9] == 1
      juvsize.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
      
      if (sizeb_used) {
        juvsizeb.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
        juvsizeb.data <- juvsizeb.data[which(!is.na(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])])),]
      }
      
      if (sizec_used) {
        juvsizec.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
        juvsizec.data <- juvsizec.data[which(!is.na(juvsizec.data[, which(names(juvsizec.data) == sizec[1])])),]
      }
    } else {
      juvsize.data <- juvobs.data
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
      
      if (sizeb_used) {
        juvsizeb.data <- juvobs.data
        juvsizeb.data <- juvsizeb.data[which(!is.na(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])])),]
      }
      
      if (sizec_used) {
        juvsizec.data <- juvobs.data
        juvsizec.data <- juvsizec.data[which(!is.na(juvsizec.data[, which(names(juvsizec.data) == sizec[1])])),]
      }
    }
    juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
    juvsize.trans <- dim(juvsize.data)[1]
    
    if (indiv_used) juvsize.sole[1] <- length(which(table(juvsize.data[, which(names(juvsize.data) == indiv)]) == 1))
    if (year_used) juvsize.sole[2] <- length(which(table(juvsize.data[, which(names(juvsize.data) == year)]) == 1))
    if (patch_used) juvsize.sole[3] <- length(which(table(juvsize.data[, which(names(juvsize.data) == patch)]) == 1))
    if (indcova_used & random.indcova) juvsize.sole[4] <- length(which(table(juvsize.data[, which(names(juvsize.data) == indcova[2])]) == 1))
    if (indcovb_used & random.indcovb) juvsize.sole[5] <- length(which(table(juvsize.data[, which(names(juvsize.data) == indcovb[2])]) == 1))
    if (indcovc_used & random.indcovc) juvsize.sole[6] <- length(which(table(juvsize.data[, which(names(juvsize.data) == indcovc[2])]) == 1))
    
    if (sizeb_used) {
      juvsizeb.ind <- length(unique(juvsizeb.data[, which(names(juvsizeb.data) == indiv)]))
      juvsizeb.trans <- dim(juvsizeb.data)[1]
      
      if (indiv_used) juvsizeb.sole[1] <- length(which(table(juvsizeb.data[, which(names(juvsizeb.data) == indiv)]) == 1))
      if (year_used) juvsizeb.sole[2] <- length(which(table(juvsizeb.data[,   which(names(juvsizeb.data) == year)]) == 1))
      if (patch_used) juvsizeb.sole[3] <- length(which(table(juvsizeb.data[,    which(names(juvsizeb.data) == patch)]) == 1))
      if (indcova_used & random.indcova) juvsizeb.sole[4] <- length(which(table(juvsizeb.data[, which(names(juvsizeb.data) == indcova[2])]) == 1))
      if (indcovb_used & random.indcovb) juvsizeb.sole[5] <- length(which(table(juvsizeb.data[, which(names(juvsizeb.data) == indcovb[2])]) == 1))
      if (indcovc_used & random.indcovc) juvsizeb.sole[6] <- length(which(table(juvsizeb.data[, which(names(juvsizeb.data) == indcovc[2])]) == 1))
      
    }
    if (sizec_used) {
      juvsizec.ind <- length(unique(juvsizec.data[, which(names(juvsizec.data) == indiv)]))
      juvsizec.trans <- dim(juvsizec.data)[1]
      
      if (indiv_used) juvsizec.sole[1] <- length(which(table(juvsizec.data[, which(names(juvsizec.data) == indiv)]) == 1))
      if (year_used) juvsizec.sole[2] <- length(which(table(juvsizec.data[, which(names(juvsizec.data) == year)]) == 1))
      if (patch_used) juvsizec.sole[3] <- length(which(table(juvsizec.data[, which(names(juvsizec.data) == patch)]) == 1))
      if (indcova_used & random.indcova) juvsizec.sole[4] <- length(which(table(juvsizec.data[, which(names(juvsizec.data) == indcova[2])]) == 1))
      if (indcovb_used & random.indcovb) juvsizec.sole[5] <- length(which(table(juvsizec.data[, which(names(juvsizec.data) == indcovb[2])]) == 1))
      if (indcovc_used & random.indcovc) juvsizec.sole[6] <- length(which(table(juvsizec.data[, which(names(juvsizec.data) == indcovc[2])]) == 1))
      
    }
    
    juvrepst.data <- juvsize.data
    juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
    juvrepst.trans <- dim(juvrepst.data)[1]
    
    if (indiv_used) juvrepst.sole <- length(which(table(juvrepst.data[, which(names(juvrepst.data) == indiv)]) == 1))
    
    data <- data[adultindivs,] #This line resets the main dataset to adults only
  }
  
  surv.data <- subset(data, data[,which(names(data) == surv[2])] == 1)
  surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
  surv.trans <- dim(surv.data)[1]
  
  if (indiv_used) surv.sole[1] <- length(which(table(surv.data[, which(names(surv.data) == indiv)]) == 1))
  if (year_used) surv.sole[2] <- length(which(table(surv.data[, which(names(surv.data) == year)]) == 1))
  if (patch_used) surv.sole[3] <- length(which(table(surv.data[, which(names(surv.data) == patch)]) == 1))
  if (indcova_used & random.indcova) surv.sole[4] <- length(which(table(surv.data[, which(names(surv.data) == indcova[2])]) == 1))
  if (indcovb_used & random.indcovb) surv.sole[5] <- length(which(table(surv.data[, which(names(surv.data) == indcovb[2])]) == 1))
  if (indcovc_used & random.indcovc) surv.sole[6] <- length(which(table(surv.data[, which(names(surv.data) == indcovc[2])]) == 1))
  
  if(any(!suppressWarnings(!is.na(as.numeric(as.character(surv.data[, which(names(surv.data) == size[1])])))))) {
    warning("Modelsearch(), flefko3(), flefko2(), and aflefko2() are made to work
      with numeric size variables. Use of categorical variables may result in
      errors and unexpected behavior.", call. = FALSE)
  }
  if (suite == "full" | suite == "main" | suite == "size") {
    if (any(is.na(surv.data[, which(names(surv.data) == size[2])]))) {
      warning("NAs in size variables may cause model selection to fail.",
        call. = FALSE)
    }
    if (sizeb_used) {
      if (any(is.na(surv.data[, which(names(surv.data) == sizeb[2])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
    }
    if (sizec_used) {
      if (any(is.na(surv.data[, which(names(surv.data) == sizec[2])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == size[3])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
      if (sizeb_used) {
        if (any(is.na(surv.data[, which(names(surv.data) == sizeb[3])]))) {
          warning("NAs in size variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
      if (sizec_used) {
        if (any(is.na(surv.data[, which(names(surv.data) == sizec[3])]))) {
          warning("NAs in size variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
    }
  }
  if (suite == "full" | suite == "main" | suite == "rep") {
    if (any(is.na(surv.data[, which(names(surv.data) == repst[2])]))) {
      warning("NAs in reproductive status variables may cause model selection to fail.",
        call. = FALSE)
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == repst[3])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.",
          call. = FALSE)
      }
    }
  }
  
  obs.data <- subset(surv.data, surv.data[, which(names(surv.data) == surv[1])] == 1)
  obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
  obs.trans <- dim(obs.data)[1]
  
  if (indiv_used) obs.sole[1] <- length(which(table(obs.data[, which(names(obs.data) == indiv)]) == 1))
  if (year_used) obs.sole[2] <- length(which(table(obs.data[, which(names(obs.data) == year)]) == 1))
  if (patch_used) obs.sole[3] <- length(which(table(obs.data[, which(names(obs.data) == patch)]) == 1))
  if (indcova_used & random.indcova) obs.sole[4] <- length(which(table(obs.data[, which(names(obs.data) == indcova[2])]) == 1))
  if (indcovb_used & random.indcovb) obs.sole[5] <- length(which(table(obs.data[, which(names(obs.data) == indcovb[2])]) == 1))
  if (indcovc_used & random.indcovc) obs.sole[6] <- length(which(table(obs.data[, which(names(obs.data) == indcovc[2])]) == 1))
  
  if (models_to_estimate[2] == 1) {
    size.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
    
    if (sizeb_used) {
      sizeb.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
      sizeb.data <- sizeb.data[which(!is.na(sizeb.data[, which(names(sizeb.data) == sizeb[1])])),]
    }
    if (sizec_used) {
      sizec.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
      sizec.data <- sizec.data[which(!is.na(sizec.data[, which(names(sizec.data) == sizec[1])])),]
    }
  } else {
    size.data <- obs.data
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
    
    if (sizeb_used) {
      sizeb.data <- obs.data
      sizeb.data <- sizeb.data[which(!is.na(sizeb.data[, which(names(sizeb.data) == sizeb[1])])),]
    }
    if (sizec_used) {
      sizec.data <- obs.data
      sizec.data <- sizec.data[which(!is.na(sizec.data[, which(names(sizec.data) == sizec[1])])),]
    }
  }
  size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
  size.trans <- dim(size.data)[1]
  
  if (indiv_used) size.sole[1] <- length(which(table(size.data[, which(names(size.data) == indiv)]) == 1))
  if (year_used) size.sole[2] <- length(which(table(size.data[, which(names(size.data) == year)]) == 1))
  if (patch_used) size.sole[3] <- length(which(table(size.data[, which(names(size.data) == patch)]) == 1))
  if (indcova_used & random.indcova) size.sole[4] <- length(which(table(size.data[, which(names(size.data) == indcova[2])]) == 1))
  if (indcovb_used & random.indcovb) size.sole[5] <- length(which(table(size.data[, which(names(size.data) == indcovb[2])]) == 1))
  if (indcovc_used & random.indcovc) size.sole[6] <- length(which(table(size.data[, which(names(size.data) == indcovc[2])]) == 1))
  
  if (sizeb_used) {
    sizeb.ind <- length(unique(sizeb.data[, which(names(sizeb.data) == indiv)]))
    sizeb.trans <- dim(sizeb.data)[1]
    
    if (indiv_used) sizeb.sole[1] <- length(which(table(sizeb.data[, which(names(sizeb.data) == indiv)]) == 1))
    if (year_used) sizeb.sole[2] <- length(which(table(sizeb.data[, which(names(sizeb.data) == year)]) == 1))
    if (patch_used) sizeb.sole[3] <- length(which(table(sizeb.data[, which(names(sizeb.data) == patch)]) == 1))
    if (indcova_used & random.indcova) sizeb.sole[4] <- length(which(table(sizeb.data[, which(names(sizeb.data) == indcova[2])]) == 1))
    if (indcovb_used & random.indcovb) sizeb.sole[5] <- length(which(table(sizeb.data[, which(names(sizeb.data) == indcovb[2])]) == 1))
    if (indcovc_used & random.indcovc) sizeb.sole[6] <- length(which(table(sizeb.data[, which(names(sizeb.data) == indcovc[2])]) == 1))
    
  }
  if (sizec_used) {
    sizec.ind <- length(unique(sizec.data[, which(names(sizec.data) == indiv)]))
    sizec.trans <- dim(sizec.data)[1]
    
    if (indiv_used) sizec.sole[1] <- length(which(table(sizec.data[, which(names(sizec.data) == indiv)]) == 1))
    if (year_used) sizec.sole[2] <- length(which(table(sizec.data[, which(names(sizec.data) == year)]) == 1))
    if (patch_used) sizec.sole[3] <- length(which(table(sizec.data[, which(names(sizec.data) == patch)]) == 1))
    if (indcova_used & random.indcova) sizec.sole[4] <- length(which(table(sizec.data[, which(names(sizec.data) == indcova[2])]) == 1))
    if (indcovb_used & random.indcovb) sizec.sole[5] <- length(which(table(sizec.data[, which(names(sizec.data) == indcovb[2])]) == 1))
    if (indcovc_used & random.indcovc) sizec.sole[6] <- length(which(table(sizec.data[, which(names(sizec.data) == indcovc[2])]) == 1))
    
  }
  
  repst.data <- size.data
  
  if (indiv_used) repst.sole[1] <- length(which(table(repst.data[, which(names(repst.data) == indiv)]) == 1))
  if (year_used) repst.sole[2] <- length(which(table(repst.data[, which(names(repst.data) == year)]) == 1))
  if (patch_used) repst.sole[3] <- length(which(table(repst.data[, which(names(repst.data) == patch)]) == 1))
  if (indcova_used & random.indcova) repst.sole[4] <- length(which(table(repst.data[, which(names(repst.data) == indcova[2])]) == 1))
  if (indcovb_used & random.indcovb) repst.sole[5] <- length(which(table(repst.data[, which(names(repst.data) == indcovb[2])]) == 1))
  if (indcovc_used & random.indcovc) repst.sole[6] <- length(which(table(repst.data[, which(names(repst.data) == indcovc[2])]) == 1))
  
  if (models_to_estimate[6] == 1) {
    if (fectime == 2) {
      fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[2])] == 1)
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[1])] == 1)
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
    }
  } else {
    fec.data <- surv.data
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
    }
  }
  fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
  fec.trans <- dim(fec.data)[1]
  
  if (indiv_used) fec.sole[1] <- length(which(table(fec.data[, which(names(fec.data) == indiv)]) == 1))
  if (year_used) fec.sole[2] <- length(which(table(fec.data[, which(names(fec.data) == year)]) == 1))
  if (patch_used) fec.sole[3] <- length(which(table(fec.data[, which(names(fec.data) == patch)]) == 1))
  if (indcova_used & random.indcova) fec.sole[4] <- length(which(table(fec.data[, which(names(fec.data) == indcova[2])]) == 1))
  if (indcovb_used & random.indcovb) fec.sole[5] <- length(which(table(fec.data[, which(names(fec.data) == indcovb[2])]) == 1))
  if (indcovc_used & random.indcovc) fec.sole[6] <- length(which(table(fec.data[, which(names(fec.data) == indcovc[2])]) == 1))
  
  if (is.element("fec", vitalrates)) {
    if (fectime == 2) {
      usedfec <- which(names(fec.data) == fec[2])
    } else if (fectime == 3) {
      usedfec <- which(names(fec.data) == fec[1])
    }
  }
  
  ran_vars_names <- c("indiv id: ", "year2: ", "patch: ", "indcova: ", "indcovb: ", "indcovc: ")
  chosen_ran_vars <- which(ran_vars != "none")
  
  # The major variable tests
  if (is.data.frame(surv.data) & is.element("surv", vitalrates)) {
    writeLines("Survival:\n")
    writeLines(paste0("  Data subset has ", dim(surv.data)[2], " variables and ",
      dim(surv.data)[1], " transitions.\n"))
    
    .intbin_check(surv.data, surv[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(surv.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", surv.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(obs.data) & is.element("obs", vitalrates)) {
    writeLines("\nObservation status:\n")
    writeLines(paste0("  Data subset has ", dim(obs.data)[2], " variables and ",
      dim(obs.data)[1], " transitions.\n"))
    
    .intbin_check(obs.data, obs[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(obs.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", obs.sole[i], ")"))
      }
    }
  }
  
  if (is.data.frame(size.data) & is.element("size", vitalrates)) {
    writeLines("\nPrimary size:\n")
    writeLines(paste0("  Data subset has ", dim(size.data)[2], " variables and ",
      dim(size.data)[1], " transitions.\n"))
    
    .sf_dist_check(size.data, size[1], size[2])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(size.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", size.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(sizeb.data) & is.element("siz", vitalrates)) {
    writeLines("\nSecondary size:\n")
    writeLines(paste0("  Data subset has ", dim(sizeb.data)[2], " variables and ",
      dim(sizeb.data)[1], " transitions.\n"))
    
    .sf_dist_check(sizeb.data, sizeb[1], sizeb[2])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(sizeb.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", sizeb.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(sizec.data) & is.element("siz", vitalrates)) {
    writeLines("\nTertiary size:\n")
    writeLines(paste0("  Data subset has ", dim(sizec.data)[2], " variables and ",
      dim(sizec.data)[1], " transitions.\n"))
    
    .sf_dist_check(sizec.data, sizec[1], sizec[2])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(sizec.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", sizec.sole[i], ")"))
      }
    }
  }
  
  if (is.data.frame(repst.data) & is.element("repst", vitalrates)) {
    writeLines("\nReproductive status:\n")
    writeLines(paste0("  Data subset has ", dim(repst.data)[2], " variables and ",
      dim(repst.data)[1], " transitions.\n"))
    
    .intbin_check(repst.data, repst[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(repst.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", repst.sole[i], ")"))
      }
    }
  }
  
  if (is.data.frame(fec.data) & is.element("fec", vitalrates)) {
    writeLines("\nFecundity:\n")
    writeLines(paste0("  Data subset has ", dim(fec.data)[2], " variables and ",
      dim(fec.data)[1], " transitions.\n"))
    
    if (fectime == 2) {
      term2 <- size[2]
    } else {
      term2 <- size[1]
    }
    .sf_dist_check(fec.data, names(fec.data)[usedfec], term2)
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(fec.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", fec.sole[i], ")"))
      }
    }
  }
  
  if (is.data.frame(juvsurv.data) & is.element("surv", vitalrates)) {
    writeLines("\nJuvenile survival:\n")
    writeLines(paste0("  Data subset has ", dim(juvsurv.data)[2], " variables and ",
      dim(juvsurv.data)[1], " transitions.\n"))
    
    .intbin_check(juvsurv.data, surv[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvsurv.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvsurv.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(juvobs.data) & is.element("obs", vitalrates)) {
    writeLines("\nJuvenile observation status:\n")
    writeLines(paste0("  Data subset has ", dim(juvobs.data)[2], " variables and ",
      dim(juvobs.data)[1], " transitions.\n"))
    
    .intbin_check(juvobs.data, obs[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvobs.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvobs.sole[i], ")"))
      }
    }
  }
  
  if (is.data.frame(juvsize.data) & is.element("size", vitalrates)) {
    writeLines("\nJuvenile primary size:\n")
    writeLines(paste0("  Data subset has ", dim(juvsize.data)[2], " variables and ",
      dim(juvsize.data)[1], " transitions.\n"))
    
    .sf_dist_check(juvsize.data, size[1], size[2])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvsize.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvsize.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(juvsizeb.data)& is.element("siz", vitalrates)) {
    writeLines("\nJuvenile secondary size:\n")
    writeLines(paste0("  Data subset has ", dim(juvsizeb.data)[2], " variables and ",
      dim(juvsizeb.data)[1], " transitions.\n"))
    
    .sf_dist_check(juvsizeb.data, sizeb[1], sizeb[2])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvsizeb.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvsizeb.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(juvsizec.data) & is.element("siz", vitalrates)) {
    writeLines("\nJuvenile tertiary size:\n")
    writeLines(paste0("  Data subset has ", dim(juvsizec.data)[2], " variables and ",
      dim(juvsizec.data)[1], " transitions.\n"))
    
    .sf_dist_check(juvsizec.data, sizec[1], sizec[2])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvsizec.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvsizec.sole[i], ")"))
      }
    }
  }
  
  if (is.data.frame(juvrepst.data) & is.element("repst", vitalrates)) {
    writeLines("\nJuvenile reproductive status:\n")
    writeLines(paste0("  Data subset has ", dim(juvrepst.data)[2], " variables and ",
      dim(juvrepst.data)[1], " transitions.\n"))
    
    .intbin_check(juvrepst.data, repst[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvrepst.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvrepst.sole[i], ")"))
      }
    }
  }
  if (is.data.frame(juvobs.data)) {
    writeLines("\nJuvenile maturity status:\n")
    writeLines(paste0("  Data subset has ", dim(juvobs.data)[2], " variables and ",
      dim(juvobs.data)[1], " transitions.\n"))
    
    .intbin_check(juvobs.data, matstat[1])
    
    if(length(chosen_ran_vars) > 0) {
      writeLines("  Numbers of categories in data subset in possible random variables:")
      
      for(i in c(1:length(chosen_ran_vars))) {
        found_uns <- length(unique(juvobs.data[,ran_vars[chosen_ran_vars[i]]]))
        writeLines(paste0("  ", ran_vars_names[chosen_ran_vars[i]], found_uns,
            "   (singleton categories: ", juvobs.sole[i], ")"))
      }
    }
  }
}

#' Create a Vector of Stage Weights for Density Dependence Calculations
#' 
#' Function \code{stage_weight()} creates a data frame summarizing the weighting
#' of individuals of different stages, stage-stage combinations, ages, or
#' age-stages relative to one another in terms of the total population size.
#' Stage weights can be thought of as Lotka-Volterra coefficients applied to
#' individuals of different stages within a single population.
#' 
#' @name stage_weight
#' 
#' @param mpm The lefkoMat object to be used in projection. Can be an example
#' MPM if function-based projection is planned.
#' @param stage2 A vector showing the name or number of a stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Abbreviations for groups of stages are also usable (see Notes).
#' This input is required for all stage-based and age-by-stage MPMs. Defaults to
#' \code{NA}.
#' @param stage1 A vector showing the name or number of a stage in occasion
#' \emph{t}-1 that should be set to a positive number of individuals in the
#' start vector. Abbreviations for groups of stages are also usable (see Notes).
#' This is only used for historical MPMs, since the rows of hMPMs correspond to
#' stage-pairs in times \emph{t} and \emph{t}-1 together. Only required for
#' historical MPMs, and will result in errors if otherwise used.
#' @param age2 A vector showing the age of each respective stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Only used for Leslie and age-by-stage MPMs. Defaults to \code{NA}.
#' @param value A vector showing the values, in order, of the number of
#' individuals set for the stage or stage-pair in question. Defaults to
#' \code{1}.
#' 
#' @return A list of class \code{lefkoEq}, with four objects, which can be
#' used as input in function \code{\link{projection3}()} and
#' \code{\link{f_projection3}()}. The last three include the \code{ahstages},
#' \code{hstages}, and \code{agestages} objects from the \code{lefkoMat} object
#' supplied in \code{mpm}. The first element in the list is a data frame with
#' the following variables:
#' 
#' \item{stage2}{Stage at occasion \emph{t}.}
#' \item{stage_id_2}{The stage number associated with \code{stage2}.}
#' \item{stage1}{Stage at occasion \emph{t}-1, if historical. Otherwise NA.}
#' \item{stage_id_1}{The stage number associated with \code{stage1}.}
#' \item{age2}{The age of individuals in \code{stage2} and, if applicable,
#' \code{stage1}. Only used in age-by-stage MPMs.}
#' \item{row_num}{A number indicating the respective starting vector element.}
#' \item{value}{Number of individuals in corresponding stage or stage-pair.}
#' 
#' @section Notes:
#' Users should generally consider which stage to set as the reference, which
#' should be designated with a value of \code{1}.
#' 
#' In some life histories, certain stages may not count against population size
#' with regards to the kind and value of density dependence decided on. An
#' example might be the dormant seed stage in some plant species. These stage
#' weight values should be set to \code{0}.
#' 
#' Entries in \code{stage2}, and \code{stage1} can include abbreviations for
#' groups of stages. Use \code{rep} if all reproductive stages are to be used,
#' \code{nrep} if all mature but non-reproductive stages are to be used,
#' \code{mat} if all mature stages are to be used, \code{immat} if all immature
#' stages are to be used, \code{prop} if all propagule stages are to be used,
#' \code{npr} if all non-propagule stages are to be used, \code{obs} if all
#' observable stages are to be used, \code{nobs} if all unobservable stages are
#' to be used, and leave empty or use \code{all} if all stages in stageframe are
#' to be used.
#' 
#' 
#' @examples
#' library(lefko3)
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
#' cycaraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
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
#' cycamatrix2r <- rlefko2(data = cycaraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cyca2_start <- start_input(cycamatrix2r, stage2 = c("SD", "P1", "P2"),
#'   value = c(500, 100, 200))
#' 
#' cyp2_dv <- density_input(cycamatrix2r, stage3 = c("SD", "P1"),
#'   stage2 = c("rep", "rep"), style = c(1, 1), alpha = c(0.5, 1.2),
#'   beta = c(1.0, 2.0), type = c(2, 1))
#' 
#' cyp_eq <- stage_weight(cycamatrix2r,
#'   stage2 = c("SD", "P1", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg"),
#'   value = c(0, 1, 1, 1, 1, 1, 1, 1, 1))
#' 
#' cyca_proj <- projection3(mpm = cycamatrix2r, start_frame = cyca2_start,
#'   density = cyp2_dv, stage_weights = cyp_eq, times = 10)
#'   
#' summary(cyca_proj)
#' 
#' @export
stage_weight <- function(mpm, stage2 = NA, stage1 = NA, age2 = NA, value = 1.0) {
  
  mpmrows <- stage2_id <- stage1_id <- start_vec <- full_length <- NULL
  
  if (all(!is(mpm, "lefkoMat"))) {
    stop("A regular lefkoMat object is required as input.", call. = FALSE)
  }
  
  if (all(is.na(stage2)) & all(is.na(age2))) {
    stop("Options stage2 and age2 cannot both be set to NA.", call. = FALSE)
  }
  if (all(is.null(stage2)) & all(is.null(age2))) {
    stop("Options stage2 and age2 cannot both be empty.", call. = FALSE)
  }
  
  if (!is.element("stage", names(mpm$ahstages))) {
    stop("Stageframe appears to be modified. Please make sure that a stage
      column exists holding stage names.", call. = FALSE)
  }
  
  if (all(is.na(mpm$hstages)) | all(is.null(mpm$hstages))) {
    historical <- FALSE
  } else {
    historical <- TRUE
  }
  if (all(is.na(mpm$agestages)) | all(is.null(mpm$agestages))) {
    agebystage <- FALSE
  } else {
    agebystage <- TRUE
  }
  
  if (historical & all(is.na(stage1))) {
    stop("Historical projection analysis requires that stage in time t-1 be
      designated for all stage pairs.", call. = FALSE)
  } else if (!historical & !all(is.na(stage1))) {
    stop("Ahistorical projection analysis cannot include designated stages in
      time t-1.", call. = FALSE)
  }
  
  full_length <- max(length(stage2), length(stage1), length(age2), length(value))
  
  if (length(value) == 1 & full_length > 1) {
    value <- rep(value, full_length)
  }
  
  if (length(stage2) == 1 & full_length > 1) {
    stage2 <- rep(stage2, full_length)
  }
  
  if (length(stage1) == 1 & full_length > 1) {
    stage1 <- rep(stage1, full_length)
  }
  
  if (length(age2) == 1 & full_length > 1) {
    age2 <- rep(age2, full_length)
  }
  
  if ((all(is.na(stage2)) | all(is.null(stage2))) & (all(is.na(age2)) | all(is.null(age2)))) {
    stop("Either stage2 or age2 must be provided.", call. = FALSE)
  }
  
  if (all(is.character(stage2))) {
    unknown_stage2 <- which(!is.element(tolower(stage2), c(tolower(mpm$ahstages$stage),
        c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "obs", "nobs"))))
    if (length(unknown_stage2) > 0) {
      stop(paste0("Unknown stage designations used in stage2: ",
        stage2[unknown_stage2]), call. = FALSE)
    }
    
    reassessed <- apply(as.matrix(c(1:length(stage2))), 1, function(X) {
      if (!is.na(stage2[X])) {
        if (is.element(stage2[X], mpm$ahstages$stage)) {
          shrubbery.small <- cbind.data.frame(stage2 = stage2[X], stage1 = stage1[X],
            age2 = age2[X], value = value[X], stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (is.element(stage2[X], as.character(mpm$ahstages$stage_id))) {
          shrubbery.small <- cbind.data.frame(stage2 = as.numeric(stage2[X]),
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "rep") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "all") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage,
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nrep") {
          shrubbery.small <- cbind.data.frame(
            stage2 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                which(mpm$ahstages$matstatus == 1))],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) =="mat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "immat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "prop") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "npr") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "obs") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nobs") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        }
      }
    })
    shrubbery <- do.call(rbind.data.frame, reassessed)
    
  } else if (all(is.numeric(stage2)) & !any(is.na(stage2))) {
    stage2_id <- stage2
    
    if (any(stage2_id > max(mpm$ahstages$stage_id)) | any(stage2_id < min(mpm$ahstages$stage_id))) {
      stop("Unknown stage2 codes used.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(stage2_id), 1, function(X) {
      return(mpm$ahstages$stage[X])
    })
    
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
    
  } else if (all(is.numeric(age2)) & !any(is.na(age2))) {
    if (!(all(is.na(stage2)))) {
      stop("Leslie MPMs should not be entered with the stage2 option set to
        values other than NA.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(age2), 1, function(X) {
      cross_ref <- which(mpm$ahstages$stage_id == X)
      return(mpm$ahstages$stage[cross_ref])
    })
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
    
  } else {
    stop("Input stage2 codes do not conform to accepted inputs.", call. = FALSE)
  }
  
  if (historical) {
    if (all(is.character(shrubbery$stage1)) & !all(is.na(shrubbery$stage1))) {
      
      unknown_stage1 <- which(!is.element(tolower(stage1), c(tolower(mpm$ahstages$stage),
          c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "obs", "nobs", "almostborn"))))
      if (length(unknown_stage1) > 0) {
      stop(paste0("Unknown stage designations used in stage1: ",
        stage1[unknown_stage1]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$stage1[X])) {
          if (is.element(shrubbery$stage1[X], mpm$ahstages$stage)) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (is.element(stage1[X], as.character(mpm$ahstages$stage_id))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = as.numeric(shrubbery$stage1[X]), age2 = age2[X],
              value = value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "rep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage, age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nrep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                  which(mpm$ahstages$matstatus == 1))],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) =="mat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "immat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "prop") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "npr") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "obs") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nobs") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$obsstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      })
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
      
    } else if (all(is.numeric(shrubbery$stage1)) & !any(is.na(shrubbery$stage1))) {
      stage1_id <- shrubbery$stage1
      
      if (any(stage1_id > max(mpm$ahstages$stage_id)) | any(stage1_id < min(mpm$ahstages$stage_id))) {
        stop("Unknown stage1 codes used.", call. = FALSE)
      }
      
      stage1 <- apply(as.matrix(stage1_id), 1, function(X) {
        return(mpm$ahstages$stage[X])
      })
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (agebystage) {
    if (all(is.character(shrubbery$age2)) & !all(is.na(shrubbery$age2))) {
      
      unknown_age2 <- which(!is.element(tolower(age2), c(tolower(mpm$agestages$age),
          "all")))
      if (length(unknown_age2) > 0) {
        stop(paste0("Unknown age designations used in age2: ",
          stage1[unknown_age2]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$age2[X])) {
          common_ages <- unique(mpm$agestages$age[which(mpm$agestages$stage == shrubbery$stage2[X])])
          
          if (is.element(shrubbery$age2[X], as.character(common_ages))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = as.numeric(shrubbery$age2[X]),
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$age2[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = common_ages,
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      }) 
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
    } else if (all(is.numeric(shrubbery$age2)) & !any(is.na(shrubbery$age2))) {
      if (any(shrubbery$age2 > max(mpm$agestages$age)) | any(shrubbery$age2 < min(mpm$agestages$age))) {
        stop("Unknown age2 values used.", call. = FALSE)
      }
      
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = shrubbery$stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (!all(is.numeric(value))) {
    stop("Object value must be composed only of valid numbers.", call. = FALSE)
  }
  
  shrubbery$stage2_id <- apply(as.matrix(shrubbery$stage2), 1, function(X) {
    return(mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)])
  })
  shrubbery$stage1_id <- apply(as.matrix(shrubbery$stage1), 1, function(X) {
    possible_option <- mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)]
    if (length(possible_option) > 0) return(possible_option) else return(NA)
  })
  
  full_length <- dim(shrubbery)[1]
  
  if (!historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$ahstages)[1]) {
      stop("This ahistorical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- shrubbery$stage2_id
    
  } else if (agebystage & !historical) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$agestages)[1]) {
      stop("This age-by-stage mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    if (any(is.na(shrubbery$age2)) | any(!is.numeric(shrubbery$age2))) {
      stop("Option age2 must include only numbers for age-by-stage MPMs.",
        call. = FALSE)
    }
    
    if (any(shrubbery$age2 < min(mpm$agestages$age)) | any(shrubbery$age2 > max(mpm$agestages$age))) {
      stop("Option age2 can only take ages shown in element agestages within the input MPM.",
        call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$agestages$stage_id == shrubbery$stage2_id[X])
      vec1 <- which(mpm$agestages$age == shrubbery$age2[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else if (historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$hstages)[1]) {
      stop("This historical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$hstages$stage_id_2 == shrubbery$stage2_id[X])
      vec1 <- which(mpm$hstages$stage_id_1 == shrubbery$stage1_id[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else {
    stop("Format of mpm not recognized.", call. = FALSE)
  }
  
  output_tab <- cbind.data.frame(shrubbery$stage2, shrubbery$stage2_id,
    shrubbery$stage1, shrubbery$stage1_id, shrubbery$age2, start_vec,
    shrubbery$value, stringsAsFactors = FALSE)
  
  names(output_tab) <- c("stage2", "stage_id_2", "stage1", "stage_id_1", "age2",
    "row_num", "value")
  
  if (!historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2")])
  } else if (historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "stage1", "stage_id_1")])
  } else if (!historical & agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "age2")])
  }
  if (dim(out_check)[1] < dim(output_tab)[1]) {
    warning("Some stages, stage-pairs, or age-stages appear to be listed
      multiple times. This may cause errors in analysis.", call. = FALSE)
  }

  class(output_tab) <- append(class(output_tab), "lefkoEq")
  
  return(output_tab)
}

