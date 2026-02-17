#' @title Historical and Ahistorical Population Projection Matrix Analysis
#' 
#' @description This package creates population matrix projection models (MPMs)
#' for use in population ecological analyses. It presents a complete working
#' environment for the construction and analysis of ALL kinds of MPMs and IPMs,
#' including age, stage, and age-by-stage versions. Its specialty is the
#' estimation of historical MPMs, which are 2d matrices comprising 3 monitoring
#' occasions (2 time steps or periods) of demographic information. The package
#' constructs both function-based and raw MPMs for both standard ahistorical
#' (i.e. 2 occasions, 1 time step) and historical analyses, has functions for
#' complex density-dependent and independent, and stochastic and cyclical,
#' projections, and also includes the automatic calculation of quality control
#' metrics throughout every step of analysis. It also includes powerful
#' functions to standardize demographic datasets.
#' 
#' @details The lefko3 package provides seven categories of functions:
#' 
#' 1. Data transformation and handling functions
#' 
#' 2. Functions determining population characteristics from vertical data
#' 
#' 3. Model building and selection
#' 
#' 4. Matrix / integral projection model creation functions
#' 
#' 5. Population dynamics analysis and projection functions
#' 
#' 6. Functions describing, summarizing, or visualizing MPMs and derived
#' structures
#' 
#' 7. Extra functions used to illustrate core theory and ideas.
#' 
#' @details lefko3 also includes example datasets complete with sample code.
#' 
#' @aliases lefko3-package
"_PACKAGE"
#' @author Richard P. Shefferson <cdorm@g.ecc.u-tokyo.ac.jp>
#' @author Johan Ehrlén
#' @references Shefferson, R.P., J. Ehrlen, and S. Kurokawa. 2021. 
#' \emph{lefko3}: analyzing individual history through size-classified matrix 
#' population models. \emph{Methods in Ecology and Evolution} 12(2): 378-382.
#' @import Rcpp
#' @importFrom glmmTMB glmmTMB nbinom2 truncated_nbinom2 truncated_poisson
#' @importFrom graphics lines
#' @importFrom grDevices palette xy.coords
#' @importFrom lme4 glmer lmer
#' @importFrom MASS glm.nb
#' @importFrom Matrix colSums
#' @importFrom methods as is
#' @importFrom MuMIn dredge
#' @importFrom pscl zeroinfl
#' @importFrom rlang expr
#' @importFrom Rcpp evalCpp
#' @importFrom stats getCall glm lm na.action na.fail na.omit rnorm sd setNames xtabs
#' @importFrom stats as.formula median pchisq poisson var logLik complete.cases
#' @importFrom VGAM posnegbinomial pospoisson vglm predictvglm
#' @useDynLib lefko3
#' @name lefko3-package
NULL
