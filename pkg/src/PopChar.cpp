#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;


// Index of functions
// 1. List sf_create  Create Stageframe for Population Matrix Projection Analysis
// 2. List actualstage3  Calculate Actual Stage, Age, Stage-Pair, or Age-Stage Distributions
// 3. DataFrame density_reassess  Check and Reorganize Density Input Table Into Usable Format
// 4. DataFrame density_input  Set Density Dependence Relationships in Matrix Elements
// 5. List supplemental  Create a Data Frame of Supplemental Data for MPM Development
// 6. void edit_lM_single  Edit Single lefkoMat Object based on Supplemental Data
// 7. List edit_lM  Edit lefkoMat or lefkoMatList Object based on Supplemental Data


//' Create Stageframe for Population Matrix Projection Analysis
//' 
//' Function \code{sf_create()} returns a data frame describing each ahistorical
//' life history stage in the life history model. This data frame can be used as 
//' input into MPM creation functions including \code{\link{flefko3}()}, 
//' \code{\link{flefko2}()}, \code{\link{aflefko2}()}, \code{\link{rlefko3}()},
//' \code{\link{rlefko2}()}, and \code{\link{arlefko2}()}, in which it
//' determines how each stage is treated during matrix estimation.
//' 
//' @name sf_create
//' 
//' @param sizes A numeric vector of the typical or representative size of each
//' life history stage. If making function-based MPMs, then this may be a
//' vector composed of the midpoints of each size bin, or simply of sizes
//' characteristic of the size bins. If denoting the boundary of an automated
//' size classification group, then should denote the absolute minimum size of
//' that group, or the absolute size of that group (see \code{Notes}).
//' @param stagenames A vector of stage names, in the same order as elements in
//' sizes. Can also be set to \code{ipm} for automated size classification (see
//' \code{Notes} section).
//' @param sizesb An optional numeric vector for a second size metric for each
//' life history stage. Only to be used if stages are defined by at least two
//' size metrics in all cases. Same issues apply as in \code{sizes}.
//' @param sizesc An optional numeric vector for a third size metric for each
//' life history stage. Only to be used if stages are defined by at least three
//' size metrics in all cases. Same issues apply as in \code{sizes}.
//' @param repstatus A vector denoting the binomial reproductive status of each
//' life history stage. Defaults to \code{1}.
//' @param obsstatus A vector denoting the binomial observation status of each
//' life history stage. Defaults to \code{1}, but may be changed for
//' unobservable stages.
//' @param propstatus A vector denoting whether each life history stage is a 
//' propagule. Such stages are generally only used in fecundity estimation. 
//' Defaults to \code{0}.
//' @param matstatus A vector denoting whether each stage is mature. Must be
//' composed of binomial values if given. Defaults to 1 for all stages defined 
//' in \code{sizes}.
//' @param immstatus A vector denoting whether each stage is immature. Must be
//' composed of binomial values if given. Defaults to the complement of vector
//' \code{matstatus}.
//' @param minage An optional vector denoting the minimum age at which a stage
//' can occur. Only used in age x stage matrix development. Defaults to
//' \code{NA}.
//' @param maxage An optional vector denoting the maximum age at which a stage
//' should occur. Only used in age x stage matrix development. Defaults to
//' \code{NA}.
//' @param indataset A vector designating which stages are found within the 
//' dataset. While \code{\link{rlefko2}()} and \code{\link{rlefko3}()} can use
//' all stages in the input dataset, \code{\link{flefko3}()} and
//' \code{\link{flefko2}()} can only handle size-classified stages with
//' non-overlapping combinations of size and status variables. Stages that do
//' not actually exist within the dataset should be marked as \code{0} in this
//' vector.
//' @param sizemin A vector giving the absolute minimum values corresponding to
//' each size in the \code{sizes} vector. Requires associated values for
//' \code{sizemax} if used. Only required if not using \code{binhalfwidth}.
//' @param sizebmin A vector giving the absolute minimum values corresponding to
//' each size in the \code{sizesb} vector. Requires associated values for
//' \code{sizebmax} if used. Only required if not using \code{binhalfwidthb}.
//' @param sizecmin A vector giving the absolute minimum values corresponding to
//' each size in the \code{sizesc} vector. Requires associated values for
//' \code{sizecmax} if used. Only required if not using \code{binhalfwidthc}.
//' @param sizemax A vector giving the absolute maximum values corresponding to
//' each size in the \code{sizes} vector. Requires associated values for
//' \code{sizemin} if used. Only required if not using \code{binhalfwidth}.
//' @param sizebmax A vector giving the absolute maximum values corresponding to
//' each size in the \code{sizesb} vector. Requires associated values for
//' \code{sizebmin} if used. Only required if not using \code{binhalfwidthb}.
//' @param sizecmax A vector giving the absolute maximum values corresponding to
//' each size in the \code{sizesc} vector. Requires associated values for
//' \code{sizecmin} if used. Only required if not using \code{binhalfwidthc}.
//' @param binhalfwidth A numeric vector giving the half-width of size bins.
//' Required if \code{sizemin} and \code{sizemax} are not used. Defaults to
//' \code{0.5} for all sizes.
//' @param binhalfwidthb A numeric vector giving the half-width of size bins
//' used for the optional second size metric. Required if \code{sizebmin} and
//' \code{sizebmax} are not used but two or three size classes are used.
//' Defaults to \code{0.5} for all sizes.
//' @param binhalfwidthc A numeric vector giving the half-width of size bins
//' used for the optional third size metric. Required if \code{sizecmin} and
//' \code{sizecmax} are not used but three size classes are used. Defaults to
//' \code{0.5} for all sizes.
//' @param group An integer vector providing information on each respective
//' stage's size classification group. If used, then function-based MPM creation
//' functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
//' \code{\link{aflefko2}()} will estimate transitions only within these groups
//' and for allowed cross-group transitions noted within the supplement table.
//' Defaults to \code{0}.
//' @param comments An optional vector of text entries holding useful text
//' descriptions of all stages.
//' @param roundsize This parameter sets the precision of size classification,
//' and equals the number of digits used in rounding sizes. Defaults to
//' \code{5}.
//' @param roundsizeb This parameter sets the precision of size classification
//' in the optional second size metric, and equals the number of digits used in
//' rounding sizes. Defaults to \code{5}.
//' @param roundsizec This parameter sets the precision of size classification
//' in the optional third size metric, and equals the number of digits used in
//' rounding sizes. Defaults to \code{5}.
//' @param ipmbins An integer giving the number of size bins to create using the
//' primary size classification variable. This number is in addition to any
//' stages that are not size classified. Defaults to \code{100}, and numbers
//' greater than this yield a warning about the loss of statistical power and
//' increasing chance of matrix over-parameterization resulting from increasing
//' numbers of stages.
//' @param ipmbinsb An optional integer giving the number of size bins to create
//' using the secondary size classification variable. This number is in addition
//' to any stages that are not size classified, as well as in addition to any
//' automated size classification using the primary and tertiary size variables.
//' Defaults to \code{NA}, and must be set to a positive integer for automated
//' size classification to progress.
//' @param ipmbinsc An optional integer giving the number of size bins to create
//' using the tertiary size classification variable. This number is in addition
//' to any stages that are not size classified, as well as in addition to any
//' automated size classification using the primary and secondary size
//' variables. Defaults to \code{NA}, and must be set to a positive integer for
//' automated size classification to progress.
//' 
//' @return A data frame of class \code{stageframe}, which includes information
//' on the stage name, size, reproductive status, observation status, propagule 
//' status, immaturity status, maturity status, presence within the core dataset, 
//' stage group classification, raw bin half-width, and the minimum, 
//' center, and maximum of each size bin, as well as its width. If minimum and
//' maximum ages were specified, then these are also included. Also includes an 
//' empty string variable that can be used to describe stages meaningfully. This
//' object can be used as the \code{stageframe} input for \code{\link{flefko3}()} 
//' \code{\link{flefko2}()}, \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.
//' 
//' Variables in this data frame include the following:
//' \item{stage}{The unique names of the stages to be analyzed.}
//' \item{size}{The typical or representative size at which each stage occurs.}
//' \item{size_b}{Size at which each stage occurs in terms of a second size
//' variable, if one exists.}
//' \item{size_c}{Size at which each stage occurs in terms of a third size
//' variable, if one exists.}
//' \item{min_age}{The minimum age at which the stage may occur.}
//' \item{max_age}{The maximum age at which the stage may occur.}
//' \item{repstatus}{A binomial variable showing whether each stage is
//' reproductive.}
//' \item{obsstatus}{A binomial variable showing whether each stage is
//' observable.}
//' \item{propstatus}{A binomial variable showing whether each stage is a
//' propagule.}
//' \item{immstatus}{A binomial variable showing whether each stage can occur as
//' immature.}
//' \item{matstatus}{A binomial variable showing whether each stage occurs in
//' maturity.}
//' \item{indataset}{A binomial variable describing whether each stage occurs in
//' the input dataset.}
//' \item{binhalfwidth_raw}{The half-width of the size bin, as input.}
//' \item{sizebin_min}{The minimum size at which the stage may occur.}
//' \item{sizebin_max}{The maximum size at which the stage may occur.}
//' \item{sizebin_center}{The midpoint of the size bin at which the stage may
//' occur.}
//' \item{sizebin_width}{The width of the size bin corresponding to the stage.}
//' \item{binhalfwidthb_raw}{The half-width of the size bin of a second size
//' variable, as input.}
//' \item{sizebinb_min}{The minimum size at which the stage may occur.}
//' \item{sizebinb_max}{The maximum size at which the stage may occur.}
//' \item{sizebinb_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a second size variable.}
//' \item{sizebinb_width}{The width of the size bin corresponding to the stage,
//' in terms of a second size variable.}
//' \item{binhalfwidthc_raw}{The half-width of the size bin of a third size
//' variable, as input.}
//' \item{sizebinc_min}{The minimum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_max}{The maximum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a third size variable.}
//' \item{sizebinc_width}{The width of the size bin corresponding to the stage,
//' in terms of a third size variable.}
//' \item{group}{An integer denoting the size classification group that the
//' stage falls within.}
//' \item{comments}{A text field for stage descriptions.}
//' 
//' @section Notes:
//' Vectors used to create a stageframe may not mix \code{NA} values with
//' non-\code{NA} values.
//' 
//' If an IPM or function-based matrix with automated size classification is
//' desired, then two stages that occur within the dataset and represent the
//' lower and upper size limits of the IPM must be marked with \code{ipm} in
//' the stagenames vector. These stages should have all characteristics other
//' than size equal, and the size input for whichever size will be classified
//' automatically must include the minimum in one stage and the maximum in the
//' other. The actual characteristics of the first stage encountered in the
//' inputs will be used as the template for the creation of these sizes. Note
//' that \code{ipm} refers to size classification with the primary size
//' variable. To automate size classification with the secondary size variable,
//' use \code{ipmb}, and to automate size classification with the tertiary size
//' variable, use \code{ipmc}. To nest automated size classifications, use 
//' \code{ipmab} for the primary and secondary size variables, \code{ipmac} for
//' the primary and tertiary size variables, \code{ipmbc} for the secondary and
//' tertiary size variables, and \code{ipmabc} for all three size variables.
//' The primary size variable can also be set with \code{ipma}.
//' 
//' If two or more groups of stages, each with its own characteristics, are to
//' be developed for an IPM or function-based MPM, then an even number of stages
//' with two stages marking the minimum and maximum size of each group should be
//' marked with the same code as given above, with all other characteristics
//' equal within each group.
//' 
//' Stage classification groups set with the \code{group} variable create zones
//' within function-based matrices in which survival transitions are estimated.
//' These groups should not be set if transitions are possible between all
//' stages regardless of group. To denote specific transitions as estimable
//' between stage groups, use the \code{\link{supplemental}()} function.
//' 
//' If importing an IPM rather than building one with \code{lefko3}: Using the
//' \code{vrm_import} approach to building function-based MPMs with provided
//' linear model slope coefficients requires careful attention to the
//' stageframe. Although no hfv data frame needs to be entered in this instance,
//' stages for which vital rates are to be estimated via linear models
//' parameterized with coefficients provided via function
//' \code{\link{vrm_import}()} should be marked as occurring within the dataset.
//' Stages for which the provided coefficients should not be used should be
//' marked as not occurring within the dataset.
//' 
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' # Cypripedium example
//' data(cypdata)
//' 
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' @export sf_create
// [[Rcpp::export]]
Rcpp::List sf_create (NumericVector sizes,
  Nullable<StringVector> stagenames = R_NilValue,
  Nullable<NumericVector> sizesb = R_NilValue,
  Nullable<NumericVector> sizesc = R_NilValue,
  Nullable<IntegerVector> repstatus = R_NilValue,
  Nullable<IntegerVector> obsstatus = R_NilValue,
  Nullable<IntegerVector> propstatus = R_NilValue,
  Nullable<IntegerVector> matstatus = R_NilValue,
  Nullable<IntegerVector> immstatus = R_NilValue,
  Nullable<NumericVector> minage = R_NilValue,
  Nullable<NumericVector> maxage = R_NilValue,
  Nullable<IntegerVector> indataset = R_NilValue,
  Nullable<NumericVector> sizemin = R_NilValue,
  Nullable<NumericVector> sizebmin = R_NilValue,
  Nullable<NumericVector> sizecmin = R_NilValue,
  Nullable<NumericVector> sizemax = R_NilValue,
  Nullable<NumericVector> sizebmax = R_NilValue,
  Nullable<NumericVector> sizecmax = R_NilValue,
  Nullable<NumericVector> binhalfwidth = R_NilValue,
  Nullable<NumericVector> binhalfwidthb = R_NilValue,
  Nullable<NumericVector> binhalfwidthc = R_NilValue,
  Nullable<IntegerVector> group = R_NilValue,
  Nullable<StringVector> comments = R_NilValue, int roundsize = 5,
  int roundsizeb = 5, int roundsizec = 5, int ipmbins = 100,
  int ipmbinsb = NA_INTEGER, int ipmbinsc = NA_INTEGER) {
  
  Rcpp::List output_longlist(29);
  Rcpp::CharacterVector varnames {"stage", "size", "size_b", "size_c",
    "min_age", "max_age", "repstatus", "obsstatus", "propstatus", "immstatus",
    "matstatus", "indataset", "binhalfwidth_raw", "sizebin_min", "sizebin_max",
    "sizebin_center", "sizebin_width", "binhalfwidthb_raw", "sizebinb_min",
    "sizebinb_max", "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw",
    "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
    "group", "comments"};
  
  int matsize = static_cast<int>(sizes.size()); // Core vector size
  int used_sizes {1};
  arma::uvec ipm_calls_a (matsize, fill::zeros);
  arma::uvec ipm_calls_b (matsize, fill::zeros);
  arma::uvec ipm_calls_c (matsize, fill::zeros);
  arma::uvec ipm_calls_ab (matsize, fill::zeros);
  arma::uvec ipm_calls_bc (matsize, fill::zeros);
  arma::uvec ipm_calls_ac (matsize, fill::zeros);
  arma::uvec ipm_calls_abc (matsize, fill::zeros);
  
  StringVector stagenames_true (matsize, NA_STRING);
  NumericVector sizesb_true (matsize, NA_REAL);
  NumericVector sizesc_true (matsize, NA_REAL);
  IntegerVector repstatus_true (matsize, 1);
  IntegerVector obsstatus_true (matsize, 1);
  IntegerVector propstatus_true (matsize, 0);
  IntegerVector matstatus_true (matsize, 1);
  IntegerVector immstatus_true (matsize, 0);
  IntegerVector indataset_true (matsize, 1);
  NumericVector minage_true (matsize, NA_REAL);
  NumericVector maxage_true (matsize, NA_REAL);
  NumericVector binhalfwidth_true (matsize, 0.5);
  NumericVector binhalfwidthb_true (matsize, NA_REAL);
  NumericVector binhalfwidthc_true (matsize, NA_REAL);
  NumericVector sizemin_true (matsize, NA_REAL);
  NumericVector sizemax_true (matsize, NA_REAL);
  NumericVector sizebmin_true (matsize, NA_REAL);
  NumericVector sizebmax_true (matsize, NA_REAL);
  NumericVector sizecmin_true (matsize, NA_REAL);
  NumericVector sizecmax_true (matsize, NA_REAL);
  NumericVector sizebin_min (matsize, NA_REAL);
  NumericVector sizebin_max (matsize, NA_REAL);
  NumericVector sizebin_center (matsize, NA_REAL);
  NumericVector sizebin_width (matsize, NA_REAL);
  NumericVector sizebinb_min (matsize, NA_REAL);
  NumericVector sizebinb_max (matsize, NA_REAL);
  NumericVector sizebinb_center (matsize, NA_REAL);
  NumericVector sizebinb_width (matsize, NA_REAL);
  NumericVector sizebinc_min (matsize, NA_REAL);
  NumericVector sizebinc_max (matsize, NA_REAL);
  NumericVector sizebinc_center (matsize, NA_REAL);
  NumericVector sizebinc_width (matsize, NA_REAL);
  IntegerVector group_true (matsize, 0);
  StringVector comments_true (matsize, "No description");
  
  bool size_manual {false};
  bool sizeb_manual {false};
  bool sizec_manual {false};
  bool binhalfwidth_used {false};
  bool binhalfwidthb_used {false};
  bool binhalfwidthc_used {false};
  bool sizemin_used {false};
  bool sizebmin_used {false};
  bool sizecmin_used {false};
  bool sizemax_used {false};
  bool sizebmax_used {false};
  bool sizecmax_used {false};
  
  if (stagenames.isNotNull()) {
    Rcpp::StringVector stagenames_thru(stagenames);
    
    if (stagenames_thru.length() == matsize) {
      stagenames_true = stagenames_thru;
      
      StringVector st_t(stagenames_thru.size());
      
      std::transform(stagenames_thru.begin(), stagenames_thru.end(), st_t.begin(), 
        make_string_transformer(tolower));
      
      for (int i = 0; i < static_cast<int>(stagenames_thru.size()); i++) {
        String check_elem = trimws(st_t(i));
        
        if (check_elem == "ipm" || check_elem == "ipma" || check_elem == "ipm_a" ||
          check_elem == "ipm1" || check_elem == "ipm_1") {
          ipm_calls_a(i) = 1;
        } else if (check_elem == "ipmb" || check_elem == "ipm_b" ||
          check_elem == "ipm2" || check_elem == "ipm_2") {
          ipm_calls_b(i) = 1;
        } else if (check_elem == "ipmc" || check_elem == "ipm_c" ||
          check_elem == "ipm3" || check_elem == "ipm_3") {
          ipm_calls_c(i) = 1;
        } else if (check_elem == "ipmab" || check_elem == "ipm_ab" ||
          check_elem == "ipm12" || check_elem == "ipm_12") {
          ipm_calls_ab(i) = 1;
        } else if (check_elem == "ipmac" || check_elem == "ipm_ac" ||
          check_elem == "ipm13" || check_elem == "ipm_13") {
          ipm_calls_ac(i) = 1;
        } else if (check_elem == "ipmbc" || check_elem == "ipm_bc" ||
          check_elem == "ipm23" || check_elem == "ipm_23") {
          ipm_calls_bc(i) = 1;
        } else if (check_elem == "ipmabc" || check_elem == "ipm_abc" ||
          check_elem == "ipm123" || check_elem == "ipm_123") {
          ipm_calls_abc(i) = 1;
        }
      }
    } else if (stagenames_thru.length() == 1) {
      
      if (!stagenames_thru.is_na(0)) {
        for (int i = 0; i < matsize; i++) {
          Rcpp::String part1(stagenames_thru(0));
          part1 += (static_cast<char>(i + 1));
          
          stagenames_true(i) = part1;
        }
      } else {
        for (int i = 0; i < matsize; i++) {
          Rcpp::String part1("Stage");
          part1 += (static_cast<char>(i + 1));
          
          stagenames_true(i) = part1;
        }
      }
    } else {
      throw Rcpp::exception("Vector stagenames should be same length as vector sizes.",
        false);
    }
  } else {
    for (int i = 0; i < matsize; i++) {
      Rcpp::String part1("Stage");
      part1 += (static_cast<char>(i + 1));
      
      stagenames_true(i) = part1;
    }
  }
  
  if (sizesb.isNotNull()) {
    Rcpp::NumericVector sizesb_thru(sizesb);
    used_sizes++;
    
    if (sizesb_thru.length() == matsize) {
      sizesb_true = sizesb_thru;
    } else if (sizesb_thru.length() == 1) {
      NumericVector try_size (matsize, sizesb_thru(0));
      sizesb_true = try_size;
    } else {
      throw Rcpp::exception("Vector sizesb should be same length as vector sizes.",
        false);
    }
  }
  if (sizesc.isNotNull()) {
    Rcpp::NumericVector sizesc_thru(sizesc);
    used_sizes++;
    
    if (sizesc_thru.length() == matsize) {
      sizesc_true = sizesc_thru;
    } else if (sizesc_thru.length() == 1) {
      NumericVector try_size (matsize, sizesc_thru(0));
      sizesc_true = try_size;
    } else {
      throw Rcpp::exception("Vector sizesc should be same length as vector sizes.",
        false);
    }
    
    if (used_sizes != 3) {
      throw Rcpp::exception("Only use vector sizesc if vector sizesb is also used.",
        false);
    }
  }
  if (minage.isNotNull()) {
    Rcpp::NumericVector minage_thru(minage);
    
    if (minage_thru.length() == matsize) {
      minage_true = minage_thru;
    } else if (minage_thru.length() == 1) {
      NumericVector try_mna (matsize, minage_thru(0));
      minage_true = try_mna;
    } else {
      throw Rcpp::exception("Vector minage should be same length as vector sizes.",
        false);
    }
  }
  if (maxage.isNotNull()) {
    Rcpp::NumericVector maxage_thru(maxage);
    
    if (maxage_thru.length() == matsize) {
      maxage_true = maxage_thru;
    } else if (maxage_thru.length() == 1) {
      NumericVector try_mxa (matsize, maxage_thru(0));
      maxage_true = try_mxa;
    } else {
      throw Rcpp::exception("Vector maxage should be same length as vector sizes.",
        false);
    }
  }
  
  if (repstatus.isNotNull()) {
    Rcpp::IntegerVector repstatus_thru(repstatus);
    
    if (repstatus_thru.length() == matsize) {
      repstatus_true = repstatus_thru;
    } else if (repstatus_thru.length() == 1) {
      IntegerVector try_rep (matsize, repstatus_thru(0));
      repstatus_true = try_rep;
    } else {
      throw Rcpp::exception("Vector repstatus should be same length as vector sizes.",
        false);
    }
    
    if (max(repstatus_true) > 1 || min(repstatus_true) < 0) {
      throw Rcpp::exception("Vector repstatus should be composed only of 0s and 1s.",
        false);
    }
  }
  if (obsstatus.isNotNull()) {
    Rcpp::IntegerVector obsstatus_thru(obsstatus);
    
    if (obsstatus_thru.length() == matsize) {
      obsstatus_true = obsstatus_thru;
    } else if (obsstatus_thru.length() == 1) {
      IntegerVector try_obs (matsize, obsstatus_thru(0));
      obsstatus_true = try_obs;
    } else {
      throw Rcpp::exception("Vector obsstatus should be same length as vector sizes.",
        false);
    }
    
    if (max(obsstatus_true) > 1 || min(obsstatus_true) < 0) {
      throw Rcpp::exception("Vector obsstatus should be composed only of 0s and 1s.",
        false);
    }
  }
  if (propstatus.isNotNull()) {
    Rcpp::IntegerVector propstatus_thru(propstatus);
    
    if (propstatus_thru.length() == matsize) {
      propstatus_true = propstatus_thru;
    } else if (propstatus_thru.length() == 1) {
      IntegerVector try_prop (matsize, propstatus_thru(0));
      propstatus_true = try_prop;
    } else {
      throw Rcpp::exception("Vector propstatus should be same length as vector sizes.",
        false);
    }
    
    if (max(propstatus_true) > 1 || min(propstatus_true) < 0) {
      throw Rcpp::exception("Vector propstatus should be composed only of 0s and 1s.",
        false);
    }
  }
  if (matstatus.isNotNull()) {
    Rcpp::IntegerVector matstatus_thru(matstatus);
    
    if (matstatus_thru.length() == matsize) {
      matstatus_true = matstatus_thru;
    } else if (matstatus_thru.length() == 1) {
      IntegerVector try_mat (matsize, matstatus_thru(0));
      matstatus_true = try_mat;
    } else {
      throw Rcpp::exception("Vector matstatus should be same length as vector sizes.",
        false);
    }
    
    if (max(matstatus_true) > 1 || min(matstatus_true) < 0) {
      throw Rcpp::exception("Vector matstatus should be composed only of 0s and 1s.",
        false);
    }
  }
  if (immstatus.isNotNull()) {
    Rcpp::IntegerVector immstatus_thru(immstatus);
    
    if (immstatus_thru.length() == matsize) {
      immstatus_true = immstatus_thru;
    } else if (immstatus_thru.length() == 1) {
      IntegerVector try_imm (matsize, immstatus_thru(0));
      immstatus_true = try_imm;
    } else {
      throw Rcpp::exception("Vector immstatus should be same length as vector sizes.",
        false);
    }
    
    if (max(immstatus_true) > 1 || min(immstatus_true) < 0) {
      throw Rcpp::exception("Vector immstatus should be composed only of 0s and 1s.",
        false);
    }
  } else {
    IntegerVector trialones(matsize, 1);
    immstatus_true = trialones - matstatus_true;
  }
  if (indataset.isNotNull()) {
    Rcpp::IntegerVector indataset_thru(indataset);
    
    if (indataset_thru.length() == matsize) {
      indataset_true = indataset_thru;
    } else if (indataset_thru.length() == 1) {
      IntegerVector try_ind (matsize, indataset_thru(0));
      indataset_true = try_ind;
    } else {
      throw Rcpp::exception("Vector indataset should be same length as vector sizes.",
        false);
    }
    
    if (max(indataset_true) > 1 || min(indataset_true) < 0) {
      throw Rcpp::exception("Vector indataset should be composed only of 0s and 1s.",
        false);
    }
  }
  
  if (sizemin.isNotNull()) {
    Rcpp::NumericVector sizemin_thru(sizemin);
    
    if (sizemin_thru.length() == matsize) {
      sizemin_true = sizemin_thru;
    } else if (sizemin_thru.length() == 1) {
      NumericVector try_szmna (matsize, sizemin_thru(0));
      sizemin_true = try_szmna;
    } else {
        String ems = "Vector sizemin should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (!sizemax.isNotNull()) {
      throw Rcpp::exception("Vector sizemax is required if sizemin if input.",
        false);
    }
    
    sizemin_used = true;
    size_manual = true;
  }
  if (sizemax.isNotNull()) {
    Rcpp::NumericVector sizemax_thru(sizemax);
    
    if (sizemax_thru.length() == matsize) {
      sizemax_true = sizemax_thru;
    } else if (sizemax_thru.length() == 1) {
      NumericVector try_szmxa (matsize, sizemax_thru(0));
      sizemax_true = try_szmxa;
    } else {
        String ems = "Vector sizemax should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (!sizemin.isNotNull()) {
      throw Rcpp::exception("Vector sizemin is required if sizemax if input.",
        false);
    }
    if (sizemax_true.length() != sizemin_true.length()) {
      throw Rcpp::exception("Vectors sizemax and sizemin should be same length.",
        false);
    }
    
    sizemax_used = true;
  }
  
  if (sizebmin.isNotNull()) {
    Rcpp::NumericVector sizebmin_thru(sizebmin);
    
    if (sizebmin_thru.length() == matsize) {
      sizebmin_true = sizebmin_thru;
    } else if (sizebmin_thru.length() == 1) {
      NumericVector try_szmnb (matsize, sizebmin_thru(0));
      sizebmin_true = try_szmnb;
    } else {
        String ems = "Vector sizebmin should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (!sizebmax.isNotNull()) {
      throw Rcpp::exception("Vector sizebmax is required if sizebmin if input.",
        false);
    }
    
    sizebmin_used = true;
    sizeb_manual = true;
  }
  if (sizebmax.isNotNull()) {
    Rcpp::NumericVector sizebmax_thru(sizebmax);
    
    if (sizebmax_thru.length() == matsize) {
      sizebmax_true = sizebmax_thru;
    } else if (sizebmax_thru.length() == 1) {
      NumericVector try_szmxb (matsize, sizebmax_thru(0));
      sizebmax_true = try_szmxb;
    } else {
        String ems = "Vector sizebmax should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (!sizebmin.isNotNull()) {
      throw Rcpp::exception("Vector sizebmin is required if sizebmax if input.",
        false);
    }
    if (sizebmax_true.length() != sizebmin_true.length()) {
      throw Rcpp::exception("Vectors sizebmax and sizebmin should be same length.",
        false);
    }
    
    sizebmax_used = true;
  }
  
  if (sizecmin.isNotNull()) {
    Rcpp::NumericVector sizecmin_thru(sizecmin);
    
    if (sizecmin_thru.length() == matsize) {
      sizecmin_true = sizecmin_thru;
    } else if (sizecmin_thru.length() == 1) {
      NumericVector try_szmnc (matsize, sizecmin_thru(0));
      sizecmin_true = try_szmnc;
    } else {
        String ems = "Vector sizecmin should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (!sizecmax.isNotNull()) {
      throw Rcpp::exception("Vector sizecmax is required if sizecmin if input.",
        false);
    }
    
    sizec_manual = true;
    sizecmin_used = true;
  }
  if (sizecmax.isNotNull()) {
    Rcpp::NumericVector sizecmax_thru(sizecmax);
    
    if (sizecmax_thru.length() == matsize) {
      sizecmax_true = sizecmax_thru;
    } else if (sizecmax_thru.length() == 1) {
      NumericVector try_szmxc (matsize, sizecmax_thru(0));
      sizecmax_true = try_szmxc;
    } else {
        String ems = "Vector sizecmax should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (!sizecmin.isNotNull()) {
      throw Rcpp::exception("Vector sizecmin is required if sizecmax if input.",
        false);
    }
    if (sizecmax_true.length() != sizecmin_true.length()) {
      throw Rcpp::exception("Vectors sizecmax and sizecmin should be same length.",
        false);
    }
    
    sizecmax_used = true;
  }
  
  if (binhalfwidth.isNotNull()) {
    Rcpp::NumericVector binhalfwidth_thru(binhalfwidth);
    
    if (!size_manual) {
      if (binhalfwidth_thru.length() == matsize) {
        binhalfwidth_true = binhalfwidth_thru;
      } else if (binhalfwidth_thru.length() == 1) {
        NumericVector try_hwa (matsize, binhalfwidth_thru(0));
        binhalfwidth_true = try_hwa;
      } else {
        String ems = "Vector binhalfwidth should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
      }
    } else {
      Rf_warningcall(R_NilValue,
        "Size half-binwidths should not be used with sizemin and sizemax.\n");
    }
    
    binhalfwidth_used = true;
  }
  if (binhalfwidthb.isNotNull()) {
    if (used_sizes == 1) {
      String ems = "Vector binhalfwidthb should only be used if multiple ";
      String ems1 = "size variables are being used for classification.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    Rcpp::NumericVector binhalfwidthb_thru(binhalfwidthb);
    
    if (!sizeb_manual) {
      if (binhalfwidthb_thru.length() == matsize) {
        binhalfwidthb_true = binhalfwidthb_thru;
      } else if (binhalfwidthb_thru.length() == 1) {
        NumericVector try_hwb (matsize, binhalfwidthb_thru(0));
        binhalfwidthb_true = try_hwb;
      } else {
        String ems = "Vector binhalfwidthb should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
      }
    } else {
      Rf_warningcall(R_NilValue,
        "Sizeb half-binwidths should not be used with sizebmin and sizebmax.\n");
    }
    
    binhalfwidthb_used = true;
  } else if (used_sizes > 1 && !sizeb_manual) {
    for (int i = 0; i < binhalfwidthb_true.length(); i++) {
      binhalfwidthb_true(i) = 0.5;
    }
  }
  if (binhalfwidthc.isNotNull()) {
    if (used_sizes < 3) {
      String ems = "Vector binhalfwidthc should only be used if three ";
      String ems1 = "size variables are being used for classification.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    Rcpp::NumericVector binhalfwidthc_thru(binhalfwidthc);
    
    if (!sizec_manual) {
      if (binhalfwidthc_thru.length() == matsize) {
        binhalfwidthc_true = binhalfwidthc_thru;
      } else if (binhalfwidthc_thru.length() == 1) {
        NumericVector try_hwc (matsize, binhalfwidthc_thru(0));
        binhalfwidthc_true = try_hwc;
      } else {
        String ems = "Vector binhalfwidthc should be same ";
        String ems1 = "length as vector sizes.";
        ems += ems1;
        
        throw Rcpp::exception(ems.get_cstring(), false);
      }
    } else {
      Rf_warningcall(R_NilValue,
        "Sizec half-binwidths should not be used with sizecmin and sizecmax.\n");
    }
    
    binhalfwidthc_used = true;
  } else if (used_sizes > 2 && !sizec_manual) {
    for (int i = 0; i < binhalfwidthc_true.length(); i++) {
      binhalfwidthc_true(i) = 0.5;
    }
  }
  
  if (group.isNotNull()) {
    Rcpp::IntegerVector group_thru(group);
    
    if (group_thru.length() == matsize) {
      group_true = group_thru;
    } else if (group_thru.length() == 1) {
      IntegerVector try_grp (matsize, group_thru(0));
      group_true = try_grp;
    } else {
      throw Rcpp::exception("Vector group should be the same length as vector sizes.",
        false);
    }
    
    if (min(group_true) < 0) {
      throw Rcpp::exception("Please use positive integers for group designations.",
        false);
    }
  }
  
  if (comments.isNotNull()) {
    Rcpp::StringVector comments_thru(comments);
    
    if (comments_thru.length() == matsize) {
      comments_true = comments_thru;
    } else if (comments_thru.length() == 1) {
      StringVector try_com (matsize, comments_thru(0));
      comments_true = try_com;
    } else {
      throw Rcpp::exception("Comments vector should be the same length as vector sizes.",
        false);
    }
  }
  
  if (ipmbins < 2) {
    throw Rcpp::exception("Option ipmbins takes only positive integers.", false);
    
  } else if (ipmbins > 100) {
    Rf_warningcall(R_NilValue,
      "High ipmbin numbers may yield overparameterized matrices.\n");
  }
  
  // Now the automated size classification processing
  int ipm_a = sum(ipm_calls_a);
  int ipm_b = sum(ipm_calls_b);
  int ipm_c = sum(ipm_calls_c);
  int ipm_ab = sum(ipm_calls_ab);
  int ipm_bc = sum(ipm_calls_bc);
  int ipm_ac = sum(ipm_calls_ac);
  int ipm_abc = sum(ipm_calls_abc);
  
  // This vector will act as an unique marker for most characteristics
  arma::ivec pair_check (matsize, fill::zeros); 
  arma::uvec entries_to_delete (1);
  int no_entries_to_delete {0};
  
  for (int i = 0; i < matsize; i++) {
    pair_check(i) = repstatus_true(i) * 1000000 + obsstatus_true(i) * 100000 + 
      propstatus_true(i) * 10000 + immstatus_true(i) * 1000 + matstatus_true(i) * 100 +
      indataset_true(i) * 10 + group_true(i);
  }
  
  // Automated size classification with sizea
  if (ipm_a > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      String ems = "Number of bins for auto size classification of ";
      String ems1 = "the primary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_a % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the ";
      String ems1 = "end size, requiring an even number of calls. Calls for auto ";
      String ems2 = "size classification must be matched and not overlapping.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_a); // Shows rows with ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    if (sizemin_used || sizemax_used || binhalfwidth_used) {
      Rf_warningcall(R_NilValue, "Values supplied in vectors sizemin, sizemax, and binhalfwidth will be overwritten for ipm bins.\n");
    }
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minima and maxima. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single ";
          String ems2 = "stages with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize = sizes(check_elems(i));
          minsize = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize = sizes(check_elems(i));
          maxsize = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbins, 0.0);
        StringVector newstagenames (ipmbins, "");
        NumericVector newsizesb (ipmbins, 0.0);
        NumericVector newsizesc (ipmbins, 0.0);
        NumericVector newminage (ipmbins, 0.0);
        NumericVector newmaxage (ipmbins, 0.0);
        IntegerVector newrepstatus (ipmbins, 0);
        IntegerVector newobsstatus (ipmbins, 0);
        IntegerVector newpropstatus (ipmbins, 0);
        IntegerVector newmatstatus (ipmbins, 0);
        IntegerVector newimmstatus (ipmbins, 0);
        IntegerVector newindataset (ipmbins, 0);
        NumericVector newbinhalfwidth (ipmbins, 0.0);
        NumericVector newbinhalfwidthb (ipmbins, 0.0);
        NumericVector newbinhalfwidthc (ipmbins, 0.0);
        IntegerVector newgroup (ipmbins, 0);
        StringVector newcomments (ipmbins, "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbins; j++) {
          newsizes(j) = minsize + (j * standard_increment) + standard_midpoint;
          newsizemin(j) = minsize + (j * standard_increment);
          newsizemax(j) = minsize + ((j + 1) * standard_increment);
          newsizesb(j) = sizesb_true(check_elems(i));
          newsizebmin(j) = sizebmin_true(check_elems(i));
          newsizebmax(j) = sizebmax_true(check_elems(i));
          newsizesc(j) = sizesc_true(check_elems(i));
          newsizecmin(j) = sizecmin_true(check_elems(i));
          newsizecmax(j) = sizecmax_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizes(j));
          newstagenames(j) = "sza_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // IPM classification with sizeb
  if (ipm_b > 0) {
    if (IntegerVector::is_na(ipmbinsb)) {
      String ems = "Number of bins for auto size classification ";
      String ems1 = "of the secondary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_b % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the end ";
      String ems1 = "size, requiring an even number of calls. Calls for auto ";
      String ems2 = "size classification must be matched and not overlap.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_b); // Shows rows with ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    if (sizebmin_used || sizebmax_used || binhalfwidthb_used) {
      Rf_warningcall(R_NilValue, "Values supplied in vectors sizebmin, sizebmax, and binhalfwidthb will be overwritten for ipm bins.\n"); // Rf_warningcall(R_NilValue, emsb.get_cstring());
    }
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }

      if (go_ahead == 1) {
        
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single ";
          String ems2 = "stages with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize = sizesb_true(check_elems(i));
          minsize = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize = sizesb_true(check_elems(i));
          maxsize = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;

        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbinsb, 0.0);
        StringVector newstagenames (ipmbinsb, "");
        NumericVector newsizesb (ipmbinsb, 0.0);
        NumericVector newsizesc (ipmbinsb, 0.0);
        NumericVector newminage (ipmbinsb, 0.0);
        NumericVector newmaxage (ipmbinsb, 0.0);
        IntegerVector newrepstatus (ipmbinsb, 0);
        IntegerVector newobsstatus (ipmbinsb, 0);
        IntegerVector newpropstatus (ipmbinsb, 0);
        IntegerVector newmatstatus (ipmbinsb, 0);
        IntegerVector newimmstatus (ipmbinsb, 0);
        IntegerVector newindataset (ipmbinsb, 0);
        NumericVector newbinhalfwidth (ipmbinsb, 0.0);
        NumericVector newbinhalfwidthb (ipmbinsb, 0.0);
        NumericVector newbinhalfwidthc (ipmbinsb, 0.0);
        IntegerVector newgroup (ipmbinsb, 0);
        StringVector newcomments (ipmbinsb, "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbinsb; j++) {
          newsizesb(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizebmin(j) = minsize + (j * standard_increment);
          newsizebmax(j) = minsize + ((j + 1) * standard_increment);
          newsizes(j) = sizes(check_elems(i));
          newsizemin(j) = sizemin_true(check_elems(i));
          newsizemax(j) = sizemax_true(check_elems(i));
          newsizesc(j) = sizesc_true(check_elems(i));
          newsizecmin(j) = sizecmin_true(check_elems(i));
          newsizecmax(j) = sizecmax_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizesb(j));
          newstagenames(j) = "szb_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // IPM classification with sizec
  if (ipm_c > 0) {
    if (IntegerVector::is_na(ipmbinsc)) {
      String ems = "Number of bins for auto size classification ";
      String ems1 = "of the tertiary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_c % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the end ";
      String ems1 = "size, requiring an even number of calls. Calls for auto ";
      String ems2 = "size classification must be matched and not overlapping.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_c); // Shows rows with ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    if (sizecmin_used || sizecmax_used || binhalfwidthc_used) {
      Rf_warningcall(R_NilValue, "Values supplied in vectors sizecmin, sizecmax, and binhalfwidthc will be overwritten for ipm bins.\n"); // .get_cstring()
    }
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single stages ";
          String ems2 = "with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize = sizesc_true(check_elems(i));
          minsize = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize = sizesc_true(check_elems(i));
          maxsize = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbinsc, 0.0);
        StringVector newstagenames (ipmbinsc, "");
        NumericVector newsizesb (ipmbinsc, 0.0);
        NumericVector newsizesc (ipmbinsc, 0.0);
        NumericVector newminage (ipmbinsc, 0.0);
        NumericVector newmaxage (ipmbinsc, 0.0);
        IntegerVector newrepstatus (ipmbinsc, 0);
        IntegerVector newobsstatus (ipmbinsc, 0);
        IntegerVector newpropstatus (ipmbinsc, 0);
        IntegerVector newmatstatus (ipmbinsc, 0);
        IntegerVector newimmstatus (ipmbinsc, 0);
        IntegerVector newindataset (ipmbinsc, 0);
        NumericVector newbinhalfwidth (ipmbinsc, 0.0);
        NumericVector newbinhalfwidthb (ipmbinsc, 0.0);
        NumericVector newbinhalfwidthc (ipmbinsc, 0.0);
        IntegerVector newgroup (ipmbinsc, 0);
        StringVector newcomments (ipmbinsc, "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbinsc; j++) {
          newsizesc(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizecmin(j) = minsize + (j * standard_increment);
          newsizecmax(j) = minsize + ((j + 1) * standard_increment);
          newsizesb(j) = sizesb_true(check_elems(i));
          newsizebmin(j) = sizebmin_true(check_elems(i));
          newsizebmax(j) = sizebmax_true(check_elems(i));
          newsizes(j) = sizes(check_elems(i));
          newsizemin(j) = sizemin_true(check_elems(i));
          newsizemax(j) = sizemax_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizesc(j));
          newstagenames(j) = "szc_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea & sizeb together
  if (ipm_ab > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      String ems = "Number of bins for auto size classification ";
      String ems1 = "of the primary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (IntegerVector::is_na(ipmbinsb)) {
      String ems = "Number of bins for auto size classification ";
      String ems1 = "of the secondary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_ab % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the end ";
      String ems1 = "size, requiring an even number of calls. Calls for auto ";
      String ems2 = "size classification must be matched and not overlapping.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_ab); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single stages ";
          String ems2 = "with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_b {0};
        double maxsize_b {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsb), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsb), "");
        NumericVector newsizesb ((ipmbins * ipmbinsb), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsb), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsb), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsb), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsb), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsb), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsb), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsb), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsb), 0);
        StringVector newcomments ((ipmbins * ipmbinsb), "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsb; k++) {
            newsizes((j * ipmbinsb) + k) = minsize_a + (j*standard_increment_a) + 
              standard_midpoint_a;
            newsizemin((j * ipmbinsb) + k) = minsize_a + (j * standard_increment_a);
            newsizemax((j * ipmbinsb) + k) = minsize_a + 
              ((j + 1) * standard_increment_a);
            newsizesb((j * ipmbinsb) + k) = minsize_b + (k*standard_increment_b) + 
              standard_midpoint_b;;
            newsizebmin((j * ipmbinsb) + k) = minsize_b + (j * standard_increment_b);
            newsizebmax((j * ipmbinsb) + k) = minsize_b + 
              ((j + 1) * standard_increment_b);
            newsizesc((j * ipmbinsb) + k) = sizesc_true(check_elems(i));
            newsizecmin((j * ipmbinsb) + k) = sizecmin_true(check_elems(i));
            newsizecmax((j * ipmbinsb) + k) = sizecmax_true(check_elems(i));
            newminage((j * ipmbinsb) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsb) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsb) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsb) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsb) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsb) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsb) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsb) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsb) + k) = standard_midpoint_a;
            newbinhalfwidthb((j * ipmbinsb) + k) = standard_midpoint_b;
            newbinhalfwidthc((j * ipmbinsb) + k) = binhalfwidthc_true(check_elems(i));
            newgroup((j * ipmbinsb) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsb) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizes((j * ipmbinsb) + k));
            std::string sizenums2 = std::to_string(newsizesb((j * ipmbinsb) + k));
            newstagenames((j * ipmbinsb) + k) = "sza_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsb) + k) += "_szb_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsb) + k) += "_";
            newstagenames((j * ipmbinsb) + k) += std::to_string(newgroup((j * ipmbinsb) + k));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea & sizec together
  if (ipm_ac > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the primary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the tertiary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_ac % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the end ";
      String ems1 = "size, requiring an even number of calls. Calls for auto ";
      String ems2 = "size classification must be matched and not overlapping.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_ac); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single stages ";
          String ems2 = "with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsc), "");
        NumericVector newsizesb ((ipmbins * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsc), 0);
        StringVector newcomments ((ipmbins * ipmbinsc), "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsc; k++) {
            newsizes((j * ipmbinsc) + k) = minsize_a + (j*standard_increment_a) + 
              standard_midpoint_a;
            newsizemin((j * ipmbinsc) + k) = minsize_a + (j * standard_increment_a);
            newsizemax((j * ipmbinsc) + k) = minsize_a + 
              ((j + 1) * standard_increment_a);
            newsizesc((j * ipmbinsc) + k) = minsize_c + (k*standard_increment_c) + 
              standard_midpoint_c;;
            newsizecmin((j * ipmbinsc) + k) = minsize_c + (j * standard_increment_c);
            newsizecmax((j * ipmbinsc) + k) = minsize_c + 
              ((j + 1) * standard_increment_c);
            newsizesb((j * ipmbinsc) + k) = sizesb_true(check_elems(i));
            newsizebmin((j * ipmbinsc) + k) = sizebmin_true(check_elems(i));
            newsizebmax((j * ipmbinsc) + k) = sizebmax_true(check_elems(i));
            newminage((j * ipmbinsc) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsc) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsc) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsc) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsc) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsc) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsc) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsc) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsc) + k) = standard_midpoint_a;
            newbinhalfwidthb((j * ipmbinsc) + k) = binhalfwidthb_true(check_elems(i));
            newbinhalfwidthc((j * ipmbinsc) + k) = standard_midpoint_c;
            newgroup((j * ipmbinsc) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsc) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizes((j * ipmbinsc) + k));
            std::string sizenums2 = std::to_string(newsizesc((j * ipmbinsc) + k));
            newstagenames((j * ipmbinsc) + k) = "sza_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_szc_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_";
            newstagenames((j * ipmbinsc) + k) += std::to_string(newgroup((j * ipmbinsc)));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizeb & sizec together
  if (ipm_bc > 0) {
    if (IntegerVector::is_na(ipmbinsb)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the primary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the primary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_bc % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the end ";
      String ems1 = "size, requiring an even number of calls. Calls for automated ";
      String ems2 = "size classification must be matched and not overlapping.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_bc); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single stages ";
          String ems2 = "with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_b {0};
        double maxsize_b {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbinsb * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbinsb * ipmbinsc), "");
        NumericVector newsizesb ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbinsb * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbinsb * ipmbinsc), 0);
        StringVector newcomments ((ipmbinsb * ipmbinsc), "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbinsb; j++) {
          for (int k = 0; k < ipmbinsc; k++) {
            newsizesb((j * ipmbinsc) + k) = minsize_b + (j*standard_increment_b) + 
              standard_midpoint_b;
            newsizebmin((j * ipmbinsc) + k) = minsize_b + (j * standard_increment_b);
            newsizebmax((j * ipmbinsc) + k) = minsize_b + 
              ((j + 1) * standard_increment_b);
            newsizesc((j * ipmbinsc) + k) = minsize_c + (k*standard_increment_c) + 
              standard_midpoint_c;;
            newsizecmin((j * ipmbinsc) + k) = minsize_c + (j * standard_increment_c);
            newsizecmax((j * ipmbinsc) + k) = minsize_c + 
              ((j + 1) * standard_increment_c);
            newsizes((j * ipmbinsc) + k) = sizes(check_elems(i));
            newsizemin((j * ipmbinsc) + k) = sizemin_true(check_elems(i));
            newsizemax((j * ipmbinsc) + k) = sizemax_true(check_elems(i));
            newminage((j * ipmbinsc) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsc) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsc) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsc) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsc) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsc) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsc) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsc) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsc) + k) = binhalfwidth_true(check_elems(i));
            newbinhalfwidthb((j * ipmbinsc) + k) = standard_midpoint_b;
            newbinhalfwidthc((j * ipmbinsc) + k) = standard_midpoint_c;
            newgroup((j * ipmbinsc) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsc) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizesb((j * ipmbinsc) + k));
            std::string sizenums2 = std::to_string(newsizesc((j * ipmbinsc) + k));
            newstagenames((j * ipmbinsc) + k) = "szb_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_szc_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_";
            newstagenames((j * ipmbinsc) + k) += std::to_string(newgroup((j * ipmbinsc) + k));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea, sizeb, & sizec together
  if (ipm_abc > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the primary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (IntegerVector::is_na(ipmbinsb)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the secondary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      String ems = "Number of bins for automated size classification ";
      String ems1 = "of the tertiary size variable has not been set.";
      ems += ems1;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    if (ipm_abc % 2 != 0) {
      String ems = "The ipm designation must specify both the start size and the end ";
      String ems1 = "size, requiring an even number of calls. Calls for automated ";
      String ems2 = "size classification must be matched and not overlapping.";
      ems += ems1;
      ems += ems2;
      
      throw Rcpp::exception(ems.get_cstring(), false);
    }
    
    arma::uvec check_elems = find(ipm_calls_abc); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = static_cast<int>(called_pair_check.n_elem);
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. More than 2 ";
          String ems2 = "stages with the same characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        } else if (match_count == 0) {
          String ems = "Stages marked ipm must have equal characteristics in pairs ";
          String ems1 = "only, corresponding to size minimum and maximum. Single stages ";
          String ems2 = "with unique characteristics cannot be marked ipm.";
          ems += ems1;
          ems += ems2;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_b {0};
        double maxsize_b {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.",
            false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = static_cast<int>(entries_to_delete.size());
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsb * ipmbinsc), "");
        NumericVector newsizesb ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsb * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsb * ipmbinsc), 0);
        StringVector newcomments ((ipmbins * ipmbinsb * ipmbinsc), "");
        
        NumericVector newsizemin (ipmbins, NA_REAL);
        NumericVector newsizebmin (ipmbins, NA_REAL);
        NumericVector newsizecmin (ipmbins, NA_REAL);
        NumericVector newsizemax (ipmbins, NA_REAL);
        NumericVector newsizebmax (ipmbins, NA_REAL);
        NumericVector newsizecmax (ipmbins, NA_REAL);
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsb; k++) {
            for (int l = 0; l < ipmbinsc; l++) {
              newsizes((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_a + 
                (j*standard_increment_a) + standard_midpoint_a;
              newsizemin((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_a + 
                (j * standard_increment_a);
              newsizemax((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_a + 
                ((j + 1) * standard_increment_a);
              newsizesb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_b + 
                (k*standard_increment_b) + standard_midpoint_b;
              newsizebmin((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_b + 
                (j * standard_increment_b);
              newsizebmax((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_b + 
                ((j + 1) * standard_increment_b);
              newsizesc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_c + 
                (l*standard_increment_c) + standard_midpoint_c;
              newsizecmin((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_c + 
                (j * standard_increment_c);
              newsizecmax((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_c + 
                ((j + 1) * standard_increment_c);
              newminage((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                minage_true(check_elems(i));
              newmaxage((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                maxage_true(check_elems(i));
              newrepstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                repstatus_true(check_elems(i));
              newobsstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                obsstatus_true(check_elems(i));
              newpropstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                propstatus_true(check_elems(i));
              newmatstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                matstatus_true(check_elems(i));
              newimmstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) =
                immstatus_true(check_elems(i));
              newindataset((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                indataset_true(check_elems(i));
              newbinhalfwidth((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                standard_midpoint_a;
              newbinhalfwidthb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                standard_midpoint_b;
              newbinhalfwidthc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                standard_midpoint_c;
              newgroup((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                group_true(check_elems(i));
              newcomments((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = 
                comments_true(check_elems(i));
              
              std::string sizenums1 = std::to_string(newsizes((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              std::string sizenums2 = std::to_string(newsizesb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              std::string sizenums3 = std::to_string(newsizesc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) =
                "sza_" + sizenums1.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) +=
                "_szb_" + sizenums2.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += 
                "_szc_" + sizenums3.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_";
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += 
                std::to_string(newgroup((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
            }
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizemin_true = concat_dbl(sizemin_true, newsizemin);
        sizemax_true = concat_dbl(sizemax_true, newsizemax);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizebmin_true = concat_dbl(sizebmin_true, newsizebmin);
        sizebmax_true = concat_dbl(sizebmax_true, newsizebmax);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        sizecmin_true = concat_dbl(sizecmin_true, newsizecmin);
        sizecmax_true = concat_dbl(sizecmax_true, newsizecmax);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Delete unneeded elements - run for loop backwards to delete elements properly
  arma::uvec unique_delete = unique(entries_to_delete);
  int no_unique_delete = static_cast<int>(unique_delete.n_elem);
  
  if (no_entries_to_delete > 0) {
    for (int i = 0; i < no_unique_delete; i++) {
      sizes.erase(unique_delete(no_unique_delete - (i + 1)));
      sizesb_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizesc_true.erase(unique_delete(no_unique_delete - (i + 1)));
      stagenames_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      minage_true.erase(unique_delete(no_unique_delete - (i + 1)));
      maxage_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      repstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      obsstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      propstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      matstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      immstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      indataset_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      binhalfwidth_true.erase(unique_delete(no_unique_delete - (i + 1)));
      binhalfwidthb_true.erase(unique_delete(no_unique_delete - (i + 1)));
      binhalfwidthc_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      sizemin_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizebmin_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizecmin_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizemax_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizebmax_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizecmax_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      group_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      comments_true.erase(unique_delete(no_unique_delete - (i + 1)));
    }
  }
  
  // Post-IPM calculations
  matsize = static_cast<int>(sizes.size()); // Redefined length
  int elems_to_grow = matsize - sizebin_min.length();
  NumericVector zeros_to_append (elems_to_grow, NA_REAL);
  
  sizebin_min = concat_dbl(sizebin_min, zeros_to_append);
  sizebin_max = concat_dbl(sizebin_max, zeros_to_append);
  sizebin_center = concat_dbl(sizebin_center, zeros_to_append);
  sizebin_width = concat_dbl(sizebin_width, zeros_to_append);
  
  sizebinb_min = concat_dbl(sizebinb_min, zeros_to_append);
  sizebinb_max = concat_dbl(sizebinb_max, zeros_to_append);
  sizebinb_center = concat_dbl(sizebinb_center, zeros_to_append);
  sizebinb_width = concat_dbl(sizebinb_width, zeros_to_append);
  
  sizebinc_min = concat_dbl(sizebinc_min, zeros_to_append);
  sizebinc_max = concat_dbl(sizebinc_max, zeros_to_append);
  sizebinc_center = concat_dbl(sizebinc_center, zeros_to_append);
  sizebinc_width = concat_dbl(sizebinc_width, zeros_to_append);
  
  for (int i = 0; i < matsize; i++) {
  
    if (!size_manual) {
      sizebin_min(i) = sizes(i) - binhalfwidth_true(i);
      sizebin_max(i) = sizes(i) + binhalfwidth_true(i);
      sizebin_center(i) = sizebin_min(i) + ((sizebin_max(i) - sizebin_min(i)) / 2);
      sizebin_width(i) = sizebin_max(i) - sizebin_min(i);
    } else {
      if (sizemin_true(i) <= sizes(i)) {
        sizebin_min(i) = sizemin_true(i);
      } else {
        throw Rcpp::exception("Some sizemin values are larger than their sizes.",
          false);
      }
      if (sizemax_true(i) >= sizes(i)) {
        sizebin_max(i) = sizemax_true(i);
      } else {
        throw Rcpp::exception("Some sizemax values are smaller than their sizes.",
          false);
      }
      sizebin_center(i) = sizebin_min(i) + ((sizebin_max(i) - sizebin_min(i))/ 2);
      sizebin_width(i) = sizebin_max(i) - sizebin_min(i);
    }
    
    if (used_sizes > 1) {
      if (!sizeb_manual) {
        sizebinb_min(i) = sizesb_true(i) - binhalfwidthb_true(i);
        sizebinb_max(i) = sizesb_true(i) + binhalfwidthb_true(i);
        sizebinb_center(i) = sizebinb_min(i) + ((sizebinb_max(i) - sizebinb_min(i))/ 2);
        sizebinb_width(i) = sizebinb_max(i) - sizebinb_min(i);
      } else {
        if (sizebmin_true(i) <= sizesb_true(i)) {
          sizebinb_min(i) = sizebmin_true(i);
        } else {
          throw Rcpp::exception("Some sizebmin values are larger than their sizes.",
            false);
        }
        if (sizebmax_true(i) >= sizesb_true(i)) {
          sizebinb_max(i) = sizebmax_true(i);
        } else {
          throw Rcpp::exception("Some sizebmax values are smaller than their sizes.",
            false);
        }
        sizebinb_center(i) = sizebinb_min(i) + ((sizebinb_max(i) - sizebinb_min(i))/ 2);
        sizebinb_width(i) = sizebinb_max(i) - sizebinb_min(i);
      }
    }
    if (used_sizes > 2) {
      if (!sizec_manual) {
        sizebinc_min(i) = sizesc_true(i) - binhalfwidthc_true(i);
        sizebinc_max(i) = sizesc_true(i) + binhalfwidthc_true(i);
        sizebinc_center(i) = sizebinc_min(i) + ((sizebinc_max(i) - sizebinc_min(i))/ 2);
        sizebinc_width(i) = sizebinc_max(i) - sizebinc_min(i);
      } else {
        if (sizecmin_true(i) <= sizesc_true(i)) {
          sizebinc_min(i) = sizecmin_true(i);
        } else {
          throw Rcpp::exception("Some sizecmin values are larger than their sizes.",
            false);
        }
        if (sizecmax_true(i) >= sizesc_true(i)) {
          sizebinc_max(i) = sizecmax_true(i);
        } else {
          throw Rcpp::exception("Some sizecmax values are smaller than their sizes.",
            false);
        }
        sizebinc_center(i) = sizebinc_min(i) + ((sizebinc_max(i) - sizebinc_min(i))/ 2);
        sizebinc_width(i) = sizebinc_max(i) - sizebinc_min(i);
      }
    }
  }
  
  sizebin_min = round(sizebin_min, roundsize);
  sizebin_max = round(sizebin_max, roundsize);
  sizebin_center = round(sizebin_center, roundsize);
  sizebin_width = round(sizebin_width, roundsize);
  
  if (used_sizes > 1) {
    sizebinb_min = round(sizebinb_min, roundsizeb);
    sizebinb_max = round(sizebinb_max, roundsizeb);
    sizebinb_center = round(sizebinb_center, roundsizeb);
    sizebinb_width = round(sizebinb_width, roundsizeb);
  }
  if (used_sizes > 2) {
    sizebinc_min = round(sizebinc_min, roundsizec);
    sizebinc_max = round(sizebinc_max, roundsizec);
    sizebinc_center = round(sizebinc_center, roundsizec);
    sizebinc_width = round(sizebinc_width, roundsizec);
  }
  
  output_longlist(0) = stagenames_true;
  output_longlist(1) = sizes;
  output_longlist(2) = sizesb_true;
  output_longlist(3) = sizesc_true;
  
  output_longlist(4) = minage_true;
  output_longlist(5) = maxage_true;
  output_longlist(6) = repstatus_true;
  output_longlist(7) = obsstatus_true;
  output_longlist(8) = propstatus_true;
  output_longlist(9) = immstatus_true;
  output_longlist(10) = matstatus_true;
  
  output_longlist(11) = indataset_true;
  
  output_longlist(12) = binhalfwidth_true;
  output_longlist(13) = sizebin_min;
  output_longlist(14) = sizebin_max;
  output_longlist(15) = sizebin_center;
  output_longlist(16) = sizebin_width;
  
  output_longlist(17) = binhalfwidthb_true;
  output_longlist(18) = sizebinb_min;
  output_longlist(19) = sizebinb_max;
  output_longlist(20) = sizebinb_center;
  output_longlist(21) = sizebinb_width;
  
  output_longlist(22) = binhalfwidthc_true;
  output_longlist(23) = sizebinc_min;
  output_longlist(24) = sizebinc_max;
  output_longlist(25) = sizebinc_center;
  output_longlist(26) = sizebinc_width;
  
  output_longlist(27) = group_true;
  output_longlist(28) = comments_true;
  
  output_longlist.attr("names") = varnames;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, matsize);
  StringVector needed_classes {"data.frame", "stageframe"};
  output_longlist.attr("class") = needed_classes; // data.frame
  
  return output_longlist;
}

//' Calculate Actual Stage, Age, Stage-Pair, or Age-Stage Distributions
//' 
//' Function \code{actualstage3()} shows the frequencies and proportions of
//' each stage, stage pair, age-stage, or age in each year.
//' 
//' @name actualstage3
//' 
//' @param data A demographic dataset in hfv format.
//' @param check_stage A logical value indicating whether to assess frequencies
//' and proportions of stages. Defaults to \code{TRUE}.
//' @param check_age A logical value indicating whether to assess frequencies and
//' proportions of ages. Defaults to \code{FALSE}.
//' @param historical A logical value indicating whether the stage structure
//' should be ahistorical (\code{FALSE}) or historical (\code{TRUE}). Defaults to
//' \code{FALSE}.
//' @param year2 A string value indicating the name of the variable coding for
//' monitoring occasion at time \emph{t}. Defaults to \code{"year2"}.
//' @param indices A vector of three strings, indicating the stage indices for
//' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
//' Defaults to \code{c("stage3index", "stage2index", "stage1index")}.
//' @param stagecol A vector of three strings, indicating the stage name columns
//' for times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
//' Defaults to \code{stagecol = c("stage3", "stage2", "stage1")}.
//' @param agecol A single string indicating the age of individuals in time
//' \emph{t}. Defaults to \code{"obsage"}.
//' @param remove_stage A string vector indicating the names of stages to remove
//' from consideration. Defaults to \code{"NotAlive"}.
//' @param t1_allow A string vector indicating which stages to be removed should
//' be allowed in the stage at time \emph{t}-1 portion of historical stage
//' pairs, if \code{historical = TRUE}. Defaults to \code{"NotAlive"}. Can also
//' be set to \code{"none"}.
//' 
//' @return A data frame with the following variables:
//' \item{rowid}{A string identifier term, equal to the monitoring occasion in
//' time \emph{t} and the stage index.}
//' \item{stageindex}{The stageframe index of the stage. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage}{The name of each stage, or \code{NA}. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage2}{The name of the stage in time \emph{t}. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage1}{The name of the stage in time \emph{t}-1, or \code{NA}. Only
//' output if \code{check_stage = TRUE}.}
//' \item{age}{The age at time \emph{t}. Only output if \code{check_age = TRUE}.}
//' \item{year2}{Monitoring occasion in time \emph{t}.}
//' \item{frequency}{The number of individuals in the respective stage and time.}
//' \item{actual_prop}{The proportion of individuals alive in time \emph{t} in
//' the respective stage.}
//' 
//' @section Notes:
//' This function produces frequencies and proportions of stages in hfv formatted
//' data using stage index variables rather than stage name variables, and so
//' requires the former. The latter is only required if the user wants to know
//' the associated stage names.
//' 
//' Frequencies and proportions will be calculated for all times, including the
//' last time, which is generally found in the \code{stage3} columns of the last
//' \code{year2} entry in object \code{data}. The default is to treat the
//' \code{year2} entry for that time as \code{max(year2) + 1}.
//' 
//' If \code{check_stage = TRUE} and \code{check_age = FALSE}, then this function
//' will assess frequencies and proportions of stages or historical stage-pairs.
//' If both \code{check_stage = TRUE} and \code{check_age = TRUE}, then this
//' function will assess frequencies and proportions of age-stages. If
//' \code{check_stage = FALSE} and \code{check_age = TRUE}, then the frequencies
//' and proportions of ages only will be assessed.
//' 
//' Note that no stageframe is required for this function to operate. Stage
//' names and their order are inferred directly from the object \code{data}.
//' 
//' @examples
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 3, 6, 11, 19.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1.5, 1.5, 3.5, 5)
//' comments <- c("Dormant seed", "1st yr protocorm", "2nd yr protocorm",
//'   "3rd yr protocorm", "Seedling", "Dormant adult",
//'   "Extra small adult (1 shoot)", "Small adult (2-4 shoots)",
//'   "Medium adult (5-7 shoots)", "Large adult (8-14 shoots)",
//'   "Extra large adult (>14 shoots)")
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset, 
//'   binhalfwidth = binvec, comments = comments)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE, age_offset = 4)
//' 
//' all_stage_props_ah <- actualstage3(cypraw_v1)
//' all_stage_props_h <- actualstage3(cypraw_v1, historical = TRUE)
//' all_stage_props_h_NANotAllow <- actualstage3(cypraw_v1, historical = TRUE,
//'   t1_allow = "none")
//' all_stage_props_as <- actualstage3(cypraw_v1, check_age = TRUE)
//' all_age_props <- actualstage3(cypraw_v1, check_stage = FALSE,
//'   check_age = TRUE)
//' 
//' @export actualstage3
// [[Rcpp::export(actualstage3)]]
List actualstage3(RObject data, bool check_stage = true, bool check_age = false,
  bool historical = false, Nullable<RObject> year2 = R_NilValue,
  Nullable<RObject> indices = R_NilValue, Nullable<RObject> stagecol = R_NilValue,
  Nullable<RObject> agecol= R_NilValue, Nullable<RObject> remove_stage = R_NilValue,
  Nullable<RObject> t1_allow = R_NilValue) {
  
  DataFrame data_;
  
  String year2_;
  String agecol_;
  StringVector remove_stage_;
  StringVector t1_allow_;
  
  int year2_num_ {0};
  int agecol_num_ {0};
  IntegerVector remove_stage_num_;
  IntegerVector remove_stage_num_t1a;
  IntegerVector t1_allow_num_;
  
  bool year2_num_yn = false;
  bool agecol_num_yn = false;
  bool remove_stage_num_yn = false;
  bool t1_allow_num_yn = false;
  
  StringVector indices_;
  StringVector stagecol_;
  
  IntegerVector indices_num_;
  IntegerVector stagecol_num_;
  
  bool indices_num_yn = false;
  bool stagecol_num_yn = false;
  
  int stage3index_ {0};
  int stage2index_ {0};
  int stage1index_ {0};
  int stage3_ {0};
  int stage2_ {0};
  int stage1_ {0};
  
  bool stages3_supplied = false;
  bool stages2_supplied = false;
  bool stages1_supplied = false;
  bool indices3_supplied = false;
  bool indices2_supplied = false;
  bool indices1_supplied = false;
  bool ages_supplied = false;
  bool years_supplied = false;
  
  if (!check_age && !check_stage) {
    throw Rcpp::exception("Options check_age and check_stage cannot both be FALSE.",
      false);
  }
  
  if (year2.isNotNull()) {
    if (is<StringVector>(year2)) {
      StringVector year2_long = as<StringVector>(year2);
      year2_ = year2_long(0);
      
    } else if (is<IntegerVector>(year2)) {
      IntegerVector year2_long = as<IntegerVector>(year2);
      year2_num_ = year2_long(0);
      year2_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object year2 must be either a variable name or a column number from object data.", 
        false);
    }
    
  } else year2_ = "year2";
  
  if (agecol.isNotNull()) {
    if (is<StringVector>(agecol)) {
      StringVector agecol_long = as<StringVector>(agecol);
      agecol_ = agecol_long(0);
      
    } else if (is<IntegerVector>(agecol)) {
      IntegerVector agecol_long = as<IntegerVector>(agecol);
      agecol_num_ = agecol_long(0);
      agecol_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object agecol must be either a variable name or a column number from object data.", 
        false);
    }
    
  } else agecol_ = "obsage";
  
  if (remove_stage.isNotNull()) {
    if (is<StringVector>(remove_stage)) {
      remove_stage_ = as<StringVector>(remove_stage);
      
      if (remove_stage_.length() == 1) {
        if (remove_stage_(0) == "none") {
          remove_stage_ = {};
          remove_stage_num_ = {};
          
        }
      }
      
    } else if (is<IntegerVector>(remove_stage)) {
      remove_stage_num_ = as<IntegerVector>(remove_stage);
      remove_stage_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object remove_stage must be either a string of variable names or a vector of column numbers from object data.", 
        false);
    }
    
  } else remove_stage_ = {"NotAlive"};
  
  if (t1_allow.isNotNull()) {
    if (is<StringVector>(t1_allow)) {
      t1_allow_ = as<StringVector>(t1_allow);
      
      if (t1_allow_.length() == 1) {
        if (t1_allow_(0) == "none") {
          t1_allow_ = {};
          t1_allow_num_ = {};
          
        }
      }
      
    } else if (is<IntegerVector>(t1_allow)) {
      t1_allow_num_ = as<IntegerVector>(t1_allow);
      t1_allow_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object t1_allow must be either a string of variable names or a vector of column numbers from object data.",
        false);
    }
  } else t1_allow_ = {"NotAlive"};
  
  if (indices.isNotNull()) {
    if (is<StringVector>(indices)) {
      indices_ = as<StringVector>(indices);
      
    } else if (is<IntegerVector>(indices)) {
      indices_num_ = as<IntegerVector>(indices);
      indices_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object indices must be a vector of variable names or column numbers from object data.", 
        false);
    }
    
  } else indices_ = {"stage3index", "stage2index", "stage1index"};
  
  if (stagecol.isNotNull()) {
    if (is<StringVector>(stagecol)) {
      stagecol_ = as<StringVector>(stagecol);
      
    } else if (is<IntegerVector>(stagecol)) {
      stagecol_num_ = as<IntegerVector>(stagecol);
      stagecol_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object stagecol must be a vector of variable names or column numbers from object data.", 
        false);
    }
    
  } else stagecol_ = {"stage3", "stage2", "stage1"};
  
  if (is<DataFrame>(data)) {
    data_ = as<DataFrame>(data);
    
  } else {
    throw Rcpp::exception("Object data must be a hfv data frame.", false);
  }
  
  // Check data frame for the correct variables, and pull out variables as vectors
  CharacterVector data_vars = data_.attr("names");
  int data_vars_num = data_vars.length();
  int data_rows = data_.nrows();
  int indices_length = indices_.length();
  int stagecol_length = stagecol_.length();
  
  if (historical) {
    if (check_age) {
      throw Rcpp::exception("Package lefko3 does not currently support historical age-by-stage analyses.",
        false);
    }
    
    if (indices_length < 3) {
      throw Rcpp::exception("Object indices must contain stage index variables for times t+1, t, and t-1 if historical = TRUE.",
        false);
    }
    
    if (stagecol_length < 3) {
      throw Rcpp::exception("Object stagecol must contain stage classifications for times t+1, t, and t-1 if historical = TRUE.",
        false);
    }
  } else {
    if (indices_length < 2) {
      throw Rcpp::exception("Object indices must contain stage index variables for times t+1 and t if historical = FALSE.",
        false);
    }
    
    if (stagecol_length < 2) {
      throw Rcpp::exception("Object stagecol must contain stage classifications for times t+1 and t if historical = FALSE.",
        false);
    }
  }
  
  // Loop to find column numbers of variables given as strings
  for (int i = 0; i < data_vars_num; i++) {
    if (!year2_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), year2_)) {
        year2_num_ = i;
        years_supplied = true;
      }
    }
    
    if (!agecol_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), agecol_)) {
        agecol_num_ = i;
        ages_supplied = true;
      }
    }
    
    if (!indices_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(0)))) {
        stage3index_ = i;
        indices3_supplied = true;
      }
      
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(1)))) {
        stage2index_ = i;
        indices2_supplied = true;
      }
      
      if (historical) {
        if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(2)))) {
          stage1index_ = i;
          indices1_supplied = true;
        }
      }
    }
    
    if (!stagecol_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(0)))) {
        stage3_ = i;
        stages3_supplied = true;
      }
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(1)))) {
        stage2_ = i;
        stages2_supplied = true;
      }
      
      if (historical) {
        if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(2)))) {
          stage1_ = i;
          stages1_supplied = true;
        }
      }
    }
  }
  
  // Check if variable column numbers actually fit within data frame provided
  if (year2_num_yn) {
    if (year2_num_ < 0 || year2_num_ >= data_vars_num) {
      throw Rcpp::exception("Variable year2 is not defined within data frame as input.",
        false);
      
    } else {
      years_supplied = true;
    }
  }
  
  if (agecol_num_yn) {
    if (agecol_num_ < 0 || agecol_num_ >= data_vars_num) {
      throw Rcpp::exception("Variable agecol is not defined within data frame as input.",
        false);
      
    } else {
      ages_supplied = true;
    }
  }
  
  if (indices_num_yn) {
    if (indices_num_(0) >= 0 && indices_num_(0) < data_vars_num) {
      stage3index_ = indices_num_(0);
      indices3_supplied = true;
    }
    
    if (indices_num_(1) >= 0 && indices_num_(1) < data_vars_num) {
      stage2index_ = indices_num_(1);
      indices2_supplied = true;
    }
    if (historical) {
      if (indices_num_(2) >= 0 && indices_num_(2) < data_vars_num) {
        stage1index_ = indices_num_(2);
        indices1_supplied = true;
      }
    }
  }
  
  if (stagecol_num_yn) {
    if (stagecol_num_(0) >= 0 && stagecol_num_(0) < data_vars_num) {
      stage3_ = stagecol_num_(0);
      stages3_supplied = true;
    }
    
    if (stagecol_num_(1) >= 0 && stagecol_num_(1) < data_vars_num) {
      stage2_ = stagecol_num_(1);
      stages2_supplied = true;
    }
    
    if (historical) {
      if (stagecol_num_(2) >= 0 && stagecol_num_(2) < data_vars_num) {
        stage1_ = stagecol_num_(2);
        stages1_supplied = true;
      }
    }
  }
  
  // New data frame variable assignments
  IntegerVector data_year2;
  if (years_supplied) data_year2 = as<IntegerVector>(data_[year2_num_]);
  arma::ivec data_year2_arma = as<arma::ivec>(data_year2);
  
  if (static_cast<int>(data_year2.length()) != data_rows) {
    throw Rcpp::exception("Object year2 does not contain a valid variable.", false);
  }
  
  IntegerVector data_agecol;
  
  if (check_age) {
    if (ages_supplied) {
      data_agecol = data_[agecol_];
      
      if (data_agecol.length() != data_rows) {
        throw Rcpp::exception("Object agecol does not contain a valid variable for age at time t.",
          false);
      }
    } else {
      throw Rcpp::exception("Object agecol must be provided if check_age = TRUE.", false);
    }
  }
  
  CharacterVector data_stage3;
  CharacterVector data_stage2;
  CharacterVector data_stage1;
  
  arma::ivec data_stage3index;
  arma::ivec data_stage2index;
  arma::ivec data_stage1index;
  
  if (check_stage) {
    if (stages3_supplied) {
      data_stage3 = as<CharacterVector>(data_[stage3_]);
      
      if (data_stage3.length() != data_rows) {
        throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t+1.",
          false);
      }
    }
    if (stages2_supplied) {
      data_stage2 = as<CharacterVector>(data_[stage2_]);
      
      if (data_stage2.length() != data_rows) {
        throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t.",
          false);
      }
    }
    
    if (indices3_supplied) {
      data_stage3index = as<arma::ivec>(data_[stage3index_]);
      
      if (static_cast<int>(data_stage3index.n_elem) != data_rows) {
        throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t+1.",
          false);
      }
    }
    if (indices2_supplied) {
      data_stage2index = as<arma::ivec>(data_[stage2index_]);
      
      if (static_cast<int>(data_stage2index.n_elem) != data_rows) {
        throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t.",
          false);
      }
    }
    
    if (!stages3_supplied && !indices3_supplied) {
      throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.",
        false);
    }
    if (!stages2_supplied && !indices2_supplied) {
      throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.",
        false);
    }
    
    if (historical) {
      if (stages1_supplied) {
        data_stage1 = data_[stage1_];
        
        if (data_stage1.length() != data_rows) {
          throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t-1.",
            false);
        }
      }
      if (indices3_supplied) {
        data_stage1index = as<arma::ivec>(data_[stage1index_]);
        
        if (static_cast<int>(data_stage1index.n_elem) != data_rows) {
          throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t-1.",
            false);
        }
      }
      
      if (!stages1_supplied && !indices1_supplied) {
        throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.",
          false);
      }
    }
  }
  
  // Now we develop vectors of all values
  IntegerVector all_year2 = sort_unique(data_year2);
  arma::ivec all_year2_arma = as<arma::ivec>(all_year2);
  int year2_num = static_cast<int>(all_year2.length());
  int years_num = year2_num + 1;
  
  IntegerVector all_years (years_num);
  for (int i = 0; i < year2_num; i++) {
    all_years(i) = all_year2(i);
  }
  all_years(year2_num) = all_years(year2_num - 1) + 1;
  arma::ivec all_years_arma = as<arma::ivec>(all_years);
  
  CharacterVector all_stage3;
  CharacterVector all_stage2;
  CharacterVector all_stage1;
  CharacterVector all_stages;
  CharacterVector all_stages_t1a;
  
  arma::ivec all_stage3index;
  arma::ivec all_stage2index;
  arma::ivec all_stage1index;
  arma::ivec all_stageindices;
  arma::ivec all_stageindices_t1a;
  
  int num_stages = 0;
  int num_stages_t1a = 0;
  
  if (check_stage) {
    if (indices3_supplied && indices2_supplied) {
      all_stage3index = unique(data_stage3index);
      all_stage2index = unique(data_stage2index);
      arma::ivec all_stage32indices = unique(join_cols(all_stage3index, all_stage2index));
      
      if (historical) {
        all_stage1index = unique(data_stage1index);
        all_stageindices = unique(join_cols(all_stage32indices, all_stage1index));
      } else {
        all_stageindices = all_stage32indices;
      }
      
      int s_length = static_cast<int>(all_stageindices.n_elem);
      // Rcout << "\n Length of all_stageindices: " << s_length << "\n";
      
      if (stages3_supplied && stages2_supplied) {
        CharacterVector temp_stages (s_length);
        for (int i = 0; i < s_length; i++) {
          arma::uvec found_stages3 = find(data_stage3index == all_stageindices(i));
          
          if (found_stages3.n_elem > 0) {
            temp_stages(i) = data_stage3(found_stages3(0));
          } else {
            arma::uvec found_stages2 = find(data_stage2index == all_stageindices(i));
            
            if (found_stages2.n_elem > 0) {
              temp_stages(i) = data_stage2(found_stages2(0));
              
            } else if (historical && stages1_supplied) {
              arma::uvec found_stages1 = find(data_stage1index == all_stageindices(i));
              
              if (found_stages1.n_elem > 0) {
                temp_stages(i) = data_stage1(found_stages1(0));
              }
            }
          }
        }
        all_stages = temp_stages;
        num_stages = all_stages.length();
        num_stages_t1a = all_stages.length();
      } else {
        num_stages = static_cast<int>(all_stageindices.n_elem);
        num_stages_t1a = static_cast<int>(all_stageindices.n_elem);
      }
    } else if (stages3_supplied && stages2_supplied) {
      all_stage3 = sort_unique(data_stage3);
      all_stage2 = sort_unique(data_stage2);
      
      CharacterVector all_stage32 = union_(all_stage3, all_stage2);
      
      if (historical) {
        all_stage1 = sort_unique(data_stage1);
        
        all_stages = union_(all_stage32, all_stage1);
        
      } else {
        all_stages = all_stage32;
      }
      
      num_stages = all_stages.length();
      num_stages_t1a = all_stages.length();
      IntegerVector temp_stageindices = seq(1, num_stages);
      all_stageindices = as<arma::ivec>(temp_stageindices);
    }
    
    if (remove_stage_.length() != 0 && !remove_stage_num_yn) {
      IntegerVector transfer_remove_stage_num (remove_stage_.length());
      IntegerVector transfer_remove_stage_num_t1a (remove_stage_.length());
      int internal_counter = 0;
      int internal_counter_t1a = 0;
      
      for (int i = 0; i < num_stages; i++) {
        for (int j = 0; j < remove_stage_.length(); j++) {
          if (all_stages(i) == remove_stage_(j)) {
            transfer_remove_stage_num(internal_counter) = i;
            remove_stage_num_yn = true;
            internal_counter++;
            
            if (historical) {
              if (t1_allow_.length() > 0) {
                
                bool found_stage = false;
                for (int k = 0; k < t1_allow_.length(); k++) {
                  if (all_stages(i) == t1_allow_(k)) {
                    found_stage = true;
                    // Rcout << "Found stage (1) " << all_stages(i) << "\n";
                    
                    transfer_remove_stage_num_t1a = shrink(transfer_remove_stage_num_t1a);
                  }
                }
                
                if (!found_stage) {
                  transfer_remove_stage_num_t1a(internal_counter_t1a) = i;
                  internal_counter_t1a++;
                }
              } else if (t1_allow_num_yn) {
                bool found_stage = false;
                for (int k = 0; k < t1_allow_num_.length(); k++) {
                  if (all_stageindices(i) == t1_allow_num_(k)) {
                    found_stage = true;
                    // Rcout << "Found stage (2) " << all_stages(i) << "\n";
                    
                    transfer_remove_stage_num_t1a = shrink(transfer_remove_stage_num_t1a);
                  }
                }
                
                if (!found_stage) {
                  transfer_remove_stage_num_t1a(internal_counter_t1a) = i;
                  internal_counter_t1a++;
                }
              }
            }
          }
        }
      }
      remove_stage_num_ = transfer_remove_stage_num;
      if (historical) {
        remove_stage_num_t1a = clone(transfer_remove_stage_num_t1a);
      } else {
        remove_stage_num_t1a = clone(transfer_remove_stage_num);
      }
      
      if (!remove_stage_num_yn) {
        Rf_warningcall(R_NilValue,
          "Stages in option remove_stage could not be found and will be ignored.");
      }
      
      IntegerVector unique_rsn = sort_unique(remove_stage_num_);
      IntegerVector unique_rsn_t1a = sort_unique(remove_stage_num_t1a);
      
      if (unique_rsn.length() < remove_stage_num_.length() || 
          unique_rsn_t1a.length() < remove_stage_num_t1a.length()) {
        Rf_warningcall(R_NilValue,
          "Some stages in option remove_stage or t1_allow are duplicates and will be ignored.");
        
        remove_stage_num_ = clone(unique_rsn);
        unique_rsn_t1a = clone(unique_rsn_t1a);
      }
      
      // Redetermination of the number of stages to be used
      IntegerVector sorted_remst = remove_stage_num_.sort();
      IntegerVector rev_sorted_remst = rev(sorted_remst);
      
      CharacterVector rep_all_stages = clone(all_stages);
      CharacterVector orig_all_stages = clone(all_stages);
      IntegerVector rep_all_stageindices = as<IntegerVector>(wrap(all_stageindices));
      IntegerVector orig_rep_all_stageindices = clone(rep_all_stageindices);
      
      for (int i = 0; i < rev_sorted_remst.length(); i++) {
        rep_all_stages.erase(rev_sorted_remst(i));
        rep_all_stageindices.erase(rev_sorted_remst(i));
      }
      all_stages = rep_all_stages;
      all_stageindices = as<arma::ivec>(rep_all_stageindices);
      num_stages = all_stages.length();
      
      if (historical) {
        IntegerVector sorted_remst_t1a = remove_stage_num_t1a.sort();
        IntegerVector rev_sorted_remst_t1a = rev(sorted_remst_t1a);
        
        CharacterVector rep_all_stages_t1a = clone(orig_all_stages);
        IntegerVector rep_all_stageindices_t1a = clone(orig_rep_all_stageindices);
        
        for (int i = 0; i < rev_sorted_remst_t1a.length(); i++) {
          rep_all_stages_t1a.erase(rev_sorted_remst_t1a(i));
          rep_all_stageindices_t1a.erase(rev_sorted_remst_t1a(i));
        }
        all_stages_t1a = rep_all_stages_t1a;
        all_stageindices_t1a = as<arma::ivec>(rep_all_stageindices_t1a);
        num_stages_t1a = all_stages_t1a.length();
      }
    } else if (historical) {
      all_stages_t1a = all_stages;
      all_stageindices_t1a = all_stageindices;
      num_stages_t1a = all_stages.length();
    }
  }
  
  IntegerVector all_ages;
  int num_ages = 0;
  
  if (check_age) {
    IntegerVector all_age2 = sort_unique(data_agecol);
    int num_age2 = all_age2.length();
    num_ages = num_age2 + 1;
    
    IntegerVector all_ages_ (num_age2 + 1);
    
    for (int i = 0; i < num_age2; i++) {
      all_ages_(i) = all_age2(i);
    }
    all_ages_(num_age2) = all_ages_(num_age2 - 1) + 1;
    
    all_ages = all_ages_;
  }
  
  // New data frame variables
  int new_df_rows = 0;
  
  if (check_stage && !check_age) {
    if (!historical) {
      new_df_rows = years_num * num_stages;
      
    } else {
      new_df_rows = years_num * num_stages * num_stages_t1a;
    }
    
  } else if (!check_stage && check_age) {
    new_df_rows = years_num * num_ages;
    
  } else if (check_stage && check_age) {
    new_df_rows = years_num * num_stages * num_ages;
  }
  
  IntegerVector new_row_names (new_df_rows);
  CharacterVector new_row_id (new_df_rows);
  IntegerVector new_stageindex (new_df_rows);
  CharacterVector new_stage (new_df_rows);
  CharacterVector new_stage2 (new_df_rows);
  CharacterVector new_stage1 (new_df_rows);
  IntegerVector new_age (new_df_rows);
  IntegerVector new_year (new_df_rows);
  IntegerVector new_frequency (new_df_rows);
  NumericVector new_actualprop (new_df_rows);
  
  if (check_stage && !check_age) {
    if (!historical) {
      IntegerVector first_years = rep(all_years(0), num_stages);
      CharacterVector first_stages2 = clone(all_stages);
      IntegerVector first_indices = as<IntegerVector>(wrap(all_stageindices));
      
      CharacterVector dud_stages (num_stages);
      for (int i = 0; i < num_stages; i++) {
        dud_stages(i) = "";
      }
      
      CharacterVector first_stages1 (num_stages);
      for (int i = 0; i < num_stages; i++) {
        first_stages1(i) = "";
      }
      
      for (int i = 1; i < years_num; i++) {
        IntegerVector next_years = rep(all_years(i), num_stages);
        first_years = concat_int(first_years, next_years);
        
        first_stages2 = concat_str(first_stages2, all_stages);
        first_stages1 = concat_str(first_stages1, dud_stages);
        
        first_indices = concat_int(first_indices, as<IntegerVector>(wrap(all_stageindices)));
      }
      new_year = first_years;
      
      new_stage2 = first_stages2;
      new_stage1 = first_stages1;
      new_stageindex = first_indices;
      
    } else {
      IntegerVector first_years = rep(all_years(0), num_stages * num_stages_t1a);
      CharacterVector first_stages2 = clone(all_stages);
      
      CharacterVector first_stages1_short (num_stages_t1a);
      CharacterVector first_stages1;
      
      for (int i = 0; i < num_stages_t1a; i++) {
        for (int j = 0; j < num_stages_t1a; j++) {
          first_stages1_short(j) = all_stages_t1a(i);
        }
        
        if (i > 0) {
          first_stages1 = concat_str(first_stages1, first_stages1_short);
          first_stages2 = concat_str(first_stages2, all_stages);
        } else {
          first_stages1 = clone(first_stages1_short);
        }
      }
      
      CharacterVector new_first_stages1 = clone(first_stages1);
      CharacterVector new_first_stages2 = clone(first_stages2);
      
      for (int i = 1; i < years_num; i++) {
        IntegerVector next_years = rep(all_years(i), num_stages * num_stages_t1a);
        first_years = concat_int(first_years, next_years);
        
        new_first_stages1 = concat_str(new_first_stages1, first_stages1);
        new_first_stages2 = concat_str(new_first_stages2, first_stages2);
      }
      
      new_year = first_years;
      
      new_stage2 = new_first_stages2;
      new_stage1 = new_first_stages1;
    }
    
  } else if (!check_stage && check_age) {
  
    IntegerVector first_years = rep(all_years(0), num_ages);
    IntegerVector first_ages = clone(all_ages);
    
    for (int i = 1; i < years_num; i++) {
      IntegerVector next_years = rep(all_years(i), num_ages);
      first_years = concat_int(first_years, next_years);
      
      first_ages = concat_int(first_ages, all_ages);
    }
    
    new_year = first_years;
    new_age = first_ages;
    
  } else if (check_stage && check_age) {
    
    CharacterVector first_stages2 = clone(all_stages);
    IntegerVector first_indices = as<IntegerVector>(wrap(all_stageindices));
    
    CharacterVector dud_stages (num_stages);
    for (int i = 0; i < num_stages; i++) {
      dud_stages(i) = "";
    }
    
    CharacterVector first_stages1 (num_stages);
    for (int i = 0; i < num_stages; i++) {
      first_stages1(i) = "";
    }
    
    IntegerVector first_ages = rep(all_ages(0), num_stages);
    for (int i = 1; i < num_ages; i++) {
      IntegerVector next_ages = rep(all_ages(i), num_stages);
      first_ages = concat_int(first_ages, next_ages);
      
      first_stages2 = concat_str(first_stages2, all_stages);
      first_stages1 = concat_str(first_stages1, dud_stages);
    }
    
    IntegerVector first_years = rep(all_years(0), num_stages * num_ages);
    IntegerVector cfirst_ages = clone(first_ages);
    CharacterVector cfirst_stages2 = clone(first_stages2);
    CharacterVector cfirst_stages1 = clone(first_stages1);
    
    for (int i = 1; i < years_num; i++) {
      IntegerVector next_years = rep(all_years(i), num_stages * num_ages);
      first_years = concat_int(first_years, next_years);
      
      first_ages = concat_int(first_ages, cfirst_ages);
      
      first_stages2 = concat_str(first_stages2, cfirst_stages2);
      first_stages1 = concat_str(first_stages1, cfirst_stages1);
      
      first_indices = concat_int(first_indices, as<IntegerVector>(wrap(all_stageindices)));
    }
    
    new_year = first_years;
    new_age = first_ages;
    new_stage2 = first_stages2;
    new_stage1 = first_stages1;
    new_stageindex = first_indices;
    
  }
  
  // Main analysis loop
  IntegerVector pop_by_years (years_num);
  
  for (int i = 0; i < new_df_rows; i++) {
    new_row_names(i) = i + 1;
    
    new_row_id(i) = new_year(i);
    new_row_id(i) += " ";
    
    // Here we set the stage and row id designations
    if (check_stage && !check_age) {
      if (!historical) {
        new_stage(i) = new_stage2(i);
        
      } else {
        new_stage(i) = new_stage2(i);
        new_stage(i) += " ";
        new_stage(i) += new_stage1(i);
        
      }
      new_row_id(i) += new_stage(i);
      
    } else if (!check_stage && check_age) {
      new_row_id(i) += new_age(i);
      
    } else if (check_stage && check_age) {
      new_stage(i) = new_stage2(i);
      
      new_row_id(i) += new_stage(i);
      new_row_id(i) += " ";
      new_row_id(i) += new_age(i);
      
    }
    
    // Now we'll find the frequencies of age-stage-year combos
    arma::uvec year_guys = find(data_year2_arma == new_year(i));
    int year_guys_length = static_cast<int>(year_guys.n_elem);
    
    if (check_stage && !check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (stringcompare_hard(as<std::string>(data_stage2(year_guys(j))), as<std::string>(new_stage2(i)))) {
            if (!historical) {
              arma::uvec year_found = find(all_years_arma == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            } else {
            
              if (stringcompare_hard(as<std::string>(data_stage1(year_guys(j))), as<std::string>(new_stage1(i)))) {
                arma::uvec year_found = find(all_years_arma == new_year(i));
                pop_by_years(year_found(0))++;
                
                new_frequency(i)++;
              }
            }
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2_arma == (new_year(i) - 1));
        int year_guys_length3 = static_cast<int>(year_guys3.n_elem);
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (stringcompare_hard(as<std::string>(data_stage3(year_guys3(j))), as<std::string>(new_stage2(i)))) {
            if (!historical) {
              arma::uvec year_found = find(all_years_arma == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            } else {
              if (stringcompare_hard(as<std::string>(data_stage2(year_guys3(j))), as<std::string>(new_stage1(i)))) {
                arma::uvec year_found = find(all_years_arma == new_year(i));
                pop_by_years(year_found(0))++;
                
                new_frequency(i)++;
              }
            }
          }
        }
      }
      
    } else if (check_stage && check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (stringcompare_hard(as<std::string>(data_stage2(year_guys(j))), as<std::string>(new_stage2(i)))) {
            if (data_agecol(year_guys(j)) == new_age(i)) {
              
              arma::uvec year_found = find(all_years_arma == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            }
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2_arma == (new_year(i) - 1));
        int year_guys_length3 = static_cast<int>(year_guys3.n_elem);
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (stringcompare_hard(as<std::string>(data_stage3(year_guys3(j))), as<std::string>(new_stage2(i)))) {
            if (data_agecol(year_guys3(j)) == new_age(i) - 1) {
              
              arma::uvec year_found = find(all_years_arma == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            }
          }
        }
      }
    } else if (!check_stage && check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (data_agecol(year_guys(j)) == new_age(i)) {
            
            arma::uvec year_found = find(all_years_arma == new_year(i));
            pop_by_years(year_found(0))++;
            
            new_frequency(i)++;
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2_arma == (new_year(i) - 1));
        int year_guys_length3 = static_cast<int>(year_guys3.n_elem);
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (data_agecol(year_guys3(j)) == new_age(i) - 1) {
            
            arma::uvec year_found = find(all_years_arma == new_year(i));
            pop_by_years(year_found(0))++;
            
            new_frequency(i)++;
          }
        }
      }
    }
  }
  
  /*
  Rcout << "Post-loop\n";
  Rcout << "\n Length of new_row_id: " << new_row_id.length() << "\n";
  Rcout << " First entry: " << new_row_id(0) << "\n";
  
  Rcout << "\n Length of new_stage: " << new_stage.length() << "\n";
  Rcout << " First entry: " << new_stage(0) << "\n";
  
  Rcout << "\n Length of new_age: " << new_age.length() << "\n";
  Rcout << " First entry: " << new_age(0) << "\n";
  
  Rcout << "\n Length of new_year: " << new_year.length() << "\n";
  Rcout << " First entry: " << new_year(0) << "\n";
  
  Rcout << "\n Length of new_frequency: " << new_frequency.length() << "\n";
  Rcout << " First entry: " << new_frequency(0) << "\n";
  
  Rcout << "\n Length of new_actualprop: " << new_actualprop.length() << "\n";
  Rcout << " First entry: " << new_actualprop(0) << "\n";
  */
  
  for (int i = 0; i < new_df_rows; i++) {
    arma::uvec year_found = find(all_years_arma == new_year(i));
    
    if (pop_by_years(year_found(0)) > 0) {
      new_actualprop(i) = new_frequency(i) / static_cast<double>(pop_by_years(year_found(0)));
    } else {
      new_actualprop(i) = 0;
    }
  }
  
  
  // Structure the output data frame
  List output;
  
  if (check_stage && !check_age) {
    List raw_output (8);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_stageindex;
    raw_output(2) = new_stage;
    raw_output(3) = new_stage2;
    raw_output(4) = new_stage1;
    raw_output(5) = new_year;
    raw_output(6) = new_frequency;
    raw_output(7) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "stageindex", "stage", "stage2",
      "stage1", "year2", "Freq", "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
    
  } else if (!check_stage && check_age) {
    List raw_output (5);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_age;
    raw_output(2) = new_year;
    raw_output(3) = new_frequency;
    raw_output(4) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "age", "year2", "Freq",
      "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
    
  } else if (check_stage && check_age) {
    List raw_output (9);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_stageindex;
    raw_output(2) = new_stage;
    raw_output(3) = new_stage2;
    raw_output(4) = new_stage1;
    raw_output(5) = new_age;
    raw_output(6) = new_year;
    raw_output(7) = new_frequency;
    raw_output(8) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "stageindex", "stage", "stage2",
      "stage1", "age", "year2", "Freq", "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
  
  }
  
  output.attr("row.names") = new_row_names;
  output.attr("class") = "data.frame";
  return output;
}

//' Check and Reorganize Density Input Table Into Usable Format
//' 
//' Function \code{density_reassess()} takes a density input table as supplied
//' by the \code{\link{density_input}()} function, and checks and rearranges it
//' into a single, complete density input table.
//' 
//' @name density_reassess
//' 
//' @param stageframe The correct stageframe, already modified by
//' \code{\link{sf_reassess}()}.
//' @param dens_inp The density input data frame as is toward the end of
//' \code{\link{density_input}()}.
//' @param agestages The agestages element from the used \code{lefkoMat} object.
//' Only needed if an age-by-stage MPM will be used.
//' @param historical A logical value denoting whether MPM is historical.
//' Defaults to \code{FALSE}.
//' @param agebystage A logical value denoting whether MPM is age-by-stage.
//' Defaults to \code{FALSE}.
//' 
//' @return A corrected density input data frame, usable in density-dependent
//' MPM creation.
//' 
//' @keywords internal
//' @noRd
Rcpp::DataFrame density_reassess(DataFrame stageframe, DataFrame dens_inp,
  Nullable<DataFrame> agestages = R_NilValue, bool historical = false,
  bool agebystage = false) {
  
  StringVector stagevec = as<StringVector>(stageframe["stage"]);
  arma::ivec stageidvec = as<arma::ivec>(stageframe["stage_id"]);
  NumericVector minagevec = as<NumericVector>(stageframe["min_age"]);
  NumericVector maxagevec = as<NumericVector>(stageframe["max_age"]);
  arma::uvec repvec = as<arma::uvec>(stageframe["repstatus"]);
  arma::uvec obsvec = as<arma::uvec>(stageframe["obsstatus"]);
  arma::uvec propvec = as<arma::uvec>(stageframe["propstatus"]);
  arma::uvec immvec = as<arma::uvec>(stageframe["immstatus"]);
  arma::uvec matvec = as<arma::uvec>(stageframe["matstatus"]);
  arma::uvec indvec = as<arma::uvec>(stageframe["indataset"]);
  arma::ivec groupvec = as<arma::ivec>(stageframe["group"]);
  
  int no_stages = static_cast<int>(indvec.n_elem);
  arma::ivec alive(no_stages, fill::ones);

  // Identify all groups in the stageframe
  arma::ivec all_groups = unique(groupvec);
  int no_groups = static_cast<int>(all_groups.n_elem);
  StringVector group_text(no_groups);
  
  for (int i = 0; i < no_groups; i++) {
    group_text(i) = "group";
    group_text(i) += std::to_string(all_groups(i));
  }
  
  // Input variables in current density input frame
  StringVector stage3_di = dens_inp["stage3"];
  StringVector stage2_di = dens_inp["stage2"];
  StringVector stage1_di = dens_inp["stage1"];
  IntegerVector age2_di = dens_inp["age2"];
  IntegerVector style_di = dens_inp["style"];
  IntegerVector time_delay_di = dens_inp["time_delay"];
  NumericVector alpha_di = dens_inp["alpha"];
  NumericVector beta_di = dens_inp["beta"];
  IntegerVector type_di = dens_inp["type"];
  IntegerVector type_t12_di = dens_inp["type_t12"];
  int di_rows = stage3_di.length(); // Size of density input frame
  
  StringVector unique_stages = unique(stagevec);
  StringVector extra_terms = {"rep", "nrep", "immat", "mat", "prop", "npr", "all", "obs", "nobs"};
  
  int no_newstages = unique_stages.length();
  int no_extraterms = extra_terms.length();
  
  StringVector all_possible_stage_terms(no_newstages + no_extraterms + no_groups);
  for (int i = 0; i < no_newstages; i++) {
    all_possible_stage_terms(i) = unique_stages(i);
  }
  for (int i = 0; i < no_extraterms; i++) {
    all_possible_stage_terms(i + no_newstages) = extra_terms(i);
  }
  for (int i = 0; i < no_groups; i++) {
    all_possible_stage_terms(i + no_newstages + no_extraterms) = group_text(i);
  }
  
  DataFrame agestages_;
  arma::uvec agestages_stageid;
  StringVector agestages_stage;
  arma::ivec agestages_age;
  int agestages_rows = 0;
  int age_min = 0;
  
  if (agestages.isNotNull()) {
    if (is<NumericVector>(agestages)) {
      if (agebystage) {
        Rf_warningcall(R_NilValue,
          "Function density_input() requires an agestages object if input MPM is age-by-stage.");
        agebystage = false;
      }
      
    } else if (is<DataFrame>(agestages)) {
      agestages_ = as<DataFrame>(agestages);
      
      if (agestages_.length() == 1) {
        agebystage = false;
        
      } else {
        agestages_stage = as<StringVector>(agestages_["stage"]);
        agestages_stageid = as<arma::uvec>(agestages_["stage_id"]);
        agestages_age = as<arma::ivec>(agestages_["age"]);
        agestages_rows = agestages_stage.length();
        age_min = agestages_age.min();
        
        agebystage = true;
      }
    } else {
      throw Rcpp::exception("Object input as agestages is not recognized.", false);
    }
  }
  
  if (historical && agebystage) {
    throw Rcpp::exception("MPMs cannot be both historical and age-by-stage.", false);
  }
  
  // Check for good entries density input data frame
  for (int i = 0; i < stage3_di.length(); i++) {
    int s3di_count {0};
    int s2di_count {0};
    int s1di_count {0};
    
    for (int j = 0; j < all_possible_stage_terms.length(); j++) {
      std::string s3used = as<std::string>(stage3_di(i));
      std::string s2used = as<std::string>(stage2_di(i));
      std::string s1used = as<std::string>(stage1_di(i));
      
      for (int k = 0; k < static_cast<int>(s3used.size()); k++) {
        s3used[k] = tolower(s3used[k]);
      }
      for (int k = 0; k < static_cast<int>(s3used.size()); k++) {
        s2used[k] = tolower(s2used[k]);
      }
      for (int k = 0; k < static_cast<int>(s3used.size()); k++) {
        s1used[k] = tolower(s1used[k]);
      }
      
      if (stringcompare_hard(s3used, as<std::string>(all_possible_stage_terms(j)))) s3di_count++;
      if (stringcompare_hard(as<std::string>(stage3_di(i)), as<std::string>(all_possible_stage_terms(j)))) s3di_count++;
      
      if (stringcompare_hard(s2used, as<std::string>(all_possible_stage_terms(j)))) s2di_count++;
      if (stringcompare_hard(as<std::string>(stage2_di(i)), as<std::string>(all_possible_stage_terms(j)))) s2di_count++;
      
      if (stringcompare_hard(s3used, "notalive")) {
        throw Rcpp::exception("Stage NotAlive is not allowed.", false);
      }
      if (stringcompare_hard(s2used, "notalive")) {
        throw Rcpp::exception("Stage NotAlive is not allowed.", false);
      }
      
      if (historical) {
        if (stringcompare_hard(s1used, as<std::string>(all_possible_stage_terms(j)))) s1di_count++;
        if (stringcompare_hard(as<std::string>(stage1_di(i)), as<std::string>(all_possible_stage_terms(j)))) s1di_count++;
        
        if (stringcompare_hard(s1used, "notalive")) {
          throw Rcpp::exception("Stage NotAlive is not allowed.", false);
        }
      } 
    }
    
    if (s3di_count == 0) {
      throw Rcpp::exception("Stage names in density input frame (stage3) must match stageframe",
        false);
    }
    if (s2di_count == 0) {
      throw Rcpp::exception("Stage names in density input frame (stage2) must match stageframe",
        false);
    }
    if (historical) {
      if (s1di_count == 0) {
        throw Rcpp::exception("Stage names in density input frame (stage1) must match stageframe",
          false);
      }
    }
  }
    
  IntegerVector s1_calls (di_rows, 1);
  IntegerVector s2_calls (di_rows, 1);
  IntegerVector s3_calls (di_rows, 1);
  IntegerVector s3_planned (di_rows, 1);
  IntegerVector s2_planned (di_rows, 1);
  IntegerVector s1_planned (di_rows, 1);
  
  IntegerVector s123_calls (di_rows, 1);
  
  List age3_calls (di_rows);
  List age2_calls (di_rows);
  
  List stageid3_calls (di_rows);
  List stageid2_calls (di_rows);
  
  arma::uvec prop_stages = find(propvec);
  arma::uvec prop0_stages = find(propvec == 0);
  arma::uvec imm_stages = find(immvec);
  arma::uvec alive_stages = find(alive);
  arma::uvec mat_stages = find(matvec);
  arma::uvec rep_stages = find(repvec);
  arma::uvec rep0_stages = find(repvec == 0);
  arma::uvec mat_rep0_stages = intersect(mat_stages, rep0_stages);
  arma::uvec obs_stages = find(obsvec);
  arma::uvec obs0_stages = find(obsvec == 0);
  arma::uvec all_stages = find(alive); // 7 "all"
  int no_current_group {0};
  
  // Now we build the expanded and edited density input frame
  for (int i = 0; i < di_rows; i++) {
    
    std::string s3used = as<std::string>(stage3_di(i));
    std::string s2used = as<std::string>(stage2_di(i));
    std::string s1used = as<std::string>(stage1_di(i));
    
    for (int j = 0; j < static_cast<int>(s3used.size()); j++) {
      s3used[j] = tolower(s3used[j]);
    }
    for (int j = 0; j < static_cast<int>(s2used.size()); j++) {
      s2used[j] = tolower(s2used[j]);
    }
    for (int j = 0; j < static_cast<int>(s1used.size()); j++) {
      s1used[j] = tolower(s1used[j]);
    }
    
    // Time t+1
    if (stringcompare_hard(s3used, "prop")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(prop_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(prop_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (static_cast<int>(prop_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(prop_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(prop_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (static_cast<int>(prop_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(prop_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "npr")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(prop0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(prop0_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (static_cast<int>(prop0_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(prop0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(prop0_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (static_cast<int>(prop0_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
         s3_calls(i) = static_cast<int>(prop0_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "immat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(imm_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(imm_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (static_cast<int>(imm_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(imm_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(imm_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (static_cast<int>(imm_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(imm_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "mat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(mat_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (static_cast<int>(mat_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (static_cast<int>(mat_stages(j)) == static_cast<int>(agestages_stageid(k) - 1)) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(mat_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
         s3_calls(i) = static_cast<int>(mat_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "rep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(rep_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(rep_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(rep_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "nrep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(mat_rep0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(mat_rep0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(mat_rep0_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "obs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(obs_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(obs_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(obs_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "nobs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(obs0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(obs0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(obs0_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s3used, "all")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(all_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < static_cast<int>(all_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = static_cast<int>(all_stages.n_elem);
      }
      
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage3_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = static_cast<int>(current_group.n_elem);
          
          if (agebystage) {
            int found_stages = 0;
            
            for (int j = 0; j < static_cast<int>(current_group.n_elem); j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (type_di(i) == 1) {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (IntegerVector::is_na(age2_di(i))) {
                      found_stages++;
                    } else if (agestages_age(k) == (age2_di(i) + 1)) {
                      found_stages++;
                    }
                  }
                } else {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (agestages_age(k) == age_min) found_stages++;
                  }
                }
              }
            }
            
            s3_calls(i) = found_stages;
            IntegerVector age3_vec(found_stages);
            IntegerVector stageid3_vec(found_stages);
            int a3_counter = 0;
            
            for (int j = 0; j < static_cast<int>(current_group.n_elem); j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (type_di(i) == 1) {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (IntegerVector::is_na(age2_di(i))) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    } else if (agestages_age(k) == (age2_di(i) + 1)) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    }
                  }
                } else {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (agestages_age(k) == age_min) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    }
                  }
                }
              }
            }
            age3_calls(i) = age3_vec;
            stageid3_calls(i) = stageid3_vec;
            
          } else {
            s3_calls(i) = no_current_group;
          }
        }
      }
    }
    
    // Correct for unnecessary rows
    int cv3_len {0};
    
    RObject check_stageid3_calls_i = as<RObject>(stageid3_calls(i));
    if (!check_stageid3_calls_i.isNULL()) {
      IntegerVector check_vec_3 = as<IntegerVector>(check_stageid3_calls_i);
      cv3_len = static_cast<int>(check_vec_3.length());
    }
    if (s3_calls(i) == 0 && cv3_len > 0) s3_calls(i) = 1;

    // Time t
    if (stringcompare_hard(s2used, "prop")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(prop_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(prop_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(prop_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "npr")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(prop0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(prop0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(prop0_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "immat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(imm_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (imm_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(imm_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (imm_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(imm_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "mat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(mat_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(mat_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(mat_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "rep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(rep_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (rep_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(rep_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (rep_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(rep_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "nrep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(mat_rep0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(mat_rep0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(mat_rep0_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "obs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(obs_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(obs_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(obs_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "nobs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(obs0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(obs0_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(obs0_stages.n_elem);
      }
      
    } else if (stringcompare_hard(s2used, "all")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < static_cast<int>(all_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (all_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < static_cast<int>(all_stages.n_elem); j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (all_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = static_cast<int>(all_stages.n_elem);
      }
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage2_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = static_cast<int>(current_group.n_elem);
          
          if (agebystage) {
            int found_stages = 0;
            
            for (int j = 0; j < static_cast<int>(current_group.n_elem); j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (current_group(j) == agestages_stageid(k) - 1) {
                  if (IntegerVector::is_na(age2_di(i))) {
                    found_stages++;
                  } else if (agestages_age(k) == age2_di(i)) {
                    found_stages++;
                  }
                }
              }
            }
            
            s2_calls(i) = found_stages;
            IntegerVector age2_vec(found_stages);
            IntegerVector stageid2_vec(found_stages);
            int a2_counter = 0;
            
            for (int j = 0; j < static_cast<int>(current_group.n_elem); j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (current_group(j) == agestages_stageid(k) - 1) {
                  if (IntegerVector::is_na(age2_di(i))) {
                    age2_vec(a2_counter) = agestages_age(k);
                    stageid2_vec(a2_counter) = agestages_stageid(k);
                    a2_counter++;
                  } else if (agestages_age(k) == age2_di(i)) {
                    age2_vec(a2_counter) = agestages_age(k);
                    stageid2_vec(a2_counter) = agestages_stageid(k);
                    a2_counter++;
                  }
                }
              }
            }
            age2_calls(i) = age2_vec;
            stageid2_calls(i) = stageid2_vec;
            
          } else {
            s2_calls(i) = no_current_group;
          }
        }
      }
    }
    if (s2_calls(i) == 0) s2_calls(i) = 1;
    
    // Time t-1
    if (stringcompare_hard(s1used, "prop")) {
      s1_calls(i) = static_cast<int>(prop_stages.n_elem);
    } else if (stringcompare_hard(s1used, "npr")) {
      s1_calls(i) = static_cast<int>(prop0_stages.n_elem);
    } else if (stringcompare_hard(s1used, "immat")) {
      s1_calls(i) = static_cast<int>(imm_stages.n_elem);
    } else if (stringcompare_hard(s1used, "mat")) {
      s1_calls(i) = static_cast<int>(mat_stages.n_elem);
    } else if (stringcompare_hard(s1used, "rep")) {
      s1_calls(i) = static_cast<int>(rep_stages.n_elem);
    } else if (stringcompare_hard(s1used, "nrep")) {
      s1_calls(i) = static_cast<int>(mat_rep0_stages.n_elem);
    } else if (stringcompare_hard(s1used, "obs")) {
      s1_calls(i) = static_cast<int>(obs_stages.n_elem);
    } else if (stringcompare_hard(s1used, "nobs")) {
      s1_calls(i) = static_cast<int>(obs0_stages.n_elem);
    } else if (stringcompare_hard(s1used, "all")) {
      s1_calls(i) = static_cast<int>(all_stages.n_elem);
    } else if (StringVector::is_na(stage1_di(i))) {
      s1_calls(i) = 1;
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage1_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = static_cast<int>(current_group.n_elem);
          
          s1_calls(i) = no_current_group;
        }
      }
    }
    if (s1_calls(i) == 0) s1_calls(i) = 1;
    
    s123_calls(i) = s3_calls(i) * s2_calls(i) * s1_calls(i);
  }
  
  NumericVector basepoints(di_rows, 0.0);
  for (int i = 0; i < (di_rows - 1); i++) {
    basepoints(i+1) = basepoints(i) + s123_calls(i);
  }
  
  // New output data frame set-up
  int newdi_rows = sum(s123_calls);
  
  StringVector stage3_newdi(newdi_rows);
  StringVector stage2_newdi(newdi_rows);
  StringVector stage1_newdi(newdi_rows);
  IntegerVector age2_newdi(newdi_rows);
  IntegerVector style_newdi(newdi_rows);
  NumericVector alpha_newdi(newdi_rows);
  NumericVector beta_newdi(newdi_rows);
  IntegerVector time_delay_newdi(newdi_rows);
  IntegerVector type_newdi(newdi_rows);
  IntegerVector type_t12_newdi(newdi_rows);
  
  int overall_counter {0};
  // int group_check {0};
  
  int group_baseline3 {0};
  int group_baseline2 {0};
  int group_baseline1 {0};
  
  int group_ratchet3 {0};
  int group_ratchet2 {0};
  int group_ratchet1 {0};
  
  int prevl3 {0};
  int prevl2 {0};
  int prevl1 {0};
  
  for (int i = 0; i < di_rows; i++) {
    overall_counter = 0;
    
    int age3_counter = 0;
    int age2_counter = 0;
    
    std::string s3used = as<std::string>(stage3_di(i));
    std::string s2used = as<std::string>(stage2_di(i));
    std::string s1used = as<std::string>(stage1_di(i));
    
    for (int j = 0; j < static_cast<int>(s3used.size()); j++) {
      s3used[j] = tolower(s3used[j]);
    }
    for (int j = 0; j < static_cast<int>(s2used.size()); j++) {
      s2used[j] = tolower(s2used[j]);
    }
    for (int j = 0; j < static_cast<int>(s1used.size()); j++) {
      s1used[j] = tolower(s1used[j]);
    }
    
    for (int j = 0; j < s1_calls(i); j++) {
      for (int k = 0; k < s2_calls(i); k++) {
        for (int l = 0; l < s3_calls(i); l++) {
          
          // Time t+1
          if (stringcompare_hard(s3used, "prop")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "npr")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "immat")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "mat")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "rep")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "nrep")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "obs")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "nobs")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "all")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(l));
            }
            
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage3_di(i) == group_text(m)) {
                if (agebystage) {
                  group_check = 1;
                  
                  IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
                  IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
                  int a3v_length = age3_vec.length();
                  
                  stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
                  
                  age3_counter++;
                  if (age3_counter == a3v_length) age3_counter = 0;
                
                } else {
                  if (l == 0) group_ratchet3 = 0;
                  if (l != prevl3 && l != 0) group_ratchet3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(groupvec == m);
                  int current_group_length = static_cast<int>(current_group.n_elem);
                  if (group_ratchet3 > (current_group_length - 1)) {
                    group_ratchet3 = 0;
                  }
                  
                  if (group_ratchet3 == 0) {
                    group_baseline3 = l;
                  }
                  
                  stage3_newdi(basepoints(i) + overall_counter) = 
                    stagevec(current_group(l - group_baseline3));
                  
                  prevl3 = l;
                  
                }
              }
            }
            if (group_check == 0) {
              stage3_newdi(basepoints(i) + overall_counter) = stage3_di(i);
            }
            
            group_check = 0;
          }
          
          // Set up of most core variables in output data frame
          age2_newdi(basepoints(i) + overall_counter) = age2_di(i);
          style_newdi(basepoints(i) + overall_counter) = style_di(i);
          alpha_newdi(basepoints(i) + overall_counter) = alpha_di(i);
          beta_newdi(basepoints(i) + overall_counter) = beta_di(i);
          time_delay_newdi(basepoints(i) + overall_counter) = time_delay_di(i);
          type_newdi(basepoints(i) + overall_counter) = type_di(i);
          type_t12_newdi(basepoints(i) + overall_counter) = type_t12_di(i);
          
          // Time t
          if (stringcompare_hard(s2used, "prop")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(k));
            }
          } else if (stringcompare_hard(s2used, "npr")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "immat")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(k));
            }
          } else if (stringcompare_hard(s2used, "mat")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(k));
            }
          } else if (stringcompare_hard(s2used, "rep")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(k));
            }
          } else if (stringcompare_hard(s2used, "nrep")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "obs")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(k));
            }
          } else if (stringcompare_hard(s2used, "nobs")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "all")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(k));
            }
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage2_di(i) == group_text(m)) {
                if (agebystage) {
                  group_check = 1;
                  
                  IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
                  IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
                  int a2v_length = age2_vec.length();
                  
                  stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
                  age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
                  
                  age2_counter++;
                  if (age2_counter == a2v_length) age2_counter = 0;
                
                } else {
                  if (k == 0) group_ratchet2 = 0;
                  if (k != prevl2 && k != 0) group_ratchet2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(groupvec == m);
                  int current_group_length = static_cast<int>(current_group.n_elem);
                  if (group_ratchet2 > (current_group_length - 1)) {
                    group_ratchet2 = 0;
                  }
                  
                  if (group_ratchet2 == 0) {
                    group_baseline2 = k;
                  }
                  
                  stage2_newdi(basepoints(i) + overall_counter) =
                    stagevec(current_group(k - group_baseline2));
                  
                  prevl2 = k;
                }
              }
            }
             
            if (group_check == 0) {
             stage2_newdi(basepoints(i) + overall_counter) = stage2_di(i);
            }
            
            group_check = 0;
          }
          
          // Time t-1
          if (stringcompare_hard(s1used, "prop")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(j));
          } else if (stringcompare_hard(s1used, "npr")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(j));
          } else if (stringcompare_hard(s1used, "immat")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(j));
          } else if (stringcompare_hard(s1used, "mat")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(j));
          } else if (stringcompare_hard(s1used, "rep")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(j));
          } else if (stringcompare_hard(s1used, "nrep")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(j));
          } else if (stringcompare_hard(s1used, "obs")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(j));
          } else if (stringcompare_hard(s1used, "nobs")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(j));
          } else if (stringcompare_hard(s1used, "all")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(j));
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage1_di(i) == group_text(m)) {
                if (j == 0) group_ratchet1 = 0;
                if (j != prevl1 && j != 0) group_ratchet1 += 1;
                
                group_check = 1;
                arma::uvec current_group = find(groupvec == m);
                int current_group_length = static_cast<int>(current_group.n_elem);
                if (group_ratchet1 > (current_group_length - 1)) {
                  group_ratchet1 = 0;
                }
                
                if (group_ratchet1 == 0) {
                  group_baseline1 = j;
                }
                
                stage1_newdi(basepoints(i) + overall_counter) =
                  stagevec(current_group(j - group_baseline1));
                
                prevl1 = j;
              }
            }
            
            if (group_check == 0) {
              stage1_newdi(basepoints(i) + overall_counter) = stage1_di(i);
            }
            
            group_check = 0;
          }
          
          overall_counter++;
        }
      }
    }
  }
  
  // Output final set-up
  Rcpp::List new_di(10);
  
  new_di(0) = stage3_newdi;
  new_di(1) = stage2_newdi;
  new_di(2) = stage1_newdi;
  new_di(3) = age2_newdi;
  new_di(4) = style_newdi;
  new_di(5) = alpha_newdi;
  new_di(6) = beta_newdi;
  new_di(7) = time_delay_newdi;
  new_di(8) = type_newdi;
  new_di(9) = type_t12_newdi;
  
  CharacterVector namevec = {"stage3", "stage2", "stage1", "age2", "style",
    "alpha", "beta", "time_delay", "type", "type_t12"};
  CharacterVector newclass = {"data.frame"};
  new_di.attr("names") = namevec;
  new_di.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newdi_rows);
  new_di.attr("class") = newclass;
  
  return new_di;
}

//' Set Density Dependence Relationships in Matrix Elements
//' 
//' Function \code{density_input()} provides all necessary data to incorporate
//' density dependence into a \code{lefkoMat} object, a list of matrices, or a
//' single matrix. Five forms of density dependence are allowed, including the
//' Ricker function, the Beverton-Holt function, the Usher function, the
//' logistic function, and the additive limit function. In each case, density
//' must have an effect with a delay of at least one time-step (see Notes). The
//' resulting data frame provides a guide for other \code{lefko3} functions to
//' modify matrix elements by density.
//'
//' @name density_input
//' 
//' @param mpm The \code{lefkoMat} object that will be subject to density
//' dependent projection.
//' @param stage3 A vector showing the name or number of the stage in occasion
//' \emph{t}+1 in the transitions to be affected by density. Abbreviations for
//' groups of stages are also usable (see Notes).
//' @param stage2 A vector showing the name or number of the stage in occasion
//' \emph{t} in the transition to be affected by density. Abbreviations for
//' groups of stages are also usable (see Notes).
//' @param stage1 A vector showing the name or number of the stage in occasion
//' \emph{t}-1 in the transition to be affected by density. Only needed if a
//' historical MPM is used. Abbreviations for groups of stages are also usable
//' (see Notes).
//' @param age2 A vector showing the age of the stage in occasion \emph{t} in the
//' transition to be affected by density. Only needed if an age-by-stage MPM is
//' used.
//' @param style A vector coding for the style of density dependence on each
//' transition subject to density dependence. Options include \code{1},
//' \code{ricker}, \code{ric}, or \code{r} for the Ricker function; \code{2},
//' \code{beverton}, \code{bev}, and \code{b} for the Beverton-Holt function;
//' \code{3}, \code{usher}, \code{ush}, and \code{u} for the Usher function;
//' \code{4}, \code{logistic}, \code{log}, and \code{l} for the logistic
//' function; and \code{5}, \code{additive}, \code{add}, and \code{a} for the
//' additive limit function. If only a single code is provided, then all noted
//' transitions are assumed to be subject to this style of density dependence.
//' Defaults to \code{ricker}.
//' @param time_delay An integer vector indicating the number of occasions back
//' on which density dependence operates. Defaults to \code{1}, and may not equal
//' any integer less than 1. If a single number is input, then all noted
//' transitions are assumed to be subject to this time delay.  Defaults to
//' \code{1}.
//' @param alpha A vector indicating the numeric values to use as the
//' alpha term in the two parameter Ricker, Beverton-Holt, or Usher function, or
//' the value of the carrying capacity \emph{K} to use in the logistic or
//' additive limit functions (see \code{Notes} section for more on this term).
//' If a single number is provided, then all noted transitions are assumed to be
//' subject to this value of alpha. Defaults to \code{1}.
//' @param beta A vector indicating the numeric values to use as the beta term in
//' the two parameter Ricker, Beverton-Holt, or Usher function, or the
//' multiplier on the previous population size in the additive limit function.
//' Also used to indicate whether to use \emph{K} as a hard limit in the
//' logistic equation (see section \code{Notes} below). If a single number is
//' provided, then all noted transitions are assumed to be subject to this value
//' of \code{beta}. Defaults to \code{1}.
//' @param type A vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition; or 
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition. If empty or not provided, then defaults to \code{1} for survival
//' transition.
//' @param type_t12 An optional vector denoting the kind of transition between
//' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
//' deVries format is desired. This should be entered as \code{1}, \code{S}, or
//' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
//' fecundity transitions. Defaults to \code{1} for survival transition, with
//' impacts only on the construction of deVries-format hMPMs.
//' 
//' @return A data frame of class \code{lefkoDens}. This object can be used as
//' input in function \code{\link{projection3}()}.
//' 
//' Variables in this object include the following:
//' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
//' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
//' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced,
//' if applicable.}
//' \item{age2}{Age at occasion \emph{t} in the transition to be replaced, if
//' applicable.}
//' \item{style}{Style of density dependence, coded as 1, 2, 3, or 4 for the
//' Ricker, Beverton-Holt, Usher, or logistic function, respectively.}
//' \item{time_delay}{The time delay on density dependence, in time steps.}
//' \item{alpha}{The value of alpha in the Ricker, Beverton-Holt, or Usher
//' function, or the value of carrying capacity, \emph{K}, in the logistic
//' function.}
//' \item{beta}{The value of beta in the Ricker, Beverton-Holt, or Usher
//' function.}
//' \item{type}{Designates whether the transition from occasion \emph{t} to
//' occasion \emph{t}+1 is a survival transition probability (1), or a fecundity
//' rate (2).}
//' \item{type_t12}{Designates whether the transition from occasion \emph{t}-1 to
//' occasion \emph{t} is a survival transition probability (1), a fecundity rate
//' (2).}
//' 
//' @section Notes:
//' This function provides inputs when density dependence is operationalized
//' directly on matrix elements. It can be used in both \code{projection3()} and
//' \code{f_projection3()}. Users wishing to modify vital rate functions by
//' density dependence functions for use in function-based projections with
//' function \code{f_projection3()} should use function \code{density_vr()} to
//' provide the correct inputs.
//' 
//' The parameters \code{alpha} and \code{beta} are applied according to the
//' two-parameter Ricker function, the two-parameter Beverton-Holt function, the
//' two-parameter Usher function, the one-parameter logistic function, or the
//' two-parameter additive limit function. Although the default is that a 1 time
//' step delay is assumed, greater time delays can be set through the
//' \code{time_delay} option.
//' 
//' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
//' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
//' are to be used, \code{nrep} if all mature but non-reproductive stages are to
//' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
//' immature stages are to be used, \code{prop} if all propagule stages are to be
//' used, \code{npr} if all non-propagule stages are to be used, \code{obs} if
//' all observable stages are to be used, \code{nobs} if all unobservable stages
//' are to be used, and leave empty or use \code{all} if all stages in stageframe
//' are to be used.
//' 
//' When using the logistic function, it is possible that the time delay used in
//' density dependent simulations will cause matrix elements to become negative.
//' To prevent this behavior, set the associated \code{beta} term to \code{1.0}.
//' Doing so will set \code{K} as the hard limit in the logistic equation,
//' essentially setting a minimum limit at \code{0} for all matrix elements
//' modified.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{projection3}()}
//' 
//' @examples
//' \donttest{
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
//'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
//'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
//'   yearcol = "year2", indivcol = "individ")
//' 
//' ehrlen3mean <- lmean(ehrlen3)
//' 
//' e3d <- density_input(ehrlen3mean, stage3 = c("Sd", "Sdl"),
//'   stage2 = c("rep", "rep"), stage1 = c("all", "all"), style = 1,
//'   time_delay = 1, alpha = 1, beta = 0, type = c(2, 2), type_t12 = c(1, 1))
//' 
//' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE, substoch = 2,
//'   density = e3d)
//' }
//' 
//' @export density_input
// [[Rcpp::export(density_input)]]
DataFrame density_input(List mpm, Nullable<RObject> stage3 = R_NilValue,
  Nullable<RObject> stage2 = R_NilValue, Nullable<RObject> stage1 = R_NilValue,
  Nullable<RObject> age2 = R_NilValue, Nullable<RObject> style = R_NilValue,
  Nullable<RObject> time_delay = R_NilValue, Nullable<RObject> alpha = R_NilValue,
  Nullable<RObject> beta = R_NilValue, Nullable<RObject> type = R_NilValue,
  Nullable<RObject> type_t12 = R_NilValue) {
  
  bool historical = false;
  bool agebystage = false;
  bool ageonly = false;
  bool stage3used = false;
  bool stage2used = false;
  
  // Check quality of mpm input
  StringVector mpm_class_vec;
  if (mpm.hasAttribute("class")) {
    mpm_class_vec = mpm.attr("class");
  } else mpm_class_vec = {"list"};
  std::string mpm_class = as<std::string>(mpm_class_vec(0));
  
  CharacterVector mpm_elems = mpm.names();
  int mpm_name_check = 0;
  for (int i = 0; i < mpm_elems.length(); i++) {
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "ahstages")) mpm_name_check++;
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "hstages")) mpm_name_check++;
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "agestages")) mpm_name_check++;
  }
  
  if (!(mpm_name_check == 3 && stringcompare_hard(mpm_class, "lefkoMat"))) {
    throw Rcpp::exception("This function requires a lefkoMat object as input.", false);
  }
  
  DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
  DataFrame hstages = as<DataFrame>(mpm["hstages"]);
  DataFrame agestages = as<DataFrame>(mpm["agestages"]);
  
  if (hstages.length() > 1) {
    historical = true;
  }
  if (agestages.length() > 1) {
    agebystage = true;
  } 
  if (stageframe.length() == 1) {
    throw Rcpp::exception("Input lefkoMat object appears to be missing a stageframe.",
      false);
  }
  
  CharacterVector ahstages_elems = stageframe.names();
  int ahs_name_check = 0;
  for (int i = 0; i < ahstages_elems.length(); i++) {
    if (stringcompare_hard(as<std::string>(ahstages_elems(i)), "stage_id")) ahs_name_check++;
    if (stringcompare_hard(as<std::string>(ahstages_elems(i)), "stage")) ahs_name_check++;
  }
  
  if (ahs_name_check != 2) {
    throw Rcpp::exception("Stageframe appears to be modified. Please make sure that element ahstages in the lefkoMat object includes both a stage column holding stage names and a stage_id column holding unique, stage identifying integers.", 
      false);
  }
  
  IntegerVector stage_id_sf = as<IntegerVector>(stageframe["stage_id"]);
  StringVector stage_sf = as<StringVector>(stageframe["stage"]);
  int no_stages = stage_sf.length();
  
  // Density input vector standardization
  StringVector stage3_names;
  StringVector stage2_names;
  StringVector stage1_names;
  
  if (stage3.isNotNull()) {
    if (is<StringVector>(stage3)) {
      StringVector stage3_names_ = as<StringVector>(stage3);
      stage3_names = stage3_names_;
      
    } else if (is<IntegerVector>(stage3)) {
      arma::ivec stage3_ids = as<arma::ivec>(stage3);
      int stage3_entries = static_cast<int>(stage3_ids.n_elem);
      
      arma::uvec bad_lows = find(stage3_ids < 1);
      arma::uvec bad_highs = find(stage3_ids > no_stages);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
      }
      
      StringVector new_stage3 (stage3_entries);
      for (int i = 0; i < stage3_entries; i++) {
        new_stage3(i) = stage_sf(stage3_ids(i) - 1);
      }
      
      stage3_names = new_stage3;
      
    } else {
      throw Rcpp::exception("Option stage3 is not a recognized input type.", false);
    }
    stage3used = true;
    
  } else {
    if (!age2.isNotNull()) {
      throw Rcpp::exception("This function requires inputs for stage2 and 3 in stage-based MPMs, and inputs for age2 for age-based MPMs. Age-stage MPMs require all three of these.", false);
    }
  }
  
  if (stage2.isNotNull()) {
    if (is<StringVector>(stage2)) {
      StringVector stage2_names_ = as<StringVector>(stage2);
      stage2_names = stage2_names_;
      
    } else if (is<IntegerVector>(stage2)) {
      arma::ivec stage2_ids = as<arma::ivec>(stage2);
      int stage2_entries = static_cast<int>(stage2_ids.n_elem);
      
      arma::uvec bad_lows = find(stage2_ids < 1);
      arma::uvec bad_highs = find(stage2_ids > no_stages);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
      }
      
      StringVector new_stage2 (stage2_entries);
      for (int i = 0; i < stage2_entries; i++) {
        new_stage2(i) = stage_sf(stage2_ids(i) - 1);
      }
      
      stage2_names = new_stage2;
      
    } else {
      throw Rcpp::exception("Option stage2 is not a recognized input type.", false);
    }
    stage2used = true;
  } else {
    if (!age2.isNotNull()) {
      throw Rcpp::exception("This function requires inputs for stage2 and 3 in stage-based MPMs, and inputs for age2 for age-based MPMs. Age-stage MPMs require all three of these.", false);
    }
  }
  
  if (stage1.isNotNull()) {
    if (is<IntegerVector>(stage1)) {
      arma::ivec stage1_ids = as<arma::ivec>(stage1);
      int stage1_entries = static_cast<int>(stage1_ids.n_elem);
      
      arma::uvec bad_lows = find(stage1_ids < 1);
      arma::uvec bad_highs = find(stage1_ids > no_stages);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
      }
      
      StringVector new_stage1 (stage1_entries);
      for (int i = 0; i < stage1_entries; i++) {
        new_stage1(i) = stage_sf(stage1_ids(i) - 1);
      }
      
      stage1_names = new_stage1;
      
    } else if (is<StringVector>(stage1)) {
      StringVector stage1_names_ = as<StringVector>(stage1);
      stage1_names = stage1_names_;
      
    } else {
      throw Rcpp::exception("Option stage1 is not a recognized input type.", false);
    }
    
    if (!historical) {
      throw Rcpp::exception("Do not include option stage1 for ahistorical MPMs.", false);
    }
    
  } else if (historical) {
    throw Rcpp::exception("Option stage1 is needed for historical MPMs.", false);
    
  } else {
    StringVector stage1_names_ = {NA_STRING};
    stage1_names = stage1_names_;
  }
  
  if (stage3_names.length() != stage2_names.length()) {
    throw Rcpp::exception("All transitions to modify require information at least for stage2 and stage3. These inputs must also be of equal length.", 
      false);
  }
  
  if (historical) {
    if (stage3_names.length() != stage1_names.length()) {
      throw Rcpp::exception("All historical transitions to modify require information for stage1, stage2, and stage3. These inputs must also be of equal length.", 
        false);
    }
  }
  
  IntegerVector age2_vec;
  IntegerVector style_vec;
  IntegerVector time_delay_vec;
  NumericVector alpha_vec;
  NumericVector beta_vec;
  IntegerVector type_vec;
  IntegerVector type_t12_vec;
  
  if (age2.isNotNull()) {
    if (!stage2used && !stage3used) ageonly = true;
    
    if (!agebystage && !ageonly) {
      throw Rcpp::exception("Do not use the age2 option unless the MPM is either purely age-based or age-by-stage.", false);
    }
    
    if (is<IntegerVector>(age2)) {
      IntegerVector age2_vec_ = as<IntegerVector>(age2);
      age2_vec = age2_vec_;
      
      arma::ivec age2_arma = as<arma::ivec>(age2_vec_);
      arma::uvec bad_lows = find(age2_arma < 0);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Negative ages are not allowed.", false);
      }
    } else if (is<NumericVector>(age2)) {
      IntegerVector age2_vec_ = as<IntegerVector>(age2);
      age2_vec = age2_vec_;
      
      arma::ivec age2_arma = as<arma::ivec>(age2_vec_);
      arma::uvec bad_lows = find(age2_arma < 0);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Negative ages are not allowed.", false);
      }
    } else {
      throw Rcpp::exception("Only integers are allowed in option age2.", false);
    }
  } else {
    IntegerVector age2_vec_ = {NA_INTEGER};
    age2_vec = age2_vec_;
  }
  
  if (ageonly) {
    DataFrame mpm_ahstages = as<DataFrame>(mpm["ahstages"]);
    IntegerVector age_order = as<IntegerVector>(mpm_ahstages["stage_id"]);
    CharacterVector age_names = as<CharacterVector>(mpm_ahstages["stage"]);
    
    IntegerVector common_ages = intersect(age2_vec, age_order);
    if (common_ages.length() != age2_vec.length()){
      throw Rcpp::exception("All ages used in option age2 must be listed explicitly as ages in the ahstages element of the input mpm.", false);
    }
    
    CharacterVector age2_names (static_cast<int>(age2_vec.length()));
    CharacterVector age3_names (static_cast<int>(age2_vec.length()));
    for (int i = 0; i < static_cast<int>(age2_vec.length()); i++) {
      age2_names(i) = age_names(age2_vec(i) - 1);
      
      if (i < (static_cast<int>(age2_vec.length())) - 1) {
        age3_names(i) = age_names(age2_vec(i) - 1);
      } else {
        age3_names(i) = age_names(age2_vec(i));
      }
    }
    stage3_names = age3_names;
    stage2_names = age2_names;
  }
  
  if (style.isNotNull()) {
    StringVector ricker_style = {"1", "ricker", "ricke", "rick", "ric", "ri", "r"};
    StringVector beverton_style = {"2", "beverton", "beverto", "bevert", "bever",
      "beve", "bev", "be", "b", "holt", "hol", "ho", "h"};
    StringVector usher_style = {"3", "usher", "ushe", "ush", "us", "u"};
    StringVector logistic_style = {"4", "logistic", "logisti", "logist", "logis",
      "logi", "log", "lo", "l"};
    StringVector recr_box_style = {"5", "additive", "additiv", "additi",
      "addit", "addi", "add", "ad", "a"};
    
    if (is<IntegerVector>(style)) {
      IntegerVector style_ = as<IntegerVector>(style);
      style_vec = style_;
      
      arma::ivec style_arma = as<arma::ivec>(style_);
      
      arma::uvec bad_lows = find(style_arma < 1);
      arma::uvec bad_highs = find(style_arma > 5);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Invalid density dependence style entered.", false);
      }
      
    } else if (is<NumericVector>(style)) {
      IntegerVector style_ = as<IntegerVector>(style);
      style_vec = style_;
      
      arma::ivec style_arma = as<arma::ivec>(style_);
      
      arma::uvec bad_lows = find(style_arma < 1);
      arma::uvec bad_highs = find(style_arma > 5);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Invalid density dependence style entered.", false);
      }
      
    } else if (is<StringVector>(style)) {
      StringVector style_stringvec = as<StringVector>(style);
      int style_elems = style_stringvec.length();
      
      IntegerVector style_vec_ (style_elems, 0);
      
      for (int i = 0; i < style_elems; i++) {
        std::string ssv = as<std::string>(style_stringvec(i));
        for (int j = 0; j < static_cast<int>(ssv.size()); j++) {
          ssv[j] = tolower(ssv[j]);
        }
        
        for (int j = 0; j < ricker_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(ricker_style(j)))) {
            style_vec_(i) = 1;
          }
        }
        
        for (int j = 0; j < beverton_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(beverton_style(j)))) {
            style_vec_(i) = 2;
          }
        }
        
        for (int j = 0; j < usher_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(usher_style(j)))) {
            style_vec_(i) = 3;
          }
        }
        
        for (int j = 0; j < logistic_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(logistic_style(j)))) {
            style_vec_(i) = 4;
          }
        }
        
        for (int j = 0; j < recr_box_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(recr_box_style(j)))) {
            style_vec_(i) = 5;
          }
        }
      }
      
      if (min(style_vec_) == 0) {
        throw Rcpp::exception("Some density dependence styles were not recognized.", false);
      }
      style_vec = style_vec_;
      
    } else {
      throw Rcpp::exception("Invalid entry for option style.", false);
      
    }
  } else {
    IntegerVector style_vec_ = {1};
    style_vec = style_vec_;
  }
  
  if (time_delay.isNotNull()) {
    if (is<IntegerVector>(time_delay)) {
      IntegerVector time_delay_vec_ = as<IntegerVector>(time_delay);
      time_delay_vec = time_delay_vec_;
      
      arma::ivec time_delay_arma = as<arma::ivec>(time_delay_vec_);
      arma::uvec bad_lows = find(time_delay_arma < 1);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Time delays less than 1 time step are not allowed.", false);
      }
    } else if (is<NumericVector>(time_delay)) {
      IntegerVector time_delay_vec_ = as<IntegerVector>(time_delay);
      time_delay_vec = time_delay_vec_;
      
      arma::ivec time_delay_arma = as<arma::ivec>(time_delay_vec_);
      arma::uvec bad_lows = find(time_delay_arma < 1);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Time delays less than 1 time step are not allowed.", false);
      }
    } else {
      throw Rcpp::exception("Only integers are allowed in option time_delay.", false);
    }
  } else {
    IntegerVector time_delay_vec_ = {1};
    time_delay_vec = time_delay_vec_;
  }
  
  if (alpha.isNotNull()) {
    if (is<NumericVector>(alpha)) {
      NumericVector alpha_vec_ = as<NumericVector>(alpha);
      alpha_vec = alpha_vec_;
      
    } else {
      throw Rcpp::exception("Only floating point decimals are allowed in option alpha.", false);
    }
    
  } else {
    NumericVector alpha_vec_ = {1.};
    alpha_vec = alpha_vec_;
  }
  
  if (beta.isNotNull()) {
    if (is<NumericVector>(beta)) {
      NumericVector beta_vec_ = as<NumericVector>(beta);
      beta_vec = beta_vec_;
      
    } else {
      throw Rcpp::exception("Only floating point decimals are allowed in option beta.", false);
    }
  } else {
    NumericVector beta_vec_ = {1.};
    beta_vec = beta_vec_;
  }
  
  if (type.isNotNull()) {
    if (is<IntegerVector>(type)) {
      IntegerVector type_vec_ = as<IntegerVector>(type);
      type_vec = type_vec_;
      
      arma::ivec type_arma = as<arma::ivec>(type_vec_);
      arma::uvec bad_lows = find(type_arma < 1);
      arma::uvec bad_highs = find(type_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Transition types may only be type 1 or 2.", false);
      }
    } else if (is<StringVector>(type)) {
      StringVector type_stringvec = as<StringVector>(type);
      int type_elems = type_stringvec.length();
      
      IntegerVector type_vec_ (type_elems, 0);
      
      for (int i = 0; i < type_elems; i++) {
        if (stringcompare_simple(as<std::string>(type_stringvec(i)), "s", true)) {
          type_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "1", true)) {
          type_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "f", true)) {
          type_vec_(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "2", true)) {
          type_vec_(i) = 2;
        } else {
          throw Rcpp::exception("Invalid entry in option type.", false);
        }
      }
      type_vec = type_vec_;
    } else if (is<NumericVector>(type)) {
      IntegerVector type_vec_ = as<IntegerVector>(type);
      type_vec = type_vec_;
      
      arma::ivec type_arma = as<arma::ivec>(type_vec_);
      arma::uvec bad_lows = find(type_arma < 1);
      arma::uvec bad_highs = find(type_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Transition types may only be type 1 or 2.", false);
      }
      
    } else {
      throw Rcpp::exception("Only integers are allowed in option type.", false);
    }
  } else {
    IntegerVector type_vec_ = {1};
    type_vec = type_vec_;
  }
  
  if (type_t12.isNotNull()) {
    if (is<IntegerVector>(type_t12)) {
      IntegerVector type_t12_vec_ = as<IntegerVector>(type_t12);
      type_t12_vec = type_t12_vec_;
      
      arma::ivec type_t12_arma = as<arma::ivec>(type_t12_vec_);
      arma::uvec bad_lows = find(type_t12_arma < 1);
      arma::uvec bad_highs = find(type_t12_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Historical transition types may only be type 1 or 2.", false);
      }
    } else if (is<NumericVector>(type_t12)) {
      IntegerVector type_t12_vec_ = as<IntegerVector>(type_t12);
      type_t12_vec = type_t12_vec_;
      
      arma::ivec type_t12_arma = as<arma::ivec>(type_t12_vec_);
      arma::uvec bad_lows = find(type_t12_arma < 1);
      arma::uvec bad_highs = find(type_t12_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Historical transition types may only be type 1 or 2.", false);
      }
    } else if (is<StringVector>(type_t12)) {
      StringVector type_t12_stringvec = as<StringVector>(type_t12);
      int type_t12_elems = type_t12_stringvec.length();
      
      IntegerVector type_t12_vec_ (type_t12_elems, 0);
      
      for (int i = 0; i < type_t12_elems; i++) {
        if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "s", true)) {
          type_t12_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "1", true)) {
          type_t12_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "f", true)) {
          type_t12_vec_(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "2", true)) {
          type_t12_vec_(i) = 2;
        } else {
          throw Rcpp::exception("Invalid entry in option type_t12.", false);
        }
      }
      type_t12_vec = type_t12_vec_;
      
    } else {
      throw Rcpp::exception("Only integers are allowed in option type_t12.", false);
    }
  } else {
    IntegerVector type_t12_vec_ = {1};
    type_t12_vec = type_t12_vec_;
  }
  
  StringVector new_stage3_names;
  StringVector new_stage2_names;
  StringVector new_stage1_names;
  IntegerVector new_age2_vec;
  IntegerVector new_style_vec;
  IntegerVector new_time_delay_vec;
  NumericVector new_alpha_vec;
  NumericVector new_beta_vec;
  IntegerVector new_type_vec;
  IntegerVector new_type_t12_vec;
  
  IntegerVector vec_lengths = {static_cast<int>(stage3_names.length()),
    static_cast<int>(stage2_names.length()), static_cast<int>(stage1_names.length()),
    static_cast<int>(age2_vec.length()), static_cast<int>(style_vec.length()),
    static_cast<int>(time_delay_vec.length()), static_cast<int>(alpha_vec.length()),
    static_cast<int>(beta_vec.length()), static_cast<int>(type_vec.length()),
    static_cast<int>(type_t12_vec.length())};
  int vec_max = max(vec_lengths);
  
  if (stage3_names.length() != vec_max) {
    if (stage3_names.length() == 1) {
      StringVector ns3n (vec_max, stage3_names(0));
      new_stage3_names = ns3n;
      
    } else {
      throw Rcpp::exception("Vector stage3 is not the correct length.", false);
    }
  } else {
    new_stage3_names = stage3_names;
  }
  
  if (stage2_names.length() != vec_max) {
    if (stage2_names.length() == 1) {
      StringVector ns2n (vec_max, stage2_names(0));
      new_stage2_names = ns2n;
      
    } else {
      throw Rcpp::exception("Vector stage2 is not the correct length.", false);
    }
  } else {
    new_stage2_names = stage2_names;
  }
  
  if (stage1_names.length() != vec_max) {
    if (stage1_names.length() == 1) {
      StringVector ns1n (vec_max, stage1_names(0));
      new_stage1_names = ns1n;
      
    } else {
      throw Rcpp::exception("Vector stage1 is not the correct length.", false);
    }
  } else {
    new_stage1_names = stage1_names;
  }
  
  if (age2_vec.length() != vec_max) {
    if (age2_vec.length() == 1) {
      IntegerVector na2v (vec_max, age2_vec(0));
      new_age2_vec = na2v;
      
    } else {
      throw Rcpp::exception("Vector age2 is not the correct length.", false);
    }
  } else {
    new_age2_vec = age2_vec;
  }
  
  if (style_vec.length() != vec_max) {
    if (style_vec.length() == 1) {
      IntegerVector ns2v (vec_max, style_vec(0));
      new_style_vec = ns2v;
      
    } else {
      throw Rcpp::exception("Vector style is not the correct length.", false);
    }
  } else {
    new_style_vec = style_vec;
  }
  
  if (time_delay_vec.length() != vec_max) {
    if (time_delay_vec.length() == 1) {
      IntegerVector ntd2v (vec_max, time_delay_vec(0));
      new_time_delay_vec = ntd2v;
      
    } else {
      throw Rcpp::exception("Vector time_delay is not the correct length.", false);
    }
  } else {
    new_time_delay_vec = time_delay_vec;
  }
  
  if (type_vec.length() != vec_max) {
    if (type_vec.length() == 1) {
      IntegerVector nt2v (vec_max, type_vec(0));
      new_type_vec = nt2v;
      
    } else {
      throw Rcpp::exception("Vector type is not the correct length.", false);
    }
  } else {
    new_type_vec = type_vec;
  }
  
  if (type_t12_vec.length() != vec_max) {
    if (type_t12_vec.length() == 1) {
      IntegerVector ntt122v (vec_max, type_t12_vec(0));
      new_type_t12_vec = ntt122v;
      
    } else {
      throw Rcpp::exception("Vector type_t12 is not the correct length.", false);
    }
  } else {
    new_type_t12_vec = type_t12_vec;
  }
  
  if (alpha_vec.length() != vec_max) {
    if (alpha_vec.length() == 1) {
      NumericVector nav (vec_max, alpha_vec(0));
      new_alpha_vec = nav;
      
    } else {
      throw Rcpp::exception("Vector alpha is not the correct length.", false);
    }
  } else {
    new_alpha_vec = alpha_vec;
  }
  
  if (beta_vec.length() != vec_max) {
    if (beta_vec.length() == 1) {
      NumericVector nbv (vec_max, beta_vec(0));
      new_beta_vec = nbv;
      
    } else {
      throw Rcpp::exception("Vector beta is not the correct length.", false);
    }
  } else {
    new_beta_vec = beta_vec;
  }
  
  DataFrame output_tab = DataFrame::create(_["stage3"] = new_stage3_names,
    _["stage2"] = new_stage2_names, _["stage1"] = new_stage1_names,
    _["age2"] = new_age2_vec, _["style"] = new_style_vec,
    _["time_delay"] = new_time_delay_vec, _["alpha"] = new_alpha_vec,
    _["beta"] = new_beta_vec, _["type"] = new_type_vec,
    _["type_t12"] = new_type_t12_vec);
  
  DataFrame output = density_reassess(stageframe, output_tab, agestages,
    historical, agebystage);
  
  StringVector needed_classes {"data.frame", "lefkoDens"};
  output.attr("class") = needed_classes;
  
  StringVector stage3_final = as<StringVector>(output["stage3"]);
  StringVector stage2_final = as<StringVector>(output["stage2"]);
  StringVector stage1_final = as<StringVector>(output["stage1"]);
  StringVector age2_final = as<StringVector>(output["age2"]);
  
  int final_rows = stage3_final.length();
  
  StringVector check_elems(final_rows);
  
  for (int i = 0; i < final_rows; i++) {
    check_elems(i) = stage3_final(i);
    check_elems(i) += " ";
    check_elems(i) += stage2_final(i);
    check_elems(i) += " ";
    if (!StringVector::is_na(stage1_final(i))) {
      check_elems(i) += stage1_final(i);
      check_elems(i) += " ";
    }
    if (!StringVector::is_na(age2_final(i))) check_elems(i) += age2_final(i);
  }
  
  StringVector unique_elems = unique(check_elems);
  if (unique_elems.length() != final_rows) {
    Rf_warningcall(R_NilValue,
      "Some transitions are listed multiple times. This may cause errors.");
  }
  
  return output;
}

//' Create a Data Frame of Supplemental Data for MPM Development
//' 
//' Function \code{supplemental()} provides all necessary supplemental data for
//' matrix estimation. It allows the establishment of proxy rates, the entry of
//' data to overwrite or offset existing rates, the identification of complete
//' reproductive transitions, and the entry of rate multipliers. The function
//' should be used to incorporate data that affects all matrices to be created.
//' To edit MPMs after creation, use \code{\link{edit_lM}()} instead.
//' 
//' @name supplemental
//' 
//' @param historical A single logical value indicating whether the MPMs
//' intended will be historical or ahistorical. Defaults to \code{TRUE}.
//' @param stagebased A single logical value indicating whether the MPM will be
//' stage-based or age-by-stage. Defaults to \code{TRUE}.
//' @param agebased A single logical value indicating whether the MPM will be
//' age-based or age-by-stage. Defaults to \code{FALSE}.
//' @param stageframe The stageframe used to produce the MPM. Required if
//' producing any stage-based or age-by-stage MPM. Must be omitted for purely
//' age-based MPMs.
//' @param stage3 String vector of stage names in occasion \emph{t}+1 in the
//' transition to be affected. Abbreviations for groups of stages are also
//' usable (see \code{Notes}). Required in all stage-based and age-by-stage
//' MPMs.
//' @param stage2 String vector of stage names in occasion \emph{t} in the
//' transition to be affected. Abbreviations for groups of stages are also
//' usable (see \code{Notes}). Required in all stage-based and age-by-stage
//' MPMs.
//' @param stage1 String vector of stage names in occasion \emph{t}-1 in the
//' transition to be affected. Only needed if a historical matrix is to be
//' produced. Abbreviations for groups of stages are also usable (see
//' \code{Notes}). Required for historical stage-based MPMs.
//' @param age2 An integer vector of the ages in occasion \emph{t} to use in
//' transitions to be affected. Required for all age- and age-by-stage MPMs.
//' @param eststage3 String vector of stage names to replace \code{stage3} in a
//' proxy transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage2 String vector of stage names to replace \code{stage2} in a
//' proxy transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage1 String vector of stage names to replace \code{stage1} in a
//' proxy historical transition. Only needed if a transition will be replaced by
//' another estimated transition, and the matrix to be estimated is historical
//' and stage-based. Stage \code{NotAlive} is also possible for raw hMPMs as a
//' means of handling the prior stage for individuals entering the population in
//' occasion \emph{t}.
//' @param estage2 Integer vector of age at time \emph{t} to replace \code{age2}
//' in a proxy transition. Only needed if a transition will be replaced by
//' another estimated transition, and only in age-based and age-by-stage MPMs.
//' @param givenrate A numeric vector of fixed rates or probabilities to replace
//' for the transition described by \code{stage3}, \code{stage2}, \code{stage1},
//' and/or \code{age2}.
//' @param offset A numeric vector of fixed numeric values to add to the
//' transitions described by \code{stage3}, \code{stage2}, \code{stage1}, and/or
//' \code{age2}.
//' @param multiplier A numeric vector of multipliers for the transition
//' described by \code{stage3}, \code{stage2}, \code{stage1}, and/or
//' \code{age2}, or for the proxy transitions described by \code{eststage3},
//' \code{eststage2}, \code{eststage1}, and/or \code{estage2}. Defaults to
//' \code{1}.
//' @param type Integer vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition;
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition; or \code{3}, \code{R}, or \code{r} for a fecundity set value /
//' general multiplier. If empty or not provided, then defaults to \code{1} for
//' survival transition.
//' @param type_t12 An optional integer vector denoting the kind of transition
//' between occasions \emph{t}-1 and \emph{t}. Only necessary if a historical
//' MPM in deVries format is desired. This should be entered as \code{1},
//' \code{S}, or \code{s} for a survival transition; or \code{2}, \code{F}, or
//' \code{f} for a fecundity transitions. Defaults to \code{1} for survival
//' transition, with impacts only on the construction of deVries-format hMPMs.
//' 
//' @return A data frame of class \code{lefkoSD}. This object can be used as
//' input in \code{\link{flefko3}()}, \code{\link{flefko2}()}, 
//' \code{\link{rlefko3}()}, \code{\link{rlefko2}()}, and 
//' \code{\link{aflefko2}()}.
//' 
//' Variables in this object include the following:
//' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be
//' replaced.}
//' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
//' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be
//' replaced.}
//' \item{age2}{Age at occasion \emph{t} in the transition to be replaced.}
//' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
//' the transition designated by \code{stage3}, \code{stage2}, and 
//' \code{stage1}.}
//' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace the
//' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
//' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
//' the transition designated by \code{stage3}, \code{stage2}, and 
//' \code{stage1}.}
//' \item{estage2}{Age at occasion \emph{t} in the transition to replace the
//' transition designated by \code{age2}.}
//' \item{givenrate}{A constant to be used as the value of the transition.}
//' \item{offset}{A constant value to be added to the transition or proxy
//' transition.}
//' \item{multiplier}{A multiplier for proxy transitions or for fecundity.}
//' \item{convtype}{Designates whether the transition from occasion \emph{t} to
//' occasion \emph{t}+1 is a survival transition probability (1), a fecundity
//' rate (2), or a fecundity multiplier (3).}
//' \item{convtype_t12}{Designates whether the transition from occasion
//' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), a
//' fecundity rate (2).}
//' 
//' @section Notes:
//' Negative values are not allowed in \code{givenrate} and \code{multiplier}
//' input, but are allowed in \code{offset}, if values are to be subtracted from
//' specific estimated transitions. Stage entries should not be used for purely
//' age-based MPMs, and age entries should not be used for purely stage-based
//' MPMs.
//' 
//' Fecundity multiplier data supplied via the \code{supplemental()} function
//' acts in the same way as non-zero entries supplied via a reproductive matrix,
//' but gets priority in all matrix creations. Thus, in cases where fecundity
//' multipliers are provided for the same function via the reproductive matrix
//' and function \code{supplemental()}, the latter is used.
//' 
//' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
//' abbreviations for groups of stages. Use \code{rep} if all reproductive
//' stages are to be used, \code{nrep} if all mature but non-reproductive stages
//' are to be used, \code{mat} if all mature stages are to be used, \code{immat}
//' if all immature stages are to be used, \code{prop} if all propagule stages
//' are to be used, \code{npr} if all non-propagule stages are to be used,
//' \code{obs} if all observable stages are to be used, \code{nobs} if all
//' unobservable stages are to be used, and leave empty or use \code{all} if all
//' stages in stageframe are to be used. Also use \code{groupX} to denote all
//' stages in group X (e.g. \code{group1} will use all stages in the respective
//' stageframe's group 1).
//' 
//' Type 3 conversions are referred to as fecundity set values, or general
//' fecundity multipliers. These set the transitions to be used as fecundity
//' transitions. Transitions set here will be interpreted as being generally
//' reproductive, meaning that the from and to stages will be used to determine
//' the general fecundity transitions to incorporate into stage-based MPMs,
//' while the age portion of the input will be used to incorporate the actual
//' multiplier(s) specified. If only stage transitions at certain ages are
//' expected to be the sole contributors to fecundity, then type 2 conversions
//' should also be included in the supplement (Type 1 and 2 conversions can be
//' purely age-specific, and do not set reproductive transitions in MPM
//' creation). For example, if all stage 2 to stage 3 transitions above age 2
//' yield fecundity, then stage 2 to stage 3 can be set to
//' \code{multiplier = 1.0} with \code{convtype = 3}, and the same transition
//' for \code{age2 = c(1, 2)} can be set to \code{multiplier = c(0, 0)}.
//' 
//' Several operations may be included per transition. Operations on the same
//' row of the resulting data frame are generally handled with given rate
//' substitutions first, then with proxy transitions, then by additive offsets,
//' and finally by multipliers. This order can be manipulated by ordering
//' operations across rows, with higher numbered rows in the data frame being
//' performed later.
//' 
//' @seealso \code{\link{edit_lM}()}
//' 
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
//'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
//'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
//'   yearcol = "year2", indivcol = "individ")
//' 
//' # Cypripedium example
//' data(cypdata)
//' 
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' 
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' @export supplemental
// [[Rcpp::export(supplemental)]]
Rcpp::List supplemental (bool historical = true, bool stagebased = true,
  bool agebased = false, Nullable<RObject> stageframe = R_NilValue,
  Nullable<RObject> stage3 = R_NilValue, Nullable<RObject> stage2 = R_NilValue,
  Nullable<RObject> stage1 = R_NilValue, Nullable<RObject> age2 = R_NilValue,
  Nullable<RObject> eststage3 = R_NilValue, Nullable<RObject> eststage2 = R_NilValue,
  Nullable<RObject> eststage1 = R_NilValue, Nullable<RObject> estage2 = R_NilValue,
  Nullable<RObject> givenrate = R_NilValue, Nullable<RObject> offset = R_NilValue,
  Nullable<RObject> multiplier = R_NilValue, Nullable<RObject> type = R_NilValue,
  Nullable<RObject> type_t12 = R_NilValue) {
  
  int wtf {-1};
  
  if (historical && stagebased && !agebased) {
    wtf = 0;
  } else if (!historical && stagebased && !agebased) {
    wtf = 1;
  } else if (!historical && stagebased && agebased) {
    wtf = 2;
  } else if (!historical && !stagebased && agebased) {
    wtf = 3;
  } else {
    throw Rcpp::exception("Unsupported MPM type.", false);
  }
  
  DataFrame stageframe_;
  int sf_yes {0};
  
  if (wtf < 3) {
    if (!stageframe.isNotNull()) {
      throw Rcpp::exception("Stageframe required for stage-based MPMs.", false);
    }
    
    if (is<DataFrame>(stageframe)) stageframe_ = as<DataFrame>(stageframe);
    StringVector sf_class = stageframe_.attr("class");
    
    String sf_error = "Please enter an object of class stageframe as input.";
    if (stageframe_.containsElementNamed("stage")) {
      sf_yes++;
    }
    if (stageframe_.containsElementNamed("min_age")) {
      sf_yes++;
    }
    if (stageframe_.containsElementNamed("max_age")) {
      sf_yes++;
    }
    if (stageframe_.containsElementNamed("group")) {
      sf_yes++;
    }
    
    for (int i = 0; i < static_cast<int>(sf_class.length()); i++) {
      if (sf_class(i) == "stageframe")  sf_yes++;
    }
    if (sf_yes < 5) throw Rcpp::exception(sf_error.get_cstring(), false);
  }
  
  StringVector stage3_;
  StringVector stage2_;
  StringVector stage1_;
  IntegerVector age2_;
  StringVector eststage3_;
  StringVector eststage2_;
  StringVector eststage1_;
  IntegerVector estage2_;
  NumericVector givenrate_;
  NumericVector offset_;
  NumericVector multiplier_;
  IntegerVector type_;
  IntegerVector type_t12_;
  
  int stage3_length {0};
  int stage2_length {0};
  int stage1_length {0};
  int age2_length {0};
  int eststage3_length {0};
  int eststage2_length {0};
  int eststage1_length {0};
  int estage2_length {0};
  int type_length {0};
  int type_t12_length {0};
  
  StringVector all_stages;
  StringVector wildcard_names = {"all", "rep", "nrep", "mat", "immat", "prop",
    "npr", "notalive", "obs", "nobs"};
  StringVector all_groups;
  
  if (wtf < 3) {
    all_stages = as<StringVector>(stageframe_["stage"]);
    
    StringVector all_groups_ = as<StringVector>(stageframe_["group"]);
    StringVector all_groups_temp (all_groups_.length());
    for (int i = 0; i < all_groups_.length(); i++) {
      String group_slot = "group";
      String group_term_added = all_groups_(i);
      group_slot += group_term_added;
      all_groups_temp(i) = group_slot;
    }
    all_groups = all_groups_temp;
  }
  
  if (stage3.isNotNull() && wtf != 3) {
    if (is<StringVector>(stage3)) {
      stage3_ = as<StringVector>(stage3);
      stage3_length = static_cast<int>(stage3_.length());
      
      for (int i = 0; i < stage3_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage3_(i)) || stage3_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage3_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage3_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage3_(i) == all_groups(j)) found_stage = true;
          }
          if (stage3_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.", false);
          }
        }
        if (!found_stage) {
          String ems = "Some entered stage names (stage3) do not ";
          String ems1 = "match expectation.";
          ems += ems1;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
      }
      
    } else {
      throw Rcpp::exception("Please enter argument stage3 as text.", 
        false);
    }
  } else if (wtf != 3) {
    throw Rcpp::exception("Stage information (stage3) for transitions is required.",
      false);
  }
  
  if (stage2.isNotNull() && wtf != 3) {
    if (is<StringVector>(stage2)) {
      stage2_ = as<StringVector>(stage2);
      stage2_length = static_cast<int>(stage2_.length());
      
      for (int i = 0; i < stage2_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage2_(i)) || stage2_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage2_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage2_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage2_(i) == all_groups(j)) found_stage = true;
          }
          if (stage2_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" can only be used in eststage1.",
              false);
          }
        }
        if (!found_stage) {
          String ems = "Some entered stage names (stage2) do not ";
          String ems1 = "match expectation.";
          ems += ems1;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
      }
      
    } else {
      throw Rcpp::exception("Please enter argument stage2 as text.",
        false);
    }
  } else if (wtf != 3) {
    throw Rcpp::exception("Stage information (stage2) for transitions is required.",
      false);
  }
  
  if (stage1.isNotNull() && wtf == 0) {
    if (is<StringVector>(stage1)) {
      stage1_ = as<StringVector>(stage1);
      stage1_length = static_cast<int>(stage1_.length());
      
      for (int i = 0; i < stage1_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage1_(i)) || stage1_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage1_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage1_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage1_(i) == all_groups(j)) found_stage = true;
          }
          if (stage1_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.", false);
          }
        }
        if (!found_stage) {
          String ems = "Some entered stage names (stage1) do not ";
          String ems1 = "match expectation.";
          ems += ems1;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
      }
      
    } else if (is<LogicalVector>(stage1)) {
      if (wtf == 0) {
        throw Rcpp::exception("Stage information (stage1) for transitions is required.",
          false);
      }
      
      StringVector stage1_temp (stage2_length, NA_STRING);
      stage1_ = stage1_temp;
      stage1_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter argument stage1 as text.", 
        false);
    }
  } else {
    if (wtf == 0) {
      throw Rcpp::exception("Stage information (stage1) for transitions is required.",
        false);
    } else if (wtf != 3) {
      StringVector stage1_temp (stage2_length, NA_STRING);
      stage1_ = stage1_temp;
      stage1_length = stage2_length;
    }
  }
  
  if (age2.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(age2) || is<NumericVector>(age2)) {
      age2_ = as<IntegerVector>(age2);
      age2_length = static_cast<int>(age2_.length());
      
      arma::ivec a2_arma = as<arma::ivec>(age2_);
      arma::uvec neg_tester = find(a2_arma < 0 && a2_arma > -2147483648);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered age2 values are negative.");
      }
    } else if (is<LogicalVector>(age2)) { 
      IntegerVector age2_temp (stage2_length, NA_INTEGER);
      age2_ = age2_temp;
      age2_length = stage2_length;
    } else {
      throw Rcpp::exception("Please enter argument age2 in integer format.",
        false);
    }
  } else {
    IntegerVector age2_temp (stage2_length, NA_INTEGER);
    age2_ = age2_temp;
    age2_length = stage2_length;
  }
  
  if (eststage3.isNotNull() && wtf != 3) {
    if (is<StringVector>(eststage3)) {
      eststage3_ = as<StringVector>(eststage3);
      eststage3_length = static_cast<int>(eststage3_.length());
      
      for (int i = 0; i < eststage3_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage3_(i)) || eststage3_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage3_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage3_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage3_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage3_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed in argument eststage1.", false);
          }
        } else found_stage = true;
        
        if (!found_stage) {
          String ems = "Some entered stage names (eststage3) do not ";
          String ems1 = "match expectation.";
          ems += ems1;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
      }
      
    } else if (is<LogicalVector>(eststage3)) {
      StringVector eststage3_temp (stage3_length, NA_STRING);
      eststage3_ = eststage3_temp;
      eststage3_length = stage3_length;
      
    } else {
      throw Rcpp::exception("Please enter argument eststage3 as text.", false);
    }
  } else if (stage3_length != 0) {
    StringVector eststage3_temp (stage3_length, NA_STRING);
    eststage3_ = eststage3_temp;
    eststage3_length = stage3_length;
  }
  
  if (eststage2.isNotNull() && wtf != 3) {
    if (is<StringVector>(eststage2)) {
      eststage2_ = as<StringVector>(eststage2);
      eststage2_length = static_cast<int>(eststage2_.length());
      
      for (int i = 0; i < eststage2_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage2_(i)) || eststage2_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage2_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage2_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage2_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage2_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed in argument eststage1.", false);
          }
        } else found_stage = true;
        
        if (!found_stage) {
          String ems = "Some entered stage names (eststage2) do not ";
          String ems1 = "match expectation.";
          ems += ems1;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
      }
      
    } else if (is<LogicalVector>(eststage2)) {
      StringVector eststage2_temp (stage2_length, NA_STRING);
      eststage2_ = eststage2_temp;
      eststage2_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter argument eststage2 as text.", false);
    }
  } else if (stage2_length != 0) {
    StringVector eststage2_temp (stage2_length, NA_STRING);
    eststage2_ = eststage2_temp;
    eststage2_length = stage2_length;
  }
  
  if (eststage1.isNotNull() && wtf == 0) {
    if (is<StringVector>(eststage1)) {
      eststage1_ = as<StringVector>(eststage1);
      eststage1_length = static_cast<int>(eststage1_.length());
      
      for (int i = 0; i < eststage1_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage1_(i)) || eststage1_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage1_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage1_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage1_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage1_(i) == "NotAlive") found_stage = true;
        } else found_stage = true;
        
        if (!found_stage) {
          String ems = "Some entered stage names (eststage1) do not ";
          String ems1 = "match expectation.";
          ems += ems1;
          
          throw Rcpp::exception(ems.get_cstring(), false);
        }
      }
      
    } else if (is<LogicalVector>(eststage1)) {
      StringVector eststage1_temp (stage2_length, NA_STRING);
      eststage1_ = eststage1_temp;
      eststage1_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter argument eststage1 as text.", false);
    }
  } else if (stage2_length != 0) {
    StringVector eststage1_temp (stage2_length, NA_STRING);
    eststage1_ = eststage1_temp;
    eststage1_length = stage2_length;
  }
  
  if (estage2.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(estage2) || is<NumericVector>(estage2)) {
      estage2_ = as<IntegerVector>(estage2);
      estage2_length = static_cast<int>(estage2_.length());
      
      arma::ivec ea2_arma = as<arma::ivec>(estage2_);
      arma::uvec neg_tester = find(ea2_arma < 0 && ea2_arma > -2147483648);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered estage2 values are negative.");
      }
    } else if (is<LogicalVector>(estage2)) {
      LogicalVector estage2_temp = as<LogicalVector>(estage2);
      estage2_length = static_cast<int>(estage2_temp.length());
      
      IntegerVector estage2_new (estage2_length, NA_INTEGER);
      estage2_ = estage2_new;
      
    } else {
      throw Rcpp::exception("Please enter argument estage2 in integer format.",
        false);
    }
  } else {
    if (stage2_length > 0) {
      IntegerVector estage2_temp (stage2_length, NA_INTEGER);
      estage2_ = estage2_temp;
      estage2_length = stage2_length;
    } else if (age2_length > 0) {
      IntegerVector estage2_temp (age2_length, NA_INTEGER);
      estage2_ = estage2_temp;
      estage2_length = age2_length;
    }
  }
  
  LefkoInputs::numeric_vectorizer(givenrate_, givenrate, "givenrate",
    stage2_length, age2_length, false, 0.);
  LefkoInputs::numeric_vectorizer(offset_, offset, "offset", stage2_length,
    age2_length, false, 0.);
  LefkoInputs::numeric_vectorizer(multiplier_, multiplier, "multiplier",
    stage2_length, age2_length, true, 1.);
  
  IntegerVector type_int_limits = {1, 2, 3};
  CharacterVector type_char_limits = {"S", "F", "R"};
  IntegerVector type12_int_limits = {1, 2};
  CharacterVector type12_char_limits = {"S", "F"};
  
  LefkoInputs::integer_char_vectorizer (type_, type, "type", stage2_length, age2_length,
    type_int_limits, type_char_limits, true, 1);
  LefkoInputs::integer_char_vectorizer (type_t12_, type_t12, "type_t12", stage2_length,
    age2_length, type12_int_limits, type12_char_limits, true, 1);
  
  type_t12_length = static_cast<int>(type_t12_.length());
  
  if (wtf < 3 && type_length != 0  && type_length != stage2_length) { 
    throw Rcpp::exception("All input vectors must be of the same length.", false);
  }
  
  if (wtf == 0) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    
    if (stage2_length != stage3_length || stage2_length != stage1_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage1_length != 0 && stage1_length != eststage1_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage2_length != 0 && stage2_length != eststage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage3_length != 0 && stage3_length != eststage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
  } else if (wtf == 1) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    
    if (stage2_length != stage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage2_length != 0 && stage2_length != eststage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage3_length != 0 && stage3_length != eststage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
  } else if (wtf == 2) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    if (age2_length == 0) throw Rcpp::exception("Age (age2) information required.", false);
    
    if (stage2_length != stage3_length || stage2_length != age2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (estage2_length != 0 && age2_length != estage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  
  } else if (wtf == 3) {
    if (age2_length == 0) throw Rcpp::exception("Age (age2) information required.", false);
    
    if (estage2_length != 0 && age2_length != estage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  }
  
  if ((type_t12_length == 0) & historical) {
    if (age2_length != 0) {
      IntegerVector cvt12_temp (age2_length, 1);
      type_t12_ = cvt12_temp;
    } else if (stage2_length != 0) {
      IntegerVector cvt12_temp (stage2_length, 1);
      type_t12_ = cvt12_temp;
    }
  }
  
  if (stage3_length == 0) { 
    StringVector s3_temp (age2_length, NA_STRING);
    stage3_ = s3_temp;
    stage3_length = age2_length;
  }
  if (stage2_length == 0) { 
    StringVector s2_temp (age2_length, NA_STRING);
    stage2_ = s2_temp;
    stage2_length = age2_length;
  }
  if (stage1_length == 0) { 
    StringVector s1_temp (stage2_length, NA_STRING);
    stage1_ = s1_temp;
    stage1_length = stage2_length;
  }
  if (eststage3_length == 0) { 
    StringVector s3_temp (stage2_length, NA_STRING);
    eststage3_ = s3_temp;
    eststage3_length = stage2_length;
  }
  if (eststage2_length == 0) { 
    StringVector s2_temp (stage2_length, NA_STRING);
    eststage2_ = s2_temp;
    eststage2_length = stage2_length;
  }
  if (eststage1_length == 0) { 
    StringVector s1_temp (stage2_length, NA_STRING);
    eststage1_ = s1_temp;
    eststage1_length = stage2_length;
  }
  
  List supplement (13);
  supplement(0) = stage3_;
  supplement(1) = stage2_;
  supplement(2) = stage1_;
  supplement(3) = age2_;
  supplement(4) = eststage3_;
  supplement(5) = eststage2_;
  supplement(6) = eststage1_;
  supplement(7) = estage2_;
  supplement(8) = givenrate_;
  supplement(9) = offset_;
  supplement(10) = multiplier_;
  supplement(11) = type_;
  supplement(12) = type_t12_;
  
  StringVector supp_names = {"stage3", "stage2", "stage1", "age2", "eststage3",
    "eststage2", "eststage1", "estage2", "givenrate", "offset", "multiplier",
    "convtype", "convtype_t12"};
  StringVector supp_class = {"data.frame", "lefkoSD"};
  
  supplement.attr("class") = supp_class;
  supplement.attr("names") = supp_names;
  supplement.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stage2_length);
  
  return supplement;
}

//' Edit Single lefkoMat Object based on Supplemental Data
//' 
//' Function \code{edit_lM_single()} edits a single existing \code{lefkoMat}
//' objects with external data supplied by the user. This is the workhorse
//' function for function \code{edit_lM}.
//' 
//' @name edit_lM_single
//' 
//' @param final_output The final list structure to be modified.
//' @param mpm The \code{lefkoMat} object to be edited.
//' @param pop A string vector denoting the populations to be edited. Defaults
//' to \code{NULL}, in which case all populations are edited.
//' @param patch A string vector denoting the patches to be edited. Defaults
//' to \code{NULL}, in which case all patches are edited.
//' @param year2 A string vector denoting the years to be edited. Defaults
//' to \code{NULL}, in which case all years are edited.
//' @param stage3 The name of the stage in occasion \emph{t}+1 in the transition
//' to be replaced. Abbreviations for groups of stages are also usable (see
//' \code{Notes}). Required in all stage-based and age-by-stage MPMs.
//' @param stage2 The name of the stage in occasion \emph{t} in the transition
//' to be replaced. Abbreviations for groups of stages are also usable (see
//' \code{Notes}). Required in all stage-based and age-by-stage MPMs.
//' @param stage1 The name of the stage in occasion \emph{t}-1 in the transition
//' to be replaced. Only needed if a historical matrix is to be produced.
//' Abbreviations for groups of stages are also usable (see \code{Notes}).
//' Required for historical stage-based MPMs.
//' @param age2 An integer vector of the ages in occasion \emph{t} to use in
//' transitions to be changed or replaced. Required for all age- and
//' age-by-stage MPMs.
//' @param eststage3 The name of the stage to replace \code{stage3} in a proxy
//' transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage2 The name of the stage to replace \code{stage2} in a proxy
//' transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage1 The name of the stage to replace \code{stage1} in a proxy
//' historical transition. Only needed if a transition will be replaced by
//' another estimated transition, and the matrix to be estimated is historical
//' and stage-based. Stage \code{NotAlive} is also possible for raw hMPMs as a
//' means of handling the prior stage for individuals entering the population in
//' occasion \emph{t}.
//' @param estage2 The age at time \emph{t} to replace \code{age2} in a proxy
//' transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in age-based and age-by-stage MPMs.
//' @param givenrate A fixed rate or probability to replace for the transition
//' described by \code{stage3}, \code{stage2}, and \code{stage1}.
//' @param offset A numeric vector of fixed numeric values to add to the
//' transitions described by \code{stage3}, \code{stage2}, \code{stage1}, and/or
//' \code{age2}.
//' @param multiplier A vector of numeric multipliers for fecundity or for proxy
//' transitions. Defaults to \code{1}.
//' @param type A vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition;
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition; or \code{3}, \code{R}, or \code{r} for a fecundity multiplier.
//' If empty or not provided, then defaults to \code{1} for survival transition.
//' @param type_t12 An optional vector denoting the kind of transition between
//' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
//' deVries format is desired. This should be entered as \code{1}, \code{S}, or
//' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
//' fecundity transitions. Defaults to \code{1} for survival transition, with
//' impacts only on the construction of deVries-format hMPMs.
//' @param target_mpm If modifying a \code{lefkoMatList} object, then this
//' allows the user to specify which MPMs to modify. To modify, enter a vector
//' with the number of each MPM to modify, or enter \code{"all"} to modify all
//' MPMs. Defaults to \code{"all"}.
//' 
//' @return This creates a modified \code{lefkoMat} object, and places it at
//' the reference for object \code{final_output}.
//' 
//' @keywords internal
//' @noRd
void edit_lM_single (List& final_output, const List& mpm, Nullable<RObject> pop = R_NilValue,
  Nullable<RObject> patch = R_NilValue, Nullable<RObject> year2 = R_NilValue,
  Nullable<RObject> stage3 = R_NilValue, Nullable<RObject> stage2 = R_NilValue,
  Nullable<RObject> stage1 = R_NilValue, Nullable<RObject> age2 = R_NilValue,
  Nullable<RObject> eststage3 = R_NilValue, Nullable<RObject> eststage2 = R_NilValue,
  Nullable<RObject> eststage1 = R_NilValue, Nullable<RObject> estage2 = R_NilValue,
  Nullable<RObject> givenrate = R_NilValue, Nullable<RObject> offset = R_NilValue,
  Nullable<RObject> multiplier = R_NilValue, Nullable<RObject> type = R_NilValue,
  Nullable<RObject> type_t12 = R_NilValue) {
  
  String mpm_error = "Please enter a lefkoMat or lefkoMatList object as input.";
  if (!mpm.containsElementNamed("ahstages")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  if (!mpm.containsElementNamed("hstages")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  if (!mpm.containsElementNamed("agestages")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  if (!mpm.containsElementNamed("labels")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  Rcpp::DataFrame ahstages = as<DataFrame>(mpm["ahstages"]);
  Rcpp::DataFrame hstages = as<DataFrame>(mpm["hstages"]);
  Rcpp::DataFrame agestages = as<DataFrame>(mpm["agestages"]);
  Rcpp::DataFrame labels = as<DataFrame>(mpm["labels"]);
  
  int wtf = LefkoUtils::whichbrew(ahstages, hstages, agestages);
  
  StringVector pop_;
  StringVector patch_;
  StringVector year2_;
  StringVector stage3_;
  StringVector stage2_;
  StringVector stage1_;
  IntegerVector age2_;
  StringVector eststage3_;
  StringVector eststage2_;
  StringVector eststage1_;
  IntegerVector estage2_;
  NumericVector givenrate_;
  NumericVector offset_;
  NumericVector multiplier_;
  IntegerVector type_;
  IntegerVector type_t12_;
  
  int stage3_length {0};
  int stage2_length {0};
  int stage1_length {0};
  int age2_length {0};
  int eststage3_length {0};
  int eststage2_length {0};
  int eststage1_length {0};
  int estage2_length {0};
  int type_length {0};
  int type_t12_length {0};
  int min_age {0};
  int max_age {0};
  
  StringVector all_stages = as<StringVector>(ahstages["stage"]);
  StringVector wildcard_names = {"all", "rep", "nrep", "mat", "immat", "prop",
    "npr", "notalive", "obs", "nobs"};
  
  StringVector all_groups_ = as<StringVector>(ahstages["group"]);
  StringVector all_groups (all_groups_.length());
  for (int i = 0; i < all_groups_.length(); i++) {
    String group_slot = "group";
    String group_term_added = all_groups_(i);
    group_slot += group_term_added;
    all_groups(i) = group_slot;
  }
  
  if (stage3.isNotNull() && wtf != 3) {
    if (is<StringVector>(stage3)) {
      stage3_ = as<StringVector>(stage3);
      stage3_length = static_cast<int>(stage3_.length());
      
      for (int i = 0; i < stage3_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage3_(i)) || stage3_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage3_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage3_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage3_(i) == all_groups(j)) found_stage = true;
          }
          if (stage3_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.", false);
          }
        }
        if (!found_stage) {
          throw Rcpp::exception("Some stage names (stage3) do not match expectation.",
            false);
        }
      }
      
    } else {
      throw Rcpp::exception("Please enter stage information (stage3) as text.", 
        false);
    }
  } else if (wtf != 3) {
    throw Rcpp::exception("Stage information (stage3) for transitions is required.",
      false);
  }
  
  if (stage2.isNotNull() && wtf != 3) {
    if (is<StringVector>(stage2)) {
      stage2_ = as<StringVector>(stage2);
      stage2_length = static_cast<int>(stage2_.length());
      
      for (int i = 0; i < stage2_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage2_(i)) || stage2_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage2_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage2_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage2_(i) == all_groups(j)) found_stage = true;
          }
          if (stage2_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.", false);
          }
        }
        if (!found_stage) {
          throw Rcpp::exception("Some stage names (stage2) do not match expectation.",
            false);
        }
      }
      
    } else {
      throw Rcpp::exception("Please enter stage information (stage2) as text.",
        false);
    }
  } else if (wtf != 3) {
    throw Rcpp::exception("Stage information (stage2) for transitions is required.",
      false);
  }
  
  if (stage1.isNotNull() && wtf == 0) {
    if (is<StringVector>(stage1)) {
      stage1_ = as<StringVector>(stage1);
      stage1_length = static_cast<int>(stage1_.length());
      
      for (int i = 0; i < stage1_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(stage1_(i)) || stage1_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (stage1_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (stage1_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (stage1_(i) == all_groups(j)) found_stage = true;
          }
          if (stage1_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.", false);
          }
        }
        if (!found_stage) {
          throw Rcpp::exception("Some stage names (stage1) do not match expectation.", false);
        }
      }
      
    } else if (is<LogicalVector>(stage1)) {
      if (wtf == 0) {
        throw Rcpp::exception("Stage information (stage1) for transitions is required.",
          false);
      }
      
      StringVector stage1_temp (stage2_length, NA_STRING);
      stage1_ = stage1_temp;
      stage1_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter stage information (stage1) as text.", 
        false);
    }
  } else {
    if (wtf == 0) {
      throw Rcpp::exception("Stage information (stage1) for transitions is required.",
        false);
    } else if (wtf != 3) {
      StringVector stage1_temp (stage2_length, NA_STRING);
      stage1_ = stage1_temp;
      stage1_length = stage2_length;
    }
  }
  
  if (age2.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(age2) || is<NumericVector>(age2)) {
      age2_ = as<IntegerVector>(age2);
      age2_length = static_cast<int>(age2_.length());
      
      arma::ivec a2_arma = as<arma::ivec>(age2_);
      arma::uvec neg_tester = find(a2_arma < 0);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered age2 values are negative.");
      }
    } else {
      throw Rcpp::exception("Please enter age information (age2) in integer format.",
        false);
    }
  } else {
    if (wtf > 1) {
      throw Rcpp::exception("Age information (age2) for transitions is required.",
        false);
    } else if (wtf != 3) {
      IntegerVector age2_temp (stage2_length, NA_INTEGER);
      age2_ = age2_temp;
      age2_length = stage2_length;
    }
  }
  
  if (pop.isNotNull()) {
    pop_ = as<StringVector>(pop);
  } else {
    if (wtf < 3) {
      StringVector pop_temp (stage2_length, NA_STRING);
      pop_ = pop_temp;
    } else {
      StringVector pop_temp (age2_length, NA_STRING);
      pop_ = pop_temp;
    }
  }
  
  if (patch.isNotNull()) {
    patch_ = as<StringVector>(patch);
  } else {
    if (wtf < 3) {
      StringVector patch_temp (stage2_length, NA_STRING);
      patch_ = patch_temp;
    } else {
      StringVector patch_temp (age2_length, NA_STRING);
      patch_ = patch_temp;
    }
  }
  
  if (year2.isNotNull()) {
    year2_ = as<StringVector>(year2);
  } else {
    if (wtf < 3) {
      StringVector year2_temp (stage2_length, NA_STRING);
      year2_ = year2_temp;
    } else {
      StringVector year2_temp (age2_length, NA_STRING);
      year2_ = year2_temp;
    }
  }
  
  if (eststage3.isNotNull() && wtf != 3) {
    if (is<StringVector>(eststage3)) {
      eststage3_ = as<StringVector>(eststage3);
      eststage3_length = static_cast<int>(eststage3_.length());
      
      for (int i = 0; i < eststage3_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage3_(i)) || eststage3_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage3_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage3_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage3_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage3_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.",
              false);
          }
        } else found_stage = true;
        
        if (!found_stage) {
          throw Rcpp::exception("Some stage names (eststage3) do not match expectation.",
            false);
        }
      }
      
    } else if (is<LogicalVector>(eststage3)) {
      StringVector eststage3_temp (stage3_length, NA_STRING);
      eststage3_ = eststage3_temp;
      eststage3_length = stage3_length;
      
    } else {
      throw Rcpp::exception("Please enter stage information (eststage3) as text.",
        false);
    }
  } else if (stage3_length != 0) {
    StringVector eststage3_temp (stage3_length, NA_STRING);
    eststage3_ = eststage3_temp;
    eststage3_length = stage3_length;
  }
  
  if (eststage2.isNotNull() && wtf != 3) {
    if (is<StringVector>(eststage2)) {
      eststage2_ = as<StringVector>(eststage2);
      eststage2_length = static_cast<int>(eststage2_.length());
      
      for (int i = 0; i < eststage2_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage2_(i)) || eststage2_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage2_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage2_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage2_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage2_(i) == "NotAlive") {
            throw Rcpp::exception("Stage \"NotAlive\" is only allowed as input in eststage1.",
              false);
          }
        } else found_stage = true;
        
        if (!found_stage) {
          throw Rcpp::exception("Some stage names (eststage2) do not match expectation.", false);
        }
      }
      
    } else if (is<LogicalVector>(eststage2)) {
      StringVector eststage2_temp (stage2_length, NA_STRING);
      eststage2_ = eststage2_temp;
      eststage2_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter stage information (eststage2) as text.",
        false);
    }
  } else if (stage2_length != 0) {
    StringVector eststage2_temp (stage2_length, NA_STRING);
    eststage2_ = eststage2_temp;
    eststage2_length = stage2_length;
  }
  
  if (eststage1.isNotNull() && wtf == 0) {
    if (is<StringVector>(eststage1)) {
      eststage1_ = as<StringVector>(eststage1);
      eststage1_length = static_cast<int>(eststage1_.length());
      
      for (int i = 0; i < eststage1_length; i++) {
        bool found_stage {false};
        
        if (!(StringVector::is_na(eststage1_(i)) || eststage1_(i) == "NA")) { 
          for (int j = 0; j < all_stages.length(); j++) {
            if (eststage1_(i) == all_stages(j)) found_stage = true;
          }
          for (int j = 0; j < wildcard_names.length(); j++) {
            if (eststage1_(i) == wildcard_names(j)) found_stage = true;
          }
          for (int j = 0; j < all_groups.length(); j++) {
            if (eststage1_(i) == all_groups(j)) found_stage = true;
          }
          if (eststage1_(i) == "NotAlive") found_stage = true;
        } else found_stage = true;
        
        if (!found_stage) {
          throw Rcpp::exception("Some stage names (eststage1) do not match expectation.",
            false);
        }
      }
      
    } else if (is<LogicalVector>(eststage1)) {
      StringVector eststage1_temp (stage2_length, NA_STRING);
      eststage1_ = eststage1_temp;
      eststage1_length = stage2_length;
      
    } else {
      throw Rcpp::exception("Please enter stage information (eststage1) as text.",
        false);
    }
  } else if (stage2_length != 0) {
    StringVector eststage1_temp (stage2_length, NA_STRING);
    eststage1_ = eststage1_temp;
    eststage1_length = stage2_length;
  }
  
  if (estage2.isNotNull() && wtf > 1) {
    if (is<IntegerVector>(estage2) || is<NumericVector>(estage2)) {
      estage2_ = as<IntegerVector>(estage2);
      estage2_length = static_cast<int>(estage2_.length());
      
      arma::ivec ea2_arma = as<arma::ivec>(estage2_);
      arma::uvec neg_tester = find(ea2_arma < 0);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered estage2 values are negative.");
      }
    } else {
      throw Rcpp::exception("Please enter age information (estage2) in integer format.",
        false);
    }
  } else {
    if (stage2_length > 0) {
      IntegerVector estage2_temp (stage2_length, NA_INTEGER);
      estage2_ = estage2_temp;
      estage2_length = stage2_length;
    } else if (age2_length > 0) {
      IntegerVector estage2_temp (age2_length, NA_INTEGER);
      estage2_ = estage2_temp;
      estage2_length = age2_length;
    }
  }
  
  if (givenrate.isNotNull()) {
    if (is<NumericVector>(givenrate)) {
      givenrate_ = as<NumericVector>(givenrate);
      
      arma::vec gvr_arma = as<arma::vec>(givenrate_);
      arma::uvec neg_tester = find(gvr_arma < 0.0);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered given rate values are negative.");
      }
    } else if (is<LogicalVector>(givenrate)) {
      if (age2_length != 0) {
        NumericVector givenrate_temp (age2_length, NA_REAL);
        givenrate_ = givenrate_temp;
      } else if (stage2_length != 0) {
        NumericVector givenrate_temp (stage2_length, NA_REAL);
        givenrate_ = givenrate_temp;
      }
      
    } else {
      throw Rcpp::exception("Please enter given rate information (givenrate) in numeric format.",
        false);
    }
  } else {
    if (age2_length != 0) {
      NumericVector givenrate_temp (age2_length, NA_REAL);
      givenrate_ = givenrate_temp;
    } else if (stage2_length != 0) {
      NumericVector givenrate_temp (stage2_length, NA_REAL);
      givenrate_ = givenrate_temp;
    }
  }
  
  if (offset.isNotNull()) {
    if (is<NumericVector>(offset)) {
      offset_ = as<NumericVector>(offset);
    } else if (is<LogicalVector>(offset)) {
      if (age2_length != 0) {
        NumericVector offset_temp (age2_length, NA_REAL);
        offset_ = offset_temp;
      } else if (stage2_length != 0) {
        NumericVector offset_temp (stage2_length, NA_REAL);
        offset_ = offset_temp;
      }
    } else {
      throw Rcpp::exception("Please enter offset information (offset) in numeric format.",
        false);
    }
  } else {
    if (age2_length != 0) {
      NumericVector offset_temp (age2_length, NA_REAL);
      offset_ = offset_temp;
    } else if (stage2_length != 0) {
      NumericVector offset_temp (stage2_length, NA_REAL);
      offset_ = offset_temp;
    }
  }
  
  if (multiplier.isNotNull()) {
    if (is<NumericVector>(multiplier)) {
      multiplier_ = as<NumericVector>(multiplier);
      
      arma::vec mpl_arma = as<arma::vec>(multiplier_);
      arma::uvec neg_tester = find(mpl_arma < 0.0);
      if (neg_tester.n_elem > 0) {
        Rf_warningcall(R_NilValue, "Some entered multiplier values are negative.");
      }
    } else if (is<LogicalVector>(multiplier)) {
      if (age2_length != 0) {
        NumericVector multiplier_temp (age2_length, 1.0);
        multiplier_ = multiplier_temp;
      } else if (stage2_length != 0) {
        NumericVector multiplier_temp (stage2_length, 1.0);
        multiplier_ = multiplier_temp;
      }
    } else {
      throw Rcpp::exception("Please enter multiplier information (multiplier) in numeric format.",
        false);
    }
  } else {
    if (age2_length != 0) {
      NumericVector multiplier_temp (age2_length, 1.0);
      multiplier_ = multiplier_temp;
    } else if (stage2_length != 0) {
      NumericVector multiplier_temp (stage2_length, 1.0);
      multiplier_ = multiplier_temp;
    }
  }
  
  if (type.isNotNull()) {
    if (is<IntegerVector>(type) || is<NumericVector>(type)) {
      type_ = as<IntegerVector>(type);
      type_length = static_cast<int>(type_.length());
      
      int check_type_min = min(type_);
      int check_type_max = max(type_);
      if (check_type_min < 1 || check_type_max > 3) {
        throw Rcpp::exception("Enter transition type (type) using only integers 1, 2, and 3.",
          false);
      }
    } else if (is<StringVector>(type)) {
      StringVector type_sv = as<StringVector>(type);
      type_length = static_cast<int>(type_sv.length());
      
      IntegerVector type_temp (type_length);
      
      for (int i = 0; i < type_length; i++) {
        if (stringcompare_simple(as<std::string>(type_sv(i)), "s", true)) {
          type_temp(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_sv(i)), "f", true)) {
          type_temp(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_sv(i)), "r", true)) {
          type_temp(i) = 3;
        } else {
          throw Rcpp::exception("Please enter transition type information (type) using only integers 1, 2, and 3.",
            false);
        }
        type_ = type_temp;
      }
    } else {
      throw Rcpp::exception("Please enter transition type information (type) in integer format.",
        false);
    }
  } else {
    throw Rcpp::exception("Information on type of transition (type) is required.",
      false);
  }
  
  if (type_t12.isNotNull()) {
    if (is<IntegerVector>(type_t12) || is<NumericVector>(type_t12)) {
      type_t12_ = as<IntegerVector>(type_t12);
      type_t12_length = static_cast<int>(type_t12_.length());
      
      int check_type_min = min(type_t12_);
      int check_type_max = max(type_t12_);
      if (check_type_min < 1 || check_type_max > 3) {
        throw Rcpp::exception("Enter transition type (type_t12) using only integers 1, 2, and 3.",
          false);
      }
    } else if (is<StringVector>(type_t12)) {
      StringVector type_t12_sv = as<StringVector>(type_t12);
      type_t12_length = static_cast<int>(type_t12_sv.length());
      
      IntegerVector type_t12_temp (type_t12_length);
      
      for (int i = 0; i < type_t12_length; i++) {
        if (stringcompare_simple(as<std::string>(type_t12_sv(i)), "s", true)) {
          type_t12_temp(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_sv(i)), "f", true)) {
          type_t12_temp(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_t12_sv(i)), "r", true)) {
          type_t12_temp(i) = 3;
        } else {
          throw Rcpp::exception("Enter historical transition type (type_t12) using only integers 1, 2, and 3.",
            false);
        }
        type_t12_ = type_t12_temp;
      }
    } else {
      throw Rcpp::exception("Enter transition type (type_t12) in integer format.",
        false);
    }
  }
  
  if (wtf < 3 && type_length != 0  && type_length != stage2_length) { 
    throw Rcpp::exception("All input vectors must be of the same length.", false);
  }
  
  if (wtf == 0) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    
    if (stage2_length != stage3_length || stage2_length != stage1_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage1_length != 0 && stage1_length != eststage1_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage2_length != 0 && stage2_length != eststage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage3_length != 0 && stage3_length != eststage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  } else if (wtf == 1) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    
    if (stage2_length != stage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage2_length != 0 && stage2_length != eststage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (eststage3_length != 0 && stage3_length != eststage3_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  } else if (wtf == 2) {
    if (stage2_length == 0) throw Rcpp::exception("Stage (stage2) information required.", false);
    if (age2_length == 0) throw Rcpp::exception("Age (age2) information required.", false);
    
    if (stage2_length != stage3_length || stage2_length != age2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
    
    if (estage2_length != 0 && age2_length != estage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  } else if (wtf == 3) {
    if (age2_length == 0) throw Rcpp::exception("Age (age2) information required.", false);
    
    if (estage2_length != 0 && age2_length != estage2_length) {
      throw Rcpp::exception("All input vectors must be of the same length.", false);
    }
  }
  
  if (type_t12_length == 0) {
    if (age2_length != 0) {
      IntegerVector cvt12_temp (age2_length, 1);
      type_t12_ = cvt12_temp;
    } else if (stage2_length != 0) {
      IntegerVector cvt12_temp (stage2_length, 1);
      type_t12_ = cvt12_temp;
    }
  }
  
  if (stage3_length == 0) { 
    StringVector s3_temp (age2_length, NA_STRING);
    stage3_ = s3_temp;
    stage3_length = age2_length;
  }
  if (stage2_length == 0) { 
    StringVector s2_temp (age2_length, NA_STRING);
    stage2_ = s2_temp;
    stage2_length = age2_length;
  }
  if (stage1_length == 0) { 
    StringVector s1_temp (stage2_length, NA_STRING);
    stage1_ = s1_temp;
    stage1_length = stage2_length;
  }
  if (eststage3_length == 0) { 
    StringVector s3_temp (stage2_length, NA_STRING);
    eststage3_ = s3_temp;
    eststage3_length = stage2_length;
  }
  if (eststage2_length == 0) { 
    StringVector s2_temp (stage2_length, NA_STRING);
    eststage2_ = s2_temp;
    eststage2_length = stage2_length;
  }
  if (eststage1_length == 0) { 
    StringVector s1_temp (stage2_length, NA_STRING);
    eststage1_ = s1_temp;
    eststage1_length = stage2_length;
  }
  
  DataFrame supplement = DataFrame::create(_["stage3"] = stage3_,
    _["stage2"] = stage2_, _["stage1"] = stage1_, _["age2"] = age2_,
    _["eststage3"] = eststage3_, _["eststage2"] = eststage2_,
    _["eststage1"] = eststage1_, _["estage2"] = estage2_,
    _["givenrate"] = givenrate_, _["offset"] = offset_,
    _["multiplier"] = multiplier_, _["convtype"] = type_ ,
    _["convtype_t12"] = type_t12_, _["pop"] = pop_, _["patch"] = patch_,
    _["year2"] = year2_);
  
  DataFrame newsupplement;
  DataFrame noage_supplement;
  
  if (wtf == 0) {
    newsupplement = LefkoMats::supp_reassess(ahstages, true, supplement);
  } else if (wtf == 1) {
    newsupplement = LefkoMats::supp_reassess(ahstages, false, supplement);
  } else if (wtf == 2) {
    noage_supplement = LefkoMats::supp_reassess(ahstages, false, supplement);
  } else {
    IntegerVector ah_minage = as<IntegerVector>(ahstages["min_age"]);
    IntegerVector ah_maxage = as<IntegerVector>(ahstages["max_age"]);
    int ahstages_length = static_cast<int>(ah_minage.length());
    
    int first_age2 = ah_minage(0);
    int final_age2 {0};
    if (IntegerVector::is_na(ah_maxage(ahstages_length - 1))) {
      final_age2 = first_age2 + ahstages_length - 1;
    } else {
      final_age2 = ah_maxage(ahstages_length - 1);
    }
    
    IntegerVector newsupp_vec (age2_length);
    for (int i = 0; i < static_cast<int>(age2_length); i++) {
      if (IntegerVector::is_na(age2_(i))) {
        if (age2_(i) < first_age2 || age2_(i) > final_age2) {
          throw Rcpp::exception("Some age2 inputs are outside of the modeled ages.", false);
        }
        newsupp_vec(i) = ahstages_length;
      } else { 
        newsupp_vec(i) = 1;
      }
    }
    int newsupp_length = sum(newsupp_vec);
    
    IntegerVector newsupp_age2 (newsupp_length);
    IntegerVector newsupp_estage2 (newsupp_length);
    NumericVector newsupp_givenrate (newsupp_length);
    NumericVector newsupp_offset (newsupp_length);
    NumericVector newsupp_multiplier (newsupp_length);
    IntegerVector newsupp_convtype (newsupp_length);
    StringVector newsupp_pop (newsupp_length);
    StringVector newsupp_patch (newsupp_length);
    StringVector newsupp_year2 (newsupp_length);
    
    int newsupp_counter {0};
    for (int i = 0; i < age2_length; i++) {
      for (int j = 0; j < newsupp_vec(i); j++) {
        if (IntegerVector::is_na(age2_(i))) {
          newsupp_age2(newsupp_counter) = first_age2 + j;
        } else {
          newsupp_age2(newsupp_counter) = age2_(i);
        }
        
        newsupp_estage2(newsupp_counter) = estage2_(i);
        newsupp_givenrate(newsupp_counter) = givenrate_(i);
        newsupp_offset(newsupp_counter) = offset_(i);
        newsupp_multiplier(newsupp_counter) = multiplier_(i);
        newsupp_convtype(newsupp_counter) = type_(i);
        newsupp_pop(newsupp_counter) = pop_(i);
        newsupp_patch(newsupp_counter) = patch_(i);
        newsupp_year2(newsupp_counter) = year2_(i);
        
        newsupp_counter++;
      }
    }
    
    DataFrame newsupp = DataFrame::create(_["age2"] = newsupp_age2,
      _["estage2"] = newsupp_estage2, _["givenrate"] = newsupp_givenrate,
      _["offset"] = newsupp_offset, _["multiplier"] = newsupp_multiplier,
      _["convtype"] = newsupp_convtype, _["pop"] = newsupp_pop,
      _["patch"] = newsupp_patch, _["year2"] = newsupp_year2);
    newsupplement = newsupp;
  }
  
  // Core element index calculation
  IntegerVector new_col_index_;
  IntegerVector new_row_index_;
  IntegerVector new_use_est_;
  IntegerVector new_est_col_index_;
  IntegerVector new_est_row_index_;
  NumericVector new_givenrate_;
  NumericVector new_offset_;
  NumericVector new_multiplier_;
  IntegerVector new_convtype_;
  
  if (wtf == 0) { // Historical
    NumericVector new_givenrate = as<NumericVector>(newsupplement["givenrate"]);
    NumericVector new_offset = as<NumericVector>(newsupplement["offset"]);
    NumericVector new_multiplier = as<NumericVector>(newsupplement["multiplier"]);
    IntegerVector new_convtype = as<IntegerVector>(newsupplement["convtype"]);
    
    StringVector new_stage3 = as<StringVector>(newsupplement["stage3"]);
    StringVector new_stage2 = as<StringVector>(newsupplement["stage2"]);
    StringVector new_stage1 = as<StringVector>(newsupplement["stage1"]);
    
    StringVector new_eststage3 = as<StringVector>(newsupplement["eststage3"]);
    StringVector new_eststage2 = as<StringVector>(newsupplement["eststage2"]);
    StringVector new_eststage1 = as<StringVector>(newsupplement["eststage1"]);
    
    int check_loop_length = static_cast<int>(new_stage3.length());
    
    arma::uvec hst_stage2_id = as<arma::uvec>(hstages["stage_id_2"]);
    arma::uvec hst_stage1_id = as<arma::uvec>(hstages["stage_id_1"]);
    IntegerVector ahst_stage_id = as<IntegerVector>(ahstages["stage_id"]);
    StringVector ahst_stage = as<StringVector>(ahstages["stage"]);
    int ahstages_length = static_cast<int>(ahst_stage.length());
    
    IntegerVector new_col_index (check_loop_length, -1);
    IntegerVector new_row_index (check_loop_length, -1);
    IntegerVector new_use_est (check_loop_length);
    IntegerVector new_est_col_index (check_loop_length, -1);
    IntegerVector new_est_row_index (check_loop_length, -1);
    
    for (int i = 0; i < check_loop_length; i++) {
      int found_stage3 {0};
      int found_stage2 {0};
      int found_stage1 {0};
      
      int found_eststage3 {0};
      int found_eststage2 {0};
      int found_eststage1 {0};
      
      for (int j = 0; j < ahstages_length; j++) {
        if (ahst_stage(j) == new_stage3(i)) found_stage3 = ahst_stage_id(j);
        if (ahst_stage(j) == new_stage2(i)) found_stage2 = ahst_stage_id(j);
        if (ahst_stage(j) == new_stage1(i)) found_stage1 = ahst_stage_id(j);
        
        if (!StringVector::is_na(new_eststage3(i)) && new_eststage3(i) != "NA") {
          new_use_est(i) = 1;
          if (ahst_stage(j) == new_eststage3(i)) found_eststage3 = ahst_stage_id(j);
          if (ahst_stage(j) == new_eststage2(i)) found_eststage2 = ahst_stage_id(j);
          if (ahst_stage(j) == new_eststage1(i)) found_eststage1 = ahst_stage_id(j);
        }
      }
      
      arma::uvec hst_found_stage3of32 = find(hst_stage2_id == found_stage3);
      arma::uvec hst_found_stage2of32 = find(hst_stage1_id == found_stage2);
      arma::uvec hst_found_stage2of21 = find(hst_stage2_id == found_stage2);
      arma::uvec hst_found_stage1of21 = find(hst_stage1_id == found_stage1);
      
      arma::uvec hst_found_stage32 = intersect(hst_found_stage3of32,
        hst_found_stage2of32);
      arma::uvec hst_found_stage21 = intersect(hst_found_stage2of21,
        hst_found_stage1of21);
      
      if (hst_found_stage32.n_elem > 0 && hst_found_stage21.n_elem > 0) {
        new_col_index(i) = static_cast<int>(hst_found_stage21(0));
        new_row_index(i) = static_cast<int>(hst_found_stage32(0));
      } else {
        Rf_warningcall(R_NilValue, "Some stage designations could not be found.");
      }
      
      if (new_use_est(i) == 1) {
        arma::uvec hst_found_eststage3of32 = find(hst_stage2_id == found_eststage3);
        arma::uvec hst_found_eststage2of32 = find(hst_stage1_id == found_eststage2);
        arma::uvec hst_found_eststage2of21 = find(hst_stage2_id == found_eststage2);
        arma::uvec hst_found_eststage1of21 = find(hst_stage1_id == found_eststage1);
        
        arma::uvec hst_found_eststage32 = intersect(hst_found_eststage3of32,
          hst_found_eststage2of32);
        arma::uvec hst_found_eststage21 = intersect(hst_found_eststage2of21,
          hst_found_eststage1of21);
        
        if (hst_found_eststage32.n_elem > 0 && hst_found_eststage21.n_elem > 0) {
          new_est_col_index(i) = static_cast<int>(hst_found_eststage21(0));
          new_est_row_index(i) = static_cast<int>(hst_found_eststage32(0));
        } else if (!StringVector::is_na(new_eststage3(i))) {
          Rf_warningcall(R_NilValue, "Some eststage designations could not be found.");
        }
      }
    }
    
    new_col_index_ = new_col_index;
    new_row_index_ = new_row_index;
    new_use_est_ = new_use_est;
    new_est_col_index_ = new_est_col_index;
    new_est_row_index_ = new_est_row_index;
    new_givenrate_ = new_givenrate;
    new_offset_ = new_offset;
    new_multiplier_ = new_multiplier;
    new_convtype_ = new_convtype;
    
  } else if (wtf == 1) { // Ahistorical
    NumericVector new_givenrate = as<NumericVector>(newsupplement["givenrate"]);
    NumericVector new_offset = as<NumericVector>(newsupplement["offset"]);
    NumericVector new_multiplier = as<NumericVector>(newsupplement["multiplier"]);
    IntegerVector new_convtype = as<IntegerVector>(newsupplement["convtype"]);
    
    StringVector new_stage3 = as<StringVector>(newsupplement["stage3"]);
    StringVector new_stage2 = as<StringVector>(newsupplement["stage2"]);
    
    StringVector new_eststage3 = as<StringVector>(newsupplement["eststage3"]);
    StringVector new_eststage2 = as<StringVector>(newsupplement["eststage2"]);
    
    int check_loop_length = static_cast<int>(new_stage3.length());
    
    IntegerVector ahst_stage_id = as<IntegerVector>(ahstages["stage_id"]);
    StringVector ahst_stage = as<StringVector>(ahstages["stage"]);
    int agestages_length = static_cast<int>(ahst_stage.length());
    
    IntegerVector new_col_index (check_loop_length, -1);
    IntegerVector new_row_index (check_loop_length, -1);
    IntegerVector new_use_est (check_loop_length);
    IntegerVector new_est_col_index (check_loop_length, -1);
    IntegerVector new_est_row_index (check_loop_length, -1);
    
    for (int i = 0; i < check_loop_length; i++) {
      int found_stage2 {0};
      int found_stage3 {0};
      bool found_stage2_log {false};
      bool found_stage3_log {false};
      
      int found_eststage2 {0};
      int found_eststage3 {0};
      bool found_eststage2_log {false};
      bool found_eststage3_log {false};
      
      for (int j = 0; j < agestages_length; j++) {
        if (ahst_stage(j) == new_stage3(i)) {
          found_stage3 = j;
          found_stage3_log = true;
        }
        if (ahst_stage(j) == new_stage2(i)) {
          found_stage2 = j;
          found_stage2_log = true;
        }
        
        if (!StringVector::is_na(new_eststage3(i)) && new_eststage3(i) != "NA") {
          new_use_est(i) = 1;
          if (ahst_stage(j) == new_eststage3(i)) {
            found_eststage3 = j;
            found_eststage3_log = true;
          }
          if (ahst_stage(j) == new_eststage2(i)) {
            found_eststage2 = j;
            found_eststage2_log = true;
          }
        }
      }
      
      if (found_stage3_log && found_stage2_log) {
        new_col_index(i) = found_stage2;
        new_row_index(i) = found_stage3;
      } else {
        Rf_warningcall(R_NilValue, "Some stage designations could not be found.");
      }
      
      if (new_use_est(i) == 1) {
        if (found_eststage3_log && found_eststage2_log) {
          new_est_col_index(i) = found_eststage2;
          new_est_row_index(i) = found_eststage3;
        } else {
          Rf_warningcall(R_NilValue, "Some eststage designations could not be found.");
        }
      }
    }
    
    new_col_index_ = new_col_index;
    new_row_index_ = new_row_index;
    new_use_est_ = new_use_est;
    new_est_col_index_ = new_est_col_index;
    new_est_row_index_ = new_est_row_index;
    new_givenrate_ = new_givenrate;
    new_offset_ = new_offset;
    new_multiplier_ = new_multiplier;
    new_convtype_ = new_convtype;
    
  } else if (wtf == 2) {
    arma::uvec agst_stage_id = as<arma::uvec>(agestages["stage_id"]);
    StringVector agst_stage = as<StringVector>(agestages["stage"]);
    arma::uvec agst_age = as<arma::uvec>(agestages["age"]);
    int agestages_length = static_cast<int>(agst_stage_id.n_elem);
    
    min_age = static_cast<int>(agst_age.min());
    max_age = static_cast<int>(agst_age.max());
    
    newsupplement = LefkoMats::age_expanded(noage_supplement, min_age, max_age);
    
    NumericVector new_givenrate = as<NumericVector>(newsupplement["givenrate"]);
    NumericVector new_offset = as<NumericVector>(newsupplement["offset"]);
    NumericVector new_multiplier = as<NumericVector>(newsupplement["multiplier"]);
    IntegerVector new_convtype = as<IntegerVector>(newsupplement["convtype"]);
    
    IntegerVector new_age2 = as<IntegerVector>(newsupplement["age2"]);
    StringVector new_stage3 = as<StringVector>(newsupplement["stage3"]);
    StringVector new_stage2 = as<StringVector>(newsupplement["stage2"]);
    
    IntegerVector new_estage2 = as<IntegerVector>(newsupplement["estage2"]);
    StringVector new_eststage3 = as<StringVector>(newsupplement["eststage3"]);
    StringVector new_eststage2 = as<StringVector>(newsupplement["eststage2"]);
    
    int check_loop_length = static_cast<int>(new_stage2.length());
    
    IntegerVector new_col_index (check_loop_length, -1);
    IntegerVector new_row_index (check_loop_length, -1);
    IntegerVector new_use_est (check_loop_length);
    IntegerVector new_est_col_index (check_loop_length, -1);
    IntegerVector new_est_row_index (check_loop_length, -1);
    
    for (int i = 0; i < check_loop_length; i++) {
      int found_age2 {0};
      int found_stage3 {0};
      int found_stage2 {0};
      
      int found_estage2 {0};
      int found_eststage3 {0};
      int found_eststage2 {0};
      
      for (int j = 0; j < agestages_length; j++) {
        if (agst_stage(j) == new_stage3(i)) found_stage3 = agst_stage_id(j);
        if (agst_stage(j) == new_stage2(i)) found_stage2 = agst_stage_id(j);
        if (static_cast<int>(agst_age(j)) == static_cast<int>(new_age2(i))) found_age2 = agst_age(j);
        
        if (!IntegerVector::is_na(new_estage2(i))) {
          new_use_est(i) = 1;
          if (agst_stage(j) == new_eststage3(i)) found_eststage3 = agst_stage_id(j);
          if (agst_stage(j) == new_eststage2(i)) found_eststage2 = agst_stage_id(j);
          if (static_cast<int>(agst_age(j)) == static_cast<int>(new_estage2(i))) found_estage2 = agst_age(j);
        }
      }
      
      arma::uvec agst_found_stage3of3a2 = find(agst_stage_id == found_stage3);
      arma::uvec agst_found_stage2of3a2;
      if (found_age2 < max_age) {
        agst_found_stage2of3a2= find(agst_age == (found_age2 + 1));
      } else {
        agst_found_stage2of3a2= find(agst_age == found_age2);
      }
      arma::uvec agst_found_stage2of2a1 = find(agst_stage_id == found_stage2);
      arma::uvec agst_found_stage1of2a1 = find(agst_age == found_age2);
      
      arma::uvec agst_found_stage3a2 = intersect(agst_found_stage3of3a2,
        agst_found_stage2of3a2);
      arma::uvec agst_found_stage2a1 = intersect(agst_found_stage2of2a1,
        agst_found_stage1of2a1);
      
      if (agst_found_stage3a2.n_elem > 0 && agst_found_stage2a1.n_elem > 0) {
        new_col_index(i) = static_cast<int>(agst_found_stage2a1(0));
        new_row_index(i) = static_cast<int>(agst_found_stage3a2(0));
      } else {
        Rf_warningcall(R_NilValue, "Some stage designations could not be found.");
      }
      
      if (new_use_est(i) == 1) {
        arma::uvec agst_found_eststage3of3a2 = find(agst_stage_id == found_eststage3);
        arma::uvec agst_found_eststage2of3a2;
        if (found_estage2 < max_age) {
          agst_found_eststage2of3a2 = find(agst_age == (found_estage2 + 1));
        } else {
          agst_found_eststage2of3a2 = find(agst_age == found_estage2);
        }
        arma::uvec agst_found_eststage2of2a1 = find(agst_stage_id == found_eststage2);
        arma::uvec agst_found_eststage1of2a1 = find(agst_age == found_estage2);
        
        arma::uvec agst_found_eststage3a2 = intersect(agst_found_eststage3of3a2,
        agst_found_eststage2of3a2);
        arma::uvec agst_found_eststage2a1 = intersect(agst_found_eststage2of2a1,
        agst_found_eststage1of2a1);
        
        if (agst_found_eststage3a2.n_elem > 0 && agst_found_eststage2a1.n_elem > 0) {
          new_est_col_index(i) = static_cast<int>(agst_found_eststage2a1(0));
          new_est_row_index(i) = static_cast<int>(agst_found_eststage3a2(0));
        } else {
          Rf_warningcall(R_NilValue, "Some eststage designations could not be found.");
        }
      }
    }
    
    new_col_index_ = new_col_index;
    new_row_index_ = new_row_index;
    new_use_est_ = new_use_est;
    new_est_col_index_ = new_est_col_index;
    new_est_row_index_ = new_est_row_index;
    new_givenrate_ = new_givenrate;
    new_offset_ = new_offset;
    new_multiplier_ = new_multiplier;
    new_convtype_ = new_convtype;
    
  } else { // Leslie matrix
    NumericVector new_givenrate = as<NumericVector>(newsupplement["givenrate"]);
    NumericVector new_offset = as<NumericVector>(newsupplement["offset"]);
    NumericVector new_multiplier = as<NumericVector>(newsupplement["multiplier"]);
    IntegerVector new_convtype = as<IntegerVector>(newsupplement["convtype"]);
    
    IntegerVector new_age2 = as<IntegerVector>(newsupplement["age2"]);
    IntegerVector new_estage2 = as<IntegerVector>(newsupplement["estage2"]);
    
    int check_loop_length = static_cast<int>(new_age2.length());
    
    IntegerVector ahst_stage_id = as<IntegerVector>(ahstages["stage_id"]);
    IntegerVector ahst_minage = as<IntegerVector>(ahstages["min_age"]);
    IntegerVector ahst_maxage = as<IntegerVector>(ahstages["max_age"]);
    int ages_length = static_cast<int>(ahst_stage_id.length());
    
    int first_age2 = ahst_minage(0);
    int final_age2 {0};
    if (IntegerVector::is_na(ahst_maxage(ages_length - 1))) {
      final_age2 = first_age2 + ages_length - 1;
    } else {
      final_age2 = ahst_maxage(ages_length - 1);
    }
    IntegerVector all_age2 = seq(first_age2, final_age2);
    arma::ivec all_age2_arma = as<arma::ivec>(all_age2);
    
    IntegerVector new_col_index (check_loop_length, -1);
    IntegerVector new_row_index (check_loop_length, -1);
    IntegerVector new_use_est (check_loop_length);
    IntegerVector new_est_col_index (check_loop_length, -1);
    IntegerVector new_est_row_index (check_loop_length, -1);
    
    for (int i = 0; i < check_loop_length; i++) {
      arma::uvec found_age_elems = find(all_age2_arma == new_age2(i));
      int found_age_elem0 = static_cast<int>(found_age_elems(0));
      int found_age2 = all_age2(found_age_elem0);
      
      new_col_index(i) = found_age_elem0;
      if (new_convtype(i) == 1) {
        if (found_age2 < final_age2) {
          new_row_index(i) = found_age_elem0 + 1;
        } else {
          new_row_index(i) = found_age_elem0;
        }
      } else {
        new_row_index(i) = 0;
      }
      
      if (!IntegerVector::is_na(new_estage2(i))) {
        arma::uvec found_estage_elems = find(all_age2_arma == new_estage2(i));
        
        if (found_estage_elems.n_elem > 0) {
          int found_estage_elem0 = static_cast<int>(found_estage_elems(0));
          int found_estage2 = all_age2(found_estage_elem0);
          
          new_est_col_index(i) = found_estage_elem0;
          if (new_convtype(i) == 1) {
            if (found_estage2 < final_age2) {
              new_est_row_index(i) = found_estage_elem0 + 1;
            } else {
              new_est_row_index(i) = found_estage_elem0;
            }
          } else {
            new_est_row_index(i) = 0;
          }
          
          new_use_est(i) = 1;
        } else {
          throw Rcpp::exception("Value entered for estage2 cannot be found.", false);
        }
      }
    }
      
    new_col_index_ = new_col_index;
    new_row_index_ = new_row_index;
    new_use_est_ = new_use_est;
    new_est_col_index_ = new_est_col_index;
    new_est_row_index_ = new_est_row_index;
    new_givenrate_ = new_givenrate;
    new_offset_ = new_offset;
    new_multiplier_ = new_multiplier;
    new_convtype_ = new_convtype;
  }
  
  // Core matrix editing
  List A_mats;
  List U_mats;
  List F_mats;
  
  bool A_used {false};
  bool U_used {false};
  bool F_used {false};
  
  bool mat_sparse {false};
  int mat_num {0};
  
  if (mpm.containsElementNamed("A")) {
    if (is<List>(mpm["A"])) {
      List new_A_mats = as<List>(mpm["A"]);
      A_mats = clone(new_A_mats);
      A_used = true;
      mat_num = static_cast<int>(A_mats.length());
      
      if (is<S4>(A_mats(0))) mat_sparse = true;
    }
    if (is<List>(mpm["U"])) {
      List new_U_mats = as<List>(mpm["U"]);
      U_mats = clone(new_U_mats);
      U_used = true;
      mat_num = static_cast<int>(U_mats.length());
      
      if (is<S4>(U_mats(0))) mat_sparse = true;
    }
    if (is<List>(mpm["F"])) {
      List new_F_mats = as<List>(mpm["F"]);
      F_mats = clone(new_F_mats);
      F_used = true;
      mat_num = static_cast<int>(F_mats.length());
      
      if (is<S4>(F_mats(0))) mat_sparse = true;
    }
  }
  
  if (!A_used && !U_used && !F_used) {
    throw Rcpp::exception("This input object does not appear to hold matrices.", false);
  }
  
  StringVector new_pop = as<StringVector>(newsupplement["pop"]);
  StringVector new_patch = as<StringVector>(newsupplement["patch"]);
  StringVector new_year2 = as<StringVector>(newsupplement["year2"]);
  
  StringVector labels_pop = labels["pop"];
  StringVector labels_patch = labels["patch"];
  StringVector labels_year2 = labels["year2"];
  
  int labels_length = static_cast<int>(labels_pop.length());
  arma::uvec main_labels_index (labels_length, fill::ones);
  
  int new_supp_length = static_cast<int>(new_col_index_.length());
  int U_adjustment {0};
  int F_adjustment {0};
  
  for (int i = 0; i < new_supp_length; i++) {
    arma::ivec chosen_pops (labels_length);
    arma::ivec chosen_patches (labels_length);
    arma::ivec chosen_year2s (labels_length);
    chosen_pops.fill(-1);
    chosen_patches.fill(-1);
    chosen_year2s.fill(-1);
    
    if (StringVector::is_na(new_pop(i)) || new_pop(i) == "NA" || new_pop(i) == "all") {
      chosen_pops = arma::conv_to<arma::ivec>::from(main_labels_index);
    } else {
      for (int j = 0; j < labels_length; j++) {
        if (new_pop(i) == labels_pop(j)) chosen_pops(j) = j;
      }
    }
    if (StringVector::is_na(new_patch(i)) || new_patch(i) == "NA" || new_patch(i) == "all") {
      chosen_patches = arma::conv_to<arma::ivec>::from(main_labels_index);
    } else {
      for (int j = 0; j < labels_length; j++) {
        if (new_patch(i) == labels_patch(j)) chosen_patches(j) = j;
      }
    }
    if (StringVector::is_na(new_year2(i)) || new_year2(i) == "NA" || new_year2(i) == "all") {
      chosen_year2s = arma::conv_to<arma::ivec>::from(main_labels_index);
    } else {
      for (int j = 0; j < labels_length; j++) {
        if (new_year2(i) == labels_year2(j)) chosen_year2s(j) = j;
      }
    }
    
    arma::uvec pop_indices = find(chosen_pops > -1);
    arma::uvec patch_indices = find(chosen_patches > -1);
    arma::uvec year2_indices = find(chosen_year2s > -1);
    
    if (pop_indices.n_elem == 0) {
      Rf_warningcall(R_NilValue,
        "Pop designations could not be found in input MPM.");
    }
    if (patch_indices.n_elem == 0) {
      Rf_warningcall(R_NilValue,
        "Patch designations could not be found in input MPM.");
    }
    if (year2_indices.n_elem == 0) {
      Rf_warningcall(R_NilValue,
        "Year2 designations could not be found in input MPM.");
    }
    
    arma::uvec pop_patch_intersect = intersect(pop_indices, patch_indices);
    arma::uvec chosen_matrices = intersect(pop_patch_intersect, year2_indices);
    int chosen_matrices_num = static_cast<int>(chosen_matrices.n_elem);
    
    for (int j = 0; j < chosen_matrices_num; j++) {
      if (U_used && F_used) {
        if (!mat_sparse) {
          arma::mat chosen_U = as<arma::mat>(U_mats(chosen_matrices(j)));
          arma::mat chosen_F = as<arma::mat>(F_mats(chosen_matrices(j)));
          
          if (new_convtype_(i) == 1 && new_col_index_(i) != -1) {
            if (!NumericVector::is_na(new_givenrate_(i))) {
              if (chosen_U(new_row_index_(i), new_col_index_(i)) == 0.0 &&
                  new_givenrate_(i) != 0.0) {
                U_adjustment += 1;
              } else if (chosen_U(new_row_index_(i), new_col_index_(i)) != 0.0 &&
                  new_givenrate_(i) == 0.0) {
                U_adjustment -= 1;
              }
              
              chosen_U(new_row_index_(i), new_col_index_(i)) = new_givenrate_(i);
            }
            if (!NumericVector::is_na(new_offset_(i))) {
              if (new_offset_(i) != 0.0) {
                chosen_U(new_row_index_(i), new_col_index_(i)) += new_offset_(i);
                U_adjustment += 1;
                
                if (chosen_U(new_row_index_(i), new_col_index_(i)) < 0.0) {
                  Rf_warningcall(R_NilValue,
                    "Some numeric offsets have produced negative matrix elements (survival-transition).");
                }
              }
            }
            if (new_use_est_(i) == 1) {
              chosen_U(new_row_index_(i), new_col_index_(i)) =
                chosen_U(new_est_row_index_(i), new_est_col_index_(i));
            }
            if (!NumericVector::is_na(new_multiplier_(i))) {
              chosen_U(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
            }
            
            U_mats(chosen_matrices(j)) = chosen_U;
            
          } else if (new_convtype_(i) == 2 && new_col_index_(i) != -1) {
            if (!NumericVector::is_na(new_givenrate_(i))) {
              if (chosen_F(new_row_index_(i), new_col_index_(i)) == 0.0 &&
                  new_givenrate_(i) != 0.0) {
                F_adjustment += 1;
              } else if (chosen_F(new_row_index_(i), new_col_index_(i)) != 0.0 &&
                  new_givenrate_(i) == 0.0) {
                F_adjustment -= 1;
              }
              chosen_F(new_row_index_(i), new_col_index_(i)) = new_givenrate_(i);
            }
            if (!NumericVector::is_na(new_offset_(i))) {
              if (new_offset_(i) != 0.0) {
                chosen_F(new_row_index_(i), new_col_index_(i)) += new_offset_(i);
                F_adjustment += 1;
                
                if (chosen_F(new_row_index_(i), new_col_index_(i)) < 0.0) {
                  Rf_warningcall(R_NilValue,
                    "Some numeric offsets have produced negative matrix elements (fecundity).");
                }
              }
            }
            if (new_use_est_(i) == 1) {
              chosen_F(new_row_index_(i), new_col_index_(i)) =
                chosen_F(new_est_row_index_(i), new_est_col_index_(i));
            }
            if (!NumericVector::is_na(new_multiplier_(i))) {
              chosen_F(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
            }
            
            F_mats(chosen_matrices(j)) = chosen_F;
            
          } else if (new_convtype_(i) == 3 && new_col_index_(i) != -1) {
            if (!NumericVector::is_na(new_multiplier_(i))) {
              
              double replacement_value = chosen_F(new_row_index_(i), new_col_index_(i)) * new_multiplier_(i);
              chosen_F(new_row_index_(i), new_col_index_(i)) = replacement_value;
            }
            
            F_mats(chosen_matrices(j)) = chosen_F;
          }
          
        } else {
          arma::sp_mat chosen_U = as<arma::sp_mat>(U_mats(chosen_matrices(j)));
          arma::sp_mat chosen_F = as<arma::sp_mat>(F_mats(chosen_matrices(j)));
          
          if (new_convtype_(i) == 1 && new_col_index_(i) != -1) {
            if (!NumericVector::is_na(new_givenrate_(i))) {
              if (chosen_U(new_row_index_(i), new_col_index_(i)) == 0.0 &&
                  new_givenrate_(i) != 0.0) {
                U_adjustment += 1;
              } else if (chosen_U(new_row_index_(i), new_col_index_(i)) != 0.0 &&
                  new_givenrate_(i) == 0.0) {
                U_adjustment -= 1;
              }
              chosen_U(new_row_index_(i), new_col_index_(i)) = new_givenrate_(i);
            }
            if (!NumericVector::is_na(new_offset_(i))) {
              if (new_offset_(i) != 0.0) {
                chosen_U(new_row_index_(i), new_col_index_(i)) += new_offset_(i);
                U_adjustment += 1;
                
                if (chosen_U(new_row_index_(i), new_col_index_(i)) < 0.0) {
                  Rf_warningcall(R_NilValue,
                    "Some numeric offsets have produced negative matrix elements (survival-transition).");
                }
              }
            }
            if (new_use_est_(i) == 1) {
              chosen_U(new_row_index_(i), new_col_index_(i)) =
                chosen_U(new_est_row_index_(i), new_est_col_index_(i));
            }
            if (!NumericVector::is_na(new_multiplier_(i))) {
              chosen_U(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
            }
            
            U_mats(chosen_matrices(j)) = chosen_U;
            
          } else if (new_convtype_(i) == 2 && new_col_index_(i) != -1) {
            if (!NumericVector::is_na(new_givenrate_(i))) {
              if (chosen_F(new_row_index_(i), new_col_index_(i)) == 0.0 &&
                  new_givenrate_(i) != 0.0) {
                F_adjustment += 1;
              } else if (chosen_F(new_row_index_(i), new_col_index_(i)) != 0.0 &&
                  new_givenrate_(i) == 0.0) {
                F_adjustment -= 1;
              }
              chosen_F(new_row_index_(i), new_col_index_(i)) = new_givenrate_(i);
            }
            if (!NumericVector::is_na(new_offset_(i))) {
              if (new_offset_(i) != 0.0) {
                chosen_F(new_row_index_(i), new_col_index_(i)) += new_offset_(i);
                F_adjustment += 1;
                
                if (chosen_F(new_row_index_(i), new_col_index_(i)) < 0.0) {
                  Rf_warningcall(R_NilValue,
                    "Some numeric offsets have produced negative matrix elements (fecundity).");
                }
              }
            }
            if (new_use_est_(i) == 1) {
              chosen_F(new_row_index_(i), new_col_index_(i)) =
                chosen_F(new_est_row_index_(i), new_est_col_index_(i));
            }
            if (!NumericVector::is_na(new_multiplier_(i))) {
              chosen_F(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
            }
            
            F_mats(chosen_matrices(j)) = chosen_F;
            
          } else if (new_convtype_(i) == 3 && new_col_index_(i) != -1) {
            if (!NumericVector::is_na(new_multiplier_(i))) {
              chosen_F(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
            }
            
            F_mats(chosen_matrices(j)) = chosen_F;
          }
        }
      } else if (A_used) {
        if (!mat_sparse) {
          arma::mat chosen_A = as<arma::mat>(A_mats(chosen_matrices(j)));
          
          if (!NumericVector::is_na(new_givenrate_(i)) && new_col_index_(i) != -1) {
            if (chosen_A(new_row_index_(i), new_col_index_(i)) == 0.0 &&
                new_givenrate_(i) != 0.0) {
              U_adjustment += 1;
            } else if (chosen_A(new_row_index_(i), new_col_index_(i)) != 0.0 &&
                new_givenrate_(i) == 0.0) {
              U_adjustment -= 1;
            }
            chosen_A(new_row_index_(i), new_col_index_(i)) = new_givenrate_(i);
          }
          if (!NumericVector::is_na(new_offset_(i))) {
            if (new_offset_(i) != 0.0) {
              chosen_A(new_row_index_(i), new_col_index_(i)) += new_offset_(i);
              U_adjustment += 1;
              
              if (chosen_A(new_row_index_(i), new_col_index_(i)) < 0.0) {
                Rf_warningcall(R_NilValue,
                  "Some numeric offsets have produced negative matrix elements.");
              }
            }
          }
          if (new_use_est_(i) == 1) {
            chosen_A(new_row_index_(i), new_col_index_(i)) =
              chosen_A(new_est_row_index_(i), new_est_col_index_(i));
          }
          if (!NumericVector::is_na(new_multiplier_(i)) && new_col_index_(i) != -1) {
            chosen_A(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
          }
          
          A_mats(chosen_matrices(j)) = chosen_A;
          
        } else {
          arma::sp_mat chosen_A = as<arma::sp_mat>(A_mats(chosen_matrices(j)));
          
          if (!NumericVector::is_na(new_givenrate_(i)) && new_col_index_(i) != -1) {
            if (chosen_A(new_row_index_(i), new_col_index_(i)) == 0.0 &&
                new_givenrate_(i) != 0.0) {
              U_adjustment += 1;
            } else if (chosen_A(new_row_index_(i), new_col_index_(i)) != 0.0 &&
                new_givenrate_(i) == 0.0) {
              U_adjustment -= 1;
            }
            chosen_A(new_row_index_(i), new_col_index_(i)) = new_givenrate_(i);
          }
          if (!NumericVector::is_na(new_offset_(i))) {
            if (new_offset_(i) != 0.0) {
              chosen_A(new_row_index_(i), new_col_index_(i)) += new_offset_(i);
              U_adjustment += 1;
              
              if (chosen_A(new_row_index_(i), new_col_index_(i)) < 0.0) {
                Rf_warningcall(R_NilValue,
                  "Some numeric offsets have produced negative matrix elements.");
              }
            }
          }
          if (new_use_est_(i) == 1) {
            chosen_A(new_row_index_(i), new_col_index_(i)) =
              chosen_A(new_est_row_index_(i), new_est_col_index_(i));
          }
          if (!NumericVector::is_na(new_multiplier_(i)) && new_col_index_(i) != -1) {
            chosen_A(new_row_index_(i), new_col_index_(i)) *= new_multiplier_(i);
          }
          
          A_mats(chosen_matrices(j)) = chosen_A;
        }
      }
    }
  }
  
  if (U_used && F_used) {
    for (int i = 0; i < mat_num; i++) {
      if (!mat_sparse) { 
        arma::mat current_U = as<arma::mat>(U_mats(i));
        arma::mat current_F = as<arma::mat>(F_mats(i));
        arma::mat current_A = current_U + current_F;
        A_mats(i) = current_A;
      } else {
        arma::sp_mat current_U = as<arma::sp_mat>(U_mats(i));
        arma::sp_mat current_F = as<arma::sp_mat>(F_mats(i));
        arma::sp_mat current_A = current_U + current_F;
        A_mats(i) = current_A;
      }
    }
  }
  
  IntegerVector dataqc;
  if (mpm.containsElementNamed("dataqc")) {
    dataqc = as<IntegerVector>(mpm["dataqc"]);
  } else {
    IntegerVector mayhem_highway = {NA_INTEGER, NA_INTEGER};
    dataqc = mayhem_highway;
  }
  
  IntegerVector matrixqc = as<IntegerVector>(mpm["matrixqc"]);
  
  int total_Us = static_cast<int>(matrixqc(0)) + U_adjustment;
  int total_Fs = static_cast<int>(matrixqc(1)) + F_adjustment;
  
  IntegerVector new_matrixqc = {total_Us, total_Fs, mat_num};
  
  List true_output;
  
  if (mpm.containsElementNamed("modelqc")) {
    DataFrame modelqc = as<DataFrame>(mpm["modelqc"]);
    
    List new_mpm (10);
    new_mpm(0) = A_mats;
    new_mpm(1) = U_mats;
    new_mpm(2) = F_mats;
    new_mpm(3) = ahstages;
    new_mpm(4) = hstages;
    new_mpm(5) = agestages;
    new_mpm(6) = labels;
    new_mpm(7) = dataqc;
    new_mpm(8) = new_matrixqc;
    new_mpm(9) = modelqc;
    
    final_output = new_mpm;
    
    StringVector new_mpm_names = {"A", "U", "F", "ahstages", "hstages",
      "agestages", "labels", "dataqc", "matrixqc", "modelqc"};
    StringVector new_mpm_class = {"lefkoMat"};
    final_output.attr("class") = new_mpm_class;
    final_output.attr("names") = new_mpm_names;
    
  } else {
    List new_mpm (9);
    new_mpm(0) = A_mats;
    new_mpm(1) = U_mats;
    new_mpm(2) = F_mats;
    new_mpm(3) = ahstages;
    new_mpm(4) = hstages;
    new_mpm(5) = agestages;
    new_mpm(6) = labels;
    new_mpm(7) = dataqc;
    new_mpm(8) = new_matrixqc;
    
    final_output = new_mpm;
    
    StringVector new_mpm_names = {"A", "U", "F", "ahstages", "hstages",
      "agestages", "labels", "dataqc", "matrixqc"};
    StringVector new_mpm_class = {"lefkoMat"};
    final_output.attr("class") = new_mpm_class;
    final_output.attr("names") = new_mpm_names;
  }
}

//' Edit lefkoMat or lefkoMatList Object based on Supplemental Data
//' 
//' Function \code{edit_lM()} edits existing \code{lefkoMat} and
//' \code{lefkoMatList} objects with external data supplied by the user. The
//' effects are similar to function \code{\link{supplemental}()}, though
//' function \code{edit_lM()} allows individuals matrices within \code{lefkoMat}
//' objects to be edited after creation, while \code{\link{supplemental}()}
//' provides external data that modifies all matrices within a \code{lefkoMat}
//' object, or within all the \code{lefkoMat} objects within a
//' \code{lefkoMatList} object.
//' 
//' @name edit_lM
//' 
//' @param mpm The \code{lefkoMat} and \code{lefkoMatList} object to be edited.
//' @param pop A string vector denoting the populations to be edited. Defaults
//' to \code{NULL}, in which case all populations are edited.
//' @param patch A string vector denoting the patches to be edited. Defaults
//' to \code{NULL}, in which case all patches are edited.
//' @param year2 A string vector denoting the years to be edited. Defaults
//' to \code{NULL}, in which case all years are edited.
//' @param stage3 The name of the stage in occasion \emph{t}+1 in the transition
//' to be replaced. Abbreviations for groups of stages are also usable (see
//' \code{Notes}). Required in all stage-based and age-by-stage MPMs.
//' @param stage2 The name of the stage in occasion \emph{t} in the transition
//' to be replaced. Abbreviations for groups of stages are also usable (see
//' \code{Notes}). Required in all stage-based and age-by-stage MPMs.
//' @param stage1 The name of the stage in occasion \emph{t}-1 in the transition
//' to be replaced. Only needed if a historical matrix is to be produced.
//' Abbreviations for groups of stages are also usable (see \code{Notes}).
//' Required for historical stage-based MPMs.
//' @param age2 An integer vector of the ages in occasion \emph{t} to use in
//' transitions to be changed or replaced. Required for all age- and
//' age-by-stage MPMs.
//' @param eststage3 The name of the stage to replace \code{stage3} in a proxy
//' transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage2 The name of the stage to replace \code{stage2} in a proxy
//' transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in stage-based and age-by-stage MPMs.
//' @param eststage1 The name of the stage to replace \code{stage1} in a proxy
//' historical transition. Only needed if a transition will be replaced by
//' another estimated transition, and the matrix to be estimated is historical
//' and stage-based. Stage \code{NotAlive} is also possible for raw hMPMs as a
//' means of handling the prior stage for individuals entering the population in
//' occasion \emph{t}.
//' @param estage2 The age at time \emph{t} to replace \code{age2} in a proxy
//' transition. Only needed if a transition will be replaced by another
//' estimated transition, and only in age-based and age-by-stage MPMs.
//' @param givenrate A fixed rate or probability to replace for the transition
//' described by \code{stage3}, \code{stage2}, and \code{stage1}.
//' @param offset A numeric vector of fixed numeric values to add to the
//' transitions described by \code{stage3}, \code{stage2}, \code{stage1}, and/or
//' \code{age2}.
//' @param multiplier A vector of numeric multipliers for fecundity or for proxy
//' transitions. Defaults to \code{1}.
//' @param type A vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition;
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition; or \code{3}, \code{R}, or \code{r} for a fecundity multiplier.
//' If empty or not provided, then defaults to \code{1} for survival transition.
//' @param type_t12 An optional vector denoting the kind of transition between
//' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
//' deVries format is desired. This should be entered as \code{1}, \code{S}, or
//' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
//' fecundity transitions. Defaults to \code{1} for survival transition, with
//' impacts only on the construction of deVries-format hMPMs.
//' @param target_mpm If modifying a \code{lefkoMatList} object, then this
//' allows the user to specify which MPMs to modify. To modify, enter a vector
//' with the number of each MPM to modify, or enter \code{"all"} to modify all
//' MPMs. Defaults to \code{"all"}.
//' 
//' @return An edited copy of the original MPM is returned, also as a
//' \code{lefkoMat} object.
//' 
//' @section Notes:
//' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
//' abbreviations for groups of stages. Use \code{rep} if all reproductive
//' stages are to be used, \code{nrep} if all mature but non-reproductive stages
//' are to be used, \code{mat} if all mature stages are to be used, \code{immat}
//' if all immature stages are to be used, \code{prop} if all propagule stages
//' are to be used, \code{npr} if all non-propagule stages are to be used,
//' \code{obs} if all observable stages are to be used, \code{nobs} if all
//' unobservable stages are to be used, and leave empty or use \code{all} if all
//' stages in stageframe are to be used. Also use \code{groupX} to denote all
//' stages in group X (e.g. \code{group1} will use all stages in the respective
//' stageframe's group 1).
//' 
//' @seealso \code{\link{supplemental}()}
//' 
//' @examples
//' data(cypdata)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   age_offset = 3, NAas0 = TRUE, NRasRep = TRUE)
//' 
//' cyp_rl <- rleslie(data = cypraw_v1, start_age = 0, last_age = 6, continue = TRUE,
//'   fecage_min = 3, year = "all", pop = NA, patch = "all", yearcol = "year2",
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' ddd1 <- edit_lM(cyp_rl, age2 = c(0, 1, 2, 3, 4, 5, 6),
//'   givenrate = c(0.25, 0.25, 0.4, 0.4, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 2000, 2000, 2000),
//'   type = c(1, 1, 1, 1, 3, 3, 3))
//'   
//' ddd1 <- edit_lM(ddd1, age2 = 6, multiplier = 1.5, type = 3, patch = "B",
//'   year2 = "2005")
//' 
//' @export edit_lM
// [[Rcpp::export(edit_lM)]]
Rcpp::List edit_lM (const RObject mpm, Nullable<RObject> pop = R_NilValue,
  Nullable<RObject> patch = R_NilValue, Nullable<RObject> year2 = R_NilValue,
  Nullable<RObject> stage3 = R_NilValue, Nullable<RObject> stage2 = R_NilValue,
  Nullable<RObject> stage1 = R_NilValue, Nullable<RObject> age2 = R_NilValue,
  Nullable<RObject> eststage3 = R_NilValue, Nullable<RObject> eststage2 = R_NilValue,
  Nullable<RObject> eststage1 = R_NilValue, Nullable<RObject> estage2 = R_NilValue,
  Nullable<RObject> givenrate = R_NilValue, Nullable<RObject> offset = R_NilValue,
  Nullable<RObject> multiplier = R_NilValue, Nullable<RObject> type = R_NilValue,
  Nullable<RObject> type_t12 = R_NilValue, Nullable<RObject> target_mpm = R_NilValue) {
  
  List mpm_list;
  List true_output;
  
  if (is<List>(mpm)) mpm_list = as<List>(mpm);
  StringVector mpm_class = mpm_list.attr("class");
  
  String mpm_error = "Please enter a lefkoMat or lefkoMatList object as input.";
  bool mpm_yes {false};
  bool mpmlist_yes {false};
  
  for (int i = 0; i < static_cast<int>(mpm_class.length()); i++) {
    if (mpm_class(i) == "lefkoMat") mpm_yes = true;
    if (mpm_class(i) == "lefkoMatList") mpmlist_yes = true;
  }
  if (!mpm_yes && !mpmlist_yes) throw Rcpp::exception(mpm_error.get_cstring(), false);
  
  
  if (mpm_yes && !mpmlist_yes) {
    edit_lM_single(true_output, mpm_list, pop, patch, year2, stage3, stage2,
      stage1, age2, eststage3, eststage2, eststage1, estage2, givenrate, offset,
      multiplier, type, type_t12);
  } else if (mpmlist_yes) {
    int mpmlist_length = static_cast<int>(mpm_list.length());
    List semifinal_output (mpmlist_length);
    
    for (int i = 0; i < mpmlist_length; i++) {
      List single_elem_output;
      List current_mpm = as<List>(mpm_list(i));
      edit_lM_single(single_elem_output, current_mpm, pop, patch, year2, stage3,
        stage2, stage1, age2, eststage3, eststage2, eststage1, estage2,
        givenrate, offset, multiplier, type, type_t12);
      semifinal_output(i) = single_elem_output;
    }
    
    CharacterVector out_class = "lefkoMatList";
    semifinal_output.attr("class") = out_class;
    true_output = semifinal_output;
  } else throw Rcpp::exception(mpm_error.get_cstring(), false);
  
  return true_output;
}


