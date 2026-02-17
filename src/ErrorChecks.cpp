#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Index of Functions
// 1. List .hoffmannofstuttgart Core Engine for cond_hmpm()
// 2. List .hoffmannofstuttgart_sp Core Engine for cond_hmpm() with Sparse Matrices
// 3. List cond_hmpm Extract Conditional Ahistorical Matrices from Historical MPM
// 4. List cond_diff Extract Conditional Ahistorical Difference Matrices

//' Core Engine for cond_hmpm()
//' 
//' Creates a list of conditional ahistorical matrices in the style noted in
//' deVries and Caswell (2018).
//' 
//' @name .hoffmannofstuttgart
//' 
//' @param mainmat Historical matrix.
//' @param indices Data frame including the stages at times \emph{t}-1,
//' \emph{t}, and \emph{t}+1, asvwell as indices corresponding to elements in
//' the main historical matrix andvthe conditional matrices to be produced.
//' @param ahstages The number of stages in the stageframe.
//' @param stageframe The original stageframe for the input matrices.
//'
//' @return A list of ahistorical matrices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.hoffmannofstuttgart)]]
Rcpp::List hoffmannofstuttgart(arma::mat& mainmat, DataFrame indices,
  int ahstages, StringVector stagenames) {
  arma::uvec stage1 = indices["stage1"];
  arma::uvec stage2 = indices["stage2"];
  arma::uvec stage3 = indices["stage3"];
  arma::ivec main_index = indices["main_index"];
  arma::ivec new_index = indices["new_index"];
  
  arma::mat newmatrix(ahstages, ahstages);
  newmatrix.zeros();
  
  Rcpp::List condlist (ahstages);
  
  arma::uvec condidx = find(stage1 == 1);
  int condlength = static_cast<int>(condidx.n_elem);
  
  for (int i = 0; i < ahstages; i++) {
    newmatrix.zeros();
    
    condidx = find(stage1 == (i+1));
    condlength = static_cast<int>(condidx.n_elem);
    
    for (int j = 0; j < condlength; j++) {
      if (mainmat(main_index(condidx(j))) > 0) newmatrix(new_index(condidx(j))) = 
        newmatrix(new_index(condidx(j))) + mainmat(main_index(condidx(j)));
    }
    
    condlist(i) = newmatrix;
  }
  condlist.names() = stagenames;
  
  return (condlist);
}

//' Core Engine for cond_hmpm() with Sparse Matrices
//' 
//' Creates a list of conditional ahistorical matrices in the style noted in
//' deVries and Caswell (2018).
//' 
//' @name .hoffmannofstuttgart_sp
//' 
//' @param mainmat Historical matrix in sparse format.
//' @param indices Data frame including the stages at times \emph{t}-1,
//' \emph{t}, and \emph{t}+1, asvwell as indices corresponding to elements in
//' the main historical matrix andvthe conditional matrices to be produced.
//' @param ahstages The number of stages in the stageframe.
//' @param stageframe The original stageframe for the input matrices.
//'
//' @return A list of ahistorical matrices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.hoffmannofstuttgart_sp)]]
Rcpp::List hoffmannofstuttgart_sp(arma::sp_mat& mainmat, DataFrame indices,
  int ahstages, StringVector stagenames) {
  arma::uvec stage1 = indices["stage1"];
  arma::uvec stage2 = indices["stage2"];
  arma::uvec stage3 = indices["stage3"];
  arma::ivec main_index = indices["main_index"];
  arma::ivec new_index = indices["new_index"];
  
  arma::mat newmatrix(ahstages, ahstages);
  newmatrix.zeros();
  
  Rcpp::List condlist (ahstages);
  
  arma::uvec condidx = find(stage1 == 1);
  int condlength = static_cast<int>(condidx.n_elem);
  
  for (int i = 0; i < ahstages; i++) {
    newmatrix.zeros();
    
    condidx = find(stage1 == (i+1));
    condlength = static_cast<int>(condidx.n_elem);
    
    for (int j = 0; j < condlength; j++) {
      if (mainmat(main_index(condidx(j))) > 0) newmatrix(new_index(condidx(j))) = 
        newmatrix(new_index(condidx(j))) + mainmat(main_index(condidx(j)));
    }
    
    condlist(i) = newmatrix;
  }
  condlist.names() = stagenames;
  
  return (condlist);
}

//' Extract Conditional Ahistorical Matrices from Historical MPM
//' 
//' Function \code{cond_hmpm()} takes historical MPMs and decomposes them into 
//' ahistorical matrices conditional upon stage in time \emph{t}-1. In effect,
//' the function takes each historical matrix within a lefkoMat object, and
//' forms one ahistorical matrix for each stage in time \emph{t}-1.
//' 
//' @name cond_hmpm
//' 
//' @param hmpm A historical matrix projection model of class \code{lefkoMat}.
//' @param matchoice A character denoting whether to use A, U, or F matrices.
//' Defaults to \code{A} matrices.
//' @param err_check A logical value denoting whether to include a data frame
//' of element equivalence from the conditional matrices to the original
//' matrices. Used only for debugging purposes. Defaults to \code{FALSE}.
//' 
//' @return A \code{lefkoCondMat} object, with the following elements:
//' 
//' \item{Mcond}{A multi-level list holding the conditional A matrices derived
//' from the input \code{lefkoMat} object. The top level of the list corresponds
//' to each historical matrix in turn, and the lower level corresponds to each
//' stage in time \emph{t}-1, with individual conditional matrices named for the
//' latter.}
//' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
//' used to create historical stage pairs.}
//' \item{ahstages}{A data frame detailing the characteristics of associated
//' ahistorical stages.}
//' \item{labels}{A data frame showing the patch and year of each input full A 
//' matrix in order.}
//' \item{err_check}{An optional data frame showing the order of used element
//' indices to create conditional matrices.}
//' 
//' @examples
//' data(cypdata)
//'  
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md",
//'   "Lg", "XLg")
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
//' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
//'     "D", "XSm", "Sm", "D", "XSm", "Sm", "mat", "mat", "mat", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
//'     "SL", "SL", "D", "XSm", "Sm", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
//'     "SL", "SL", "SL", "SL", "SL", "SL", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
//'     "mat", "mat", "mat", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", "D", "XSm", "Sm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
//'     NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
//'     NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"),
//'   supplement = cypsupp3r, yearcol = "year2", patchcol = "patchid",
//'   indivcol = "individ")
//' 
//' cypcondmats <- cond_hmpm(cypmatrix3r)
//' summary(cypcondmats)
//' 
//' @export cond_hmpm
// [[Rcpp::export]]
Rcpp::List cond_hmpm(List hmpm, Nullable<CharacterVector> matchoice = R_NilValue,
  Nullable<LogicalVector> err_check = R_NilValue) {
  
  CharacterVector usedmats("a");
  int iusedmats {1};
  LogicalVector err_v = {false};
  bool err_c = false;
  
  if (matchoice.isNotNull()) {
    usedmats = matchoice;
    
    if (usedmats(0) == "A" || usedmats(0) == "a") {
      iusedmats = 1;
    } else if (usedmats(0) == "U" || usedmats(0) == "u") {
      iusedmats = 2;
    } else if (usedmats(0) == "F" || usedmats(0) == "f") {
      iusedmats = 3;
    } else {
      Rcpp::Rcout << "Choice of matrix not recognized. Using A matrices.";
      iusedmats = 1;
    }
  }
  
  if (err_check.isNotNull()) {
    err_v = err_check;
    if (err_v(0)) {
      err_c = true;
    }
  }
  
  List amats = hmpm["A"];
  List umats = hmpm["U"];
  List fmats = hmpm["F"];
  List hstages = hmpm["hstages"];
  List stageframe = hmpm["ahstages"];
  List labels = hmpm["labels"];
  int numofmats = amats.length();
  
  arma::uvec hstage1 = hstages["stage_id_1"];
  arma::uvec hstage2 = hstages["stage_id_2"];
  arma::uvec ahstages = stageframe["stage_id"];
  StringVector stagenames = stageframe["stage"];
  StringVector matnames(numofmats);
  
  for (int i = 0; i < numofmats; i++) {
    matnames(i) = i + 1;
  }
  
  int ahmpm_rows = static_cast<int>(ahstages.n_elem);
  int hmpm_rows = static_cast<int>(hstage1.n_elem);
  int hmpm_elems = 2 * ahmpm_rows * ahmpm_rows * ahmpm_rows;
  
  int format_int {0};
  if (stagenames(ahmpm_rows - 1) == "AlmostBorn") format_int = 1;
  
  arma::uvec stage1(hmpm_elems);
  arma::uvec stage2(hmpm_elems);
  arma::uvec stage3(hmpm_elems);
  arma::ivec main_index(hmpm_elems);
  arma::ivec new_index(hmpm_elems);
  stage1.zeros();
  stage2.zeros();
  stage3.zeros();
  main_index.fill(-1);
  new_index.fill(-1);
  
  int stage1_proxy {0};
  int stage2o_proxy {0};
  int stage2n_proxy {0};
  int stage3_proxy {0};
  
  int counter = 0;
  
  // Create index of matrix elements in original matrices & conditional matrices
  // Must be run on hstages inputs to handle reduced matrices properly
  if (format_int == 0) {
    // Ehrlen-format hMPMs
    for (int prior = 0; prior < hmpm_rows; prior++) {
      for (int post = 0; post < hmpm_rows; post++) {
        stage1_proxy = hstage1[prior];
        stage2o_proxy = hstage2[prior];
        stage2n_proxy = hstage1[post];
        stage3_proxy = hstage2[post];
        
        if (stage2o_proxy == stage2n_proxy) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2n_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
          
          counter++;
        }
      }
    }
  } else {
    // deVries-format hMPMs
    for (int prior = 0; prior < hmpm_rows; prior++) {
      for (int post = 0; post < hmpm_rows; post++) {
        stage1_proxy = hstage1[prior];
        stage2o_proxy = hstage2[prior];
        stage2n_proxy = hstage1[post];
        stage3_proxy = hstage2[post];
        
        if (stage2o_proxy == stage2n_proxy) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2o_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
          
          counter++;
        } else if (stage2n_proxy == ahmpm_rows && stage1_proxy < ahmpm_rows) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2o_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
              
          counter++;
        }
      }
    }
  }
  
  arma::uvec elements_to_use = find(stage1 > 0);
  stage1 = stage1.elem(elements_to_use);
  stage2 = stage2.elem(elements_to_use);
  stage3 = stage3.elem(elements_to_use);
  main_index = main_index.elem(elements_to_use);
  new_index = new_index.elem(elements_to_use);
  
  DataFrame loveontherocks = DataFrame::create(Named("stage1") = stage1,
    _["stage2"] = stage2, _["stage3"] = stage3, _["main_index"] = main_index,
    _["new_index"] = new_index);
  
  // Creates conditional matrices in list form
  List allout (numofmats);
  List mats1;
  
  if (iusedmats == 3) {
    for (int i = 0; i < numofmats; i++) {
      if (is<NumericMatrix>(fmats(i))) {
        arma::mat chosen_mat = as<arma::mat>(fmats(i));
        mats1 = hoffmannofstuttgart(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      } else { 
        arma::sp_mat chosen_mat = as<arma::sp_mat>(fmats(i));
        mats1 = hoffmannofstuttgart_sp(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      }
      allout(i) = mats1;
    }
  } else if (iusedmats == 2) {
    for (int i = 0; i < numofmats; i++) {
      if (is<NumericMatrix>(umats(i))) {
        arma::mat chosen_mat = as<arma::mat>(umats(i));
        mats1 = hoffmannofstuttgart(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      } else { 
        arma::sp_mat chosen_mat = as<arma::sp_mat>(umats(i));
        mats1 = hoffmannofstuttgart_sp(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      }
      allout(i) = mats1;
    }
  } else {
    for (int i = 0; i < numofmats; i++) {
      if (is<NumericMatrix>(amats(i))) {
        arma::mat chosen_mat = as<arma::mat>(amats(i));
        mats1 = hoffmannofstuttgart(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      } else { 
        arma::sp_mat chosen_mat = as<arma::sp_mat>(amats(i));
        mats1 = hoffmannofstuttgart_sp(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      }
      allout(i) = mats1;
    }
  }
  allout.names() = matnames;
  
  int panama_parts = 4;
  if (err_c) panama_parts = 5;
  
  List panama (panama_parts);
  
  panama(0) = allout;
  panama(1) = hstages;
  panama(2) = stageframe;
  panama(3) = labels;
  
  if (err_c) {
    panama(4) = loveontherocks;
    
    StringVector pan_1 = {"Mcond", "hstages", "ahstages", "labels", "err_check"};
    panama.names() = pan_1;
    
  } else {
    StringVector pan_2 = {"Mcond", "hstages", "ahstages", "labels"};
    panama.names() = pan_2;
  }
  
  panama.attr("class") = "lefkoCondMat";
  
  return panama;
}

//' Extract Conditional Ahistorical Difference Matrices
//' 
//' Function \code{cond_diff()} takes a set of historical difference matrices
//' resulting from function \code{\link{diff_lM}()} and decomposes them into 
//' ahistorical difference matrices conditional upon stage in time \emph{t}-1.
//' 
//' @name cond_diff
//' 
//' @param lDiff An object of class \code{lefkoDiff}.
//' @param ref Choice of mpm to use as reference. Defaults to \code{1}, which
//' means that the \code{ahstages}, \code{hstages}, and \code{labels} elements
//' for mpm1 will be used for all calculations. Only \code{1} amd \code{2} are
//' possible inputs.
//' @param matchoice A character denoting whether to use A, U, or F matrices.
//' Defaults to \code{A} matrices.
//' @param err_check A logical value denoting whether to include a data frame
//' of element equivalence from the conditional matrices to the original
//' matrices. Used only for debugging purposes. Defaults to \code{FALSE}.
//' 
//' @return A \code{lefkoCondDiff} object, with the following elements:
//' 
//' \item{Mcond}{A multi-level list holding the conditional matrices derived
//' from the input \code{lefkoDiff} object. The top level of the list
//' corresponds to each historical difference matrix in turn, and the lower
//' level corresponds to each stage in time \emph{t}-1, with individual
//' conditional matrices named for the latter.}
//' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
//' used to create historical stage pairs.}
//' \item{ahstages}{A data frame detailing the characteristics of associated
//' ahistorical stages.}
//' \item{labels}{A data frame showing the patch and year of each input full A 
//' matrix in order.}
//' \item{err_check}{An optional data frame showing the order of used element
//' indices to create conditional matrices.}
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
//'   NRasRep = TRUE)
//' 
//' seeds_per_pod <- 5000
//' 
//' cypsupp2_raw <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
//'     "XSm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
//'   givenrate = c(0.03, 0.15, 0.1, 0.1, 0.1, 0.05, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, (0.5 * seeds_per_pod),
//'     (0.5 * seeds_per_pod)),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' cypsupp3_raw <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
//'     "D", "XSm", "Sm", "D", "XSm", "Sm", "mat", "mat", "mat", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
//'     "SL", "SL", "D", "XSm", "Sm", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
//'     "SL", "SL", "SL", "SL", "SL", "SL", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
//'     "mat", "mat", "mat", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", "D", "XSm", "Sm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
//'     NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
//'     NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix2rp <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw,
//'   year = "all", patch = "all", stages = c("stage3", "stage2"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2_raw, 
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw,
//'   year = "all", stages = c("stage3", "stage2"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2_raw, 
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cypmatrix3rp <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"), 
//'   size = c("size3added", "size2added", "size1added"), supplement = cypsupp3_raw, 
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
//'   year = "all", stages = c("stage3", "stage2", "stage1"), 
//'   size = c("size3added", "size2added", "size1added"), supplement = cypsupp3_raw, 
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cypmatrix2r_3 <- hist_null(cypmatrix2r)
//' cypmatrix2r_3 <- delete_lM(cypmatrix2r_3, year = 2004)
//' diff_r <- diff_lM(cypmatrix3r, cypmatrix2r_3)
//' 
//' cypmatrix2rp_3 <- hist_null(cypmatrix2rp)
//' cypmatrix2rp_3 <- delete_lM(cypmatrix2rp_3, year = 2004)
//' diff_rp <- diff_lM(cypmatrix3rp, cypmatrix2rp_3)
//' 
//' condr1 <- cond_diff(diff_r, ref = 1)
//' condr2 <- cond_diff(diff_r, ref = 2)
//' 
//' condrp1 <- cond_diff(diff_rp, matchoice = "U", ref = 1)
//' condrp2 <- cond_diff(diff_rp, matchoice = "F", ref = 2)
//' 
//' @export cond_diff
// [[Rcpp::export]]
Rcpp::List cond_diff(List lDiff, int ref = 1,
  Nullable<CharacterVector> matchoice = R_NilValue,
  Nullable<LogicalVector> err_check = R_NilValue) {
  
  CharacterVector usedmats("a");
  int iusedmats {1};
  LogicalVector err_v = {false};
  bool err_c = false;
  
  if (ref != 1 && ref != 2) {
    throw Rcpp::exception("Option ref must equal 1 or 2.");
  }
  
  if (matchoice.isNotNull()) {
    usedmats = matchoice;
    
    if (usedmats(0) == "A" || usedmats(0) == "a") {
      iusedmats = 1;
    } else if (usedmats(0) == "U" || usedmats(0) == "u") {
      iusedmats = 2;
    } else if (usedmats(0) == "F" || usedmats(0) == "f") {
      iusedmats = 3;
    } else {
      Rcpp::Rcout << "Choice of matrix not recognized. Using A matrices.\n";
      iusedmats = 1;
    }
  }
  
  if (err_check.isNotNull()) {
    err_v = err_check;
    if (err_v(0)) {
      err_c = true;
    }
  }
  
  List amats = lDiff["A"];
  List umats = lDiff["U"];
  List fmats = lDiff["F"];
  List hstages1 = lDiff["hstages1"];
  List hstages2 = lDiff["hstages2"];
  List stageframe1 = lDiff["ahstages1"];
  List stageframe2 = lDiff["ahstages2"];
  List labels1 = lDiff["labels1"];
  List labels2 = lDiff["labels2"];
  int numofmats = amats.length();
  
  arma::uvec hstage1 = hstages1["stage_id_1"];
  arma::uvec hstage2 = hstages1["stage_id_2"];
  arma::uvec ahstages = stageframe1["stage_id"];
  StringVector stagenames = stageframe1["stage"];
  
  List hstages = hstages1;
  List stageframe = stageframe1;
  List labels = labels1;
  
  if (ref == 2) {
    hstage1 = as<arma::uvec>(hstages2["stage_id_1"]);
    hstage2 = as<arma::uvec>(hstages2["stage_id_2"]);
    ahstages = as<arma::uvec>(stageframe2["stage_id"]);
    stagenames = stageframe2["stage"];
    
    hstages = hstages2;
    stageframe = stageframe2;
    labels = labels2;
  }
  
  StringVector matnames(numofmats);
  
  for (int i = 0; i < numofmats; i++) {
    matnames(i) = i + 1;
  }
  
  int ahmpm_rows = static_cast<int>(ahstages.n_elem);
  int hmpm_rows = static_cast<int>(hstage1.n_elem);
  int hmpm_elems = 2 * ahmpm_rows * ahmpm_rows * ahmpm_rows;
  
  int format_int {0};
  if (stagenames(ahmpm_rows - 1) == "AlmostBorn") format_int = 1;
  
  arma::uvec stage1(hmpm_elems);
  arma::uvec stage2(hmpm_elems);
  arma::uvec stage3(hmpm_elems);
  arma::ivec main_index(hmpm_elems);
  arma::ivec new_index(hmpm_elems);
  stage1.zeros();
  stage2.zeros();
  stage3.zeros();
  main_index.fill(-1);
  new_index.fill(-1);
  
  int stage1_proxy {0};
  int stage2o_proxy {0};
  int stage2n_proxy {0};
  int stage3_proxy {0};
  
  int counter = 0;
  
  // Create index of matrix elements in original matrices & conditional matrices
  // Run on hstages inputs to handle reduced matrices properly
  if (format_int == 0) {
    // Ehrlen-format hMPMs
    for (int prior = 0; prior < hmpm_rows; prior++) {
      for (int post = 0; post < hmpm_rows; post++) {
        stage1_proxy = hstage1[prior];
        stage2o_proxy = hstage2[prior];
        stage2n_proxy = hstage1[post];
        stage3_proxy = hstage2[post];
        
        if (stage2o_proxy == stage2n_proxy) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2n_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
          
          counter++;
        }
      }
    }
  } else {
    // deVries-format hMPMs
    for (int prior = 0; prior < hmpm_rows; prior++) {
      for (int post = 0; post < hmpm_rows; post++) {
        stage1_proxy = hstage1[prior];
        stage2o_proxy = hstage2[prior];
        stage2n_proxy = hstage1[post];
        stage3_proxy = hstage2[post];
        
        if (stage2o_proxy == stage2n_proxy) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2o_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
          
          counter++;
        } else if (stage2n_proxy == ahmpm_rows && stage1_proxy < ahmpm_rows) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2o_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
              
          counter++;
        }
      }
    }
  }
  
  arma::uvec elements_to_use = find(stage1 > 0);
  stage1 = stage1.elem(elements_to_use);
  stage2 = stage2.elem(elements_to_use);
  stage3 = stage3.elem(elements_to_use);
  main_index = main_index.elem(elements_to_use);
  new_index = new_index.elem(elements_to_use);
  
  DataFrame loveontherocks = DataFrame::create(Named("stage1") = stage1,
    _["stage2"] = stage2, _["stage3"] = stage3, _["main_index"] = main_index,
    _["new_index"] = new_index);
  
  // Create conditional matrices in list form
  List allout (numofmats);
  List mats1;
  
  if (iusedmats == 3) {
    for (int i = 0; i < numofmats; i++) {
      if (is<NumericMatrix>(fmats(i))) {
        arma::mat chosen_mat = as<arma::mat>(fmats(i));
        mats1 = hoffmannofstuttgart(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      } else { 
        arma::sp_mat chosen_mat = as<arma::sp_mat>(fmats(i));
        mats1 = hoffmannofstuttgart_sp(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      }
      allout(i) = mats1;
    }
  } else if (iusedmats == 2) {
    for (int i = 0; i < numofmats; i++) {
      if (is<NumericMatrix>(umats(i))) {
        arma::mat chosen_mat = as<arma::mat>(umats(i));
        mats1 = hoffmannofstuttgart(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      } else { 
        arma::sp_mat chosen_mat = as<arma::sp_mat>(umats(i));
        mats1 = hoffmannofstuttgart_sp(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      }
      allout(i) = mats1;
    }
  } else {
    for (int i = 0; i < numofmats; i++) {
      if (is<NumericMatrix>(amats(i))) {
        arma::mat chosen_mat = as<arma::mat>(amats(i));
        mats1 = hoffmannofstuttgart(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      } else { 
        arma::sp_mat chosen_mat = as<arma::sp_mat>(amats(i));
        mats1 = hoffmannofstuttgart_sp(chosen_mat, loveontherocks, ahmpm_rows, stagenames);
      }
      allout(i) = mats1;
    }
  }
  allout.names() = matnames;
  
  int panama_parts = 4;
  if (err_c) panama_parts = 5;
  
  List panama (panama_parts);
  
  panama(0) = allout;
  panama(1) = hstages;
  panama(2) = stageframe;
  panama(3) = labels;
  
  if (err_c) {
    panama(4) = loveontherocks;
    
    StringVector pan_1 = {"Mcond", "hstages", "ahstages", "labels", "err_check"};
    panama.names() = pan_1;
    
  } else {
    StringVector pan_2 = {"Mcond", "hstages", "ahstages", "labels"};
    panama.names() = pan_2;
  }
  
  panama.attr("class") = "lefkoCondDiffMat";
  
  return panama;
}

