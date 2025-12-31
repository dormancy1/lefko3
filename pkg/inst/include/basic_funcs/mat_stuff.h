#ifndef LEFKOUTILS_mat_stuff_H
#define LEFKOUTILS_mat_stuff_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// Function index:
// 1. arma::uvec spmat_index  Create Element Index Meeting Condition for Sparse Matrix
// 2. arma::uvec general_index  Create General Element Index for Any lefkoMat Matrix
// 3. List decomp3  Full Eigen Analysis of a Single Dense Matrix
// 4. List decomp3sp  Full Eigen Analysis of a Single Sparse Matrix
// 5. List decomp3sp_inp  Full Eigen Analysis of a Single Sparse Matrix, with Sparse Input
// 6. arma::mat ovreplace  Re-index Projection Matrix On Basis of Overwrite Table
// 7. DataFrame sf_core  Creates Base Skeleton Stageframe
// 8. DataFrame paramnames_skeleton  Base Skeleton Data Frame for Paramnames Objects
// 9. List turbogeodiesel  Estimates Mean LefkoMat Object for Historical MPM
// 10. List geodiesel  Estimates Mean LefkoMat Object for Ahistorical MPM
// 11. int supp_decision1  Create Skeleton Plan of Expanded Supplemental Table
// 12. String supp_decision2  Decide on Stage for Each Entry in Supplemental Table
// 13. DataFrame supp_reassess  Expand Supplemental Table Given User Input
// 14. DataFrame age_expanded  Expand Supplemental Table by Age Inputs
// 15. void hst_maker  Creates hstages Data Frames
// 16. DataFrame age_maker  Creates agestages Data Frames
// 17. List theoldpizzle  Create Element Index for Matrix Estimation
// 18. List sf_reassess_internal Standardize Stageframe For MPM Analysis
// 19. List sf_leslie Create Stageframe for Population Matrix Projection Analysis


namespace LefkoMats {
  //' Create Element Index Meeting Condition for Sparse Matrix
  //' 
  //' This function takes a single sparse matrix (dgCMatrix) and creates a
  //' vector of indices meeting the condition that elements must be greater
  //' than a tolerance threshold.
  //' 
  //' @name spmat_index
  //' 
  //' @param M The sparse matrix of interest.
  //' @param tol The tolerance threshold.
  //' 
  //' @return An arma::uvec vector giving the indices of elements greater than
  //' the threshold tolerance.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::uvec spmat_index (arma::sp_mat& M, double tol) {
    int mat_dim = M.n_cols;
    arma::sp_mat::iterator it_start = M.begin();
    arma::sp_mat::iterator it_end = M.end();
    int n = std::distance(it_start, it_end);
    
    arma::uvec it_elems (n, fill::zeros);
    
    unsigned int found_count {0};
    for (sp_mat::const_iterator it = it_start; it != it_end; ++it) {
      if (*it > tol) {
        unsigned int it_row = it.row();
        unsigned int it_col = it.col();
        it_elems(found_count) = (it_col * mat_dim) + it_row;
        
        found_count++;
      }
    }
    
    arma::uvec new_index (found_count, fill::zeros);
    for (unsigned int i = 0; i < found_count; i++) {
      new_index(i) = it_elems(i);
    }
    
    return new_index;
  }

  //' Create General Element Index for Any lefkoMat Matrix
  //' 
  //' This function creates a general element index by taking a list of matrices
  //' and finding all non-zero values, which it then return C++ indices of
  //' within an arma::uvec object.
  //' 
  //' @name general_index
  //' 
  //' @param mats A list of matrices.
  //' @param tol A tolerance limit to use for assessing whether an element is
  //' equal to zero. Defaults to \code{1e-30}.
  //' @param use_tol A logical value indicating whether to use the tolerance
  //' limit set in the \code{tol} option. Defaults to \code{FALSE}.
  //' 
  //' @return An arma::uvec object listing indices of non-zero values in order.
  //' This object is essentially a union of all non-zero indices across all
  //' matrices in the list.
  //' 
  //' Note that if \code{use_tol = TRUE}, then only values greater than the
  //' tolerance limit are indexed.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::uvec general_index (Rcpp::List mats, double tol = 1e-30,
    bool use_tol = false) {
    int mat_length = mats.length();
    arma::uvec torture_chamber;
    if (!use_tol) tol = 0.0;
    
    for (int i = 0; i < mat_length; i++) {
      arma::uvec iron_maiden;
      if (is<S4>(mats(i))) {
        sp_mat rel_mat = as<arma::sp_mat>(mats(i));
        iron_maiden = LefkoMats::spmat_index(rel_mat, tol);
      } else {
        iron_maiden = find(as<arma::mat>(mats(i)) > tol);
      }
      
      if (i == 0) {
        torture_chamber = iron_maiden;
      } else {
        torture_chamber = unique(join_cols(torture_chamber, iron_maiden));
      }
    }
    
    return torture_chamber;
  }
  
  //' Full Eigen Analysis of a Single Dense Matrix
  //' 
  //' Function \code{decomp3()} returns all eigenvalues, right eigenvectors, and
  //' left eigenvectors estimated for a matrix by the \code{eig_gen}() function
  //' in the C++ Armadillo library. Works with dense matrices.
  //' 
  //' @name decomp3
  //' 
  //' @param Amat A population projection matrix of class \code{matrix}.
  //'
  //' @return This function returns all estimated eigenvalues, right
  //' eigenvectors, and left eigenvectors of a single matrix. This output is
  //' provided as a list with three parts, named appropriately.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List decomp3 (arma::mat Amat) {
    arma::cx_vec Aeigval;
    arma::cx_mat Aeigvecl;
    arma::cx_mat Aeigvecr;
    
    eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
    
    List output = List::create(Named("eigenvalues") = Aeigval,
      _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
    
    return output;
  }
  
  //' Full Eigen Analysis of a Single Sparse Matrix
  //' 
  //' Function \code{decomp3sp()} returns all eigenvalues, right eigenvectors, and
  //' left eigenvectors estimated for a matrix by the \code{eigs_gen}() function
  //' in the C++ Armadillo library. Works with sparse matrices.
  //' 
  //' @name decomp3sp
  //' 
  //' @param Amat A population projection matrix of class \code{matrix}.
  //'
  //' @return This function returns all estimated eigenvalues, right
  //' eigenvectors, and left eigenvectors of a single matrix. This output is
  //' provided as a list with three parts, named appropriately.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List decomp3sp (arma::mat Amat) {
    arma::sp_mat spAmat(Amat);
    arma::sp_mat t_spAmat = spAmat.t();
    arma::cx_vec Aeigval;
    arma::cx_vec Aeigvall;
    arma::cx_mat Aeigvecl;
    arma::cx_mat Aeigvecr;
    
    eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
    eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
    
    List output = List::create(Named("eigenvalues") = Aeigval,
      _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
    
    return output;
  }
  
  //' Full Eigen Analysis of a Single Sparse Matrix, with Sparse Input
  //' 
  //' \code{decomp3sp_inp()} returns all eigenvalues, right eigenvectors, and left
  //' eigenvectors estimated for a matrix by the \code{eigs_gen}() function
  //' in the C++ Armadillo library. Works with sparse matrices.
  //' 
  //' @name decomp3sp_inp
  //' 
  //' @param Amat A population projection matrix of class \code{matrix}.
  //'
  //' @return This function returns all estimated eigenvalues, right
  //' eigenvectors, and left eigenvectors of a single matrix. This output is
  //' provided as a list with three parts, named appropriately.
  //' 
  //' @section Notes:
  //' This function works slightly differently from function \code{decomp3sp()} in
  //' that the latter function requires a sparse matrix provided in dense format,
  //' while this function requires a sparse matrix in sparse format.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List decomp3sp_inp (arma::sp_mat spAmat) {
    arma::sp_mat t_spAmat = spAmat.t();
    arma::cx_vec Aeigval;
    arma::cx_vec Aeigvall;
    arma::cx_mat Aeigvecl;
    arma::cx_mat Aeigvecr;
    
    eigs_gen(Aeigval, Aeigvecr, spAmat, 1);
    eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
    
    List output = List::create(Named("eigenvalues") = Aeigval,
      _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
    
    return output;
  }
  
  //' Re-index Projection Matrix On Basis of Overwrite Table
  //' 
  //' Function \code{ovreplace()} takes matrix indices provided by functions
  //' \code{\link{rlefko3}()}, \code{\link{rlefko2}()}, \code{\link{flefko3}()},
  //' \code{\link{flefko2}()}, and \code{\link{aflefko2}()} and updates them with
  //' information provided in the overwrite table used as input in that function.
  //' 
  //' @name ovreplace
  //' 
  //' @param allst321 Vector containing the original element-by-element matrix
  //' index.
  //' @param idx321old Vector containing the indices of matrix elements to be
  //' updated.
  //' @param idx321new Vector containing the replacement matrix element indices.
  //' @param convtype Vector denoting survival transition (1), fecundity (2), or
  //' fecundity multiplier (3).
  //' @param eststag3 Vector of new stages in time \emph{t}+1.
  //' @param gvnrate Vector of replacement transition values.
  //' @param offset Vector of additive offset values for transitions.
  //' @param multipl Vector of fecundity multipliers.
  //' 
  //' @return A matrix with 9 columns. Column 1 is the given rate for a survival
  //' transition. Column 2 is the proxy transition to be used to estimate that
  //' transition. Column 3 is the given rate for a fecundity transition. Column
  //' 4 is the proxy transition to be used to estimate that transition. Column 5
  //' is a vector of fecundity multipliers, in cases where no given rate or
  //' proxy is to be used but fecundity is to be multiplied by some value.
  //' Column 6 is a vector of survival transition multipliers. Column 7 is a
  //' vector of fecundity transition multipliers. Column 8 is a vector of
  //' additive offsets to survival transitions. Column 9 is a vector of additive
  //' offsets to fecundity transitions.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::mat ovreplace (arma::vec allst321, arma::vec idx321old,
    arma::vec idx321new, arma::vec convtype, arma::vec eststag3, 
    arma::vec gvnrate, arma::vec offset, arma::vec multipl) {
    
    int n = static_cast<int>(idx321new.n_elem);
    
    arma::mat replacements(allst321.n_elem, 9);
    replacements.fill(-1.0);
    replacements.col(7).fill(0.);
    replacements.col(8).fill(0.);
    
    for (int i = 0; i < n; i++) {
      arma::uvec correctplace = find(allst321 == idx321old(i));
      
      int m = static_cast<int>(correctplace.n_elem); 
      
      for (int j = 0; j < m; j++) {
        if (convtype(i) == 1.) {
          if (gvnrate(i) >= 0.) {replacements(correctplace(j), 0) = gvnrate(i);}
          if (eststag3(i) != -1 && idx321new(i) >= 0) {replacements(correctplace(j), 1) = idx321new(i);}
          if (offset(i) != 0.) {replacements(correctplace(j), 7) = offset(i);}
        }
        
        if (convtype(i) == 2.) {
          if (gvnrate(i) >= 0.) {replacements(correctplace(j), 2) = gvnrate(i);}
          if (eststag3(i) != -1 && idx321new(i) >= 0) {replacements(correctplace(j), 3) = idx321new(i);}
          if (offset(i) != 0.) {replacements(correctplace(j), 8) = offset(i);}
        }
        
        if (convtype(i) == 3.) {
          replacements(correctplace(j), 4) = multipl(i);
        } else if (convtype(i) == 1.0) {
          replacements(correctplace(j), 5) = multipl(i);
        } else if (convtype(i) == 2.0) {
          replacements(correctplace(j), 6) = multipl(i);
        }
      }
    }
    
    return replacements;
  }
  
  //' Creates Base Skeleton Stageframe
  //' 
  //' Function \code{sf_core()} creates a skeleton stageframe.
  //' 
  //' @name sf_core
  //' 
  //' @param num_stages An integer denoting the number of stages to include.
  //' @param reassessed A logical value indicating whether to create a
  //' processed stageframe, or an original stageframe.
  //' @param small A logical value indicating whether to create a small
  //' stageframe composed of only five variables, or a full stageframe.
  //' 
  //' @return A data frame with empty values for all variables except
  //' \code{"stage"}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::DataFrame sf_core (int num_stages, bool reassessed = false,
    bool small = false) {
    
    if (num_stages < 1) {
      throw Rcpp::exception("Stageframe cannot be made for fewer than 1 stage.", false);
    }
    
    List output;
    
    StringVector stage (num_stages);
    NumericVector original_size (num_stages);
    IntegerVector entrystage (num_stages);
    IntegerVector repstatus (num_stages);
    NumericVector zero_double (num_stages);
    IntegerVector zero_int (num_stages);
    IntegerVector one_int (num_stages, 1);
    
    for (int i = 0; i < num_stages; i++) {
      String new_stage = "stage ";
      new_stage += String((i + 1));
      
      stage(i) = new_stage;
      
      original_size(i) = static_cast<double>(i + 1);
    }
    
    if (small) {
      IntegerVector stage_id (num_stages);
      for (int i = 0; i < num_stages; i++) {
        stage_id(i) = i+1;
      }
      
      List output_new (5);
      output_new(0) = stage_id;
      output_new(1) = stage;
      output_new(2) = original_size;
      output_new(3) = repstatus;
      output_new(4) = entrystage;
      
      output = output_new;
      Rcpp::CharacterVector varnames = {"stage_id", "stage", "original_size", 
        "repstatus", "entrystage"};
      output.attr("names") = varnames;
      
    } else {
      StringVector comments (num_stages);
      NumericVector bin_half (num_stages, 0.5);
      NumericVector bin_full (num_stages, 1.0);
      
      for (int i = 0; i < num_stages; i++) {
        String new_comment = "stage ";
        new_comment += String((i + 1));
        new_comment += " comment";
        
        comments(i) = new_comment;
      }
      
      if (reassessed) {
        IntegerVector stage_id (num_stages);
        for (int i = 0; i < num_stages; i++) {
          stage_id(i) = i+1;
        }
        
        List output_new (33);
        output_new(0) = stage_id;
        output_new(1) = stage;
        output_new(2) = original_size;
        output_new(3) = clone(zero_double);
        output_new(4) = clone(zero_double);
        output_new(5) = clone(zero_int);
        output_new(6) = clone(zero_int);
        output_new(7) = repstatus;
        output_new(8) = clone(one_int);
        output_new(9) = clone(zero_int);
        
        output_new(10) = clone(zero_int);
        output_new(11) = clone(zero_int);
        output_new(12) = entrystage;
        output_new(13) = clone(one_int);
        output_new(14) = bin_half;
        output_new(15) = clone(zero_double);
        output_new(16) = clone(zero_double);
        output_new(17) = clone(original_size);
        output_new(18) = bin_full;
        output_new(19) = clone(zero_double);
        
        output_new(20) = clone(zero_double);
        output_new(21) = clone(zero_double);
        output_new(22) = clone(zero_double);
        output_new(23) = clone(zero_double);
        output_new(24) = clone(zero_double);
        output_new(25) = clone(zero_double);
        output_new(26) = clone(zero_double);
        output_new(27) = clone(zero_double);
        output_new(28) = clone(zero_double);
        output_new(29) = clone(zero_int);
        
        output_new(30) = comments;
        output_new(31) = clone(one_int);
        output_new(32) = clone(zero_int);
        
        output = output_new;
        Rcpp::CharacterVector varnames = {"stage_id", "stage", "original_size",
          "original_size_b", "original_size_c", "min_age", "max_age",
          "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus",
          "entrystage", "indataset", "binhalfwidth_raw", "sizebin_min",
          "sizebin_max", "sizebin_center", "sizebin_width", "binhalfwidthb_raw",
          "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
          "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max",
          "sizebinc_center", "sizebinc_width", "group", "comments", "alive",
          "almostborn"};
        output.attr("names") = varnames;
        
      } else {
        List output_new (29);
        output_new(0) = stage;
        output_new(1) = original_size;
        output_new(2) = clone(zero_double);
        output_new(3) = clone(zero_double);
        output_new(4) = clone(zero_int);
        output_new(5) = clone(zero_int);
        output_new(6) = repstatus;
        output_new(7) = clone(one_int);
        output_new(8) = clone(zero_int);
        output_new(9) = clone(zero_int);
        
        output_new(10) = clone(zero_int);
        output_new(11) = clone(one_int);
        output_new(12) = bin_half;
        output_new(13) = clone(zero_double);
        output_new(14) = clone(zero_double);
        output_new(15) = clone(original_size);
        output_new(16) = bin_full;
        output_new(17) = clone(zero_double);
        output_new(18) = clone(zero_double);
        output_new(19) = clone(zero_double);
        
        output_new(20) = clone(zero_double);
        output_new(21) = clone(zero_double);
        output_new(22) = clone(zero_double);
        output_new(23) = clone(zero_double);
        output_new(24) = clone(zero_double);
        output_new(25) = clone(zero_double);
        output_new(26) = clone(zero_double);
        output_new(27) = clone(zero_int);
        output_new(28) = comments;
        
        output = output_new;
        Rcpp::CharacterVector varnames = {"stage", "size", "size_b", "size_c",
          "min_age", "max_age", "repstatus", "obsstatus", "propstatus",
          "immstatus", "matstatus", "indataset", "binhalfwidth_raw",
          "sizebin_min", "sizebin_max", "sizebin_center", "sizebin_width",
          "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
          "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw",
          "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
          "group", "comments"};
        output.attr("names") = varnames;
      }
    }
    
    CharacterVector newclasses = {"data.frame", "stageframe"};
    output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, num_stages);
    output.attr("class") = newclasses;
    
    return output;
  }
  
  //' Base Skeleton Data Frame for Paramnames Objects
  //' 
  //' Creates a simple skeleton \code{paramnames} object that can be entered as
  //' input in functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
  //' \code{\link{aflefko2}()}.
  //' 
  //' @name paramnames_skeleton
  //' 
  //' @param name_terms A logical value indicating whether to start each variable
  //' name as \code{none} if \code{FALSE}, or as the default \code{modelparams}
  //' name if \code{TRUE}. Defaults to \code{FALSE}.
  //' 
  //' @return A three column data frame, of which the first describes the
  //' parameters in reasonably plain English, the second gives the name of the
  //' parameter within the MPM generating functions, and the third is to be
  //' edited with the names of the variables as they appear in the models.
  //' 
  //' @section Notes:
  //' The third column in the resulting object should be edited with the names only
  //' of those variables actually used in vital rate modeling. This
  //' \code{paramnames} object should apply to all models used in a single MPM
  //' building exercise. So, for example, if the models used include random terms,
  //' then they should all have the same random terms. Fixed terms can vary,
  //' however.
  //' 
  //' @keywords internal
  //' @noRd
  inline DataFrame paramnames_skeleton (bool name_terms = false) {
    
    CharacterVector parameter_names = {"time t", "individual", "patch",
      "alive in time t+1", "observed in time t+1", "sizea in time t+1",
      "sizeb in time t+1", "sizec in time t+1", "reproductive status in time t+1",
      "fecundity in time t+1", "fecundity in time t", "sizea in time t",
      "sizea in time t-1", "sizeb in time t", "sizeb in time t-1",
      "sizec in time t", "sizec in time t-1", "reproductive status in time t",
      "reproductive status in time t-1", "maturity status in time t+1",
      "maturity status in time t", "age in time t", "density in time t",
      "individual covariate a in time t", "individual covariate a in time t-1",
      "individual covariate b in time t", "individual covariate b in time t-1",
      "individual covariate c in time t", "individual covariate c in time t-1",
      "stage group in time t", "stage group in time t-1",
      "annual covariate a in time t", "annual covariate a in time t-1",
      "annual covariate b in time t", "annual covariate b in time t-1",
      "annual covariate c in time t", "annual covariate c in time t-1"};
    
    CharacterVector mainparams = {"year2", "individ", "patch", "surv3", "obs3",
      "size3", "sizeb3", "sizec3", "repst3", "fec3", "fec2", "size2", "size1",
      "sizeb2", "sizeb1", "sizec2", "sizec1", "repst2", "repst1", "matst3",
      "matst2", "age", "density", "indcova2", "indcova1", "indcovb2", "indcovb1",
      "indcovc2", "indcovc1", "group2", "group1", "annucova2", "annucova1",
      "annucovb2", "annucovb1", "annucovc2", "annucovc1"};
    
    CharacterVector modelparams = {"none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none"};
    
    CharacterVector modelparams_alt = {"year2", "individ", "patchid", "alive3",
      "obsstatus3", "sizea3", "sizeb3", "sizec3", "repstatus3", "feca3",
      "feca2", "sizea2", "sizea1", "sizeb2", "sizeb1", "sizec2", "sizec1",
      "repstatus2", "repstatus1", "matstatus3", "matstatus2", "obsage",
      "density", "indcova2", "indcova1", "indcovb2", "indcovb1", "indcovc2",
      "indcovc1", "group2", "group1", "annucova2", "annucova1", "annucovb2",
      "annucovb1", "annucovc2", "annucovc1"};
    
    if (name_terms) modelparams = modelparams_alt;
    
    DataFrame output = DataFrame::create(_["parameter_names"] = parameter_names,
      _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
    return output;
  }
  
  //' Estimates Mean LefkoMat Object for Historical MPM
  //' 
  //' Function \code{turbogeodiesel()} estimates mean historical population
  //' projection matrices, treating the mean as element-wise arithmetic.
  //' 
  //' @name turbogeodiesel
  //' 
  //' @param loy A data frame denoting the population, patch, and occasion
  //' designation for each matrix. Includes a total of 9 variables.
  //' @param Umats A matrix with all U matrices turned into columns.
  //' @param Fmats A matrix with all F matrices turned into columns.
  //' @param hstages This is the \code{hstages} object held by \code{mats}.
  //' @param agestages This is the \code{agestages} object held by \code{mats}.
  //' @param stages This is the core stageframe held by \code{mats}, equivalent to
  //' \code{ahstages}.
  //' @param patchmats A logical value stating whether to estimate patch-level
  //' means.
  //' @param popmats A logical value stating whether to estimate population-level
  //' means.
  //' @param mat_input A logical value indicating whether the input matrix class
  //' is a standard \code{NumericMatrix}.
  //' @param sparse_switch A number denoting whether to force simple matrices to
  //' be projected as sparse (\code{1}) or not (\code{0}).
  //' 
  //' @return A list using the structure of a lefkoMat object.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List turbogeodiesel (DataFrame& loy, List Umats, List Fmats,
    DataFrame hstages, DataFrame agestages, DataFrame stages, bool patchmats,
    bool popmats, bool mat_input, int sparse_switch) {
    
    StringVector pops = as<StringVector>(loy["pop"]);
    arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
    StringVector patches = as<StringVector>(loy["patch"]);
    arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
    arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
    arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
    arma::uvec uniquepops = unique(pop_num);
    arma::uvec uniquepoppatches = unique(poppatchc);
    int loydim = pops.length();
    int numofpops = static_cast<int>(uniquepops.n_elem);
    int numofpatches = static_cast<int>(uniquepoppatches.n_elem);
    
    if (numofpatches == 1) popmats = false;
    
    StringVector poporderlong(loydim);
    arma::uvec poporderlong_num(loydim);
    StringVector patchorderlong(loydim);
    arma::uvec annmatriceslong(loydim);
    arma::uvec meanassign(loydim);
    poporderlong_num.zeros();
    annmatriceslong.zeros();
    meanassign.zeros();
    
    pop_num = pop_num + 1;
    poppatchc = poppatchc + 1;
    
    poporderlong(0) = pops(0);
    poporderlong_num(0) = pop_num(0);
    patchorderlong(0) = patches(0);
    annmatriceslong(0) = 1;
    meanassign(0) = 1;
    
    int counter {0};
    
    StringVector uniquepops_str(numofpops);
    uniquepops_str(0) = pops(0);
    int popcounter {0};
    
    // Needed matrix number and order determination
    if (loydim > 1) {
      for (int i = 1; i < loydim; i++) {
        if (poppatchc(i) != poppatchc(i-1)) {
          counter++;
          poporderlong(counter) = pops(i);
          poporderlong_num(counter) = pop_num(i);
          patchorderlong(counter) = patches(i);
          annmatriceslong(counter) = 1;
          meanassign(i) = meanassign(i-1) + 1;
          
          if (pop_num(i) != pop_num(i-1)) {
            popcounter += 1;
            uniquepops_str(popcounter) = pops(i);
          }
          
        } else {
          annmatriceslong(counter) = annmatriceslong(counter) + 1;
          meanassign(i) = meanassign(i-1);
        }
      }
    }
    
    arma::uvec toestimate = find(poporderlong_num);
    int popcount = static_cast<int>(toestimate.n_elem);
    
    int totalmatrices = static_cast<int>(toestimate.n_elem) + numofpops;
    
    if (patchmats && !popmats) {
      totalmatrices = static_cast<int>(toestimate.n_elem);
    } else if (!patchmats && popmats) {
      totalmatrices = numofpops;
    }
    
    arma::uvec poporder = poporderlong_num.elem(toestimate);
    arma::uvec patchorder = poppatchc.elem(toestimate);
    arma::uvec annmatrices = annmatriceslong.elem(toestimate);
    
    StringVector poporder_str(popcount);
    StringVector patchorder_str(popcount);
    
    for (int i = 0; i < popcount; i++) {
      poporder_str(i) = pops(toestimate(i));
      patchorder_str(i) = patches(toestimate(i));
    }
    
    // Predicts which elements will be targeted for arithmetic mean estimation
    int format_int {0};
    arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
    StringVector stagenames = as<StringVector>(stages["stage"]);
    int numstages = static_cast<int>(astages.n_elem);
    
    if (stagenames(numstages - 1) == "AlmostBorn") format_int = 1;
    
    arma::uvec hstage3in = as<arma::uvec>(hstages["stage_id_2"]);
    arma::uvec hstage2nin = as<arma::uvec>(hstages["stage_id_1"]);
    int numhstages = static_cast<int>(hstage3in.n_elem);
    
    int predictedsize = 2 * numstages * numstages * numstages;
    
    arma::uvec hsindexl(predictedsize);
    hsindexl.zeros();
    
    counter = 0;
    
    if (format_int == 0) { // Ehrlen format
      for (int i1 = 0; i1 < numhstages; i1++) {
        for (int i2 = 0; i2 < numhstages; i2++) {
          if (hstage3in(i1) == hstage2nin(i2)) {
            hsindexl(counter) = (i1 * numhstages) + i2;
            counter++;
          }
        }
      }
    } else { // deVries format
      for (int i1 = 0; i1 < numhstages; i1++) {
        for (int i2 = 0; i2 < numhstages; i2++) {
          if (hstage3in(i1) == hstage2nin(i2)) {
            hsindexl(counter) = (i1 * numhstages) + i2;
            counter++;
            
          } else if (static_cast<int>(hstage2nin(i2)) == numstages || 
              static_cast<int>(hstage3in(i1)) == numstages) {
            hsindexl(counter) = (i1 * numhstages) + i2;
            counter++;
          }
        }
      }
    }
    
    arma::uvec hsgood = find(hsindexl);
    arma::uvec hsindex = hsindexl.elem(hsgood);
    arma::uvec zerovec(1, fill::zeros);
    arma::uvec allindices = join_cols(zerovec, hsindex);
    
    // Build U & F matrices of element-wise arithmetic means
    // Each column holds predicted non-zero elements of each mean matrix
    // Each matrix is a column vector within the overall matrix
    int core_elem = counter;
    
    arma::mat umatvec(core_elem, totalmatrices, fill::zeros);
    arma::mat fmatvec(core_elem, totalmatrices, fill::zeros);
    
    int patchchoice {0};
    int popchoice {0};
    
    pop_num = pop_num - 1;
    poppatchc = poppatchc - 1;
    
    for (int i = 0; i < loydim; i ++) {
      if (patchmats) {
        patchchoice = poppatchc(i);
        
        if (mat_input) {
          arma::mat Umats_invaded = as<arma::mat>(Umats(i));
          arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
          
          umatvec.col(patchchoice) = umatvec.col(patchchoice) +
            (Umats_invaded.elem(allindices) / yearsinpatch(i));
          fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
            (Fmats_invaded.elem(allindices) / yearsinpatch(i));
        } else {
          arma::sp_mat Umats_invaded = as<arma::sp_mat>(Umats(i));
          arma::sp_mat Fmats_invaded = as<arma::sp_mat>(Fmats(i));
          
          arma::vec Umats_invaded_allindices (allindices.n_elem, fill::zeros);
          arma::vec Fmats_invaded_allindices (allindices.n_elem, fill::zeros);
          
          for (int j = 0; j < static_cast<int>(allindices.n_elem); j++) {
            Umats_invaded_allindices(j) = Umats_invaded(allindices(j));
            Fmats_invaded_allindices(j) = Fmats_invaded(allindices(j));
          }
        
          umatvec.col(patchchoice) = umatvec.col(patchchoice) +
            (Umats_invaded_allindices / yearsinpatch(i));
          fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
            (Fmats_invaded_allindices / yearsinpatch(i));
        }
        
        if (popmats) {
          if (patchmats) {
            popchoice = numofpatches + pop_num(i);
          } else {
            popchoice = pop_num(i);
          }
          
          if (mat_input) {
            arma::mat Umats_invaded = as<arma::mat>(Umats(i));
            arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
            
            umatvec.col(popchoice) = umatvec.col(popchoice) +
              (Umats_invaded.elem(allindices) / (yearsinpatch(i) * patchesinpop(i)));
            fmatvec.col(popchoice) = fmatvec.col(popchoice) +
              (Fmats_invaded.elem(allindices) / (yearsinpatch(i) * patchesinpop(i)));
          } else {
            arma::sp_mat Umats_invaded = as<arma::sp_mat>(Umats(i));
            arma::sp_mat Fmats_invaded = as<arma::sp_mat>(Fmats(i));
            
            arma::vec Umats_invaded_allindices (allindices.n_elem, fill::zeros);
            arma::vec Fmats_invaded_allindices (allindices.n_elem, fill::zeros);
            
            for (int j = 0; j < static_cast<int>(allindices.n_elem); j++) {
              Umats_invaded_allindices(j) = Umats_invaded(allindices(j));
              Fmats_invaded_allindices(j) = Fmats_invaded(allindices(j));
            }
          
            umatvec.col(popchoice) = umatvec.col(popchoice) +
              (Umats_invaded_allindices / (yearsinpatch(i) * patchesinpop(i)));
            fmatvec.col(popchoice) = fmatvec.col(popchoice) +
              (Fmats_invaded_allindices / (yearsinpatch(i) * patchesinpop(i)));
          }
        }
      }
    }
    
    // New labels data frame
    int cheatsheetlength {1};
    
    if (patchmats && !popmats) {
      cheatsheetlength = numofpatches;
      
    } else if (!patchmats && popmats) {
      cheatsheetlength = numofpops;
      
    } else if (popmats && patchmats) {
      cheatsheetlength = numofpops + numofpatches;
    }
    
    StringVector poporder_redone(cheatsheetlength);
    StringVector patchorder_redone(cheatsheetlength);
    
    int new_labels_counter {0};
    
    if (patchmats) {
      for (int i = 0; i < numofpatches; i++) {
        poporder_redone(i) = poporderlong(i);
        patchorder_redone(i) = patchorderlong(i);
        
        new_labels_counter++;
      }
    }
    
    if (popmats) {
      for (int i = 0; i < numofpops; i++) {
        poporder_redone(new_labels_counter) = uniquepops_str(i);
        patchorder_redone(new_labels_counter) = "0";
        
        new_labels_counter++;
      }
    }
    
    DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone,
      _["patch"] = patchorder_redone);
    
    // Create the main list objects holding matrices
    List U (totalmatrices);
    List F (totalmatrices);
    List A (totalmatrices);
    int totalutrans {0};
    int totalftrans {0};
    
    if (mat_input && sparse_switch == 0) {
      for (int i = 0; i < totalmatrices; i++) {
        arma::mat umat_base(numhstages, numhstages, fill::zeros);
        arma::mat fmat_base(numhstages, numhstages, fill::zeros);
        
        umat_base.elem(allindices) = umatvec.col(i);
        fmat_base.elem(allindices) = fmatvec.col(i);
        arma::mat amat_base = umat_base + fmat_base;
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = static_cast<int>(utrans.n_elem);
      totalftrans = static_cast<int>(ftrans.n_elem);
      
    } else {
      for (int i = 0; i < totalmatrices; i++) {
        arma::sp_mat umat_base(numhstages, numhstages);
        arma::sp_mat fmat_base(numhstages, numhstages);
        
        for (int j = 0; j < static_cast<int>(allindices.n_elem); j++) {
          umat_base(allindices(j)) = umatvec(j, i);
          fmat_base(allindices(j)) = fmatvec(j, i);
        }
        arma::sp_mat amat_base = umat_base + fmat_base;
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = static_cast<int>(utrans.n_elem);
      totalftrans = static_cast<int>(ftrans.n_elem);
    }
    
    NumericVector matrixqc(3);
    matrixqc(0) = totalutrans; // summed number of non-zero u transitions
    matrixqc(1) = totalftrans; // summed number of non-zero f transitions
    matrixqc(2) = totalmatrices;
    
    // Final output
    List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F,
      _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
    
    return output;
  }
  
  //' Estimates Mean LefkoMat Object for Ahistorical MPM
  //' 
  //' Function \code{geodiesel()} estimates mean ahistorical population
  //' projection matrices, treating the mean as element-wise arithmetic. The
  //' function can handle both normal ahistorical MPMs and age x stage ahistorical
  //' MPMs.
  //' 
  //' @name geodiesel
  //' 
  //' @param loy A data frame denoting the population, patch, and occasion
  //' designation of each matrix. Includes a total of 9 variables.
  //' @param Umats A matrix with all U matrices turned into columns.
  //' @param Fmats A matrix with all F matrices turned into columns.
  //' @param agestages This is the \code{agestages} object held by \code{mats}.
  //' @param stages This is the core stageframe held by \code{mats}, equivalent to
  //' \code{ahstages}.
  //' @param patchmats A logical value stating whether to estimate patch-level
  //' means.
  //' @param popmats A logical value stating whether to estimate population-level
  //' means.
  //' @param mat_input A logical value indicating whether the input matrix class
  //' is a standard \code{NumericMatrix}.
  //' @param sparse_switch A number denoting whether to force simple matrices to
  //' be projected as sparse (\code{1}) or not (\code{0}).
  //' 
  //' @return A list using the structure of a LefkoMat object.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List geodiesel (DataFrame& loy, List Umats, List Fmats,
    DataFrame agestages, DataFrame stages, bool patchmats, bool popmats,
    bool mat_input, int sparse_switch) {
    
    StringVector pops = as<StringVector>(loy["pop"]);
    arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
    StringVector patches = as<StringVector>(loy["patch"]);
    arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
    arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
    arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
    arma::uvec uniquepops = unique(pop_num);
    arma::uvec uniquepoppatches = unique(poppatchc);
    int loydim = pops.length();
    int numofpops = static_cast<int>(uniquepops.n_elem);
    int numofpatches = static_cast<int>(uniquepoppatches.n_elem);
    
    if (numofpatches == 1) popmats = false;
    
    StringVector poporderlong(loydim);
    arma::uvec poporderlong_num(loydim);
    StringVector patchorderlong(loydim);
    arma::uvec annmatriceslong(loydim);
    arma::uvec meanassign(loydim);
    poporderlong_num.zeros();
    annmatriceslong.zeros();
    meanassign.zeros();
    
    pop_num = pop_num + 1;
    poppatchc = poppatchc + 1;
    
    poporderlong(0) = pops(0);
    poporderlong_num(0) = pop_num(0);
    patchorderlong(0) = patches(0);
    annmatriceslong(0) = 1;
    meanassign(0) = 1;
    
    int counter {0};
    
    StringVector uniquepops_str(numofpops);
    uniquepops_str(0) = pops(0);
    int popcounter {0};
    
    // Needed matrix number and order determination
    if (loydim > 1) {
      for (int i = 1; i < loydim; i++) {
        if (poppatchc(i) != poppatchc(i-1)) {
          counter++;
          poporderlong(counter) = pops(i);
          poporderlong_num(counter) = pop_num(i);
          patchorderlong(counter) = patches(i);
          annmatriceslong(counter) = 1;
          meanassign(i) = meanassign(i-1) + 1;
          
          if (pop_num(i) != pop_num(i-1)) {
            popcounter += 1;
            uniquepops_str(popcounter) = pops(i);
          }
          
        } else {
          annmatriceslong(counter) = annmatriceslong(counter) + 1;
          meanassign(i) = meanassign(i-1);
        }
      }
    }
    
    arma::uvec toestimate = find(poporderlong_num);
    int popcount = static_cast<int>(toestimate.n_elem);
    
    int totalmatrices = static_cast<int>(toestimate.n_elem) + numofpops;
    
    if (patchmats && !popmats) {
      totalmatrices = static_cast<int>(toestimate.n_elem);
    } else if (!patchmats && popmats) {
      totalmatrices = numofpops;
    }
    
    arma::uvec poporder = poporderlong_num.elem(toestimate);
    arma::uvec patchorder = poppatchc.elem(toestimate);
    arma::uvec annmatrices = annmatriceslong.elem(toestimate);
    
    StringVector poporder_str(popcount);
    StringVector patchorder_str(popcount);
    
    for (int i = 0; i < popcount; i++) {
      poporder_str(i) = pops(toestimate(i));
      patchorder_str(i) = patches(toestimate(i));
    }
    
    // Predicts which elements will be targeted for arithmetic mean estimation
    arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
    int initialstages = static_cast<int>(astages.n_elem);
    
    // Test for the presence of ages, and determine matrix dimensions
    int colsused {0};
    
    if (mat_input && sparse_switch == 0) {
      arma::mat initUmat = as<arma::mat>(Umats(0));
      colsused = initUmat.n_cols;
    } else if (mat_input) {
      arma::sp_mat initUmat(as<arma::mat>(Umats(0)));
      colsused = initUmat.n_cols;
    } else {
      arma::sp_mat initUmat = as<arma::sp_mat>(Umats(0));
      colsused = initUmat.n_cols;
    }
    
    int agemultiplier = colsused / initialstages;
    
    int numstages = static_cast<int>(astages.n_elem) * agemultiplier;
    
    // Build U & F matrices of element-wise arithmetic means
    // Each column holds predicted non-zero elements of each mean matrix
    // Each matrix is a column vector within the overall matrix
    int core_elem = numstages * numstages;
    
    arma::mat umatvec;
    arma::mat fmatvec;
    arma::mat amatvec;
    arma::sp_mat umatvec_sp;
    arma::sp_mat fmatvec_sp;
    arma::sp_mat amatvec_sp;
    
    if (mat_input && sparse_switch == 0) {
      arma::mat umatvec_(core_elem, totalmatrices);
      arma::mat fmatvec_(core_elem, totalmatrices);
      umatvec_.zeros();
      fmatvec_.zeros();
      
      int patchchoice {0};
      int popchoice {0};
      
      pop_num = pop_num - 1;
      poppatchc = poppatchc - 1;
      
      for (int i = 0; i < loydim; i ++) {
        if (patchmats) {
          patchchoice = poppatchc(i);
          
          arma::mat Umats_invaded = as<arma::mat>(Umats(i));
          arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
          
          umatvec_.col(patchchoice) = umatvec_.col(patchchoice) +
            (arma::vectorise(Umats_invaded) / yearsinpatch(i));
          fmatvec_.col(patchchoice) = fmatvec_.col(patchchoice) +
            (arma::vectorise(Fmats_invaded) / yearsinpatch(i));
        }
        
        if (popmats) {
          if (patchmats) {
            popchoice = numofpatches + pop_num(i);
          } else {
            popchoice = pop_num(i);
          }
          
          arma::mat Umats_invaded = as<arma::mat>(Umats(i));
          arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
          
          umatvec_.col(popchoice) = umatvec_.col(popchoice) +
            (arma::vectorise(Umats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
          fmatvec_.col(popchoice) = fmatvec_.col(popchoice) +
            (arma::vectorise(Fmats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
        }
      }
      umatvec = umatvec_;
      fmatvec = fmatvec_;
      amatvec = umatvec + fmatvec;
      
    } else {
      arma::sp_mat umatvec_sp_(core_elem, totalmatrices);
      arma::sp_mat fmatvec_sp_(core_elem, totalmatrices);
      
      int patchchoice {0};
      int popchoice {0};
      
      pop_num = pop_num - 1;
      poppatchc = poppatchc - 1;
      
      for (int i = 0; i < loydim; i ++) {
        if (patchmats) {
          patchchoice = poppatchc(i);
          
          arma::sp_mat Umats_invaded;
          arma::sp_mat Fmats_invaded;
          
          if (mat_input) { 
            arma::sp_mat Umats_invaded_(as<arma::mat>(Umats(i)));
            arma::sp_mat Fmats_invaded_(as<arma::mat>(Fmats(i)));
            Umats_invaded = Umats_invaded_;
            Fmats_invaded = Fmats_invaded_;
          } else { 
            Umats_invaded = as<arma::sp_mat>(Umats(i));
            Fmats_invaded = as<arma::sp_mat>(Fmats(i));
          }
          
          umatvec_sp_.col(patchchoice) = umatvec_sp_.col(patchchoice) +
            (arma::vectorise(Umats_invaded) / yearsinpatch(i));
          fmatvec_sp_.col(patchchoice) = fmatvec_sp_.col(patchchoice) +
            (arma::vectorise(Fmats_invaded) / yearsinpatch(i));
        }
        
        if (popmats) {
          if (patchmats) {
            popchoice = numofpatches + pop_num(i);
          } else {
            popchoice = pop_num(i);
          }
          
          arma::sp_mat Umats_invaded;
          arma::sp_mat Fmats_invaded;
          
          if (mat_input) {
            arma::sp_mat Umats_invaded_(as<arma::mat>(Umats(i)));
            arma::sp_mat Fmats_invaded_(as<arma::mat>(Fmats(i)));
            Umats_invaded = Umats_invaded_;
            Fmats_invaded = Fmats_invaded_;
          } else {
            Umats_invaded = as<arma::sp_mat>(Umats(i));
            Fmats_invaded = as<arma::sp_mat>(Fmats(i));
          }
          
          umatvec_sp_.col(popchoice) = umatvec_sp_.col(popchoice) +
            (arma::vectorise(Umats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
          fmatvec_sp_.col(popchoice) = fmatvec_sp_.col(popchoice) +
            (arma::vectorise(Fmats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
        }
      }
      
      umatvec_sp = umatvec_sp_;
      fmatvec_sp = fmatvec_sp_;
      amatvec_sp = umatvec_sp + fmatvec_sp;
    }
    
    // New labels data frame
    int cheatsheetlength {1};
    
    if (patchmats && !popmats) {
      cheatsheetlength = numofpatches;
      
    } else if (!patchmats && popmats) {
      cheatsheetlength = numofpops;
      
    } else if (popmats && patchmats) {
      cheatsheetlength = numofpops + numofpatches;
    }
    StringVector poporder_redone(cheatsheetlength);
    StringVector patchorder_redone(cheatsheetlength);
    
    int new_labels_counter {0};
    
    if (patchmats) {
      for (int i = 0; i < numofpatches; i++) {
        poporder_redone(i) = poporderlong(i);
        patchorder_redone(i) = patchorderlong(i);
        
        new_labels_counter++;
      }
    }
    
    if (popmats) {
      for (int i = 0; i < numofpops; i++) {
        poporder_redone(new_labels_counter) = uniquepops_str(i);
        patchorder_redone(new_labels_counter) = "0";
        
        new_labels_counter++;
      }
    }
    
    DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, 
      _["patch"] = patchorder_redone);
    
    // Create the main list objects to hold the matrices
    List U (totalmatrices);
    List F (totalmatrices);
    List A (totalmatrices);
    int totalutrans {0};
    int totalftrans {0};
    
    if (mat_input && sparse_switch == 0) {
      for (int i = 0; i < totalmatrices; i++) {
        arma::mat umat_base = umatvec.col(i);
        arma::mat fmat_base = fmatvec.col(i);
        arma::mat amat_base = amatvec.col(i);
        
        umat_base.reshape(numstages, numstages);
        fmat_base.reshape(numstages, numstages);
        amat_base.reshape(numstages, numstages);
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = static_cast<int>(utrans.n_elem);
      totalftrans = static_cast<int>(ftrans.n_elem);
      
    } else {
      for (int i = 0; i < totalmatrices; i++) {
        arma::sp_mat umat_base = umatvec_sp.col(i);
        arma::sp_mat fmat_base = fmatvec_sp.col(i);
        arma::sp_mat amat_base = amatvec_sp.col(i);
        
        umat_base.reshape(numstages, numstages);
        fmat_base.reshape(numstages, numstages);
        amat_base.reshape(numstages, numstages);
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = static_cast<int>(utrans.n_elem);
      totalftrans = static_cast<int>(ftrans.n_elem);
    }
    
    NumericVector matrixqc(3);
    matrixqc(0) = totalutrans; // summed number of U transitions
    matrixqc(1) = totalftrans; // summed number of F transitions
    matrixqc(2) = totalmatrices;
    
    // Final output
    List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F,
      _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
    
    return output;
  }
  
  //' Create Skeleton Plan of Expanded Supplemental Table
  //' 
  //' This function is used to take a supplemental table input and expand it
  //' given shorthand codes that users may have used.
  //' 
  //' @name supp_decision1
  //' 
  //' @param base_check The string input from the user's supplemental table.
  //' @param np_s Number of propagule stages.
  //' @param np0_s Number of non-propagule stages.
  //' @param ni_s Number of immature stages.
  //' @param nm_s Number of mature stages.
  //' @param nr_s Number of reproductive stages.
  //' @param nmr0_s Number of mature, non-reproductive stages.
  //' @param no_s Number of observable stages.
  //' @param no0_s Number of unobservable stages.
  //' @param a_s Total number of stages.
  //' @param no_groups Total number of stage groups.
  //' @param newgroupvec An integer vector giving the group number for each stage.
  //' @param group_text A string vector giving the names of all stage groups.
  //' 
  //' @return This function returns a single integer value corresponding to the
  //' number of stages to include for the given code in the expanded supplemental
  //' table.
  //' 
  //' @keywords internal
  //' @noRd
  inline int supp_decision1 (std::string base_check, int np_s, int np0_s,
    int ni_s, int nm_s, int nr_s, int nmr0_s, int no_s, int no0_s, int a_s,
    int no_groups, arma::ivec newgroupvec, StringVector group_text) {
    
    int decided {0};
    int no_current_group {0};
    
    if (base_check == "prop") {
      decided = np_s;
    } else if (base_check == "npr") {
      decided = np0_s;
    } else if (base_check == "immat") {
      decided = ni_s;
    } else if (base_check == "mat") {
      decided = nm_s;
    } else if (base_check == "rep") {
      decided = nr_s;
    } else if (base_check == "nrep") {
      decided = nmr0_s;
    } else if (base_check == "obs") {
      decided = no_s;
    } else if (base_check == "nobs") {
      decided = no0_s;
    } else if (base_check == "all") {
      decided = a_s;
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (base_check == as<std::string>(group_text(j))) {
          arma::uvec current_group = find(newgroupvec == j);
          no_current_group = static_cast<int>(current_group.n_elem);
          
          decided = no_current_group;
        }
      }
    }
    if (decided == 0) decided = 1;
    
    return decided;
  }
  
  //' Decide on Stage for Each Entry in Supplemental Table
  //' 
  //' This function is used to take a supplemental table input and expand it
  //' given shorthand codes that users may have used. It provides the actual stage
  //' designations given the numbers of stages decided on via function
  //' \code{supp_decision1()}.
  //' 
  //' @name supp_decision2
  //' 
  //' @param base_check The string input from the user's supplemental table.
  //' @param newprop_stages Integer designations of all propagule stages.
  //' @param newprop0_stages Integer designations of all non-propagule stages.
  //' @param newimm_stages Integer designations of all immature stages.
  //' @param newmat_stages Integer designations of all mature stages.
  //' @param newrep_stages Integer designations of all reproductive stages.
  //' @param newmat_rep0_stages Integer designations of all mature,
  //' non-reproductive stages.
  //' @param newobs_stages Integer designations of all observable stages.
  //' @param newobs0_stages Integer designations of all unobservable stages.
  //' @param all_stages Integer designations of all stages.
  //' @param no_groups Total number of stage groups.
  //' @param newgroupvec An integer vector giving the group number for each stage.
  //' @param group_text A string vector giving the names of all stage groups.
  //' @param stagevec A string vector with the names of all stages.
  //' @param counter An integer giving the correct position within vectors for
  //' stage designations.
  //' @param group_check An integer used to decide whether to continue checking
  //' group deisngation.
  //' @param group_ratchet An integer giving a cut-off number for group
  //' designation, input as a pointer to allow processing across supplemental
  //' table lines.
  //' @param group_baseline An integer used in calculating the correct stage
  //' designation given a group code.
  //' @param prevl An integer allowing the counter value in the preceding time
  //' step to be kept in memory.
  //' 
  //' @return This function returns a single string value corresponding to the
  //' correct stage to include for the given code in the input supplement table.
  //' 
  //' @keywords internal
  //' @noRd
  inline String supp_decision2 (std::string base_check,
    arma::uvec newprop_stages, arma::uvec newprop0_stages,
    arma::uvec newimm_stages, arma::uvec newmat_stages,
    arma::uvec newrep_stages, arma::uvec newmat_rep0_stages,
    arma::uvec newobs_stages, arma::uvec newobs0_stages, arma::uvec all_stages,
    int no_groups, arma::ivec newgroupvec, StringVector group_text,
    StringVector stagevec, int counter, int group_check, int& group_ratchet,
    int& group_baseline, int& prevl) {
    
    String decided;
    
    if (base_check == "prop") {
      decided = stagevec(newprop_stages(counter));
    } else if (base_check == "npr") {
      decided = stagevec(newprop0_stages(counter));
    } else if (base_check == "immat") {
      decided = stagevec(newimm_stages(counter));
    } else if (base_check == "mat") {
      decided = stagevec(newmat_stages(counter));
    } else if (base_check == "rep") {
      decided = stagevec(newrep_stages(counter));
    } else if (base_check == "nrep") {
      decided = stagevec(newmat_rep0_stages(counter));
    } else if (base_check == "obs") {
      decided = stagevec(newobs_stages(counter));
    } else if (base_check == "nobs") {
      decided = stagevec(newobs0_stages(counter));
    } else if (base_check == "all") {
      decided = stagevec(all_stages(counter));
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (base_check == as<std::string>(group_text(j))) {
          if (counter == 0) group_ratchet = 0;
          if (counter != prevl && counter != 0) group_ratchet += 1;
          
          group_check = 1;
          arma::uvec current_group = find(newgroupvec == j);
          int current_group_length = static_cast<int>(current_group.n_elem);
          if (group_ratchet > (current_group_length - 1)) {
            group_ratchet = 0;
          }
          
          if (group_ratchet == 0) {
            group_baseline = counter;
          }
          
          decided = stagevec(current_group(counter - group_baseline));
          
          prevl = counter;
        }
      }
      
      if (group_check == 0) {
        decided = base_check;
      }
      
      group_check = 0;
    }
    if (decided == 0) decided = 1;
    
    return decided;
  }
  
  //' Expand Supplemental Table Given User Input
  //' 
  //' The function takes a supplemental table as input and produces an edited and
  //' expanded version for calculation.
  //' 
  //' @name supp_reassess
  //' 
  //' @param stageframe The stageframe used for assessment.
  //' @param historical A logical value indicating whether the MPM should be
  //' historical (\code{TRUE}) or not (\code{FALSE}).
  //' @param supplement The user's supplemental table, if used.
  //' @param overwrite The user's overwrite table, if used.
  //' 
  //' @return This function returns a new data frame that acts as the expanded
  //' supplemental table without shorthand codes. It uses only stage
  //' designations from the stageframe used.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::DataFrame supp_reassess (Rcpp::DataFrame stageframe,
    bool historical, Nullable<DataFrame> supplement = R_NilValue,
    Nullable<DataFrame> overwrite = R_NilValue) {
    
    StringVector stagevec = as<StringVector>(stageframe["stage"]);
    
    arma::ivec groupvec = as<arma::ivec>(stageframe["group"]);
    int stageframe_length {static_cast<int>(stagevec.length())};
    IntegerVector stage_id = seq(1, stageframe_length);
    
    // Identify all groups
    arma::ivec all_groups = unique(groupvec);
    int no_groups {static_cast<int>(all_groups.n_elem)};
    StringVector group_text(no_groups);
    
    for (int i = 0; i < no_groups; i++) {
      group_text(i) = "group";
      group_text(i) += std::to_string(all_groups(i));
    }
    
    StringVector stage3_supp;
    StringVector stage2_supp;
    StringVector stage1_supp;
    IntegerVector age2_supp;
    StringVector eststage3_supp;
    StringVector eststage2_supp;
    StringVector eststage1_supp;
    IntegerVector estage2_supp;
    NumericVector givenrate_supp;
    NumericVector offset_supp;
    NumericVector multiplier_supp;
    IntegerVector convtype_supp;
    IntegerVector convtype_t12_supp;
    
    StringVector pop_supp;
    StringVector patch_supp;
    StringVector year2_supp;
    
    int supp_rows {0};
    
    if (supplement.isNotNull()) {
      Rcpp::DataFrame supplement_true(supplement);
      
      stage3_supp = as<StringVector>(supplement_true["stage3"]);
      stage2_supp = as<StringVector>(supplement_true["stage2"]);
      stage1_supp = as<StringVector>(supplement_true["stage1"]);
      eststage3_supp = as<StringVector>(supplement_true["eststage3"]);
      eststage2_supp = as<StringVector>(supplement_true["eststage2"]);
      eststage1_supp = as<StringVector>(supplement_true["eststage1"]);
      givenrate_supp = as<NumericVector>(supplement_true["givenrate"]);
      offset_supp = as<NumericVector>(supplement_true["offset"]);
      multiplier_supp = as<NumericVector>(supplement_true["multiplier"]);
      convtype_supp = as<IntegerVector>(supplement_true["convtype"]);
      convtype_t12_supp = as<IntegerVector>(supplement_true["convtype_t12"]);
      supp_rows = stage3_supp.length();
      
      if (supplement_true.containsElementNamed("age2")) {
        age2_supp = as<IntegerVector>(supplement_true["age2"]);
        estage2_supp = as<IntegerVector>(supplement_true["estage2"]);
      } else {
        IntegerVector age2_supp_ (stage2_supp.length(), NA_INTEGER);
        IntegerVector estage2_supp_ (eststage2_supp.length(), NA_INTEGER);
        
        age2_supp = age2_supp_;
        estage2_supp = estage2_supp_;
      }
      
      if (supplement_true.containsElementNamed("pop")) {
        pop_supp = as<StringVector>(supplement_true["pop"]);
      } else {
        StringVector pop_supp_ (stage2_supp.length(), NA_STRING);
        pop_supp = pop_supp_;
      }
      if (supplement_true.containsElementNamed("patch")) {
        patch_supp = as<StringVector>(supplement_true["patch"]);
      } else {
        StringVector patch_supp_ (stage2_supp.length(), NA_STRING);
        patch_supp = patch_supp_;
      }
      if (supplement_true.containsElementNamed("year2")) {
        year2_supp = as<StringVector>(supplement_true["year2"]);
      } else {
        StringVector year2_supp_ (stage2_supp.length(), NA_STRING);
        year2_supp = year2_supp_;
      }
      
    } else if (overwrite.isNotNull()) {
      Rcpp::DataFrame supplement_true(supplement);
      
      stage3_supp = as<StringVector>(supplement_true["stage3"]);
      stage2_supp = as<StringVector>(supplement_true["stage2"]);
      stage1_supp = as<StringVector>(supplement_true["stage1"]);
      eststage3_supp = as<StringVector>(supplement_true["eststage3"]);
      eststage2_supp = as<StringVector>(supplement_true["eststage2"]);
      eststage1_supp = as<StringVector>(supplement_true["eststage1"]);
      givenrate_supp = as<NumericVector>(supplement_true["givenrate"]);
      convtype_supp = as<IntegerVector>(supplement_true["convtype"]);
      convtype_t12_supp = as<IntegerVector>(supplement_true["convtype_t12"]);
      
      supp_rows = givenrate_supp.length();
      offset_supp = Rcpp::NumericVector::create(0., supp_rows);
      multiplier_supp = Rcpp::NumericVector::create(1., supp_rows);
      
      IntegerVector age2_supp_ (stage2_supp.length(), NA_INTEGER);
      IntegerVector estage2_supp_ (eststage2_supp.length(), NA_INTEGER);
      
      age2_supp = age2_supp_;
      estage2_supp = estage2_supp_;
      
      StringVector pop_supp_ (stage2_supp.length(), NA_STRING);
      StringVector patch_supp_ (stage2_supp.length(), NA_STRING);
      StringVector year2_supp_ (stage2_supp.length(), NA_STRING);
      
      pop_supp = pop_supp_;
      patch_supp = patch_supp_;
      year2_supp = year2_supp_;
      
    } else {
      throw Rcpp::exception("No supplement provided.", false);
    }
    
    StringVector unique_stages = unique(stagevec);
    StringVector extra_terms = {"rep", "nrep", "immat", "mat", "prop", "npr", "all", "obs", "nobs"};
    
    int no_newstages {static_cast<int>(unique_stages.length())};
    int no_extraterms {static_cast<int>(extra_terms.length())};
    
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
    
    // Check entries in supplement / overwrite table
    for (int i = 0; i < static_cast<int>(stage3_supp.length()); i++) {
      int s3supp_count {0};
      int s2supp_count {0};
      int s1supp_count {0};
      
      bool ests3_used {false};
      bool ests2_used {false};
      bool ests1_used {false};
      
      for (int j = 0; j < static_cast<int>(all_possible_stage_terms.length()); j++) {
        if (stage3_supp(i) == all_possible_stage_terms(j)) s3supp_count++;
        if (stage2_supp(i) == all_possible_stage_terms(j)) s2supp_count++;
        
        if (!StringVector::is_na(eststage3_supp(i))) {
          ests3_used = true;
        }
        if (!StringVector::is_na(eststage2_supp(i))) {
          ests2_used = true;
        }
        
        if (historical) {
          if (stage1_supp(i) == all_possible_stage_terms(j)) s1supp_count++;
          
          if (!StringVector::is_na(eststage1_supp(i))) {
            ests1_used = true;
          }
        } 
      }
      
      if (s3supp_count == 0) {
        String eat_my_shorts = "Stage names in supplement or overwrite table ";
        String eat_my_shorts1 = "(stage3) must match stageframe.";
        eat_my_shorts += eat_my_shorts1;
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
      if (s2supp_count == 0) {
        String eat_my_shorts = "Stage names in supplement or overwrite table ";
        String eat_my_shorts1 = "(stage2) must match stageframe.";
        eat_my_shorts += eat_my_shorts1;
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
      if (ests3_used) {
        if (s3supp_count == 0) {
          String eat_my_shorts = "Stage names in supplement or overwrite table ";
          String eat_my_shorts1 = "(eststage3) must match stageframe.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      if (ests2_used) {
        if (s2supp_count == 0) {
          String eat_my_shorts = "Stage names in supplement or overwrite table ";
          String eat_my_shorts1 = "(eststage2) must match stageframe.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      if (historical) {
        if (s1supp_count == 0) {
          String eat_my_shorts = "Stage names in supplement or overwrite table ";
          String eat_my_shorts1 = "(stage1) must match stageframe.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        if (ests1_used) {
          if (s1supp_count == 0) {
            String eat_my_shorts = "Stage names in supplement or overwrite table ";
            String eat_my_shorts1 = "(eststage1) must match stageframe.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
      }
    }
    
    IntegerVector s1_calls (supp_rows, 1);
    IntegerVector s2_calls (supp_rows, 1);
    IntegerVector s3_calls (supp_rows, 1);
    IntegerVector ests1_calls (supp_rows, 1);
    IntegerVector ests2_calls (supp_rows, 1);
    IntegerVector ests3_calls (supp_rows, 1);
    IntegerVector s3_planned (supp_rows, 1);
    IntegerVector s2_planned (supp_rows, 1);
    IntegerVector s1_planned (supp_rows, 1);
    
    IntegerVector s123_calls (supp_rows, 1);
    
    // Create indices for edited supplement/overwrite table
    arma::uvec alive;
    if (stageframe.containsElementNamed("alive")) {
      arma::uvec alive_temp = as<arma::uvec>(stageframe["alive"]);
      alive = alive_temp;
    } else {
      arma::uvec alive_temp (stageframe_length, fill::ones);
      alive = alive_temp;
    }
    arma::uvec repvec = as<arma::uvec>(stageframe["repstatus"]);
    arma::uvec obsvec = as<arma::uvec>(stageframe["obsstatus"]);
    arma::uvec propvec = as<arma::uvec>(stageframe["propstatus"]);
    arma::uvec immvec = as<arma::uvec>(stageframe["immstatus"]);
    arma::uvec matvec = as<arma::uvec>(stageframe["matstatus"]);
    arma::uvec indvec = as<arma::uvec>(stageframe["indataset"]);
    
    arma::uvec newprop_stages = find(propvec);
    arma::uvec newprop0_stages = find(propvec == 0);
    arma::uvec newimm_stages = find(immvec);
    arma::uvec newalive_stages = find(alive);
    arma::uvec newmat_stages1 = find(matvec);
    arma::uvec newmat_stages = intersect(newalive_stages, newmat_stages1);
    arma::uvec newrep_stages = find(repvec);
    arma::uvec newrep0_stages = find(repvec == 0);
    arma::uvec newmat_rep0_stages = intersect(newmat_stages, newrep0_stages);
    arma::uvec newobs_stages = find(obsvec);
    arma::uvec newobs0_stages = find(obsvec == 0);
    arma::uvec all_stages = find(alive);
    
    int np_s = static_cast<int>(newprop_stages.n_elem);
    int np0_s = static_cast<int>(newprop0_stages.n_elem);
    int ni_s = static_cast<int>(newimm_stages.n_elem);
    int nm_s = static_cast<int>(newmat_stages.n_elem);
    int nr_s = static_cast<int>(newrep_stages.n_elem);
    int nmr0_s = static_cast<int>(newmat_rep0_stages.n_elem);
    int no_s = static_cast<int>(newobs_stages.n_elem);
    int no0_s = static_cast<int>(newobs0_stages.n_elem);
    int a_s = static_cast<int>(all_stages.n_elem);
    
    // Build expanded supplement table
    for (int i = 0; i < supp_rows; i++) {
      s3_calls(i) = supp_decision1(as<std::string>(stage3_supp(i)), np_s, np0_s,
        ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups, groupvec,
        group_text);
      
      ests3_calls(i) = supp_decision1(as<std::string>(eststage3_supp(i)), np_s,
        np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups, groupvec,
        group_text);
      
      s2_calls(i) = supp_decision1(as<std::string>(stage2_supp(i)), np_s, np0_s,
        ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups, groupvec,
        group_text);
      
      ests2_calls(i) = supp_decision1(as<std::string>(eststage2_supp(i)), np_s,
        np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups, groupvec,
        group_text);
      
      s1_calls(i) = supp_decision1(as<std::string>(stage1_supp(i)), np_s, np0_s,
        ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups, groupvec,
        group_text);
      
      ests1_calls(i) = supp_decision1(as<std::string>(eststage1_supp(i)), np_s,
        np0_s, ni_s, nm_s, nr_s, nmr0_s, no_s, no0_s, a_s, no_groups, groupvec,
        group_text);
      
      String eat_my_shorts_gse = "If stage group shorthand is used to designate ";
      String eat_my_shorts1_gse = "both a transition and a proxy, then the ";
      String eat_my_shorts2_gse = "shorthand group must be the same in both cases.";
      eat_my_shorts_gse += eat_my_shorts1_gse;
      eat_my_shorts_gse += eat_my_shorts2_gse;
      
      if (!StringVector::is_na(eststage3_supp(i))) {
        if (eststage3_supp(i) != stage3_supp(i)) {
          if (s3_calls(i) == 1 && ests3_calls(i) > 1) {
            s3_planned(i) = ests3_calls(i);
          } else if (s3_calls(i) > 1 && ests3_calls(i) > 1) {
            throw Rcpp::exception(eat_my_shorts_gse.get_cstring(), false);
          }
        } else {
          s3_planned(i) = s3_calls(i);
        }
      } else {
        s3_planned(i) = s3_calls(i);
      }
      
      if (!StringVector::is_na(eststage2_supp(i))) {
        if (eststage2_supp(i) != stage2_supp(i)) {
          if (s2_calls(i) == 1 && ests2_calls(i) > 1) {
            s2_planned(i) = ests2_calls(i);
          } else if (s2_calls(i) > 1 && ests2_calls(i) > 1) {
            throw Rcpp::exception(eat_my_shorts_gse.get_cstring(), false);
          }
        } else {
          s2_planned(i) = s2_calls(i);
        }
      } else {
        s2_planned(i) = s2_calls(i);
      }
      
      if (!StringVector::is_na(eststage1_supp(i))) {
        if (historical && eststage1_supp(i) != stage1_supp(i)) {
          if (s1_calls(i) == 1 && ests1_calls(i) > 1) {
            s1_planned(i) = ests1_calls(i);
          } else if (s1_calls(i) > 1 && ests1_calls(i) > 1) {
            throw Rcpp::exception(eat_my_shorts_gse.get_cstring(), false);
          }
        } else if (historical) {
          s1_planned(i) = s1_calls(i);
        } else if (!historical) {
          s1_planned(i) = 1;
        }
      } else {
        s1_planned(i) = s1_calls(i);
      }
      
      s123_calls(i) = s3_planned(i) * s2_planned(i) * s1_planned(i);
    }
    
    NumericVector basepoints(supp_rows, 0.0);
    for (int i = 0; i < (supp_rows - 1); i++) {
      basepoints(i+1) = basepoints(i) + s123_calls(i);
    }
    
    int newsupp_rows = sum(s123_calls);
    
    StringVector stage3_newsupp (newsupp_rows);
    StringVector stage2_newsupp (newsupp_rows);
    StringVector stage1_newsupp (newsupp_rows);
    IntegerVector age2_newsupp (newsupp_rows);
    StringVector eststage3_newsupp (newsupp_rows);
    StringVector eststage2_newsupp (newsupp_rows);
    StringVector eststage1_newsupp (newsupp_rows);
    IntegerVector estage2_newsupp (newsupp_rows);
    NumericVector givenrate_newsupp (newsupp_rows);
    NumericVector offset_newsupp (newsupp_rows);
    NumericVector multiplier_newsupp (newsupp_rows);
    IntegerVector convtype_newsupp (newsupp_rows);
    IntegerVector convtype_t12_newsupp (newsupp_rows);
    
    StringVector pop_newsupp (newsupp_rows);
    StringVector patch_newsupp (newsupp_rows);
    StringVector year2_newsupp (newsupp_rows);
    
    int overall_counter {0};
    int group_check {0};
    
    int group_baseline3 {0};
    int group_baseline2 {0};
    int group_baseline1 {0};
    int group_baselinee3 {0};
    int group_baselinee2 {0};
    int group_baselinee1 {0};
    
    int group_ratchet3 {0};
    int group_ratchet2 {0};
    int group_ratchet1 {0};
    int group_ratchete3 {0};
    int group_ratchete2 {0};
    int group_ratchete1 {0};
    
    int prevl3 {0};
    int prevl2 {0};
    int prevl1 {0};
    int prevle3 {0};
    int prevle2 {0};
    int prevle1 {0};
    
    for (int i = 0; i < supp_rows; i++) {
      overall_counter = 0;
      for (int j = 0; j < s1_planned(i); j++) {
        for (int k = 0; k < s2_planned(i); k++) {
          for (int l = 0; l < s3_planned(i); l++) {
            stage3_newsupp(basepoints(i) + overall_counter) =
              supp_decision2(as<std::string>(stage3_supp(i)), newprop_stages,
                newprop0_stages, newimm_stages, newmat_stages, newrep_stages,
                newmat_rep0_stages, newobs_stages, newobs0_stages, all_stages,
                no_groups, groupvec, group_text, stagevec, l, group_check,
                group_ratchet3, group_baseline3, prevl3);
            
            givenrate_newsupp(basepoints(i) + overall_counter) = givenrate_supp(i);
            offset_newsupp(basepoints(i) + overall_counter) = offset_supp(i);
            multiplier_newsupp(basepoints(i) + overall_counter) = multiplier_supp(i);
            convtype_newsupp(basepoints(i) + overall_counter) = convtype_supp(i);
            convtype_t12_newsupp(basepoints(i) + overall_counter) = convtype_t12_supp(i);
            
            eststage3_newsupp(basepoints(i) + overall_counter) =
              supp_decision2(as<std::string>(eststage3_supp(i)), newprop_stages,
                newprop0_stages, newimm_stages, newmat_stages, newrep_stages,
                newmat_rep0_stages, newobs_stages, newobs0_stages, all_stages,
                no_groups, groupvec, group_text, stagevec, l, group_check,
                group_ratchete3, group_baselinee3, prevle3);
            
            stage2_newsupp(basepoints(i) + overall_counter) =
              supp_decision2(as<std::string>(stage2_supp(i)), newprop_stages,
                newprop0_stages, newimm_stages, newmat_stages, newrep_stages,
                newmat_rep0_stages, newobs_stages, newobs0_stages, all_stages,
                no_groups, groupvec, group_text, stagevec, k, group_check,
                group_ratchet2, group_baseline2, prevl2);
            
            eststage2_newsupp(basepoints(i) + overall_counter) =
              supp_decision2(as<std::string>(eststage2_supp(i)), newprop_stages,
                newprop0_stages, newimm_stages, newmat_stages, newrep_stages,
                newmat_rep0_stages, newobs_stages, newobs0_stages, all_stages,
                no_groups, groupvec, group_text, stagevec, k, group_check,
                group_ratchete2, group_baselinee2, prevle2);
            
            stage1_newsupp(basepoints(i) + overall_counter) =
              supp_decision2(as<std::string>(stage1_supp(i)), newprop_stages,
                newprop0_stages, newimm_stages, newmat_stages, newrep_stages,
                newmat_rep0_stages, newobs_stages, newobs0_stages, all_stages,
                no_groups, groupvec, group_text, stagevec, j, group_check,
                group_ratchet1, group_baseline1, prevl1);
            
            eststage1_newsupp(basepoints(i) + overall_counter) =
              supp_decision2(as<std::string>(eststage1_supp(i)), newprop_stages,
                newprop0_stages, newimm_stages, newmat_stages, newrep_stages,
                newmat_rep0_stages, newobs_stages, newobs0_stages, all_stages,
                no_groups, groupvec, group_text, stagevec, j, group_check,
                group_ratchete1, group_baselinee1, prevle1);
            
            age2_newsupp(basepoints(i) + overall_counter) = age2_supp(i);
            estage2_newsupp(basepoints(i) + overall_counter) = estage2_supp(i);
            
            pop_newsupp(basepoints(i) + overall_counter) = pop_supp(i);
            patch_newsupp(basepoints(i) + overall_counter) = patch_supp(i);
            year2_newsupp(basepoints(i) + overall_counter) = year2_supp(i);
            
            overall_counter++;
          }
        }
      }
    }
    
    Rcpp::List newsupplement(16);
    
    newsupplement(0) = stage3_newsupp;
    newsupplement(1) = stage2_newsupp;
    newsupplement(2) = stage1_newsupp;
    newsupplement(3) = age2_newsupp;
    newsupplement(4) = eststage3_newsupp;
    newsupplement(5) = eststage2_newsupp;
    newsupplement(6) = eststage1_newsupp;
    newsupplement(7) = estage2_newsupp;
    newsupplement(8) = givenrate_newsupp;
    newsupplement(9) = offset_newsupp;
    newsupplement(10) = multiplier_newsupp;
    newsupplement(11) = convtype_newsupp;
    newsupplement(12) = convtype_t12_newsupp;
    newsupplement(13) = pop_newsupp;
    newsupplement(14) = patch_newsupp;
    newsupplement(15) = year2_newsupp;
    
    CharacterVector su_namevec = {"stage3", "stage2", "stage1", "age2",
      "eststage3", "eststage2", "eststage1", "estage2", "givenrate", "offset",
      "multiplier", "convtype", "convtype_t12", "pop", "patch", "year2"};
    CharacterVector su_newclasses = {"data.frame", "lefkoSD"};
    newsupplement.attr("names") = su_namevec;
    newsupplement.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newsupp_rows);
    newsupplement.attr("class") = su_newclasses;
    
    return newsupplement;
  }
  
  //' Expand Supplemental Table by Age Inputs
  //' 
  //' The function takes an already expanded supplemental table as input and
  //' produces an edited version expanded by age. This is used in editing
  //' age-by-stage MPMs.
  //' 
  //' @name age_expanded
  //' 
  //' @param supplement The expanded supplemental table.
  //' @param minage An integer denoting the minimum age used in the MPM.
  //' @param maxage An integer denoting the maximum age used in the MPM.
  //' 
  //' @return This function returns a new data frame that acts as the expanded
  //' supplemental table without shorthand codes and with ages explicit.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::DataFrame age_expanded (DataFrame supplement, int minage, int maxage) {
    
    StringVector supp_stage3 = as<StringVector>(supplement["stage3"]);
    StringVector supp_stage2 = as<StringVector>(supplement["stage2"]);
    StringVector supp_stage1 = as<StringVector>(supplement["stage1"]);
    IntegerVector supp_age2 = as<IntegerVector>(supplement["age2"]);
    StringVector supp_eststage3 = as<StringVector>(supplement["eststage3"]);
    StringVector supp_eststage2 = as<StringVector>(supplement["eststage2"]);
    StringVector supp_eststage1 = as<StringVector>(supplement["eststage1"]);
    IntegerVector supp_estage2 = as<IntegerVector>(supplement["estage2"]);
    NumericVector supp_givenrate = as<NumericVector>(supplement["givenrate"]);
    NumericVector supp_offset = as<NumericVector>(supplement["offset"]);
    NumericVector supp_multiplier = as<NumericVector>(supplement["multiplier"]);
    IntegerVector supp_convtype = as<IntegerVector>(supplement["convtype"]);
    IntegerVector supp_convtype_t12 = as<IntegerVector>(supplement["convtype_t12"]);
    
    StringVector supp_pop;
    StringVector supp_patch;
    StringVector supp_year2;
    
    if (supplement.containsElementNamed("pop")) {
      supp_pop = as<StringVector>(supplement["pop"]);
      supp_patch = as<StringVector>(supplement["patch"]);
      supp_year2 = as<StringVector>(supplement["year2"]);
    } else {
      int new_string_length = static_cast<int>(supp_age2.length());
      StringVector new1s(new_string_length);
      for (int i = 0; i < new_string_length; i++) {
        new1s(i) = NA_STRING;
      }
      
      supp_pop = new1s;
      supp_patch = clone(new1s);
      supp_year2 = clone(new1s);
    }
    
    int base_length = static_cast<int>(supp_stage2.length());
    arma::uvec estimated_elements (base_length, fill::zeros);
    
    int total_ages = maxage - minage + 1; // Originally maxage - minage + 1
    IntegerVector total_ages_seq = seq(minage, maxage); // Originally seq(minage, maxage)
    
    for (int i = 0; i < base_length; i++) {
      if (IntegerVector::is_na(supp_age2(i))) {
        estimated_elements(i) = total_ages;
      } else {
        estimated_elements(i) = 1;
      }
    }
    
    int new_length = static_cast<int>(accu(estimated_elements));
    
    StringVector new_supp_stage3 (new_length);
    StringVector new_supp_stage2 (new_length);
    StringVector new_supp_stage1 (new_length);
    IntegerVector new_supp_age2 (new_length);
    StringVector new_supp_eststage3 (new_length);
    StringVector new_supp_eststage2 (new_length);
    StringVector new_supp_eststage1 (new_length);
    IntegerVector new_supp_estage2 (new_length);
    NumericVector new_supp_givenrate (new_length);
    NumericVector new_supp_offset (new_length);
    NumericVector new_supp_multiplier (new_length);
    IntegerVector new_supp_convtype (new_length);
    IntegerVector new_supp_convtype_t12 (new_length);
    StringVector new_supp_pop (new_length);
    StringVector new_supp_patch (new_length);
    StringVector new_supp_year2 (new_length);
    
    int cookie_counter {0};
    for (int i = 0; i < base_length; i++) {
      for (int j = 0; j < static_cast<int>(estimated_elements(i)); j++) {
        new_supp_stage3(cookie_counter) = supp_stage3(i);
        new_supp_stage2(cookie_counter) = supp_stage2(i);
        new_supp_stage1(cookie_counter) = supp_stage1(i);
        
        if (estimated_elements(i) > 1) {
          new_supp_age2(cookie_counter) = total_ages_seq(j);
        } else {
          new_supp_age2(cookie_counter) = supp_age2(i);
        }
        
        new_supp_eststage3(cookie_counter) = supp_eststage3(i);
        new_supp_eststage2(cookie_counter) = supp_eststage2(i);
        new_supp_eststage1(cookie_counter) = supp_eststage1(i);
        
        if (estimated_elements(i) > 1) {
          new_supp_estage2(cookie_counter) = total_ages_seq(j);
        } else {
          new_supp_estage2(cookie_counter) = supp_estage2(i);
        }
        
        new_supp_givenrate(cookie_counter) = supp_givenrate(i);
        new_supp_offset(cookie_counter) = supp_offset(i);
        new_supp_multiplier(cookie_counter) = supp_multiplier(i);
        new_supp_convtype(cookie_counter) = supp_convtype(i);
        new_supp_convtype_t12(cookie_counter) = supp_convtype_t12(i);
        
        new_supp_pop(cookie_counter) = supp_pop(i);
        new_supp_patch(cookie_counter) = supp_patch(i);
        new_supp_year2(cookie_counter) = supp_year2(i);
        
        cookie_counter++;
      }
    }
    
    Rcpp::DataFrame new_supplement = DataFrame::create(_["stage3"] = new_supp_stage3,
      _["stage2"] = new_supp_stage2, _["stage1"] = new_supp_stage1,
      _["age2"] = new_supp_age2, _["eststage3"] = new_supp_eststage3,
      _["eststage2"] = new_supp_eststage2, _["eststage1"] = new_supp_eststage1,
      _["estage2"] = new_supp_estage2, _["givenrate"] = new_supp_givenrate,
      _["offset"] = new_supp_offset, _["multiplier"] = new_supp_multiplier,
      _["convtype"] = new_supp_convtype,
      _["convtype_t12"] = new_supp_convtype_t12, _["pop"] = new_supp_pop,
      _["patch"] = new_supp_patch, _["year2"] = new_supp_year2);
    
    return new_supplement;
  }
  
  //' Create hstages Index Object
  //' 
  //' Function \code{hst_maker()} creates \code{hstages} index data frames from
  //' \code{stageframe} inputs.
  //' 
  //' @name hst_maker
  //' 
  //' @param sframe The ahistorical stageframe used in MPM development.
  //' @param format An integer value creating an Ehrlen format data frame if
  //' \code{format = 1}, and a deVries format data frame if \code{format = 2}.
  //' Defaults to \code{format = 1}.
  //' 
  //' @return Creates a data frame by reference with the following columns:
  //' \item{stage_id_2}{Integer index of stage in time \emph{t}+1.}
  //' \item{stage_id_1}{Integer index of stage in time \emph{t}.}
  //' \item{stage_2}{String name of stage in time \emph{t}+1.}
  //' \item{stage_1}{String name of stage in time \emph{t}.}
  //' 
  //' @keywords internal
  //' @noRd
  inline void hst_maker (DataFrame& hstages_out, const DataFrame& sframe,
    int format = 1) {
    
    StringVector stage_name = as<StringVector>(sframe["stage"]);
    int true_stages = stage_name.length();
    
    IntegerVector stage_id = seq(1, true_stages);
    int correction {0};
    if (format == 2) correction = 1;
    
    int h_stages = (true_stages - correction) * true_stages;
    
    IntegerVector stage_id_2 (h_stages);
    IntegerVector stage_id_1 (h_stages);
    StringVector stage_2 (h_stages);
    StringVector stage_1 (h_stages);
    
    int current_elem {0};
    for (int s1 = 0; s1 < true_stages; s1++) {
      for (int s2 = 0; s2 < (true_stages - correction); s2++) {
        stage_id_2[current_elem] = stage_id[s2];
        stage_id_1[current_elem] = stage_id[s1];
        stage_2[current_elem] = stage_name[s2];
        stage_1[current_elem] = stage_name[s1];
        current_elem++;
      }
    }
    
    DataFrame output = DataFrame::create(_["stage_id_2"] = stage_id_2,
      _["stage_id_1"] = stage_id_1, _["stage_2"] = stage_2,
      _["stage_1"] = stage_1);
    
    hstages_out = output;
  }
  
  //' Create agestages Index Object
  //' 
  //' Function \code{age_maker()} creates \code{agestages} index data frames from
  //' \code{stageframe} inputs.
  //' 
  //' @name age_maker
  //' 
  //' @param sframe The ahistorical stageframe used in MPM development.
  //' 
  //' @return A data frame with the following columns:
  //' \item{stage_id}{Integer index of stage.}
  //' \item{stage}{String name of stage.}
  //' \item{age}{The age of stage in current time.}
  //' 
  //' @keywords internal
  //' @noRd
  inline DataFrame age_maker (const DataFrame& sframe, int start_age, int last_age) {
    StringVector stage_name = as<StringVector>(sframe["stage"]);
    int true_stages = stage_name.length();
    
    IntegerVector stage_id = seq(1, true_stages);
    IntegerVector all_ages = seq(start_age, last_age);
    int num_ages = all_ages.length();
    int age_stages = true_stages * num_ages;
    
    IntegerVector stage_id_new (age_stages);
    StringVector stage_new (age_stages);
    IntegerVector age_new (age_stages);
    
    for (int s1 = 0; s1 < num_ages; s1++) {
      for (int s2 = 0; s2 < true_stages; s2++) {
        int current_elem = (s1 * true_stages) + s2;
        
        stage_id_new[current_elem] = stage_id[s2];
        stage_new[current_elem] = stage_name[s2];
        age_new[current_elem] = all_ages[s1];
      }
    }
    
    DataFrame output = DataFrame::create(_["stage_id"] = stage_id_new,
      _["stage"] = stage_new, _["age"] = age_new);
    
    return output;
  }
  
  //' Create Element Index for Matrix Estimation
  //' 
  //' Function \code{theoldpizzle()} creates a data frame object used by 
  //' functions \code{specialpatrolgroup()}, \code{normalpatrolgroup()},
  //' \code{subvertedpatrolgroup()}, and \code{jerzeibalowski()} to estimate
  //' raw and function-based matrices.
  //' 
  //' @name theoldpizzle
  //'
  //' @param StageFrame The stageframe object identifying the life history model
  //' being operationalized.
  //' @param OverWrite The supplement or overwrite table used in analysis, as
  //' modified by \code{sf_reassess()}.
  //' @param repmatrix The reproductive matrix used in analysis.
  //' @param firstage The first age to be used in the analysis. Should typically
  //' be \code{0} for pre-breeding and \code{1} for post-breeding life history
  //' models. If not building age-by-stage MPMs, then should be set to \code{0}.
  //' @param finalage The final age to be used in analysis. If not building
  //' age-by-stage MPMs, then should be set to \code{0}.
  //' @param format Indicates whether historical matrices should be in (\code{1})
  //' Ehrlen or (\code{2}) deVries format.
  //' @param style The style of analysis, where \code{0} is historical, \code{1}
  //' is ahistorical, and \code{2} is age-by-stage.
  //' @param cont Denotes whether age-by-stage matrix continues past the final
  //' age.
  //' @param filter An integer denoting whether to filter the output data frame to
  //' eliminate unusable rows, and if so, how to do so. Possible values: \code{0}:
  //' no filtering, \code{1}: filter out rows with \code{index321 == -1}, and
  //' \code{2}: filter out rows with \code{aliveandequal == -1}.
  //' 
  //' @return The output is a large data frame describing every element to be
  //' estimated in matrices.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List theoldpizzle(const DataFrame& StageFrame, const DataFrame& OverWrite,
    const arma::mat& repmatrix, int firstage, int finalage, int format, int style,
    int cont, int filter) {
    
    StringVector ovstage3 = as<StringVector>(OverWrite["stage3"]);
    StringVector ovstage2 = as<StringVector>(OverWrite["stage2"]);
    StringVector ovstage1 = as<StringVector>(OverWrite["stage1"]);
    StringVector oveststage3 = as<StringVector>(OverWrite["eststage3"]);
    StringVector oveststage2 = as<StringVector>(OverWrite["eststage2"]);
    StringVector oveststage1 = as<StringVector>(OverWrite["eststage1"]);
    arma::vec ovgivenrate = as<arma::vec>(OverWrite["givenrate"]);
    arma::vec ovoffset = as<arma::vec>(OverWrite["offset"]);
    arma::vec ovmultiplier = as<arma::vec>(OverWrite["multiplier"]);
    arma::vec ovconvtype = as<arma::vec>(OverWrite["convtype"]);
    arma::vec ovconvt12 = as<arma::vec>(OverWrite["convtype_t12"]);
    int ovrows = static_cast<int>(ovconvtype.n_elem);
    
    IntegerVector ovage2;
    IntegerVector ovestage2;
    if (OverWrite.containsElementNamed("age2")) {
      ovage2 = as<IntegerVector>(OverWrite["age2"]);
      ovestage2 = as<IntegerVector>(OverWrite["estage2"]);
    }
    
    int totalages = (finalage - firstage) + 1;
    
    arma::vec ovindex3 (ovrows); // Originally (ovrows * totalages)
    arma::vec ovindex2 (ovrows);
    arma::vec ovindex1 (ovrows);
    arma::vec ovnew3 (ovrows);
    arma::vec ovnew2 (ovrows);
    arma::vec ovnew1 (ovrows);
    arma::vec ovindexold321 (ovrows);
    arma::vec ovindexnew321 (ovrows);
    arma::vec ovnewgivenrate (ovrows);
    arma::vec ovnewoffset (ovrows, fill::zeros);
    arma::vec ovnewmultiplier (ovrows, fill::zeros);
    arma::vec ovconvtypeage (ovrows);
    ovindex3.fill(-1.0);
    ovindex2.fill(-1.0);
    ovindex1.fill(-1.0);
    ovnew3.fill(-1.0);
    ovnew2.fill(-1.0);
    ovnew1.fill(-1.0);
    ovindexold321.fill(-1.0);
    ovindexnew321.fill(-1.0);
    ovnewgivenrate.fill(-1.0);
    ovconvtypeage.fill(-1.0);
    
    arma::ivec newstageid = as<arma::ivec>(StageFrame["stage_id"]);
    StringVector origstageid = as<StringVector>(StageFrame["stage"]);
    arma::vec binsizectr = as<arma::vec>(StageFrame["sizebin_center"]);
    arma::vec repstatus = as<arma::vec>(StageFrame["repstatus"]);
    arma::vec obsstatus = as<arma::vec>(StageFrame["obsstatus"]);
    arma::vec immstatus = as<arma::vec>(StageFrame["immstatus"]);
    arma::vec matstatus = as<arma::vec>(StageFrame["matstatus"]);
    arma::vec indata = as<arma::vec>(StageFrame["indataset"]);
    arma::vec binsizewidth = as<arma::vec>(StageFrame["sizebin_width"]);
    arma::vec alive = as<arma::vec>(StageFrame["alive"]);
    arma::vec minage = as<arma::vec>(StageFrame["min_age"]);
    arma::vec maxage = as<arma::vec>(StageFrame["max_age"]);
    arma::vec group = as<arma::vec>(StageFrame["group"]);
    arma::vec almostborn = as<arma::vec>(StageFrame["almostborn"]);
    
    arma::vec binsizebctr = as<arma::vec>(StageFrame["sizebinb_center"]);
    arma::vec binsizecctr = as<arma::vec>(StageFrame["sizebinc_center"]);
    arma::vec binsizebwidth = as<arma::vec>(StageFrame["sizebinb_width"]);
    arma::vec binsizecwidth = as<arma::vec>(StageFrame["sizebinc_width"]);
    
    // Determine length of matrix map data frame
    int nostages = static_cast<int>(newstageid.n_elem);
    int nostages_nodead = nostages - 1;
    int nostages_nounborn = nostages;
    int nostages_nodead_nounborn = nostages_nodead;
    int prior_stage = -1;
    arma::vec ovrepentry_prior(nostages, fill::zeros);
    IntegerVector stageorder = seq(1, nostages);
    int totallength {0};
    
    if (style == 2) {
      totallength = (nostages * nostages * totalages * totalages);
    } else if (style == 1) {
      totallength = (nostages * nostages_nodead);
    } else {
      if (format == 2) {
        nostages_nodead_nounborn = nostages - 2;
        prior_stage = nostages_nodead_nounborn;
        nostages_nounborn = nostages - 1;
        totallength = (2 * nostages_nodead_nounborn * nostages_nounborn * nostages_nounborn);
      } else {
        totallength = (nostages * (nostages_nodead * nostages_nodead));
      }
    }
    
    // Set up repmatrix
    int reprows = repmatrix.n_rows;
    int repmattype = 0;
    
    if (reprows == (nostages - 1) || reprows == (nostages - 2)) {
      repmattype = 1; // repmatrix is ahistorical
    } else if (reprows == ((nostages - 1) * (nostages - 1)) || 
        reprows == ((nostages - 2) * (nostages - 2))) {
      repmattype = 2; // repmatrix is historical
    }
    
    // Set up vectors that will be put together into matrix map data frame
    arma::ivec stage3 (totallength, fill::zeros);
    arma::ivec stage2n (totallength, fill::zeros);
    arma::ivec stage2o (totallength, fill::zeros);
    arma::ivec stage1 (totallength, fill::zeros);
    
    arma::ivec stageorder3 (totallength, fill::zeros);
    arma::ivec stageorder2n (totallength, fill::zeros);
    arma::ivec stageorder2o (totallength, fill::zeros);
    arma::ivec stageorder1 (totallength, fill::zeros);
    
    arma::vec size3 (totallength, fill::zeros);
    arma::vec size2n (totallength, fill::zeros);
    arma::vec size2o (totallength, fill::zeros);
    arma::vec size1 (totallength, fill::zeros);
    
    arma::vec sizeb3 (totallength, fill::zeros);
    arma::vec sizeb2n (totallength, fill::zeros);
    arma::vec sizeb2o (totallength, fill::zeros);
    arma::vec sizeb1 (totallength, fill::zeros);
    
    arma::vec sizec3 (totallength, fill::zeros);
    arma::vec sizec2n (totallength, fill::zeros);
    arma::vec sizec2o (totallength, fill::zeros);
    arma::vec sizec1 (totallength, fill::zeros);
    
    arma::vec obs3 (totallength, fill::zeros);
    arma::vec obs2n (totallength, fill::zeros);
    arma::vec obs2o (totallength, fill::zeros);
    arma::vec obs1 (totallength, fill::zeros);
    
    arma::vec rep3 (totallength, fill::zeros);
    arma::vec rep2n (totallength, fill::zeros);
    arma::vec rep2o (totallength, fill::zeros);
    arma::vec rep1 (totallength, fill::zeros);
    
    arma::vec mat3 (totallength, fill::zeros);
    arma::vec mat2n (totallength, fill::zeros);
    arma::vec mat2o (totallength, fill::zeros);
    arma::vec mat1 (totallength, fill::zeros);
    
    arma::vec imm3 (totallength, fill::zeros);
    arma::vec imm2n (totallength, fill::zeros);
    arma::vec imm2o (totallength, fill::zeros);
    arma::vec imm1 (totallength, fill::zeros);
    
    arma::vec repentry3 (totallength, fill::zeros);
    arma::vec repentry2o (totallength, fill::zeros);
    arma::vec almostborn1 (totallength, fill::zeros);
    
    arma::vec binwidth (totallength, fill::zeros);
    arma::vec binbwidth (totallength, fill::zeros);
    arma::vec bincwidth (totallength, fill::zeros);
    
    arma::vec indata3 (totallength, fill::zeros);
    arma::vec indata2n (totallength, fill::zeros);
    arma::vec indata2o (totallength, fill::zeros);
    arma::vec indata1 (totallength, fill::zeros);
    
    arma::vec minage3 (totallength, fill::zeros);
    arma::vec minage2 (totallength, fill::zeros);
    arma::vec maxage3 (totallength, fill::zeros);
    arma::vec maxage2 (totallength, fill::zeros);
    
    arma::vec grp3 (totallength, fill::zeros);
    arma::vec grp2n (totallength, fill::zeros);
    arma::vec grp2o (totallength, fill::zeros);
    arma::vec grp1 (totallength, fill::zeros);
    
    arma::vec actualage (totallength, fill::zeros);
    arma::vec index321 (totallength); // No death transitions
    arma::vec index321d (totallength); // Death transitions included
    arma::vec index21 (totallength);
    arma::vec indatalong (totallength, fill::zeros);
    arma::vec aliveequal (totallength);
    arma::vec included (totallength, fill::zeros);
    index321.fill(-1.0);
    index321d.fill(-1.0);
    index21.fill(-1.0);
    aliveequal.fill(-1.0);
    
    arma::mat asadditions (totallength, 5, fill::zeros);
    
    arma::vec ovgivent (totallength);
    arma::vec ovgivenf (totallength);
    arma::vec ovestt (totallength);
    arma::vec ovestf (totallength);
    arma::vec ovoffsett (totallength, fill::zeros);
    arma::vec ovoffsetf (totallength, fill::zeros);
    arma::vec ovrepentry (totallength, fill::zeros);
    arma::vec ovsurvmult (totallength, fill::ones);
    arma::vec ovfecmult (totallength, fill::ones);
    ovgivent.fill(-1.0);
    ovestt.fill(-1.0);
    ovgivenf.fill(-1.0);
    ovestf.fill(-1.0);
    
    int repm_elem {-1};
    double deadandnasty {0};
    long long int currentindex {0};
    
    // Change stage names to stage numbers per input stageframe for styles 0 and 1
    if (style < 2) {
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
            for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
              if (ovstage3(i) == origstageid(j)) {
                ovindex3(i) = newstageid(j);
              }
              
              if (ovstage2(i) == origstageid(j)) {
                ovindex2(i) = newstageid(j);
              }
              
              if (ovstage1(i) == origstageid(j)) {
                ovindex1(i) = newstageid(j);
              }
              
              if (oveststage3(i) == origstageid(j)) {
                ovnew3(i) = newstageid(j);
              }
              
              if (oveststage2(i) == origstageid(j)) {
                ovnew2(i) = newstageid(j);
              }
              
              if (oveststage1(i) == origstageid(j)) {
                ovnew1(i) = newstageid(j);
              }
            } // j for loop
          } // i for loop
        } // ovrows if statement
      }
    } // style if statement
    
    // Main data frame creation loops
    if (style == 0 && format == 2) { // Historical MPM deVries format
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          for (int i = 0; i < ovrows; i++) {  // Loop across overwrite rows
            if (ovconvtype(i) > 1.0) { // Catches all changes to fecundity and reproductive multipliers
              ovindexold321(i) = (ovindex3(i) - 1) + (prior_stage * nostages) + 
                ((ovindex2(i) - 1) * nostages * nostages) + 
                ((ovindex1(i) - 1) * nostages * nostages * nostages);
                
              ovindexnew321(i) = (ovnew3(i) - 1) + (prior_stage * nostages) + 
                ((ovnew2(i) - 1) * nostages * nostages) + 
                ((ovnew1(i) - 1) * nostages * nostages * nostages);
            } else if (ovconvt12(i) == 2.0) { // Catches all survival terms with historical reproduction events
              ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages) + 
                ((ovindex2(i) - 1) * nostages * nostages) + 
                (prior_stage * nostages * nostages * nostages);
                
              ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages) + 
                ((ovnew2(i) - 1) * nostages * nostages) + 
                (prior_stage * nostages * nostages * nostages);
            } else { // Full survival transitions
              ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages) + 
                ((ovindex2(i) - 1) * nostages * nostages) + 
                ((ovindex1(i) - 1) * nostages * nostages * nostages);
                
              ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages) + 
                ((ovnew2(i) - 1) * nostages * nostages) + 
                ((ovnew1(i) - 1) * nostages * nostages * nostages);
            }
            if (ovindexold321(i) < 0.0) ovindexold321(i) = -1.0;
            if (ovindexnew321(i) < 0.0) ovindexnew321(i) = -1.0;
            
            if (!NumericVector::is_na(ovgivenrate(i))) {
              ovnewgivenrate(i) = ovgivenrate(i);
            }
            if (!NumericVector::is_na(ovoffset(i))) {
              if (ovoffset(i) != 0.) ovnewoffset(i) = ovoffset(i);
            }
            if (NumericVector::is_na(ovmultiplier(i))) {
              ovmultiplier(i) = 1;
            }
            ovnewmultiplier(i) = ovmultiplier(i);
            
            if (ovconvtype(i) == 3.0) {
              for (int j = 0; j < nostages; j++) {
                if (origstageid(j) == ovstage3(i)) ovrepentry_prior(j) = 1.0;
              }
            }
          } // i for loop
        } // ovrows if statement
      }
      
      arma::uvec marked_for_repentry (nostages, fill::zeros); // Only in deVries format
      
      for (int time1 = 0; time1 < nostages_nodead; time1++) {
        for (int time2o = 0; time2o < nostages_nodead_nounborn; time2o++) {
          for (int time2n = 0; time2n < nostages; time2n++) {
            for (int time3 = 0; time3 < nostages; time3++) {
              
              if (time3 != prior_stage) {
                if (time2n == time2o || time2n == prior_stage){
                  
                  included(currentindex) = 1.0;
                  
                  stageorder3(currentindex) = stageorder(time3);
                  stageorder2n(currentindex) = stageorder(time2n);
                  stageorder2o(currentindex) = stageorder(time2o);
                  stageorder1(currentindex) = stageorder(time1);
                  
                  stage3(currentindex) = newstageid(time3);
                  stage2n(currentindex) = newstageid(time2n);
                  stage2o(currentindex) = newstageid(time2o);
                  stage1(currentindex) = newstageid(time1);
                  
                  size3(currentindex) = binsizectr(time3);
                  size2n(currentindex) = binsizectr(time2n);
                  size2o(currentindex) = binsizectr(time2o);
                  size1(currentindex) = binsizectr(time1);
                  
                  sizeb3(currentindex) = binsizebctr(time3);
                  sizeb2n(currentindex) = binsizebctr(time2n);
                  sizeb2o(currentindex) = binsizebctr(time2o);
                  sizeb1(currentindex) = binsizebctr(time1);
                  
                  if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0.0;
                  
                  sizec3(currentindex) = binsizecctr(time3);
                  sizec2n(currentindex) = binsizecctr(time2n);
                  sizec2o(currentindex) = binsizecctr(time2o);
                  sizec1(currentindex) = binsizecctr(time1);
                  
                  if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0.0;
                  
                  obs3(currentindex) = obsstatus(time3);
                  obs2n(currentindex) = obsstatus(time2n);
                  obs2o(currentindex) = obsstatus(time2o);
                  obs1(currentindex) = obsstatus(time1);
                  
                  rep3(currentindex) = repstatus(time3);
                  rep2n(currentindex) = repstatus(time2n);
                  rep2o(currentindex) = repstatus(time2o);
                  rep1(currentindex) = repstatus(time1);
                  
                  mat3(currentindex) = matstatus(time3);
                  mat2n(currentindex) = matstatus(time2n);
                  mat2o(currentindex) = matstatus(time2o);
                  mat1(currentindex) = matstatus(time1);
                  
                  imm3(currentindex) = immstatus(time3);
                  imm2n(currentindex) = immstatus(time2n);
                  imm2o(currentindex) = immstatus(time2o);
                  imm1(currentindex) = immstatus(time1);
                  
                  // Fill in repentry info from repmatrix
                  if (time2n == prior_stage && time3 < prior_stage && time2o < prior_stage) {
                    if (repmattype == 1) { 
                      repm_elem = (time3 + (time2o * nostages_nodead_nounborn));
                    } else if (repmattype == 2) {
                      repm_elem = time3 + (time2o * nostages_nodead_nounborn) + 
                        (time2o * nostages_nodead_nounborn * nostages_nodead_nounborn) +
                        (time1 * nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn);
                    } else repm_elem = -1;
                    
                    if (repmatrix(repm_elem) > 0.0) {
                      repentry3(currentindex) = repmatrix(repm_elem);
                      if (repentry3(currentindex) == 0.0 && ovrepentry_prior(time3) != 0.0) {
                        repentry3(currentindex) = 1.0;
                        marked_for_repentry(stage3(currentindex)) = 1;
                      } 
                    }
                  } else repentry3(currentindex) = 0.0;
                  
                  almostborn1(currentindex) = almostborn(time1);
                  
                  indata3(currentindex) = indata(time3);
                  indata2n(currentindex) = indata(time2n);
                  indata2o(currentindex) = indata(time2o);
                  indata1(currentindex) = indata(time1);
                  
                  binwidth(currentindex) = binsizewidth(time3);
                  binbwidth(currentindex) = binsizebwidth(time3);
                  bincwidth(currentindex) = binsizecwidth(time3);
                  
                  if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
                  if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                  
                  minage3(currentindex) = minage(time3);
                  minage2(currentindex) = minage(time2o);
                  maxage3(currentindex) = maxage(time3);
                  maxage2(currentindex) = maxage(time2o);
                  actualage(currentindex) = 0.0;
                  
                  grp3(currentindex) = group(time3);
                  grp2n(currentindex) = group(time2n);
                  grp2o(currentindex) = group(time2o);
                  grp1(currentindex) = group(time1);
                  
                  if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
                    deadandnasty = 1.0;
                  } else if (stage2o(currentindex) == nostages || stage1(currentindex) == nostages) {
                    deadandnasty = 1.0;
                  } else {
                    deadandnasty = 0.0;
                  }
                  
                  // Required for proper fecundity estimation in rlefko3
                  index321d(currentindex) = (stage3(currentindex) - 1) + 
                    ((stage2n(currentindex) - 1) * nostages) + 
                    ((stage2o(currentindex) - 1) * nostages * nostages) + 
                    ((stage1(currentindex) - 1) * nostages * nostages * nostages);
                  
                  if (deadandnasty == 0.0) {
                    // Next index variable gives element in the final matrix
                    aliveequal(currentindex) = (stageorder3(currentindex) - 1) + 
                      ((stageorder2n(currentindex) - 1) * nostages_nodead_nounborn) + 
                      ((stageorder2o(currentindex) - 1) * nostages_nodead * nostages_nodead_nounborn) + 
                      ((stageorder1(currentindex) - 1) * nostages_nodead_nounborn * 
                        nostages_nodead * nostages_nodead_nounborn);
                    
                    // Next two index variables used by ovreplace
                    index321(currentindex) = (stage3(currentindex) - 1) + 
                      ((stage2n(currentindex) - 1) * nostages) + 
                      ((stage2o(currentindex) - 1) * nostages * nostages) + 
                      ((stage1(currentindex) - 1) * nostages * nostages * nostages);
                      
                    index21(currentindex) = (stage2o(currentindex) - 1) + 
                      ((stage1(currentindex) - 1) * nostages);
                  }
                  
                  indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
                    indata2o(currentindex) * indata1(currentindex);
                  
                  currentindex += 1;
                } // if (time2n == tim2o || time2n == prior_stage) statement
              } // if (time3n != dead_stage) statement
            } // time3 loop
          } // time2n loop
        } // time2o loop
      } // time1 loop 
      
      // Edit data frame to make sure that almostborn situations in time 1
      // lead to estimated elements only if a repentry stage occurs in time 2
      arma::uvec marked_only = find(marked_for_repentry);
      if (marked_only.n_elem > 0) {
        for (int i = 0; i < static_cast<int>(marked_only.n_elem); i++) {
          arma::uvec total_indices_to_change = find(stage2o == marked_only(i));
          
          if (total_indices_to_change.n_elem > 0) {
            for (int j = 0; j < static_cast<int>(total_indices_to_change.n_elem); j++) {
              repentry2o(total_indices_to_change(j)) = 1;
            }
          }
        }
      }
      
      arma::uvec alm_only = find(almostborn1);
      if (alm_only.n_elem > 0) {
        for (int i = 0; i < static_cast<int>(alm_only.n_elem); i++) {
          if (repentry2o(alm_only(i)) < 1.0) {
            index321(alm_only(i)) = -1.0;
          }
        }
      }
      
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          asadditions = LefkoMats::ovreplace(index321, ovindexold321,
            ovindexnew321, ovconvtype, ovnew3, ovnewgivenrate, ovnewoffset,
            ovnewmultiplier);
          
          ovgivent = asadditions.col(0);
          ovestt = asadditions.col(1);
          ovgivenf = asadditions.col(2);
          ovestf = asadditions.col(3);
          
          ovrepentry = asadditions.col(4);
          ovsurvmult = asadditions.col(5);
          ovfecmult = asadditions.col(6);
          ovoffsett = asadditions.col(7);
          ovoffsetf = asadditions.col(8);
          
          arma::uvec workedupindex = find(ovrepentry > 0.0);
          int changedreps = static_cast<int>(workedupindex.n_elem);
          
          if (changedreps > 0) {
            for (int i = 0; i < changedreps; i++) {
              repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
            }
          }
        } // ovreplace if statement
      }
    } else if (style == 0 && format == 1) { // Historical MPM Ehrlen format
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
            ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages_nodead_nounborn) + 
              ((ovindex2(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
              ((ovindex1(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                nostages_nodead_nounborn);
              
            ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages_nodead) + 
              ((ovnew2(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
              ((ovnew1(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                nostages_nodead_nounborn);
            
            if (ovindexold321(i) < 0) ovindexold321(i) = -1.0;
            if (ovindexnew321(i) < 0) ovindexnew321(i) = -1.0;
            
            if (!NumericVector::is_na(ovgivenrate(i))) {
              ovnewgivenrate(i) = ovgivenrate(i);
            }
            if (!NumericVector::is_na(ovoffset(i))) {
              if (ovoffset(i) != 0.) ovnewoffset(i) = ovoffset(i);
            }
            if (NumericVector::is_na(ovmultiplier(i))) {
              ovmultiplier(i) = 1.0;
            }
            ovnewmultiplier(i) = ovmultiplier(i);
          } // i for loop
        } // ovrows if statement
      }
      
      for (int time1 = 0; time1 < nostages_nodead; time1++) {
        for (int time2o = 0; time2o < nostages_nodead; time2o++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            included(currentindex) = 1.0;
            
            stageorder3(currentindex) = stageorder(time3);
            stageorder2n(currentindex) = stageorder(time2o);
            stageorder2o(currentindex) = stageorder(time2o);
            stageorder1(currentindex) = stageorder(time1);
                  
            stage3(currentindex) = newstageid(time3);
            stage2n(currentindex) = newstageid(time2o);
            stage2o(currentindex) = newstageid(time2o);
            stage1(currentindex) = newstageid(time1);
            
            size3(currentindex) = binsizectr(time3);
            size2n(currentindex) = binsizectr(time2o);
            size2o(currentindex) = binsizectr(time2o);
            size1(currentindex) = binsizectr(time1);
            
            sizeb3(currentindex) = binsizebctr(time3);
            sizeb2n(currentindex) = binsizebctr(time2o);
            sizeb2o(currentindex) = binsizebctr(time2o);
            sizeb1(currentindex) = binsizebctr(time1);
            
            if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0.0;
                  
            sizec3(currentindex) = binsizecctr(time3);
            sizec2n(currentindex) = binsizecctr(time2o);
            sizec2o(currentindex) = binsizecctr(time2o);
            sizec1(currentindex) = binsizecctr(time1);
            
            if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
            if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
            if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
            if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0.0;
            
            obs3(currentindex) = obsstatus(time3);
            obs2n(currentindex) = obsstatus(time2o);
            obs2o(currentindex) = obsstatus(time2o);
            obs1(currentindex) = obsstatus(time1);
            
            rep3(currentindex) = repstatus(time3);
            rep2n(currentindex) = repstatus(time2o);
            rep2o(currentindex) = repstatus(time2o);
            rep1(currentindex) = repstatus(time1);
            
            mat3(currentindex) = matstatus(time3);
            mat2n(currentindex) = matstatus(time2o);
            mat2o(currentindex) = matstatus(time2o);
            mat1(currentindex) = matstatus(time1);
            
            imm3(currentindex) = immstatus(time3);
            imm2n(currentindex) = immstatus(time2o);
            imm2o(currentindex) = immstatus(time2o);
            imm1(currentindex) = immstatus(time1);
            
            // Determine repentry3 on basis of input repmatrix
            if (time3 < nostages_nodead_nounborn) {
              if (repmattype == 1) {
                repm_elem = time3 + (time2o * nostages_nodead_nounborn);
              } else if (repmattype == 2) {
                repm_elem = time3 + (time2o * nostages_nodead_nounborn) + 
                  (time2o * nostages_nodead_nounborn * nostages_nodead_nounborn) +
                  (time1 * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                    nostages_nodead_nounborn);
              } else {
                repm_elem = -1;
              }
            }
            
            if(repm_elem > -1) {
              if (repmatrix(repm_elem) > 0.0) {
                repentry3(currentindex) = repmatrix(repm_elem);
              }
            }
            
            if (time3 < nostages_nodead_nounborn) {
              if (repmattype == 1) { // Ahistorical repmatrix
                repentry3(currentindex) = repmatrix((time3 + (nostages_nodead_nounborn * time2o)));
              } else if (repmattype == 2) {  // Historical repmatrix
                repentry3(currentindex) = repmatrix((time3 + (nostages_nodead_nounborn * time2o)) + 
                  ((nostages_nodead_nounborn * nostages_nodead_nounborn * time2o)) +
                  (nostages_nodead_nounborn * nostages_nodead_nounborn * 
                    nostages_nodead_nounborn * time1));
              }
            } else {
              repentry3(currentindex) = 0.0;
            }
            
            indata3(currentindex) = indata(time3);
            indata2n(currentindex) = indata(time2o);
            indata2o(currentindex) = indata(time2o);
            indata1(currentindex) = indata(time1);
            
            binwidth(currentindex) = binsizewidth(time3);
            binbwidth(currentindex) = binsizebwidth(time3);
            bincwidth(currentindex) = binsizecwidth(time3);
            
            if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
            if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
            
            minage3(currentindex) = minage(time3);
            minage2(currentindex) = minage(time2o);
            maxage3(currentindex) = maxage(time3);
            maxage2(currentindex) = maxage(time2o);
            actualage(currentindex) = 0.0;
            
            grp3(currentindex) = group(time3);
            grp2n(currentindex) = group(time2o);
            grp2o(currentindex) = group(time2o);
            grp1(currentindex) = group(time1);
            
            if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
              deadandnasty = 1.0;
            } else if (stage2o(currentindex) == nostages || stage1(currentindex) == nostages) {
              deadandnasty = 1.0;
            } else {
              deadandnasty = 0.0;
            }
            
            // Required for proper fecundity estimation in rlefko3
            index321d(currentindex) = (stage3(currentindex) - 1) + 
              ((stage2n(currentindex) - 1) * nostages_nodead_nounborn) + 
              ((stage2n(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
              ((stage1(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                nostages_nodead_nounborn);
            
            if (deadandnasty == 0.0) {
              aliveequal(currentindex) = (stageorder3(currentindex) - 1) + ((stageorder2n(currentindex) - 1) * 
                  (nostages - 1)) + ((stageorder2o(currentindex) - 1) * (nostages - 1) * (nostages - 1)) + 
                ((stageorder1(currentindex) - 1) * (nostages - 1) * (nostages - 1) * (nostages - 1));
              
              index321(currentindex) = (stage3(currentindex) - 1) + 
                ((stage2n(currentindex) - 1) * nostages_nodead_nounborn) + 
                ((stage2n(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
                ((stage1(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                  nostages_nodead_nounborn);
              index21(currentindex) = (stage2n(currentindex) - 1) + ((stage1(currentindex) - 1) * nostages);
            }
            
            indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
              indata2o(currentindex) * indata1(currentindex);
            
            currentindex += 1;
          } // time3 loop
        } // time2o loop
      } // time1 loop 
      
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          asadditions = LefkoMats::ovreplace(index321, ovindexold321,
            ovindexnew321, ovconvtype, ovnew3, ovnewgivenrate, ovnewoffset,
            ovnewmultiplier);
          
          ovgivent = asadditions.col(0);
          ovestt = asadditions.col(1);
          ovgivenf = asadditions.col(2);
          ovestf = asadditions.col(3);
          ovrepentry = asadditions.col(4);
          ovsurvmult = asadditions.col(5);
          ovfecmult = asadditions.col(6);
          ovoffsett = asadditions.col(7);
          ovoffsetf = asadditions.col(8);
          
          arma::uvec workedupindex = find(ovrepentry > 0.0);
          int changedreps = static_cast<int>(workedupindex.n_elem);
          
          if (changedreps > 0) {
            for (int i = 0; i < changedreps; i++) {
              repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
            }
          }
        } // ovreplace if statement
      }
    } else if (style == 1) { // Ahistorical case
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
          
            ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages);
            ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages);
            
            if (ovindexold321(i) < 0) ovindexold321(i) = -1.0;
            if (ovindexnew321(i) < 0) ovindexnew321(i) = -1.0;
            
            if (!NumericVector::is_na(ovgivenrate(i))) {
              ovnewgivenrate(i) = ovgivenrate(i);
            }
            if (!NumericVector::is_na(ovoffset(i))) {
              if (ovoffset(i) != 0.) ovnewoffset(i) = ovoffset(i);
            }
            if (NumericVector::is_na(ovmultiplier(i))) {
              ovmultiplier(i) = 1;
            }
            ovnewmultiplier(i) = ovmultiplier(i);
          } // i for loop
        } // ovrows if statement
      }
      
      for (int time2n = 0; time2n < nostages_nodead; time2n++) {
        for (int time3 = 0; time3 < nostages; time3++) {
          stageorder3(currentindex) = stageorder(time3);
          stageorder2n(currentindex) = stageorder(time2n);
          stageorder2o(currentindex) = stageorder(time2n);
          stageorder1(currentindex) = 0;
                  
          stage3(currentindex) = newstageid(time3);
          stage2n(currentindex) = newstageid(time2n);
          stage2o(currentindex) = newstageid(time2n);
          stage1(currentindex) = 0;
          
          size3(currentindex) = binsizectr(time3);
          size2n(currentindex) = binsizectr(time2n);
          size2o(currentindex) = binsizectr(time2n);
          size1(currentindex) = 0.0;
          
          sizeb3(currentindex) = binsizebctr(time3);
          sizeb2n(currentindex) = binsizebctr(time2n);
          sizeb2o(currentindex) = binsizebctr(time2n);
          sizeb1(currentindex) = 0.0;
          
          if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
          if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
          if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
          
          sizec3(currentindex) = binsizecctr(time3);
          sizec2n(currentindex) = binsizecctr(time2n);
          sizec2o(currentindex) = binsizecctr(time2n);
          sizec1(currentindex) = 0.0;
          
          if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
          if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
          if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
          
          obs3(currentindex) = obsstatus(time3);
          obs2n(currentindex) = obsstatus(time2n);
          obs2o(currentindex) = obsstatus(time2n);
          obs1(currentindex) = 0.0;
          
          rep3(currentindex) = repstatus(time3);
          rep2n(currentindex) = repstatus(time2n);
          rep2o(currentindex) = repstatus(time2n);
          rep1(currentindex) = 0.0;
          
          mat3(currentindex) = matstatus(time3);
          mat2n(currentindex) = matstatus(time2n);
          mat2o(currentindex) = matstatus(time2n);
          mat1(currentindex) = 0.0;
          
          imm3(currentindex) = immstatus(time3);
          imm2n(currentindex) = immstatus(time2n);
          imm2o(currentindex) = immstatus(time2n);
          imm1(currentindex) = 0.0;
          
          if (time3 < nostages_nodead) {
            repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
          } else {
            repentry3(currentindex) = 0.0;
          }
          
          indata3(currentindex) = indata(time3);
          indata2n(currentindex) = indata(time2n);
          indata2o(currentindex) = indata(time2n);
          indata1(currentindex) = 1.0;
          
          binwidth(currentindex) = binsizewidth(time3);
          binbwidth(currentindex) = binsizebwidth(time3);
          bincwidth(currentindex) = binsizecwidth(time3);
          
          if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
          if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
          
          minage3(currentindex) = minage(time3);
          minage2(currentindex) = minage(time2n);
          maxage3(currentindex) = maxage(time3);
          maxage2(currentindex) = maxage(time2n);
          actualage(currentindex) = 0.0;
          
          grp3(currentindex) = group(time3);
          grp2n(currentindex) = group(time2n);
          grp2o(currentindex) = group(time2n);
          grp1(currentindex) = 0.0;
          
          if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
            deadandnasty = 1.0;
          } else {
            deadandnasty = 0.0;
          }
          
          if (deadandnasty == 0.0) {
            aliveequal(currentindex) = (stageorder3(currentindex) - 1) + 
              ((stageorder2n(currentindex) - 1) * nostages_nodead);
            
            index321(currentindex) = (stage3(currentindex) - 1) + 
              ((stage2n(currentindex) - 1) * nostages);
            index21(currentindex) = (stage2n(currentindex) - 1);
          }
          
          indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
            indata2o(currentindex);
            
          currentindex += 1;
        } // time3 loop
      } // time2n loop
      
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          asadditions = LefkoMats::ovreplace(index321, ovindexold321, ovindexnew321,
            ovconvtype, ovnew3, ovnewgivenrate, ovnewoffset, ovnewmultiplier);
          
          ovgivent = asadditions.col(0);
          ovestt = asadditions.col(1);
          ovgivenf = asadditions.col(2);
          ovestf = asadditions.col(3);
          ovrepentry = asadditions.col(4);
          ovsurvmult = asadditions.col(5);
          ovfecmult = asadditions.col(6);
          ovoffsett = asadditions.col(7);
          ovoffsetf = asadditions.col(8);
          
          arma::uvec workedupindex = find(ovrepentry > 0.0);
          int changedreps = static_cast<int>(workedupindex.n_elem);
          
          if (changedreps > 0) {
            for (int i = 0; i < changedreps; i++) {
              repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
            }
          }
        } // ovreplace if statement
      }
    } else if (style == 2) { // Age-by-stage case
      int age3 {firstage};
      
      for (int time3 = 0; time3 < nostages; time3++) {
        if (NumericVector::is_na(maxage(time3))) {
          maxage(time3) = finalage + cont;
        }
      }
      
      // Sets up overwrite tables
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
          // First set of loops establishes a number of indices
          for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
            int age2 = ovage2(i);
            
            for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
              ovconvtypeage(i) = ovconvtype(i);
                
              if (age2 < totalages) {
                if (ovconvtype(i) == 1.0) {
                  age3 = age2 + 1;
                } else {
                  age3 = firstage;
                }
                
                if (ovstage3(i) == origstageid(j)) {
                  ovindex3(i) = j; // newstageid(j) - 1.0
                }
                
                if (ovstage2(i) == origstageid(j)) {
                  ovindex2(i) = j; // newstageid(j) - 1.0
                }
                
                if (oveststage3(i) == origstageid(j)) {
                  ovnew3(i) = j; // newstageid(j) - 1.0
                }
                
                if (oveststage2(i) == origstageid(j)) {
                  ovnew2(i) = j; // newstageid(j) - 1.0
                }
                
                if (ovindex3(i) != -1.0 && ovindex2(i) != -1.0) {
                  ovindexold321(i) = ovindex3(i) + ((age3 - firstage) * nostages) +
                    (ovindex2(i) * nostages * totalages) + 
                    ((age2 - firstage) * nostages * nostages * totalages);
                }
                
                if (ovnew3(i) != -1.0 && ovnew2(i) != -1.0) {
                  if (!IntegerVector::is_na(ovestage2(i)) && ovestage2(i) != -1) {
                    int newage2 = ovestage2(i);
                    int newage3 = newage2 + 1;
                    
                    ovindexnew321(i) = ovnew3(i) + ((newage3 - firstage) * nostages) +
                      (ovnew2(i) * nostages * totalages) +
                      ((newage2 - firstage) * nostages * nostages * totalages);
                  } else {
                    ovindexnew321(i) = ovnew3(i) + ((age3 - firstage) * nostages) +
                      (ovnew2(i) * nostages * totalages) +
                      ((age2 - firstage) * nostages * nostages * totalages);
                  }
                }
                
                if (!NumericVector::is_na(ovgivenrate(i))) {
                  ovnewgivenrate(i) = ovgivenrate(i);
                }
                if (!NumericVector::is_na(ovoffset(i))) {
                  if (ovoffset(i) != 0.) ovnewoffset(i) = ovoffset(i);
                }
                if (NumericVector::is_na(ovmultiplier(i))) ovmultiplier(i) = 1.0;
                
                ovnewmultiplier(i) = ovmultiplier(i);
              } else {
                if (ovconvtype(i) == 1.0) {
                  age3 = age2;
                } else {
                  age3 = firstage;
                }
                
                if (ovstage3(i) == origstageid(j)) {
                  ovindex3(i) = j; // newstageid(j) - 1.0
                }
                
                if (ovstage2(i) == origstageid(j)) {
                  ovindex2(i) = j; // newstageid(j) - 1.0
                }
                
                if (oveststage3(i) == origstageid(j)) {
                  ovnew3(i) = j; // newstageid(j) - 1.0
                }
                
                if (oveststage2(i) == origstageid(j)) {
                  ovnew2(i) = j; // newstageid(j) - 1.0
                }
                
                if (ovindex3(i) != -1.0 && ovindex2(i) != -1.0) {
                  ovindexold321(i) = ovindex3(i) + ((age3 - firstage) * nostages) +
                    (ovindex2(i) * nostages * totalages) +
                    ((age2 - firstage) * nostages * nostages * totalages);
                }
                
                if (ovnew3(i) != -1.0 && ovnew2(i) != -1.0) {
                  if (!IntegerVector::is_na(ovestage2(i)) && ovestage2(i) != -1) {
                    int newage2 = ovestage2(i);
                    int newage3 = newage2 + 1;
                    
                    ovindexnew321(i) = ovnew3(i) + ((newage3 - firstage) * nostages) +
                      (ovnew2(i) * nostages * totalages) +
                      ((newage2 - firstage) * nostages * nostages * totalages);
                  } else {
                    ovindexnew321(i) = ovnew3(i) + ((age3 - firstage) * nostages) +
                      (ovnew2(i) * nostages * totalages) +
                      ((age2 - firstage) * nostages * nostages * totalages);
                  }
                }
                if (!NumericVector::is_na(ovgivenrate(i))) {
                  ovnewgivenrate(i) = ovgivenrate(i);
                }
                if (!NumericVector::is_na(ovoffset(i))) {
                  if (ovoffset(i) != 0.) ovnewoffset(i) = ovoffset(i);
                }
                if (NumericVector::is_na(ovmultiplier(i))) ovmultiplier(i) = 1.0;
                
                ovnewmultiplier(i) = ovmultiplier(i);
              }
            } // j for loop
            
          if (ovindexold321(i) < 0) ovindexold321(i) = -1.0;
          if (ovindexnew321(i) < 0) ovindexnew321(i) = -1.0;
            
          } // i for loop
        } // ovrows if statement
      }
      for (int age2 = firstage; age2 <= finalage; age2++) {
        if (age2 < finalage) { // First loop takes care of age transitions
          for (int time2n = 0; time2n < nostages; time2n++) {
            for (int time3 = 0; time3 < nostages; time3++) {
              
              // First survival
              age3 = age2 + 1;
              currentindex = time3 + ((age3 - firstage) * nostages) + 
                (time2n * nostages * totalages) +
                ((age2 - firstage) * nostages * nostages * totalages);
              
              stage3(currentindex) = newstageid(time3);
              stage2n(currentindex) = newstageid(time2n);
              stage2o(currentindex) = newstageid(time2n);
              stage1(currentindex) = 0;
              
              size3(currentindex) = binsizectr(time3);
              size2n(currentindex) = binsizectr(time2n);
              size2o(currentindex) = binsizectr(time2n);
              size1(currentindex) = 0.0;
              
              sizeb3(currentindex) = binsizebctr(time3);
              sizeb2n(currentindex) = binsizebctr(time2n);
              sizeb2o(currentindex) = binsizebctr(time2n);
              sizeb1(currentindex) = 0.0;
              
              if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
              if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
              if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
              
              sizec3(currentindex) = binsizecctr(time3);
              sizec2n(currentindex) = binsizecctr(time2n);
              sizec2o(currentindex) = binsizecctr(time2n);
              sizec1(currentindex) = 0.0;
              
              if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
              if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
              if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
              
              obs3(currentindex) = obsstatus(time3);
              obs2n(currentindex) = obsstatus(time2n);
              obs2o(currentindex) = obsstatus(time2n);
              obs1(currentindex) = 0.0;
              
              rep3(currentindex) = repstatus(time3);
              rep2n(currentindex) = repstatus(time2n);
              rep2o(currentindex) = repstatus(time2n);
              rep1(currentindex) = 0.0;
              
              mat3(currentindex) = matstatus(time3);
              mat2n(currentindex) = matstatus(time2n);
              mat2o(currentindex) = matstatus(time2n);
              mat1(currentindex) = 0.0;
              
              imm3(currentindex) = immstatus(time3);
              imm2n(currentindex) = immstatus(time2n);
              imm2o(currentindex) = immstatus(time2n);
              imm1(currentindex) = 0.0;
              
              repentry3(currentindex) = 0.0;
              
              indata3(currentindex) = indata(time3);
              indata2n(currentindex) = indata(time2n);
              indata2o(currentindex) = indata(time2n);
              indata1(currentindex) = 0.0;
              
              binwidth(currentindex) = binsizewidth(time3);
              binbwidth(currentindex) = binsizebwidth(time3);
              bincwidth(currentindex) = binsizecwidth(time3);
              
              if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
              if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
              
              minage3(currentindex) = minage(time3);
              minage2(currentindex) = minage(time2n);
              maxage3(currentindex) = maxage(time3);
              maxage2(currentindex) = maxage(time2n);
              actualage(currentindex) = age2;
              
              grp3(currentindex) = group(time3);
              grp2n(currentindex) = group(time2n);
              grp2o(currentindex) = group(time2n);
              grp1(currentindex) = 0.0;
              
              // Indexer order: (1st # age blocks) + (1st # stage cols) +
              // (1st # age rows) + stage in time 3
              index321(currentindex) = currentindex;
              index21(currentindex) = time2n + ((age2 - firstage) * nostages);
              indatalong(currentindex) = 1.0;
              
              // Identify elements with non-zero entries by element number in final matrix
              if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
                if (age2 >= minage2(currentindex) && age2 < maxage3(currentindex)) { 
                  
                  // Survival transitions
                  aliveequal(currentindex) =
                    ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                    (time2n * (nostages - 1) * totalages) +
                    ((age3 - firstage) * (nostages - 1)) + time3;
                }
              }
              
              if (time3 < nostages_nodead && time2n < nostages_nodead) {
                
                if (repmatrix((time3 + (nostages_nodead * time2n))) > 0.0) {
                  
                  // Now fecundity
                  age3 = firstage;
                  currentindex = time3 + ((age3 - firstage) * nostages) + 
                    (time2n * nostages * totalages) +
                    ((age2 - firstage) * nostages * nostages * totalages);
                  
                  stage3(currentindex) = newstageid(time3);
                  stage2n(currentindex) = newstageid(time2n);
                  stage2o(currentindex) = newstageid(time2n);
                  stage1(currentindex) = 0.0;
                  
                  size3(currentindex) = binsizectr(time3);
                  size2n(currentindex) = binsizectr(time2n);
                  size2o(currentindex) = binsizectr(time2n);
                  size1(currentindex) = 0.0;
                  
                  sizeb3(currentindex) = binsizebctr(time3);
                  sizeb2n(currentindex) = binsizebctr(time2n);
                  sizeb2o(currentindex) = binsizebctr(time2n);
                  sizeb1(currentindex) = 0.0;
                  
                  if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
                  
                  sizec3(currentindex) = binsizecctr(time3);
                  sizec2n(currentindex) = binsizecctr(time2n);
                  sizec2o(currentindex) = binsizecctr(time2n);
                  sizec1(currentindex) = 0.0;
                  
                  if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                  
                  obs3(currentindex) = obsstatus(time3);
                  obs2n(currentindex) = obsstatus(time2n);
                  obs2o(currentindex) = obsstatus(time2n);
                  obs1(currentindex) = 0.0;
                  
                  rep3(currentindex) = repstatus(time3);
                  rep2n(currentindex) = repstatus(time2n);
                  rep2o(currentindex) = repstatus(time2n);
                  rep1(currentindex) = 0.0;
                  
                  mat3(currentindex) = matstatus(time3);
                  mat2n(currentindex) = matstatus(time2n);
                  mat2o(currentindex) = matstatus(time2n);
                  mat1(currentindex) = 0.0;
                  
                  imm3(currentindex) = immstatus(time3);
                  imm2n(currentindex) = immstatus(time2n);
                  imm2o(currentindex) = immstatus(time2n);
                  imm1(currentindex) = 0.0;
                  
                  if (rep2n(currentindex) > 0.0 && time3 < nostages_nodead && time2n < nostages_nodead) {
                    repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
                  } else repentry3(currentindex) = 0.0;
                  
                  indata3(currentindex) = indata(time3);
                  indata2n(currentindex) = indata(time2n);
                  indata2o(currentindex) = indata(time2n);
                  indata1(currentindex) = 0.0;
                  
                  binwidth(currentindex) = binsizewidth(time3);
                  binbwidth(currentindex) = binsizebwidth(time3);
                  bincwidth(currentindex) = binsizecwidth(time3);
                  
                  if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
                  if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                  
                  minage3(currentindex) = minage(time3);
                  minage2(currentindex) = minage(time2n);
                  maxage3(currentindex) = maxage(time3);
                  maxage2(currentindex) = maxage(time2n);
                  actualage(currentindex) = age2;
                  
                  grp3(currentindex) = group(time3);
                  grp2n(currentindex) = group(time2n);
                  grp2o(currentindex) = group(time2n);
                  grp1(currentindex) = 0.0;
                  
                  // The next indexer includes the following order: (1st # of age blocks) + 
                  // (1st # of stage cols) + (1st # of age rows) + stage in time 3
                  index321(currentindex) = currentindex;
                  index21(currentindex) = time2n + ((age2 - firstage) * nostages);
                  indatalong(currentindex) = 1.0;
                  
                  // This section identifies elements with non-zero entries by their
                  // element number in the final matrix
                  if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
                    if (age2 >= minage2(currentindex) && age2 <= maxage2(currentindex)) { 
                      
                      // Fecundity transitions
                      aliveequal(currentindex) = 
                        ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                        (time2n * (nostages - 1) * totalages) +
                        ((age3 - firstage) * (nostages - 1)) + time3;
                    }
                  } // if statement leading to aliveequal assignment
                } // if statement yielding fecundity estimation
              } // if statement checking time3 and time2n
            } // time3 loop
          } // time2n loop
        } else if (cont == 1) { // Self-loop on final age, if organism can live past final age
          for (int time2n = 0; time2n < nostages; time2n++) {
            for (int time3 = 0; time3 < nostages; time3++) {
              
              // First survival
              age3 = age2;
              currentindex = time3 + ((age3 - firstage) * nostages) + 
                (time2n * nostages * totalages) +
                ((age2 - firstage) * nostages * nostages * totalages);
              
              stage3(currentindex) = newstageid(time3);
              stage2n(currentindex) = newstageid(time2n);
              stage2o(currentindex) = newstageid(time2n);
              stage1(currentindex) = 0.0;
              
              size3(currentindex) = binsizectr(time3);
              size2n(currentindex) = binsizectr(time2n);
              size2o(currentindex) = binsizectr(time2n);
              size1(currentindex) = 0.0;
              
              sizeb3(currentindex) = binsizebctr(time3);
              sizeb2n(currentindex) = binsizebctr(time2n);
              sizeb2o(currentindex) = binsizebctr(time2n);
              sizeb1(currentindex) = 0.0;
              
              if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
              if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
              if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
              
              sizec3(currentindex) = binsizecctr(time3);
              sizec2n(currentindex) = binsizecctr(time2n);
              sizec2o(currentindex) = binsizecctr(time2n);
              sizec1(currentindex) = 0.0;
              
              if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
              if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
              if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                  
              obs3(currentindex) = obsstatus(time3);
              obs2n(currentindex) = obsstatus(time2n);
              obs2o(currentindex) = obsstatus(time2n);
              obs1(currentindex) = 0.0;
              
              rep3(currentindex) = repstatus(time3);
              rep2n(currentindex) = repstatus(time2n);
              rep2o(currentindex) = repstatus(time2n);
              rep1(currentindex) = 0.0;
              
              mat3(currentindex) = matstatus(time3);
              mat2n(currentindex) = matstatus(time2n);
              mat2o(currentindex) = matstatus(time2n);
              mat1(currentindex) = 0.0;
              
              imm3(currentindex) = immstatus(time3);
              imm2n(currentindex) = immstatus(time2n);
              imm2o(currentindex) = immstatus(time2n);
              imm1(currentindex) = 0.0;
              
              repentry3(currentindex) = 0.0;
              
              indata3(currentindex) = indata(time3);
              indata2n(currentindex) = indata(time2n);
              indata2o(currentindex) = indata(time2n);
              indata1(currentindex) = 0.0;
              
              binwidth(currentindex) = binsizewidth(time3);
              binbwidth(currentindex) = binsizebwidth(time3);
              bincwidth(currentindex) = binsizecwidth(time3);
              
              if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
              if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                  
              minage3(currentindex) = minage(time3);
              minage2(currentindex) = minage(time2n);
              maxage3(currentindex) = maxage(time3);
              maxage2(currentindex) = maxage(time2n);
              actualage(currentindex) = age2;
              
              grp3(currentindex) = group(time3);
              grp2n(currentindex) = group(time2n);
              grp2o(currentindex) = group(time2n);
              grp1(currentindex) = 0.0;
              
              // Indexer order: (1st # age blocks) + (1st # stage cols) +
              // (1st # age rows) + stage in time 3
              index321(currentindex) = currentindex;
              index21(currentindex) = time2n + ((age2 - firstage) * nostages);
              indatalong(currentindex) = 1;
              
              // Identify elements with non-zero entries by element number in final matrix
              if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
                if (age2 >= minage2(currentindex) && age2 < maxage3(currentindex)) { 
  
                  // Survival transitions
                  aliveequal(currentindex) = 
                    ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                    (time2n * (nostages - 1) * totalages) +
                    ((age3 - firstage) * (nostages - 1)) + time3;
                }
              }
              
              if (time3 < nostages_nodead && time2n < nostages_nodead) {
                if (repmatrix((time3 + (nostages_nodead * time2n))) > 0.0) {
                  
                  // Now fecundity
                  age3 = firstage;
                  currentindex = time3 + ((age3 - firstage) * nostages) + 
                    (time2n * nostages * totalages) +
                    ((age2 - firstage) * nostages * nostages * totalages);
                  
                  stage3(currentindex) = newstageid(time3);
                  stage2n(currentindex) = newstageid(time2n);
                  stage2o(currentindex) = newstageid(time2n);
                  stage1(currentindex) = 0.0;
                  
                  size3(currentindex) = binsizectr(time3);
                  size2n(currentindex) = binsizectr(time2n);
                  size2o(currentindex) = binsizectr(time2n);
                  size1(currentindex) = 0.0;
                  
                  sizeb3(currentindex) = binsizebctr(time3);
                  sizeb2n(currentindex) = binsizebctr(time2n);
                  sizeb2o(currentindex) = binsizebctr(time2n);
                  sizeb1(currentindex) = 0.0;
                  
                  if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
                  if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
                  
                  sizec3(currentindex) = binsizecctr(time3);
                  sizec2n(currentindex) = binsizecctr(time2n);
                  sizec2o(currentindex) = binsizecctr(time2n);
                  sizec1(currentindex) = 0.0;
                  
                  if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
                  if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                  
                  obs3(currentindex) = obsstatus(time3);
                  obs2n(currentindex) = obsstatus(time2n);
                  obs2o(currentindex) = obsstatus(time2n);
                  obs1(currentindex) = 0.0;
                  
                  rep3(currentindex) = repstatus(time3);
                  rep2n(currentindex) = repstatus(time2n);
                  rep2o(currentindex) = repstatus(time2n);
                  rep1(currentindex) = 0.0;
                  
                  mat3(currentindex) = matstatus(time3);
                  mat2n(currentindex) = matstatus(time2n);
                  mat2o(currentindex) = matstatus(time2n);
                  mat1(currentindex) = 0.0;
                  
                  imm3(currentindex) = immstatus(time3);
                  imm2n(currentindex) = immstatus(time2n);
                  imm2o(currentindex) = immstatus(time2n);
                  imm1(currentindex) = 0.0;
                  
                  if (rep2n(currentindex) == 1) {
                    repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
                  } else repentry3(currentindex) = 0.0;
                  
                  indata3(currentindex) = indata(time3);
                  indata2n(currentindex) = indata(time2n);
                  indata2o(currentindex) = indata(time2n);
                  indata1(currentindex) = 0.0;
                  
                  binwidth(currentindex) = binsizewidth(time3);
                  binbwidth(currentindex) = binsizebwidth(time3);
                  bincwidth(currentindex) = binsizecwidth(time3);
                  
                  if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
                  if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                  
                  minage3(currentindex) = minage(time3);
                  minage2(currentindex) = minage(time2n);
                  maxage3(currentindex) = maxage(time3);
                  maxage2(currentindex) = maxage(time2n);
                  actualage(currentindex) = age2;
                  grp3(currentindex) = group(time3);
                  grp2n(currentindex) = group(time2n);
                  grp2o(currentindex) = group(time2n);
                  grp1(currentindex) = 0.0;
                  
                  // Indexer order: (1st # age blocks) + (1st # stage cols) + 
                  // (1st # age rows) + stage in time 3
                  index321(currentindex) = currentindex;
                  index21(currentindex) = time2n + ((age2 - firstage) * nostages);
                  indatalong(currentindex) = 1.0;
                  
                  // Identify elements with non-zero entries by element number in final matrix
                  if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
                    if (age2 >= minage2(currentindex) && age2 <= maxage2(currentindex)) { 
                      
                      // Fecundity transitions
                      aliveequal(currentindex) =
                        ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                        (time2n * (nostages - 1) * totalages) +
                        ((age3 - firstage) * (nostages - 1)) + time3;
                    }
                  } // if statement leading to aliveequal assignment
                } // if statement yielding fecundity estimation
              } // if statement checking time3 and time2n
            } // time3 loop
          } // time2n loop
        }// if-else statement
      } // age2 loop
      
      if (ovrows > 0) {
        if (ovrows > 1 || ovconvtype(0) != -1.0) {
        
          asadditions = LefkoMats::ovreplace(index321, ovindexold321, ovindexnew321,
            ovconvtypeage, ovnew3, ovnewgivenrate, ovnewoffset, ovnewmultiplier);
          
          ovgivent = asadditions.col(0);
          ovestt = asadditions.col(1);
          ovgivenf = asadditions.col(2);
          ovestf = asadditions.col(3);
          ovrepentry = asadditions.col(4);
          ovsurvmult = asadditions.col(5);
          ovfecmult = asadditions.col(6);
          ovoffsett = asadditions.col(7);
          ovoffsetf = asadditions.col(8);
          
          arma::uvec workedupindex = find(ovrepentry > 0.0);
          int changedreps = static_cast<int>(workedupindex.n_elem);
          
          if (changedreps > 0) {
            for (int i = 0; i < changedreps; i++) {
              repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
            }
          }
        } // ovreplace if statement
      }
    } // Age-by-stage loop (style = 2)
    
    // Output formatting
    Rcpp::List output_longlist(62);
    int stage3_length = 0;
    
    arma::uvec used_indices;
    
    if (filter == 1) {
      used_indices = find(index321 != -1.0);
    } else if (filter == 2) {
      used_indices = find(aliveequal != -1.0);
    }
    
    if (filter > 0) {
      int new_length = static_cast<int>(used_indices.n_elem);
      stage3_length = new_length;
      
      IntegerVector stage3_new(new_length);
      IntegerVector stage2n_new(new_length);
      IntegerVector stage2o_new(new_length);
      IntegerVector stage1_new(new_length);
      
      NumericVector size3_new(new_length);
      NumericVector size2n_new(new_length);
      NumericVector size2o_new(new_length);
      NumericVector size1_new(new_length);
      
      NumericVector sizeb3_new(new_length);
      NumericVector sizeb2n_new(new_length);
      NumericVector sizeb2o_new(new_length);
      NumericVector sizeb1_new(new_length);
      
      NumericVector sizec3_new(new_length);
      NumericVector sizec2n_new(new_length);
      NumericVector sizec2o_new(new_length);
      NumericVector sizec1_new(new_length);
      
      NumericVector obs3_new(new_length);
      NumericVector obs2n_new(new_length);
      NumericVector obs2o_new(new_length);
      NumericVector obs1_new(new_length);
      
      NumericVector rep3_new(new_length);
      NumericVector rep2n_new(new_length);
      NumericVector rep2o_new(new_length);
      NumericVector rep1_new(new_length);
      
      NumericVector mat3_new(new_length);
      NumericVector mat2n_new(new_length);
      NumericVector mat2o_new(new_length);
      NumericVector mat1_new(new_length);
      
      NumericVector imm3_new(new_length);
      NumericVector imm2n_new(new_length);
      NumericVector imm2o_new(new_length);
      NumericVector imm1_new(new_length);
      
      NumericVector repentry3_new(new_length);
      NumericVector indata3_new(new_length);
      NumericVector indata2n_new(new_length);
      NumericVector indata2o_new(new_length);
    
      NumericVector indata1_new(new_length);
      NumericVector binwidth_new(new_length);
      NumericVector binbwidth_new(new_length);
      NumericVector bincwidth_new(new_length);
      
      NumericVector minage3_new(new_length);
      NumericVector minage2_new(new_length);
      NumericVector maxage3_new(new_length);
      NumericVector maxage2_new(new_length);
      NumericVector actualage_new(new_length);
      
      NumericVector grp3_new(new_length);
      NumericVector grp2n_new(new_length);
      NumericVector grp2o_new(new_length);
      NumericVector grp1_new(new_length);
      
      NumericVector indatalong_new(new_length);
      NumericVector ovgivent_new(new_length);
      NumericVector ovgivenf_new(new_length);
      
      NumericVector ovoffsett_new(new_length);
      NumericVector ovoffsetf_new(new_length);
      
      NumericVector ovestt_new(new_length);
      NumericVector ovestf_new(new_length);
      
      NumericVector ovsurvmult_new(new_length);
      NumericVector ovfecmult_new(new_length);
      
      NumericVector aliveequal_new(new_length);
      NumericVector index321_new(new_length);
      NumericVector index321d_new(new_length);
      NumericVector index21_new(new_length);
      
      for (int i = 0; i < new_length; i++) {
        stage3_new(i) = stage3(used_indices(i));
        stage2n_new(i) = stage2n(used_indices(i));
        stage2o_new(i) = stage2o(used_indices(i));
        stage1_new(i) = stage1(used_indices(i));
        
        size3_new(i) = size3(used_indices(i));
        size2n_new(i) = size2n(used_indices(i));
        size2o_new(i) = size2o(used_indices(i));
        size1_new(i) = size1(used_indices(i));
        
        sizeb3_new(i) = sizeb3(used_indices(i));
        sizeb2n_new(i) = sizeb2n(used_indices(i));
        sizeb2o_new(i) = sizeb2o(used_indices(i));
        sizeb1_new(i) = sizeb1(used_indices(i));
        
        sizec3_new(i) = sizec3(used_indices(i));
        sizec2n_new(i) = sizec2n(used_indices(i));
        sizec2o_new(i) = sizec2o(used_indices(i));
        sizec1_new(i) = sizec1(used_indices(i));
        
        obs3_new(i) = obs3(used_indices(i));
        obs2n_new(i) = obs2n(used_indices(i));
        obs2o_new(i) = obs2o(used_indices(i));
        obs1_new(i) = obs1(used_indices(i));
        
        rep3_new(i) = rep3(used_indices(i));
        rep2n_new(i) = rep2n(used_indices(i));
        rep2o_new(i) = rep2o(used_indices(i));
        rep1_new(i) = rep1(used_indices(i));
        
        mat3_new(i) = mat3(used_indices(i));
        mat2n_new(i) = mat2n(used_indices(i));
        mat2o_new(i) = mat2o(used_indices(i));
        mat1_new(i) = mat1(used_indices(i));
        
        imm3_new(i) = imm3(used_indices(i));
        imm2n_new(i) = imm2n(used_indices(i));
        imm2o_new(i) = imm2o(used_indices(i));
        imm1_new(i) = imm1(used_indices(i));
        
        repentry3_new(i) = repentry3(used_indices(i));
        indata3_new(i) = indata3(used_indices(i));
        indata2n_new(i) = indata2n(used_indices(i));
        indata2o_new(i) = indata2o(used_indices(i));
      
        indata1_new(i) = indata1(used_indices(i));
        binwidth_new(i) = binwidth(used_indices(i));
        binbwidth_new(i) = binbwidth(used_indices(i));
        bincwidth_new(i) = bincwidth(used_indices(i));
        
        minage3_new(i) = minage3(used_indices(i));
        minage2_new(i) = minage2(used_indices(i));
        maxage3_new(i) = maxage3(used_indices(i));
        maxage2_new(i) = maxage2(used_indices(i));
        actualage_new(i) = actualage(used_indices(i));
        
        grp3_new(i) = grp3(used_indices(i));
        grp2n_new(i) = grp2n(used_indices(i));
        grp2o_new(i) = grp2o(used_indices(i));
        grp1_new(i) = grp1(used_indices(i));
        
        indatalong_new(i) = indatalong(used_indices(i));
        ovgivent_new(i) = ovgivent(used_indices(i));
        ovestt_new(i) = ovestt(used_indices(i));
        ovoffsett_new(i) = ovoffsett(used_indices(i));
        
        ovgivenf_new(i) = ovgivenf(used_indices(i));
        ovestf_new(i) = ovestf(used_indices(i));
        ovoffsetf_new(i) = ovoffsetf(used_indices(i));
        
        ovsurvmult_new(i) = ovsurvmult(used_indices(i));
        ovfecmult_new(i) = ovfecmult(used_indices(i));
        
        aliveequal_new(i) = aliveequal(used_indices(i));
        index321_new(i) = index321(used_indices(i));
        index321d_new(i) = index321d(used_indices(i));
        index21_new(i) = index21(used_indices(i));
      }
      
      output_longlist(0) = stage3_new;
      output_longlist(1) = stage2n_new;
      output_longlist(2) = stage2o_new;
      output_longlist(3) = stage1_new;
      output_longlist(4) = size3_new;
      output_longlist(5) = size2n_new;
      output_longlist(6) = size2o_new;
      output_longlist(7) = size1_new;
      output_longlist(8) = sizeb3_new;
      output_longlist(9) = sizeb2n_new;
      
      output_longlist(10) = sizeb2o_new;
      output_longlist(11) = sizeb1_new;
      output_longlist(12) = sizec3_new;
      output_longlist(13) = sizec2n_new;
      output_longlist(14) = sizec2o_new;
      output_longlist(15) = sizec1_new;
      output_longlist(16) = obs3_new;
      output_longlist(17) = obs2n_new;
      output_longlist(18) = obs2o_new;
      output_longlist(19) = obs1_new;
      
      output_longlist(20) = rep3_new;
      output_longlist(21) = rep2n_new;
      output_longlist(22) = rep2o_new;
      output_longlist(23) = rep1_new;
      output_longlist(24) = mat3_new;
      output_longlist(25) = mat2n_new;
      output_longlist(26) = mat2o_new;
      output_longlist(27) = mat1_new;
      output_longlist(28) = imm3_new;
      output_longlist(29) = imm2n_new;
      
      output_longlist(30) = imm2o_new;
      output_longlist(31) = imm1_new;
      output_longlist(32) = repentry3_new;
      output_longlist(33) = indata3_new;
      output_longlist(34) = indata2n_new;
      output_longlist(35) = indata2o_new;
      output_longlist(36) = indata1_new;
      output_longlist(37) = binwidth_new;
      output_longlist(38) = binbwidth_new;
      output_longlist(39) = bincwidth_new;
      
      output_longlist(40) = minage3_new;
      output_longlist(41) = minage2_new;
      output_longlist(42) = maxage3_new;
      output_longlist(43) = maxage2_new;
      output_longlist(44) = actualage_new;
      
      output_longlist(45) = grp3_new;
      output_longlist(46) = grp2n_new;
      output_longlist(47) = grp2o_new;
      output_longlist(48) = grp1_new;
      
      output_longlist(49) = indatalong_new;
      output_longlist(50) = ovgivent_new;
      output_longlist(51) = ovestt_new;
      output_longlist(52) = ovgivenf_new;
      output_longlist(53) = ovestf_new;
      output_longlist(54) = ovsurvmult_new;
      output_longlist(55) = ovfecmult_new;
      
      output_longlist(56) = aliveequal_new;
      output_longlist(57) = index321_new;
      output_longlist(58) = index321d_new;
      output_longlist(59) = index21_new;
      
      output_longlist(60) = ovoffsett_new;
      output_longlist(61) = ovoffsetf_new;
      
    } else {
      stage3_length = static_cast<int>(stage3.n_elem);
      
      output_longlist(0) = Rcpp::IntegerVector(stage3.begin(), stage3.end());
      output_longlist(1) = Rcpp::IntegerVector(stage2n.begin(), stage2n.end());
      output_longlist(2) = Rcpp::IntegerVector(stage2o.begin(), stage2o.end());
      output_longlist(3) = Rcpp::IntegerVector(stage1.begin(), stage1.end());
      output_longlist(4) = Rcpp::NumericVector(size3.begin(), size3.end());
      output_longlist(5) = Rcpp::NumericVector(size2n.begin(), size2n.end());
      output_longlist(6) = Rcpp::NumericVector(size2o.begin(), size2o.end());
      output_longlist(7) = Rcpp::NumericVector(size1.begin(), size1.end());
      output_longlist(8) = Rcpp::NumericVector(sizeb3.begin(), sizeb3.end());
      output_longlist(9) = Rcpp::NumericVector(sizeb2n.begin(), sizeb2n.end());
      
      output_longlist(10) = Rcpp::NumericVector(sizeb2o.begin(), sizeb2o.end());
      output_longlist(11) = Rcpp::NumericVector(sizeb1.begin(), sizeb1.end());
      output_longlist(12) = Rcpp::NumericVector(sizec3.begin(), sizec3.end());
      output_longlist(13) = Rcpp::NumericVector(sizec2n.begin(), sizec2n.end());
      output_longlist(14) = Rcpp::NumericVector(sizec2o.begin(), sizec2o.end());
      output_longlist(15) = Rcpp::NumericVector(sizec1.begin(), sizec1.end());
      output_longlist(16) = Rcpp::NumericVector(obs3.begin(), obs3.end());
      output_longlist(17) = Rcpp::NumericVector(obs2n.begin(), obs2n.end());
      output_longlist(18) = Rcpp::NumericVector(obs2o.begin(), obs2o.end());
      output_longlist(19) = Rcpp::NumericVector(obs1.begin(), obs1.end());
      
      output_longlist(20) = Rcpp::NumericVector(rep3.begin(), rep3.end());
      output_longlist(21) = Rcpp::NumericVector(rep2n.begin(), rep2n.end());
      output_longlist(22) = Rcpp::NumericVector(rep2o.begin(), rep2o.end());
      output_longlist(23) = Rcpp::NumericVector(rep1.begin(), rep1.end());
      output_longlist(24) = Rcpp::NumericVector(mat3.begin(), mat3.end());
      output_longlist(25) = Rcpp::NumericVector(mat2n.begin(), mat2n.end());
      output_longlist(26) = Rcpp::NumericVector(mat2o.begin(), mat2o.end());
      output_longlist(27) = Rcpp::NumericVector(mat1.begin(), mat1.end());
      output_longlist(28) = Rcpp::NumericVector(imm3.begin(), imm3.end());
      output_longlist(29) = Rcpp::NumericVector(imm2n.begin(), imm2n.end());
      
      output_longlist(30) = Rcpp::NumericVector(imm2o.begin(), imm2o.end());
      output_longlist(31) = Rcpp::NumericVector(imm1.begin(), imm1.end());
      output_longlist(32) = Rcpp::NumericVector(repentry3.begin(), repentry3.end());
      output_longlist(33) = Rcpp::NumericVector(indata3.begin(), indata3.end());
      output_longlist(34) = Rcpp::NumericVector(indata2n.begin(), indata2n.end());
      output_longlist(35) = Rcpp::NumericVector(indata2o.begin(), indata2o.end());
      output_longlist(36) = Rcpp::NumericVector(indata1.begin(), indata1.end());
      output_longlist(37) = Rcpp::NumericVector(binwidth.begin(), binwidth.end());
      output_longlist(38) = Rcpp::NumericVector(binbwidth.begin(), binbwidth.end());
      output_longlist(39) = Rcpp::NumericVector(bincwidth.begin(), bincwidth.end());
      
      output_longlist(40) = Rcpp::NumericVector(minage3.begin(), minage3.end());
      output_longlist(41) = Rcpp::NumericVector(minage2.begin(), minage2.end());
      output_longlist(42) = Rcpp::NumericVector(maxage3.begin(), maxage3.end());
      output_longlist(43) = Rcpp::NumericVector(maxage2.begin(), maxage2.end());
      output_longlist(44) = Rcpp::NumericVector(actualage.begin(), actualage.end());
      
      output_longlist(45) = Rcpp::NumericVector(grp3.begin(), grp3.end());
      output_longlist(46) = Rcpp::NumericVector(grp2n.begin(), grp2n.end());
      output_longlist(47) = Rcpp::NumericVector(grp2o.begin(), grp2o.end());
      output_longlist(48) = Rcpp::NumericVector(grp1.begin(), grp1.end());
      
      output_longlist(49) = Rcpp::NumericVector(indatalong.begin(), indatalong.end());
      output_longlist(50) = Rcpp::NumericVector(ovgivent.begin(), ovgivent.end());
      output_longlist(51) = Rcpp::NumericVector(ovestt.begin(), ovestt.end());
      output_longlist(52) = Rcpp::NumericVector(ovgivenf.begin(), ovgivenf.end());
      output_longlist(53) = Rcpp::NumericVector(ovestf.begin(), ovestf.end());
      output_longlist(54) = Rcpp::NumericVector(ovsurvmult.begin(), ovsurvmult.end());
      output_longlist(55) = Rcpp::NumericVector(ovfecmult.begin(), ovfecmult.end());
      
      output_longlist(56) = Rcpp::NumericVector(aliveequal.begin(), aliveequal.end());
      output_longlist(57) = Rcpp::NumericVector(index321.begin(), index321.end());
      output_longlist(58) = Rcpp::NumericVector(index321d.begin(), index321d.end());
      output_longlist(59) = Rcpp::NumericVector(index21.begin(), index21.end());
      
      output_longlist(60) = Rcpp::NumericVector(ovoffsett.begin(), ovoffsett.end());
      output_longlist(61) = Rcpp::NumericVector(ovoffsetf.begin(), ovoffsetf.end());
    }
    
    CharacterVector namevec = {"stage3", "stage2n", "stage2o", "stage1", "size3",
      "size2n", "size2o", "size1", "sizeb3", "sizeb2n", "sizeb2o", "sizeb1", 
      "sizec3", "sizec2n", "sizec2o", "sizec1", "obs3", "obs2n", "obs2o", "obs1",
      "rep3", "rep2n", "rep2o", "rep1", "mat3", "mat2n", "mat2o", "mat1", "imm3",
      "imm2n", "imm2o", "imm1", "repentry3", "indata3", "indata2n", "indata2o",
      "indata1", "binwidth", "binbwidth", "bincwidth", "minage3", "minage2",
      "maxage3", "maxage2", "actualage", "group3", "group2n", "group2o", "group1",
      "indata", "ovgiven_t", "ovest_t", "ovgiven_f", "ovest_f", "ovsurvmult",
      "ovfecmult", "aliveandequal", "index321", "index321d", "index21",
      "ovoffset_t", "ovoffset_f"};
    output_longlist.attr("names") = namevec;
    output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stage3_length);
    output_longlist.attr("class") = "data.frame";
    
    return output_longlist;
  }
  
  //' Standardize Stageframe For MPM Analysis
  //' 
  //' Function \code{sf_reassess_internal()} takes a stageframe as input, and uses
  //' information supplied there and through the supplement, reproduction and
  //' overwrite tables to rearrange this into a format usable by the matrix
  //' creation functions, \code{mpm_create()}, \code{flefko3()},
  //' \code{flefko2()}, \code{aflefko2()}, \code{rlefko3()}, and \code{rlefko2()}.
  //' 
  //' @name sf_reassess_internal
  //' 
  //' @param stageframe The original stageframe.
  //' @param supplement The original supplemental data input (class
  //' \code{lefkoSD}). Can also equal NA.
  //' @param overwrite An overwrite table.
  //' @param repmatrix The original reproduction matrix. Can also equal \code{NA},
  //' \code{0}, or \code{NULL} (the last value by default).
  //' @param agemat A Boolean value indicating whether MPM is age-by-stage.
  //' @param historical A Boolean value indicating whether MPM is historical.
  //' @param format An integer indicating whether matrices will be in Ehrlen
  //' format (if set to 1), or deVries format (if set to 2). Setting to deVries
  //' format adds one extra stage to account for the prior status of newborns.
  //' @param import_lM A Boolean value indicating whether the lefkoMat object
  //' being dealt with was imported. If \code{TRUE}, then will skip warnings
  //' related to the lack of supplements being input (since supplements cannot
  //' be input via the create_lM() function and its relatives). Defaults to
  //' \code{FALSE}.
  //' 
  //' @return This function returns a list with a modified \code{stageframe}
  //' usable in MPM construction, an associated \code{repmatrix}, and a general
  //' \code{supplement} table that takes over for any input \code{supplement} or
  //' \code{overwrite} table. Note that if a \code{supplement} is provided and a
  //' \code{repmatrix} is not, or if \code{repmatrix} is set to 0, then it will be
  //' assumed that a \code{repmatrix} should not be used.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List sf_reassess_internal(const DataFrame& stageframe,
    Nullable<DataFrame> supplement = R_NilValue,
    Nullable<DataFrame> overwrite = R_NilValue,
    Nullable<NumericMatrix> repmatrix = R_NilValue,
    bool agemat = false, bool historical = false, int format = 1,
    bool import_lM = false) {
    
    bool supp_provided = false;
    
    Rcpp::DataFrame supplement_true;
    arma::mat repmatrix_true;
    
    StringVector stagevec = as<StringVector>(stageframe["stage"]);
    
    NumericVector origsizevec;
    NumericVector origsizebvec;
    NumericVector origsizecvec;
    if (stageframe.containsElementNamed("size")) {
      origsizevec = as<NumericVector>(stageframe["size"]);
      origsizebvec = as<NumericVector>(stageframe["size_b"]);
      origsizecvec = as<NumericVector>(stageframe["size_c"]);
    } else {
      origsizevec = as<NumericVector>(stageframe["original_size"]);
      origsizebvec = as<NumericVector>(stageframe["original_size_b"]);
      origsizecvec = as<NumericVector>(stageframe["original_size_c"]);
    }
    
    NumericVector minagevec = as<NumericVector>(stageframe["min_age"]);
    NumericVector maxagevec = as<NumericVector>(stageframe["max_age"]);
    arma::uvec repvec = as<arma::uvec>(stageframe["repstatus"]);
    arma::uvec obsvec = as<arma::uvec>(stageframe["obsstatus"]);
    arma::uvec propvec = as<arma::uvec>(stageframe["propstatus"]);
    arma::uvec immvec = as<arma::uvec>(stageframe["immstatus"]);
    arma::uvec matvec = as<arma::uvec>(stageframe["matstatus"]);
    arma::uvec indvec = as<arma::uvec>(stageframe["indataset"]);
    NumericVector binvec = as<NumericVector>(stageframe["binhalfwidth_raw"]);
    NumericVector binbvec = as<NumericVector>(stageframe["binhalfwidthb_raw"]);
    NumericVector bincvec = as<NumericVector>(stageframe["binhalfwidthc_raw"]);
    NumericVector sizeminvec = as<NumericVector>(stageframe["sizebin_min"]);
    NumericVector sizemaxvec = as<NumericVector>(stageframe["sizebin_max"]);
    NumericVector sizectrvec = as<NumericVector>(stageframe["sizebin_center"]);
    NumericVector sizewidthvec = as<NumericVector>(stageframe["sizebin_width"]);
    NumericVector sizebminvec = as<NumericVector>(stageframe["sizebinb_min"]);
    NumericVector sizebmaxvec = as<NumericVector>(stageframe["sizebinb_max"]);
    NumericVector sizebctrvec = as<NumericVector>(stageframe["sizebinb_center"]);
    NumericVector sizebwidthvec = as<NumericVector>(stageframe["sizebinb_width"]);
    NumericVector sizecminvec = as<NumericVector>(stageframe["sizebinc_min"]);
    NumericVector sizecmaxvec = as<NumericVector>(stageframe["sizebinc_max"]);
    NumericVector sizecctrvec = as<NumericVector>(stageframe["sizebinc_center"]);
    NumericVector sizecwidthvec = as<NumericVector>(stageframe["sizebinc_width"]);
    arma::ivec groupvec = as<arma::ivec>(stageframe["group"]);
    StringVector comvec = as<StringVector>(stageframe["comments"]);
    
    StringVector stage3_supp;
    StringVector stage2_supp;
    StringVector stage1_supp;
    NumericVector age2_supp;
    StringVector eststage3_supp;
    StringVector eststage2_supp;
    StringVector eststage1_supp;
    NumericVector estage2_supp;
    NumericVector givenrate_supp;
    NumericVector offset_supp;
    NumericVector multiplier_supp;
    IntegerVector convtype_supp;
    IntegerVector convtype_t12_supp;
    
    int stageframe_length {static_cast<int>(repvec.n_elem)};
    arma::uvec repentryvec(stageframe_length, fill::zeros);
    IntegerVector stage_id = seq(1, stageframe_length);
    
    // Identify all groups in stageframe
    arma::ivec all_groups = unique(groupvec);
    int no_groups {static_cast<int>(all_groups.n_elem)};
    StringVector group_text(no_groups);
    
    for (int i = 0; i < no_groups; i++) {
      group_text(i) = "group";
      group_text(i) += std::to_string(all_groups(i));
    }
    
    if (supplement.isNotNull()) {
      supp_provided = true;
      
      Rcpp::DataFrame supplement_thru(supplement);
      supplement_true = supplement_thru;
      
      stage3_supp = supplement_true["stage3"];
      stage2_supp = supplement_true["stage2"];
      multiplier_supp = supplement_true["multiplier"];
      offset_supp = supplement_true["offset"];
      convtype_supp = supplement_true["convtype"];
    }
    
    if (repmatrix.isNotNull()) {
      NumericMatrix repmatrix_thru(repmatrix);
      repmatrix_true = as<arma::mat>(repmatrix_thru);
      
      arma::rowvec rep_sums = sum(repmatrix_true, 0);
      arma::vec rep_rowsums = sum(repmatrix_true, 1);
      double repmat_sum {static_cast<double>(sum(rep_sums))};
      
      if (static_cast<int>(rep_sums.n_elem) != stageframe_length) {
        throw Rcpp::exception("Argument repmatrix must be a square matrix.");
      }
      
      if (repmat_sum > 0) {
        for (int i = 0; i < static_cast<int>(rep_sums.n_elem); i++) {
          if (rep_sums(i) > 0) {
            repvec(i) = 1;
          } else {
            repvec(i) = 0;
          }
          
          if (rep_rowsums(i) > 0) {
            repentryvec(i) = 1;
          } else {
            repentryvec(i) = 0;
          }
        }
      } else {
        repentryvec = propvec + immvec;
        int rev_checksum {static_cast<int>(sum(repentryvec))};
        
        if (rev_checksum == 0) {
          repentryvec(0) = 1;
          Rf_warningcall(R_NilValue,
            "Reproductive entry stages unknown. Assuming stage 1 is entry stage.");
        }
        
        arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
          
        arma::uvec repentry_calls = find(repentryvec);
        arma::uvec rep_calls = find(repvec);
        
        if (repentry_calls.n_elem > 0 && rep_calls.n_elem > 0) {
          for (int i = 0; i < static_cast<int>(repentry_calls.n_elem); i++) {
            for (int j = 0; j < static_cast<int>(rep_calls.n_elem); j++) {
              token_mat(repentry_calls(i), rep_calls(j)) = 1;
            }
          }
        }
        repmatrix_true = token_mat;
      }
    } else if (supp_provided) {
      arma::ivec cv_supp_arma = as<arma::ivec>(convtype_supp);
      arma::uvec mult_elems = find(cv_supp_arma == 3);
      
      int fec_mults {static_cast<int>(mult_elems.n_elem)};
      for (int i = 0; i < fec_mults; i++) {
        for (int j = 0; j < stageframe_length; j++) {
          if (stage3_supp(mult_elems(i)) == stagevec(j)) {
            repentryvec(j) = 1;
          }
        }
        
        if (stage3_supp(mult_elems(i)) == "prop") {
          arma::uvec propvec_1elems = find(propvec);
          
          int propvec_no1s = static_cast<int>(propvec_1elems.n_elem);
          for (int j = 0; j < propvec_no1s; j++) {
            repentryvec(propvec_1elems(j)) = 1;
          }
        }
        
        if (stage3_supp(mult_elems(i)) == "immat") {
          arma::uvec immvec_1elems = find(immvec);
          
          int immvec_no1s = static_cast<int>(immvec_1elems.n_elem);
          for (int j = 0; j < immvec_no1s; j++) {
            repentryvec(immvec_1elems(j)) = 1;
          }
        }
        
        for (int j = 0; j < no_groups; j++) {
          if (stage3_supp(mult_elems(i)) == group_text(j)) {
            arma::uvec group_elems = find(groupvec == all_groups(j));
            
            int group_noelems = static_cast<int>(group_elems.n_elem);
            for (int k = 0; k < group_noelems; k++) {
              repentryvec(group_elems(k)) = 1;
            }
          }
        }
      }
      
      if (agemat) {
        arma::ivec needed_reprods(stageframe_length, fill::zeros);
        arma::vec needed_mults(stageframe_length, fill::zeros);
        
        for (int i = 0; i < fec_mults; i++) {
          for (int j = 0; j < stageframe_length; j++) {
            if (stage2_supp(mult_elems(i)) == stagevec(j)) {
              needed_reprods(j) = 1;
              needed_mults(j) = 1.0; // multiplier_supp(mult_elems(i));
            }
            
            if (stage2_supp(mult_elems(i)) == "rep") {
              arma::uvec repvec_1elems = find(repvec);
              
              int repvec_no1s = static_cast<int>(repvec_1elems.n_elem);
              for (int k = 0; k < repvec_no1s; k++) {
                needed_reprods(repvec_1elems(k)) = 1;
                needed_mults(repvec_1elems(k)) = 1.0; // multiplier_supp(mult_elems(i));
              }
            }
            
            for (int k = 0; k < no_groups; k++) {
              if (stage2_supp(mult_elems(i)) == group_text(k)) {
                arma::uvec group_elems = find(groupvec == all_groups(k));
                
                int group_noelems = static_cast<int>(group_elems.n_elem);
                for (int l = 0; l < group_noelems; l++) {
                  needed_reprods(group_elems(l)) = 1;
                  needed_mults(group_elems(l)) = 1.0; // multiplier_supp(mult_elems(i));
                }
              }
            }
          }
        }
        
        arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
        
        arma::uvec repentry_calls = find(repentryvec);
        arma::uvec neededrep_calls = find(needed_reprods);
        
        if (repentry_calls.n_elem > 0 && neededrep_calls.n_elem > 0) {
          for (int i = 0; i < static_cast<int>(repentry_calls.n_elem); i++) {
            for (int j = 0; j < static_cast<int>(neededrep_calls.n_elem); j++) {
              token_mat(repentry_calls(i), neededrep_calls(j)) = needed_mults(neededrep_calls(j));
            }
          }
        }
        repmatrix_true = token_mat;
        
      } else {
        arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
        repmatrix_true = token_mat;
      }
    } else {
      int count_of_offspring_stages {0};
      
      for (int i = 0; i < stageframe_length; i++) {
        if (propvec(i) > 0) {
          repentryvec(i) = 1;
          count_of_offspring_stages++;
        } else if (immvec(i) > 0) {
          repentryvec(i) = 1;
          count_of_offspring_stages++;
        }
      }
      
      int rev_checksum = sum(repentryvec);
      if(rev_checksum == 0) {
        repentryvec(0) = 1;
        Rf_warningcall(R_NilValue,
          "Reproductive entry stages unknown. Assuming stage 1 is entry stage.");
      } else if (!import_lM) {
        Rf_warningcall(R_NilValue,
          "No supplement provided. Assuming fecundity yields all propagule and immature stages.");
      }
      
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
        
      arma::uvec repentry_calls = find(repentryvec);
      arma::uvec rep_calls = find(repvec);
      
      if (repentry_calls.n_elem > 0 && rep_calls.n_elem > 0) {
        for (int i = 0; i < static_cast<int>(repentry_calls.n_elem); i++) {
          for (int j = 0; j < static_cast<int>(rep_calls.n_elem); j++) {
            token_mat(repentry_calls(i), rep_calls(j)) = 1;
          }
        }
      }
      repmatrix_true = token_mat;
      
      // Build new supplement for this case
      StringVector novel_stage3 (count_of_offspring_stages);
      StringVector novel_stage2 (count_of_offspring_stages);
      StringVector novel_stage1 (count_of_offspring_stages);
      IntegerVector novel_convtype_t12 (count_of_offspring_stages);
      NumericVector novel_age2 (count_of_offspring_stages);
      StringVector novel_eststage3 (count_of_offspring_stages);
      StringVector novel_eststage2 (count_of_offspring_stages);
      StringVector novel_eststage1 (count_of_offspring_stages);
      NumericVector novel_estage2 (count_of_offspring_stages);
      NumericVector novel_givenrate (count_of_offspring_stages);
      NumericVector novel_offset (count_of_offspring_stages);
      NumericVector novel_multiplier (count_of_offspring_stages);
      IntegerVector novel_convtype (count_of_offspring_stages);
      
      int place_counter {0};
      
      for (int i = 0; i < stageframe_length; i++) { 
        if (repentryvec(i) == 1) {
          novel_stage3(place_counter) = String(stagevec(i));
          novel_stage2(place_counter) = "rep";
          
          if (historical) {
            novel_stage1(place_counter) = "all";
            novel_convtype_t12(place_counter) = 1;
          } else {
            novel_stage1(place_counter) = NA_STRING;
            novel_convtype_t12(place_counter) = NA_INTEGER;
          }
          
          novel_age2(place_counter) = NA_REAL;
        
          novel_eststage3(place_counter) = NA_STRING;
          novel_eststage2(place_counter) = NA_STRING;
          novel_eststage1(place_counter) = NA_STRING;
          novel_estage2(place_counter) = NA_REAL;
          
          novel_givenrate(place_counter) = NA_REAL;
          novel_offset(place_counter) = 0.;
          novel_multiplier(place_counter) = 1.;
          novel_convtype(place_counter) = 3;
          
          place_counter++;
        }
      }
      
      Rcpp::DataFrame supplement_thru = DataFrame::create(_["stage3"] = novel_stage3,
        _["stage2"] = novel_stage2, _["stage1"] = novel_stage1, _["age2"] = novel_age2,
        _["eststage3"] =  novel_eststage3, _["eststage2"] =  novel_eststage2,
        _["eststage1"] = novel_eststage1, _["estage2"] = novel_estage2,
        _["givenrate"] = novel_givenrate, _["offset"] = novel_offset,
        _["multiplier"] = novel_multiplier, _["convtype"] = novel_convtype,
        _["convtype_t12"] = novel_convtype_t12);
      
      supplement_true = supplement_thru;
      
      stage3_supp = novel_stage3;
      stage2_supp = novel_stage2;
      stage1_supp = novel_stage1;
      age2_supp = novel_age2;
      eststage3_supp = novel_eststage3;
      eststage2_supp = novel_eststage2;
      eststage1_supp = novel_eststage1;
      estage2_supp = novel_estage2;
      givenrate_supp = novel_givenrate;
      offset_supp = novel_offset;
      multiplier_supp = novel_multiplier;
      convtype_supp = novel_convtype;
      convtype_t12_supp = novel_convtype_t12;
      
      supp_provided = true;
    }
    
    // Reorder stageframe
    arma::uvec prop_stages = find(propvec);
    arma::uvec prop0_stages = find(propvec == 0);
    arma::uvec imm_stages = find(immvec);
    arma::uvec imm0_stages = find(immvec == 0);
    arma::uvec p0_im_stages = intersect(prop0_stages, imm_stages);
    
    arma::uvec mat_stages = find(matvec);
    arma::uvec mat_imm0_stages = intersect(mat_stages, imm0_stages);
    arma::uvec rep_stages = find(repvec);
    arma::uvec rep0_stages = find(repvec == 0);
    arma::uvec rep0_mat_imm0_stages = intersect(mat_imm0_stages, rep0_stages);
    arma::uvec mat_rep0_stages = intersect(mat_stages, rep0_stages);
    arma::uvec neworder(stageframe_length, fill::zeros);
    arma::uvec tracked_elems(stageframe_length, fill::ones);
    
    int no_prop_stages {static_cast<int>(prop_stages.n_elem)};
    int no_p0_im_stages {static_cast<int>(p0_im_stages.n_elem)};
    int no_r0_im_mt_stages {static_cast<int>(rep0_mat_imm0_stages.n_elem)};
    
    int counter {0};
    for (int j = 0; j < no_prop_stages; j++) {
      neworder(counter) = prop_stages(j);
      tracked_elems(prop_stages(j)) = 0;
      counter++;
    }
    for (int j = 0; j < no_p0_im_stages; j++) {
      neworder(counter) = p0_im_stages(j);
      tracked_elems(p0_im_stages(j)) = 0;
      counter++;
    }
    for (int j = 0; j < no_r0_im_mt_stages; j++) {
      neworder(counter) = rep0_mat_imm0_stages(j);
      tracked_elems(rep0_mat_imm0_stages(j)) = 0;
      counter++;
    }
    
    arma::uvec remaining_elems = find(tracked_elems);
    int no_remaining_elems = static_cast<int>(remaining_elems.n_elem);
    for (int j = 0; j < no_remaining_elems; j++) {
      neworder(counter) = remaining_elems(j);
      counter++;
    }
    
    StringVector newstagevec(stageframe_length + format);
    IntegerVector newstageidvec(stageframe_length + format);
    NumericVector neworigsizevec(stageframe_length + format);
    NumericVector neworigsizebvec(stageframe_length + format);
    NumericVector neworigsizecvec(stageframe_length + format);
    NumericVector newminagevec(stageframe_length + format);
    NumericVector newmaxagevec(stageframe_length + format);
    arma::uvec newrepvec((stageframe_length + format), fill::zeros);
    arma::uvec newobsvec((stageframe_length + format), fill::zeros);
    arma::uvec newpropvec((stageframe_length + format), fill::zeros);
    arma::uvec newimmvec((stageframe_length + format), fill::zeros);
    arma::uvec newmatvec((stageframe_length + format), fill::zeros);
    arma::uvec newrepentryvec((stageframe_length + format), fill::zeros);
    arma::uvec newindvec((stageframe_length + format), fill::zeros);
    NumericVector newbinvec((stageframe_length + format));
    NumericVector newbinbvec((stageframe_length + format));
    NumericVector newbincvec((stageframe_length + format));
    NumericVector newsizeminvec((stageframe_length + format));
    NumericVector newsizemaxvec((stageframe_length + format));
    NumericVector newsizectrvec((stageframe_length + format));
    NumericVector newsizewidthvec((stageframe_length + format));
    NumericVector newsizebminvec((stageframe_length + format));
    NumericVector newsizebmaxvec((stageframe_length + format));
    NumericVector newsizebctrvec((stageframe_length + format));
    NumericVector newsizebwidthvec((stageframe_length + format));
    NumericVector newsizecminvec((stageframe_length + format));
    NumericVector newsizecmaxvec((stageframe_length + format));
    NumericVector newsizecctrvec((stageframe_length + format));
    NumericVector newsizecwidthvec((stageframe_length + format));
    arma::ivec newgroupvec((stageframe_length + format), fill::zeros);
    StringVector newcomvec((stageframe_length + format));
    arma::uvec newalive((stageframe_length + format), fill::zeros);
    NumericVector newalmostbornvec((stageframe_length + format));
    
    arma::mat repmat1(stageframe_length, stageframe_length);
    arma::mat repmat2(stageframe_length, stageframe_length);
    
    for (int i = 0; i < stageframe_length; i++) {
      newstagevec(i) = stagevec(neworder(i));
      newstageidvec(i) = stage_id(neworder(i));
      neworigsizevec(i) = origsizevec(neworder(i));
      neworigsizebvec(i) = origsizebvec(neworder(i));
      neworigsizecvec(i) = origsizecvec(neworder(i));
      newminagevec(i) = minagevec(neworder(i));
      if (NumericVector::is_na(newminagevec(i))) newminagevec(i) = 0.0;
      newmaxagevec(i) = maxagevec(neworder(i));
      newrepvec(i) = repvec(neworder(i));
      newobsvec(i) = obsvec(neworder(i));
      newpropvec(i) = propvec(neworder(i));
      newimmvec(i) = immvec(neworder(i));
      newmatvec(i) = matvec(neworder(i));
      newrepentryvec(i) = repentryvec(neworder(i));
      newindvec(i) = indvec(neworder(i));
      newbinvec(i) = binvec(neworder(i));
      newbinbvec(i) = binbvec(neworder(i));
      newbincvec(i) = bincvec(neworder(i));
      newsizeminvec(i) = sizeminvec(neworder(i));
      newsizemaxvec(i) = sizemaxvec(neworder(i));
      newsizectrvec(i) = sizectrvec(neworder(i));
      newsizewidthvec(i) = sizewidthvec(neworder(i));
      newsizebminvec(i) = sizebminvec(neworder(i));
      newsizebmaxvec(i) = sizebmaxvec(neworder(i));
      newsizebctrvec(i) = sizebctrvec(neworder(i));
      newsizebwidthvec(i) = sizebwidthvec(neworder(i));
      newsizecminvec(i) = sizecminvec(neworder(i));
      newsizecmaxvec(i) = sizecmaxvec(neworder(i));
      newsizecctrvec(i) = sizecctrvec(neworder(i));
      newsizecwidthvec(i) = sizecwidthvec(neworder(i));
      newgroupvec(i) = groupvec(neworder(i));
      newcomvec(i) = comvec(neworder(i));
      newalive(i) = 1;
      
      repmat1.col(i) = repmatrix_true.col(neworder(i));
    }
    
    for (int i = 0; i < stageframe_length; i++) {
      repmat2.row(i) = repmat1.row(neworder(i));
    }
    
    if (format == 2) {
      newstagevec(stageframe_length) = "AlmostBorn";
      newstageidvec(stageframe_length) = stageframe_length + 1;
      neworigsizevec(stageframe_length) = 0.0;
      neworigsizebvec(stageframe_length) = 0.0;
      neworigsizecvec(stageframe_length) = 0.0;
      newrepvec(stageframe_length) = 0;
      newobsvec(stageframe_length) = 1;
      newpropvec(stageframe_length) = 0;
      newimmvec(stageframe_length) = 1;
      newmatvec(stageframe_length) = 0;
      newrepentryvec(stageframe_length) = 0;
      newindvec(stageframe_length) = 1;
      newbinvec(stageframe_length) = 0.0;
      newbinbvec(stageframe_length) = 0.0;
      newbincvec(stageframe_length) = 0.0;
      newsizeminvec(stageframe_length) = 0.0;
      newsizemaxvec(stageframe_length) = 0.0;
      newsizectrvec(stageframe_length) = 0.0;
      newsizewidthvec(stageframe_length) = 0.0;
      newsizebminvec(stageframe_length) = 0.0;
      newsizebmaxvec(stageframe_length) = 0.0;
      newsizebctrvec(stageframe_length) = 0.0;
      newsizebwidthvec(stageframe_length) = 0.0;
      newsizecminvec(stageframe_length) = 0.0;
      newsizecmaxvec(stageframe_length) = 0.0;
      newsizecctrvec(stageframe_length) = 0.0;
      newsizecwidthvec(stageframe_length) = 0.0;
      newgroupvec(stageframe_length) = no_groups + 1;
      newcomvec(stageframe_length) = "Almost born (t-1)";
      newalive(stageframe_length) = 0;
      newalmostbornvec(stageframe_length) = 1.0;
      
      if (agemat) {
        newminagevec(stageframe_length) = 0.0;
        newmaxagevec(stageframe_length) = 0.0;
      } else {
        newminagevec(stageframe_length) = NA_REAL;
        newmaxagevec(stageframe_length) = NA_REAL;
      }
    }
    
    newstagevec(stageframe_length + (format - 1)) = "Dead";
    newstageidvec(stageframe_length + (format - 1)) = stageframe_length + 1 + (format - 1);
    neworigsizevec(stageframe_length + (format - 1)) = 0.0;
    neworigsizebvec(stageframe_length + (format - 1)) = 0.0;
    neworigsizecvec(stageframe_length + (format - 1)) = 0.0;
    newrepvec(stageframe_length + (format - 1)) = 0;
    newobsvec(stageframe_length + (format - 1)) = 1;
    newpropvec(stageframe_length + (format - 1)) = 0;
    newimmvec(stageframe_length + (format - 1)) = 0;
    newmatvec(stageframe_length + (format - 1)) = 1;
    newrepentryvec(stageframe_length + (format - 1)) = 0;
    newindvec(stageframe_length + (format - 1)) = 1;
    newbinvec(stageframe_length + (format - 1)) = 0.0;
    newbinbvec(stageframe_length + (format - 1)) = 0.0;
    newbincvec(stageframe_length + (format - 1)) = 0.0;
    newsizeminvec(stageframe_length + (format - 1)) = 0.0;
    newsizemaxvec(stageframe_length + (format - 1)) = 0.0;
    newsizectrvec(stageframe_length + (format - 1)) = 0.0;
    newsizewidthvec(stageframe_length + (format - 1)) = 0.0;
    newsizebminvec(stageframe_length + (format - 1)) = 0.0;
    newsizebmaxvec(stageframe_length + (format - 1)) = 0.0;
    newsizebctrvec(stageframe_length + (format - 1)) = 0.0;
    newsizebwidthvec(stageframe_length + (format - 1)) = 0.0;
    newsizecminvec(stageframe_length + (format - 1)) = 0.0;
    newsizecmaxvec(stageframe_length + (format - 1)) = 0.0;
    newsizecctrvec(stageframe_length + (format - 1)) = 0.0;
    newsizecwidthvec(stageframe_length + (format - 1)) = 0.0;
    newgroupvec(stageframe_length + (format - 1)) = no_groups + 1;
    newcomvec(stageframe_length + (format - 1)) = "Dead";
    newalive(stageframe_length + (format - 1)) = 0;
    
    if (agemat) {
      newminagevec(stageframe_length + (format - 1)) = 0.0;
      newmaxagevec(stageframe_length + (format - 1)) = 0.0;
    } else {
      newminagevec(stageframe_length + (format - 1)) = NA_REAL;
      newmaxagevec(stageframe_length + (format - 1)) = NA_REAL;
    }
    
    StringVector newstagevec_check = unique(newstagevec);
    if (newstagevec_check.length() < newstagevec.length()) {
      throw Rcpp::exception("All stage names must be unique.", false);
    }
    
    Rcpp::List newstageframe(33);
    
    newstageframe(0) = newstageidvec;
    newstageframe(1) = newstagevec;
    newstageframe(2) = neworigsizevec;
    newstageframe(3) = neworigsizebvec;
    newstageframe(4) = neworigsizecvec;
    newstageframe(5) = newminagevec;
    newstageframe(6) = newmaxagevec;
    newstageframe(7) = Rcpp::IntegerVector(newrepvec.begin(), newrepvec.end());
    newstageframe(8) = Rcpp::IntegerVector(newobsvec.begin(), newobsvec.end());
    newstageframe(9) = Rcpp::IntegerVector(newpropvec.begin(), newpropvec.end());
    newstageframe(10) = Rcpp::IntegerVector(newimmvec.begin(), newimmvec.end());
    newstageframe(11) = Rcpp::IntegerVector(newmatvec.begin(), newmatvec.end());
    newstageframe(12) = Rcpp::IntegerVector(newrepentryvec.begin(), newrepentryvec.end());
    newstageframe(13) = Rcpp::IntegerVector(newindvec.begin(), newindvec.end());
    newstageframe(14) = newbinvec;
    newstageframe(15) = newsizeminvec;
    newstageframe(16) = newsizemaxvec;
    newstageframe(17) = newsizectrvec;
    newstageframe(18) = newsizewidthvec;
    newstageframe(19) = newbinbvec;
    newstageframe(20) = newsizebminvec;
    newstageframe(21) = newsizebmaxvec;
    newstageframe(22) = newsizebctrvec;
    newstageframe(23) = newsizebwidthvec;
    newstageframe(24) = newbincvec;
    newstageframe(25) = newsizecminvec;
    newstageframe(26) = newsizecmaxvec;
    newstageframe(27) = newsizecctrvec;
    newstageframe(28) = newsizecwidthvec;
    newstageframe(29) = Rcpp::IntegerVector(newgroupvec.begin(), newgroupvec.end());
    newstageframe(30) = newcomvec;
    newstageframe(31) = Rcpp::IntegerVector(newalive.begin(), newalive.end());
    newstageframe(32) = newalmostbornvec; // For use in deVries-format hMPMs
    
    CharacterVector namevec = {"stage_id", "stage", "original_size",
      "original_size_b", "original_size_c", "min_age", "max_age", "repstatus",
      "obsstatus", "propstatus", "immstatus", "matstatus", "entrystage",
      "indataset", "binhalfwidth_raw", "sizebin_min", "sizebin_max",
      "sizebin_center", "sizebin_width", "binhalfwidthb_raw", "sizebinb_min",
      "sizebinb_max", "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw",
      "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
      "group", "comments", "alive", "almostborn"};
    CharacterVector newclasses = {"data.frame", "stageframe"};
    newstageframe.attr("names") = namevec;
    newstageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
        (stageframe_length + format));
    newstageframe.attr("class") = newclasses;
    
    if (supp_provided) {
      DataFrame newsupplement = supp_reassess(newstageframe, historical,
        supplement_true, overwrite);
      
      return Rcpp::List::create(_["stageframe"] = newstageframe,
        _["repmatrix"] = repmat2, _["ovtable"] = newsupplement);
    } else {
      return Rcpp::List::create(_["stageframe"] = newstageframe,
        _["repmatrix"] = repmat2, _["ovtable"] = NULL);
    }
  }
  
  //' Create Stageframe for Population Matrix Projection Analysis
  //' 
  //' Function \code{sf_leslie()} returns a data frame describing each age in a
  //' Leslie MPM in terms of ahistorical stage information. This function is
  //' internal to \code{rleslie()} and \code{fleslie()}.
  //' 
  //' @name sf_leslie
  //' 
  //' @param min_age The first age to include in the matrix.
  //' @param max_age The maximum age to include in the matrix.
  //' @param min_fecage The first age in which reproduction is possible.
  //' @param max_fecage The final age in which reproduction is possible.
  //' @param cont A logical value indicating whether survival continues past the
  //' last described age.
  //' 
  //' @return A data frame of class \code{stageframe}, which includes information
  //' on the stage name, size, reproductive status, observation status, propagule 
  //' status, immaturity status, maturity status, presence within the core dataset, 
  //' stage group classification, raw bin half-width, and the minimum, 
  //' center, and maximum of each size bin, as well as its width. If minimum and
  //' maximum ages were specified, then these are also included. Also includes an 
  //' empty string variable that can be used to describe stages meaningfully.
  //' 
  //' Variables in this data frame include the following:
  //' \item{stage_id}{An unique integer representing each age, in order.}
  //' \item{stage}{The unique names of the ages to be analyzed.}
  //' \item{original_size}{The typical or representative size at which each stage
  //' occurs. Since ages are not characterized by size, this is generally
  //' \code{NA}.}
  //' \item{original_size_b}{Size at which each stage occurs in terms of a second
  //' size variable, if one exists. In Leslie MPMs, generally \code{NA}.}
  //' \item{original_size_c}{Size at which each stage occurs in terms of a third
  //' size variable, if one exists. In Leslie MPMs, generally \code{NA}.}
  //' \item{min_age}{The minimum age at which the stage may occur. In Leslie MPMs,
  //' defaults to the current age.}
  //' \item{max_age}{The maximum age at which the stage may occur. In Leslie MPMs,
  //' will generally equal the current age or \code{NA}, depending on whether
  //' individuals are allowed to remain at the maximum age.}
  //' \item{repstatus}{A binomial variable showing whether each age is
  //' reproductive.}
  //' \item{obsstatus}{A binomial variable showing whether each age is
  //' observable.}
  //' \item{propstatus}{A binomial variable showing whether each age is a
  //' propagule.}
  //' \item{immstatus}{A binomial variable showing whether each age can occur as
  //' immature.}
  //' \item{matstatus}{A binomial variable showing whether each age occurs in
  //' maturity.}
  //' \item{entrystage}{A binomial variable showing whether each age is an entry
  //' stage. In Leslie MPMs, only the first stage is set to \code{1}, while all
  //' others are set to \code{0}.}
  //' \item{indataset}{A binomial variable describing whether each age occurs in
  //' the input dataset.}
  //' \item{binhalfwidth_raw}{The half-width of the size bin, as input.}
  //' \item{sizebin_min}{The minimum primary size at which the age may occur.}
  //' \item{sizebin_max}{The maximum primary size at which the age may occur.}
  //' \item{sizebin_center}{The midpoint of the primary size bin at which the age
  //' may occur.}
  //' \item{sizebin_width}{The width of the primary size bin corresponding to the
  //' age.}
  //' \item{binhalfwidthb_raw}{The half-width of the size bin of a second size
  //' variable, as input.}
  //' \item{sizebinb_min}{The minimum secondary size at which the age may occur.}
  //' \item{sizebinb_max}{The maximum secondary size at which the age may occur.}
  //' \item{sizebinb_center}{The midpoint of the secondary size bin at which the
  //' age may occur.}
  //' \item{sizebinb_width}{The width of the secondary size bin corresponding to
  //' the age.}
  //' \item{binhalfwidthc_raw}{The half-width of the size bin of a third size
  //' variable, as input.}
  //' \item{sizebinc_min}{The minimum tertiary size at which the age may occur.}
  //' \item{sizebinc_max}{The maximum tertiary size at which the age may occur.}
  //' \item{sizebinc_center}{The midpoint of the tertiary size bin at which the
  //' age may occur.}
  //' \item{sizebinc_width}{The width of the tertiary size bin corresponding to
  //' the age.}
  //' \item{group}{An integer denoting the size classification group that the
  //' age falls within.}
  //' \item{comments}{A text field for stage descriptions.}
  //' \item{alive}{An integer vector denoting whether the age is alive. Defaults
  //' to \code{1} for all ages.}
  //' \item{almost_born}{An integer vector denoting whether the age corresponds to
  //' the prior stage of a newly produced individual in a historical model. In
  //' Leslie MPMs, defaults to \code{0}.}
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List sf_leslie (int min_age, int max_age, int min_fecage,
    int max_fecage, bool cont) {
    
    const int age_range {max_age - min_age};
    const int total_ages {age_range + 1};
    if (age_range < 1) {
      throw Rcpp::exception("Leslie MPMs require at least 2 ages.", false);
    }
    
    Rcpp::List output_longlist(33);
    Rcpp::CharacterVector varnames {"stage_id", "stage", "original_size",
      "original_size_b", "original_size_c", "min_age", "max_age", "repstatus",
      "obsstatus", "propstatus", "immstatus", "matstatus", "entrystage",
      "indataset", "binhalfwidth_raw", "sizebin_min", "sizebin_max",
      "sizebin_center", "sizebin_width", "binhalfwidthb_raw", "sizebinb_min",
      "sizebinb_max", "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw",
      "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
      "group", "comments", "alive", "almostborn"};
    
    int matsize {total_ages};
    
    IntegerVector stage_id (matsize, 1);
    StringVector agenames_true (matsize, NA_STRING);
    NumericVector size_true (matsize, NA_REAL);
    NumericVector sizesb_true (matsize, NA_REAL);
    NumericVector sizesc_true (matsize, NA_REAL);
    IntegerVector repstatus_true (matsize, 0);
    IntegerVector obsstatus_true (matsize, 1);
    IntegerVector propstatus_true (matsize, 0);
    IntegerVector matstatus_true (matsize, 0);
    IntegerVector immstatus_true (matsize, 1);
    IntegerVector entrystage_true (matsize, 0);
    IntegerVector indataset_true (matsize, 1);
    NumericVector minage_true (matsize, 0.0);
    NumericVector maxage_true (matsize, 0.0);
    NumericVector binhalfwidth_true (matsize, NA_REAL);
    NumericVector binhalfwidthb_true (matsize, NA_REAL);
    NumericVector binhalfwidthc_true (matsize, NA_REAL);
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
    IntegerVector alive_true (matsize, 1);
    IntegerVector almost_born (matsize, 0);
    
    entrystage_true(0) = 1;
    
    for (int i = 0; i < total_ages; i++) {
      stage_id(i) = i + 1;
      Rcpp::String part1 {"Age"};
      part1 += (static_cast<char>(i + min_age));
      agenames_true(i) = part1;
      
      if ((i + min_age) >= min_fecage) {
        repstatus_true(i) = 1;
        matstatus_true(i) = 1;
        immstatus_true(i) = 0;
      }
      
      minage_true(i) = static_cast<double>(i + min_age);
      if ((i + min_age < max_age) || !cont) {
        maxage_true(i) = static_cast<double>(i + min_age);
      } else {
        maxage_true(i) = NA_REAL;
      }
    }
    
    output_longlist(0) = stage_id;
    output_longlist(1) = agenames_true;
    output_longlist(2) = size_true;
    output_longlist(3) = sizesb_true;
    output_longlist(4) = sizesc_true;
    
    output_longlist(5) = minage_true;
    output_longlist(6) = maxage_true;
    output_longlist(7) = repstatus_true;
    output_longlist(8) = obsstatus_true;
    output_longlist(9) = propstatus_true;
    output_longlist(10) = immstatus_true;
    output_longlist(11) = matstatus_true;
    output_longlist(12) = entrystage_true;
    
    output_longlist(13) = indataset_true;
    
    output_longlist(14) = binhalfwidth_true;
    output_longlist(15) = sizebin_min;
    output_longlist(16) = sizebin_max;
    output_longlist(17) = sizebin_center;
    output_longlist(18) = sizebin_width;
    
    output_longlist(19) = binhalfwidthb_true;
    output_longlist(20) = sizebinb_min;
    output_longlist(21) = sizebinb_max;
    output_longlist(22) = sizebinb_center;
    output_longlist(23) = sizebinb_width;
    
    output_longlist(24) = binhalfwidthc_true;
    output_longlist(25) = sizebinc_min;
    output_longlist(26) = sizebinc_max;
    output_longlist(27) = sizebinc_center;
    output_longlist(28) = sizebinc_width;
    
    output_longlist(29) = group_true;
    output_longlist(30) = comments_true;
    output_longlist(31) = alive_true;
    output_longlist(32) = almost_born;
    
    output_longlist.attr("names") = varnames;
    output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, matsize);
    StringVector needed_classes {"data.frame", "stageframe"};
    output_longlist.attr("class") = needed_classes;
    
    return output_longlist;
  }

}
#endif