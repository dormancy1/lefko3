#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;


// Index of Functions
// 
// 1. DataFrame .bambi3  Creates Size Index for Elasticity Summaries of hMPMs
// 2. DataFrame .bambi2  Creates Size Index for Elasticity Summaries of ahMPMs
// 3. List .demolition4  Sum Positive and Negative LTRE or Elasticity Contributions
// 4. List .demolition3  Creates Summary Data for Elasticity Matrix Inputs
// 5. List .demolition3sp  Creates Summary Data for Elasticity Matrix Inputs
// 6. RObject lambda3  Estimate Actual or Deterministic Population Growth Rate
// 7. DataFrame matrix_interp  Arranges Matrix Elements in Order of Magnitude for Interpretation
// 8. List append_lP  Append Projections Into New lefkoProj Object


//' Creates Size Index for Elasticity Summaries of hMPMs
//' 
//' Function \code{bambi3()} creates an index of estimable elements in
//' historical matrices, and details the kind of transition that it is.
//' 
//' @name .bambi3
//' 
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' 
//' @return A data frame with the following elements:
//' \item{index}{Vector index of matrix element in C++ terms.}
//' \item{transition}{Category of transition.}
//' \item{size3}{Size in occasion \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in occasion \emph{t}+1.}
//' \item{entrystatus3}{Entry status in occasion \emph{t}+1.}
//' \item{size2}{Size in occasion \emph{t}.}
//' \item{repstatus2}{Reproductive status in occasion \emph{t}.}
//' \item{entrystatus2}{Entry status in occasion \emph{t}.}
//' \item{size1}{Size in occasion \emph{t}-1.}
//' \item{repstatus1}{Reproductive status in occasion \emph{t}-1.}
//' \item{entrystatus1}{Entry status in occasion \emph{t}-1.}
//'
//' The kind of transitions conforms to the following code: \code{10}: full
//' stasis, \code{11}: stasis to growth, \code{12}: full growth, \code{13}:
//' growth to stasis, \code{14}: stasis to shrinkage, \code{15}: full shrinkage,
//' \code{16}: shrinkage to stasis, \code{17}: growth to shrinkage, \code{18}:
//' shrinkage to growth, \code{20}: stasis to fecundity, \code{21}: growth to
//' fecundity, \code{22}: shrinkage to fecundity, \code{23}: fecundity to
//' stasis, \code{24}: fecundity to growth, \code{25}: fecundity to shrinkage,
//' \code{26}: fecundity to fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.bambi3)]]
Rcpp::DataFrame bambi3(const DataFrame& stages, const DataFrame& hstages) {
  
  StringVector stagenames = as<StringVector>(stages["stage"]);
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  arma::vec sizes = as<arma::vec>(stages["original_size"]);
  arma::uvec repstatus = as<arma::uvec>(stages["repstatus"]);
  arma::uvec entrystage = as<arma::uvec>(stages["entrystage"]);
  int numstages = static_cast<int>(astages.n_elem);
  
  arma::uvec hstage3in = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec hstage2nin = as<arma::uvec>(hstages["stage_id_1"]);
  int numhstages = static_cast<int>(hstage3in.n_elem);
  
  hstage3in = hstage3in - 1;
  hstage2nin = hstage2nin - 1;
  int predictedsize = numstages * numstages * numstages;
  
  arma::ivec hsindexl(predictedsize);
  arma::uvec transition_type(predictedsize);
  hsindexl.fill(-1);
  transition_type.zeros();
  
  arma::vec size1(predictedsize);
  arma::vec size2(predictedsize);
  arma::vec size3(predictedsize);
  size1.fill(-1);
  size2.fill(-1);
  size3.fill(-1);
  
  arma::uvec repstatus1(predictedsize);
  arma::uvec repstatus2(predictedsize);
  arma::uvec repstatus3(predictedsize);
  repstatus1.zeros();
  repstatus2.zeros();
  repstatus3.zeros();
  
  arma::uvec entrystatus1(predictedsize);
  arma::uvec entrystatus2(predictedsize);
  arma::uvec entrystatus3(predictedsize);
  entrystatus1.zeros();
  entrystatus2.zeros();
  entrystatus3.zeros();
  
  StringVector longnames3(predictedsize);
  StringVector longnames2(predictedsize);
  StringVector longnames1(predictedsize);
  
  int counter = 0;
  
  for (int i1 = 0; i1 < numhstages; i1++) {
    for (int i2 = 0; i2 < numhstages; i2++) {
      if (hstage3in(i1) == (hstage2nin(i2))) {
        hsindexl(counter) = (i1 * numhstages) + i2;
        
        unsigned int stage1 = hstage2nin(i1);
        arma::uvec ahst_rows_1 = find(astages == (stage1 + 1));
        unsigned int ahst_row_1 = ahst_rows_1(0);
        longnames1(counter) = stagenames(ahst_row_1);
        size1(counter) = sizes(ahst_row_1);
        repstatus1(counter) = repstatus(ahst_row_1);
        entrystatus1(counter) = entrystage(ahst_row_1);
        
        unsigned int stage2 = hstage2nin(i2);
        arma::uvec ahst_rows_2 = find(astages == (stage2 + 1));
        unsigned int ahst_row_2 = ahst_rows_2(0);
        longnames2(counter) = stagenames(ahst_row_2);
        size2(counter) = sizes(ahst_row_2);
        repstatus2(counter) = repstatus(ahst_row_2);
        entrystatus2(counter) = entrystage(ahst_row_2);
        
        unsigned int stage3 = hstage3in(i2);
        arma::uvec ahst_rows_3 = find(astages == (stage3 + 1));
        unsigned int ahst_row_3 = ahst_rows_3(0);
        longnames3(counter) = stagenames(ahst_row_3);
        size3(counter) = sizes(ahst_row_3);
        repstatus3(counter) = repstatus(ahst_row_3);
        entrystatus3(counter) = entrystage(ahst_row_3);
        
        if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
          if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
            transition_type(counter) = 26; // Fecundity to fecundity
          } else if (size2(counter) == size1(counter)) {
            if (repstatus2(counter) > repstatus1(counter) ||
                entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 21; // Growth to fecundity
            } else if (repstatus2(counter) < repstatus1(counter) ||
              entrystatus2(counter) > entrystatus1(counter)) {
              transition_type(counter) = 22; // Shrinkage to fecundity
            } else {
              transition_type(counter) = 20; // Stasis to fecundity
            }
          } else if (size2(counter) > size1(counter)) {
            transition_type(counter) = 21; // Growth to fecundity
          } else if (size2(counter) < size1(counter)) {
            transition_type(counter) = 22; // Shrinkage to fecundity
          }
        } else if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
          if (size3(counter) == size2(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 24; // Fecundity to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 25; // Fecundity to shrinkage
            } else {
              transition_type(counter) = 23; // Fecundity to stasis
            }
          } else if (size3(counter) > size2(counter)) {
            transition_type(counter) = 24; // Fecundity to growth
          } else if (size3(counter) < size2(counter)) {
            transition_type(counter) = 25; // Fecundity to shrinkage
          }
        } else if (size3(counter) == size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 12; // Full growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 17; // Growth to shrinkage
            } else {
              transition_type(counter) = 13; // Growth to stasis
            }
          } else if (repstatus2(counter) < repstatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 18; // Shrinkage to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 15; // Full shrinkage
            } else {
              transition_type(counter) = 16; // Shrinkage to stasis
            }
          } else {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 11; // Stasis to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 14; // Stasis to shrinkage
            } else {
              transition_type(counter) = 10; // Full stasis
            }
          }
        } else if (size3(counter) > size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus2(counter) < repstatus1(counter) ||
            entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else {
            transition_type(counter) = 11; // Stasis to growth
          }
        } else if (size3(counter) > size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 12; // Full growth
        } else if (size3(counter) == size2(counter) && size2(counter) > size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) ||
            entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus3(counter) < repstatus2(counter) ||
            entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else {
            transition_type(counter) = 13; // Growth to stasis
          }
        } else if (size3(counter) < size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else if (repstatus2(counter) < repstatus1(counter) ||
            entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 14; // Stasis to shrinkage
          }
        } else if (size3(counter) < size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 15; // Full shrinkage
        } else if (size3(counter) == size2(counter) && size2(counter) < size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) ||
            entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else if (repstatus3(counter) < repstatus2(counter) ||
            entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 16; // Shrinkage to stasis
          }
        } else if (size3(counter) < size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 17; // Growth to shrinkage
        } else if (size3(counter) > size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 18; // Shrinkage to growth
        }
        counter++;
      }
    }
  }
  
  StringVector names3(counter);
  StringVector names2(counter);
  StringVector names1(counter);
  for (int i = 0; i < counter; i++) {
    names3(i) = longnames3(i);
    names2(i) = longnames2(i);
    names1(i) = longnames1(i);
  }
  
  arma::uvec targetindices = find(hsindexl > -1);
  arma::ivec hsindex = hsindexl.elem(targetindices);
  arma::uvec t_type = transition_type.elem(targetindices);
  arma::vec size3c = size3.elem(targetindices);
  arma::vec size2c = size2.elem(targetindices);
  arma::vec size1c = size1.elem(targetindices);
  arma::uvec r_status3 = repstatus3.elem(targetindices);
  arma::uvec r_status2 = repstatus2.elem(targetindices);
  arma::uvec r_status1 = repstatus1.elem(targetindices);
  arma::uvec e_status3 = entrystatus3.elem(targetindices);
  arma::uvec e_status2 = entrystatus2.elem(targetindices);
  arma::uvec e_status1 = entrystatus1.elem(targetindices);
  
  DataFrame output = DataFrame::create(Named("index") = hsindex,
    _["transition"] = t_type, _["stage3"] = names3, _["size3"] = size3c,
    _["repstatus3"] = r_status3, _["entrystatus3"] = e_status3,
    _["stage2"] = names2, _["size2"] = size2c, _["repstatus2"] = r_status2,
    _["entrystatus2"] = e_status2, _["stage1"] = names1, _["size1"] = size1c,
    _["repstatus1"] = r_status1, _["entrystatus1"] = e_status1);
  
  return output;
}

//' Creates Size Index for Elasticity Summaries of ahMPMs
//' 
//' Function \code{bambi2()} creates an index of estimable elements in
//' ahistorical matrices, and details the kind of transition that it is.
//' 
//' @name .bambi2
//' 
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' 
//' @return A data frame with the following elements:
//' \item{index}{Vector index of matrix element in C++ terms.}
//' \item{transition}{Category of transition.}
//' \item{stage3}{Stage in occasion \emph{t}+1.}
//' \item{size3}{Size in occasion \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in occasion \emph{t}+1.}
//' \item{entrystatus3}{Entry status in occasion \emph{t}+1.}
//' \item{stage2}{Stage in occasion \emph{t}.}
//' \item{size2}{Size in occasion \emph{t}.}
//' \item{repstatus2}{Reproductive status in occasion \emph{t}.}
//' \item{entrystatus2}{Entry status in occasion \emph{t}.}
//'
//' The kind of transitions conforms to the following code: \code{1}: stasis, 
//' \code{2}: growth, \code{3}: shrinkage, \code{4}: fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.bambi2)]]
Rcpp::DataFrame bambi2(const DataFrame& stages) {
  
  StringVector stagenames = as<StringVector>(stages["stage"]);
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  arma::vec sizes = as<arma::vec>(stages["original_size"]);
  arma::uvec repstatus = as<arma::uvec>(stages["repstatus"]);
  arma::uvec entrystage = as<arma::uvec>(stages["entrystage"]);
  int numstages = static_cast<int>(astages.n_elem);
  astages = astages - 1;
  int predictedsize = numstages * numstages;
  
  arma::ivec ahsindexl(predictedsize);
  arma::uvec transition_type(predictedsize);
  ahsindexl.fill(-1);
  transition_type.zeros();
  
  StringVector longstages3(predictedsize);
  StringVector longstages2(predictedsize);
  
  arma::vec size2(predictedsize);
  arma::vec size3(predictedsize);
  size2.fill(-1);
  size3.fill(-1);
  
  arma::uvec repstatus2(predictedsize);
  arma::uvec repstatus3(predictedsize);
  repstatus2.zeros();
  repstatus3.zeros();
  
  arma::uvec entrystatus2(predictedsize);
  arma::uvec entrystatus3(predictedsize);
  entrystatus2.zeros();
  entrystatus3.zeros();
  
  int counter = 0;
  
  for (int i1 = 0; i1 < numstages; i1++) {
    for (int i2 = 0; i2 < numstages; i2++) {
      ahsindexl(counter) = (i1 * numstages) + i2;
      
      int stage2 = astages(i1);
      arma::uvec ahst_rows_2 = find(astages == stage2);
      int ahst_row_2 = ahst_rows_2(0);
      longstages2(counter) = stagenames(ahst_row_2);
      size2(counter) = sizes(ahst_row_2);
      repstatus2(counter) = repstatus(ahst_row_2);
      entrystatus2(counter) = entrystage(ahst_row_2);
      
      int stage3 = astages(i2);
      arma::uvec ahst_rows_3 = find(astages == stage3);
      int ahst_row_3 = ahst_rows_3(0);
      longstages3(counter) = stagenames(ahst_row_3);
      size3(counter) = sizes(ahst_row_3);
      repstatus3(counter) = repstatus(ahst_row_3);
      entrystatus3(counter) = entrystage(ahst_row_3);
      
      if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
        transition_type(counter) = 4; // Fecundity
      } else if (size3(counter) == size2(counter)) {
        if (repstatus3(counter) > repstatus2(counter) ||
          entrystatus3(counter) < entrystatus2(counter)) {
          transition_type(counter) = 2; // Growth
        } else if (repstatus3(counter) < repstatus2(counter) ||
          entrystatus3(counter) > entrystatus2(counter)) {
          transition_type(counter) = 3; // Shrinkage
        } else {
          transition_type(counter) = 1; // Stasis
        }
      } else if (size3(counter) > size2(counter)) {
        transition_type(counter) = 2; // Growth
      } else if (size3(counter) < size2(counter)) {
        transition_type(counter) = 3; // Shrinkage
      }
      
      counter++;
    }
  }
  
  arma::uvec targetindices = find(ahsindexl > -1);
  arma::ivec ahsindex = ahsindexl.elem(targetindices);
  
  DataFrame output = DataFrame::create(Named("index") = ahsindex,
    _["transition"] = transition_type, _["stage3"] = longstages3,
    _["size3"] = size3, _["repstatus3"] = repstatus3,
    _["entrystatus3"] = entrystatus3, _["stage2"] = longstages2,
    _["size2"] = size2, _["repstatus2"] = repstatus2,
    _["entrystatus2"] = entrystatus2);
  
  return output;
}

//' Sum Positive and Negative LTRE or Elasticity Contributions
//' 
//' @name .demolition4
//' 
//' Function \code{demolition4()} takes \code{lefkoElas} and \code{lefkoLTRE}
//' inputs and returns a data frame summarizing their positive and negative
//' contributions.
//' 
//' @param cmats Any \code{lefkoElas} or \code{lefkoLTRE} object.
//' 
//' @return A data frame with three columns per contribution type, corresponding
//' to summed positive, negative, and total contributions, and rows
//' corresponding to each input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.demolition4)]]
Rcpp::List demolition4 (List cmats) {
  bool ltre_input {false};
  bool sparse_input {false};
  
  StringVector cmats_class = as<StringVector>(cmats.attr("class"));
  for (int i = 0; i < cmats_class.length(); i++) {
    if (stringcompare_hard(as<std::string>(cmats_class(i)), "lefkoLTRE")) {
      ltre_input = true;
    } else if (stringcompare_hard(as<std::string>(cmats_class(i)), "lefkoElas")) {
      ltre_input = false;
    } else {
      throw Rcpp::exception("Unrecognized input in summary.lefkoLTRE().", false);
    }
  }
  
  List output_df;
  StringVector varnames_;
  IntegerVector rowlabels;
  
  if (ltre_input) { 
    List cont_0 = as<List>(cmats(0));
    int conts_num = cmats.length() - 4;
    if (conts_num > 4) conts_num -= 2;
    int mats_num = cont_0.length();
    
    if (is<S4>(cont_0(0))) sparse_input = true;
    
    List temp_df ((conts_num * 4) + 1);
    rowlabels = seq(1, mats_num);
    temp_df(0) = rowlabels;
    StringVector varnames ((conts_num * 4) + 1);
    varnames(0) = "matrix";
    
    for (int i = 0; i < conts_num; i++) {
      List current_cont = as<List>(cmats(i));
      
      NumericVector pos_cont (mats_num);
      NumericVector neg_cont (mats_num);
      NumericVector abs_cont (mats_num);
      NumericVector tot_cont (mats_num);
      
      for (int j = 0; j < mats_num; j++) {
        if (!sparse_input) { 
          arma::mat current_mat = as<arma::mat>(current_cont(j));
          arma::uvec pos_vec = find(current_mat > 0.0);
          arma::uvec neg_vec = find(current_mat < 0.0);
          
          double pos_sum = accu(current_mat.elem(pos_vec));
          double neg_sum = accu(current_mat.elem(neg_vec));
          
          pos_cont(j) = pos_sum;
          neg_cont(j) = neg_sum;
          abs_cont(j) = pos_sum + (-1.0 * neg_sum);
          tot_cont(j) = pos_sum + neg_sum;
          
        } else {
          arma::sp_mat current_mat = as<arma::sp_mat>(current_cont(j));
          
          double pos_sum {0.0};
          double neg_sum {0.0};
          
          for (int k = 0; k < static_cast<int>(current_mat.n_elem); k++) {
            if (current_mat(k) > 0.0) {
              pos_sum += static_cast<double>(current_mat(k));
            } else if (current_mat(k) < 0.0) {
              neg_sum += static_cast<double>(current_mat(k));
            }
          }
          
          pos_cont(j) = pos_sum;
          neg_cont(j) = neg_sum;
          abs_cont(j) = pos_sum + (-1.0 * neg_sum);
          tot_cont(j) = pos_sum + neg_sum;
          
        }
      }
      
      if (i == 0) {
        varnames(1) = "means_positive";
        varnames(2) = "means_negative";
        varnames(3) = "means_abs_sum";
        varnames(4) = "means_total";
      }
      
      if (i == 1) {
        if (conts_num == 2) {
          varnames(5) = "sd_positive";
          varnames(6) = "sd_negative";
          varnames(7) = "sd_abs_sum";
          varnames(8) = "sd_total";
        } else {
          varnames(5) = "elas_positive";
          varnames(6) = "elas_negative";
          varnames(7) = "elas_abs_sum";
          varnames(8) = "elas_total";
        }
      }
      
      if (i == 2) {
        varnames(9) = "cv_positive";
        varnames(10) = "cv_negative";
        varnames(11) = "cv_abs_sum";
        varnames(12) = "cv_total";
      }
      
      if (i == 3) {
        varnames(13) = "cor_positive";
        varnames(14) = "cor_negative";
        varnames(15) = "cor_abs_sum";
        varnames(16) = "cor_total";
      }
      
      temp_df(1 + (4 * i) + 0) = pos_cont;
      temp_df(1 + (4 * i) + 1) = neg_cont;
      temp_df(1 + (4 * i) + 2) = abs_cont;
      temp_df(1 + (4 * i) + 3) = tot_cont;
    }
    output_df = temp_df;
    varnames_ = varnames;
  }
  
  StringVector rowlabels_str = as<StringVector>(rowlabels);
  StringVector out_class = {"data.frame"};
  output_df.attr("class") = out_class;
  output_df.attr("names") = varnames_;
  output_df.attr("row.names") = rowlabels_str;
  
  return output_df;
}

//' Creates Summary Data for Elasticity Matrix Inputs
//' 
//' Function \code{demolition3()} sums elasticity values from elasticity
//' matrices, and LTRE contributions from LTRE and sLTRE matrices, according to
//' the categories developed by functions \code{bambi2()} and \code{bambi3()}.
//' It requires \code{matrix} class inputs.
//' 
//' @name .demolition3
//' 
//' @param e_amat A single elasticity, LTRE, or sLTRE matrix of class
//' \code{matrix}.
//' @param bambesque This is the output from \code{bambi2()} or \code{bambi3()}
//' corresponding to the current lefkoMat object. The format is a data frame
//' giving the indices and characteristics of all predicted potential non-zero
//' elements in the supplied matrix.
//' @param amat_ The A matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions.
//' @param fmat_ The F matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions.
//' 
//' @return A list with two data frames, one showing the summed elasticities for
//' the historical matrix supplied (if supplied), and the other showing the
//' ahistorical summary of the historical matrix or the summed elasticities of
//' a supplied ahistorical elasticity matrix. Also includes sums of only the
//' positive elements and only the negative elements, in all cases.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.demolition3)]]
Rcpp::List demolition3(const arma::mat& e_amat, const DataFrame& bambesque,
  Nullable<Rcpp::NumericMatrix> amat_ = R_NilValue,
  Nullable<Rcpp::NumericMatrix> fmat_ = R_NilValue) {
  
  arma::uvec eindices = as<arma::uvec>(bambesque["index"]);
  arma::uvec categories = as<arma::uvec>(bambesque["transition"]);
  
  int e_amatsize = static_cast<int>(e_amat.n_elem);
  int e_amatrows = static_cast<int>(e_amat.n_rows);
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  arma::mat amat;
  arma::mat fmat;
  
  if (maxelem > e_amatsize) {
    throw Rcpp::exception("Supplied info does not correspond to input matrices.",
      false);
  }
  
  if (amat_.isNotNull() && fmat_.isNotNull()) {
    amat = Rcpp::as<arma::mat>(amat_);
    fmat = Rcpp::as<arma::mat>(fmat_);
    
  } else {
    arma::mat amat1(e_amatrows, e_amatrows, fill::ones);
    arma::mat fmat1(e_amatrows, e_amatrows, fill::zeros);
    amat = amat1;
    fmat = fmat1;
    
    arma::uvec fec_trans4 = find(categories == 4);
    if (fec_trans4.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans4.n_elem); i ++) {
        fmat(eindices(fec_trans4(i))) = 1;
      }
    }
    
    arma::uvec fec_trans20 = find(categories == 20);
    if (fec_trans20.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans20.n_elem); i ++) {
        fmat(eindices(fec_trans20(i))) = 1;
      }
    }
    
    arma::uvec fec_trans21 = find(categories == 21);
    if (fec_trans21.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans21.n_elem); i ++) {
        fmat(eindices(fec_trans21(i))) = 1;
      }
    }
    
    arma::uvec fec_trans22 = find(categories == 22);
    if (fec_trans22.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans22.n_elem); i ++) {
        fmat(eindices(fec_trans22(i))) = 1;
      }
    }
    
    arma::uvec fec_trans26 = find(categories == 26);
    if (fec_trans26.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans26.n_elem); i ++) {
        fmat(eindices(fec_trans26(i))) = 1;
      }
    }
  }
  
  DataFrame histout;
  DataFrame ahistout;
  
  StringVector histcats {"Full stasis", "Stasis to growth", "Full growth",
    "Growth to stasis", "Stasis to shrinkage", "Full shrinkage",
    "Shrinkage to stasis", "Growth to shrinkage", "Shrinkage to growth",
    "Stasis to fecundity", "Growth to fecundity", "Shrinkage to fecundity",
    "Fecundity to stasis", "Fecundity to growth", "Fecundity to shrinkage",
    "Fecundity to fecundity"};
  arma::uvec histcatnums {10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23,
    24, 25, 26};
  arma::vec histsums(16);
  arma::vec histpos(16);
  arma::vec histneg(16);
  arma::vec hc_ahistsums(4);
  arma::vec hc_ahistpos(4);
  arma::vec hc_ahistneg(4);
  histsums.zeros();
  histpos.zeros();
  histneg.zeros();
  hc_ahistsums.zeros();
  hc_ahistpos.zeros();
  hc_ahistneg.zeros();
  
  StringVector ahistcats {"Stasis", "Growth", "Shrinkage", "Fecundity"};
  arma::uvec ahistcatnums {1, 2, 3, 4};
  arma::vec ahistsums(4);
  arma::vec ahistpos(4);
  arma::vec ahistneg(4);
  ahistsums.zeros();
  ahistpos.zeros();
  ahistneg.zeros();
  
  if (minindex > 9) {
    // Object minindex will only be > 9 if historical
    arma::vec size3 = as<arma::vec>(bambesque["size3"]);
    arma::vec size2 = as<arma::vec>(bambesque["size2"]);
    arma::vec size1 = as<arma::vec>(bambesque["size1"]);
    
    arma::uvec repstatus3 = as<arma::uvec>(bambesque["repstatus3"]);
    arma::uvec repstatus2 = as<arma::uvec>(bambesque["repstatus2"]);
    arma::uvec repstatus1 = as<arma::uvec>(bambesque["repstatus1"]);
    
    arma::uvec entrystatus3 = as<arma::uvec>(bambesque["entrystatus3"]);
    arma::uvec entrystatus2 = as<arma::uvec>(bambesque["entrystatus2"]);
    arma::uvec entrystatus1 = as<arma::uvec>(bambesque["entrystatus1"]);
    
    for (int i = 0; i < 16; i++) {
      arma::uvec currentguys = find(categories == histcatnums(i));
      int currentguysnem = static_cast<int>(currentguys.n_elem);
      
      if (histcatnums(i) == 20 || histcatnums(i) == 21 || histcatnums(i) == 22 ||
        histcatnums(i) == 26) { // Fecundity transitions
      
        // Splits transitions that are actually combos of fecundity and survival
        for (int j = 0; j < currentguysnem; j++) {
          int this_guy = eindices(currentguys(j));
          
          hc_ahistsums(3) += (e_amat(this_guy));
          histsums(i) += (e_amat(this_guy));
          
          if (e_amat(this_guy) > 0) {
            hc_ahistpos(3) += (e_amat(this_guy));
            histpos(i) += (e_amat(this_guy));
          } else if (e_amat(this_guy) < 0) {
            hc_ahistneg(3) += (e_amat(this_guy));
            histneg(i) += (e_amat(this_guy));
          }
        }
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 ||
          histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(2) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(2) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(2) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 ||
          histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(0) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(0) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(0) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 ||
          histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(1) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(1) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(1) += sum(all_es.elem(all_es_neg));
        }
      }
    }
    
    histout = DataFrame::create(Named("category") = histcats,
      _["elas"] = histsums, _["elas_pos"] = histpos, _["elas_neg"] = histneg);
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = hc_ahistsums, _["elas_pos"] = hc_ahistpos,
      _["elas_neg"] = hc_ahistneg);
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      arma::vec all_es = e_amat.elem(eindices(currentguys));
      arma::uvec all_es_pos = find(all_es > 0.0);
      arma::uvec all_es_neg = find(all_es < 0.0);
      int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
      int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
      
      double getoutofdodge = sum(all_es);
      ahistsums(i) += getoutofdodge;
      if (all_es_pos_num > 0) {
        ahistpos(i) += sum(all_es.elem(all_es_pos));
      }
      if (all_es_neg_num > 0) {
        ahistneg(i) += sum(all_es.elem(all_es_neg));
      }
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = ahistsums, _["elas_pos"] = ahistpos,
      _["elas_neg"] = ahistneg);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}

//' Creates Summary Data for Elasticity Matrix Inputs
//' 
//' Function \code{demolition3sp()} sums elasticity values from elasticity
//' matrices, and LTRE contributions from LTRE and sLTRE matrices, according to
//' the categories developed by functions \code{bambi2()} and \code{bambi3()}.
//' It requires \code{dgCMatrix} class inputs.
//' 
//' @name .demolition3sp
//' 
//' @param e_amat A single elasticity, LTRE, or sLTRE matrix of class
//' \code{dgCMatrix}.
//' @param bambesque This is the output from \code{bambi2()} or \code{bambi3()}
//' corresponding to the current lefkoMat object. The format is a data frame
//' giving the indices and characteristics of all predicted potential non-zero
//' elements in the supplied matrix.
//' @param amat_ The A matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions. Must also be of class \code{dgCMatrix}.
//' @param fmat_ The F matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions. Must also be of class \code{dgCMatrix}.
//' 
//' @return A list with two data frames, one showing the summed elasticities for
//' the historical matrix supplied (if supplied), and the other showing the
//' ahistorical summary of the historical matrix or the summed elasticities of
//' a supplied ahistorical elasticity matrix. Also includes sums of only the
//' positive elements and only the negative elements, in all cases.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.demolition3sp)]]
Rcpp::List demolition3sp(const arma::sp_mat& e_amat, const DataFrame& bambesque,
  Nullable<arma::sp_mat> amat_ = R_NilValue,
  Nullable<arma::sp_mat> fmat_ = R_NilValue) {
  
  arma::uvec eindices = as<arma::uvec>(bambesque["index"]);
  arma::uvec categories = as<arma::uvec>(bambesque["transition"]);
  
  int e_amatsize = static_cast<int>(e_amat.n_elem);
  int e_amatrows = static_cast<int>(e_amat.n_rows);
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  arma::sp_mat amat;
  arma::sp_mat fmat;
  
  if (maxelem > e_amatsize) {
    throw Rcpp::exception("Supplied info does not correspond to input matrices.",
      false);
  }
  
  if (amat_.isNotNull() && fmat_.isNotNull()) {
    arma::sp_mat amat_temp = as<arma::sp_mat>(amat_);
    arma::sp_mat fmat_temp = as<arma::sp_mat>(fmat_);
    
    amat = amat_temp;
    fmat = fmat_temp;
    
  } else {
    arma::sp_mat fmat1(e_amatrows, e_amatrows);
    fmat = fmat1;
    
    arma::uvec fec_trans4 = find(categories == 4);
    if (fec_trans4.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans4.n_elem); i ++) {
        fmat(eindices(fec_trans4(i))) = 1;
      }
    }
    arma::uvec fec_trans20 = find(categories == 20);
    if (fec_trans20.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans20.n_elem); i ++) {
        fmat(eindices(fec_trans20(i))) = 1;
      }
    }
    arma::uvec fec_trans21 = find(categories == 21);
    if (fec_trans21.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans21.n_elem); i ++) {
        fmat(eindices(fec_trans21(i))) = 1;
      }
    }
    arma::uvec fec_trans22 = find(categories == 22);
    if (fec_trans22.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans22.n_elem); i ++) {
        fmat(eindices(fec_trans22(i))) = 1;
      }
    }
    arma::uvec fec_trans26 = find(categories == 26);
    if (fec_trans26.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans26.n_elem); i ++) {
        fmat(eindices(fec_trans26(i))) = 1;
      }
    }
  }
  
  DataFrame histout;
  DataFrame ahistout;
  
  StringVector histcats {"Full stasis", "Stasis to growth", "Full growth",
    "Growth to stasis", "Stasis to shrinkage", "Full shrinkage",
    "Shrinkage to stasis", "Growth to shrinkage", "Shrinkage to growth",
    "Stasis to fecundity", "Growth to fecundity", "Shrinkage to fecundity",
    "Fecundity to stasis", "Fecundity to growth", "Fecundity to shrinkage",
    "Fecundity to fecundity"};
  arma::uvec histcatnums {10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23,
    24, 25, 26};
  arma::vec histsums(16);
  arma::vec histpos(16);
  arma::vec histneg(16);
  arma::vec hc_ahistsums(4);
  arma::vec hc_ahistpos(4);
  arma::vec hc_ahistneg(4);
  histsums.zeros();
  histpos.zeros();
  histneg.zeros();
  hc_ahistsums.zeros();
  hc_ahistpos.zeros();
  hc_ahistneg.zeros();
  
  StringVector ahistcats {"Stasis", "Growth", "Shrinkage", "Fecundity"};
  arma::uvec ahistcatnums {1, 2, 3, 4};
  arma::vec ahistsums(4);
  arma::vec ahistpos(4);
  arma::vec ahistneg(4);
  ahistsums.zeros();
  ahistpos.zeros();
  ahistneg.zeros();
  
  if (minindex > 9) {
    // Object minindex will only be > 9 if historical
    arma::vec size3 = as<arma::vec>(bambesque["size3"]);
    arma::vec size2 = as<arma::vec>(bambesque["size2"]);
    arma::vec size1 = as<arma::vec>(bambesque["size1"]);
    
    arma::uvec repstatus3 = as<arma::uvec>(bambesque["repstatus3"]);
    arma::uvec repstatus2 = as<arma::uvec>(bambesque["repstatus2"]);
    arma::uvec repstatus1 = as<arma::uvec>(bambesque["repstatus1"]);
    
    arma::uvec entrystatus3 = as<arma::uvec>(bambesque["entrystatus3"]);
    arma::uvec entrystatus2 = as<arma::uvec>(bambesque["entrystatus2"]);
    arma::uvec entrystatus1 = as<arma::uvec>(bambesque["entrystatus1"]);
    
    for (int i = 0; i < 16; i++) {
      arma::uvec currentguys = find(categories == histcatnums(i));
      int currentguysnem = static_cast<int>(currentguys.n_elem);
      
      if (histcatnums(i) == 20 || histcatnums(i) == 21 || histcatnums(i) == 22 ||
        histcatnums(i) == 26) { // Fecundity transitions
      
        // Splits transitions that are actually combos of fecundity and survival
        for (int j = 0; j < currentguysnem; j++) {
          int this_guy = eindices(currentguys(j));
          
          hc_ahistsums(3) += (e_amat(this_guy));
          histsums(i) += (e_amat(this_guy));
          
          if (e_amat(this_guy) > 0) {
            hc_ahistpos(3) += (e_amat(this_guy));
            histpos(i) += (e_amat(this_guy));
          } else if (e_amat(this_guy) < 0) {
            hc_ahistneg(3) += (e_amat(this_guy));
            histneg(i) += (e_amat(this_guy));
          }
        }
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 ||
          histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
        for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
          all_es(m) = e_amat(eindices(currentguys(m)));
        }
        
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(2) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(2) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(2) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 ||
          histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
        for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
          all_es(m) = e_amat(eindices(currentguys(m)));
        }
        
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(0) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(0) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(0) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 ||
          histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
        for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
          all_es(m) = e_amat(eindices(currentguys(m)));
        }
        
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(1) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(1) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(1) += sum(all_es.elem(all_es_neg));
        }
      }
    }
    
    histout = DataFrame::create(Named("category") = histcats,
      _["elas"] = histsums, _["elas_pos"] = histpos, _["elas_neg"] = histneg);
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = hc_ahistsums, _["elas_pos"] = hc_ahistpos,
      _["elas_neg"] = hc_ahistneg);
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      
      arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
      for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
        all_es(m) = e_amat(eindices(currentguys(m)));
      }
        
      arma::uvec all_es_pos = find(all_es > 0.0);
      arma::uvec all_es_neg = find(all_es < 0.0);
      int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
      int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
      
      double getoutofdodge = sum(all_es);
      ahistsums(i) += getoutofdodge;
      if (all_es_pos_num > 0) {
        ahistpos(i) += sum(all_es.elem(all_es_pos));
      }
      if (all_es_neg_num > 0) {
        ahistneg(i) += sum(all_es.elem(all_es_neg));
      }
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = ahistsums, _["elas_pos"] = ahistpos,
      _["elas_neg"] = ahistneg);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}

//' Estimate Actual or Deterministic Population Growth Rate
//' 
//' Function \code{lambda3()} is a generic function that returns the dominant
//' eigenvalue of a matrix, set of dominant eigenvalues of a set of matrices,
//' set of dominant eigenvalues for a \code{lefkoMat} object, or actual
//' \eqn{\lambda} in each year in a \code{lefkoProj} object. It can handle
//' large and sparse matrices supplied as \code{lefkoMat} objects or as
//' individual matrices, and can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices, and general
//' projetions.
//' 
//' @param mpm A \code{lefkoMat} object, a \code{lefkoMatList} object, a list of
//' projection matrices, a \code{lefkoProj} object, or a single projection
//' matrix in either standard or sparse format.
//' @param force_sparse A logical value or string detailing whether to force
//' sparse matrix encoding for simple matrix input. Defaults to \code{"auto"},
//' which only forces sparse matrix coding if simple matrices are input that are
//' both sparse (i.e, percentage of matrix elements that are non-zero <= 50%)
//' and have more than 20 rows. Can also be set to \code{"yes"}, \code{"no"},
//' \code{TRUE}, or \code{FALSE}. Note that sparse matrix coding is always used
//' for \code{lefkoMat} objects with matrices in sparse format (class
//' \code{dgCMatrix}). Ignored with \code{lefkoProj} objects.
//' 
//' @return The value returned depends on the class of the \code{mpm} argument.
//' If a \code{lefkoMat} object is provided, then this function will return the
//' \code{labels} data frame with a new column named \code{lambda} showing the
//' dominant eigenvalues for each matrix. If a \code{lefkoMatList} object is
//' provided, then a two element list in which the first element is a vector
//' composed of the mean lambda values of each \code{lefkoMat} element within
//' the list is provided, and the second element is a list of \code{lefkoMat}
//' lambda summaries as previously described. If a list of matrices is provided,
//' then this function will produce a numeric vector with the dominant
//' eigenvalues provided in order of matrix. If a single matrix is provided,
//' then this function will return the dominant eigenvalue of that matrix. Only
//' the largest real parts of the eigenvalues are returned.
//' 
//' If a \code{lefkoProj} object is provided, then the output consists of a list
//' with three elements. The second and third elements are lists of matrices
//' with each lower-level list elements corresponding to \code{labels} rows,
//' and matrices within these lists showing the actual \eqn{\lambda} and
//' \code{log} \eqn{\lambda} for each consecutive year or time index (columns)
//' within each replicate (row).
//' 
//' @seealso \code{\link{slambda3}()}
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
//' ehrlen3mean <- lmean(ehrlen3)
//' lambda3(ehrlen3mean)
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
//' # Here we use supplemental() to provide overwrite and reproductive info
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
//' lambda3(cypmatrix2r)
//' 
//' @export lambda3
// [[Rcpp::export(lambda3)]]
Rcpp::RObject lambda3(RObject& mpm, Nullable<RObject> force_sparse = R_NilValue) {
  
  RObject output;
  
  int sparse_check {2};
  bool matrix_class_input {false};
  bool lefkoMatList_input {false};
  //bool lefkoMat_input {false};
  
  if (force_sparse.isNotNull()) {
    if (is<LogicalVector>(force_sparse)) {
      LogicalVector sparse_ = as<LogicalVector>(force_sparse);
      
      if (static_cast<bool>(sparse_(0))) {
        sparse_check = 1; // Forced sparse matrix
      } else {
        sparse_check = 0; // Not forced sparse matrix
      }
    } else if (is<StringVector>(force_sparse)) {
      StringVector sparse_ = as<StringVector>(force_sparse);
      
      if (stringcompare_simple(String(sparse_(0)), "au", true)) {
        sparse_check = 2; // Auto decision
        
      } else if (stringcompare_simple(String(sparse_(0)), "y", true) ||
          stringcompare_simple(String(sparse_(0)), "t", true)) {
        sparse_check = 1; // Forced sparse matrix
        
      } else if (stringcompare_simple(String(sparse_(0)), "n", true) ||
          stringcompare_simple(String(sparse_(0)), "f", true)) {
        sparse_check = 0; // Not forced sparse matrix
      } else {  
        throw Rcpp::exception("Value entered for argument sparse not understood.",
          false);
      }
    } else {
      throw Rcpp::exception("Value entered for argument sparse not understood.",
        false);
    }
  }
  
  if (is<List>(mpm)) {
    List mpm_ = as<List>(mpm);
    
    CharacterVector mpm_class = as<CharacterVector>(mpm_.attr("class"));
    for (int i = 0; i < static_cast<int>(mpm_class.length()); i++) {
      if (stringcompare_hard(as<std::string>(mpm_class(i)), "lefkoMatList")) lefkoMatList_input = true;
    }
    
    if (lefkoMatList_input) {
      int lefkoMatList_length = static_cast<int>(mpm_.length());
      NumericVector all_lambdas (lefkoMatList_length);
      List catch_all (lefkoMatList_length);
      
      for (int i = 0; i < lefkoMatList_length; i++) {
        List current_mpm = as<List>(mpm_(i));
        CharacterVector mpm_names;
        if (current_mpm.hasAttribute("names")) mpm_names = current_mpm.attr("names");
        int no_mpm_names = mpm_names.length();
        
        bool A_check {false};
        bool ps_check {false};
        bool labels_check {false};
        for (int i = 0; i < no_mpm_names; i++) {
          if (stringcompare_simple(as<std::string>(mpm_names(i)), "A", false)) A_check = true;
          if (stringcompare_simple(as<std::string>(mpm_names(i)), "labels", false)) labels_check = true;
          if (stringcompare_simple(as<std::string>(mpm_names(i)), "pop_size", false)) ps_check = true;
        }
        
        if (!labels_check) {
          // List of matrices input
          if (is<NumericMatrix>(current_mpm(0))) { 
            matrix_class_input = true;
          } else if (is<S4>(current_mpm(0))) { 
            matrix_class_input = false;
          } else {
            throw Rcpp::exception("Object mpm list structure is not recognized.", false);
          }
          
          int no_matrices = current_mpm.length();
          
          if (matrix_class_input && sparse_check == 2) {
            arma::mat a1 = current_mpm[0];
            int mat_rows = a1.n_rows;
            int mat_cols = a1.n_cols;
            int total_elems = mat_rows * mat_cols;
            
            arma::uvec nonzeros = find(a1);
            int no_nonzeros = static_cast<int>(nonzeros.n_elem);
            
            double density = static_cast<double>(no_nonzeros) / static_cast<double>(total_elems);
            
            if (density <= 0.5 && total_elems > 399) {
              sparse_check = 1;
            } else {
              sparse_check = 0;
            }
          }
          
          NumericVector lambda_prog (no_matrices);
          
          for (int i = 0; i < no_matrices; i++) {
            if (sparse_check == 0 && matrix_class_input) {
              arma::mat Amat = as<arma::mat>(current_mpm(i));
              
              arma::cx_vec Aeigval;
              arma::cx_mat Aeigvecl;
              arma::cx_mat Aeigvecr;
              
              eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
              
              arma::vec all_eigenvalues = real(Aeigval);
              double maxval = max(all_eigenvalues);
              
              arma::uvec max_elems = find(all_eigenvalues == maxval);
              if (max_elems.n_elem == 0) {
                throw Rcpp::exception("Eigenanalysis failed.", false);
              }
              
              arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
              arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
              
              if (pos_max_eigvals.n_elem > 0) {
                lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
              } else {
                lambda_prog(i) = 0.0;
                Rf_warningcall(R_NilValue,
                  "A matrix with an eigenvalue of 0 has been detected.");
              }
            } else {
              arma::sp_mat spAmat;
              
              if (matrix_class_input) { 
                arma::sp_mat spAmat_(as<arma::mat>(current_mpm(i)));
                spAmat = spAmat_;
              } else { 
                spAmat = as<arma::sp_mat>(current_mpm(i));
              }
              
              arma::cx_vec Aeigval;
              arma::cx_mat Aeigvecr;
              
              eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
              
              arma::vec all_eigenvalues = real(Aeigval);
              double maxval = max(all_eigenvalues);
              
              arma::uvec max_elems = find(all_eigenvalues == maxval);
              if (max_elems.n_elem == 0) {
                throw Rcpp::exception("Eigen analysis failed.", false);
              }
              arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
              arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
              
              if (pos_max_eigvals.n_elem > 0) {
                lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
              } else {
                lambda_prog(i) = 0.0;
                Rf_warningcall(R_NilValue,
                  "A matrix with an eigenvalue of 0 has been detected.");
              }
            }
          }
          catch_all(i) = lambda_prog;
          all_lambdas(i) = mean(lambda_prog);
        } else {
          
          DataFrame labels = as<DataFrame>(current_mpm["labels"]);
          CharacterVector l_pop = labels["pop"];
          CharacterVector l_patch = labels["patch"];
          int l_length = static_cast<int>(labels.length());
          int num_poppatchyears = static_cast<int>(l_patch.length());
          
          if (A_check) {
            // lefkoMat input
            List A_list = current_mpm["A"];
            if (is<NumericMatrix>(A_list(0))) { 
              matrix_class_input = true;
            } else if (is<S4>(A_list(0))) { 
              matrix_class_input = false;
            } else {
              throw Rcpp::exception("Object mpm does not appear to contain matrices.", false);
            }
            
            int no_matrices = A_list.length();
            
            if (matrix_class_input && sparse_check == 2) {
              arma::mat a1 = A_list[0];
              int mat_rows = a1.n_rows;
              int mat_cols = a1.n_cols;
              int total_elems = mat_rows * mat_cols;
              
              arma::uvec nonzeros = find(a1);
              int no_nonzeros = static_cast<int>(nonzeros.n_elem);
              
              double density = static_cast<double>(no_nonzeros) /
                static_cast<double>(total_elems);
              
              if (density <= 0.5 && total_elems > 399) {
                sparse_check = 1;
              } else {
                sparse_check = 0;
              }
            }
            
            NumericVector lambda_prog (no_matrices);
            
            for (int i = 0; i < no_matrices; i++) {
              if (sparse_check == 0 && matrix_class_input) {
                arma::mat Amat = as<arma::mat>(A_list(i));
                
                arma::cx_vec Aeigval;
                arma::cx_mat Aeigvecl;
                arma::cx_mat Aeigvecr;
                
                eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
                
                arma::vec all_eigenvalues = real(Aeigval);
                double maxval = max(all_eigenvalues);
                
                arma::uvec max_elems = find(all_eigenvalues == maxval);
                if (max_elems.n_elem == 0) {
                  throw Rcpp::exception("Eigenanalysis failed.", false);
                }
                
                arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
                arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
                
                if (pos_max_eigvals.n_elem > 0) {
                  lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
                } else {
                  lambda_prog(i) = 0.0;
                  Rf_warningcall(R_NilValue,
                    "A matrix with an eigenvalue of 0 has been detected.");
                }
              } else {
                arma::sp_mat spAmat;
                
                if (matrix_class_input) {
                  arma::sp_mat spAmat_(as<arma::mat>(A_list(i)));
                  spAmat = spAmat_;
                } else { 
                  spAmat = as<arma::sp_mat>(A_list(i));
                }
                
                arma::cx_vec Aeigval;
                arma::cx_mat Aeigvecr;
                
                eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
                
                arma::vec all_eigenvalues = real(Aeigval);
                double maxval = max(all_eigenvalues);
                
                arma::uvec max_elems = find(all_eigenvalues == maxval);
                if (max_elems.n_elem == 0) {
                  throw Rcpp::exception("Eigen analysis failed.", false);
                }
                arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
                arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
                
                if (pos_max_eigvals.n_elem > 0) {
                  lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
                } else {
                  lambda_prog(i) = 0.0;
                  Rf_warningcall(R_NilValue,
                    "A matrix with an eigenvalue of 0 has been detected.");
                }
              }
            }
            
            DataFrame new_out;
            
            if (l_length == 3) {
              CharacterVector l_year2 = labels["year2"];
              
              new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
                _["year2"] = l_year2, _["lambda"] = lambda_prog);
            } else {
              new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
                _["lambda"] = lambda_prog);
            }
            catch_all(i) = new_out;
            all_lambdas(i) = mean(lambda_prog);
          } else if (ps_check) {
            // lefkoProj input
            List ps_list = current_mpm["pop_size"];
            int ps_list_length = static_cast<int>(ps_list.length());
            
            if (ps_list_length != num_poppatchyears) {
              throw Rcpp::exception("Structure of lefkoProj object appears broken.", false);
            }
            
            Rcpp::List fine_lambdas (ps_list_length);
            Rcpp::List fine_log_lambdas (ps_list_length);
            Rcpp::List lambda_summary (ps_list_length);
            
            // Summary labels construction
            int full_summary_labels_rows {0};
            for (int i = 0; i < ps_list_length; i++) {
              NumericMatrix current_ps_mat = as<NumericMatrix>(ps_list[i]);
              int ps_mat_rows = static_cast<int>(current_ps_mat.nrow());
              
              full_summary_labels_rows += ps_mat_rows;
            }
            
            IntegerVector sumlab_rownum = seq(1, full_summary_labels_rows);
            StringVector sumlab_pop (full_summary_labels_rows);
            StringVector sumlab_patch (full_summary_labels_rows);
            IntegerVector sumlab_repl (full_summary_labels_rows);
            NumericVector summary_mean_lambdas (full_summary_labels_rows);
            NumericVector summary_sd_lambdas (full_summary_labels_rows);
            NumericVector summary_mean_log_lambdas (full_summary_labels_rows);
            NumericVector summary_sd_log_lambdas (full_summary_labels_rows);
            IntegerVector summary_lambda_zeros (full_summary_labels_rows);
            IntegerVector summary_lambda_nonzeros (full_summary_labels_rows);
            
            int sumlab_rowtracker {0};
            
            for (int i = 0; i < ps_list_length; i++) {
              NumericMatrix current_ps_mat = as<NumericMatrix>(ps_list[i]);
              int ps_mat_rows = static_cast<int>(current_ps_mat.nrow());
              int ps_mat_cols = static_cast<int>(current_ps_mat.ncol());
              
              NumericMatrix lambda_mat (ps_mat_rows, (ps_mat_cols - 1));
              NumericMatrix log_lambda_mat (ps_mat_rows, (ps_mat_cols - 1));
              IntegerVector real_num_counts (ps_mat_rows);
              
              for (int repl = 0; repl < ps_mat_rows; repl++) {
                int current_real_num_count {0};
                int current_real_zero_count {0};
                int current_real_nonzero_count {0};
                
                
                for (int time = 1; time < ps_mat_cols; time++) {
                  if (!NumericVector::is_na(current_ps_mat(repl, (time - 1))) &&
                      !NumericVector::is_na(current_ps_mat(repl, time))) {
                    
                    double lam_numer = current_ps_mat(repl, time);
                    double lam_denom = current_ps_mat(repl, (time - 1));
                    
                    if (lam_denom > 0.0) {
                      lambda_mat(repl, (time - 1)) = lam_numer / lam_denom;
                      log_lambda_mat(repl, (time - 1)) = log(lambda_mat(repl, (time - 1)));
                      
                      current_real_num_count++;
                      current_real_nonzero_count++;
                    } else if (lam_denom == 0.0) {
                      lambda_mat(repl, (time - 1)) = R_NaN;
                      log_lambda_mat(repl, (time - 1)) = R_NaN;
                      
                      current_real_zero_count++;
                    } else if (lam_numer == 0.0) {
                      lambda_mat(repl, (time - 1)) = 0.0;
                      log_lambda_mat(repl, (time - 1)) = R_NaN;
                      
                      current_real_num_count++;
                      current_real_zero_count++;
                    }
                  } else {
                    lambda_mat(repl, (time - 1)) = NA_REAL;
                    log_lambda_mat(repl, (time - 1)) = NA_REAL;
                  }
                }
                real_num_counts(repl) = current_real_num_count;
                summary_lambda_zeros(sumlab_rowtracker) = current_real_zero_count;
                summary_lambda_nonzeros(sumlab_rowtracker) = current_real_nonzero_count;
                
                NumericVector core_row = lambda_mat(repl, _);
                NumericVector core_log_row = log_lambda_mat(repl, _);
                
                double s1 {0.0};
                double s2 {0.0};
                double s1_log {0.0};
                double s2_log {0.0};
                for (int j = 0; j < current_real_nonzero_count; j++) {
                  s1 += core_row(j) ;
                  s2 += (core_row(j) * core_row(j)) ;
                  
                  if (current_real_zero_count == 0) {
                    s1_log += core_log_row(j);
                    s2_log += (core_log_row(j) * core_log_row(j));
                  } 
                }
                
                double used_denominator = static_cast<double>(current_real_nonzero_count);
                summary_mean_lambdas(sumlab_rowtracker) = s1 / (used_denominator);
                summary_mean_log_lambdas(sumlab_rowtracker) = s1_log / (used_denominator);
                summary_sd_lambdas(sumlab_rowtracker) = sqrt((used_denominator) / (used_denominator - 1.0)) *
                    sqrt((s2 / used_denominator) - ((s1 / used_denominator) * (s1 / used_denominator)));
                summary_sd_log_lambdas(sumlab_rowtracker) = sqrt((used_denominator) / (used_denominator - 1.0)) *
                    sqrt((s2_log / used_denominator) - ((s1_log / used_denominator) * (s1_log / used_denominator)));
                
                if (current_real_zero_count > 0) summary_mean_log_lambdas(sumlab_rowtracker) = R_NaN;
                
                sumlab_repl(sumlab_rowtracker) = repl;
                sumlab_pop(sumlab_rowtracker) = l_pop(i);
                sumlab_patch(sumlab_rowtracker) = l_patch(i);
                
                sumlab_rowtracker++;
              }
              
              fine_lambdas(i) = lambda_mat;
              fine_log_lambdas(i) = log_lambda_mat;
            }
            
            DataFrame summary_out;
            summary_out = DataFrame::create(_["index"] = sumlab_rownum,
              _["pop"] = sumlab_pop, _["patch"] = sumlab_patch,
              _["lambda_mean"] = summary_mean_lambdas, _["lambda_sd"] = summary_sd_lambdas,
              _["log_lambda_mean"] = summary_mean_log_lambdas,
              _["log_lambda_sd"] = summary_sd_log_lambdas,
              _["lambda_zeros"] = summary_lambda_zeros);
            
            List new_out = List::create(_["summary"] = summary_out,
              _["lambda"] = fine_lambdas, _["log_lambda"] = fine_log_lambdas);
            catch_all(i) = new_out;
            all_lambdas(i) = mean(summary_mean_lambdas);
            
          } else {
            throw Rcpp::exception("Object mpm does not appear to be an appropriate MPM, list of matrices, lefkoProj object, or matrix.",
              false);
          }
        }
      }
      
      List general_output (2);
      general_output(0) = all_lambdas;
      general_output(1) = catch_all;
      CharacterVector go_names = {"lambdas", "summaries"};
      general_output.attr("names") = go_names;
      output = general_output;
    } else {
      CharacterVector mpm_names;
      if (mpm_.hasAttribute("names")) mpm_names = mpm_.attr("names");
      int no_mpm_names = mpm_names.length();
      
      bool A_check {false};
      bool ps_check {false};
      bool labels_check {false};
      for (int i = 0; i < no_mpm_names; i++) {
        if (stringcompare_simple(as<std::string>(mpm_names(i)), "A", false)) A_check = true;
        if (stringcompare_simple(as<std::string>(mpm_names(i)), "labels", false)) labels_check = true;
        if (stringcompare_simple(as<std::string>(mpm_names(i)), "pop_size", false)) ps_check = true;
      }
      
      if (!labels_check) {
        // List of matrices input
        if (is<NumericMatrix>(mpm_(0))) { 
          matrix_class_input = true;
        } else if (is<S4>(mpm_(0))) { 
          matrix_class_input = false;
        } else {
          throw Rcpp::exception("Object mpm list structure is not recognized.", false);
        }
        
        int no_matrices = mpm_.length();
        
        if (matrix_class_input && sparse_check == 2) {
          arma::mat a1 = mpm_[0];
          int mat_rows = a1.n_rows;
          int mat_cols = a1.n_cols;
          int total_elems = mat_rows * mat_cols;
          
          arma::uvec nonzeros = find(a1);
          int no_nonzeros = static_cast<int>(nonzeros.n_elem);
          
          double density = static_cast<double>(no_nonzeros) / static_cast<double>(total_elems);
          
          if (density <= 0.5 && total_elems > 399) {
            sparse_check = 1;
          } else {
            sparse_check = 0;
          }
        }
        
        NumericVector lambda_prog (no_matrices);
        
        for (int i = 0; i < no_matrices; i++) {
          if (sparse_check == 0 && matrix_class_input) {
            arma::mat Amat = as<arma::mat>(mpm_(i));
            
            arma::cx_vec Aeigval;
            arma::cx_mat Aeigvecl;
            arma::cx_mat Aeigvecr;
            
            eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
            
            arma::vec all_eigenvalues = real(Aeigval);
            double maxval = max(all_eigenvalues);
            
            arma::uvec max_elems = find(all_eigenvalues == maxval);
            if (max_elems.n_elem == 0) {
              throw Rcpp::exception("Eigenanalysis failed.", false);
            }
            
            arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
            arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
            
            if (pos_max_eigvals.n_elem > 0) {
              lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
            } else {
              lambda_prog(i) = 0.0;
              Rf_warningcall(R_NilValue,
                "A matrix with an eigenvalue of 0 has been detected.");
            }
          } else {
            arma::sp_mat spAmat;
            
            if (matrix_class_input) { 
              arma::sp_mat spAmat_(as<arma::mat>(mpm_(i)));
              spAmat = spAmat_;
            } else { 
              spAmat = as<arma::sp_mat>(mpm_(i));
            }
            
            arma::cx_vec Aeigval;
            arma::cx_mat Aeigvecr;
            
            eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
            
            arma::vec all_eigenvalues = real(Aeigval);
            double maxval = max(all_eigenvalues);
            
            arma::uvec max_elems = find(all_eigenvalues == maxval);
            if (max_elems.n_elem == 0) {
              throw Rcpp::exception("Eigen analysis failed.", false);
            }
            arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
            arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
            
            if (pos_max_eigvals.n_elem > 0) {
              lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
            } else {
              lambda_prog(i) = 0.0;
              Rf_warningcall(R_NilValue,
                "A matrix with an eigenvalue of 0 has been detected.");
            }
          }
        }
        output = lambda_prog;
        
      } else {
        
        DataFrame labels = as<DataFrame>(mpm_["labels"]);
        CharacterVector l_pop = labels["pop"];
        CharacterVector l_patch = labels["patch"];
        int l_length = static_cast<int>(labels.length());
        int num_poppatchyears = static_cast<int>(l_patch.length());
        
        if (A_check) {
          // lefkoMat input
          List A_list = mpm_["A"];
          if (is<NumericMatrix>(A_list(0))) { 
            matrix_class_input = true;
          } else if (is<S4>(A_list(0))) { 
            matrix_class_input = false;
          } else {
            throw Rcpp::exception("Object mpm does not appear to contain matrices.", false);
          }
          
          int no_matrices = A_list.length();
          
          if (matrix_class_input && sparse_check == 2) {
            arma::mat a1 = A_list[0];
            int mat_rows = a1.n_rows;
            int mat_cols = a1.n_cols;
            int total_elems = mat_rows * mat_cols;
            
            arma::uvec nonzeros = find(a1);
            int no_nonzeros = static_cast<int>(nonzeros.n_elem);
            
            double density = static_cast<double>(no_nonzeros) /
              static_cast<double>(total_elems);
            
            if (density <= 0.5 && total_elems > 399) {
              sparse_check = 1;
            } else {
              sparse_check = 0;
            }
          }
          
          NumericVector lambda_prog (no_matrices);
          
          for (int i = 0; i < no_matrices; i++) {
            if (sparse_check == 0 && matrix_class_input) {
              arma::mat Amat = as<arma::mat>(A_list(i));
              
              arma::cx_vec Aeigval;
              arma::cx_mat Aeigvecl;
              arma::cx_mat Aeigvecr;
              
              eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
              
              arma::vec all_eigenvalues = real(Aeigval);
              double maxval = max(all_eigenvalues);
              
              arma::uvec max_elems = find(all_eigenvalues == maxval);
              if (max_elems.n_elem == 0) {
                throw Rcpp::exception("Eigenanalysis failed.", false);
              }
              
              arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
              arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
              
              if (pos_max_eigvals.n_elem > 0) {
                lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
              } else {
                lambda_prog(i) = 0.0;
                Rf_warningcall(R_NilValue,
                  "A matrix with an eigenvalue of 0 has been detected.");
              }
            } else {
              arma::sp_mat spAmat;
              
              if (matrix_class_input) {
                arma::sp_mat spAmat_(as<arma::mat>(A_list(i)));
                spAmat = spAmat_;
              } else { 
                spAmat = as<arma::sp_mat>(A_list(i));
              }
              
              arma::cx_vec Aeigval;
              arma::cx_mat Aeigvecr;
              
              eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
              
              arma::vec all_eigenvalues = real(Aeigval);
              double maxval = max(all_eigenvalues);
              
              arma::uvec max_elems = find(all_eigenvalues == maxval);
              if (max_elems.n_elem == 0) {
                throw Rcpp::exception("Eigen analysis failed.", false);
              }
              arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
              arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
              
              if (pos_max_eigvals.n_elem > 0) {
                lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
              } else {
                lambda_prog(i) = 0.0;
                Rf_warningcall(R_NilValue,
                  "A matrix with an eigenvalue of 0 has been detected.");
              }
            }
          }
          
          DataFrame new_out;
          
          if (l_length == 3) {
            CharacterVector l_year2 = labels["year2"];
            
            new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
              _["year2"] = l_year2, _["lambda"] = lambda_prog);
          } else {
            new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
              _["lambda"] = lambda_prog);
          }
          output = new_out;
          
        } else if (ps_check) {
          // lefkoProj input
          List ps_list = mpm_["pop_size"];
          int ps_list_length = static_cast<int>(ps_list.length());
          
          if (ps_list_length != num_poppatchyears) {
            throw Rcpp::exception("Structure of lefkoProj object appears broken.", false);
          }
          
          Rcpp::List fine_lambdas (ps_list_length);
          Rcpp::List fine_log_lambdas (ps_list_length);
          Rcpp::List lambda_summary (ps_list_length);
          
          // Summary labels construction
          int full_summary_labels_rows {0};
          for (int i = 0; i < ps_list_length; i++) {
            NumericMatrix current_ps_mat = as<NumericMatrix>(ps_list[i]);
            int ps_mat_rows = static_cast<int>(current_ps_mat.nrow());
            
            full_summary_labels_rows += ps_mat_rows;
          }
          
          IntegerVector sumlab_rownum = seq(1, full_summary_labels_rows);
          StringVector sumlab_pop (full_summary_labels_rows);
          StringVector sumlab_patch (full_summary_labels_rows);
          IntegerVector sumlab_repl (full_summary_labels_rows);
          NumericVector summary_mean_lambdas (full_summary_labels_rows);
          NumericVector summary_sd_lambdas (full_summary_labels_rows);
          NumericVector summary_mean_log_lambdas (full_summary_labels_rows);
          NumericVector summary_sd_log_lambdas (full_summary_labels_rows);
          IntegerVector summary_lambda_zeros (full_summary_labels_rows);
          IntegerVector summary_lambda_nonzeros (full_summary_labels_rows);
          
          int sumlab_rowtracker {0};
          
          for (int i = 0; i < ps_list_length; i++) {
            NumericMatrix current_ps_mat = as<NumericMatrix>(ps_list[i]);
            int ps_mat_rows = static_cast<int>(current_ps_mat.nrow());
            int ps_mat_cols = static_cast<int>(current_ps_mat.ncol());
            
            NumericMatrix lambda_mat (ps_mat_rows, (ps_mat_cols - 1));
            NumericMatrix log_lambda_mat (ps_mat_rows, (ps_mat_cols - 1));
            IntegerVector real_num_counts (ps_mat_rows);
            
            for (int repl = 0; repl < ps_mat_rows; repl++) {
              int current_real_num_count {0};
              int current_real_zero_count {0};
              int current_real_nonzero_count {0};
              
              
              for (int time = 1; time < ps_mat_cols; time++) {
                if (!NumericVector::is_na(current_ps_mat(repl, (time - 1))) &&
                    !NumericVector::is_na(current_ps_mat(repl, time))) {
                  
                  double lam_numer = current_ps_mat(repl, time);
                  double lam_denom = current_ps_mat(repl, (time - 1));
                  
                  if (lam_denom > 0.0) {
                    lambda_mat(repl, (time - 1)) = lam_numer / lam_denom;
                    log_lambda_mat(repl, (time - 1)) = log(lambda_mat(repl, (time - 1)));
                    
                    current_real_num_count++;
                    current_real_nonzero_count++;
                  } else if (lam_denom == 0.0) {
                    lambda_mat(repl, (time - 1)) = R_NaN;
                    log_lambda_mat(repl, (time - 1)) = R_NaN;
                    
                    current_real_zero_count++;
                  } else if (lam_numer == 0.0) {
                    lambda_mat(repl, (time - 1)) = 0.0;
                    log_lambda_mat(repl, (time - 1)) = R_NaN;
                    
                    current_real_num_count++;
                    current_real_zero_count++;
                  }
                } else {
                  lambda_mat(repl, (time - 1)) = NA_REAL;
                  log_lambda_mat(repl, (time - 1)) = NA_REAL;
                }
              }
              real_num_counts(repl) = current_real_num_count;
              summary_lambda_zeros(sumlab_rowtracker) = current_real_zero_count;
              summary_lambda_nonzeros(sumlab_rowtracker) = current_real_nonzero_count;
              
              NumericVector core_row = lambda_mat(repl, _);
              NumericVector core_log_row = log_lambda_mat(repl, _);
              
              double s1 {0.0};
              double s2 {0.0};
              double s1_log {0.0};
              double s2_log {0.0};
              for (int j = 0; j < current_real_nonzero_count; j++) {
                s1 += core_row(j) ;
                s2 += (core_row(j) * core_row(j)) ;
                
                if (current_real_zero_count == 0) {
                  s1_log += core_log_row(j);
                  s2_log += (core_log_row(j) * core_log_row(j));
                } 
              }
              
              double used_denominator = static_cast<double>(current_real_nonzero_count);
              summary_mean_lambdas(sumlab_rowtracker) = s1 / (used_denominator);
              summary_mean_log_lambdas(sumlab_rowtracker) = s1_log / (used_denominator);
              summary_sd_lambdas(sumlab_rowtracker) = sqrt((used_denominator) / (used_denominator - 1.0)) *
                  sqrt((s2 / used_denominator) - ((s1 / used_denominator) * (s1 / used_denominator)));
              summary_sd_log_lambdas(sumlab_rowtracker) = sqrt((used_denominator) / (used_denominator - 1.0)) *
                  sqrt((s2_log / used_denominator) - ((s1_log / used_denominator) * (s1_log / used_denominator)));
              
              if (current_real_zero_count > 0) summary_mean_log_lambdas(sumlab_rowtracker) = R_NaN;
              
              sumlab_repl(sumlab_rowtracker) = repl;
              sumlab_pop(sumlab_rowtracker) = l_pop(i);
              sumlab_patch(sumlab_rowtracker) = l_patch(i);
              
              sumlab_rowtracker++;
            }
            
            fine_lambdas(i) = lambda_mat;
            fine_log_lambdas(i) = log_lambda_mat;
          }
          
          DataFrame summary_out;
          summary_out = DataFrame::create(_["index"] = sumlab_rownum,
            _["pop"] = sumlab_pop, _["patch"] = sumlab_patch,
            _["lambda_mean"] = summary_mean_lambdas, _["lambda_sd"] = summary_sd_lambdas,
            _["log_lambda_mean"] = summary_mean_log_lambdas,
            _["log_lambda_sd"] = summary_sd_log_lambdas,
            _["lambda_zeros"] = summary_lambda_zeros);
          
          List new_out = List::create(_["summary"] = summary_out,
            _["lambda"] = fine_lambdas, _["log_lambda"] = fine_log_lambdas);
          output = new_out;
          
        } else {
          throw Rcpp::exception("Object mpm does not appear to be an appropriate MPM, list of matrices, lefkoProj object, or matrix.",
            false);
        }
      }
      
    }
    
  } else if(is<NumericMatrix>(mpm)) {
    // Single matrix input
    arma::mat mpm_ = as<arma::mat>(mpm);
    
    if (sparse_check == 2) {
      int mat_rows = mpm_.n_rows;
      int mat_cols = mpm_.n_cols;
      int total_elems = mat_rows * mat_cols;
      
      arma::uvec nonzeros = find(mpm_);
      int no_nonzeros = static_cast<int>(nonzeros.n_elem);
      
      double density = static_cast<double>(no_nonzeros) /
        static_cast<double>(total_elems);
      
      if (density <= 0.5 && total_elems > 399) {
        sparse_check = 1;
      } else {
        sparse_check = 0;
      }
    }
    
    NumericVector lambda_prog (1);
    if (sparse_check == 0) {
      arma::cx_vec Aeigval;
      arma::cx_mat Aeigvecl;
      arma::cx_mat Aeigvecr;
      
      eig_gen(Aeigval, Aeigvecl, Aeigvecr, mpm_);
      
      arma::vec all_eigenvalues = real(Aeigval);
      
      double maxval = max(all_eigenvalues);
      
      arma::uvec max_elems = find(all_eigenvalues == maxval);
      if (max_elems.n_elem == 0) {
        throw Rcpp::exception("Eigen analysis failed.", false);
      }
      
      arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
      arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
      
      if (pos_max_eigvals.n_elem > 0) {
        lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
      } else {
        lambda_prog(0) = 0.0;
        Rf_warningcall(R_NilValue,
          "A matrix with an eigenvalue of 0 has been detected.");
      }
    } else {
      arma::sp_mat spAmat(mpm_);
      
      arma::cx_vec Aeigval;
      arma::cx_mat Aeigvecr;
      
      eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
      
      arma::vec all_eigenvalues = real(Aeigval);
      
      double maxval = max(all_eigenvalues);
      
      arma::uvec max_elems = find(all_eigenvalues == maxval);
      if (max_elems.n_elem == 0) {
        throw Rcpp::exception("Eigen analysis failed.", false);
      }
      
      arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
      arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
      
      if (pos_max_eigvals.n_elem > 0) {
        lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
      } else {
        lambda_prog(0) = 0.0;
        Rf_warningcall(R_NilValue,
          "A matrix with an eigenvalue of 0 has been detected.");
      }
    }
    
    output = lambda_prog;
    
  } else if(is<S4>(mpm)) {
    // Single sparse matrix input
    arma::sp_mat spAmat = as<arma::sp_mat>(mpm);
    
    NumericVector lambda_prog (1);
    arma::cx_vec Aeigval;
    arma::cx_mat Aeigvecr;
    
    eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
    
    arma::vec all_eigenvalues = real(Aeigval);
    
    double maxval = max(all_eigenvalues);
    
    arma::uvec max_elems = find(all_eigenvalues == maxval);
    if (max_elems.n_elem == 0) {
      throw Rcpp::exception("Eigen analysis failed.", false);
    }
    
    arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
    arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
    
    if (pos_max_eigvals.n_elem > 0) {
      lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
    } else {
      lambda_prog(0) = 0.0;
      Rf_warningcall(R_NilValue,
        "A matrix with an eigenvalue of 0 has been detected.");
    }
    
    output = lambda_prog;
    
  } else {
    String eat_my_shorts = "Object mpm does not appear to be an appropriate ";
    eat_my_shorts += "MPM, list of matrices, lefkoProj object, or matrix.";
    
    throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
  }
  
  return output;
}

//' Arranges Matrix Elements in Order of Magnitude for Interpretation
//' 
//' Function \code{matrix_interp} summarizes matrices from \code{lefkoMat},
//' \code{lefkoSens}, \code{lefkoElas}, and \code{lefkoLTRE} objects in terms
//' of the magnitudes of their elements. It can also create ordered summaries of
//' standard matrices and sparse matrices.
//' 
//' @name matrix_interp
//' 
//' @param object A list object in one of \code{lefko3}'s output formats, or a
//' standard matrix or sparse matrix in \code{dgCMatrix} format. Standard
//' \code{lefko3} output formats include \code{lefkoMat}, \code{lefkoSens},
//' \code{lefkoElas}, and \code{lefkoLTRE} objects.
//' @param mat_chosen The number of the matrix to assess, within the appropriate
//' matrix list. See \code{Notes} for further details.
//' @param part An integer noting whether to provide assessments of which of the
//' main types of matrices to analyze. In a standard \code{lefkoMat} object, the
//' integers \code{1}, \code{2}, and \code{3} correspond to the \code{A},
//' \code{U}, and \code{F} lists, respectively. In \code{lefkoSens} and
//' \code{lefkoElas} objects, the integers \code{1} and \code{2} correspond to
//' the ahistorical matrix sets and the historical matrix sets, respectively.
//' In deterministic and stochastic \code{lefkoLTRE} objects, the integers
//' \code{1} and \code{2} correspond to the \code{cont_mean} and \code{cont_sd}
//' lists, respectively.
//' @param type An integer corresponding to the type of order summary, including
//' most to least positive (\code{1}), most to least negative (\code{2}), and
//' greatest to lowest absolute magnitude (\code{3}). Defaults to type \code{3}.
//' 
//' @return A data frame arranging all elements in the matrix chosen from
//' greatest and smallest. This can be a data frame of only positive elements,
//' of only negative elements, or all elements in order of absolute magnitude.
//' 
//' @section Notes:
//' Argument \code{mat_chosen} refers to the number of the matrix within the
//' list that it is held in. For example, if the function is applied to the
//' \code{cont_sd} portion of a stochastic LTRE, and there are four LTRE
//' matrices within that list element corresponding to three patch LTRE matrices
//' and one overall population-level LTRE matrix, then setting this value to
//' \code{4} would focus the function on the overall population-level LTRE
//' matrix associated with contributions of the standard deviations of elements.
//' This argument should be left blank if a standard matrix or sparse matrix is
//' input.
//' 
//' Huge sparse matrices may take more time to process than small, dense
//' matrices.
//' 
//' @examples
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
//' aaa <- ltre3(cypmatrix2r, stochastic = TRUE)
//' 
//' matrix_interp(aaa, mat_chosen = 1, part = 2, type = 3)
//' 
//' @export matrix_interp
// [[Rcpp::export(matrix_interp)]]
Rcpp::DataFrame matrix_interp (RObject object, int mat_chosen = 1,
  int part = 1, int type = 3) {
  
  Rcpp::NumericMatrix mat;
  arma::sp_mat smat;
  Rcpp::List object_list;
  bool sparse_mat {false};
  bool list_only {false};
  int list_type {0}; // 1 = lefkoMat, 2 = lefkoSens, 3 = lefkoElas, 4 = lefkoLTRE, 0 = other
  
  CharacterVector object_class;
  int class_num {0};
  
  Rcpp::DataFrame stageframe; 
  Rcpp::DataFrame hstages;
  Rcpp::DataFrame agestages;
  
  if (Rf_isMatrix(object)) {
    mat = as<NumericMatrix>(object);
    
    object_class = {"matrix", "array"};
    class_num = object_class.length();
    
  } else if (object.isS4()) {
    sparse_mat = true;
    
    S4 test_S4 = as<S4>(object);
    CharacterVector test_S4_class = test_S4.attr("class");
    bool found_dgc {false};
    for (int i_c = 0; i_c < static_cast<int>(test_S4_class.length()); i_c++) {
      if (stringcompare_hard(as<std::string>(test_S4_class(i_c)), "dgCMatrix")) found_dgc = true;
    }
    if (!found_dgc) throw Rcpp::exception("Object entered is of unknown class.", false);
    
    smat = as<arma::sp_mat>(object);
    
    object_class = {"dgCMatrix"};
    class_num = object_class.length();
    
  } else if (is<List>(object)) {
    list_only = true;
    
    object_list = as<List>(object);
    
    if (!object_list.hasAttribute("class")) {
      throw Rcpp::exception("Entered list object not recognized.", false);
    }
    
    object_class = object_list.attr("class");
    class_num = object_class.length();
    
    for (int i = 0; i < class_num; i++) {
      if (stringcompare_hard(as<std::string>(object_class(i)), "lefkoLTRE")) {
        stageframe = as<DataFrame>(object_list["ahstages"]); 
        hstages = as<DataFrame>(object_list["hstages"]);
        agestages = as<DataFrame>(object_list["agestages"]);
        
        list_type = 4;
        
        if (part == 1) {
          Rcpp::List mat_list = object_list["cont_mean"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        } else if (part == 2) {
          Rcpp::List mat_list = object_list["cont_sd"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        }
      } else if (stringcompare_hard(as<std::string>(object_class(i)), "lefkoElas")) {
        stageframe = as<DataFrame>(object_list["ahstages"]); 
        hstages = as<DataFrame>(object_list["hstages"]);
        agestages = as<DataFrame>(object_list["agestages"]);
        
        list_type = 3;
        
        if (part == 1) {
          Rcpp::List mat_list = object_list["ah_elasmats"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        } else if (part == 2) {
          Rcpp::List mat_list = object_list["h_elasmats"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        }
      } else if (stringcompare_hard(as<std::string>(object_class(i)), "lefkoSens")) {
        stageframe = as<DataFrame>(object_list["ahstages"]); 
        hstages = as<DataFrame>(object_list["hstages"]);
        agestages = as<DataFrame>(object_list["agestages"]);
        
        list_type = 2;
        
        if (part == 1) {
          Rcpp::List mat_list = object_list["ah_sensmats"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        } else if (part == 2) {
          Rcpp::List mat_list = object_list["h_sensmats"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        }
      } else if (stringcompare_hard(as<std::string>(object_class(i)), "lefkoMat")) {
        stageframe = as<DataFrame>(object_list["ahstages"]); 
        hstages = as<DataFrame>(object_list["hstages"]);
        agestages = as<DataFrame>(object_list["agestages"]);
        
        list_type = 1;
        
        if (part == 1) {
          Rcpp::List mat_list = object_list["A"];
          
          RObject the_chosen_one = as<RObject>(mat_list[mat_chosen - 1]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        } else if (part == 2) {
          Rcpp::List mat_list = object_list["U"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        } else if (part == 3) {
          Rcpp::List mat_list = object_list["F"];
          
          RObject the_chosen_one = as<RObject>(mat_list[(mat_chosen - 1)]);
          
          if (is<S4>(the_chosen_one)) {
            smat = as<arma::sp_mat>(the_chosen_one);
            sparse_mat = true;
          } else {
            mat = as<NumericMatrix>(mat_list[mat_chosen - 1]);
          }
        }
      } else {
        throw Rcpp::exception("Argument object is not of a usable class.", false);
      }
    }
    
  } else {
    throw Rcpp::exception("Argument object is not of a usable class.", false);
  }
  
  int mat_rows {0};
  int mat_cols {0};
  int smat_points {0};
  int mat_negs {0};
  int mat_poss {0};
  
  arma::umat smat_loc_mat;
  arma::vec smat_values;
  
  if (!sparse_mat) {
    mat_rows = mat.nrow();
    mat_cols = mat.ncol();
  } else {
    mat_rows = smat.n_rows;
    mat_cols = smat.n_cols;
  }
  
  if (!sparse_mat) {
    for (int i_row = 0; i_row < mat_rows; i_row++) {
      Rcpp::checkUserInterrupt();
      for (int j_col = 0; j_col < mat_cols; j_col++) {
        double test_value = static_cast<double>(mat(i_row, j_col));
        
        if (test_value > 0.0) {
          mat_poss++;
        } else if (test_value < 0.0) {
          mat_negs++;
        }
      }
    }
  } else {
    // Sparse matrix element counter
    // Based on code from post of Coatless Professor at
    // https://stackoverflow.com/questions/40222092/access-and-modify-the-non-zero-elements-of-sparse-matrix-of-class-armasp-mat-u
    arma::sp_mat::const_iterator smat_start = smat.begin();
    arma::sp_mat::const_iterator smat_end = smat.end();
    
    smat_points = std::distance(smat_start, smat_end);
    if (smat_points <= 0) {
      throw Rcpp::exception("Sparse matrix contains no non-zero elements.",
        false);
    }
    
    arma::umat smat_locs (2, smat_points);
    arma::uvec temp (2);
    arma::vec temp_values (smat_points);
    
    arma::sp_mat::const_iterator it = smat_start;
    
    for (int i = 0; i < smat_points; ++i) {
      temp(0) = it.row();
      temp(1) = it.col();
      smat_locs.col(i) = temp;
      
      double smat_entry = static_cast<double>(smat(it.row(), it.col()));
      temp_values(i) = smat_entry;
      
      if (smat_entry > 0.0) {
        mat_poss++;
      } else if (smat_entry < 0.0) {
        mat_negs++;
      }
      ++it;
    }
    
    smat_loc_mat = smat_locs;
    smat_values = temp_values;
  }
  
  int used_elem_nums {0};
  
  if (type == 1) {
    used_elem_nums = mat_poss;
  } else if (type == 2) {
    used_elem_nums = mat_negs;
  } else {
    used_elem_nums = mat_poss + mat_negs;
  }
  
  // Appropriate col & vec indices
  int hstages_cols {0};
  int agestages_cols {0};
  
  Rcpp::StringVector main_stage_vector;
  
  if (list_only) {
    hstages_cols = hstages.length();
    agestages_cols = agestages.length();
    
    if (agestages_cols > 1) {
      StringVector all_stages = agestages["stage"];
      StringVector all_ages = agestages["age"];
      int all_age_stages_num = all_stages.length();
      
      StringVector new_age_stages (all_age_stages_num);
      
      for (int i = 0; i < all_age_stages_num; i++) {
        String new_jack = all_stages(i);
        new_jack += "  age ";
        new_jack += all_ages(i);
        
        new_age_stages(i) = new_jack;
      }
      
      main_stage_vector = new_age_stages;
    } else if (hstages_cols > 1) {
      if (list_type == 1 || list_type == 4 || (list_type != 1 && part == 2)) {
        StringVector all_stage1 = hstages["stage_1"];
        StringVector all_stage2 = hstages["stage_2"];
        int all_stage_pairs_num = all_stage1.length();
        
        StringVector new_paired_stages (all_stage_pairs_num);
        
        for (int i = 0; i < all_stage_pairs_num; i++) {
          String new_jack = "st2: ";
          new_jack += all_stage2(i);
          new_jack += "  st1: ";
          new_jack += all_stage1(i);
          
          new_paired_stages(i) = new_jack;
        }
        
        main_stage_vector = new_paired_stages;
      } else {
        StringVector all_stages = stageframe["stage"];
        
        main_stage_vector = all_stages;
      }
    } else {
      StringVector all_stages = stageframe["stage"];
      
      main_stage_vector = all_stages;
    }
  }
  
  // Creation of major vectors
  Rcpp::StringVector found_from_stage (used_elem_nums);
  Rcpp::StringVector found_to_stage (used_elem_nums);
  Rcpp::NumericVector found_elems_current_values (used_elem_nums);
  Rcpp::NumericVector found_elems_test_values (used_elem_nums);
  Rcpp::IntegerVector found_elems_cols (used_elem_nums);
  Rcpp::IntegerVector found_elems_rows (used_elem_nums);
  Rcpp::IntegerVector found_elems_indices (used_elem_nums);
  
  if (!sparse_mat) {
    if (used_elem_nums > 0) {
      int sorted_elem_counter {0};
      for (int j_col = 0; j_col < mat_cols; j_col++) {
        Rcpp::checkUserInterrupt();
        for (int i_row = 0; i_row < mat_rows; i_row++) {
          double current_value = static_cast<double>(mat(i_row, j_col));
          double test_value = current_value;
          
          if (test_value < 0.0 && type != 1) test_value = test_value * -1.0;
          
          if ((test_value > 0.0 && type != 2) || (current_value < 0.0 && type == 2)) {
            if (sorted_elem_counter == 0) {
              found_elems_current_values(0) = current_value;
              found_elems_test_values(0) = test_value;
              
              if (list_only) {
                found_from_stage(0) = main_stage_vector(j_col);
                found_to_stage(0) = main_stage_vector(i_row);
              }
              
              found_elems_cols(0) = j_col + 1;
              found_elems_rows(0) = i_row + 1;
              found_elems_indices(0) = (j_col * mat_rows) + i_row + 1;
              
            } else {
              bool found_slot {false};
              int new_slot {0};
              
              for (int k = 0; k < used_elem_nums; k++) {
                if (test_value >= found_elems_test_values(k) && !found_slot) {
                  new_slot = k;
                  found_slot = true;
                }
              }
              
              if (sorted_elem_counter == 1 && !found_slot &&
                  sorted_elem_counter < used_elem_nums) {
                new_slot = 1;
                found_slot = true;
              }
              
              if (found_slot) {
                for (int k = (used_elem_nums - 2); k >= new_slot ; k--) {
                  found_elems_test_values(k+1) = found_elems_test_values(k);
                  found_elems_current_values(k+1) = found_elems_current_values(k);
                  
                  if (list_only) {
                    found_from_stage(k+1) = found_from_stage(k);
                    found_to_stage(k+1) = found_to_stage(k);
                  }
                  
                  found_elems_cols(k+1) = found_elems_cols(k);
                  found_elems_rows(k+1) = found_elems_rows(k);
                  found_elems_indices(k+1) = found_elems_indices(k);
                }
                
                found_elems_test_values(new_slot) = test_value;
                found_elems_current_values(new_slot) = current_value;
                
                if (list_only) {
                  found_from_stage(new_slot) = main_stage_vector(j_col);
                  found_to_stage(new_slot) = main_stage_vector(i_row);
                }
                
                found_elems_cols(new_slot) = j_col + 1;
                found_elems_rows(new_slot) = i_row + 1;
                found_elems_indices(new_slot) = (j_col * mat_rows) + i_row + 1;
              }
            }
            sorted_elem_counter++;
          }
        }
      }
    }
  } else {
    if (used_elem_nums > 0) {
      int sorted_elem_counter {0};
      for (int point_track = 0; point_track < smat_points; point_track++) {
        Rcpp::checkUserInterrupt();
        
        double current_value = static_cast<double>(smat_values(point_track));
        double test_value = current_value;
        
        if (test_value < 0.0 && type != 1) test_value = test_value * -1.0;
        
        if ((test_value > 0.0 && type != 2) || (current_value < 0.0 && type == 2)) {
          if (sorted_elem_counter == 0) {
            found_elems_current_values(0) = current_value;
            found_elems_test_values(0) = test_value;
            
            int current_row = smat_loc_mat(0, point_track);
            int current_col = smat_loc_mat(1, point_track);
            
            if (list_only) {
              found_from_stage(0) = main_stage_vector(current_col);
              found_to_stage(0) = main_stage_vector(current_row );
            }
            
            found_elems_cols(0) = current_col + 1;
            found_elems_rows(0) = current_row + 1;
            found_elems_indices(0) = (current_col * mat_rows) + current_row + 1;
            
          } else {
            bool found_slot {false};
            int new_slot {0};
            
            for (int k = 0; k < used_elem_nums; k++) {
              if (test_value >= found_elems_test_values(k) && !found_slot) {
                new_slot = k;
                found_slot = true;
              }
            }
            
            if (sorted_elem_counter == 1 && !found_slot &&
                sorted_elem_counter < used_elem_nums) {
              new_slot = 1;
              found_slot = true;
            }
            
            if (found_slot) {
              for (int k = (used_elem_nums - 2); k >= new_slot ; k--) {
                found_elems_test_values(k+1) = found_elems_test_values(k);
                found_elems_current_values(k+1) = found_elems_current_values(k);
                
                if (list_only) {
                  found_from_stage(k+1) = found_from_stage(k);
                  found_to_stage(k+1) = found_to_stage(k);
                }
                
                found_elems_cols(k+1) = found_elems_cols(k);
                found_elems_rows(k+1) = found_elems_rows(k);
                found_elems_indices(k+1) = found_elems_indices(k);
              }
              
              found_elems_test_values(new_slot) = test_value;
              found_elems_current_values(new_slot) = current_value;
              
              int current_row = smat_loc_mat(0, point_track);
              int current_col = smat_loc_mat(1, point_track);
              
              if (list_only) {
                found_from_stage(new_slot) = main_stage_vector(current_col);
                found_to_stage(new_slot) = main_stage_vector(current_row);
              }
              
              found_elems_cols(new_slot) = current_col + 1;
              found_elems_rows(new_slot) = current_row + 1;
              found_elems_indices(new_slot) = (current_col * mat_rows) + current_row + 1;
            }
          }
          sorted_elem_counter++;
        }
      }
    }
  }
  
  Rcpp::DataFrame out_table;
  
  if (list_only) {
    out_table = DataFrame::create(_["index"] = found_elems_indices,
      _["from_stage"] = found_from_stage, _["to_stage"] = found_to_stage,
      _["column"] = found_elems_cols, _["row"] = found_elems_rows,
      _["value"] = found_elems_current_values);
  } else {
    out_table = DataFrame::create(_["index"] = found_elems_indices,
      _["column"] = found_elems_cols, _["row"] = found_elems_rows,
      _["value"] = found_elems_current_values);
  }
  
  return(out_table);
}

//' Append Projections Into New lefkoProj Object
//' 
//' Function \code{append_lP()} combines two population projections. It takes
//' two \code{lefkoProj} objects and appends them into a new \code{lefkoProj}
//' object.
//' 
//' @name append_lP
//' 
//' @param proj1 A \code{lefkoProj} object.
//' @param proj2 A second \code{lefkoProj} object, based on the same stageframe
//' as \code{proj1}.
//' 
//' @return A list of class \code{lefkoProj}, which always includes the first
//' three elements of the following, and also includes the remaining elements
//' below when a \code{lefkoMat} object is used as input:
//' \item{projection}{A list of lists of matrices showing the total number of
//' individuals per stage per occasion. The first list corresponds to each
//' pop-patch followed by each population (this top-level list is a single
//' element in \code{f_projection3()}). The inner list corresponds to
//' replicates within each pop-patch or population.}
//' \item{stage_dist}{A list of lists of the actual stage distribution in each
//' occasion in each replicate in each pop-patch or population.}
//' \item{rep_value}{A list of lists of the actual reproductive value in each
//' occasion in each replicate in each pop-patch or population.}
//' \item{pop_size}{A list of matrices showing the total population size in each
//' occasion per replicate (row within data frame) per pop-patch or population
//' (list element). \code{NA} values will result if projections with different
//' numbers of time steps are appended.}
//' \item{labels}{A data frame showing the order of populations and patches in
//' item \code{projection}.}
//' \item{ahstages}{The original stageframe used in the study.}
//' \item{hstages}{A data frame showing the order of historical stage pairs.}
//' \item{agestages}{A data frame showing the order of age-stage pairs.}
//' \item{labels}{A short data frame indicating the population (always \code{1}),
//' and patch (either the numeric index of the single chosen patch, or \code{1}
//' in all other cases). Any pop-patches having the same designation across the
//' two input projections will be appended together.}
//' \item{control}{A data frame showing the number of replicates and time steps
//' corresponding to each set of projections, where each set corresponds to a
//' pop-patch within the labels object of each input projection.}
//' \item{density}{The data frame input under the density option. Only provided
//' if input by the user for at least one of the two projections. Output as a
//' nested list corresponding to each pop-patch - replicate.}
//' \item{density_vr}{The data frame input under the density_vr option. Only
//' provided if input by the user for at least one of the two projections.
//' Output as a nested list corresponding to each pop-patch - replicate.}
//' 
//' @section Notes:
//' \code{lefkoProj} objects resulting from previous appends can also be
//' appended.
//' 
//' @seealso \code{\link{projection3}()}
//' 
//' @examples
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
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.1, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3), stageframe = cypframe_raw,
//'   historical = FALSE)
//' 
//' cypmatrix2r_AB <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = c("A", "B"), stages = c("stage3", "stage2"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2",  patchcol = "patchid", indivcol = "individ")
//' 
//' cypmatrix2r_AC <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = c("A", "C"), stages = c("stage3", "stage2"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2",  patchcol = "patchid", indivcol = "individ")
//' 
//' cypproj1 <- projection3(cypmatrix2r_AB, nreps = 5, times = 15,
//'   stochastic = TRUE)
//' cypproj2 <- projection3(cypmatrix2r_AC, nreps = 10, times = 20,
//'   stochastic = TRUE)
//' cypproj3 <- append_lP(cypproj1, cypproj2)
//' 
//' @export append_lP
// [[Rcpp::export(append_lP)]]
Rcpp::List append_lP(Nullable<RObject> proj1 = R_NilValue,
  Nullable<RObject> proj2 = R_NilValue) {
  
  Rcpp::List proj1_list;
  Rcpp::List proj2_list;
  
  CharacterVector proj1_class;
  CharacterVector proj2_class;
  
  if (proj1.isNotNull()) {
    if (is<List>(proj1)) {
      proj1_list = as<List>(proj1);
      proj1_class = proj1_list.attr("class");
      
      bool lP_found {false};
      for (int i = 0; i < proj1_class.length(); i++) {
        if (stringcompare_simple(String(proj1_class(i)), "lefkoProj")) lP_found = true;
      }
      
      if (!lP_found) pop_error("proj1", "a lefkoProj object", "", 1);
    } else {
      pop_error("proj1", "a lefkoProj object", "", 1);
    }
  } else {
    pop_error("proj1", "a lefkoProj object", "", 1);
  }
  
  if (proj2.isNotNull()) {
    if (is<List>(proj2)) {
      proj2_list = as<List>(proj2);
      proj2_class = proj2_list.attr("class");
      
      bool lP_found {false};
      for (int i = 0; i < proj2_class.length(); i++) {
        if (stringcompare_simple(String(proj2_class(i)), "lefkoProj")) lP_found = true;
      }
      
      if (!lP_found) pop_error("proj2", "a lefkoProj object", "", 1);
    } else {
      pop_error("proj2", "a lefkoProj object", "", 1);
    }
  } else {
    pop_error("proj2", "a lefkoProj object", "", 1);
  }
  
  bool proj1_lP_check {false};
  bool proj2_lP_check {false};
  
  for (int i = 0; i < static_cast<int>(proj1_class.length()); i++) {
    if (stringcompare_simple(String(proj1_class(i)), "lefkoProj")) {
      proj1_lP_check = true;
    }
  }
  
  for (int i = 0; i < static_cast<int>(proj2_class.length()); i++) {
    if (stringcompare_simple(String(proj2_class(i)), "lefkoProj")) {
      proj2_lP_check = true;
    }
  }
  
  if (!proj1_lP_check || !proj2_lP_check) {
    throw Rcpp::exception("Only lefkoProj objects are allowed in arguments proj1 and proj2.", false);
  }
  
  CharacterVector proj1_element_names = proj1_list.names();
  List proj1_projection;
  List proj1_stage_dist;
  List proj1_rep_value;
  List proj1_pop_size;
  DataFrame proj1_labels;
  DataFrame proj1_ahstages;
  DataFrame proj1_hstages;
  DataFrame proj1_agestages;
  
  bool proj1_projection_check {false};
  bool proj1_stage_dist_check {false};
  bool proj1_rep_value_check {false};
  bool proj1_pop_size_check {false};
  bool proj1_labels_check {false};
  bool proj1_ahstages_check {false};
  bool proj1_hstages_check {false};
  bool proj1_agestages_check {false};
  bool proj1_control_check {false};
  bool proj1_density_check {false};
  bool proj1_density_vr_check {false};
  
  for (int i = 0; i < static_cast<int>(proj1_element_names.length()); i++) {
    if (stringcompare_simple(String(proj1_element_names(i)), "projection")) {
      proj1_projection_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "stage_dist")) {
      proj1_stage_dist_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "rep_value")) {
      proj1_rep_value_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "pop_size")) {
      proj1_pop_size_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "labels")) {
      proj1_labels_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "ahstages")) {
      proj1_ahstages_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "hstages")) {
      proj1_hstages_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "agestages")) {
      proj1_agestages_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "control")) {
      proj1_control_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "density")) {
      proj1_density_check = true;
    }
    if (stringcompare_simple(String(proj1_element_names(i)), "density_vr")) {
      proj1_density_vr_check = true;
    }
  }
  
  CharacterVector proj2_element_names = proj2_list.names();
  List proj2_projection;
  List proj2_stage_dist;
  List proj2_rep_value;
  List proj2_pop_size;
  DataFrame proj2_labels;
  
  bool proj2_projection_check {false};
  bool proj2_stage_dist_check {false};
  bool proj2_rep_value_check {false};
  bool proj2_pop_size_check {false};
  bool proj2_labels_check {false};
  bool proj2_ahstages_check {false};
  bool proj2_hstages_check {false};
  bool proj2_agestages_check {false};
  bool proj2_control_check {false};
  bool proj2_density_check {false};
  bool proj2_density_vr_check {false};
  
  for (int i = 0; i < static_cast<int>(proj2_element_names.length()); i++) {
    if (stringcompare_simple(String(proj2_element_names(i)), "projection")) {
      proj2_projection_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "stage_dist")) {
      proj2_stage_dist_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "rep_value")) {
      proj2_rep_value_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "pop_size")) {
      proj2_pop_size_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "labels")) {
      proj2_labels_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "ahstages")) {
      proj2_ahstages_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "hstages")) {
      proj2_hstages_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "agestages")) {
      proj2_agestages_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "control")) {
      proj2_control_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "density")) {
      proj2_density_check = true;
    }
    if (stringcompare_simple(String(proj2_element_names(i)), "density_vr")) {
      proj2_density_vr_check = true;
    }
  }
  
  if (!proj1_projection_check || !proj2_projection_check) {
    throw Rcpp::exception("Only lefkoProj objects are allowed in arguments proj1 and proj2.",
      false);
  }
  if (!proj1_rep_value_check || !proj2_rep_value_check) {
    throw Rcpp::exception("All lefkoProj objects must have rep_value elements.",
      false);
  }
  if (!proj1_stage_dist_check || !proj2_stage_dist_check) {
    throw Rcpp::exception("All lefkoProj objects must have stage_dist elements.",
      false);
  }
  if (!proj1_pop_size_check || !proj2_pop_size_check) {
    throw Rcpp::exception("All lefkoProj objects must have pop_size elements.",
      false);
  }
  if (!proj1_labels_check || !proj2_labels_check) {
    throw Rcpp::exception("All lefkoProj objects must have labels elements.",
      false);
  }
  if (!proj1_ahstages_check || !proj2_ahstages_check) {
    throw Rcpp::exception("All lefkoProj objects must have ahstages elements.",
      false);
  }
  if (!proj1_hstages_check || !proj2_hstages_check) {
    throw Rcpp::exception("All lefkoProj objects must have hstages elements.",
      false);
  }
  if (!proj1_agestages_check || !proj2_agestages_check) {
    throw Rcpp::exception("All lefkoProj objects must have agestages elements.",
      false);
  }
  if (!proj1_control_check || !proj2_control_check) {
    throw Rcpp::exception("All lefkoProj objects must have control elements.",
      false);
  }
  
  
  // Check ahstages, hstages, agestages
  proj1_ahstages = as<DataFrame>(proj1_list["ahstages"]);
  
  
  { // ahstages
    DataFrame proj2_ahstages = as<DataFrame>(proj2_list["ahstages"]);
    
    if (!LefkoUtils::df_compare(proj1_ahstages, proj2_ahstages)) {
      throw Rcpp::exception("Input projections cannot use different life history models.", 
        false);
    }
    
    // hstages
    RObject proj1_hstages_ro = as<RObject>(proj1_list["hstages"]);
    RObject proj2_hstages_ro = as<RObject>(proj2_list["hstages"]);
    
    if (is<DataFrame>(proj1_hstages_ro)) {
      if (!is<DataFrame>(proj2_hstages_ro)) {
        throw Rcpp::exception("One projection seems to be historical while the other does not.", 
          false);
      }
      
      proj1_hstages = as<DataFrame>(proj1_list["hstages"]);
      DataFrame proj2_hstages = as<DataFrame>(proj2_list["hstages"]);
      
      if (!LefkoUtils::df_compare(proj1_hstages, proj2_hstages)) {
        throw Rcpp::exception("Input projections cannot use different hstages objects.",
          false);
      }
    }
    
    // agestages
    RObject proj1_agestages_ro = as<RObject>(proj1_list["agestages"]);
    RObject proj2_agestages_ro = as<RObject>(proj2_list["agestages"]);
    
    if (is<DataFrame>(proj1_agestages_ro)) {
      if (!is<DataFrame>(proj2_agestages_ro)) {
        throw Rcpp::exception("One projection seems to be age-by-stage while the other does not.", 
          false);
      }
      
      proj1_agestages = as<DataFrame>(proj1_list["agestages"]);
      DataFrame proj2_agestages = as<DataFrame>(proj2_list["agestages"]);
      
      if (!LefkoUtils::df_compare(proj1_agestages, proj2_agestages)) {
        throw Rcpp::exception("Input projections cannot use different agestages objects.",
          false);
      }
    }
  }
  
  // Projection append routine
  proj1_projection = as<List>(proj1_list["projection"]);
  proj2_projection = as<List>(proj2_list["projection"]);
  
  IntegerVector proj1_control_int;
  IntegerVector proj2_control_int;
  IntegerMatrix proj1_control_mat;
  IntegerMatrix proj2_control_mat;
  IntegerVector proj1_control_index;
  IntegerVector proj2_control_index;
  
  bool proj1_control_int_log {false};
  bool proj2_control_int_log {false};
  IntegerVector proj1_reps;
  IntegerVector proj1_times;
  IntegerVector proj2_reps;
  IntegerVector proj2_times;
  {
    RObject proj1_control = as<RObject>(proj1_list["control"]);
    RObject proj2_control = as<RObject>(proj2_list["control"]);
    
    if (is<IntegerVector>(proj1_control) || is<NumericVector>(proj1_control)) {
      if (Rf_isMatrix(proj1_control)) { 
        proj1_control_mat = IntegerMatrix(proj1_control);
        
        proj1_control_index = proj1_control_mat(_, 0);
        proj1_reps = proj1_control_mat(_, 1);
        proj1_times = proj1_control_mat(_, 2);
      } else {
        proj1_control_int = IntegerVector(proj1_control);
        
        IntegerVector p1_listing (1);
        IntegerVector p1_times (1);
        IntegerVector p1_reps (1);
        
        p1_listing(0) = 1;
        p1_reps(0) = proj1_control_int(0);
        p1_times(0) = proj1_control_int(1);
        
        proj1_control_index = p1_listing;
        proj1_times = p1_times;
        proj1_reps = p1_reps;
        
        proj1_control_int_log = true;
      }
    }
    
    if (is<IntegerVector>(proj2_control) || is<NumericVector>(proj2_control)) {
      if (Rf_isMatrix(proj2_control)) {
        proj2_control_mat = IntegerMatrix(proj2_control);
        
        proj2_control_index = proj2_control_mat(_, 0);
        proj2_reps = proj2_control_mat(_, 1);
        proj2_times = proj2_control_mat(_, 2);
      } else {
        proj2_control_int = IntegerVector(proj2_control);
        
        IntegerVector p2_listing (1);
        IntegerVector p2_times (1);
        IntegerVector p2_reps (1);
        
        p2_listing(0) = 1;
        p2_reps(0) = proj2_control_int(0);
        p2_times(0) = proj2_control_int(1);
        
        proj2_control_index = p2_listing;
        proj2_times = p2_times;
        proj2_reps = p2_reps;
        
        proj2_control_int_log = true;
      }
    }
  }
  
  proj1_rep_value = as<List>(proj1_list["rep_value"]);
  proj2_rep_value = as<List>(proj2_list["rep_value"]);
  IntegerVector proj1_rep_value_lengths (static_cast<int>(proj1_projection.length()));
  IntegerVector proj2_rep_value_lengths (static_cast<int>(proj2_projection.length()));
  LogicalVector proj1_rep_value_list_check (static_cast<int>(proj1_projection.length()));
  LogicalVector proj2_rep_value_list_check (static_cast<int>(proj2_projection.length()));
  
  for (int i = 0; i < static_cast<int>(proj1_projection.length()); i++) {
    RObject p1_rv_i = as<RObject>(proj1_rep_value(i));
    if (is<List>(p1_rv_i)) {
      List p1_rv_i_list = as<List>(p1_rv_i);
      
      proj1_rep_value_list_check(i) = true;
      
      int p1_rv_i_list_length = static_cast<int>(p1_rv_i_list.length());
      proj1_rep_value_lengths(i) = p1_rv_i_list_length;
    } else if (is<NumericVector>(p1_rv_i)) {
      proj1_rep_value_list_check(i) = false;
      proj1_rep_value_lengths(i) = 1;
    }
  }
  for (int i = 0; i < static_cast<int>(proj2_projection.length()); i++) {
    RObject p2_rv_i = as<RObject>(proj2_rep_value(i));
    if (is<List>(p2_rv_i)) {
      List p2_rv_i_list = as<List>(p2_rv_i);
      
      proj2_rep_value_list_check(i) = true;
      
      int p2_rv_i_list_length = static_cast<int>(p2_rv_i_list.length());
      proj2_rep_value_lengths(i) = p2_rv_i_list_length;
      
    } else if (is<NumericVector>(p2_rv_i)) {
      proj2_rep_value_list_check(i) = false;
      proj2_rep_value_lengths(i) = 1;
    }
  }
  
  proj1_stage_dist = as<List>(proj1_list["stage_dist"]);
  proj2_stage_dist = as<List>(proj2_list["stage_dist"]);
  IntegerVector proj1_stage_dist_lengths (static_cast<int>(proj1_projection.length()));
  IntegerVector proj2_stage_dist_lengths (static_cast<int>(proj2_projection.length()));
  LogicalVector proj1_stage_dist_list_check (static_cast<int>(proj1_projection.length()));
  LogicalVector proj2_stage_dist_list_check (static_cast<int>(proj2_projection.length()));
  
  for (int i = 0; i < static_cast<int>(proj1_projection.length()); i++) {
    RObject p1_sd_i = as<RObject>(proj1_stage_dist(i));
    if (is<List>(p1_sd_i)) {
      List p1_sd_i_list = as<List>(p1_sd_i);
      
      proj1_stage_dist_list_check(i) = true;
      
      int p1_sd_i_list_length = static_cast<int>(p1_sd_i_list.length());
      proj1_stage_dist_lengths(i) = p1_sd_i_list_length;
    } else if (is<NumericVector>(p1_sd_i)) {
      proj1_stage_dist_list_check(i) = false;
      proj1_stage_dist_lengths(i) = 1;
    }
  }
  for (int i = 0; i < static_cast<int>(proj2_projection.length()); i++) {
    RObject p2_sd_i = as<RObject>(proj2_stage_dist(i));
    if (is<List>(p2_sd_i)) {
      List p2_sd_i_list = as<List>(p2_sd_i);
      
      proj2_stage_dist_list_check(i) = true;
      
      int p2_sd_i_list_length = static_cast<int>(p2_sd_i_list.length());
      proj2_stage_dist_lengths(i) = p2_sd_i_list_length;
      
    } else if (is<NumericVector>(p2_sd_i)) {
      proj2_stage_dist_list_check(i) = false;
      proj2_stage_dist_lengths(i) = 1;
    }
  }
  
  IntegerVector pop1_size_cols (static_cast<int>(proj1_projection.length()));
  IntegerVector pop2_size_cols (static_cast<int>(proj2_projection.length()));
  IntegerVector pop1_size_rows (static_cast<int>(proj1_projection.length()));
  IntegerVector pop2_size_rows (static_cast<int>(proj2_projection.length()));
  proj1_pop_size = as<List>(proj1_list["pop_size"]);
  proj2_pop_size = as<List>(proj2_list["pop_size"]);
  
  for (int i = 0; i < proj1_projection.length(); i++) {
    NumericMatrix pop_size_mat = as<NumericMatrix>(proj1_pop_size(i));
    int found_pop_size_cols = pop_size_mat.ncol();
    int found_pop_size_rows = pop_size_mat.nrow();
    
    pop1_size_cols(i) = found_pop_size_cols;
    pop1_size_rows(i) = found_pop_size_rows;
  }
  for (int i = 0; i < proj2_projection.length(); i++) {
    NumericMatrix pop_size_mat = as<NumericMatrix>(proj2_pop_size(i));
    unsigned int found_pop_size_cols = pop_size_mat.ncol();
    unsigned int found_pop_size_rows = pop_size_mat.nrow();
    
    pop2_size_cols(i) = found_pop_size_cols;
    pop2_size_rows(i) = found_pop_size_rows;
  }
  
  // Tests of density and density_vr elements
  int proj1_density_type {0}; // 0 = NULL, 1 = DataFrame, 2 = List
  int proj2_density_type {0}; // 0 = NULL, 1 = DataFrame, 2 = List
  int proj1_density_vr_type {0}; // 0 = NULL, 1 = DataFrame, 2 = List
  int proj2_density_vr_type {0}; // 0 = NULL, 1 = DataFrame, 2 = List
  
  if (proj1_density_check) {
    RObject proj1_density_ro = as<RObject>(proj1_list["density"]);
    
    if (is<DataFrame>(proj1_density_ro)) {
      proj1_density_type = 1;
    } else if (is<List>(proj1_density_ro)) {
      proj1_density_type = 2;
      List proj1_density_li = as<List>(proj1_list["density"]);
    }
  }
  if (proj2_density_check) {
    RObject proj2_density_ro = as<RObject>(proj2_list["density"]);
    
    if (is<DataFrame>(proj2_density_ro)) {
      proj2_density_type = 1;
    } else if (is<List>(proj2_density_ro)) {
      proj2_density_type = 2;
      List proj2_density_li = as<List>(proj2_list["density"]);
    }
  }
  
  if (proj1_density_vr_check) {
    RObject proj1_density_vr_ro = as<RObject>(proj1_list["density_vr"]);
    
    if (is<DataFrame>(proj1_density_vr_ro)) {
      proj1_density_vr_type = 1;
    } else if (is<List>(proj1_density_vr_ro)) {
      proj1_density_vr_type = 2;
      List proj1_density_vr_li = as<List>(proj1_list["density_vr"]);
    }
  }
  if (proj2_density_vr_check) {
    RObject proj2_density_vr_ro = as<RObject>(proj2_list["density_vr"]);
    
    if (is<DataFrame>(proj2_density_vr_ro)) {
      proj2_density_vr_type = 1;
    } else if (is<List>(proj2_density_vr_ro)) {
      proj2_density_vr_type = 2;
      List proj2_density_vr_li = as<List>(proj2_list["density_vr"]);
    }
  }
  
  proj1_labels = as<DataFrame>(proj1_list["labels"]);
  StringVector proj1_labels_pop = as<StringVector>(proj1_labels["pop"]);
  StringVector proj1_labels_patch = as<StringVector>(proj1_labels["patch"]);
  
  proj2_labels = as<DataFrame>(proj2_list["labels"]);
  StringVector proj2_labels_pop = as<StringVector>(proj2_labels["pop"]);
  StringVector proj2_labels_patch = as<StringVector>(proj2_labels["patch"]);
  
  int base_length_projection = static_cast<int>(proj1_labels_pop.length());
  int potential_to_add = static_cast<int>(proj2_labels_pop.length());
  
  IntegerVector proj2_labels_equivalence = rep(-1, potential_to_add);
  IntegerVector proj2_labels_allassigned = rep(-1, potential_to_add);
  
  for (int i = 0; i < static_cast<int>(proj2_labels_pop.length()); i++) {
    for (int j = 0; j < static_cast<int>(proj1_labels_pop.length()); j++) {
      if (stringcompare_simple(String(proj2_labels_pop(i)), String(proj1_labels_pop(j)))) {
        if (stringcompare_simple(String(proj2_labels_patch(i)), String(proj1_labels_patch(j)))) {
          proj2_labels_equivalence(i) = j;
          proj2_labels_allassigned(i) = j;
        }
      }
    }
  }
  
  int added_length {0};
  int max_found_equivalence = max(proj2_labels_equivalence);
  int found_in_proj1 {0};
  
  for (int i = 0; i < potential_to_add; i++) {
    if (proj2_labels_equivalence(i) == -1) {
      added_length++;
      proj2_labels_allassigned(i) = max_found_equivalence + added_length;
    } else {
      found_in_proj1++;
    }
  }
  
  IntegerVector proj1_control_breaks;
  IntegerVector proj2_control_breaks;
  IntegerVector proj1_control_lengths;
  IntegerVector proj2_control_lengths;
  
  if (proj1_control_int_log) {
    IntegerVector p1_control_breaks (proj1_labels_pop.length());
    IntegerVector p1_control_lengths (proj1_labels_pop.length());
    p1_control_breaks.fill(0);
    p1_control_lengths.fill(1);
    
    proj1_control_breaks = p1_control_breaks;
    proj1_control_lengths = p1_control_lengths;
  } else {
    arma::ivec proj1_control_index_arma = as<arma::ivec>(proj1_control_index);
    
    IntegerVector p1_control_breaks (static_cast<int>(max(proj1_control_index)));
    IntegerVector p1_control_lengths (static_cast<int>(max(proj1_control_index)));
    
    for (int i = 0; i < static_cast<int>(max(proj1_control_index)); i++) {
      arma::uvec found_indices = find(proj1_control_index_arma == (i + 1));
      p1_control_breaks(i) = min(found_indices);
      p1_control_lengths(i) = found_indices.n_elem;
    }
    
    proj1_control_breaks = p1_control_breaks;
    proj1_control_lengths = p1_control_lengths;
  }
  
  if (proj2_control_int_log) {
    IntegerVector p2_control_breaks (proj2_labels_pop.length());
    IntegerVector p2_control_lengths (proj2_labels_pop.length());
    p2_control_breaks.fill(0);
    p2_control_lengths.fill(1);
    
    proj2_control_breaks = p2_control_breaks;
    proj2_control_lengths = p2_control_lengths;
  } else {
    arma::ivec proj2_control_index_arma = as<arma::ivec>(proj2_control_index);
    
    IntegerVector p2_control_breaks (static_cast<int>(max(proj2_control_index)));
    IntegerVector p2_control_lengths (static_cast<int>(max(proj2_control_index)));
    
    for (int i = 0; i < static_cast<int>(max(proj2_control_index)); i++) {
      arma::uvec found_indices = find(proj2_control_index_arma == (i + 1));
      p2_control_breaks(i) = min(found_indices);
      p2_control_lengths(i) = found_indices.n_elem;
    }
    
    proj2_control_breaks = p2_control_breaks;
    proj2_control_lengths = p2_control_lengths;
  }
  
  int new_length = base_length_projection + added_length;
  List new_proj_projection (new_length);
  List new_pop_size (new_length);
  List new_rep_value (new_length);
  List new_stage_dist (new_length);
  List new_density (new_length);
  List new_density_vr (new_length);
  
  int new_control_rows {0};
  
  if (proj1_control_int_log && proj2_control_int_log) { 
    new_control_rows = base_length_projection + found_in_proj1 + added_length;
  } else if (proj1_control_int_log) {
    new_control_rows = base_length_projection + 
      static_cast<int>(proj2_control_index.length());
  } else if (proj2_control_int_log) {
    new_control_rows = static_cast<int>(proj1_control_index.length()) + 
      added_length;
  } else {
    new_control_rows = static_cast<int>(proj1_control_index.length()) + 
      static_cast<int>(proj2_control_index.length());
  }
  
  IntegerMatrix new_control (new_control_rows, 3);
  int control_row_tracker {0};
  
  StringVector new_labels_pop (base_length_projection + added_length);
  StringVector new_labels_patch (base_length_projection + added_length);
  
  for (int i = 0; i < new_length; i++) {
    List new_proj_projection_i;
    
    if (i < base_length_projection) {
      List proj1_projection_i = as<List>(proj1_projection(i));
      int p1pi_length = static_cast<int>(proj1_projection_i.length());
      bool found_check {false};
      
      for (int j = 0; j < static_cast<int>(proj1_control_lengths(i)); j++) {
        new_control(control_row_tracker, 0) = i + 1;
        new_control(control_row_tracker, 1) = proj1_reps(proj1_control_breaks(i) + j);
        new_control(control_row_tracker, 2) = proj1_times(proj1_control_breaks(i) + j);
        control_row_tracker++;
      }
      
      new_labels_pop(i) = proj1_labels_pop(i);
      new_labels_patch(i) = proj1_labels_patch(i);
      
      for (int j = 0; j < potential_to_add; j++) {
        if (proj2_labels_equivalence(j) == i) {
          found_check = true;
          
          List proj2_projection_j = as<List>(proj2_projection(j));
          int p2pj_length = static_cast<int>(proj2_projection_j.length());
          
          int total_new_projection_length = p1pi_length + p2pj_length;
          
          List temp_new (total_new_projection_length);
          
          for (int k = 0; k < p1pi_length; k++) {
            temp_new(k) = proj1_projection_i(k);
          }
          for (int k = 0; k < p2pj_length; k++) {
            temp_new(k + p1pi_length) = proj2_projection_j(k);
          }
          new_proj_projection_i = temp_new;
          
          if (static_cast<bool>(proj1_rep_value_list_check[i]) &&
            static_cast<bool>(proj2_rep_value_list_check[j])) {
            
            List proj1_rep_value_i = as<List>(proj1_rep_value(i));
            List proj2_rep_value_j = as<List>(proj2_rep_value(j));
            
            List temp_new_rv (total_new_projection_length);
            for (int k = 0; k < p1pi_length; k++) {
              temp_new_rv(k) = proj1_rep_value_i(k);
            }
            for (int k = 0; k < p2pj_length; k++) {
              temp_new_rv(k + p1pi_length) = proj2_rep_value_j(k);
            }
            new_rep_value(i) = temp_new_rv;
          } else {
            NumericVector new_0 = {0.};
            
            new_rep_value(i) = new_0;
          }
          
          if (static_cast<bool>(proj1_stage_dist_list_check[i]) &&
            static_cast<bool>(proj2_stage_dist_list_check[j])) {
            
            List proj1_stage_dist_i = as<List>(proj1_stage_dist(i));
            List proj2_stage_dist_j = as<List>(proj2_stage_dist(j));
            
            List temp_new_sd (total_new_projection_length);
            for (int k = 0; k < p1pi_length; k++) {
              temp_new_sd(k) = proj1_stage_dist_i(k);
            }
            for (int k = 0; k < p2pj_length; k++) {
              temp_new_sd(k + p1pi_length) = proj2_stage_dist_j(k);
            }
            new_stage_dist(i) = temp_new_sd;
          } else {
            NumericVector new_0 = {0.};
            
            new_stage_dist(i) = new_0;
          }
          
          // density
          if (proj1_density_check) { 
            if (proj1_density_type == 1) {
              DataFrame proj1_density_temp = as<DataFrame>(proj1_list["density"]);
              
              List temp_new_de (total_new_projection_length);
              for (int k = 0; k < p1pi_length; k++) {
                temp_new_de(k) = proj1_density_temp;
              }
              
              if (proj2_density_type == 1) { 
                DataFrame proj2_density_temp = as<DataFrame>(proj2_list["density"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_de(k + p1pi_length) = proj2_density_temp;
                }
              } else if (proj2_density_type == 2) {
                List proj2_density_temp = as<List>(proj2_list["density"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_de(k + p1pi_length) = proj2_density_temp(j);
                }
              } else {
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_de(k + p1pi_length) = R_NilValue;
                }
              }
              
              new_density(i) = temp_new_de;
            } else if (proj1_density_type == 2) {
              List proj1_density_temp = as<List>(proj1_list["density"]);
              
              List temp_new_de (total_new_projection_length);
              for (int k = 0; k < p1pi_length; k++) {
                temp_new_de(k) = proj1_density_temp(i);
              }
              
              if (proj2_density_type == 1) { 
                DataFrame proj2_density_temp = as<DataFrame>(proj2_list["density"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_de(k + p1pi_length) = proj2_density_temp;
                }
              } else if (proj2_density_type == 2) {
                List proj2_density_temp = as<List>(proj2_list["density"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_de(k + p1pi_length) = proj2_density_temp(j);
                }
              } else {
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_de(k + p1pi_length) = R_NilValue;
                }
              }
              
              new_density(i) = temp_new_de;
            }
          } else if (proj2_density_check) { 
            List temp_new_de (total_new_projection_length);
            for (int k = 0; k < p1pi_length; k++) {
              temp_new_de(k) = R_NilValue;
            }
            
            if (proj2_density_type == 1) { 
              DataFrame proj2_density_temp = as<DataFrame>(proj2_list["density"]);
            
              for (int k = 0; k < p2pj_length; k++) {
                temp_new_de(k + p1pi_length) = proj2_density_temp;
              }
            } else if (proj2_density_type == 2) {
              List proj2_density_temp = as<List>(proj2_list["density"]);
            
              for (int k = 0; k < p2pj_length; k++) {
                temp_new_de(k + p1pi_length) = proj2_density_temp(j);
              }
            } else {
              for (int k = 0; k < p2pj_length; k++) {
                temp_new_de(k + p1pi_length) = R_NilValue;
              }
            }
            
            new_density(i) = temp_new_de;
          }
          
          // density_vr
          if (proj1_density_vr_check) { 
            if (proj1_density_vr_type == 1) {
              DataFrame proj1_density_vr_temp = as<DataFrame>(proj1_list["density_vr"]);
              
              List temp_new_devr (total_new_projection_length);
              for (int k = 0; k < p1pi_length; k++) {
                temp_new_devr(k) = proj1_density_vr_temp;
              }
              
              if (proj2_density_vr_type == 1) { 
                DataFrame proj2_density_vr_temp = as<DataFrame>(proj2_list["density_vr"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_devr(k + p1pi_length) = proj2_density_vr_temp;
                }
              } else if (proj2_density_vr_type == 2) {
                List proj2_density_vr_temp = as<List>(proj2_list["density_vr"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_devr(k + p1pi_length) = proj2_density_vr_temp(j);
                }
              } else {
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_devr(k + p1pi_length) = R_NilValue;
                }
              }
              
              new_density_vr(i) = temp_new_devr;
            } else if (proj1_density_vr_type == 2) {
              List proj1_density_vr_temp = as<List>(proj1_list["density_vr"]);
              
              List temp_new_devr (total_new_projection_length);
              for (int k = 0; k < p1pi_length; k++) {
                temp_new_devr(k) = proj1_density_vr_temp(i);
              }
              
              if (proj2_density_vr_type == 1) { 
                DataFrame proj2_density_vr_temp = as<DataFrame>(proj2_list["density_vr"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_devr(k + p1pi_length) = proj2_density_vr_temp;
                }
              } else if (proj2_density_vr_type == 2) {
                List proj2_density_vr_temp = as<List>(proj2_list["density_vr"]);
              
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_devr(k + p1pi_length) = proj2_density_vr_temp(j);
                }
              } else {
                for (int k = 0; k < p2pj_length; k++) {
                  temp_new_devr(k + p1pi_length) = R_NilValue;
                }
              }
              
              new_density_vr(i) = temp_new_devr;
            }
          } else if (proj2_density_vr_check) { 
            List temp_new_devr (total_new_projection_length);
            for (int k = 0; k < p1pi_length; k++) {
              temp_new_devr(k) = R_NilValue;
            }
            
            if (proj2_density_vr_type == 1) { 
              DataFrame proj2_density_vr_temp = as<DataFrame>(proj2_list["density_vr"]);
            
              for (int k = 0; k < p2pj_length; k++) {
                temp_new_devr(k + p1pi_length) = proj2_density_vr_temp;
              }
            } else if (proj2_density_vr_type == 2) {
              List proj2_density_vr_temp = as<List>(proj2_list["density_vr"]);
            
              for (int k = 0; k < p2pj_length; k++) {
                temp_new_devr(k + p1pi_length) = proj2_density_vr_temp(j);
              }
            } else {
              for (int k = 0; k < p2pj_length; k++) {
                temp_new_devr(k + p1pi_length) = R_NilValue;
              }
            }
            
            new_density_vr(i) = temp_new_devr;
          }
          
          for (int l = 0; l < static_cast<int>(proj2_control_lengths(j)); l++) {
            new_control(control_row_tracker, 0) = i + 1;
            new_control(control_row_tracker, 1) = proj2_reps(proj2_control_breaks(j) + l);
            new_control(control_row_tracker, 2) = proj2_times(proj2_control_breaks(j) + l);
            control_row_tracker++;
          }
          
          // Merge pop_size matrices
          int new_pop_size_mat_rows = static_cast<int>(pop1_size_rows(i)) + static_cast<int>(pop2_size_rows(j));
          int new_pop_size_mat_cols = static_cast<int>(pop1_size_cols(i));
          if (pop2_size_cols(j) > pop1_size_cols(i)) new_pop_size_mat_cols = 
            static_cast<int>(pop2_size_cols(j));
          
          NumericMatrix old_pop1_size_mat = as<NumericMatrix>(proj1_pop_size[i]);
          NumericMatrix old_pop2_size_mat = as<NumericMatrix>(proj2_pop_size[j]);
          NumericMatrix new_pop_size_mat (new_pop_size_mat_rows, new_pop_size_mat_cols);
          new_pop_size_mat.fill(NA_REAL);
          
          for (int k = 0; k < pop1_size_rows(i); k++) {
            for (int l = 0; l < pop1_size_cols(i); l++) {
              new_pop_size_mat(k, l) = old_pop1_size_mat(k, l);
            }
          }
          for (int k = 0; k < pop2_size_rows(j); k++) {
            for (int l = 0; l < pop2_size_cols(j); l++) {
              new_pop_size_mat((k + pop1_size_rows(j)), l) = old_pop2_size_mat(k, l);
            }
          }
          
          new_pop_size(i) = new_pop_size_mat;
          
          break;
        }
      }
      
      if (!found_check) {
        new_proj_projection_i = proj1_projection_i;
        new_pop_size(i) = proj1_pop_size(i);
        
        if (static_cast<bool>(proj1_rep_value_list_check(i))) {
          
          List proj1_rep_value_i = as<List>(proj1_rep_value[i]);
          new_rep_value(i) = proj1_rep_value_i;
        } else {
          NumericVector new_0 = {0.};
          new_rep_value(i) = new_0;
        }
        
        if (static_cast<bool>(proj1_stage_dist_list_check(i))) {
          
          List proj1_stage_dist_i = as<List>(proj1_stage_dist[i]);
          new_stage_dist(i) = proj1_stage_dist_i;
        } else {
          NumericVector new_0 = {0.};
          new_stage_dist(i) = new_0;
        }
        
        if (proj1_density_check) {
          if (proj1_density_type == 1) {
            DataFrame proj1_density_temp = as<DataFrame>(proj1_list["density"]);
            new_density(i) = proj1_density_temp;
          } else if (proj1_density_type == 2) {
            List proj1_density_temp = as<List>(proj1_list["density"]);
            new_density(i) = proj1_density_temp(i);
          }
        } else if (proj2_density_check) {
          new_density(i) = R_NilValue;
        }
      }
    } else {
      for (int j = 0; j < potential_to_add; j++) {
        if (proj2_labels_allassigned(j) == i) {
          List proj2_projection_i = as<List>(proj2_projection(j));
          new_proj_projection_i = proj2_projection_i;
          
          if (static_cast<bool>(proj2_rep_value_list_check(j))) {
            
            List proj2_rep_value_i = as<List>(proj2_rep_value[j]);
            new_rep_value(i) = proj2_rep_value_i;
          } else {
            NumericVector new_0 = {0.};
            new_rep_value(i) = new_0;
          }
          
          if (static_cast<bool>(proj2_stage_dist_list_check(j))) {
            
            List proj2_stage_dist_i = as<List>(proj2_stage_dist[j]);
            new_stage_dist(i) = proj2_stage_dist_i;
          } else {
            NumericVector new_0 = {0.};
            new_stage_dist(i) = new_0;
          }
          
          if (proj2_density_check) {
            if (proj2_density_type == 1) {
              DataFrame proj2_density_temp = as<DataFrame>(proj2_list["density"]);
              new_density(i) = proj2_density_temp;
            } else if (proj2_density_type == 2) {
              List proj2_density_temp = as<List>(proj2_list["density"]);
              new_density(i) = proj2_density_temp(j);
            }
          } else if (proj1_density_check) {
            new_density(i) = R_NilValue;
          }
          
          for (int l = 0; l < static_cast<int>(proj2_control_lengths(j)); l++) {
            new_control(control_row_tracker, 0) = base_length_projection + j;
            new_control(control_row_tracker, 1) = proj2_reps(proj2_control_breaks(j) + l);
            new_control(control_row_tracker, 2) = proj2_times(proj2_control_breaks(j) + l);
            control_row_tracker++;
          }
          
          int new_label_position {-1};
          for (int k = 0; k < proj2_labels_allassigned.length(); k++) {
            if (proj2_labels_allassigned(k) == i) new_label_position = k;
          }
          
          new_labels_pop(i) = proj2_labels_pop(new_label_position);
          new_labels_patch(i) = proj2_labels_patch(new_label_position);
          
          new_pop_size(i) = proj2_pop_size(new_label_position);
          
          break;
        }
      }
    }
    
    new_proj_projection(i) = new_proj_projection_i;
  }
  
  DataFrame new_labels = DataFrame::create(_["pop"] = new_labels_pop,
    _["patch"] = new_labels_patch);
  
  List new_output (11);
  new_output(0) = new_proj_projection;
  new_output(1) = new_stage_dist;
  new_output(2) = new_rep_value;
  new_output(3) = new_pop_size;
  new_output(4) = new_labels;
  new_output(5) = proj1_ahstages;
  new_output(6) = proj1_hstages;
  new_output(7) = proj1_agestages;
  new_output(8) = new_control;
  new_output(9) = new_density;
  new_output(10) = new_density_vr;
  
  CharacterVector new_output_names = {"projection", "stage_dist", "rep_value",
    "pop_size", "labels", "ahstages", "hstages", "agestages", "control",
    "density", "density_vr"};
  new_output.attr("names") = new_output_names;
  new_output.attr("class") = "lefkoProj";
  
  return(new_output);
}


