// [[Rcpp::depends(RcppArmadillo)]]
#define BOOST_DISABLE_ASSERTS

#include <RcppArmadilloExtensions/sample.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;



// Index of functions
// 
// 1. specialpatrolgroup() - Calculates matrix transitions in raw historical matrices
// 2. normalpatrolgroup() - Calculates matrix transitions in raw ahistorical martrices
// 3. minorpatrolgroup() - Calculates matrix transitions in raw Leslie MPMs
// 4. subvertedpatrolgroup() - Calculates matrix transitions in raw age-by-stage matrices
// 5. raymccooney() - Takes various vital rate models and other parameters and coordinates them as input into the function-based matrix estimation functions
// 6. mothermccooney() - Takes various vital rate models and other parameters and coordinates them as input into function fleslie()
// 7. f_projection3() - Develops and projects function-based matrix models
// 8. mpm_create() - Core workhorse function that creates all flavors of MPM in lefko3
// 9. .ss3matrix() - Returns stable stage distribution for a dense or sparse matrix
// 10. .ssmatrix_sp() - Returns stable stage distribution for a sparse matrix
// 11. .rv3matrix() - Returns reproductive values in a dense or sparse matrix
// 12. .rv3matrix_sp() - Returns reproductive values in a sparse matrix
// 13. .sens3matrix() - Returns sensitivity of lambda to each element of sparse or dense matrix
// 14. .sens3sp_matrix() - Returns sensitivity of lambda to each element in a sparse matrix, with output in sparse format
// 15. .sens3matrix_spinp() - Returns sensitivity of lambda to each element in a sparse matrix, with output in dense matrix format
// 16. .sens3hlefko() - Returns sensitivity of lambda to each historical stage-pair, and associated life stage
// 17. .sens3hlefko_sp() - Returns sensitivity of lambda to each historical stage-pair, and associated life stage, with input in sparse format
// 18. .elas3matrix() - Returns elasticity of lambda to each element in dense or sparse matrix
// 19. .ekas3sp_matrix() - Returns elasticity of lambda to each element in sparse matrix, in sparse output
// 20. .elas3hlefko() - Returns elasticity of lambda to each historical stage-pair, and each associated life stage
// 21. .elas3sp_hlefko() - Returns elasticity of lambda to each historical stage-pair, and each associated life stage, with sparse input
// 22. .proj3() - Cure functiuon running matrix projections used in other functions in lefko3
// 23. .proj3sp() - Core function running sparse matrix projections used in other functions in lefko3
// 24. .proj3dens() - Core function running density-dependent projections used in other functions in lefko3
// 25. projection3() - Runs projection simulations with lefkoMat objects
// 26. slambda3() - Estimates stochastic population growth rate in lefkoMat objects and other MPMs
// 27. .stoch_senselas() - Estimates sensitivity and elasticity of matrix elements to a
// 28. .ltre3matrix() - Returns one-way fixed deterministic LTRE matrix
// 29. .sltre3matrix() - Returns one-way stochastic LTRE matrices
// 30. .snaltre3matrix() - Returns one-way small noise approximation LTRE matrices
// 31. markov_run() - Creates vector of randomly sampled times

// Matrix Estimation

//' Estimate All Elements of Raw Historical Matrix
//' 
//' Function \code{specialpatrolgroup()} swiftly calculates matrix transitions
//' in raw historical matrices, and serves as the core workhorse function behind
//' \code{rlefko3()}.
//' 
//' @name specialpatrolgroup
//' 
//' @param sge9l The Allstages data frame developed for \code{rlefko3()}
//' covering stage pairs across times \emph{t}+1, \emph{t} and \emph{t}-1.
//' Generally termed \code{stageexpansion9}.
//' @param sge3index21 Integer index vector of stages in times \emph{t}-1 and
//' \emph{t}, from \code{stageexpansion3}.
//' @param MainData The demographic dataset modified to hold \code{usedfec}
//' columns.
//' @param StageFrame The full stageframe for the analysis.
//' @param format Indicates whether to output Ehrlen-format hMPMs (\code{1}) or
//' deVries-format hMPMs (\code{2}).
//' @param err_switch A logical value. If set to \code{TRUE}, then will also
//' output \code{probsrates} and \code{stage2fec}. Defaults to \code{FALSE}.
//' @param loypop A string vector giving the order of populations in the list of
//' years.
//' @param loypatch A string vector giving the order of patches in the list of
//' years.
//' @param loyyear2 A string vector giving the order of years at time \emph{t}
//' in the list of years.
//' @param yearorder The integer year order corresponding to \code{loyyear2}.
//' @param pop_var_int The variable number coding for population in the main
//' data set.
//' @param patch_var_int The variable number coding for patch in the main data
//' set.
//' @param year_var_int The variable number coding for year in time \emph{t} in
//' the main data set.
//' @param loy_pop_used A logical value indicating whether the population
//' variable is to be used.
//' @param loy_patch_used A logical value indicating whether the patch variable
//' is to be used.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' @param sparse If \code{TRUE}, then output will be in sparse matrix format.
//' Defaults to \code{FALSE}.
//' 
//' @return If \code{err_switch = FALSE}, then will return a list of three
//' matrices, including the survival-transition (\code{U}) matrix, the fecundity
//' matrix (\code{F}), and the sum (\code{A}) matrix, with the \code{A} matrix
//' first, followed by the \code{ahstages}, \code{agestages}, \code{hstages},
//' and \code{labels} data frames, and the \code{matrixqc} and \code{dataqc}
//' vectors. If \code{err_switch = TRUE}, then will also output four further
//' elements: \code{probsrates_all} and \code{stage2fec_all}. The former is a
//' matrix composed of the following vectors in order: \code{sge9index321},
//' which gives the historical index number for each transition possible and
//' in order (from \code{stageexpansion9}); \code{aliveandequal}, which gives
//' the element index in the matrix associated with that transition (or
//' \code{-1} if the transition is impossible); \code{probsrates0}, which gives
//' the total number of individuals counted for a particular historical
//' transition; \code{probsrates1}, which gives the total number of individuals
//' associated with a particular paired stage in times \emph{t}-1 and \emph{t};
//' \code{probsrates2}, which gives the total number of individuals for each
//' paired stage in times \emph{t}-1 and \emph{t} that survive into time
//' \emph{t}+1; \code{probsrates3}, which gives the total fecundity found for
//' the historical transition; \code{sge9fec32}, which gives the predicted
//' reproductive status of the historical transition as given in the
//' \code{repmatrix}; \code{probsrates0p}, which is \code{probsrates0} but
//' assuming a prior stage; \code{probsrates1p}, which is \code{probsrates1} but
//' assuming a prior stage; \code{probsrates2p}, which is \code{probsrates2} but
//' assuming a prior stage; and \code{probsrates3p}, which is \code{probsrates3}
//' but assuming a prior stage. Element \code{stage2fec_all} is a matrix
//' composed of two vectors: \code{stage2fec} is the fecundity associated with
//' each paired historical stage, and \code{stage2fecp} is the equivalent
//' assuming a prior stage. The final two elements are \code{dataindex321_prior}
//' and the edited dataset.
//' 
//' @keywords internal
//' @noRd
List specialpatrolgroup(const DataFrame& sge9l, const arma::ivec& sge3index21,
  const DataFrame& MainData, const DataFrame& StageFrame, int format,
  int err_switch, StringVector loypop, StringVector loypatch,
  StringVector loyyear2, IntegerVector yearorder, const int pop_var_int,
  const int patch_var_int, const int year_var_int, const bool loy_pop_used,
  const bool loy_patch_used, bool simplicity = false, bool sparse = false) {
  
  arma::ivec sge9stage3 = as<arma::ivec>(sge9l["stage3"]);
  arma::vec sge9fec32 = as<arma::vec>(sge9l["repentry3"]);
  arma::uvec sge9rep2 = as<arma::uvec>(sge9l["rep2o"]);
  arma::vec sge9ovgivent = as<arma::vec>(sge9l["ovgiven_t"]);
  arma::vec sge9ovgivenf = as<arma::vec>(sge9l["ovgiven_f"]);
  arma::vec sge9ovoffsett = as<arma::vec>(sge9l["ovoffset_t"]);
  arma::vec sge9ovoffsetf = as<arma::vec>(sge9l["ovoffset_f"]);
  arma::ivec sge9ovestt = as<arma::ivec>(sge9l["ovest_t"]);
  arma::ivec sge9ovestf = as<arma::ivec>(sge9l["ovest_f"]);
  arma::vec sge9ovsurvmult = as<arma::vec>(sge9l["ovsurvmult"]);
  arma::vec sge9ovfecmult = as<arma::vec>(sge9l["ovfecmult"]);
  arma::ivec sge9index321 = as<arma::ivec>(sge9l["index321"]);
  arma::ivec sge9index321d = as<arma::ivec>(sge9l["index321d"]);
  arma::ivec sge9index21 = as<arma::ivec>(sge9l["index21"]);
  arma::ivec aliveandequal = as<arma::ivec>(sge9l["aliveandequal"]);
  
  arma::ivec dataindex321 = as<arma::ivec>(MainData["index321"]);
  arma::ivec dataindex21 = as<arma::ivec>(MainData["pairindex21"]);
  arma::uvec dataalive3 = as<arma::uvec>(MainData["alive3"]);
  arma::vec datausedfec = as<arma::vec>(MainData["usedfec"]);
  arma::ivec dataindex2 = as<arma::ivec>(MainData["index2"]);
  arma::ivec dataindex1 = as<arma::ivec>(MainData["index1"]);
  
  int nostages = static_cast<int>(StageFrame.nrows());
  int no2stages = nostages - 1;
  int noelems = static_cast<int>(sge9index321.n_elem);
  if (format == 2) no2stages = no2stages - 1;
  int matrixdim = (nostages - 1) * no2stages;
  bool small_subset {false};
  
  arma::uvec all_repentries = find(sge9fec32 > 0.0);
  arma::ivec all_entry_stages = arma::unique(sge9stage3(all_repentries));
  int aes_count = static_cast<int>(all_entry_stages.n_elem);
  int n = static_cast<int>(dataindex321.n_elem);
  
  arma::mat dataindex321_prior(n, aes_count);
  dataindex321_prior.fill(-1);
  
  // Alternative index for fecundity under deVries format
  if (format == 2) {
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < aes_count; j++) {
        dataindex321_prior(k, j) = (all_entry_stages(j) - 1) + ((nostages - 2) * nostages) + 
          ((dataindex2(k) - 1) * nostages * nostages) + 
          ((dataindex1(k) - 1) * nostages * nostages * nostages);
      }
    }
  }
  
  int loy_length = yearorder.length();
  List A_output (loy_length);
  List U_output (loy_length);
  List F_output (loy_length);
  List conc_err (loy_length);
  List s2f_err (loy_length);
  List di321p_err (loy_length);
  
  StringVector data_pop_;
  StringVector data_patch_;
  
  if (loy_pop_used) data_pop_ = as<StringVector>(MainData[pop_var_int]);
  if (loy_patch_used) data_patch_ = as<StringVector>(MainData[patch_var_int]);
  StringVector data_year_ = as<StringVector>(MainData[year_var_int]);
  
  arma::uvec ovgiventind = find(sge9ovgivent != -1.);
  arma::uvec ovgivenfind = find(sge9ovgivenf != -1.);
  int ovgtn = static_cast<int>(ovgiventind.n_elem);
  int ovgfn = static_cast<int>(ovgivenfind.n_elem);
  
  arma::uvec ovoffsettind = find(sge9ovoffsett != 0.);
  arma::uvec ovoffsetfind = find(sge9ovoffsetf != 0.);
  int ovotn = static_cast<int>(ovoffsettind.n_elem);
  int ovofn = static_cast<int>(ovoffsetfind.n_elem);
  
  arma::uvec ovesttind = find(sge9ovestt != -1);
  arma::uvec ovestfind = find(sge9ovestf != -1);
  int ovestn = static_cast<int>(ovesttind.n_elem);
  int ovesfn = static_cast<int>(ovestfind.n_elem);
  
  for (int i = 0; i < loy_length; i++) {
    arma::uvec data_current_index = as<arma::uvec>(LefkoUtils::index_l3(data_year_, loyyear2(i)));
    
    if (loy_pop_used) {
      arma::uvec data_current_pop = as<arma::uvec>(LefkoUtils::index_l3(data_pop_, loypop(i)));
      data_current_index = intersect(data_current_index, data_current_pop);
    }
    
    if (loy_patch_used) {
      arma::uvec data_current_patch = as<arma::uvec>(LefkoUtils::index_l3(data_patch_, loypatch(i)));
      data_current_index = intersect(data_current_index, data_current_patch);
    }
    
    int data_subset_rows = static_cast<int>(data_current_index.n_elem);
    if (data_subset_rows < 10 && !small_subset) {
      small_subset = true;
      Rf_warningcall(R_NilValue,
        "Extremely small data subsets used to populate matrices.");
    }
    if (data_subset_rows == 0) continue;
    
    arma::ivec dataindex321i = dataindex321.elem(data_current_index);
    arma::ivec dataindex21i = dataindex21.elem(data_current_index);
    arma::uvec dataalive3i = dataalive3.elem(data_current_index);
    arma::vec datausedfeci = datausedfec.elem(data_current_index);
    arma::ivec dataindex2i = dataindex2.elem(data_current_index);
    arma::ivec dataindex1i = dataindex1.elem(data_current_index);
    arma::mat dataindex321_priori = dataindex321_prior.rows(data_current_index);
    
    n = static_cast<int>(dataindex321i.n_elem);
    
    arma::vec probsrates0(noelems, fill::zeros); // 1st vec: # indivs (3 trans)
    arma::vec probsrates1(noelems, fill::zeros); // 2nd vec: total indivs pair stage
    arma::vec probsrates2(noelems, fill::zeros); // 3rd vec: total indivs alive pair stage
    arma::vec probsrates3(noelems, fill::zeros); // 4th vec: total fec pair stage
    
    arma::mat stage2fec(static_cast<int>(sge3index21.n_elem), 3, fill::zeros); // col1: #inds, col2: #alive, col3: sum fec
    
    // Develop prior stage
    arma::vec probsrates0p(noelems, fill::zeros); // 1st vec: # indivs (3 trans)
    arma::vec probsrates1p(noelems, fill::zeros); // 2nd vec: total indivs pair stage
    arma::vec probsrates2p(noelems, fill::zeros); // 3rd vec: total indivs alive pair stage
    arma::vec probsrates3p(noelems, fill::zeros); // 4th vec: total fec pair stage
  
    arma::mat stage2fecp(static_cast<int>(sge3index21.n_elem), 3, fill::zeros); // col1: #inds, col2: #alive, col3: sum fec
    
    // Initialize final matrices
    arma::mat tmatrix;
    arma::mat fmatrix;
    arma::sp_mat tmatrix_sp;
    arma::sp_mat fmatrix_sp;
    
    if (sparse) { 
      arma::sp_mat tmatrix_chuck (matrixdim, matrixdim); // Main output U matrix
      arma::sp_mat fmatrix_chuck (matrixdim, matrixdim); // Main output F matrix
      
      tmatrix_sp = tmatrix_chuck;
      fmatrix_sp = fmatrix_chuck;
    } else {
      arma::mat tmatrix_chuck (matrixdim, matrixdim, fill::zeros);
      arma::mat fmatrix_chuck (matrixdim, matrixdim, fill::zeros);
      
      tmatrix = tmatrix_chuck;
      fmatrix = fmatrix_chuck;
    }
    
    // Loop counts individuals going through transitions, sums fecundities,
    // adds that info to the 3-trans and 2-trans tables
    for (int j = 0; j < n; j++) {
      // Survival portion and main individual counter
      arma::uvec choiceelement = find(sge9index321 == dataindex321i(j));
      
      arma::uvec choiceelement_sge3_vec = find(sge3index21 == dataindex21i(j));
      int choiceelement_sge3 = static_cast<int>(choiceelement_sge3_vec(0));
      
      // Indiv sum with particular transition
      stage2fec(choiceelement_sge3, 0) = stage2fec(choiceelement_sge3, 0) + 1.0; 
      
      if (choiceelement.n_elem > 0) {
        // Indiv sum with particular transition
        probsrates0(static_cast<int>(choiceelement(0))) = probsrates0(static_cast<int>(choiceelement(0))) + 1.0; 
        
        if (dataalive3i(j) > 0) {
          stage2fec(choiceelement_sge3, 1) = stage2fec(choiceelement_sge3, 1) + 1.0;
        }
      }
      
      // Fecundity sums
      arma::uvec choiceelement_f = find(sge9index321d == dataindex321i(j));
      if (choiceelement_f.n_elem > 0) {
        if (NumericVector::is_na(datausedfeci(j))) datausedfeci(j) = 0.0;
        stage2fec(choiceelement_sge3, 2) = stage2fec(choiceelement_sge3, 2) + datausedfeci(j);
      }
      
      if (format == 2) {
        for (int k = 0; k < aes_count; k++) {
          // Survival portion
          arma::uvec choiceelementp = find(sge9index321 == dataindex321_priori(j, k));
          
          arma::uvec choiceelement_sge3_vecp = find(sge3index21 == dataindex21i(j));
          int choiceelement_sge3p = static_cast<int>(choiceelement_sge3_vecp(0));
          
          stage2fecp(choiceelement_sge3p, 0) = stage2fecp(choiceelement_sge3p, 0) + 1.0;
        
          if (choiceelementp.n_elem > 0) {
            probsrates0p(choiceelementp(0)) = probsrates0p(choiceelementp(0)) + 1.0;
            
            if (dataalive3i(j) > 0) {
              stage2fecp(choiceelement_sge3p, 1) = stage2fecp(choiceelement_sge3p, 1) + 1.0;
            }
          }
          
          arma::uvec choiceelementp_f = find(sge9index321d == dataindex321_priori(j, k));
          if (choiceelementp_f.n_elem > 0) {
            stage2fecp(choiceelement_sge3p, 2) = stage2fecp(choiceelement_sge3p, 2) + datausedfeci(j);
          }
        }
      }
    }
    
    // Core data to estimate matrix elements
    for (int j = 0; j < noelems; j++) {
      int baseindex21 = sge9index21(j);
      
      if (baseindex21 > -1) {
        arma::uvec coreelementsforchoice = find(sge3index21 == baseindex21);
        unsigned int thechosenone = coreelementsforchoice(0);
        
        probsrates1(j) = stage2fec(thechosenone, 0);
        probsrates2(j) = stage2fec(thechosenone, 1);
        probsrates3(j) = stage2fec(thechosenone, 2);
        
        if (format == 2) {
          probsrates1p(j) = stage2fecp(thechosenone, 0);
          probsrates2p(j) = stage2fecp(thechosenone, 1);
          probsrates3p(j) = stage2fecp(thechosenone, 2);
        }
      }
    }
    
    // Create matrices
    for (int elem3 = 0; elem3 < noelems; elem3++) {
      if (aliveandequal(elem3) != -1) {
        if (sge9ovsurvmult(elem3) < 0.0) sge9ovsurvmult(elem3) = 1.0;
        if (!sparse) {
          tmatrix(aliveandequal(elem3)) = probsrates0(elem3)* sge9ovsurvmult(elem3) /
            probsrates1(elem3); // Survival
        } else {
          tmatrix_sp(aliveandequal(elem3)) = probsrates0(elem3)* sge9ovsurvmult(elem3) /
            probsrates1(elem3); // Survival
        }
        
        // Fecundity
        if (sge9ovfecmult(elem3) < 0.0) sge9ovfecmult(elem3) = 1.0;
        if (format == 2) {
          if (!sparse) {
            fmatrix(aliveandequal(elem3)) = (sge9fec32(elem3)) *
              static_cast<double>(sge9rep2(elem3)) * probsrates3p(elem3) *
              sge9ovfecmult(elem3) / probsrates1p(elem3);
          } else {
            fmatrix_sp(aliveandequal(elem3)) = (sge9fec32(elem3)) *
              static_cast<double>(sge9rep2(elem3)) * probsrates3p(elem3) *
              sge9ovfecmult(elem3) / probsrates1p(elem3);
          }
        } else {
          if (!sparse) {
            fmatrix(aliveandequal(elem3)) = (sge9fec32(elem3)) *
              static_cast<double>(sge9rep2(elem3)) * probsrates3(elem3) *
              sge9ovfecmult(elem3) / probsrates1(elem3);
          } else {
            fmatrix_sp(aliveandequal(elem3)) = (sge9fec32(elem3)) *
              static_cast<double>(sge9rep2(elem3)) * probsrates3(elem3) *
              sge9ovfecmult(elem3) / probsrates1(elem3);
          }
        }
      }
    }
    
    // Correct transitions for given stuff
    if (ovgtn > 0) {
      for (int j = 0; j < ovgtn; j++) {
        int matrixelement2 = aliveandequal(ovgiventind(j));
        if (!sparse) {
          tmatrix(matrixelement2) = sge9ovgivent(ovgiventind(j));
        } else {
          tmatrix_sp(matrixelement2) = sge9ovgivent(ovgiventind(j));
        }
      }
    }
    
    if (ovgfn > 0) {
      for (int j = 0; j < ovgfn; j++) {
        int matrixelement2 = aliveandequal(ovgivenfind(j));
        if (!sparse) {
          fmatrix(matrixelement2) = sge9ovgivenf(ovgivenfind(j));
        } else {
          fmatrix_sp(matrixelement2) = sge9ovgivenf(ovgivenfind(j));
        }
      }
    }
    
    // Correct transitions for additive offsets
    if (ovotn > 0) {
      for (int j = 0; j < ovotn; j++) {
        int matrixelement2 = aliveandequal(ovoffsettind(j));
        if (NumericVector::is_na(tmatrix(matrixelement2)) ||
          Rcpp::traits::is_nan<REALSXP>(tmatrix(matrixelement2))) {tmatrix(matrixelement2) = 0.;}
        if (!sparse) {
          tmatrix(matrixelement2) = tmatrix(matrixelement2) + sge9ovoffsett(ovoffsettind(j));
        } else {
          tmatrix_sp(matrixelement2) = tmatrix_sp(matrixelement2) + sge9ovoffsett(ovoffsettind(j));
        }
      }
    }
    
    if (ovofn > 0) {
      for (int j = 0; j < ovofn; j++) {
        int matrixelement2 = aliveandequal(ovoffsetfind(j));
        if (NumericVector::is_na(fmatrix(matrixelement2)) ||
          Rcpp::traits::is_nan<REALSXP>(fmatrix(matrixelement2))) {fmatrix(matrixelement2) = 0.;}
        if (!sparse) {
          fmatrix(matrixelement2) = fmatrix(matrixelement2) + sge9ovoffsetf(ovoffsetfind(j));
        } else {
          fmatrix_sp(matrixelement2) = fmatrix_sp(matrixelement2) + sge9ovoffsetf(ovoffsetfind(j));
        }
      }
    }
    
    // Replace transitions for proxy values as given in overwrite table  
    if (ovestn > 0) {
      for (int j = 0; j < ovestn; j++) {
        arma::uvec replacement = find(sge9index321 == sge9ovestt(ovesttind(j)));
        
        if (replacement.n_elem > 0) {
          double correction = sge9ovsurvmult(ovesttind(j));
          if (correction == -1.0) correction = 1.0;
          if (!sparse) {
            tmatrix(aliveandequal(ovesttind(j))) = tmatrix(aliveandequal(replacement(0))) *
              correction;
          } else {
            tmatrix_sp(aliveandequal(ovesttind(j))) = tmatrix_sp(aliveandequal(replacement(0))) *
              correction;
          }
        }
      }
    }
    
    if (ovesfn > 0) {
      for (int j = 0; j < ovesfn; j++) {
        arma::uvec replacement = find(sge9index321 == sge9ovestf(ovestfind(j)));
        
        if (replacement.n_elem > 0) {
          double correction = sge9ovfecmult(ovesttind(j));
          if (correction == -1.0) correction = 1.0;
          if (!sparse) {
            fmatrix(aliveandequal(ovestfind(j))) = fmatrix(aliveandequal(replacement(0))) *
              correction;
          } else {
            fmatrix_sp(aliveandequal(ovestfind(j))) = fmatrix_sp(aliveandequal(replacement(0))) *
              correction;
          }
        }
      }
    }
    
    if (!sparse) {
      tmatrix(find_nonfinite(tmatrix)).zeros();
      fmatrix(find_nonfinite(fmatrix)).zeros();
      
      U_output(i) = tmatrix;
      F_output(i) = fmatrix;
      
      if (!simplicity) {
        arma::mat amatrix = tmatrix + fmatrix;
        A_output(i) = amatrix;
      }
    } else {
      arma::uvec tm_sp_nan = find_nonfinite(tmatrix_sp);
      arma::uvec fm_sp_nan = find_nonfinite(fmatrix_sp);
      
      for (int nan_count = 0; nan_count < static_cast<int>(tm_sp_nan.n_elem); nan_count++) { 
        tmatrix_sp(tm_sp_nan(nan_count)) = 0.0;
      }
      for (int nan_count = 0; nan_count < static_cast<int>(fm_sp_nan.n_elem); nan_count++) { 
        fmatrix_sp(fm_sp_nan(nan_count)) = 0.0;
      }
      
      U_output(i) = tmatrix_sp;
      F_output(i) = fmatrix_sp;
      
      if (!simplicity) {
        arma::sp_mat amatrix_sp = tmatrix_sp + fmatrix_sp; // Create the A matrix
        A_output(i) = amatrix_sp;
      }
    }
    
    if (err_switch) {
      DataFrame concatenated_crap = DataFrame::create(_["sge9index321"] = sge9index321,
        _["aliveandequal"] = aliveandequal, _["probsrates0"] = probsrates0,
        _["probsrates1"] = probsrates1, _["probsrates2"] = probsrates2,
        _["probsrates3"] = probsrates3, _["sge9fec32"] = sge9fec32,
        _["probsrates0p"] = probsrates0p, _["probsrates1p"] = probsrates1p,
        _["probsrates2p"] = probsrates2p, _["probsrates3p"] = probsrates3p);
      
      arma::mat s2f = arma::join_horiz(stage2fec, stage2fecp);
      
      conc_err(i) = concatenated_crap;
      s2f_err(i) = s2f;
      di321p_err(i) = dataindex321_prior;
    }
  }
  
  List final_output;
  
  if (err_switch) {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue,
      _["probsrates_all"] = conc_err, _["s2f"] = s2f_err,
      _["dataprior"] = di321p_err, _["data"] = MainData);
    final_output = out_dude;
    
  } else {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue);
    final_output = out_dude;
  }
  
  return final_output;
}

//' Estimate All Elements of Raw Ahistorical Population Projection Matrix
//' 
//' Function \code{normalpatrolgroup()} swiftly calculates matrix transitions
//' in raw ahistorical matrices, and serves as the core workhorse function
//' behind \code{rlefko2()}.
//' 
//' @name normalpatrolgroup
//' 
//' @param sge3 The Allstages data frame developed for \code{rlefko2()} covering
//' stage pairs across times \emph{t}+1 and \emph{t}. Generally termed
//' \code{stageexpansion3}.
//' @param sge2stage2 An integer index vector giving the stage in time \emph{t},
//' from \code{stageexpansion2}.
//' @param MainData The demographic dataset modified to hold \code{usedfec} and
//' \code{usedstage} columns.
//' @param StageFrame The full stageframe for the analysis.
//' @param err_switch A logical value. If set to \code{TRUE}, then will also
//' output \code{probsrates}, \code{stage2fec}, and the edited dataset.
//' @param loypop A string vector giving the order of populations in the list of
//' years.
//' @param loypatch A string vector giving the order of patches in the list of
//' years.
//' @param loyyear2 A string vector giving the order of years at time \emph{t}
//' in the list of years.
//' @param yearorder The integer year order corresponding to \code{loyyear2}.
//' @param pop_var_int The variable number coding for population in the main
//' data set.
//' @param patch_var_int The variable number coding for patch in the main data
//' set.
//' @param year_var_int The variable number coding for year in time \emph{t} in
//' the main data set.
//' @param loy_pop_used A logical value indicating whether the population
//' variable is to be used.
//' @param loy_patch_used A logical value indicating whether the patch variable
//' is to be used.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' @param sparse If \code{TRUE}, then output will be in sparse matrix format.
//' Defaults to \code{FALSE}.
//' 
//' @return In the standard output, a list of three lists, called \code{A},
//' \code{U}, and \code{F}, each containing A, U, or F matrices, respectively.
//' Further information on the output is provided in the documentation for
//' \code{specialpatrolgroup()}.
//' 
//' @keywords internal
//' @noRd
List normalpatrolgroup(const DataFrame& sge3, const arma::ivec& sge2stage2, 
  const DataFrame& MainData, const DataFrame& StageFrame, int err_switch,
  StringVector loypop, StringVector loypatch, StringVector loyyear2,
  IntegerVector yearorder, const int pop_var_int, const int patch_var_int,
  const int year_var_int, const bool loy_pop_used, const bool loy_patch_used,
  bool simplicity = false, bool sparse = false) {
  
  arma::vec sge3fec32 = as<arma::vec>(sge3["repentry3"]);
  arma::uvec sge3rep2 = as<arma::uvec>(sge3["rep2n"]);
  arma::vec sge3ovgivent = as<arma::vec>(sge3["ovgiven_t"]);
  arma::ivec sge3ovgivenf = as<arma::ivec>(sge3["ovgiven_f"]);
  arma::vec sge3ovoffsett = as<arma::vec>(sge3["ovoffset_t"]);
  arma::vec sge3ovoffsetf = as<arma::vec>(sge3["ovoffset_f"]);
  arma::ivec sge3ovestt = as<arma::ivec>(sge3["ovest_t"]);
  arma::ivec sge3ovestf = as<arma::ivec>(sge3["ovest_f"]);
  arma::vec sge3ovsurvmult = as<arma::vec>(sge3["ovsurvmult"]);
  arma::vec sge3ovfecmult = as<arma::vec>(sge3["ovfecmult"]);
  arma::ivec sge3index32 = as<arma::ivec>(sge3["index321"]);
  arma::ivec aliveandequal = as<arma::ivec>(sge3["aliveandequal"]);
  
  arma::ivec dataindex32 = as<arma::ivec>(MainData["index32"]);
  arma::ivec dataindex2 = as<arma::ivec>(MainData["index2"]);
  arma::uvec dataalive3 = as<arma::uvec>(MainData["alive3"]);
  arma::vec datausedfec = as<arma::vec>(MainData["usedfec"]);
  
  int nostages = static_cast<int>(StageFrame.nrows());
  int no2stages = static_cast<int>(sge2stage2.n_elem) - 1; // Removes dead stage
  int noelems = static_cast<int>(sge3index32.n_elem);
  bool small_subset {false};
  
  int loy_length = yearorder.length();
  List A_output (loy_length);
  List U_output (loy_length);
  List F_output (loy_length);
  List conc_err (loy_length);
  List s2f_err (loy_length);

  StringVector data_pop_;
  StringVector data_patch_;
  
  if (loy_pop_used) data_pop_ = as<StringVector>(MainData[pop_var_int]);
  if (loy_patch_used) data_patch_ = as<StringVector>(MainData[patch_var_int]);
  StringVector data_year_ = as<StringVector>(MainData[year_var_int]);
  
  arma::uvec ovgiventind = find(sge3ovgivent != -1.);
  arma::uvec ovgivenfind = find(sge3ovgivenf != -1.);
  int ovgtn = static_cast<int>(ovgiventind.n_elem);
  int ovgfn = static_cast<int>(ovgivenfind.n_elem);
  
  arma::uvec ovoffsettind = find(sge3ovoffsett != 0.);
  arma::uvec ovoffsetfind = find(sge3ovoffsetf != 0.);
  int ovotn = static_cast<int>(ovoffsettind.n_elem);
  int ovofn = static_cast<int>(ovoffsetfind.n_elem);
  
  arma::uvec ovesttind = find(sge3ovestt != -1);
  arma::uvec ovestfind = find(sge3ovestf != -1);
  int ovestn = static_cast<int>(ovesttind.n_elem);
  int ovesfn = static_cast<int>(ovestfind.n_elem);
  
  for (int i = 0; i < loy_length; i++) {
    arma::uvec data_current_index = as<arma::uvec>(LefkoUtils::index_l3(data_year_, loyyear2(i)));
    
    if (loy_pop_used) {
      arma::uvec data_current_pop = as<arma::uvec>(LefkoUtils::index_l3(data_pop_, loypop(i)));
      data_current_index = intersect(data_current_index, data_current_pop);
    }
    
    if (loy_patch_used) {
      arma::uvec data_current_patch = as<arma::uvec>(LefkoUtils::index_l3(data_patch_, loypatch(i)));
      data_current_index = intersect(data_current_index, data_current_patch);
    }
    
    int data_subset_rows = static_cast<int>(data_current_index.n_elem);
    if (data_subset_rows < 10 && !small_subset) {
      small_subset = true;
      Rf_warningcall(R_NilValue,
        "Extremely small data subsets used to populate matrices.");
    }
    if (data_subset_rows == 0) continue;
    
    arma::ivec dataindex32i = dataindex32.elem(data_current_index);
    arma::ivec dataindex2i = dataindex2.elem(data_current_index);
    arma::uvec dataalive3i = dataalive3.elem(data_current_index);
    arma::vec datausedfeci = datausedfec.elem(data_current_index);
    
    dataindex32i.resize(data_subset_rows);
    dataindex2i.resize(data_subset_rows);
    dataalive3i.resize(data_subset_rows);
    datausedfeci.resize(data_subset_rows);
    
    int n = static_cast<int>(dataindex32i.n_elem);
    
    arma::vec probsrates0(noelems, fill::zeros); // # indivs (3 trans)
    arma::vec probsrates1(noelems, fill::zeros); // total indivs pair stage
    arma::vec probsrates2(noelems, fill::zeros); // total indivs alive pair stage
    arma::vec probsrates3(noelems, fill::zeros); // total fec pair stage
    
    // col1: # inds total, col2: no alive, col3: sum fec
    arma::mat stage2fec(no2stages, 3, fill::zeros);
    
    arma::mat tmatrix;
    arma::mat fmatrix;
    arma::sp_mat tmatrix_sp;
    arma::sp_mat fmatrix_sp;
    
    if (!sparse) {
      arma::mat tmatrix_chuck((nostages-1), (nostages-1), fill::zeros);
      arma::mat fmatrix_chuck((nostages-1), (nostages-1), fill::zeros);
      
      tmatrix = tmatrix_chuck;
      fmatrix = fmatrix_chuck;
    } else { 
      arma::sp_mat tmatrix_chuck((nostages-1), (nostages-1));
      arma::sp_mat fmatrix_chuck((nostages-1), (nostages-1));
      
      tmatrix_sp = tmatrix_chuck;
      fmatrix_sp = fmatrix_chuck;
    }
    
    // Count individuals through transitions, sum fecundities,
    // add 3-trans and 2-trans tables
    for (int j = 0; j < n; j++) { 
      // Sum all individuals with particular transition
      arma::uvec chosen_sge3index32vec = find(sge3index32 == dataindex32i(j));
      
      int chosen_sge3index32 {0};
      if (chosen_sge3index32vec.n_elem != 0) {
        chosen_sge3index32 = static_cast<int>(chosen_sge3index32vec(0));
        
        probsrates0(chosen_sge3index32) = probsrates0(chosen_sge3index32) + 1.0;
      }
      
      // Sum all individuals with particular transition
      if (dataindex2i(j) >= no2stages) {
        throw Rcpp::exception("Current stageframe does not hold all trait combinations in dataset.", false);
      }
      stage2fec((dataindex2i(j)), 0) = stage2fec((dataindex2i(j)), 0) + 1.0; 
      if (dataalive3i(j) > 0) {
        stage2fec((dataindex2i(j)), 1) = stage2fec((dataindex2i(j)), 1) + 1.0;
      }
      
      if (NumericVector::is_na(datausedfeci(j))) datausedfeci(j) = 0.0;
      stage2fec((dataindex2i(j)), 2) = stage2fec((dataindex2i(j)), 2) + datausedfeci(j);
    }
    
    // Populate vectors of individuals by stage in time t
    for (int k = 0; k < no2stages; k++) {
      unsigned int foradding = ((sge2stage2(k) - 1) * nostages);
      
      for (int j = 0; j < nostages; j++) {
        unsigned int entry = foradding + j;
        
        probsrates1(entry) = stage2fec(k, 0);
        probsrates2(entry) = stage2fec(k, 1);
        probsrates3(entry) = stage2fec(k, 2);
      }
    }
    
    // Populate main U and F matrices
    for (int elem3 = 0; elem3 < noelems; elem3++) {
      if (aliveandequal(elem3) != -1) {
        
        // Leave NaNs when 0 individuals summed through in probsrates1
        if (sge3ovsurvmult(elem3) < 0.0) sge3ovsurvmult(elem3) = 1.0;
        if (!sparse) {
          tmatrix(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
            probsrates1(elem3);
        } else {
          tmatrix_sp(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
            probsrates1(elem3);
        }
          
        if (sge3ovfecmult(elem3) < 0.0) sge3ovfecmult(elem3) = 1.0;
        if (!sparse) {
          fmatrix(aliveandequal(elem3)) = sge3fec32(elem3) *
            static_cast<double>(sge3rep2(elem3)) * probsrates3(elem3) *
            sge3ovfecmult(elem3) / probsrates1(elem3);
        } else {
          fmatrix_sp(aliveandequal(elem3)) = sge3fec32(elem3) *
            static_cast<double>(sge3rep2(elem3)) * probsrates3(elem3) *
            sge3ovfecmult(elem3) / probsrates1(elem3);
        }
      }
    }
    
    // Correct for transitions given in overwrite table
    if (ovgtn > 0) {
      for (int j = 0; j < ovgtn; j++) {
        int matrixelement2 = aliveandequal(ovgiventind(j));
        if (!sparse) {
          tmatrix(matrixelement2) = sge3ovgivent(ovgiventind(j));
        } else {
          tmatrix_sp(matrixelement2) = sge3ovgivent(ovgiventind(j));
        }
      }
    }
    
    if (ovgfn > 0) {
      for (int j = 0; j < ovgfn; j++) {
        int matrixelement2 = aliveandequal(ovgivenfind(j));
        if (!sparse) {
          fmatrix(matrixelement2) = sge3ovgivenf(ovgivenfind(j));
        } else {
          fmatrix_sp(matrixelement2) = sge3ovgivenf(ovgivenfind(j));
        }
      }
    }
    
    // Correct for additive offsets to transitions
    if (ovotn > 0) {
      for (int j = 0; j < ovotn; j++) {
        int matrixelement2 = aliveandequal(ovoffsettind(j));
        if (NumericVector::is_na(tmatrix(matrixelement2)) ||
          Rcpp::traits::is_nan<REALSXP>(tmatrix(matrixelement2))) {tmatrix(matrixelement2) = 0.;}
        if (!sparse) {
          tmatrix(matrixelement2) = tmatrix(matrixelement2) + sge3ovoffsett(ovoffsettind(j));
        } else {
          tmatrix_sp(matrixelement2) = tmatrix_sp(matrixelement2) + sge3ovoffsett(ovoffsettind(j));
        }
      }
    }
    
    if (ovofn > 0) {
      for (int j = 0; j < ovofn; j++) {
        int matrixelement2 = aliveandequal(ovoffsetfind(j));
        if (NumericVector::is_na(fmatrix(matrixelement2)) ||
          Rcpp::traits::is_nan<REALSXP>(fmatrix(matrixelement2))) {fmatrix(matrixelement2) = 0.;}
        if (!sparse) {
          fmatrix(matrixelement2) = fmatrix(matrixelement2) + sge3ovoffsetf(ovoffsetfind(j));
        } else {
          fmatrix_sp(matrixelement2) = fmatrix_sp(matrixelement2) + sge3ovoffsetf(ovoffsetfind(j));
        }
      }
    }
    
    // Replace transitions with proxies as in overwrite table
    if (ovestn > 0) {
      for (int j = 0; j < ovestn; j++) {
        arma::uvec replacement = find(sge3index32 == sge3ovestt(ovesttind(j)));
        
        double correction = sge3ovsurvmult(ovesttind(j));
        if (correction == -1.0) correction = 1.0;
        if (!sparse) {
          tmatrix(aliveandequal(ovesttind(j))) = tmatrix(aliveandequal(replacement(0))) *
            correction;
        } else {
          tmatrix_sp(aliveandequal(ovesttind(j))) = tmatrix_sp(aliveandequal(replacement(0))) *
            correction;
        }
      }
    }
    
    if (ovesfn > 0) {
      for (int j = 0; j < ovesfn; j++) {
        arma::uvec replacement = find(sge3index32 == sge3ovestf(ovestfind(j)));
        
        double correction = sge3ovfecmult(ovesttind(j));
        if (correction == -1.0) correction = 1.0;
        if (!sparse) {
          fmatrix(aliveandequal(ovestfind(j))) = fmatrix(aliveandequal(replacement(0))) *
            correction;
        } else {
          fmatrix_sp(aliveandequal(ovestfind(j))) = fmatrix_sp(aliveandequal(replacement(0))) *
            correction;
        }
      }
    }
    
    if (!sparse) {
      tmatrix(find_nonfinite(tmatrix)).zeros();
      fmatrix(find_nonfinite(fmatrix)).zeros();
      
      U_output(i) = tmatrix;
      F_output(i) = fmatrix;
      
      if (!simplicity) {
        arma::mat amatrix = tmatrix + fmatrix;
        A_output(i) = amatrix;
      }
    } else {
      arma::uvec tm_sp_nan = find_nonfinite(tmatrix_sp);
      arma::uvec fm_sp_nan = find_nonfinite(fmatrix_sp);
      
      for (int nan_count = 0; nan_count < static_cast<int>(tm_sp_nan.n_elem); nan_count++) { 
        tmatrix_sp(tm_sp_nan(nan_count)) = 0.0;
      }
      for (int nan_count = 0; nan_count < static_cast<int>(fm_sp_nan.n_elem); nan_count++) { 
        fmatrix_sp(fm_sp_nan(nan_count)) = 0.0;
      }
      
      U_output(i) = tmatrix_sp;
      F_output(i) = fmatrix_sp;
      
      if (!simplicity) {
        arma::sp_mat amatrix_sp = tmatrix_sp + fmatrix_sp;
        A_output(i) = amatrix_sp;
      }
    }
    
    if (err_switch) {
      DataFrame concatenated_crap = DataFrame::create(_["sge3index32"] = sge3index32,
        _["aliveandequal"] = aliveandequal, _["probsrates0"] = probsrates0,
        _["probsrates1"] = probsrates1, _["probsrates2"] = probsrates2,
        _["probsrates3"] = probsrates3, _["sge3fec32"] = sge3fec32);
      
      conc_err(i) = concatenated_crap;
      s2f_err(i) = stage2fec;
    }
  }
  
  List final_output;
  
  if (err_switch) {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue,
      _["probsrates_all"] = conc_err, _["s2f"] = s2f_err, _["data"] = MainData);
    final_output = out_dude;
    
  } else {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue);
    final_output = out_dude;
  }
  
  return final_output;
}

//' Estimate All Elements of Raw Ahistorical Population Projection Matrix
//' 
//' Function \code{minorpatrolgroup()} swiftly calculates matrix transitions
//' in raw Leslie MPMs, and is used internally in \code{rleslie()}.
//' 
//' @name minorpatrolgroup
//' 
//' @param MainData The demographic dataset modified internally to have needed
//' variables for living status, reproduction status, and fecundity.
//' @param StageFrame The full stageframe for the analysis.
//' @param supplement A supplement table as produced by function
//' \code{supplemental()} and edited and age-expanded by pre-processing within
//' function \code{mpm_create()}.
//' @param start_age An integer denoting the first age incorporated in the MPM.
//' @param last_age An integer denoting the last age incorporated in the MPM,
//' not including ages set to equal the last estimated age.
//' @param cont Should a self-loop transition be estimated for the final age.
//' @param fec_mod A multiplier on raw fecundity to estimate true fecundity.
//' @param err_switch A logical value. If set to \code{TRUE}, then will also
//' output \code{probsrates} and \code{stage2fec}.
//' @param loypop A string vector giving the order of populations in the list of
//' years.
//' @param loypatch A string vector giving the order of patches in the list of
//' years.
//' @param loyyear2 A string vector giving the order of years at time \emph{t}
//' in the list of years.
//' @param yearorder The integer year order corresponding to \code{loyyear2}.
//' @param pop_var_int The variable number coding for population in the main
//' data set.
//' @param patch_var_int The variable number coding for patch in the main data
//' set.
//' @param year_var_int The variable number coding for year in time \emph{t} in
//' the main data set.
//' @param loy_pop_used A logical value indicating whether the population
//' variable is to be used.
//' @param loy_patch_used A logical value indicating whether the patch variable
//' is to be used.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' @param sparse If \code{TRUE}, then output will be in sparse matrix format.
//' Defaults to \code{FALSE}.
//' 
//' @return List of three matrices, including the survival-transition (\code{U})
//' matrix, the fecundity matrix (\code{F}), and the sum (\code{A}) matrix, with
//' the \code{A} matrix first. Further information on outputs is provided in the
//' documentation for \code{specialpatrolgroup()}.
//' 
//' @keywords internal
//' @noRd
Rcpp::List minorpatrolgroup(const DataFrame& MainData,
  const DataFrame& StageFrame, const DataFrame& supplement, int start_age,
  int last_age, bool cont, double fec_mod, int err_switch, StringVector loypop,
  StringVector loypatch, StringVector loyyear2, IntegerVector yearorder,
  const int pop_var_int, const int patch_var_int, const int year_var_int,
  const bool loy_pop_used, const bool loy_patch_used, bool simplicity = false,
  bool sparse = false) {
  
  arma::ivec dataage = as<arma::ivec>(MainData["usedage"]);
  arma::uvec dataalive2 = as<arma::uvec>(MainData["alive2"]);
  arma::uvec dataalive3 = as<arma::uvec>(MainData["alive3"]);
  arma::vec datausedfec = as<arma::vec>(MainData["usedfec"]);
  
  IntegerVector sf_minage = StageFrame["min_age"];
  IntegerVector sf_repstatus = StageFrame["repstatus"];
  int noages = sf_minage.length();
  bool small_subset {false};
  
  IntegerVector ov_age2;
  IntegerVector ov_estage2;
  NumericVector ov_givenrate;
  NumericVector ov_offset;
  NumericVector ov_multiplier;
  IntegerVector ov_convtype;
  int supp_length {0};
  
  if (supplement.containsElementNamed("age2")) {
    ov_age2 = as<IntegerVector>(supplement["age2"]);
    ov_estage2 = as<IntegerVector>(supplement["estage2"]);
    ov_givenrate = as<NumericVector>(supplement["givenrate"]);
    ov_offset = as<NumericVector>(supplement["offset"]);
    ov_multiplier = as<NumericVector>(supplement["multiplier"]);
    ov_convtype = as<IntegerVector>(supplement["convtype"]);
    
    supp_length = static_cast<int>(ov_givenrate.length());
    
    for (int i = 0; i < supp_length; i++) {
      ov_age2(i) = ov_age2(i) - start_age;
      
      if (!IntegerVector::is_na(ov_estage2(i))) {
        ov_estage2(i) = ov_estage2(i) - start_age;
      }
    }
  }
  
  int loy_length = yearorder.length();
  List A_output (loy_length);
  List U_output (loy_length);
  List F_output (loy_length);
  List conc_err (loy_length);
  
  StringVector data_pop_;
  StringVector data_patch_;
  
  if (loy_pop_used) data_pop_ = as<StringVector>(MainData[pop_var_int]);
  if (loy_patch_used) data_patch_ = as<StringVector>(MainData[patch_var_int]);
  StringVector data_year_ = as<StringVector>(MainData[year_var_int]);
  
  for (int i = 0; i < loy_length; i++) {
    arma::uvec data_current_index = as<arma::uvec>(LefkoUtils::index_l3(data_year_, loyyear2(i)));
    
    if (loy_pop_used) {
      arma::uvec data_current_pop = as<arma::uvec>(LefkoUtils::index_l3(data_pop_, loypop(i)));
      data_current_index = intersect(data_current_index, data_current_pop);
    }
    
    if (loy_patch_used) {
      arma::uvec data_current_patch = as<arma::uvec>(LefkoUtils::index_l3(data_patch_, loypatch(i)));
      data_current_index = intersect(data_current_index, data_current_patch);
    }
    
    int data_subset_rows = static_cast<int>(data_current_index.n_elem);
    if (data_subset_rows < 10 && !small_subset) {
      small_subset = true;
      Rf_warningcall(R_NilValue,
        "Extremely small data subsets used to populate matrices.");
    }
    if (data_subset_rows == 0) continue;
    
    arma::ivec dataagei = dataage.elem(data_current_index);
    arma::uvec dataalive2i = dataalive2.elem(data_current_index);
    arma::uvec dataalive3i = dataalive3.elem(data_current_index);
    arma::vec datausedfeci = datausedfec.elem(data_current_index);
    
    dataagei.resize(data_subset_rows);
    dataalive2i.resize(data_subset_rows);
    dataalive3i.resize(data_subset_rows);
    datausedfeci.resize(data_subset_rows);
    
    arma::vec probsrates0(noages, fill::zeros); // total indivs in t
    arma::vec probsrates1(noages, fill::zeros); // total indivs alive in t+1
    arma::vec probsrates2(noages, fill::zeros); // total fec
    
    arma::mat tmatrix;
    arma::mat fmatrix;
    arma::sp_mat tmatrix_sp;
    arma::sp_mat fmatrix_sp;
    
    if (!sparse) {
      arma::mat tmatrix_chuck (noages, noages, fill::zeros);
      arma::mat fmatrix_chuck (noages, noages, fill::zeros);
      
      tmatrix = tmatrix_chuck;
      fmatrix = fmatrix_chuck;
    } else {
      arma::sp_mat tmatrix_chuck (noages, noages);
      arma::sp_mat fmatrix_chuck (noages, noages);
      
      tmatrix_sp = tmatrix_chuck;
      fmatrix_sp = fmatrix_chuck;
    }
    
    // Count individuals in transitions, calculate survival and fecundity
    arma::uvec data_allalivei = find(dataalive2i);
    int survsum {0};
    double fecsum {0.0};
    
    for (int k = 0; k < noages; k++) {
      arma::uvec data_indices;
      arma::uvec aget_alive;
      int num_aget_alive {0};
      
      if (k < (noages - 1) || !cont) {
        data_indices = find(dataagei == sf_minage(k));
      } else { 
        data_indices = find(dataagei >= sf_minage(k));
      }
      
      aget_alive = intersect(data_allalivei, data_indices);
      num_aget_alive = static_cast<int>(aget_alive.n_elem);
      
      if (num_aget_alive > 0) {
        for (int j = 0; j < num_aget_alive; j++) {
          if (dataalive3i(aget_alive(j)) > 0) survsum++;
          fecsum = fecsum + datausedfeci(aget_alive(j));
        }
      }
      
      probsrates0(k) = num_aget_alive;
      if (num_aget_alive > 0) {
        probsrates1(k) = static_cast<double>(survsum) / static_cast<double>(num_aget_alive);
        if (sf_repstatus(k) > 0) {
          probsrates2(k) = fec_mod * fecsum / static_cast<double>(num_aget_alive);
        }
      } else {
        probsrates1(k) = 0.0;
        probsrates2(k) = 0.0;
      }
      
      if (!sparse) {
        if (k < (noages - 1)) {
          tmatrix(k+1, k) = probsrates1(k);
        } else if (cont) {
          tmatrix(k, k) = probsrates1(k);
        }
        
        fmatrix(0, k) = probsrates2(k);
      } else {
        if (k < (noages - 1)) {
          tmatrix_sp(k+1, k) = probsrates1(k);
        } else if (cont) {
          tmatrix_sp(k, k) = probsrates1(k);
        }
        
        fmatrix_sp(0, k) = probsrates2(k);
      }
      
      survsum = 0;
      fecsum = 0.0;
    }
    
    // Supplement replacement portion
    int target_col {0};
    int target_row {0};
    int proxy_col {0};
    int proxy_row {0};
    
    for (int l = 0; l < supp_length; l++) {
      target_col = ov_age2(l);
      if (target_col >= (noages - 1) && cont) {
        target_col = noages - 1;
      }
      
      if (ov_convtype(l) == 1) {
        if (target_col >= (noages - 1) && cont) {
          target_row = target_col;
        } else {
          target_row = target_col + 1;
        }
        
        if (!NumericVector::is_na(ov_givenrate(l))) {
          tmatrix(target_row, target_col) = ov_givenrate(l);
        }
        if (!NumericVector::is_na(ov_offset(l))) {
          if (NumericVector::is_na(tmatrix(target_row, target_col)) ||
              Rcpp::traits::is_nan<REALSXP>(tmatrix(target_row, target_col))) {
            tmatrix(target_row, target_col) = 0.;
          }
          if (ov_offset(l) != 0.) tmatrix(target_row, target_col) = 
            tmatrix(target_row, target_col) + ov_offset(l);
        }
        if (!IntegerVector::is_na(ov_estage2(l))) {
          proxy_col = ov_estage2(l);
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_col = noages - 1;
          }
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_row = proxy_col; // Used to be target_col
          } else {
            proxy_row = proxy_col + 1;
          }
          
          tmatrix(target_row, target_col) = tmatrix(proxy_row, proxy_col);
        }
        if (!NumericVector::is_na(ov_multiplier(l))) {
          tmatrix(target_row, target_col) *= ov_multiplier(l);
        }
      } else if (ov_convtype(l) == 2) {
        target_row = 0;
        
        if (!NumericVector::is_na(ov_givenrate(l))) {
          fmatrix(target_row, target_col) = ov_givenrate(l);
        }
        if (!NumericVector::is_na(ov_offset(l))) {
          if (NumericVector::is_na(fmatrix(target_row, target_col)) ||
              Rcpp::traits::is_nan<REALSXP>(fmatrix(target_row, target_col))) {
            fmatrix(target_row, target_col) = 0.;
          }
          if (ov_offset(l) != 0.) fmatrix(target_row, target_col) = 
            fmatrix(target_row, target_col) + ov_offset(l);
        }
        if (!IntegerVector::is_na(ov_estage2(l))) {
          proxy_col = ov_estage2(l);
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_col = noages - 1;
          }
          
          proxy_row = 0;
          
          fmatrix(target_row, target_col) = fmatrix(proxy_row, proxy_col);
        }
        if (!NumericVector::is_na(ov_multiplier(l))) {
          fmatrix(target_row, target_col) *= ov_multiplier(l);
        }
      } else {
        target_row = 0;
        
        if (!NumericVector::is_na(ov_multiplier(l))) {
          fmatrix(target_row, target_col) *= ov_multiplier(l);
        }
      }
    }
    
    if (!sparse) {
      tmatrix(find_nonfinite(tmatrix)).zeros();
      fmatrix(find_nonfinite(fmatrix)).zeros();
      
      U_output(i) = tmatrix;
      F_output(i) = fmatrix;
      
      if (!simplicity) {
        arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
        A_output(i) = amatrix;
      }
    } else {
      arma::uvec tm_sp_nan = find_nonfinite(tmatrix_sp);
      arma::uvec fm_sp_nan = find_nonfinite(fmatrix_sp);
      
      for (int nan_count = 0; nan_count < static_cast<int>(tm_sp_nan.n_elem); nan_count++) { 
        tmatrix_sp(tm_sp_nan(nan_count)) = 0.0;
      }
      for (int nan_count = 0; nan_count < static_cast<int>(fm_sp_nan.n_elem); nan_count++) { 
        fmatrix_sp(fm_sp_nan(nan_count)) = 0.0;
      }
      
      U_output(i) = tmatrix_sp;
      F_output(i) = fmatrix_sp;
      
      if (!simplicity) {
        arma::sp_mat amatrix_sp = tmatrix_sp + fmatrix_sp;
        A_output(i) = amatrix_sp;
      }
    }
    
    if (err_switch) {
      DataFrame concatenated_crap = DataFrame::create(_["probsrates0"] = probsrates0,
        _["probsrates1"] = probsrates1, _["probsrates2"] = probsrates2);
      
      conc_err(i) = concatenated_crap;
    }
  }
  
  List final_output;
  
  if (err_switch) {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue,
      _["probsrates_all"] = conc_err, _["data"] = MainData);
    final_output = out_dude;
    
  } else {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue);
    final_output = out_dude;
  }
  
  return final_output;
}

//' Estimate All Elements of Raw Age-By-Stage Population Projection Matrix
//' 
//' Function \code{subvertedpatrolgroup()} swiftly calculates matrix
//' transitions in raw age-by-stage matrices, and serves as the core workhorse
//' function behind \code{arlefko2()}.
//' 
//' @name subvertedpatrolgroup
//' 
//' @param sge3 The Allstages data frame developed for \code{rlefko2()} covering
//' stage pairs across times \emph{t}+1 and \emph{t}. Generally termed
//' \code{stageexpansion3}.
//' @param sge2index21 An integer index vector of stage in times \emph{t}-1 and
//' \emph{t} from \code{stageexpansion2}.
//' @param MainData The demographic dataset modified to hold \code{usedfec} and
//' \code{usedstage} columns.
//' @param StageFrame The full stageframe for the analysis.
//' @param firstage The first true age to start the matrix with.
//' @param finalage The last true age to estimate.
//' @param cont A logical value indicating whether to lump survival past the
//' last age into a final age transition set on the supermatrix diagonal.
//' @param err_switch A logical value. If set to \code{TRUE}, then will also
//' output \code{probsrates} and \code{stage2fec}.
//' @param loypop A string vector giving the order of populations in the list of
//' years.
//' @param loypatch A string vector giving the order of patches in the list of
//' years.
//' @param loyyear2 A string vector giving the order of years at time \emph{t}
//' in the list of years.
//' @param yearorder The integer year order corresponding to \code{loyyear2}.
//' @param pop_var_int The variable number coding for population in the main
//' data set.
//' @param patch_var_int The variable number coding for patch in the main data
//' set.
//' @param year_var_int The variable number coding for year in time \emph{t} in
//' the main data set.
//' @param loy_pop_used A logical value indicating whether the population
//' variable is to be used.
//' @param loy_patch_used A logical value indicating whether the patch variable
//' is to be used.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' @param sparse If \code{TRUE}, then outputs matrices in sparse format.
//' Defaults to \code{FALSE}.
//' 
//' @return In the standard output, a list of three lists, called \code{A},
//' \code{U}, and \code{F}, each containing A, U, or F matrices, respectively.
//' Further information on outputs is provided in the documentation for
//' \code{specialpatrolgroup()}.
//' 
//' @keywords internal
//' @noRd
List subvertedpatrolgroup(const DataFrame& sge3, const arma::ivec& sge2index21,
  const DataFrame& MainData, const DataFrame& StageFrame, int firstage, int finalage,
  bool cont, int err_switch, StringVector loypop, StringVector loypatch,
  StringVector loyyear2, IntegerVector yearorder, const int pop_var_int,
  const int patch_var_int, const int year_var_int, const bool loy_pop_used,
  const bool loy_patch_used, bool simplicity = false, bool sparse = false) {
  
  arma::vec sge3fec32 = as<arma::vec>(sge3["repentry3"]);
  arma::uvec sge3rep2 = as<arma::uvec>(sge3["rep2n"]);
  arma::vec sge3ovgivent = as<arma::vec>(sge3["ovgiven_t"]);
  arma::vec sge3ovgivenf = as<arma::vec>(sge3["ovgiven_f"]);
  arma::vec sge3ovoffsett = as<arma::vec>(sge3["ovoffset_t"]);
  arma::vec sge3ovoffsetf = as<arma::vec>(sge3["ovoffset_f"]);
  arma::ivec sge3ovestt = as<arma::ivec>(sge3["ovest_t"]);
  arma::ivec sge3ovestf = as<arma::ivec>(sge3["ovest_f"]);
  arma::vec sge3ovsurvmult = as<arma::vec>(sge3["ovsurvmult"]);
  arma::vec sge3ovfecmult = as<arma::vec>(sge3["ovfecmult"]);
  arma::ivec sge3index321 = as<arma::ivec>(sge3["index321"]);
  arma::ivec sge3index21 = as<arma::ivec>(sge3["index21"]);
  arma::ivec sge3index2 = as<arma::ivec>(sge3["stage2n"]);
  arma::ivec aliveandequal = as<arma::ivec>(sge3["aliveandequal"]);
  
  arma::ivec dataindex321 = as<arma::ivec>(MainData["index321"]);
  arma::ivec dataindex21 = as<arma::ivec>(MainData["index21"]);
  arma::uvec dataalive3 = as<arma::uvec>(MainData["alive3"]);
  arma::vec datausedfec = as<arma::vec>(MainData["usedfec"]);
  
  int totalages = finalage - firstage + 1;
  int nostages = static_cast<int>(StageFrame.nrows());
  int no21stages = static_cast<int>(sge2index21.n_elem); // Includes Dead stage in every age
  int noelems = static_cast<int>(sge3index321.n_elem);
  bool small_subset {false};
  
  int loy_length = yearorder.length();
  List A_output (loy_length);
  List U_output (loy_length);
  List F_output (loy_length);
  List conc_err (loy_length);
  List s2f_err (loy_length);
  
  StringVector data_pop_;
  StringVector data_patch_;
  
  if (loy_pop_used) data_pop_ = as<StringVector>(MainData[pop_var_int]);
  if (loy_patch_used) data_patch_ = as<StringVector>(MainData[patch_var_int]);
  StringVector data_year_ = as<StringVector>(MainData[year_var_int]);
  
  arma::uvec ovgiventind = find(sge3ovgivent != -1.);
  arma::uvec ovgivenfind = find(sge3ovgivenf != -1.);
  int ovgtn = static_cast<int>(ovgiventind.n_elem);
  int ovgfn = static_cast<int>(ovgivenfind.n_elem);
  
  arma::uvec ovoffsettind = find(sge3ovoffsett != 0.);
  arma::uvec ovoffsetfind = find(sge3ovoffsetf != 0.);
  int ovotn = static_cast<int>(ovoffsettind.n_elem);
  int ovofn = static_cast<int>(ovoffsetfind.n_elem);
  
  arma::uvec ovesttind = find(sge3ovestt != -1);
  arma::uvec ovestfind = find(sge3ovestf != -1);
  int ovestn = static_cast<int>(ovesttind.n_elem);
  int ovesfn = static_cast<int>(ovestfind.n_elem);
  
  for (int i = 0; i < loy_length; i++) {
    arma::uvec data_current_index = as<arma::uvec>(LefkoUtils::index_l3(data_year_, loyyear2(i)));
    
    if (loy_pop_used) {
      arma::uvec data_current_pop = as<arma::uvec>(LefkoUtils::index_l3(data_pop_, loypop(i)));
      data_current_index = intersect(data_current_index, data_current_pop);
    }
    
    if (loy_patch_used) {
      arma::uvec data_current_patch = as<arma::uvec>(LefkoUtils::index_l3(data_patch_, loypatch(i)));
      data_current_index = intersect(data_current_index, data_current_patch);
    }
    
    int data_subset_rows = static_cast<int>(data_current_index.n_elem);
    if (data_subset_rows < 10 && !small_subset) {
      small_subset = true;
      Rf_warningcall(R_NilValue,
        "Extremely small data subsets used to populate matrices.");
    }
    if (data_subset_rows == 0) continue;
  
    arma::ivec dataindex321i = dataindex321.elem(data_current_index);
    arma::ivec dataindex21i = dataindex21.elem(data_current_index);
    arma::uvec dataalive3i = dataalive3.elem(data_current_index);
    arma::vec datausedfeci = datausedfec.elem(data_current_index);
    
    dataindex321i.resize(data_subset_rows);
    dataindex21i.resize(data_subset_rows);
    dataalive3i.resize(data_subset_rows);
    datausedfeci.resize(data_subset_rows);
    
    int n = static_cast<int>(dataindex321i.n_elem);
    unsigned int the_chosen_bun {0};
    
    arma::vec probsrates0(noelems, fill::zeros); // # indivs (3 trans)
    arma::vec probsrates1(noelems, fill::zeros); // total indivs for pair stage
    arma::vec probsrates2(noelems, fill::zeros); // total indivs alive for pair stage
    arma::vec probsrates3(noelems, fill::zeros); // total fec for pair stage
    
    // col0 = #indivs total, col2 = #indivs alive, col3 = sum fec
    arma::mat stage21fec(no21stages, 3, fill::zeros);
    
    arma::mat tmatrix;
    arma::mat fmatrix;
    arma::sp_mat tmatrix_sp;
    arma::sp_mat fmatrix_sp;
    
    if (!sparse) {
      arma::mat tmatrix_chuck(((nostages-1) * totalages), ((nostages-1) * totalages), fill::zeros);
      arma::mat fmatrix_chuck(((nostages-1) * totalages), ((nostages-1) * totalages), fill::zeros);
      
      tmatrix = tmatrix_chuck;
      fmatrix = fmatrix_chuck;
    } else {
      arma::sp_mat tmatrix_chuck(((nostages-1) * totalages), ((nostages-1) * totalages));
      arma::sp_mat fmatrix_chuck(((nostages-1) * totalages), ((nostages-1) * totalages));
      
      tmatrix_sp = tmatrix_chuck;
      fmatrix_sp = fmatrix_chuck;
    }
    
    // Count individuals going through transitions, sum their fecundities,
    // add that info to the 3-trans and 2-trans tables
    for (int j = 0; j < n; j++) { 
      // Sum all individuals with particular transition
      arma::uvec chosen_index321_vec = find(sge3index321 == dataindex321i(j));
      
      if (chosen_index321_vec.n_elem > 0) {
        int chosen_index321 = chosen_index321_vec(0);
        probsrates0(chosen_index321) = probsrates0(chosen_index321) + 1.0; 
      }
      
      // Sum all individuals with particular transition
      arma::uvec chosen_index21_vec = find(sge2index21 == dataindex21i(j));
      if (chosen_index21_vec.n_elem > 0) {
        int chosen_index21 = chosen_index21_vec(0);
        
        stage21fec(chosen_index21, 0) = stage21fec(chosen_index21, 0) + 1.0; 
        if (dataalive3i(j) > 0) {
          stage21fec(chosen_index21, 1) = stage21fec(chosen_index21, 1) + 1.0;
        }
        
        stage21fec(chosen_index21, 2) = stage21fec(chosen_index21, 2) + datausedfeci(j);
      }
    }
    
    // Populate vectors of individuals by stage in time t
    for (int j = 0; j < noelems; j++) {
      arma::uvec classy_aks = find(sge2index21 == sge3index21(j));
      
      if (classy_aks.n_elem > 0) {
        the_chosen_bun = classy_aks(0);
        
        probsrates1(j) = stage21fec(the_chosen_bun, 0);
        probsrates2(j) = stage21fec(the_chosen_bun, 1);
        probsrates3(j) = stage21fec(the_chosen_bun, 2);
      }
    }
    
    // Populate main U and F matrices
    for (int elem3 = 0; elem3 < noelems; elem3++) {
      if (aliveandequal(elem3) != -1) {
        
        // Leave NaNs in matrices when 0 individuals are summed in probsrates1
        if (sge3ovsurvmult(elem3) < 0.0) sge3ovsurvmult(elem3) = 1.0;
        if (!sparse) {
          tmatrix(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
            probsrates1(elem3);
        } else {
          tmatrix_sp(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
            probsrates1(elem3);
        }
          
        if (sge3ovfecmult(elem3) < 0.0) sge3ovfecmult(elem3) = 1.0;
        if (!sparse) {
          fmatrix(aliveandequal(elem3)) = static_cast<double>(sge3fec32(elem3)) *   
            static_cast<double>(sge3rep2(elem3)) * probsrates3(elem3) *
            sge3ovfecmult(elem3) / probsrates1(elem3);
        } else {
          fmatrix_sp(aliveandequal(elem3)) = static_cast<double>(sge3fec32(elem3)) *   
            static_cast<double>(sge3rep2(elem3)) * probsrates3(elem3) *
            sge3ovfecmult(elem3) / probsrates1(elem3);
        }
      }
    }
    
    // Correct for transitions given in overwrite table
    if (ovgtn > 0) {
      for (int j = 0; j < ovgtn; j++) {
        int matrixelement2 = aliveandequal(ovgiventind(j));
        if (matrixelement2 != -1) {
          if (!sparse) {
            tmatrix(matrixelement2) = sge3ovgivent(ovgiventind(j));
          } else {
            tmatrix_sp(matrixelement2) = sge3ovgivent(ovgiventind(j));
          }
        }
      }
    }
    
    if (ovgfn > 0) {
      for (int j = 0; j < ovgfn; j++) {
        int matrixelement2 = aliveandequal(ovgivenfind(j));
        if (matrixelement2 != -1) {
          if (!sparse) {
            fmatrix(matrixelement2) = sge3ovgivenf(ovgivenfind(j));
          } else {
            fmatrix_sp(matrixelement2) = sge3ovgivenf(ovgivenfind(j));
          }
        }
      }
    }
    
    // Correct for additive offsets given in overwrite table
    if (ovotn > 0) {
      for (int j = 0; j < ovotn; j++) {
        int matrixelement2 = aliveandequal(ovoffsettind(j));
        if (NumericVector::is_na(tmatrix(matrixelement2)) ||
          Rcpp::traits::is_nan<REALSXP>(tmatrix(matrixelement2))) {tmatrix(matrixelement2) = 0.;}
        if (matrixelement2 != -1) {
          if (!sparse) {
            tmatrix(matrixelement2) = tmatrix(matrixelement2) + sge3ovoffsett(ovoffsettind(j));
          } else {
            tmatrix_sp(matrixelement2) = tmatrix_sp(matrixelement2) + sge3ovoffsett(ovoffsettind(j));
          }
        }
      }
    }
    
    if (ovofn > 0) {
      for (int j = 0; j < ovofn; j++) {
        int matrixelement2 = aliveandequal(ovoffsetfind(j));
        if (NumericVector::is_na(fmatrix(matrixelement2)) ||
          Rcpp::traits::is_nan<REALSXP>(fmatrix(matrixelement2))) {fmatrix(matrixelement2) = 0.;}
        if (matrixelement2 != -1) {
          if (!sparse) {
            fmatrix(matrixelement2) = fmatrix(matrixelement2) + sge3ovoffsetf(ovoffsetfind(j));
          } else {
            fmatrix_sp(matrixelement2) = fmatrix_sp(matrixelement2) + sge3ovoffsetf(ovoffsetfind(j));
          }
        }
      }
    }
    
    // Replace transitions with proxies as given in overwrite table
    if (ovestn > 0) {
      for (int j = 0; j < ovestn; j++) {
        arma::uvec replacement = find(sge3index321 == sge3ovestt(ovesttind(j)));
        if (replacement.n_elem > 0) {
          if (aliveandequal(ovesttind(j)) != -1 && aliveandequal(replacement(0)) != -1) {
            if (!sparse) {
              tmatrix(aliveandequal(ovesttind(j))) = tmatrix(aliveandequal(replacement(0)));
            } else {
              tmatrix_sp(aliveandequal(ovesttind(j))) = tmatrix_sp(aliveandequal(replacement(0)));
            }
          }
        }
      }
    }
    
    if (ovesfn > 0) {
      for (int j = 0; j < ovesfn; j++) {
        arma::uvec replacement = find(sge3index321 == sge3ovestf(ovestfind(j)));
        if (replacement.n_elem > 0) {
          if (aliveandequal(ovestfind(j)) != -1 && aliveandequal(replacement(0)) != -1) {
            if (!sparse) {
              fmatrix(aliveandequal(ovestfind(j))) = fmatrix(aliveandequal(replacement(0)));
            } else {
              fmatrix_sp(aliveandequal(ovestfind(j))) = fmatrix_sp(aliveandequal(replacement(0)));
            }
          }
        }
      }
    }
    
    if (!sparse) {
      tmatrix(find_nonfinite(tmatrix)).zeros();
      fmatrix(find_nonfinite(fmatrix)).zeros();
      
      U_output(i) = tmatrix;
      F_output(i) = fmatrix;
      
      if (!simplicity) {
        arma::mat amatrix = tmatrix + fmatrix; // Create A matrix
        A_output(i) = amatrix;
      }
    } else {
      arma::uvec tm_sp_nan = find_nonfinite(tmatrix_sp);
      arma::uvec fm_sp_nan = find_nonfinite(fmatrix_sp);
      
      for (int nan_count = 0; nan_count < static_cast<int>(tm_sp_nan.n_elem); nan_count++) { 
        tmatrix_sp(tm_sp_nan(nan_count)) = 0.0;
      }
      for (int nan_count = 0; nan_count < static_cast<int>(fm_sp_nan.n_elem); nan_count++) { 
        fmatrix_sp(fm_sp_nan(nan_count)) = 0.0;
      }
      
      U_output(i) = tmatrix_sp;
      F_output(i) = fmatrix_sp;
      
      if (!simplicity) {
        arma::sp_mat amatrix_sp = tmatrix_sp + fmatrix_sp; // Create A matrix
        A_output(i) = amatrix_sp;
      }
    }
    
    if (err_switch) {
      DataFrame concatenated_crap = DataFrame::create(_["sge3index321"] = sge3index321,
        _["aliveandequal"] = aliveandequal, _["probsrates0"] = probsrates0,
        _["probsrates1"] = probsrates1, _["probsrates2"] = probsrates2,
        _["probsrates3"] = probsrates3, _["sge3fec32"] = sge3fec32);
      
      conc_err(i) = concatenated_crap;
      s2f_err(i) = stage21fec;
    }
  }
  
  List final_output;
  
  if (err_switch) {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue,
      _["probsrates_all"] = conc_err, _["s2f"] = s2f_err, _["data"] = MainData);
    final_output = out_dude;
    
  } else {
    List out_dude = List::create(Named("A") = A_output, _["U"] = U_output,
      _["F"] = F_output, _["ahstages"] = R_NilValue, _["agestages"] = R_NilValue,
      _["hstages"] = R_NilValue, _["labels"] = R_NilValue,
      _["matrixqc"] = R_NilValue, _["dataqc"] = R_NilValue);
    final_output = out_dude;
  }
  
  return final_output;
}

//' Key Function Passing Models and Other Parameters to Matrix Estimators
//' 
//' Function \code{raymccooney()} takes the various vital rate models and other
//' parameters and coordinates them as input into the function-based matrix
//' estimation functions.
//' 
//' @name raymccooney
//' 
//' @param listofyears A data frame where the rows designate the exact order of
//' years and patches to produce matrices for.
//' @param modelsuite An object of class \code{lefkoMod}, a similarly structured
//' list object, or a \code{vrm_input} object. All 14 vital rate models and the
//' \code{paramnames} data frame are required if not using a \code{vrm_input}
//' object.
//' @param mainyears A string vector of all times at time \emph{t}.
//' @param mainpatches A string vector of patch names.
//' @param maingroups A string vector of stage group names.
//' @param mainindcova Typically a string vector of individual covariate
//' category names.
//' @param mainindcovb Typically a string vector of individual covariate
//' category names.
//' @param mainindcovc Typically a string vector of individual covariate
//' category names.
//' @param StageFrame The stageframe object identifying the life history model
//' being operationalized.
//' @param OverWrite The supplement or overwrite table used in analysis, as
//' modified by \code{.sf_reassess()} (which technically uses function
//' \code{supp_reassess()}). Must be processed via \code{.supp_reassess} rather
//' than being a raw overwrite or supplement table.
//' @param repmatrix The reproductive matrix used in analysis.
//' @param f2_inda_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_annua_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_annua_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_annub_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_annub_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_annuc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_annuc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_inda_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{a} at each time \emph{t} to
//' be used in analysis.
//' @param f1_inda_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{a} at each time \emph{t-1}
//' to be used in analysis.
//' @param f2_indb_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{b} at each time \emph{t} to
//' be used in analysis.
//' @param f1_indb_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{b} at each time \emph{t-1}
//' to be used in analysis.
//' @param f2_indc_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{c} at each time \emph{t} to
//' be used in analysis.
//' @param f1_indc_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{c} at each time \emph{t-1}
//' to be used in analysis.
//' @param r2_inda A string vector of length equal to the number of years,
//' holding categories of individual random covariate \code{a} at each time
//' \emph{t} to be used in analysis.
//' @param r1_inda A string vector of length equal to the number of years,
//' holding categories of individual random covariate \code{a} at each time
//' \emph{t-1} to be used in analysis.
//' @param r2_indb A string vector of length equal to the number of years,
//' holding categories of individual random covariate \code{b} at each time
//' \emph{t} to be used in analysis.
//' @param r1_indb A string vector of length equal to the number of years,
//' holding categories of individual random covariate \code{b} at each time
//' \emph{t-1} to be used in analysis.
//' @param r2_indc A string vector of length equal to the number of years,
//' holding categories of individual random covariate \code{c} at each time
//' \emph{t} to be used in analysis.
//' @param r1_indc A string vector of length equal to the number of years,
//' holding categories of individual random covariate \code{c} at each time
//' \emph{t-1} to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' juvenile reproductive status, and juvenile maturity status.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param firstage The first age to be used in the analysis. Should typically
//' be \code{0} for pre-breeding and \code{1} for post-breeding life history
//' models. If not building age-by-stage MPMs, then should be set to \code{0}.
//' @param finalage The final age to be used in analysis. If not building
//' age-by-stage MPMs, then should be set to \code{0}.
//' @param format Indicates whether historical matrices should be in (1) Ehrlen
//' or (2) deVries format.
//' @param style The style of analysis, where 0 is historical, 1 is ahistorical,
//' and 2 is age-by-stage.
//' @param cont Denotes whether age-by-stage matrix continues past the final
//' age.
//' @param filter An integer denoting whether to filter the DataFrame to
//' eliminate unusable rows, and if so, how to do so. Possible values: \code{0}:
//' no filtering, \code{1}: filter out rows with \code{index321 == -1}, and
//' \code{2}: filter out rows with \code{aliveandequal == -1}.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0. Defaults to \code{FALSE}.
//' @param nodata A logical value indicating whether the modelsuite contains
//' all parameter coefficients and no hfv dataset is provided (\code{TRUE}), or
//' whether an hfv dataset and a true modelsuite are provided (\code{FALSE}).
//' Defaults to \code{FALSE}.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param cdf A logical value indicating whether to estimate size transitions
//' using the cumulative density function in cases with continuous
//' distributions. Defaults to \code{TRUE}, with the midpoint method used if
//' \code{FALSE}.
//' @param err_check If \code{TRUE}, then also output objects \code{prob_out}
//' and \code{allstages} for error checking purposes.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' @param sparse If \code{TRUE}, then outputs matrices in sparse format.
//' Defaults to \code{FALSE}.
//' 
//' @return A list with with up to 6 elements. In order: \code{A}: a list of A
//' matrices, or a list of \code{NULL} values if \code{simplicity = TRUE};
//' \code{U}: a list of U matrices, in the same order as \code{A}; \code{F}:
//' a list of F matrices, in the same order as \code{A}; \code{prob_out}: a list
//' of error-checking conditional probability matrices, or a list of \code{NULL}
//' values if \code{err_check = FALSE}; \code{allstages}: a data frame showing
//' the used values of all variables used in transition calculations; and
//' \code{proxies}: a list of model proxies output during function-based model
//' processing.
//' 
//' @keywords internal
//' @noRd
List raymccooney(const DataFrame& listofyears, const List& modelsuite,
  const CharacterVector& mainyears, const CharacterVector& mainpatches,
  RObject maingroups, RObject mainindcova, RObject mainindcovb,
  RObject mainindcovc, const DataFrame& StageFrame, const DataFrame& OverWrite,
  const arma::mat& repmatrix, NumericVector f2_inda_num,
  NumericVector f1_inda_num, NumericVector f2_indb_num,
  NumericVector f1_indb_num, NumericVector f2_indc_num,
  NumericVector f1_indc_num, NumericVector f2_annua_num,
  NumericVector f1_annua_num, NumericVector f2_annub_num,
  NumericVector f1_annub_num, NumericVector f2_annuc_num,
  NumericVector f1_annuc_num, StringVector f2_inda_cat,
  StringVector f1_inda_cat, StringVector f2_indb_cat, StringVector f1_indb_cat,
  StringVector f2_indc_cat, StringVector f1_indc_cat, StringVector r2_inda,
  StringVector r1_inda, StringVector r2_indb, StringVector r1_indb,
  StringVector r2_indc, StringVector r1_indc, const NumericVector& dev_terms,
  double dens, double fecmod, int firstage, int finalage, int format, int style,
  int cont, int filter, bool negfec = false, bool nodata = false,
  double exp_tol = 700.0, double theta_tol = 1e8, bool cdf = true,
  bool err_check = false, bool simplicity = false, bool sparse = false) {
  
  // Dud dens_vr inputs
  Rcpp::DataFrame dvr_frame;
  List err_check_proxies;
  
  // listofyears import and settings
  IntegerVector years = listofyears["yearorder"];
  IntegerVector patches = listofyears["patchorder"];
  int loy_length = years.length();
  
  int matrixformat {0};
  if (style == 2) {
    matrixformat = 4;
  } else if (style == 1) {
    matrixformat = 3;
  } else if (style == 0) {
    if (format == 1) {
      matrixformat = 1;
    } else if (format == 2) {
      matrixformat = 2;
    } else pop_error("format", "", "", 5);
  } else pop_error("style", "", "", 5);
  
  Rcpp::DataFrame allstages = theoldpizzle(StageFrame, OverWrite, repmatrix,
    firstage, finalage, format, style, cont, filter);
  
  NumericVector size3 = allstages["size3"];
  NumericVector size2n = allstages["size2n"];
  NumericVector size2o = allstages["size2o"];
  NumericVector sizeb3 = allstages["sizeb3"];
  NumericVector sizeb2n = allstages["sizeb2n"];
  NumericVector sizeb2o = allstages["sizeb2o"];
  NumericVector sizec3 = allstages["sizec3"];
  NumericVector sizec2n = allstages["sizec2n"];
  NumericVector sizec2o = allstages["sizec2o"];
  
  NumericVector maxveca = {max(size3), max(size2n), max(size2o)}; // What about NAs?
  NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
  NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
  
  double maxsize = max(maxveca); // What about NAs?
  double maxsizeb = max(maxvecb);
  double maxsizec = max(maxvecc);
  
  // Move model summaries to appropriate RObjects
  RObject surv_model;
  RObject obs_model;
  RObject size_model;
  RObject sizeb_model;
  RObject sizec_model;
  RObject repst_model;
  RObject fec_model;
  RObject jsurv_model;
  RObject jobs_model;
  RObject jsize_model;
  RObject jsizeb_model;
  RObject jsizec_model;
  RObject jrepst_model;
  RObject jmatst_model;
  DataFrame paramnames;
  
  if (nodata) {
    DataFrame vrm_frame = as<DataFrame>(modelsuite["vrm_frame"]);
    DataFrame year_frame = as<DataFrame>(modelsuite["year_frame"]);
    DataFrame patch_frame = as<DataFrame>(modelsuite["patch_frame"]);
    DataFrame group2_frame = as<DataFrame>(modelsuite["group2_frame"]);
    DataFrame group1_frame = as<DataFrame>(modelsuite["group1_frame"]);
    DataFrame dist_frame = as<DataFrame>(modelsuite["dist_frame"]);
    NumericVector st_frame = as<NumericVector>(modelsuite["st_frame"]);
    
    CharacterVector main_effect_1 = as<CharacterVector>(vrm_frame["main_effect_1"]);
    CharacterVector effects_names = clone(main_effect_1);
    
    CharacterVector main_effect_2;
    if (main_effect_1.length() > 23) {
      main_effect_2 = as<CharacterVector>(vrm_frame["main_effect_2"]);
      
      for (int i = 0; i < main_effect_1.length(); i++) {
        if (i > 16) {
          effects_names(i) += ":";
          effects_names(i) += main_effect_2(i);
        }
      }
    }
    
    CharacterVector year_names = as<CharacterVector>(year_frame["years"]);
    CharacterVector patch_names = as<CharacterVector>(patch_frame["patches"]);
    CharacterVector group_names = as<CharacterVector>(group2_frame["groups"]);
    
    bool zi_yn = false;
    int vrm_length = vrm_frame.length();
    
    NumericVector surv_num = as<NumericVector>(vrm_frame["surv"]);
    NumericVector obs_num = as<NumericVector>(vrm_frame["obs"]);
    NumericVector sizea_num = as<NumericVector>(vrm_frame["sizea"]);
    NumericVector sizeb_num = as<NumericVector>(vrm_frame["sizeb"]);
    NumericVector sizec_num = as<NumericVector>(vrm_frame["sizec"]);
    NumericVector repst_num = as<NumericVector>(vrm_frame["repst"]);
    NumericVector fec_num = as<NumericVector>(vrm_frame["fec"]);
    NumericVector jsurv_num = as<NumericVector>(vrm_frame["jsurv"]);
    NumericVector jobs_num = as<NumericVector>(vrm_frame["jobs"]);
    NumericVector jsizea_num = as<NumericVector>(vrm_frame["jsizea"]);
    NumericVector jsizeb_num = as<NumericVector>(vrm_frame["jsizeb"]);
    NumericVector jsizec_num = as<NumericVector>(vrm_frame["jsizec"]);
    NumericVector jrepst_num = as<NumericVector>(vrm_frame["jrepst"]);
    NumericVector jmatst_num = as<NumericVector>(vrm_frame["jmatst"]);
    
    NumericVector surv_year = as<NumericVector>(year_frame["surv"]);
    NumericVector obs_year = as<NumericVector>(year_frame["obs"]);
    NumericVector sizea_year = as<NumericVector>(year_frame["sizea"]);
    NumericVector sizeb_year = as<NumericVector>(year_frame["sizeb"]);
    NumericVector sizec_year = as<NumericVector>(year_frame["sizec"]);
    NumericVector repst_year = as<NumericVector>(year_frame["repst"]);
    NumericVector fec_year = as<NumericVector>(year_frame["fec"]);
    NumericVector jsurv_year = as<NumericVector>(year_frame["jsurv"]);
    NumericVector jobs_year = as<NumericVector>(year_frame["jobs"]);
    NumericVector jsizea_year = as<NumericVector>(year_frame["jsizea"]);
    NumericVector jsizeb_year = as<NumericVector>(year_frame["jsizeb"]);
    NumericVector jsizec_year = as<NumericVector>(year_frame["jsizec"]);
    NumericVector jrepst_year = as<NumericVector>(year_frame["jrepst"]);
    NumericVector jmatst_year = as<NumericVector>(year_frame["jmatst"]);
    
    NumericVector surv_patch = as<NumericVector>(patch_frame["surv"]);
    NumericVector obs_patch = as<NumericVector>(patch_frame["obs"]);
    NumericVector sizea_patch = as<NumericVector>(patch_frame["sizea"]);
    NumericVector sizeb_patch = as<NumericVector>(patch_frame["sizeb"]);
    NumericVector sizec_patch = as<NumericVector>(patch_frame["sizec"]);
    NumericVector repst_patch = as<NumericVector>(patch_frame["repst"]);
    NumericVector fec_patch = as<NumericVector>(patch_frame["fec"]);
    NumericVector jsurv_patch = as<NumericVector>(patch_frame["jsurv"]);
    NumericVector jobs_patch = as<NumericVector>(patch_frame["jobs"]);
    NumericVector jsizea_patch = as<NumericVector>(patch_frame["jsizea"]);
    NumericVector jsizeb_patch = as<NumericVector>(patch_frame["jsizeb"]);
    NumericVector jsizec_patch = as<NumericVector>(patch_frame["jsizec"]);
    NumericVector jrepst_patch = as<NumericVector>(patch_frame["jrepst"]);
    NumericVector jmatst_patch = as<NumericVector>(patch_frame["jmatst"]);
    
    NumericVector surv_group2 = as<NumericVector>(group2_frame["surv"]);
    NumericVector obs_group2 = as<NumericVector>(group2_frame["obs"]);
    NumericVector sizea_group2 = as<NumericVector>(group2_frame["sizea"]);
    NumericVector sizeb_group2 = as<NumericVector>(group2_frame["sizeb"]);
    NumericVector sizec_group2 = as<NumericVector>(group2_frame["sizec"]);
    NumericVector repst_group2 = as<NumericVector>(group2_frame["repst"]);
    NumericVector fec_group2 = as<NumericVector>(group2_frame["fec"]);
    NumericVector jsurv_group2 = as<NumericVector>(group2_frame["jsurv"]);
    NumericVector jobs_group2 = as<NumericVector>(group2_frame["jobs"]);
    NumericVector jsizea_group2 = as<NumericVector>(group2_frame["jsizea"]);
    NumericVector jsizeb_group2 = as<NumericVector>(group2_frame["jsizeb"]);
    NumericVector jsizec_group2 = as<NumericVector>(group2_frame["jsizec"]);
    NumericVector jrepst_group2 = as<NumericVector>(group2_frame["jrepst"]);
    NumericVector jmatst_group2 = as<NumericVector>(group2_frame["jmatst"]);
    
    NumericVector surv_group1 = as<NumericVector>(group1_frame["surv"]);
    NumericVector obs_group1 = as<NumericVector>(group1_frame["obs"]);
    NumericVector sizea_group1 = as<NumericVector>(group1_frame["sizea"]);
    NumericVector sizeb_group1 = as<NumericVector>(group1_frame["sizeb"]);
    NumericVector sizec_group1 = as<NumericVector>(group1_frame["sizec"]);
    NumericVector repst_group1 = as<NumericVector>(group1_frame["repst"]);
    NumericVector fec_group1 = as<NumericVector>(group1_frame["fec"]);
    NumericVector jsurv_group1 = as<NumericVector>(group1_frame["jsurv"]);
    NumericVector jobs_group1 = as<NumericVector>(group1_frame["jobs"]);
    NumericVector jsizea_group1 = as<NumericVector>(group1_frame["jsizea"]);
    NumericVector jsizeb_group1 = as<NumericVector>(group1_frame["jsizeb"]);
    NumericVector jsizec_group1 = as<NumericVector>(group1_frame["jsizec"]);
    NumericVector jrepst_group1 = as<NumericVector>(group1_frame["jrepst"]);
    NumericVector jmatst_group1 = as<NumericVector>(group1_frame["jmatst"]);
    
    StringVector distribs = as<StringVector>(dist_frame["dist"]);
    String surv_dist = distribs(0);
    String obs_dist = distribs(1);
    String sizea_dist = distribs(2);
    String sizeb_dist = distribs(3);
    String sizec_dist = distribs(4);
    String repst_dist = distribs(5);
    String fec_dist = distribs(6);
    String jsurv_dist = distribs(7);
    String jobs_dist = distribs(8);
    String jsizea_dist = distribs(9);
    String jsizeb_dist = distribs(10);
    String jsizec_dist = distribs(11);
    String jrepst_dist = distribs(12);
    String jmatst_dist = distribs(13);
    
    double sizea_st = st_frame(2);
    double sizeb_st = st_frame(3);
    double sizec_st = st_frame(4);
    double fec_st = st_frame(6);
    double jsizea_st = st_frame(9);
    double jsizeb_st = st_frame(10);
    double jsizec_st = st_frame(11);
    
    NumericVector sizea_zi;
    NumericVector sizeb_zi;
    NumericVector sizec_zi;
    NumericVector fec_zi;
    NumericVector jsizea_zi;
    NumericVector jsizeb_zi;
    NumericVector jsizec_zi;
    
    NumericVector year_sizea_zi;
    NumericVector year_sizeb_zi;
    NumericVector year_sizec_zi;
    NumericVector year_fec_zi;
    NumericVector year_jsizea_zi;
    NumericVector year_jsizeb_zi;
    NumericVector year_jsizec_zi;
    
    NumericVector patch_sizea_zi;
    NumericVector patch_sizeb_zi;
    NumericVector patch_sizec_zi;
    NumericVector patch_fec_zi;
    NumericVector patch_jsizea_zi;
    NumericVector patch_jsizeb_zi;
    NumericVector patch_jsizec_zi;
    
    NumericVector group2_sizea_zi;
    NumericVector group2_sizeb_zi;
    NumericVector group2_sizec_zi;
    NumericVector group2_fec_zi;
    NumericVector group2_jsizea_zi;
    NumericVector group2_jsizeb_zi;
    NumericVector group2_jsizec_zi;
    
    NumericVector group1_sizea_zi;
    NumericVector group1_sizeb_zi;
    NumericVector group1_sizec_zi;
    NumericVector group1_fec_zi;
    NumericVector group1_jsizea_zi;
    NumericVector group1_jsizeb_zi;
    NumericVector group1_jsizec_zi;
    
    NumericVector dud_zi;
    
    if (vrm_length > 16) {
      zi_yn = true;
      
      sizea_zi = as<NumericVector>(vrm_frame["sizea_zi"]);
      sizeb_zi = as<NumericVector>(vrm_frame["sizeb_zi"]);
      sizec_zi = as<NumericVector>(vrm_frame["sizec_zi"]);
      fec_zi = as<NumericVector>(vrm_frame["fec_zi"]);
      jsizea_zi = as<NumericVector>(vrm_frame["jsizea_zi"]);
      jsizeb_zi = as<NumericVector>(vrm_frame["jsizeb_zi"]);
      jsizec_zi = as<NumericVector>(vrm_frame["jsizec_zi"]);
      
      year_sizea_zi = as<NumericVector>(year_frame["sizea_zi"]);
      year_sizeb_zi = as<NumericVector>(year_frame["sizeb_zi"]);
      year_sizec_zi = as<NumericVector>(year_frame["sizec_zi"]);
      year_fec_zi = as<NumericVector>(year_frame["fec_zi"]);
      year_jsizea_zi = as<NumericVector>(year_frame["jsizea_zi"]);
      year_jsizeb_zi = as<NumericVector>(year_frame["jsizeb_zi"]);
      year_jsizec_zi = as<NumericVector>(year_frame["jsizec_zi"]);
      
      patch_sizea_zi = as<NumericVector>(patch_frame["sizea_zi"]);
      patch_sizeb_zi = as<NumericVector>(patch_frame["sizeb_zi"]);
      patch_sizec_zi = as<NumericVector>(patch_frame["sizec_zi"]);
      patch_fec_zi = as<NumericVector>(patch_frame["fec_zi"]);
      patch_jsizea_zi = as<NumericVector>(patch_frame["jsizea_zi"]);
      patch_jsizeb_zi = as<NumericVector>(patch_frame["jsizeb_zi"]);
      patch_jsizec_zi = as<NumericVector>(patch_frame["jsizec_zi"]);
      
      group2_sizea_zi = as<NumericVector>(group2_frame["sizea_zi"]);
      group2_sizeb_zi = as<NumericVector>(group2_frame["sizeb_zi"]);
      group2_sizec_zi = as<NumericVector>(group2_frame["sizec_zi"]);
      group2_fec_zi = as<NumericVector>(group2_frame["fec_zi"]);
      group2_jsizea_zi = as<NumericVector>(group2_frame["jsizea_zi"]);
      group2_jsizeb_zi = as<NumericVector>(group2_frame["jsizeb_zi"]);
      group2_jsizec_zi = as<NumericVector>(group2_frame["jsizec_zi"]);
      
      group1_sizea_zi = as<NumericVector>(group1_frame["sizea_zi"]);
      group1_sizeb_zi = as<NumericVector>(group1_frame["sizeb_zi"]);
      group1_sizec_zi = as<NumericVector>(group1_frame["sizec_zi"]);
      group1_fec_zi = as<NumericVector>(group1_frame["fec_zi"]);
      group1_jsizea_zi = as<NumericVector>(group1_frame["jsizea_zi"]);
      group1_jsizeb_zi = as<NumericVector>(group1_frame["jsizeb_zi"]);
      group1_jsizec_zi = as<NumericVector>(group1_frame["jsizec_zi"]);
    }
    
    CharacterVector indcova_names;
    CharacterVector indcovb_names;
    CharacterVector indcovc_names;
    
    NumericVector surv_indcova2;
    NumericVector surv_indcovb2;
    NumericVector surv_indcovc2;
    NumericVector obs_indcova2;
    NumericVector obs_indcovb2;
    NumericVector obs_indcovc2;
    NumericVector sizea_indcova2;
    NumericVector sizea_indcovb2;
    NumericVector sizea_indcovc2;
    NumericVector sizeb_indcova2;
    NumericVector sizeb_indcovb2;
    NumericVector sizeb_indcovc2;
    NumericVector sizec_indcova2;
    NumericVector sizec_indcovb2;
    NumericVector sizec_indcovc2;
    NumericVector repst_indcova2;
    NumericVector repst_indcovb2;
    NumericVector repst_indcovc2;
    NumericVector fec_indcova2;
    NumericVector fec_indcovb2;
    NumericVector fec_indcovc2;
    NumericVector jsurv_indcova2;
    NumericVector jsurv_indcovb2;
    NumericVector jsurv_indcovc2;
    NumericVector jobs_indcova2;
    NumericVector jobs_indcovb2;
    NumericVector jobs_indcovc2;
    NumericVector jsizea_indcova2;
    NumericVector jsizea_indcovb2;
    NumericVector jsizea_indcovc2;
    NumericVector jsizeb_indcova2;
    NumericVector jsizeb_indcovb2;
    NumericVector jsizeb_indcovc2;
    NumericVector jsizec_indcova2;
    NumericVector jsizec_indcovb2;
    NumericVector jsizec_indcovc2;
    NumericVector jrepst_indcova2;
    NumericVector jrepst_indcovb2;
    NumericVector jrepst_indcovc2;
    NumericVector jmatst_indcova2;
    NumericVector jmatst_indcovb2;
    NumericVector jmatst_indcovc2;
    
    NumericVector sizea_indcova2_zi;
    NumericVector sizea_indcovb2_zi;
    NumericVector sizea_indcovc2_zi;
    NumericVector sizeb_indcova2_zi;
    NumericVector sizeb_indcovb2_zi;
    NumericVector sizeb_indcovc2_zi;
    NumericVector sizec_indcova2_zi;
    NumericVector sizec_indcovb2_zi;
    NumericVector sizec_indcovc2_zi;
    NumericVector fec_indcova2_zi;
    NumericVector fec_indcovb2_zi;
    NumericVector fec_indcovc2_zi;
    NumericVector jsizea_indcova2_zi;
    NumericVector jsizea_indcovb2_zi;
    NumericVector jsizea_indcovc2_zi;
    NumericVector jsizeb_indcova2_zi;
    NumericVector jsizeb_indcovb2_zi;
    NumericVector jsizeb_indcovc2_zi;
    NumericVector jsizec_indcova2_zi;
    NumericVector jsizec_indcovb2_zi;
    NumericVector jsizec_indcovc2_zi;
    
    NumericVector surv_indcova1;
    NumericVector surv_indcovb1;
    NumericVector surv_indcovc1;
    NumericVector obs_indcova1;
    NumericVector obs_indcovb1;
    NumericVector obs_indcovc1;
    NumericVector sizea_indcova1;
    NumericVector sizea_indcovb1;
    NumericVector sizea_indcovc1;
    NumericVector sizeb_indcova1;
    NumericVector sizeb_indcovb1;
    NumericVector sizeb_indcovc1;
    NumericVector sizec_indcova1;
    NumericVector sizec_indcovb1;
    NumericVector sizec_indcovc1;
    NumericVector repst_indcova1;
    NumericVector repst_indcovb1;
    NumericVector repst_indcovc1;
    NumericVector fec_indcova1;
    NumericVector fec_indcovb1;
    NumericVector fec_indcovc1;
    NumericVector jsurv_indcova1;
    NumericVector jsurv_indcovb1;
    NumericVector jsurv_indcovc1;
    NumericVector jobs_indcova1;
    NumericVector jobs_indcovb1;
    NumericVector jobs_indcovc1;
    NumericVector jsizea_indcova1;
    NumericVector jsizea_indcovb1;
    NumericVector jsizea_indcovc1;
    NumericVector jsizeb_indcova1;
    NumericVector jsizeb_indcovb1;
    NumericVector jsizeb_indcovc1;
    NumericVector jsizec_indcova1;
    NumericVector jsizec_indcovb1;
    NumericVector jsizec_indcovc1;
    NumericVector jrepst_indcova1;
    NumericVector jrepst_indcovb1;
    NumericVector jrepst_indcovc1;
    NumericVector jmatst_indcova1;
    NumericVector jmatst_indcovb1;
    NumericVector jmatst_indcovc1;
    
    NumericVector sizea_indcova1_zi;
    NumericVector sizea_indcovb1_zi;
    NumericVector sizea_indcovc1_zi;
    NumericVector sizeb_indcova1_zi;
    NumericVector sizeb_indcovb1_zi;
    NumericVector sizeb_indcovc1_zi;
    NumericVector sizec_indcova1_zi;
    NumericVector sizec_indcovb1_zi;
    NumericVector sizec_indcovc1_zi;
    NumericVector fec_indcova1_zi;
    NumericVector fec_indcovb1_zi;
    NumericVector fec_indcovc1_zi;
    NumericVector jsizea_indcova1_zi;
    NumericVector jsizea_indcovb1_zi;
    NumericVector jsizea_indcovc1_zi;
    NumericVector jsizeb_indcova1_zi;
    NumericVector jsizeb_indcovb1_zi;
    NumericVector jsizeb_indcovc1_zi;
    NumericVector jsizec_indcova1_zi;
    NumericVector jsizec_indcovb1_zi;
    NumericVector jsizec_indcovc1_zi;
    
    NumericVector surv_annucova2;
    NumericVector surv_annucovb2;
    NumericVector surv_annucovc2;
    NumericVector obs_annucova2;
    NumericVector obs_annucovb2;
    NumericVector obs_annucovc2;
    NumericVector sizea_annucova2;
    NumericVector sizea_annucovb2;
    NumericVector sizea_annucovc2;
    NumericVector sizeb_annucova2;
    NumericVector sizeb_annucovb2;
    NumericVector sizeb_annucovc2;
    NumericVector sizec_annucova2;
    NumericVector sizec_annucovb2;
    NumericVector sizec_annucovc2;
    NumericVector repst_annucova2;
    NumericVector repst_annucovb2;
    NumericVector repst_annucovc2;
    NumericVector fec_annucova2;
    NumericVector fec_annucovb2;
    NumericVector fec_annucovc2;
    NumericVector jsurv_annucova2;
    NumericVector jsurv_annucovb2;
    NumericVector jsurv_annucovc2;
    NumericVector jobs_annucova2;
    NumericVector jobs_annucovb2;
    NumericVector jobs_annucovc2;
    NumericVector jsizea_annucova2;
    NumericVector jsizea_annucovb2;
    NumericVector jsizea_annucovc2;
    NumericVector jsizeb_annucova2;
    NumericVector jsizeb_annucovb2;
    NumericVector jsizeb_annucovc2;
    NumericVector jsizec_annucova2;
    NumericVector jsizec_annucovb2;
    NumericVector jsizec_annucovc2;
    NumericVector jrepst_annucova2;
    NumericVector jrepst_annucovb2;
    NumericVector jrepst_annucovc2;
    NumericVector jmatst_annucova2;
    NumericVector jmatst_annucovb2;
    NumericVector jmatst_annucovc2;
    
    NumericVector sizea_annucova2_zi;
    NumericVector sizea_annucovb2_zi;
    NumericVector sizea_annucovc2_zi;
    NumericVector sizeb_annucova2_zi;
    NumericVector sizeb_annucovb2_zi;
    NumericVector sizeb_annucovc2_zi;
    NumericVector sizec_annucova2_zi;
    NumericVector sizec_annucovb2_zi;
    NumericVector sizec_annucovc2_zi;
    NumericVector fec_annucova2_zi;
    NumericVector fec_annucovb2_zi;
    NumericVector fec_annucovc2_zi;
    NumericVector jsizea_annucova2_zi;
    NumericVector jsizea_annucovb2_zi;
    NumericVector jsizea_annucovc2_zi;
    NumericVector jsizeb_annucova2_zi;
    NumericVector jsizeb_annucovb2_zi;
    NumericVector jsizeb_annucovc2_zi;
    NumericVector jsizec_annucova2_zi;
    NumericVector jsizec_annucovb2_zi;
    NumericVector jsizec_annucovc2_zi;
    
    NumericVector surv_annucova1;
    NumericVector surv_annucovb1;
    NumericVector surv_annucovc1;
    NumericVector obs_annucova1;
    NumericVector obs_annucovb1;
    NumericVector obs_annucovc1;
    NumericVector sizea_annucova1;
    NumericVector sizea_annucovb1;
    NumericVector sizea_annucovc1;
    NumericVector sizeb_annucova1;
    NumericVector sizeb_annucovb1;
    NumericVector sizeb_annucovc1;
    NumericVector sizec_annucova1;
    NumericVector sizec_annucovb1;
    NumericVector sizec_annucovc1;
    NumericVector repst_annucova1;
    NumericVector repst_annucovb1;
    NumericVector repst_annucovc1;
    NumericVector fec_annucova1;
    NumericVector fec_annucovb1;
    NumericVector fec_annucovc1;
    NumericVector jsurv_annucova1;
    NumericVector jsurv_annucovb1;
    NumericVector jsurv_annucovc1;
    NumericVector jobs_annucova1;
    NumericVector jobs_annucovb1;
    NumericVector jobs_annucovc1;
    NumericVector jsizea_annucova1;
    NumericVector jsizea_annucovb1;
    NumericVector jsizea_annucovc1;
    NumericVector jsizeb_annucova1;
    NumericVector jsizeb_annucovb1;
    NumericVector jsizeb_annucovc1;
    NumericVector jsizec_annucova1;
    NumericVector jsizec_annucovb1;
    NumericVector jsizec_annucovc1;
    NumericVector jrepst_annucova1;
    NumericVector jrepst_annucovb1;
    NumericVector jrepst_annucovc1;
    NumericVector jmatst_annucova1;
    NumericVector jmatst_annucovb1;
    NumericVector jmatst_annucovc1;
    
    NumericVector sizea_annucova1_zi;
    NumericVector sizea_annucovb1_zi;
    NumericVector sizea_annucovc1_zi;
    NumericVector sizeb_annucova1_zi;
    NumericVector sizeb_annucovb1_zi;
    NumericVector sizeb_annucovc1_zi;
    NumericVector sizec_annucova1_zi;
    NumericVector sizec_annucovb1_zi;
    NumericVector sizec_annucovc1_zi;
    NumericVector fec_annucova1_zi;
    NumericVector fec_annucovb1_zi;
    NumericVector fec_annucovc1_zi;
    NumericVector jsizea_annucova1_zi;
    NumericVector jsizea_annucovb1_zi;
    NumericVector jsizea_annucovc1_zi;
    NumericVector jsizeb_annucova1_zi;
    NumericVector jsizeb_annucovb1_zi;
    NumericVector jsizeb_annucovc1_zi;
    NumericVector jsizec_annucova1_zi;
    NumericVector jsizec_annucovb1_zi;
    NumericVector jsizec_annucovc1_zi;
    
    int modelsuite_length = modelsuite.length();
    CharacterVector modelsuite_names = modelsuite.attr("names");
    
    for (int i = 0; i < modelsuite_length; i++) {
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova2_frame")) {
        DataFrame indcova2_frame = as<DataFrame>(modelsuite["indcova2_frame"]);
        
        indcova_names = indcova2_frame["indcova"];
        
        surv_indcova2 = indcova2_frame["surv"];
        obs_indcova2 = indcova2_frame["obs"];
        sizea_indcova2 = indcova2_frame["sizea"];
        sizeb_indcova2 = indcova2_frame["sizeb"];
        sizec_indcova2 = indcova2_frame["sizec"];
        repst_indcova2 = indcova2_frame["repst"];
        fec_indcova2 = indcova2_frame["fec"];
        
        jsurv_indcova2 = indcova2_frame["jsurv"];
        jobs_indcova2 = indcova2_frame["jobs"];
        jsizea_indcova2 = indcova2_frame["jsizea"];
        jsizeb_indcova2 = indcova2_frame["jsizeb"];
        jsizec_indcova2 = indcova2_frame["jsizec"];
        jrepst_indcova2 = indcova2_frame["jrepst"];
        jmatst_indcova2 = indcova2_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcova2_zi = indcova2_frame["sizea_zi"];
          sizeb_indcova2_zi = indcova2_frame["sizeb_zi"];
          sizec_indcova2_zi = indcova2_frame["sizec_zi"];
          fec_indcova2_zi = indcova2_frame["fec_zi"];
          jsizea_indcova2_zi = indcova2_frame["jsizea_zi"];
          jsizeb_indcova2_zi = indcova2_frame["jsizeb_zi"];
          jsizec_indcova2_zi = indcova2_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova1_frame")) {
        DataFrame indcova1_frame = as<DataFrame>(modelsuite["indcova1_frame"]);
        
        indcova_names = indcova1_frame["indcova"];
        
        surv_indcova1 = indcova1_frame["surv"];
        obs_indcova1 = indcova1_frame["obs"];
        sizea_indcova1 = indcova1_frame["sizea"];
        sizeb_indcova1 = indcova1_frame["sizeb"];
        sizec_indcova1 = indcova1_frame["sizec"];
        repst_indcova1 = indcova1_frame["repst"];
        fec_indcova1 = indcova1_frame["fec"];
        
        jsurv_indcova1 = indcova1_frame["jsurv"];
        jobs_indcova1 = indcova1_frame["jobs"];
        jsizea_indcova1 = indcova1_frame["jsizea"];
        jsizeb_indcova1 = indcova1_frame["jsizeb"];
        jsizec_indcova1 = indcova1_frame["jsizec"];
        jrepst_indcova1 = indcova1_frame["jrepst"];
        jmatst_indcova1 = indcova1_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcova1_zi = indcova1_frame["sizea_zi"];
          sizeb_indcova1_zi = indcova1_frame["sizeb_zi"];
          sizec_indcova1_zi = indcova1_frame["sizec_zi"];
          fec_indcova1_zi = indcova1_frame["fec_zi"];
          jsizea_indcova1_zi = indcova1_frame["jsizea_zi"];
          jsizeb_indcova1_zi = indcova1_frame["jsizeb_zi"];
          jsizec_indcova1_zi = indcova1_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb2_frame")) {
        DataFrame indcovb2_frame = as<DataFrame>(modelsuite["indcovb2_frame"]);
        
        indcovb_names = indcovb2_frame["indcovb"];
        
        surv_indcovb2 = indcovb2_frame["surv"];
        obs_indcovb2 = indcovb2_frame["obs"];
        sizea_indcovb2 = indcovb2_frame["sizea"];
        sizeb_indcovb2 = indcovb2_frame["sizeb"];
        sizec_indcovb2 = indcovb2_frame["sizec"];
        repst_indcovb2 = indcovb2_frame["repst"];
        fec_indcovb2 = indcovb2_frame["fec"];
        
        jsurv_indcovb2 = indcovb2_frame["jsurv"];
        jobs_indcovb2 = indcovb2_frame["jobs"];
        jsizea_indcovb2 = indcovb2_frame["jsizea"];
        jsizeb_indcovb2 = indcovb2_frame["jsizeb"];
        jsizec_indcovb2 = indcovb2_frame["jsizec"];
        jrepst_indcovb2 = indcovb2_frame["jrepst"];
        jmatst_indcovb2 = indcovb2_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovb2_zi = indcovb2_frame["sizea_zi"];
          sizeb_indcovb2_zi = indcovb2_frame["sizeb_zi"];
          sizec_indcovb2_zi = indcovb2_frame["sizec_zi"];
          fec_indcovb2_zi = indcovb2_frame["fec_zi"];
          jsizea_indcovb2_zi = indcovb2_frame["jsizea_zi"];
          jsizeb_indcovb2_zi = indcovb2_frame["jsizeb_zi"];
          jsizec_indcovb2_zi = indcovb2_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb1_frame")) {
        DataFrame indcovb1_frame = as<DataFrame>(modelsuite["indcovb1_frame"]);
        
        indcovb_names = indcovb1_frame["indcovb"];
        
        surv_indcovb1 = indcovb1_frame["surv"];
        obs_indcovb1 = indcovb1_frame["obs"];
        sizea_indcovb1 = indcovb1_frame["sizea"];
        sizeb_indcovb1 = indcovb1_frame["sizeb"];
        sizec_indcovb1 = indcovb1_frame["sizec"];
        repst_indcovb1 = indcovb1_frame["repst"];
        fec_indcovb1 = indcovb1_frame["fec"];
        
        jsurv_indcovb1 = indcovb1_frame["jsurv"];
        jobs_indcovb1 = indcovb1_frame["jobs"];
        jsizea_indcovb1 = indcovb1_frame["jsizea"];
        jsizeb_indcovb1 = indcovb1_frame["jsizeb"];
        jsizec_indcovb1 = indcovb1_frame["jsizec"];
        jrepst_indcovb1 = indcovb1_frame["jrepst"];
        jmatst_indcovb1 = indcovb1_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovb1_zi = indcovb1_frame["sizea_zi"];
          sizeb_indcovb1_zi = indcovb1_frame["sizeb_zi"];
          sizec_indcovb1_zi = indcovb1_frame["sizec_zi"];
          fec_indcovb1_zi = indcovb1_frame["fec_zi"];
          jsizea_indcovb1_zi = indcovb1_frame["jsizea_zi"];
          jsizeb_indcovb1_zi = indcovb1_frame["jsizeb_zi"];
          jsizec_indcovb1_zi = indcovb1_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc2_frame")) {
        DataFrame indcovc2_frame = as<DataFrame>(modelsuite["indcovc2_frame"]);
        
        indcovc_names = indcovc2_frame["indcovc"];
        
        surv_indcovc2 = indcovc2_frame["surv"];
        obs_indcovc2 = indcovc2_frame["obs"];
        sizea_indcovc2 = indcovc2_frame["sizea"];
        sizeb_indcovc2 = indcovc2_frame["sizeb"];
        sizec_indcovc2 = indcovc2_frame["sizec"];
        repst_indcovc2 = indcovc2_frame["repst"];
        fec_indcovc2 = indcovc2_frame["fec"];
        
        jsurv_indcovc2 = indcovc2_frame["jsurv"];
        jobs_indcovc2 = indcovc2_frame["jobs"];
        jsizea_indcovc2 = indcovc2_frame["jsizea"];
        jsizeb_indcovc2 = indcovc2_frame["jsizeb"];
        jsizec_indcovc2 = indcovc2_frame["jsizec"];
        jrepst_indcovc2 = indcovc2_frame["jrepst"];
        jmatst_indcovc2 = indcovc2_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovc2_zi = indcovc2_frame["sizea_zi"];
          sizeb_indcovc2_zi = indcovc2_frame["sizeb_zi"];
          sizec_indcovc2_zi = indcovc2_frame["sizec_zi"];
          fec_indcovc2_zi = indcovc2_frame["fec_zi"];
          jsizea_indcovc2_zi = indcovc2_frame["jsizea_zi"];
          jsizeb_indcovc2_zi = indcovc2_frame["jsizeb_zi"];
          jsizec_indcovc2_zi = indcovc2_frame["jsizec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc1_frame")) {
        DataFrame indcovc1_frame = as<DataFrame>(modelsuite["indcovc1_frame"]);
        
        indcovc_names = indcovc1_frame["indcovc"];
        
        surv_indcovc1 = indcovc1_frame["surv"];
        obs_indcovc1 = indcovc1_frame["obs"];
        sizea_indcovc1 = indcovc1_frame["sizea"];
        sizeb_indcovc1 = indcovc1_frame["sizeb"];
        sizec_indcovc1 = indcovc1_frame["sizec"];
        repst_indcovc1 = indcovc1_frame["repst"];
        fec_indcovc1 = indcovc1_frame["fec"];
        
        jsurv_indcovc1 = indcovc1_frame["jsurv"];
        jobs_indcovc1 = indcovc1_frame["jobs"];
        jsizea_indcovc1 = indcovc1_frame["jsizea"];
        jsizeb_indcovc1 = indcovc1_frame["jsizeb"];
        jsizec_indcovc1 = indcovc1_frame["jsizec"];
        jrepst_indcovc1 = indcovc1_frame["jrepst"];
        jmatst_indcovc1 = indcovc1_frame["jmatst"];
        
        if (zi_yn) {
          sizea_indcovc1_zi = indcovc1_frame["sizea_zi"];
          sizeb_indcovc1_zi = indcovc1_frame["sizeb_zi"];
          sizec_indcovc1_zi = indcovc1_frame["sizec_zi"];
          fec_indcovc1_zi = indcovc1_frame["fec_zi"];
          jsizea_indcovc1_zi = indcovc1_frame["jsizea_zi"];
          jsizeb_indcovc1_zi = indcovc1_frame["jsizeb_zi"];
          jsizec_indcovc1_zi = indcovc1_frame["jsizec_zi"];
        }
      }
    }
    
    CharacterVector list_names = {"fixed_slopes", "year_slopes", "patch_slopes",
      "group2_slopes", "dist", "zi", "fixed_zi", "year_zi", "patch_zi",
      "group2_zi", "indcova_names", "indcova2_slopes", "indcova2_zi",
      "indcovb_names", "indcovb2_slopes", "indcovb2_zi", "indcovc_names",
      "indcovc2_slopes", "indcovc2_zi", "year_names", "patch_names",
      "group_names", "main_effect_1", "main_effect_2", "sigma_theta",
      "effects_names", "group1_slopes", "group1_zi", "indcova1_slopes",
      "indcovb1_slopes", "indcovc1_slopes", "indcova1_zi", "indcovb1_zi",
      "indcovc1_zi"};
      
    List surv_list(34);
    surv_list(0) = surv_num;
    surv_list(1) = surv_year;
    surv_list(2) = surv_patch;
    surv_list(3) = surv_group2;
    surv_list(4) = surv_dist;
    surv_list(5) = false;
    surv_list(6) = dud_zi;
    surv_list(7) = dud_zi;
    surv_list(8) = dud_zi;
    surv_list(9) = dud_zi;
    surv_list(10) = indcova_names;
    surv_list(11) = surv_indcova2;
    surv_list(12) = dud_zi;
    surv_list(13) = indcovb_names;
    surv_list(14) = surv_indcovb2;
    surv_list(15) = dud_zi;
    surv_list(16) = indcovc_names;
    surv_list(17) = surv_indcovc2;
    surv_list(18) = dud_zi;
    surv_list(19) = year_names;
    surv_list(20) = patch_names;
    surv_list(21) = group_names;
    surv_list(22) = main_effect_1;
    surv_list(23) = main_effect_2;
    surv_list(24) = 1.0;
    surv_list(25) = effects_names;
    surv_list(26) = surv_group1;
    surv_list(27) = dud_zi;
    surv_list(28) = surv_indcova1;
    surv_list(29) = surv_indcovb1;
    surv_list(30) = surv_indcovc1;
    surv_list(31) = dud_zi;
    surv_list(32) = dud_zi;
    surv_list(33) = dud_zi;
    
    List obs_list(34);
    obs_list(0) = obs_num;
    obs_list(1) = obs_year;
    obs_list(2) = obs_patch;
    obs_list(3) = obs_group2;
    obs_list(4) = obs_dist;
    obs_list(5) = false;
    obs_list(6) = dud_zi;
    obs_list(7) = dud_zi;
    obs_list(8) = dud_zi;
    obs_list(9) = dud_zi;
    obs_list(10) = indcova_names;
    obs_list(11) = obs_indcova2;
    obs_list(12) = dud_zi;
    obs_list(13) = indcovb_names;
    obs_list(14) = obs_indcovb2;
    obs_list(15) = dud_zi;
    obs_list(16) = indcovc_names;
    obs_list(17) = obs_indcovc2;
    obs_list(18) = dud_zi;
    obs_list(19) = year_names;
    obs_list(20) = patch_names;
    obs_list(21) = group_names;
    obs_list(22) = main_effect_1;
    obs_list(23) = main_effect_2;
    obs_list(24) = 1.0;
    obs_list(25) = effects_names;
    obs_list(26) = obs_group1;
    obs_list(27) = dud_zi;
    obs_list(28) = obs_indcova1;
    obs_list(29) = obs_indcovb1;
    obs_list(30) = obs_indcovc1;
    obs_list(31) = dud_zi;
    obs_list(32) = dud_zi;
    obs_list(33) = dud_zi;
    
    List sizea_list(34);
    sizea_list(0) = sizea_num;
    sizea_list(1) = sizea_year;
    sizea_list(2) = sizea_patch;
    sizea_list(3) = sizea_group2;
    sizea_list(4) = sizea_dist;
    sizea_list(5) = zi_yn;
    sizea_list(6) = sizea_zi;
    sizea_list(7) = year_sizea_zi;
    sizea_list(8) = patch_sizea_zi;
    sizea_list(9) = group2_sizea_zi;
    sizea_list(10) = indcova_names;
    sizea_list(11) = sizea_indcova2;
    sizea_list(12) = sizea_indcova2_zi;
    sizea_list(13) = indcovb_names;
    sizea_list(14) = sizea_indcovb2;
    sizea_list(15) = sizea_indcovb2_zi;
    sizea_list(16) = indcovc_names;
    sizea_list(17) = sizea_indcovc2;
    sizea_list(18) = sizea_indcovc2_zi;
    sizea_list(19) = year_names;
    sizea_list(20) = patch_names;
    sizea_list(21) = group_names;
    sizea_list(22) = main_effect_1;
    sizea_list(23) = main_effect_2;
    sizea_list(24) = sizea_st;
    sizea_list(25) = effects_names;
    sizea_list(26) = sizea_group1;
    sizea_list(27) = group1_sizea_zi;
    sizea_list(28) = sizea_indcova1;
    sizea_list(29) = sizea_indcovb1;
    sizea_list(30) = sizea_indcovc1;
    sizea_list(31) = sizea_indcova1_zi;
    sizea_list(32) = sizea_indcovb1_zi;
    sizea_list(33) = sizea_indcovc1_zi;
    
    List sizeb_list(34);
    sizeb_list(0) = sizeb_num;
    sizeb_list(1) = sizeb_year;
    sizeb_list(2) = sizeb_patch;
    sizeb_list(3) = sizeb_group2;
    sizeb_list(4) = sizeb_dist;
    sizeb_list(5) = zi_yn;
    sizeb_list(6) = sizeb_zi;
    sizeb_list(7) = year_sizeb_zi;
    sizeb_list(8) = patch_sizeb_zi;
    sizeb_list(9) = group2_sizeb_zi;
    sizeb_list(10) = indcova_names;
    sizeb_list(11) = sizeb_indcova2;
    sizeb_list(12) = sizeb_indcova2_zi;
    sizeb_list(13) = indcovb_names;
    sizeb_list(14) = sizeb_indcovb2;
    sizeb_list(15) = sizeb_indcovb2_zi;
    sizeb_list(16) = indcovc_names;
    sizeb_list(17) = sizeb_indcovc2;
    sizeb_list(18) = sizeb_indcovc2_zi;
    sizeb_list(19) = year_names;
    sizeb_list(20) = patch_names;
    sizeb_list(21) = group_names;
    sizeb_list(22) = main_effect_1;
    sizeb_list(23) = main_effect_2;
    sizeb_list(24) = sizeb_st;
    sizeb_list(25) = effects_names;
    sizeb_list(26) = sizeb_group1;
    sizeb_list(27) = group1_sizeb_zi;
    sizeb_list(28) = sizeb_indcova1;
    sizeb_list(29) = sizeb_indcovb1;
    sizeb_list(30) = sizeb_indcovc1;
    sizeb_list(31) = sizeb_indcova1_zi;
    sizeb_list(32) = sizeb_indcovb1_zi;
    sizeb_list(33) = sizeb_indcovc1_zi;
    
    List sizec_list(34);
    sizec_list(0) = sizec_num;
    sizec_list(1) = sizec_year;
    sizec_list(2) = sizec_patch;
    sizec_list(3) = sizec_group2;
    sizec_list(4) = sizec_dist;
    sizec_list(5) = zi_yn;
    sizec_list(6) = sizec_zi;
    sizec_list(7) = year_sizec_zi;
    sizec_list(8) = patch_sizec_zi;
    sizec_list(9) = group2_sizec_zi;
    sizec_list(10) = indcova_names;
    sizec_list(11) = sizec_indcova2;
    sizec_list(12) = sizec_indcova2_zi;
    sizec_list(13) = indcovb_names;
    sizec_list(14) = sizec_indcovb2;
    sizec_list(15) = sizec_indcovb2_zi;
    sizec_list(16) = indcovc_names;
    sizec_list(17) = sizec_indcovc2;
    sizec_list(18) = sizec_indcovc2_zi;
    sizec_list(19) = year_names;
    sizec_list(20) = patch_names;
    sizec_list(21) = group_names;
    sizec_list(22) = main_effect_1;
    sizec_list(23) = main_effect_2;
    sizec_list(24) = sizec_st;
    sizec_list(25) = effects_names;
    sizec_list(26) = sizec_group1;
    sizec_list(27) = group1_sizec_zi;
    sizec_list(28) = sizec_indcova1;
    sizec_list(29) = sizec_indcovb1;
    sizec_list(30) = sizec_indcovc1;
    sizec_list(31) = sizec_indcova1_zi;
    sizec_list(32) = sizec_indcovb1_zi;
    sizec_list(33) = sizec_indcovc1_zi;
    
    List repst_list(34);
    repst_list(0) = repst_num;
    repst_list(1) = repst_year;
    repst_list(2) = repst_patch;
    repst_list(3) = repst_group2;
    repst_list(4) = repst_dist;
    repst_list(5) = false;
    repst_list(6) = dud_zi;
    repst_list(7) = dud_zi;
    repst_list(8) = dud_zi;
    repst_list(9) = dud_zi;
    repst_list(10) = indcova_names;
    repst_list(11) = repst_indcova2;
    repst_list(12) = dud_zi;
    repst_list(13) = indcovb_names;
    repst_list(14) = repst_indcovb2;
    repst_list(15) = dud_zi;
    repst_list(16) = indcovc_names;
    repst_list(17) = repst_indcovc2;
    repst_list(18) = dud_zi;
    repst_list(19) = year_names;
    repst_list(20) = patch_names;
    repst_list(21) = group_names;
    repst_list(22) = main_effect_1;
    repst_list(23) = main_effect_2;
    repst_list(24) = 1.0;
    repst_list(25) = effects_names;
    repst_list(26) = repst_group1;
    repst_list(27) = dud_zi;
    repst_list(28) = repst_indcova1;
    repst_list(29) = repst_indcovb1;
    repst_list(30) = repst_indcovc1;
    repst_list(31) = dud_zi;
    repst_list(32) = dud_zi;
    repst_list(33) = dud_zi;
    
    List fec_list(34);
    fec_list(0) = fec_num;
    fec_list(1) = fec_year;
    fec_list(2) = fec_patch;
    fec_list(3) = fec_group2;
    fec_list(4) = fec_dist;
    fec_list(5) = zi_yn;
    fec_list(6) = fec_zi;
    fec_list(7) = year_fec_zi;
    fec_list(8) = patch_fec_zi;
    fec_list(9) = group2_fec_zi;
    fec_list(10) = indcova_names;
    fec_list(11) = fec_indcova2;
    fec_list(12) = fec_indcova2_zi;
    fec_list(13) = indcovb_names;
    fec_list(14) = fec_indcovb2;
    fec_list(15) = fec_indcovb2_zi;
    fec_list(16) = indcovc_names;
    fec_list(17) = fec_indcovc2;
    fec_list(18) = fec_indcovc2_zi;
    fec_list(19) = year_names;
    fec_list(20) = patch_names;
    fec_list(21) = group_names;
    fec_list(22) = main_effect_1;
    fec_list(23) = main_effect_2;
    fec_list(24) = fec_st;
    fec_list(25) = effects_names;
    fec_list(26) = fec_group1;
    fec_list(27) = group1_fec_zi;
    fec_list(28) = fec_indcova1;
    fec_list(29) = fec_indcovb1;
    fec_list(30) = fec_indcovc1;
    fec_list(31) = fec_indcova1_zi;
    fec_list(32) = fec_indcovb1_zi;
    fec_list(33) = fec_indcovc1_zi;
    
    List jsurv_list(34);
    jsurv_list(0) = jsurv_num;
    jsurv_list(1) = jsurv_year;
    jsurv_list(2) = jsurv_patch;
    jsurv_list(3) = jsurv_group2;
    jsurv_list(4) = jsurv_dist;
    jsurv_list(5) = false;
    jsurv_list(6) = dud_zi;
    jsurv_list(7) = dud_zi;
    jsurv_list(8) = dud_zi;
    jsurv_list(9) = dud_zi;
    jsurv_list(10) = indcova_names;
    jsurv_list(11) = jsurv_indcova2;
    jsurv_list(12) = dud_zi;
    jsurv_list(13) = indcovb_names;
    jsurv_list(14) = jsurv_indcovb2;
    jsurv_list(15) = dud_zi;
    jsurv_list(16) = indcovc_names;
    jsurv_list(17) = jsurv_indcovc2;
    jsurv_list(18) = dud_zi;
    jsurv_list(19) = year_names;
    jsurv_list(20) = patch_names;
    jsurv_list(21) = group_names;
    jsurv_list(22) = main_effect_1;
    jsurv_list(23) = main_effect_2;
    jsurv_list(24) = 1.0;
    jsurv_list(25) = effects_names;
    jsurv_list(26) = jsurv_group1;
    jsurv_list(27) = dud_zi;
    jsurv_list(28) = jsurv_indcova1;
    jsurv_list(29) = jsurv_indcovb1;
    jsurv_list(30) = jsurv_indcovc1;
    jsurv_list(31) = dud_zi;
    jsurv_list(32) = dud_zi;
    jsurv_list(33) = dud_zi;
    
    List jobs_list(34);
    jobs_list(0) = jobs_num;
    jobs_list(1) = jobs_year;
    jobs_list(2) = jobs_patch;
    jobs_list(3) = jobs_group2;
    jobs_list(4) = jobs_dist;
    jobs_list(5) = false;
    jobs_list(6) = dud_zi;
    jobs_list(7) = dud_zi;
    jobs_list(8) = dud_zi;
    jobs_list(9) = dud_zi;
    jobs_list(10) = indcova_names;
    jobs_list(11) = jobs_indcova2;
    jobs_list(12) = dud_zi;
    jobs_list(13) = indcovb_names;
    jobs_list(14) = jobs_indcovb2;
    jobs_list(15) = dud_zi;
    jobs_list(16) = indcovc_names;
    jobs_list(17) = jobs_indcovc2;
    jobs_list(18) = dud_zi;
    jobs_list(19) = year_names;
    jobs_list(20) = patch_names;
    jobs_list(21) = group_names;
    jobs_list(22) = main_effect_1;
    jobs_list(23) = main_effect_2;
    jobs_list(24) = 1.0;
    jobs_list(25) = effects_names;
    jobs_list(26) = jobs_group1;
    jobs_list(27) = dud_zi;
    jobs_list(28) = jobs_indcova1;
    jobs_list(29) = jobs_indcovb1;
    jobs_list(30) = jobs_indcovc1;
    jobs_list(31) = dud_zi;
    jobs_list(32) = dud_zi;
    jobs_list(33) = dud_zi;
    
    List jsizea_list(34);
    jsizea_list(0) = jsizea_num;
    jsizea_list(1) = jsizea_year;
    jsizea_list(2) = jsizea_patch;
    jsizea_list(3) = jsizea_group2;
    jsizea_list(4) = jsizea_dist;
    jsizea_list(5) = zi_yn;
    jsizea_list(6) = jsizea_zi;
    jsizea_list(7) = year_jsizea_zi;
    jsizea_list(8) = patch_jsizea_zi;
    jsizea_list(9) = group2_jsizea_zi;
    jsizea_list(10) = indcova_names;
    jsizea_list(11) = jsizea_indcova2;
    jsizea_list(12) = jsizea_indcova2_zi;
    jsizea_list(13) = indcovb_names;
    jsizea_list(14) = jsizea_indcovb2;
    jsizea_list(15) = jsizea_indcovb2_zi;
    jsizea_list(16) = indcovc_names;
    jsizea_list(17) = jsizea_indcovc2;
    jsizea_list(18) = jsizea_indcovc2_zi;
    jsizea_list(19) = year_names;
    jsizea_list(20) = patch_names;
    jsizea_list(21) = group_names;
    jsizea_list(22) = main_effect_1;
    jsizea_list(23) = main_effect_2;
    jsizea_list(24) = jsizea_st;
    jsizea_list(25) = effects_names;
    jsizea_list(26) = jsizea_group1;
    jsizea_list(27) = group1_jsizea_zi;
    jsizea_list(28) = jsizea_indcova1;
    jsizea_list(29) = jsizea_indcovb1;
    jsizea_list(30) = jsizea_indcovc1;
    jsizea_list(31) = jsizea_indcova1_zi;
    jsizea_list(32) = jsizea_indcovb1_zi;
    jsizea_list(33) = jsizea_indcovc1_zi;
    
    List jsizeb_list(34);
    jsizeb_list(0) = jsizeb_num;
    jsizeb_list(1) = jsizeb_year;
    jsizeb_list(2) = jsizeb_patch;
    jsizeb_list(3) = jsizeb_group2;
    jsizeb_list(4) = jsizeb_dist;
    jsizeb_list(5) = zi_yn;
    jsizeb_list(6) = jsizeb_zi;
    jsizeb_list(7) = year_jsizeb_zi;
    jsizeb_list(8) = patch_jsizeb_zi;
    jsizeb_list(9) = group2_jsizeb_zi;
    jsizeb_list(10) = indcova_names;
    jsizeb_list(11) = jsizeb_indcova2;
    jsizeb_list(12) = jsizeb_indcova2_zi;
    jsizeb_list(13) = indcovb_names;
    jsizeb_list(14) = jsizeb_indcovb2;
    jsizeb_list(15) = jsizeb_indcovb2_zi;
    jsizeb_list(16) = indcovc_names;
    jsizeb_list(17) = jsizeb_indcovc2;
    jsizeb_list(18) = jsizeb_indcovc2_zi;
    jsizeb_list(19) = year_names;
    jsizeb_list(20) = patch_names;
    jsizeb_list(21) = group_names;
    jsizeb_list(22) = main_effect_1;
    jsizeb_list(23) = main_effect_2;
    jsizeb_list(24) = jsizeb_st;
    jsizeb_list(25) = effects_names;
    jsizeb_list(26) = jsizeb_group1;
    jsizeb_list(27) = group1_jsizeb_zi;
    jsizeb_list(28) = jsizeb_indcova1;
    jsizeb_list(29) = jsizeb_indcovb1;
    jsizeb_list(30) = jsizeb_indcovc1;
    jsizeb_list(31) = jsizeb_indcova1_zi;
    jsizeb_list(32) = jsizeb_indcovb1_zi;
    jsizeb_list(33) = jsizeb_indcovc1_zi;
    
    List jsizec_list(34);
    jsizec_list(0) = jsizec_num;
    jsizec_list(1) = jsizec_year;
    jsizec_list(2) = jsizec_patch;
    jsizec_list(3) = jsizec_group2;
    jsizec_list(4) = jsizec_dist;
    jsizec_list(5) = zi_yn;
    jsizec_list(6) = jsizec_zi;
    jsizec_list(7) = year_jsizec_zi;
    jsizec_list(8) = patch_jsizec_zi;
    jsizec_list(9) = group2_jsizec_zi;
    jsizec_list(10) = indcova_names;
    jsizec_list(11) = jsizec_indcova2;
    jsizec_list(12) = jsizec_indcova2_zi;
    jsizec_list(13) = indcovb_names;
    jsizec_list(14) = jsizec_indcovb2;
    jsizec_list(15) = jsizec_indcovb2_zi;
    jsizec_list(16) = indcovc_names;
    jsizec_list(17) = jsizec_indcovc2;
    jsizec_list(18) = jsizec_indcovc2_zi;
    jsizec_list(19) = year_names;
    jsizec_list(20) = patch_names;
    jsizec_list(21) = group_names;
    jsizec_list(22) = main_effect_1;
    jsizec_list(23) = main_effect_2;
    jsizec_list(24) = jsizec_st;
    jsizec_list(25) = effects_names;
    jsizec_list(26) = jsizec_group1;
    jsizec_list(27) = group1_jsizec_zi;
    jsizec_list(28) = jsizec_indcova1;
    jsizec_list(29) = jsizec_indcovb1;
    jsizec_list(30) = jsizec_indcovc1;
    jsizec_list(31) = jsizec_indcova1_zi;
    jsizec_list(32) = jsizec_indcovb1_zi;
    jsizec_list(33) = jsizec_indcovc1_zi;
    
    List jrepst_list(34);
    jrepst_list(0) = jrepst_num;
    jrepst_list(1) = jrepst_year;
    jrepst_list(2) = jrepst_patch;
    jrepst_list(3) = jrepst_group2;
    jrepst_list(4) = jrepst_dist;
    jrepst_list(5) = false;
    jrepst_list(6) = dud_zi;
    jrepst_list(7) = dud_zi;
    jrepst_list(8) = dud_zi;
    jrepst_list(9) = dud_zi;
    jrepst_list(10) = indcova_names;
    jrepst_list(11) = jrepst_indcova2;
    jrepst_list(12) = dud_zi;
    jrepst_list(13) = indcovb_names;
    jrepst_list(14) = jrepst_indcovb2;
    jrepst_list(15) = dud_zi;
    jrepst_list(16) = indcovc_names;
    jrepst_list(17) = jrepst_indcovc2;
    jrepst_list(18) = dud_zi;
    jrepst_list(19) = year_names;
    jrepst_list(20) = patch_names;
    jrepst_list(21) = group_names;
    jrepst_list(22) = main_effect_1;
    jrepst_list(23) = main_effect_2;
    jrepst_list(24) = 1.0;
    jrepst_list(25) = effects_names;
    jrepst_list(26) = jrepst_group1;
    jrepst_list(27) = dud_zi;
    jrepst_list(28) = jrepst_indcova1;
    jrepst_list(29) = jrepst_indcovb1;
    jrepst_list(30) = jrepst_indcovc1;
    jrepst_list(31) = dud_zi;
    jrepst_list(32) = dud_zi;
    jrepst_list(33) = dud_zi;
    
    List jmatst_list(34);
    jmatst_list(0) = jmatst_num;
    jmatst_list(1) = jmatst_year;
    jmatst_list(2) = jmatst_patch;
    jmatst_list(3) = jmatst_group2;
    jmatst_list(4) = jmatst_dist;
    jmatst_list(5) = false;
    jmatst_list(6) = dud_zi;
    jmatst_list(7) = dud_zi;
    jmatst_list(8) = dud_zi;
    jmatst_list(9) = dud_zi;
    jmatst_list(10) = indcova_names;
    jmatst_list(11) = jmatst_indcova2;
    jmatst_list(12) = dud_zi;
    jmatst_list(13) = indcovb_names;
    jmatst_list(14) = jmatst_indcovb2;
    jmatst_list(15) = dud_zi;
    jmatst_list(16) = indcovc_names;
    jmatst_list(17) = jmatst_indcovc2;
    jmatst_list(18) = dud_zi;
    jmatst_list(19) = year_names;
    jmatst_list(20) = patch_names;
    jmatst_list(21) = group_names;
    jmatst_list(22) = main_effect_1;
    jmatst_list(23) = main_effect_2;
    jmatst_list(24) = 1.0;
    jmatst_list(25) = effects_names;
    jmatst_list(26) = jmatst_group1;
    jmatst_list(27) = dud_zi;
    jmatst_list(28) = jmatst_indcova1;
    jmatst_list(29) = jmatst_indcovb1;
    jmatst_list(30) = jmatst_indcovc1;
    jmatst_list(31) = dud_zi;
    jmatst_list(32) = dud_zi;
    jmatst_list(33) = dud_zi;
    
    surv_model = surv_list;
    obs_model = obs_list;
    size_model = sizea_list;
    sizeb_model = sizeb_list;
    sizec_model = sizec_list;
    repst_model = repst_list;
    fec_model = fec_list;
    
    jsurv_model = jsurv_list;
    jobs_model = jobs_list;
    jsize_model = jsizea_list;
    jsizeb_model = jsizeb_list;
    jsizec_model = jsizec_list;
    jrepst_model = jrepst_list;
    jmatst_model = jmatst_list;
    
    surv_model.attr("names") = list_names;
    obs_model.attr("names") = list_names;
    size_model.attr("names") = list_names;
    sizeb_model.attr("names") = list_names;
    sizec_model.attr("names") = list_names;
    repst_model.attr("names") = list_names;
    fec_model.attr("names") = list_names;
    jsurv_model.attr("names") = list_names;
    jobs_model.attr("names") = list_names;
    jsize_model.attr("names") = list_names;
    jsizeb_model.attr("names") = list_names;
    jsizec_model.attr("names") = list_names;
    jrepst_model.attr("names") = list_names;
    jmatst_model.attr("names") = list_names;
    
    paramnames = as<DataFrame>(modelsuite["paramnames"]);
    CharacterVector modelparams = as<CharacterVector>(paramnames["modelparams"]);
    CharacterVector mainparams = as<CharacterVector>(paramnames["mainparams"]);
    CharacterVector parameter_names = as<CharacterVector>(paramnames["parameter_names"]);
    
    bool current_check = false;
    for (int i = 0; i < modelparams.length(); i++) {
      for (int j = 0; j < 17; j++) {
        current_check = stringcompare_hard(as<std::string>(mainparams(i)), as<std::string>(main_effect_1(j)));
        
        if (current_check) modelparams(i) = main_effect_1(j);
      }
    }
    
    paramnames = DataFrame::create(_["parameter_names"] = parameter_names,
      _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
  } else {
    // Standard lefkoMod modelsuite input
    surv_model = modelsuite["survival_model"];
    obs_model = modelsuite["observation_model"];
    size_model = modelsuite["size_model"];
    sizeb_model = modelsuite["sizeb_model"];
    sizec_model = modelsuite["sizec_model"];
    repst_model = modelsuite["repstatus_model"];
    fec_model = modelsuite["fecundity_model"];
    
    jsurv_model = modelsuite["juv_survival_model"];
    jobs_model = modelsuite["juv_observation_model"];
    jsize_model = modelsuite["juv_size_model"];
    jsizeb_model = modelsuite["juv_sizeb_model"];
    jsizec_model = modelsuite["juv_sizec_model"];
    jrepst_model = modelsuite["juv_reproduction_model"];
    jmatst_model = modelsuite["juv_maturity_model"];
    
    paramnames = as<DataFrame>(modelsuite["paramnames"]);
  }
  
  List surv_proxy = LefkoUtils::modelextract(surv_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List obs_proxy = LefkoUtils::modelextract(obs_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List size_proxy = LefkoUtils::modelextract(size_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List sizeb_proxy = LefkoUtils::modelextract(sizeb_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List sizec_proxy = LefkoUtils::modelextract(sizec_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List repst_proxy = LefkoUtils::modelextract(repst_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List fec_proxy = LefkoUtils::modelextract(fec_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  
  List jsurv_proxy = LefkoUtils::modelextract(jsurv_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List jobs_proxy = LefkoUtils::modelextract(jobs_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List jsize_proxy = LefkoUtils::modelextract(jsize_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List jsizeb_proxy = LefkoUtils::modelextract(jsizeb_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List jsizec_proxy = LefkoUtils::modelextract(jsizec_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List jrepst_proxy = LefkoUtils::modelextract(jrepst_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List jmatst_proxy = LefkoUtils::modelextract(jmatst_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  
  if (err_check) {
    List proxy_out_list (14);
    proxy_out_list(0) = surv_proxy;
    proxy_out_list(1) = obs_proxy;
    proxy_out_list(2) = size_proxy;
    proxy_out_list(3) = sizeb_proxy;
    proxy_out_list(4) = sizec_proxy;
    proxy_out_list(5) = repst_proxy;
    proxy_out_list(6) = fec_proxy;
    proxy_out_list(7) = jsurv_proxy;
    proxy_out_list(8) = jobs_proxy;
    proxy_out_list(9) = jsize_proxy;
    proxy_out_list(10) = jsizeb_proxy;
    proxy_out_list(11) = jsizec_proxy;
    proxy_out_list(12) = jrepst_proxy;
    proxy_out_list(13) = jmatst_proxy;
    
    CharacterVector pol_names = {"surv_proxy", "obs_proxy", "size_proxy",
      "sizeb_proxy", "sizec_proxy", "repst_proxy", "fec_proxy", "jsurv_proxy",
      "jobs_proxy", "jsize_proxy", "jsizeb_proxy", "jsizec_proxy",
      "jrepst_proxy", "jmatst_proxy"};
    proxy_out_list.attr("names") = pol_names;
    err_check_proxies = proxy_out_list;
  }
  
  // lefkoMat structure
  List A_mats(loy_length);
  List F_mats(loy_length);
  List U_mats(loy_length);
  List out_mats(loy_length);
  
  int yearnumber {0};
  int patchnumber {0};
  
  LogicalVector dvr_yn = {false, false, false, false, false, false, false,
    false, false, false, false, false, false, false};
  IntegerVector dvr_style = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector dvr_alpha = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  NumericVector dvr_beta = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  NumericVector dvr_dens = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  
  for (int i = 0; i < loy_length; i++) {
    Rcpp::checkUserInterrupt();
    
    yearnumber = years(i);
    patchnumber = patches(i);
    
    List madsexmadrigal_oneyear = jerzeibalowski(allstages, StageFrame,
      matrixformat, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
      repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy,
      jsizeb_proxy, jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_num,
      f1_inda_num, f2_indb_num, f1_indb_num, f2_indc_num, f1_indc_num,
      f2_annua_num, f1_annua_num, f2_annub_num, f1_annub_num, f2_annuc_num,
      f1_annuc_num, f2_inda_cat, f1_inda_cat, f2_indb_cat, f1_indb_cat,
      f2_indc_cat, f1_indc_cat, r2_inda, r1_inda, r2_indb, r1_indb,
      r2_indc, r1_indc, dev_terms, false, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
      dvr_dens, dens, fecmod, maxsize, maxsizeb, maxsizec, firstage, finalage,
      negfec, yearnumber, patchnumber, exp_tol, theta_tol, cdf, err_check,
      simplicity, sparse, false);
    
    if (!simplicity) A_mats(i) = madsexmadrigal_oneyear["A"];
    F_mats(i) = madsexmadrigal_oneyear["F"];
    U_mats(i) = madsexmadrigal_oneyear["U"];
    if (err_check ) out_mats(i) = madsexmadrigal_oneyear["out"];
  }
  
  List output;
  
  if (simplicity && err_check) {
    output = List::create(_["U"] = U_mats, _["F"] = F_mats, _["prob_out"] = out_mats,
      _["allstages"] = allstages, _["proxies"] = err_check_proxies);
  } else if (simplicity) {
    output = List::create(_["U"] = U_mats, _["F"] = F_mats);
  } else if (err_check) {
    output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats,
      _["prob_out"] = out_mats, _["allstages"] = allstages,
      _["proxies"] = err_check_proxies);
  } else {
    output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats);
  }
  
  return output;
}

//' Function Passing Models and Other Parameters to Leslie Matrix Estimator
//' 
//' This function takes the various vital rate models and other parameters and
//' coordinates them as input into function \code{fleslie()}.
//' 
//' @name mothermccooney
//' 
//' @param listofyears A data frame where the rows designate the exact order of
//' years and patches to produce matrices for.
//' @param modelsuite An object of class \code{lefkoMod}, a similarly structured
//' list object, or a \code{vrm_input} object. Survival model, fecundity model,
//' and the \code{paramnames} data frame are required if not using a
//' \code{vrm_input} object.
//' @param actualages An integer vector of all actual ages to be included in the
//' matrices, in order.
//' @param mainyears A numeric vector of all times at time \emph{t}.
//' @param mainpatches A string vector of patch names.
//' @param maingroups A string vector of stage group names.
//' @param mainindcova Typically a string vector of individual covariate
//' category names.
//' @param mainindcovb Typically a string vector of individual covariate
//' category names.
//' @param mainindcovc Typically a string vector of individual covariate
//' category names.
//' @param ageframe The modified stageframe used in matrix calculations.
//' @param supplement The supplement table used in analysis, as modified by
//' \code{age_expanded()} and other pre-MPM processing.
//' @param f2_inda_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f1_annua_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_annub_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_annub_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_annuc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_annuc_num A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of annual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_inda_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{a} at each time \emph{t} to
//' be used in analysis.
//' @param f1_inda_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{a} at each time \emph{t-1}
//' to be used in analysis.
//' @param f2_indb_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{b} at each time \emph{t} to
//' be used in analysis.
//' @param f1_indb_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{b} at each time \emph{t-1}
//' to be used in analysis.
//' @param f2_indc_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{c} at each time \emph{t} to
//' be used in analysis.
//' @param f1_indc_cat A string vector of length equal to the number of years,
//' holding categories of individual covariate \code{c} at each time \emph{t-1}
//' to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, and fecundity.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param finalage The final age to be used in analysis.
//' @param cont Denotes whether age-by-stage matrix continues past the final
//' age.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0. Defaults to \code{FALSE}.
//' @param nodata A logical value indicating whether the modelsuite contains
//' all parameter coefficients and no hfv dataset is provided (\code{TRUE}), or
//' whether an hfv dataset and a true modelsuite are provided (\code{FALSE}).
//' Defaults to \code{FALSE}.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param err_check If \code{TRUE}, then also output objects \code{prob_out}
//' and \code{allstages} for error checking purposes.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' @param sparse If \code{TRUE}, then outputs matrices in sparse format.
//' Defaults to \code{FALSE}.
//' 
//' @return A list with with up to 5 elements. In order: \code{A}: a list of A
//' matrices, or a list of \code{NULL} values if \code{simplicity = TRUE};
//' \code{U}: a list of U matrices, in the same order as \code{A}; \code{F}:
//' a list of F matrices, in the same order as \code{A}; \code{prob_out}: a list
//' of error-checking conditional probability matrices, or a list of \code{NULL}
//' values if \code{err_check = FALSE}; and \code{allstages}: a data frame
//' showing the used values of all variables used in transition calculations.
//' 
//' @keywords internal
//' @noRd
List mothermccooney(const DataFrame& listofyears, const List& modelsuite,
  const IntegerVector& actualages, const CharacterVector& mainyears,
  const CharacterVector& mainpatches, RObject maingroups, RObject mainindcova,
  RObject mainindcovb, RObject mainindcovc, const DataFrame& ageframe,
  const DataFrame& supplement, NumericVector f2_inda_num,
  NumericVector f1_inda_num, NumericVector f2_indb_num,
  NumericVector f1_indb_num, NumericVector f2_indc_num,
  NumericVector f1_indc_num, NumericVector f2_annua_num,
  NumericVector f1_annua_num, NumericVector f2_annub_num,
  NumericVector f1_annub_num, NumericVector f2_annuc_num,
  NumericVector f1_annuc_num, StringVector f2_inda_cat,
  StringVector f1_inda_cat, StringVector f2_indb_cat,
  StringVector f1_indb_cat, StringVector f2_indc_cat, StringVector f1_indc_cat,
  StringVector r2_inda, StringVector r1_inda, StringVector r2_indb,
  StringVector r1_indb, StringVector r2_indc, StringVector r1_indc,
  const NumericVector& dev_terms, double dens, double fecmod, int finalage,
  int cont, bool negfec = false, bool nodata = false, double exp_tol = 700.0,
  double theta_tol = 1e8, bool err_check = false,
  bool simplicity = false, bool sparse = false) {
  
  // Dud dens_vr inputs
  Rcpp::DataFrame dvr_frame;
  List err_check_proxies;
  
  // listofyears import and settings
  IntegerVector years = listofyears["yearorder"];
  IntegerVector patches = listofyears["patchorder"];
  int loy_length = years.length();
  
  // Deviation terms
  double surv_dev = dev_terms(0);
  double fec_dev = dev_terms(1);
  
  // Move model summaries to appropriate RObjects
  RObject surv_model;
  RObject fec_model;
  DataFrame paramnames;
  
  if (nodata) {
    DataFrame vrm_frame = as<DataFrame>(modelsuite["vrm_frame"]);
    DataFrame year_frame = as<DataFrame>(modelsuite["year_frame"]);
    DataFrame patch_frame = as<DataFrame>(modelsuite["patch_frame"]);
    DataFrame group2_frame = as<DataFrame>(modelsuite["group2_frame"]);
    DataFrame group1_frame = as<DataFrame>(modelsuite["group1_frame"]);
    DataFrame dist_frame = as<DataFrame>(modelsuite["dist_frame"]);
    NumericVector st_frame = as<NumericVector>(modelsuite["st_frame"]);
    
    CharacterVector main_effect_1 = as<CharacterVector>(vrm_frame["main_effect_1"]);
    CharacterVector effects_names = clone(main_effect_1);
    
    CharacterVector main_effect_2;
    if (main_effect_1.length() > 23) {
      main_effect_2 = as<CharacterVector>(vrm_frame["main_effect_2"]);
      
      for (int i = 0; i < main_effect_1.length(); i++) {
        if (i > 16) {
          effects_names(i) += ":";
          effects_names(i) += main_effect_2(i);
        }
      }
    }
    
    CharacterVector year_names = as<CharacterVector>(year_frame["years"]);
    CharacterVector patch_names = as<CharacterVector>(patch_frame["patches"]);
    CharacterVector group_names = as<CharacterVector>(group2_frame["groups"]);
    
    bool zi_yn = false;
    
    int vrm_length = vrm_frame.length();
    
    NumericVector surv_num = as<NumericVector>(vrm_frame["surv"]);
    NumericVector fec_num = as<NumericVector>(vrm_frame["fec"]);
    
    NumericVector surv_year = as<NumericVector>(year_frame["surv"]);
    NumericVector fec_year = as<NumericVector>(year_frame["fec"]);
    NumericVector surv_patch = as<NumericVector>(patch_frame["surv"]);
    NumericVector fec_patch = as<NumericVector>(patch_frame["fec"]);
    NumericVector surv_group2 = as<NumericVector>(group2_frame["surv"]);
    NumericVector fec_group2 = as<NumericVector>(group2_frame["fec"]);
    NumericVector surv_group1 = as<NumericVector>(group1_frame["surv"]);
    NumericVector fec_group1 = as<NumericVector>(group1_frame["fec"]);
    
    StringVector distribs = as<StringVector>(dist_frame["dist"]);
    String surv_dist = distribs(0);
    String fec_dist = distribs(6);
    
    double fec_st = st_frame(6);
    
    NumericVector fec_zi;
    NumericVector year_fec_zi;
    NumericVector patch_fec_zi;
    NumericVector group2_fec_zi;
    NumericVector group1_fec_zi;
    
    NumericVector dud_zi;
    
    if (vrm_length > 16) {
      zi_yn = true;
      
      fec_zi = as<NumericVector>(vrm_frame["fec_zi"]);
      year_fec_zi = as<NumericVector>(year_frame["fec_zi"]);
      patch_fec_zi = as<NumericVector>(patch_frame["fec_zi"]);
      group2_fec_zi = as<NumericVector>(group2_frame["fec_zi"]);
      group1_fec_zi = as<NumericVector>(group1_frame["fec_zi"]);
    }
    
    CharacterVector indcova_names;
    CharacterVector indcovb_names;
    CharacterVector indcovc_names;
    
    NumericVector surv_indcova2;
    NumericVector surv_indcovb2;
    NumericVector surv_indcovc2;
    
    NumericVector fec_indcova2;
    NumericVector fec_indcovb2;
    NumericVector fec_indcovc2;
    
    NumericVector fec_indcova2_zi;
    NumericVector fec_indcovb2_zi;
    NumericVector fec_indcovc2_zi;
    
    NumericVector surv_indcova1;
    NumericVector surv_indcovb1;
    NumericVector surv_indcovc1;
    
    NumericVector fec_indcova1;
    NumericVector fec_indcovb1;
    NumericVector fec_indcovc1;
    
    NumericVector fec_indcova1_zi;
    NumericVector fec_indcovb1_zi;
    NumericVector fec_indcovc1_zi;
    
    NumericVector surv_annucova2;
    NumericVector surv_annucovb2;
    NumericVector surv_annucovc2;
    
    NumericVector fec_annucova2;
    NumericVector fec_annucovb2;
    NumericVector fec_annucovc2;
    
    NumericVector fec_annucova2_zi;
    NumericVector fec_annucovb2_zi;
    NumericVector fec_annucovc2_zi;
    
    NumericVector surv_annucova1;
    NumericVector surv_annucovb1;
    NumericVector surv_annucovc1;
    
    NumericVector fec_annucova1;
    NumericVector fec_annucovb1;
    NumericVector fec_annucovc1;
    
    NumericVector fec_annucova1_zi;
    NumericVector fec_annucovb1_zi;
    NumericVector fec_annucovc1_zi;
    
    int modelsuite_length = modelsuite.length();
    CharacterVector modelsuite_names = modelsuite.attr("names");
    
    for (int i = 0; i < modelsuite_length; i++) {
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova2_frame")) {
        DataFrame indcova2_frame = as<DataFrame>(modelsuite["indcova2_frame"]);
        
        indcova_names = indcova2_frame["indcova"];
        surv_indcova2 = indcova2_frame["surv"];
        fec_indcova2 = indcova2_frame["fec"];
        
        if (zi_yn) {
          fec_indcova2_zi = indcova2_frame["fec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova1_frame")) {
        DataFrame indcova1_frame = as<DataFrame>(modelsuite["indcova1_frame"]);
        
        indcova_names = indcova1_frame["indcova"];
        surv_indcova1 = indcova1_frame["surv"];
        fec_indcova1 = indcova1_frame["fec"];
        
        if (zi_yn) {
          fec_indcova1_zi = indcova1_frame["fec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb2_frame")) {
        DataFrame indcovb2_frame = as<DataFrame>(modelsuite["indcovb2_frame"]);
        
        indcovb_names = indcovb2_frame["indcovb"];
        surv_indcovb2 = indcovb2_frame["surv"];
        fec_indcovb2 = indcovb2_frame["fec"];
        
        if (zi_yn) {
          fec_indcovb2_zi = indcovb2_frame["fec_zi"];
        }
      }
            
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb1_frame")) {
        DataFrame indcovb1_frame = as<DataFrame>(modelsuite["indcovb1_frame"]);
        
        indcovb_names = indcovb1_frame["indcovb"];
        surv_indcovb1 = indcovb1_frame["surv"];
        fec_indcovb1 = indcovb1_frame["fec"];
        
        if (zi_yn) {
          fec_indcovb1_zi = indcovb1_frame["fec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc2_frame")) {
        DataFrame indcovc2_frame = as<DataFrame>(modelsuite["indcovc2_frame"]);
        
        indcovc_names = indcovc2_frame["indcovc"];
        surv_indcovc2 = indcovc2_frame["surv"];
        fec_indcovc2 = indcovc2_frame["fec"];
        
        if (zi_yn) {
          fec_indcovc2_zi = indcovc2_frame["fec_zi"];
        }
      }
      
      if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc1_frame")) {
        DataFrame indcovc1_frame = as<DataFrame>(modelsuite["indcovc1_frame"]);
        
        indcovc_names = indcovc1_frame["indcovc"];
        surv_indcovc1 = indcovc1_frame["surv"];
        fec_indcovc1 = indcovc1_frame["fec"];
        
        if (zi_yn) {
          fec_indcovc1_zi = indcovc1_frame["fec_zi"];
        }
      }
    }
    
    CharacterVector list_names = {"fixed_slopes", "year_slopes", "patch_slopes",
      "group2_slopes", "dist", "zi", "fixed_zi", "year_zi", "patch_zi",
      "group2_zi", "indcova_names", "indcova2_slopes", "indcova2_zi",
      "indcovb_names", "indcovb2_slopes", "indcovb2_zi", "indcovc_names",
      "indcovc2_slopes", "indcovc2_zi", "year_names", "patch_names",
      "group_names", "main_effect_1", "main_effect_2", "sigma_theta",
      "effects_names", "group1_slopes", "group1_zi", "indcova1_slopes",
      "indcovb1_slopes", "indcovc1_slopes", "indcova1_zi", "indcovb1_zi",
      "indcovc1_zi"};
      
    List surv_list(34);
    surv_list(0) = surv_num;
    surv_list(1) = surv_year;
    surv_list(2) = surv_patch;
    surv_list(3) = surv_group2;
    surv_list(4) = surv_dist;
    surv_list(5) = false;
    surv_list(6) = dud_zi;
    surv_list(7) = dud_zi;
    surv_list(8) = dud_zi;
    surv_list(9) = dud_zi;
    surv_list(10) = indcova_names;
    surv_list(11) = surv_indcova2;
    surv_list(12) = dud_zi;
    surv_list(13) = indcovb_names;
    surv_list(14) = surv_indcovb2;
    surv_list(15) = dud_zi;
    surv_list(16) = indcovc_names;
    surv_list(17) = surv_indcovc2;
    surv_list(18) = dud_zi;
    surv_list(19) = year_names;
    surv_list(20) = patch_names;
    surv_list(21) = group_names;
    surv_list(22) = main_effect_1;
    surv_list(23) = main_effect_2;
    surv_list(24) = 1.0;
    surv_list(25) = effects_names;
    surv_list(26) = surv_group1;
    surv_list(27) = dud_zi;
    surv_list(28) = surv_indcova1;
    surv_list(29) = surv_indcovb1;
    surv_list(30) = surv_indcovc1;
    surv_list(31) = dud_zi;
    surv_list(32) = dud_zi;
    surv_list(33) = dud_zi;
    
    List fec_list(34);
    fec_list(0) = fec_num;
    fec_list(1) = fec_year;
    fec_list(2) = fec_patch;
    fec_list(3) = fec_group2;
    fec_list(4) = fec_dist;
    fec_list(5) = zi_yn;
    fec_list(6) = fec_zi;
    fec_list(7) = year_fec_zi;
    fec_list(8) = patch_fec_zi;
    fec_list(9) = group2_fec_zi;
    fec_list(10) = indcova_names;
    fec_list(11) = fec_indcova2;
    fec_list(12) = fec_indcova2_zi;
    fec_list(13) = indcovb_names;
    fec_list(14) = fec_indcovb2;
    fec_list(15) = fec_indcovb2_zi;
    fec_list(16) = indcovc_names;
    fec_list(17) = fec_indcovc2;
    fec_list(18) = fec_indcovc2_zi;
    fec_list(19) = year_names;
    fec_list(20) = patch_names;
    fec_list(21) = group_names;
    fec_list(22) = main_effect_1;
    fec_list(23) = main_effect_2;
    fec_list(24) = fec_st;
    fec_list(25) = effects_names;
    fec_list(26) = fec_group1;
    fec_list(27) = group1_fec_zi;
    fec_list(28) = fec_indcova1;
    fec_list(29) = fec_indcovb1;
    fec_list(30) = fec_indcovc1;
    fec_list(31) = fec_indcova1_zi;
    fec_list(32) = fec_indcovb1_zi;
    fec_list(33) = fec_indcovc1_zi;
    
    surv_model = surv_list;
    fec_model = fec_list;
    
    surv_model.attr("names") = list_names;
    fec_model.attr("names") = list_names;
    
    paramnames = as<DataFrame>(modelsuite["paramnames"]);
    CharacterVector modelparams = as<CharacterVector>(paramnames["modelparams"]);
    CharacterVector mainparams = as<CharacterVector>(paramnames["mainparams"]);
    CharacterVector parameter_names = as<CharacterVector>(paramnames["parameter_names"]);
    
    bool current_check = false;
    for (int i = 0; i < modelparams.length(); i++) {
      for (int j = 0; j < 17; j++) {
        current_check = stringcompare_hard(as<std::string>(mainparams(i)),
          as<std::string>(main_effect_1(j)));
        
        if (current_check) modelparams(i) = main_effect_1(j);
      }
    }
    
    paramnames = DataFrame::create(_["parameter_names"] = parameter_names,
      _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
  } else {
    surv_model = modelsuite["survival_model"];
    fec_model = modelsuite["fecundity_model"];
    paramnames = as<DataFrame>(modelsuite["paramnames"]);
  }
  
  List surv_proxy = LefkoUtils::modelextract(surv_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  List fec_proxy = LefkoUtils::modelextract(fec_model, paramnames, mainyears,
    mainpatches, maingroups, mainindcova, mainindcovb, mainindcovc, nodata);
  
  if (err_check) {
    List proxy_out_list (2);
    proxy_out_list(0) = surv_proxy;
    proxy_out_list(1) = fec_proxy;
    
    CharacterVector pol_names = {"surv_proxy", "fec_proxy"};
    proxy_out_list.attr("names") = pol_names;
    err_check_proxies = proxy_out_list;
  }
  
  // Create matrices and order them within correct list structure
  List A_mats(loy_length);
  List F_mats(loy_length);
  List U_mats(loy_length);
  
  int yearnumber {0};
  int patchnumber {0};
  
  LogicalVector dvr_yn = {false, false, false, false, false, false, false, false,
    false, false, false, false, false, false};
  IntegerVector dvr_style = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector dvr_alpha = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  NumericVector dvr_beta = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  NumericVector dvr_dens = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  
  for (int i = 0; i < loy_length; i++) {
    yearnumber = years(i);
    patchnumber = patches(i);
    
    List madsexmadrigal_oneyear = motherbalowski(actualages, ageframe,
      surv_proxy, fec_proxy, f2_inda_num, f1_inda_num, f2_indb_num, f1_indb_num,
      f2_indc_num, f1_indc_num, f2_annua_num, f1_annua_num, f2_annub_num,
      f1_annub_num, f2_annuc_num, f1_annuc_num, f2_inda_cat, f1_inda_cat,
      f2_indb_cat, f1_indb_cat, f2_indc_cat, f1_indc_cat, r2_inda, r1_inda,
      r2_indb, r1_indb, r2_indc, r1_indc, surv_dev, fec_dev, dens, fecmod,
      finalage, negfec, yearnumber, patchnumber, false, dvr_yn, dvr_style,
      dvr_alpha, dvr_beta, dvr_dens, exp_tol, theta_tol, simplicity, sparse,
      supplement);
    
    if (!simplicity) A_mats(i) = madsexmadrigal_oneyear["A"];
    F_mats(i) = madsexmadrigal_oneyear["F"];
    U_mats(i) = madsexmadrigal_oneyear["U"];
  }
  
  List output;
  
  if (simplicity) {
    if (err_check) {
      output = List::create(_["U"] = U_mats, _["F"] = F_mats,
        _["proxies"] = err_check_proxies);
    } else {
      output = List::create(_["U"] = U_mats, _["F"] = F_mats);
    }
  } else {
    if (err_check) {
      output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats,
        _["proxies"] = err_check_proxies);
    } else {
      output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats);
    }
  }
  return output;
}

//' Project Function-based Matrix Projection Model
//' 
//' Function \code{f_projection3()} develops and projects function-based matrix
//' models. Unlike \code{\link{projection3}()}, which uses matrices provided as
//' input via already created \code{lefkoMat} objects, function
//' \code{f_projection3()} creates matrices at each time step from vital rate
//' models and parameter inputs provided. Projections may be stochastic or not,
//' and may be density dependent in either case. Also handles replication.
//' 
//' @name f_projection3
//' 
//' @param data The historical vertical demographic data frame used to estimate
//' vital rates (class \code{hfvdata}), which is required to initialize times and
//' patches properly. Variable names should correspond to the naming conventions
//' in \code{\link{verticalize3}()} and \code{\link{historicalize3}()}.
//' @param format An integer indicating the kind of function-based MPM to create.
//' Possible choices include: \code{1}, Ehrlen-format historical MPM; \code{2},
//' deVries-format historical MPM; \code{3}, ahistorical MPM; \code{4},
//' age-by-stage MPM; and \code{5}, Leslie (age-based) MPM.
//' @param prebreeding A logical value indicating whether the life history model
//' is a pre-breeding model. Only used in Leslie and age-by-stage MPMs. Defaults
//' to \code{TRUE}.
//' @param start_age The age from which to start the matrix. Defaults to
//' \code{NA}, in which case age \code{1} is used if \code{prebreeding = TRUE},
//' and age \code{0} is used if \code{prebreeding = FALSE}.
//' @param last_age The final age to use in the matrix. Defaults to \code{NA}, in
//' which case the highest age in the dataset is used.
//' @param fecage_min The minimum age at which reproduction is possible. Defaults
//' to \code{NA}, which is interpreted to mean that fecundity should be assessed
//' starting in the minimum age observed in the dataset.
//' @param fecage_max The maximum age at which reproduction is possible. Defaults
//' to \code{NA}, which is interpreted to mean that fecundity should be assessed
//' until the final observed age.
//' @param cont A logical value designating whether to allow continued survival
//' of individuals past the final age noted in the stageframe, using the 
//' demographic characteristics of the final age. Defaults to \code{TRUE}.
//' @param stochastic A logical value denoting whether to conduct a stochastic
//' projection or a deterministic / cyclical projection.
//' @param standardize A logical value denoting whether to re-standardize the
//' population size to \code{1.0} at each occasion. Used in density-independent
//' simulations in which it is more important to know the general trend in
//' population growth than the explicit growth rate. Defaults to \code{FALSE}.
//' @param growthonly A logical value indicating whether to produce only the
//' projected population size at each occasion (\code{TRUE}), or also to produce
//' vectors showing the stage distribution at each occasion (\code{FALSE}).
//' Defaults to \code{TRUE}.
//' @param repvalue A logical value indicating whether to calculate reproductive
//' value vectors at each time step. Can only be set to \code{TRUE} if 
//' \code{growthonly = FALSE}. Setting to \code{TRUE} may dramatically increase
//' the duration of calculations. Defaults to \code{FALSE}.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion to the nearest
//' integer. Defaults to \code{TRUE}.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent and density
//' independent simulations. Defaults to \code{0}, which does not enforce
//' substochasticity. Alternatively, \code{1} forces all survival-transition
//' elements to range from 0.0 to 1.0, and forces fecundity to be non-negative;
//' and \code{2} forces all column rows in the survival-transition matrices to
//' total no more than 1.0, in addition to the actions outlined for option
//' \code{1}. Both settings \code{1} and \code{2} change negative fecundity
//' elements to \code{0.0}.
//' @param ipm_cdf A logical value indicating whether to estimate size
//' transitions using the cumulative density function in cases with continuous
//' distributions. Defaults to \code{TRUE}, with the midpoint method used if
//' \code{FALSE}.
//' @param nreps The number of replicate projections. Defaults to \code{1}.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' \code{10000}.
//' @param repmod A scalar multiplier of fecundity. Defaults to \code{1}.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param theta_tol A numeric value used to indicate a maximum value to theta as
//' used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' @param random_inda A logical value denoting whether to treat individual
//' covariate \code{a} as a random, categorical variable. Otherwise is treated as
//' a fixed, numeric variable. Defaults to \code{FALSE}.
//' @param random_indb A logical value denoting whether to treat individual
//' covariate \code{b} as a random, categorical variable. Otherwise is treated as
//' a fixed, numeric variable. Defaults to \code{FALSE}.
//' @param random_indc A logical value denoting whether to treat individual
//' covariate \code{c} as a random, categorical variable. Otherwise is treated as
//' a fixed, numeric variable. Defaults to \code{FALSE}.
//' @param err_check A logical value indicating whether to append extra output
//' for debugging purposes. Defaults to \code{FALSE}.
//' @param quiet A logical value indicating whether warning messages should be
//' suppressed. Defaults to \code{FALSE}.
//' @param stageframe An object of class \code{stageframe}. These objects are
//' generated by function \code{\link{sf_create}()}, and include information on
//' the size, observation status, propagule status, reproduction status,
//' immaturity status, maturity status, stage group, size bin widths, and other
//' key characteristics of each ahistorical stage. Required for all MPM formats
//' except Leslie MPMs.
//' @param supplement An optional data frame of class \code{lefkoSD} that
//' provides supplemental data that should be incorporated into the MPM. Three
//' kinds of data may be integrated this way: transitions to be estimated via the
//' use of proxy transitions, transition overwrites from the literature or
//' supplemental studies, and transition multipliers for survival and fecundity.
//' This data frame should be produced using the \code{\link{supplemental}()}
//' function. Can be used in place of or in addition to an overwrite table (see 
//' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
//' below).
//' @param repmatrix An optional reproduction matrix. This matrix is composed
//' mostly of \code{0}s, with non-zero entries acting as element identifiers and
//' multipliers for fecundity (with \code{1} equaling full fecundity). If left
//' blank, and no \code{supplement} is provided, then \code{flefko3()} will
//' assume that all stages marked as reproductive produce offspring at 1x that of
//' estimated fecundity, and that offspring production will yield the first stage
//' noted as propagule or immature. May be the dimensions of either a historical
//' or an ahistorical matrix. If the latter, then all stages will be used in
//' occasion \emph{t}-1 for each suggested ahistorical transition.
//' @param overwrite An optional data frame developed with the
//' \code{\link{overwrite}()} function describing transitions to be overwritten
//' either with given values or with other estimated transitions. Note that this
//' function supplements overwrite data provided in \code{supplement}.
//' @param modelsuite A \code{lefkoMod} object, at minimum with all required
//' best-fit vital rate models and a \code{paramnames} data frame, and following
//' the naming conventions used in this package. If given, then
//' \code{surv_model}, \code{obs_model}, \code{size_model}, \code{sizeb_model},
//' \code{sizec_model}, \code{repst_model}, \code{fec_model}, \code{jsurv_model},
//' \code{jobs_model}, \code{jsize_model}, \code{jsizeb_model},
//' \code{jsizec_model}, \code{jrepst_model}, \code{jmatst_model},
//' \code{paramnames}, \code{yearcol}, and \code{patchcol} are not required.
//' Although this is optional input, it is recommended, and without it separate
//' vital rate model inputs (named \code{XX_model}) are required.
//' @param paramnames A data frame with three columns, the first describing all
//' terms used in linear modeling, the second (must be called \code{mainparams})
//' giving the general model terms that will be used in matrix creation, and the
//' third showing the equivalent terms used in modeling (must be named
//' \code{modelparams}). Function \code{\link{create_pm}()} can be used to
//' create a skeleton \code{paramnames} object, which can then be edited. Only
//' required if \code{modelsuite} is not supplied.
//' @param year Either a single integer value corresponding to the year to
//' project, or a vector of \code{times} elements with the year to use at each
//' time step. Defaults to \code{NA}, in which the first year in the set of years
//' in the dataset is projected. If a vector shorter than \code{times} is
//' supplied, then this vector will be cycled.
//' @param patch A value of \code{NA}, a single string value corresponding to the
//' patch to project, or a vector of \code{times} elements with the patch to use
//' at each time step. If a vector shorter than \code{times} is supplied, then
//' this vector will be cycled. Note that this function currently does not
//' handle multiple projections for different patches in the same run.
//' @param sp_density Either a single numeric value of spatial density to use in
//' vital rate models in all time steps, or a vector of \code{times} elements of
//' such numeric values. Defaults to \code{NA}.
//' @param ind_terms An optional data frame with 3 columns and \code{times} rows
//' giving the values of individual covariates a, b, and c, respectively, for
//' each projected time. Unused terms must be set to \code{0} (use of \code{NA}
//' will produce errors.)
//' @param ann_terms An optional data frame with 3 columns and \code{times} rows
//' giving the values of annual covariates a, b, and c, respectively, for each
//' projected time. Unused terms must be set to \code{0} (use of \code{NA} will
//' produce errors.)
//' @param dev_terms An optional data frame with 14 columns and \code{times}
//' rows showing the values of the deviation terms to be added to each linear
//' vital rate. The column order should be: 1: survival, 2: observation, 3:
//' primary size, 4: secondary size, 5: tertiary size, 6: reproduction, 7:
//' fecundity, 8: juvenile survival, 9: juvenile observation, 10: juvenile
//' primary size, 11: juvenile secondary size, 12: juvenile tertiary size, 13:
//' juvenile reproduction, and 14: juvenile maturity transition.  Unused terms
//' must be set to \code{0} (use of \code{NA} will produce errors.)
//' @param surv_model A linear model predicting survival probability. This can 
//' be a model of class \code{glm} or \code{glmer}, and requires a predicted
//' binomial variable under a logit link. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param obs_model A linear model predicting sprouting or observation
//' probability. This can be a model of class \code{glm} or \code{glmer}, and
//' requires a predicted binomial variable under a logit link. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param size_model A linear model predicting primary size. This can be a model
//' of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
//' \code{vglm}, \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param sizeb_model A linear model predicting secondary size. This can be a
//' model of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
//' \code{vglm}, \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param sizec_model A linear model predicting tertiary size. This can be a
//' model of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
//' \code{vglm}, \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param repst_model A linear model predicting reproduction probability. This 
//' can be a model of class \code{glm} or \code{glmer}, and requires a predicted
//' binomial variable under a logit link. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param fec_model A linear model predicting fecundity. This can be a model of
//' class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl}, \code{vglm},
//' \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is provided. This
//' model must have been developed in a modeling exercise testing the impacts of
//' occasions \emph{t} and \emph{t}-1.
//' @param jsurv_model A linear model predicting juvenile survival probability.
//' This can be a model of class \code{glm} or \code{glmer}, and requires a
//' predicted binomial variable under a logit link. Ignored if \code{modelsuite}
//' is provided. This model must have been developed in a modeling exercise
//' testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jobs_model A linear model predicting juvenile sprouting or observation
//' probability. This can be a model of class \code{glm} or \code{glmer}, and
//' requires a predicted binomial variable under a logit link. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jsize_model A linear model predicting juvenile primary size. This
//' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
//' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jsizeb_model A linear model predicting juvenile secondary size. This
//' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
//' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jsizec_model A linear model predicting juvenile tertiary size. This
//' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
//' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jrepst_model A linear model predicting reproduction probability of a 
//' mature individual that was immature in time \emph{t}. This can be a model
//' of class \code{glm} or \code{glmer}, and requires a predicted binomial
//' variable under a logit link. Ignored if \code{modelsuite} is provided. This
//' model must have been developed in a modeling exercise testing the impacts of
//' occasions \emph{t} and \emph{t}-1.
//' @param jmatst_model A linear model predicting maturity probability of an 
//' individual that was immature in time \emph{t}. This can be a model of class
//' \code{glm} or \code{glmer}, and requires a predicted binomial variable under
//' a logit link. Ignored if \code{modelsuite} is provided. This model must have
//' been developed in a modeling exercise testing the impacts of occasions
//' \emph{t} and \emph{t}-1.
//' @param start_vec An optional numeric vector denoting the starting stage
//' distribution for the projection. Defaults to a single individual of each
//' stage.
//' @param start_frame An optional data frame characterizing stages, age-stages,
//' or stage-pairs that should be set to non-zero values in the starting vector,
//' and what those values should be. Can only be used with \code{lefkoMat}
//' objects.
//' @param tweights An optional numeric vector or matrix denoting the
//' probabilities of choosing each matrix in a stochastic projection. If a
//' matrix is input, then a first-order Markovian environment is assumed, in
//' which the probability of choosing a specific annual matrix depends on which
//' annual matrix is currently chosen. If a vector is input, then the choice of
//' annual matrix is assumed to be independent of the current matrix. Defaults
//' to equal weighting among matrices.
//' @param density An optional data frame describing the matrix elements that
//' will be subject to density dependence, and the exact kind of density
//' dependence that they will be subject to. The data frame used should be an
//' object of class \code{lefkoDens}, which is the output from function
//' \code{\link{density_input}()}.
//' @param density_vr An optional data frame describing density dependence
//' relationships in vital rates, if such relationships are to be assumed. The
//' data frame must be of class \code{lefkoDensVR}, which is the output from the
//' function \code{\link{density_vr}()}.
//' @param sparse A text string indicating whether to use sparse matrix encoding
//' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
//' \code{"auto"}, in which case sparse matrix encoding is used with square
//' matrices with at least 50 rows and no more than 50\% of elements with values
//' greater than zero. Can also be entered as a logical value if forced sparse
//' (\code{TRUE}) or forced dense (\code{FALSE}) projection is desired.
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
//' occasion in each replicate in each pop-patch or population. The list
//' structure is the same as in \code{\link{projection3}()}.}
//' \item{rep_value}{A list of lists of the actual reproductive value in each
//' occasion in each replicate in each pop-patch or population. The list
//' structure is the same as in \code{\link{projection3}()}.}
//' \item{pop_size}{A list of matrices showing the total population size in each
//' occasion per replicate (row within matrix) per pop-patch or population
//' (list element). Only a single pop-patch or population is allowed in
//' \code{f_projection3()}.}
//' \item{labels}{A data frame showing the order of populations and patches in
//' item \code{projection}.}
//' \item{ahstages}{The original stageframe used in the study.}
//' \item{hstages}{A data frame showing the order of historical stage pairs.}
//' \item{agestages}{A data frame showing the order of age-stage pairs.}
//' \item{labels}{A short data frame indicating the population (always \code{1}),
//' and patch (either the numeric index of the single chosen patch, or \code{1}
//' in all other cases).}
//' \item{control}{A short vector indicating the number of replicates and the
//' number of occasions projected per replicate.}
//' \item{density}{The data frame input under the density option. Only provided
//' if input by the user.}
//' \item{density_vr}{The data frame input under the density_vr option. Only
//' provided if input by the user.}
//' 
//' @section Notes:
//' Population projection can be a very time-consuming activity, and it is most
//' time-consuming when matrices need to be created at each time step. We have
//' created this function to work as quickly as possible, but some options will
//' slow analysis. First, the \code{err_check} option should always be set to
//' \code{FALSE}, as the added created output will not only slow the analysis
//' down but also potentially crash the memory if matrices are large enough.
//' Second, the \code{repvalue} option should be set to \code{FALSE} unless
//' reproductive values are genuinely needed, since this step requires
//' concurrent backward projection and so in some cases may double total run
//' time. Next, if the only needed data is the total population size and
//' age/stage structure at each time step, then setting \code{growthonly = TRUE}
//' will yield the quickest possible run time. Finally, the default behavior of
//' the function is to round down fractional values of individuals, and to stop
//' running projections (replicates) when the population drops to 0. Setting
//' \code{integeronly = FALSE} will have the impact of increasing runtime,
//' potentially dramatically, since the population can reach a point in which
//' the population size is extremely small but not equal to 0.
//' 
//' Projections with large matrices may take a long time to run. To assess the
//' likely running time, try using a low number of iterations on a single
//' replicate first. For example, set \code{nreps = 1} and \code{times = 10} for
//' a trial run. If a full run is set and takes too long, press the STOP button
//' in RStudio to cancel the projection run, or click \code{esc}.
//' 
//' This function currently allows three forms of density dependence. The first
//' modifies matrix elements on the basis of the input provided in option
//' \code{density}, and so alters matrix elements once the matrix has already
//' been created. The second form alters the vital rates estimated, and so
//' estimates matrix elements using vital rate values already modified by
//' density. This second form uses the input provided in option
//' \code{density_vr}. These two forms of density dependence utilize the
//' projected population size at some time to make these alterations. The third
//' form of density dependence also alters the vital rates, but using spatial
//' density supplied via option \code{sp_density} and only in vital rates in
//' which spatial density is included as a fixed factor in the associated
//' vital rate model.
//' 
//' When running density dependent simulations involving user-set exponents,
//' such as the beta term in the Ricker function and both the alpha and beta
//' terms in the Usher function, values above or below the computer limits may
//' cause unpredictable behavior. Noted odd behavior includes sudden shifts in
//' population size to negative values. This function produces warnings when
//' such values are used, and the values used for warnings may be reset with the
//' \code{exp_tol} term. In addition, this function resets beta values for the
//' Ricker function automatically to positive or negative \code{exp_tol}, giving
//' a warning when doing so.
//' 
//' Consistently positive population growth can quickly lead to population size
//' numbers larger than can be handled computationally. In that circumstance, a
//' continuously rising population size will suddenly become \code{NaN} for the
//' remainder of the projection.
//' 
//' This function does not reduce the dimensionality of matrices developed for
//' projection.
//' 
//' Speed can sometimes be increased by shifting from automatic sparse matrix
//' determination to forced dense or sparse matrix projection. This will most
//' likely occur when matrices have between 30 and 300 rows and columns.
//' Defaults work best when matrices are very small and dense, or very large and
//' sparse.
//' 
//' Some issues may arise in first-order Markovian stochastic projections if
//' the \code{year} argument is used. Use the matrix input in the
//' \code{tweights} argument to eliminate any years from consideration that are
//' not needed.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{density_input}()}
//' @seealso \code{\link{density_vr}()}
//' @seealso \code{\link{projection3}()}
//' @seealso \code{\link{flefko3}()}
//' @seealso \code{\link{flefko2}()}
//' @seealso \code{\link{aflefko2}()}
//' @seealso \code{\link{fleslie}()}
//' @seealso \code{\link{append_lP}()}
//' @seealso \code{\link{summary.lefkoProj}()}
//' @seealso \code{\link{plot.lefkoProj}()}
//' 
//' @examples
//' \donttest{
//' data(lathyrus)
//' 
//' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
//'   9)
//' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
//'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
//'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
//' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
//'   0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
//'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
//' 
//' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
//'   propstatus = propvector)
//' 
//' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
//'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988", 
//'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
//' 
//' lathvertln$feca2 <- round(lathvertln$feca2)
//' lathvertln$feca1 <- round(lathvertln$feca1)
//' lathvertln$feca3 <- round(lathvertln$feca3)
//' 
//' lathvertln_adults <- subset(lathvertln, stage2index > 2)
//' surv_model <- glm(alive3 ~ sizea2 + sizea1 + as.factor(patchid) +
//'   as.factor(year2), data = lathvertln_adults, family = "binomial")
//' 
//' obs_data <- subset(lathvertln_adults, alive3 == 1)
//' obs_model <- glm(obsstatus3 ~ as.factor(patchid), data = obs_data,
//'   family = "binomial")
//' 
//' size_data <- subset(obs_data, obsstatus3 == 1)
//' siz_model <- lm(sizea3 ~ sizea2 + sizea1 + repstatus1 + as.factor(patchid) +
//'   as.factor(year2), data = size_data)
//' 
//' reps_model <- glm(repstatus3 ~ sizea2 + sizea1 + as.factor(patchid) +
//'   as.factor(year2), data = size_data, family = "binomial")
//' 
//' fec_data <- subset(lathvertln_adults, repstatus2 == 1)
//' fec_model <- glm(feca2 ~ sizea2 + sizea1 + repstatus1 + as.factor(patchid),
//'   data = fec_data, family = "poisson")
//' 
//' lathvertln_juvs <- subset(lathvertln, stage2index < 3)
//' jsurv_model <- glm(alive3 ~ as.factor(patchid), data = lathvertln_juvs,
//'   family = "binomial")
//' 
//' jobs_data <- subset(lathvertln_juvs, alive3 == 1)
//' jobs_model <- glm(obsstatus3 ~ 1, family = "binomial", data = jobs_data)
//' 
//' jsize_data <- subset(jobs_data, obsstatus3 == 1)
//' jsiz_model <- lm(sizea3 ~ as.factor(year2), data = jsize_data)
//' 
//' jrepst_model <- 0
//' jmatst_model <- 1
//' 
//' mod_params <- create_pm(name_terms = TRUE)
//' mod_params$modelparams[3] <- "patchid"
//' mod_params$modelparams[4] <- "alive3"
//' mod_params$modelparams[5] <- "obsstatus3"
//' mod_params$modelparams[6] <- "sizea3"
//' mod_params$modelparams[9] <- "repstatus3"
//' mod_params$modelparams[11] <- "feca2"
//' mod_params$modelparams[12] <- "sizea2"
//' mod_params$modelparams[13] <- "sizea1"
//' mod_params$modelparams[18] <- "repstatus2"
//' mod_params$modelparams[19] <- "repstatus1"
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "mat", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "Sdl", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "Sd", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, "mat", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, "Sdl", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, "Sdl", NA, NA),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframeln, historical = TRUE)
//' 
//' # While we do not use MPMs to initialize f_projections3(), we do use MPMs to
//' # initialize functions start_input() and density_input().
//' lathmat3ln <- flefko3(year = "all", patch = "all", data = lathvertln,
//'   stageframe = lathframeln, supplement = lathsupp3, paramnames = mod_params,
//'   surv_model = surv_model, obs_model = obs_model, size_model = siz_model,
//'   repst_model = reps_model, fec_model = fec_model, jsurv_model = jsurv_model,
//'   jobs_model = jobs_model, jsize_model = jsiz_model,
//'   jrepst_model = jrepst_model, jmatst_model = jmatst_model, reduce = FALSE)
//' 
//' e3m_sv <- start_input(lathmat3ln, stage2 = "Sd", stage1 = "Sd", value = 1000)
//' 
//' dyn7 <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
//'   FALSE, FALSE, FALSE, FALSE, FALSE)
//' dst7 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' dal7 <- c(0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' dbe7 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' 
//' e3d_vr <- density_vr(density_yn = dyn7, style = dst7, alpha = dal7,
//'   beta = dbe7)
//' 
//' trial7_dvr_1 <- f_projection3(format = 1, data = lathvertln, supplement = lathsupp3,
//'   paramnames = mod_params, stageframe = lathframeln, nreps = 2,
//'   surv_model = surv_model, obs_model = obs_model, size_model = siz_model,
//'   repst_model = reps_model, fec_model = fec_model, jsurv_model = jsurv_model,
//'   jobs_model = jobs_model, jsize_model = jsiz_model,
//'   jrepst_model = jrepst_model, jmatst_model = jmatst_model,
//'   times = 100, stochastic = TRUE, standardize = FALSE, growthonly = TRUE,
//'   substoch = 0, sp_density = 0, start_frame = e3m_sv, density_vr = e3d_vr)
//' }
//' 
//' @export f_projection3
// [[Rcpp::export(f_projection3)]]
Rcpp::List f_projection3(int format, bool prebreeding = true, int start_age = NA_INTEGER,
  int last_age = NA_INTEGER, int fecage_min = NA_INTEGER, int fecage_max = NA_INTEGER,
  bool cont = true, bool stochastic = false, bool standardize = false,
  bool growthonly = true, bool repvalue = false, bool integeronly = true,
  int substoch = 0, bool ipm_cdf = true, int nreps = 1, int times = 10000,
  double repmod = 1.0, double exp_tol = 700.0, double theta_tol = 1e8,
  bool random_inda = false, bool random_indb = false, bool random_indc = false,
  bool err_check = false, bool quiet = false, Nullable<DataFrame> data = R_NilValue,
  Nullable<DataFrame> stageframe = R_NilValue, Nullable<DataFrame> supplement = R_NilValue,
  Nullable<NumericMatrix> repmatrix = R_NilValue, Nullable<DataFrame> overwrite = R_NilValue,
  Nullable<List> modelsuite = R_NilValue, Nullable<DataFrame> paramnames = R_NilValue,
  Nullable<NumericVector> year = R_NilValue, Nullable<CharacterVector> patch = R_NilValue,
  Nullable<NumericVector> sp_density = R_NilValue, Nullable<RObject> ind_terms = R_NilValue,
  Nullable<RObject> ann_terms = R_NilValue, Nullable<RObject> dev_terms = R_NilValue,
  Nullable<RObject> surv_model = R_NilValue, Nullable<RObject> obs_model = R_NilValue,
  Nullable<RObject> size_model = R_NilValue, Nullable<RObject> sizeb_model = R_NilValue,
  Nullable<RObject> sizec_model = R_NilValue, Nullable<RObject> repst_model = R_NilValue,
  Nullable<RObject> fec_model = R_NilValue, Nullable<RObject> jsurv_model = R_NilValue,
  Nullable<RObject> jobs_model = R_NilValue, Nullable<RObject> jsize_model = R_NilValue,
  Nullable<RObject> jsizeb_model = R_NilValue, Nullable<RObject> jsizec_model = R_NilValue,
  Nullable<RObject> jrepst_model = R_NilValue, Nullable<RObject> jmatst_model = R_NilValue,
  Nullable<NumericVector> start_vec = R_NilValue, Nullable<RObject> start_frame = R_NilValue,
  Nullable<RObject> tweights = R_NilValue, Nullable<RObject> density = R_NilValue,
  Nullable<RObject> density_vr = R_NilValue, Nullable<RObject> sparse = R_NilValue) {
  
  bool assume_markov {false};
  
  if (format < 1 || format > 5) pop_error("format", "", "", 5);
  if (substoch < 0 || substoch > 2) pop_error("substoch", "either 0, 1, or 2.", "", 1);
  
  if (nreps < 1) pop_error("nreps", "a positive integer", "", 1);
  if (times < 1) pop_error("times", "a positive integer", "", 1);
  if (format < 4 && (!IntegerVector::is_na(start_age) || !IntegerVector::is_na(last_age))) {
    if (!quiet) {
      Rf_warningcall(R_NilValue,
        "Start and final ages cannot be used with matrix formats 1-3. Resetting these parameters.");
    }
    start_age = 0;  
    last_age = 0;
  }
  if (growthonly && repvalue) {
    if (!quiet) {
      Rf_warningcall(R_NilValue,
        "Option repvalue cannot be TRUE if growthonly is TRUE. Resetting repvalue to FALSE.");
    }
    repvalue = false;
  }
  
  int sparse_switch {0};
  bool sparse_auto {true};
  bool sparse_bool = false;
  
  if (sparse.isNotNull()) {
    LefkoInputs::yesnoauto_to_logic(as<RObject>(sparse), "sparse", sparse_bool, sparse_auto);
    if (sparse_bool) sparse_switch = 1;
  }
  
  if (sparse_auto) {
    if (format < 3 || format == 4) {
      sparse_switch = 1;
    } else sparse_switch = 0;
  }
  
  // Vital rate models
  List msuite;
  RObject surmodl;
  RObject obsmodl;
  RObject sizmodl;
  RObject sibmodl;
  RObject sicmodl;
  RObject repmodl;
  RObject fecmodl;
  RObject jsurmodl;
  RObject jobsmodl;
  RObject jsizmodl;
  RObject jsibmodl;
  RObject jsicmodl;
  RObject jrepmodl;
  RObject jmatmodl;
  DataFrame pmnames;
  bool nodata = false; // Tag for vrm_input
  
  int modelcheck {0};
  NumericVector model1 = {1.0};
  NumericVector model0 = {0.0};
  CharacterVector modelnone {"none"};
  bool pmn_provided = false;
  
  if (modelsuite.isNotNull()) {
    Rcpp::List ms_intermediate(modelsuite);
    String ms_class = ms_intermediate.attr("class");
    msuite = ms_intermediate;
    
    if (stringcompare_simple(ms_class, "lefkoMo", false)) {
      surmodl = msuite["survival_model"];
      obsmodl = msuite["observation_model"];
      sizmodl = msuite["size_model"];
      sibmodl = msuite["sizeb_model"];
      sicmodl = msuite["sizec_model"];
      repmodl = msuite["repstatus_model"];
      fecmodl = msuite["fecundity_model"];
      jsurmodl = msuite["juv_survival_model"];
      jobsmodl = msuite["juv_observation_model"];
      jsizmodl = msuite["juv_size_model"];
      jsibmodl = msuite["juv_sizeb_model"];
      jsicmodl = msuite["juv_sizec_model"];
      jrepmodl = msuite["juv_reproduction_model"];
      jmatmodl = msuite["juv_maturity_model"];
      
      pmnames = as<DataFrame>(msuite["paramnames"]);
      
      pmn_provided = true;
      modelcheck++;
    } else if (stringcompare_simple(ms_class, "vrm", false)) {
      nodata = true;
      
      DataFrame vrm_frame = as<DataFrame>(msuite["vrm_frame"]);
      DataFrame year_frame = as<DataFrame>(msuite["year_frame"]);
      DataFrame patch_frame = as<DataFrame>(msuite["patch_frame"]);
      DataFrame group2_frame = as<DataFrame>(msuite["group2_frame"]);
      DataFrame group1_frame = as<DataFrame>(msuite["group1_frame"]);
      DataFrame dist_frame = as<DataFrame>(msuite["dist_frame"]);
      NumericVector st_frame = as<NumericVector>(msuite["st_frame"]);
      
      CharacterVector main_effect_1 = as<CharacterVector>(vrm_frame["main_effect_1"]);
      CharacterVector effects_names = clone(main_effect_1);
      
      CharacterVector main_effect_2;
      if (main_effect_1.length() > 23) {
        main_effect_2 = as<CharacterVector>(vrm_frame["main_effect_2"]);
        
        for (int i = 0; i < main_effect_1.length(); i++) {
          if (i > 16) {
            effects_names(i) += ":";
            effects_names(i) += main_effect_2(i);
          }
        }
      }
      
      CharacterVector year_names = as<CharacterVector>(year_frame["years"]);
      CharacterVector patch_names = as<CharacterVector>(patch_frame["patches"]);
      CharacterVector group_names = as<CharacterVector>(group2_frame["groups"]);
      
      bool zi_yn = false;
      
      int vrm_length = vrm_frame.length();
      
      NumericVector surv_num = as<NumericVector>(vrm_frame["surv"]);
      NumericVector obs_num = as<NumericVector>(vrm_frame["obs"]);
      NumericVector sizea_num = as<NumericVector>(vrm_frame["sizea"]);
      NumericVector sizeb_num = as<NumericVector>(vrm_frame["sizeb"]);
      NumericVector sizec_num = as<NumericVector>(vrm_frame["sizec"]);
      NumericVector repst_num = as<NumericVector>(vrm_frame["repst"]);
      NumericVector fec_num = as<NumericVector>(vrm_frame["fec"]);
      NumericVector jsurv_num = as<NumericVector>(vrm_frame["jsurv"]);
      NumericVector jobs_num = as<NumericVector>(vrm_frame["jobs"]);
      NumericVector jsizea_num = as<NumericVector>(vrm_frame["jsizea"]);
      NumericVector jsizeb_num = as<NumericVector>(vrm_frame["jsizeb"]);
      NumericVector jsizec_num = as<NumericVector>(vrm_frame["jsizec"]);
      NumericVector jrepst_num = as<NumericVector>(vrm_frame["jrepst"]);
      NumericVector jmatst_num = as<NumericVector>(vrm_frame["jmatst"]);
      
      NumericVector surv_year = as<NumericVector>(year_frame["surv"]);
      NumericVector obs_year = as<NumericVector>(year_frame["obs"]);
      NumericVector sizea_year = as<NumericVector>(year_frame["sizea"]);
      NumericVector sizeb_year = as<NumericVector>(year_frame["sizeb"]);
      NumericVector sizec_year = as<NumericVector>(year_frame["sizec"]);
      NumericVector repst_year = as<NumericVector>(year_frame["repst"]);
      NumericVector fec_year = as<NumericVector>(year_frame["fec"]);
      NumericVector jsurv_year = as<NumericVector>(year_frame["jsurv"]);
      NumericVector jobs_year = as<NumericVector>(year_frame["jobs"]);
      NumericVector jsizea_year = as<NumericVector>(year_frame["jsizea"]);
      NumericVector jsizeb_year = as<NumericVector>(year_frame["jsizeb"]);
      NumericVector jsizec_year = as<NumericVector>(year_frame["jsizec"]);
      NumericVector jrepst_year = as<NumericVector>(year_frame["jrepst"]);
      NumericVector jmatst_year = as<NumericVector>(year_frame["jmatst"]);
      
      NumericVector surv_patch = as<NumericVector>(patch_frame["surv"]);
      NumericVector obs_patch = as<NumericVector>(patch_frame["obs"]);
      NumericVector sizea_patch = as<NumericVector>(patch_frame["sizea"]);
      NumericVector sizeb_patch = as<NumericVector>(patch_frame["sizeb"]);
      NumericVector sizec_patch = as<NumericVector>(patch_frame["sizec"]);
      NumericVector repst_patch = as<NumericVector>(patch_frame["repst"]);
      NumericVector fec_patch = as<NumericVector>(patch_frame["fec"]);
      NumericVector jsurv_patch = as<NumericVector>(patch_frame["jsurv"]);
      NumericVector jobs_patch = as<NumericVector>(patch_frame["jobs"]);
      NumericVector jsizea_patch = as<NumericVector>(patch_frame["jsizea"]);
      NumericVector jsizeb_patch = as<NumericVector>(patch_frame["jsizeb"]);
      NumericVector jsizec_patch = as<NumericVector>(patch_frame["jsizec"]);
      NumericVector jrepst_patch = as<NumericVector>(patch_frame["jrepst"]);
      NumericVector jmatst_patch = as<NumericVector>(patch_frame["jmatst"]);
      
      NumericVector surv_group2 = as<NumericVector>(group2_frame["surv"]);
      NumericVector obs_group2 = as<NumericVector>(group2_frame["obs"]);
      NumericVector sizea_group2 = as<NumericVector>(group2_frame["sizea"]);
      NumericVector sizeb_group2 = as<NumericVector>(group2_frame["sizeb"]);
      NumericVector sizec_group2 = as<NumericVector>(group2_frame["sizec"]);
      NumericVector repst_group2 = as<NumericVector>(group2_frame["repst"]);
      NumericVector fec_group2 = as<NumericVector>(group2_frame["fec"]);
      NumericVector jsurv_group2 = as<NumericVector>(group2_frame["jsurv"]);
      NumericVector jobs_group2 = as<NumericVector>(group2_frame["jobs"]);
      NumericVector jsizea_group2 = as<NumericVector>(group2_frame["jsizea"]);
      NumericVector jsizeb_group2 = as<NumericVector>(group2_frame["jsizeb"]);
      NumericVector jsizec_group2 = as<NumericVector>(group2_frame["jsizec"]);
      NumericVector jrepst_group2 = as<NumericVector>(group2_frame["jrepst"]);
      NumericVector jmatst_group2 = as<NumericVector>(group2_frame["jmatst"]);
      
      NumericVector surv_group1 = as<NumericVector>(group1_frame["surv"]);
      NumericVector obs_group1 = as<NumericVector>(group1_frame["obs"]);
      NumericVector sizea_group1 = as<NumericVector>(group1_frame["sizea"]);
      NumericVector sizeb_group1 = as<NumericVector>(group1_frame["sizeb"]);
      NumericVector sizec_group1 = as<NumericVector>(group1_frame["sizec"]);
      NumericVector repst_group1 = as<NumericVector>(group1_frame["repst"]);
      NumericVector fec_group1 = as<NumericVector>(group1_frame["fec"]);
      NumericVector jsurv_group1 = as<NumericVector>(group1_frame["jsurv"]);
      NumericVector jobs_group1 = as<NumericVector>(group1_frame["jobs"]);
      NumericVector jsizea_group1 = as<NumericVector>(group1_frame["jsizea"]);
      NumericVector jsizeb_group1 = as<NumericVector>(group1_frame["jsizeb"]);
      NumericVector jsizec_group1 = as<NumericVector>(group1_frame["jsizec"]);
      NumericVector jrepst_group1 = as<NumericVector>(group1_frame["jrepst"]);
      NumericVector jmatst_group1 = as<NumericVector>(group1_frame["jmatst"]);
      
      StringVector distribs = as<StringVector>(dist_frame["dist"]);
      String surv_dist = distribs(0);
      String obs_dist = distribs(1);
      String sizea_dist = distribs(2);
      String sizeb_dist = distribs(3);
      String sizec_dist = distribs(4);
      String repst_dist = distribs(5);
      String fec_dist = distribs(6);
      String jsurv_dist = distribs(7);
      String jobs_dist = distribs(8);
      String jsizea_dist = distribs(9);
      String jsizeb_dist = distribs(10);
      String jsizec_dist = distribs(11);
      String jrepst_dist = distribs(12);
      String jmatst_dist = distribs(13);
      
      double sizea_st = st_frame(2);
      double sizeb_st = st_frame(3);
      double sizec_st = st_frame(4);
      double fec_st = st_frame(6);
      double jsizea_st = st_frame(9);
      double jsizeb_st = st_frame(10);
      double jsizec_st = st_frame(11);
      
      NumericVector sizea_zi;
      NumericVector sizeb_zi;
      NumericVector sizec_zi;
      NumericVector fec_zi;
      NumericVector jsizea_zi;
      NumericVector jsizeb_zi;
      NumericVector jsizec_zi;
      
      NumericVector year_sizea_zi;
      NumericVector year_sizeb_zi;
      NumericVector year_sizec_zi;
      NumericVector year_fec_zi;
      NumericVector year_jsizea_zi;
      NumericVector year_jsizeb_zi;
      NumericVector year_jsizec_zi;
      
      NumericVector patch_sizea_zi;
      NumericVector patch_sizeb_zi;
      NumericVector patch_sizec_zi;
      NumericVector patch_fec_zi;
      NumericVector patch_jsizea_zi;
      NumericVector patch_jsizeb_zi;
      NumericVector patch_jsizec_zi;
      
      NumericVector group2_sizea_zi;
      NumericVector group2_sizeb_zi;
      NumericVector group2_sizec_zi;
      NumericVector group2_fec_zi;
      NumericVector group2_jsizea_zi;
      NumericVector group2_jsizeb_zi;
      NumericVector group2_jsizec_zi;
      
      NumericVector group1_sizea_zi;
      NumericVector group1_sizeb_zi;
      NumericVector group1_sizec_zi;
      NumericVector group1_fec_zi;
      NumericVector group1_jsizea_zi;
      NumericVector group1_jsizeb_zi;
      NumericVector group1_jsizec_zi;
      
      NumericVector dud_zi;
      
      if (vrm_length > 16) {
        zi_yn = true;
        
        sizea_zi = as<NumericVector>(vrm_frame["sizea_zi"]);
        sizeb_zi = as<NumericVector>(vrm_frame["sizeb_zi"]);
        sizec_zi = as<NumericVector>(vrm_frame["sizec_zi"]);
        fec_zi = as<NumericVector>(vrm_frame["fec_zi"]);
        jsizea_zi = as<NumericVector>(vrm_frame["jsizea_zi"]);
        jsizeb_zi = as<NumericVector>(vrm_frame["jsizeb_zi"]);
        jsizec_zi = as<NumericVector>(vrm_frame["jsizec_zi"]);
        
        year_sizea_zi = as<NumericVector>(year_frame["sizea_zi"]);
        year_sizeb_zi = as<NumericVector>(year_frame["sizeb_zi"]);
        year_sizec_zi = as<NumericVector>(year_frame["sizec_zi"]);
        year_fec_zi = as<NumericVector>(year_frame["fec_zi"]);
        year_jsizea_zi = as<NumericVector>(year_frame["jsizea_zi"]);
        year_jsizeb_zi = as<NumericVector>(year_frame["jsizeb_zi"]);
        year_jsizec_zi = as<NumericVector>(year_frame["jsizec_zi"]);
        
        patch_sizea_zi = as<NumericVector>(patch_frame["sizea_zi"]);
        patch_sizeb_zi = as<NumericVector>(patch_frame["sizeb_zi"]);
        patch_sizec_zi = as<NumericVector>(patch_frame["sizec_zi"]);
        patch_fec_zi = as<NumericVector>(patch_frame["fec_zi"]);
        patch_jsizea_zi = as<NumericVector>(patch_frame["jsizea_zi"]);
        patch_jsizeb_zi = as<NumericVector>(patch_frame["jsizeb_zi"]);
        patch_jsizec_zi = as<NumericVector>(patch_frame["jsizec_zi"]);
        
        group2_sizea_zi = as<NumericVector>(group2_frame["sizea_zi"]);
        group2_sizeb_zi = as<NumericVector>(group2_frame["sizeb_zi"]);
        group2_sizec_zi = as<NumericVector>(group2_frame["sizec_zi"]);
        group2_fec_zi = as<NumericVector>(group2_frame["fec_zi"]);
        group2_jsizea_zi = as<NumericVector>(group2_frame["jsizea_zi"]);
        group2_jsizeb_zi = as<NumericVector>(group2_frame["jsizeb_zi"]);
        group2_jsizec_zi = as<NumericVector>(group2_frame["jsizec_zi"]);
        
        group1_sizea_zi = as<NumericVector>(group1_frame["sizea_zi"]);
        group1_sizeb_zi = as<NumericVector>(group1_frame["sizeb_zi"]);
        group1_sizec_zi = as<NumericVector>(group1_frame["sizec_zi"]);
        group1_fec_zi = as<NumericVector>(group1_frame["fec_zi"]);
        group1_jsizea_zi = as<NumericVector>(group1_frame["jsizea_zi"]);
        group1_jsizeb_zi = as<NumericVector>(group1_frame["jsizeb_zi"]);
        group1_jsizec_zi = as<NumericVector>(group1_frame["jsizec_zi"]);
      }
      
      CharacterVector indcova_names;
      CharacterVector indcovb_names;
      CharacterVector indcovc_names;
      
      NumericVector surv_indcova2;
      NumericVector surv_indcovb2;
      NumericVector surv_indcovc2;
      NumericVector obs_indcova2;
      NumericVector obs_indcovb2;
      NumericVector obs_indcovc2;
      NumericVector sizea_indcova2;
      NumericVector sizea_indcovb2;
      NumericVector sizea_indcovc2;
      NumericVector sizeb_indcova2;
      NumericVector sizeb_indcovb2;
      NumericVector sizeb_indcovc2;
      NumericVector sizec_indcova2;
      NumericVector sizec_indcovb2;
      NumericVector sizec_indcovc2;
      NumericVector repst_indcova2;
      NumericVector repst_indcovb2;
      NumericVector repst_indcovc2;
      NumericVector fec_indcova2;
      NumericVector fec_indcovb2;
      NumericVector fec_indcovc2;
      NumericVector jsurv_indcova2;
      NumericVector jsurv_indcovb2;
      NumericVector jsurv_indcovc2;
      NumericVector jobs_indcova2;
      NumericVector jobs_indcovb2;
      NumericVector jobs_indcovc2;
      NumericVector jsizea_indcova2;
      NumericVector jsizea_indcovb2;
      NumericVector jsizea_indcovc2;
      NumericVector jsizeb_indcova2;
      NumericVector jsizeb_indcovb2;
      NumericVector jsizeb_indcovc2;
      NumericVector jsizec_indcova2;
      NumericVector jsizec_indcovb2;
      NumericVector jsizec_indcovc2;
      NumericVector jrepst_indcova2;
      NumericVector jrepst_indcovb2;
      NumericVector jrepst_indcovc2;
      NumericVector jmatst_indcova2;
      NumericVector jmatst_indcovb2;
      NumericVector jmatst_indcovc2;
      
      NumericVector sizea_indcova2_zi;
      NumericVector sizea_indcovb2_zi;
      NumericVector sizea_indcovc2_zi;
      NumericVector sizeb_indcova2_zi;
      NumericVector sizeb_indcovb2_zi;
      NumericVector sizeb_indcovc2_zi;
      NumericVector sizec_indcova2_zi;
      NumericVector sizec_indcovb2_zi;
      NumericVector sizec_indcovc2_zi;
      NumericVector fec_indcova2_zi;
      NumericVector fec_indcovb2_zi;
      NumericVector fec_indcovc2_zi;
      NumericVector jsizea_indcova2_zi;
      NumericVector jsizea_indcovb2_zi;
      NumericVector jsizea_indcovc2_zi;
      NumericVector jsizeb_indcova2_zi;
      NumericVector jsizeb_indcovb2_zi;
      NumericVector jsizeb_indcovc2_zi;
      NumericVector jsizec_indcova2_zi;
      NumericVector jsizec_indcovb2_zi;
      NumericVector jsizec_indcovc2_zi;
      
      NumericVector surv_indcova1;
      NumericVector surv_indcovb1;
      NumericVector surv_indcovc1;
      NumericVector obs_indcova1;
      NumericVector obs_indcovb1;
      NumericVector obs_indcovc1;
      NumericVector sizea_indcova1;
      NumericVector sizea_indcovb1;
      NumericVector sizea_indcovc1;
      NumericVector sizeb_indcova1;
      NumericVector sizeb_indcovb1;
      NumericVector sizeb_indcovc1;
      NumericVector sizec_indcova1;
      NumericVector sizec_indcovb1;
      NumericVector sizec_indcovc1;
      NumericVector repst_indcova1;
      NumericVector repst_indcovb1;
      NumericVector repst_indcovc1;
      NumericVector fec_indcova1;
      NumericVector fec_indcovb1;
      NumericVector fec_indcovc1;
      NumericVector jsurv_indcova1;
      NumericVector jsurv_indcovb1;
      NumericVector jsurv_indcovc1;
      NumericVector jobs_indcova1;
      NumericVector jobs_indcovb1;
      NumericVector jobs_indcovc1;
      NumericVector jsizea_indcova1;
      NumericVector jsizea_indcovb1;
      NumericVector jsizea_indcovc1;
      NumericVector jsizeb_indcova1;
      NumericVector jsizeb_indcovb1;
      NumericVector jsizeb_indcovc1;
      NumericVector jsizec_indcova1;
      NumericVector jsizec_indcovb1;
      NumericVector jsizec_indcovc1;
      NumericVector jrepst_indcova1;
      NumericVector jrepst_indcovb1;
      NumericVector jrepst_indcovc1;
      NumericVector jmatst_indcova1;
      NumericVector jmatst_indcovb1;
      NumericVector jmatst_indcovc1;
      
      NumericVector sizea_indcova1_zi;
      NumericVector sizea_indcovb1_zi;
      NumericVector sizea_indcovc1_zi;
      NumericVector sizeb_indcova1_zi;
      NumericVector sizeb_indcovb1_zi;
      NumericVector sizeb_indcovc1_zi;
      NumericVector sizec_indcova1_zi;
      NumericVector sizec_indcovb1_zi;
      NumericVector sizec_indcovc1_zi;
      NumericVector fec_indcova1_zi;
      NumericVector fec_indcovb1_zi;
      NumericVector fec_indcovc1_zi;
      NumericVector jsizea_indcova1_zi;
      NumericVector jsizea_indcovb1_zi;
      NumericVector jsizea_indcovc1_zi;
      NumericVector jsizeb_indcova1_zi;
      NumericVector jsizeb_indcovb1_zi;
      NumericVector jsizeb_indcovc1_zi;
      NumericVector jsizec_indcova1_zi;
      NumericVector jsizec_indcovb1_zi;
      NumericVector jsizec_indcovc1_zi;
      
      int modelsuite_length = msuite.length();
      CharacterVector modelsuite_names = msuite.attr("names");
      
      for (int i = 0; i < modelsuite_length; i++) {
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova2_frame")) {
          DataFrame indcova2_frame = as<DataFrame>(msuite["indcova2_frame"]);
          
          indcova_names = indcova2_frame["indcova"];
          
          surv_indcova2 = indcova2_frame["surv"];
          obs_indcova2 = indcova2_frame["obs"];
          sizea_indcova2 = indcova2_frame["sizea"];
          sizeb_indcova2 = indcova2_frame["sizeb"];
          sizec_indcova2 = indcova2_frame["sizec"];
          repst_indcova2 = indcova2_frame["repst"];
          fec_indcova2 = indcova2_frame["fec"];
          
          jsurv_indcova2 = indcova2_frame["jsurv"];
          jobs_indcova2 = indcova2_frame["jobs"];
          jsizea_indcova2 = indcova2_frame["jsizea"];
          jsizeb_indcova2 = indcova2_frame["jsizeb"];
          jsizec_indcova2 = indcova2_frame["jsizec"];
          jrepst_indcova2 = indcova2_frame["jrepst"];
          jmatst_indcova2 = indcova2_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcova2_zi = indcova2_frame["sizea_zi"];
            sizeb_indcova2_zi = indcova2_frame["sizeb_zi"];
            sizec_indcova2_zi = indcova2_frame["sizec_zi"];
            fec_indcova2_zi = indcova2_frame["fec_zi"];
            jsizea_indcova2_zi = indcova2_frame["jsizea_zi"];
            jsizeb_indcova2_zi = indcova2_frame["jsizeb_zi"];
            jsizec_indcova2_zi = indcova2_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcova1_frame")) {
          DataFrame indcova1_frame = as<DataFrame>(msuite["indcova1_frame"]);
          
          indcova_names = indcova1_frame["indcova"];
          
          surv_indcova1 = indcova1_frame["surv"];
          obs_indcova1 = indcova1_frame["obs"];
          sizea_indcova1 = indcova1_frame["sizea"];
          sizeb_indcova1 = indcova1_frame["sizeb"];
          sizec_indcova1 = indcova1_frame["sizec"];
          repst_indcova1 = indcova1_frame["repst"];
          fec_indcova1 = indcova1_frame["fec"];
          
          jsurv_indcova1 = indcova1_frame["jsurv"];
          jobs_indcova1 = indcova1_frame["jobs"];
          jsizea_indcova1 = indcova1_frame["jsizea"];
          jsizeb_indcova1 = indcova1_frame["jsizeb"];
          jsizec_indcova1 = indcova1_frame["jsizec"];
          jrepst_indcova1 = indcova1_frame["jrepst"];
          jmatst_indcova1 = indcova1_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcova1_zi = indcova1_frame["sizea_zi"];
            sizeb_indcova1_zi = indcova1_frame["sizeb_zi"];
            sizec_indcova1_zi = indcova1_frame["sizec_zi"];
            fec_indcova1_zi = indcova1_frame["fec_zi"];
            jsizea_indcova1_zi = indcova1_frame["jsizea_zi"];
            jsizeb_indcova1_zi = indcova1_frame["jsizeb_zi"];
            jsizec_indcova1_zi = indcova1_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb2_frame")) {
          DataFrame indcovb2_frame = as<DataFrame>(msuite["indcovb2_frame"]);
          
          indcovb_names = indcovb2_frame["indcovb"];
          
          surv_indcovb2 = indcovb2_frame["surv"];
          obs_indcovb2 = indcovb2_frame["obs"];
          sizea_indcovb2 = indcovb2_frame["sizea"];
          sizeb_indcovb2 = indcovb2_frame["sizeb"];
          sizec_indcovb2 = indcovb2_frame["sizec"];
          repst_indcovb2 = indcovb2_frame["repst"];
          fec_indcovb2 = indcovb2_frame["fec"];
          
          jsurv_indcovb2 = indcovb2_frame["jsurv"];
          jobs_indcovb2 = indcovb2_frame["jobs"];
          jsizea_indcovb2 = indcovb2_frame["jsizea"];
          jsizeb_indcovb2 = indcovb2_frame["jsizeb"];
          jsizec_indcovb2 = indcovb2_frame["jsizec"];
          jrepst_indcovb2 = indcovb2_frame["jrepst"];
          jmatst_indcovb2 = indcovb2_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovb2_zi = indcovb2_frame["sizea_zi"];
            sizeb_indcovb2_zi = indcovb2_frame["sizeb_zi"];
            sizec_indcovb2_zi = indcovb2_frame["sizec_zi"];
            fec_indcovb2_zi = indcovb2_frame["fec_zi"];
            jsizea_indcovb2_zi = indcovb2_frame["jsizea_zi"];
            jsizeb_indcovb2_zi = indcovb2_frame["jsizeb_zi"];
            jsizec_indcovb2_zi = indcovb2_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovb1_frame")) {
          DataFrame indcovb1_frame = as<DataFrame>(msuite["indcovb1_frame"]);
          
          indcovb_names = indcovb1_frame["indcovb"];
          
          surv_indcovb1 = indcovb1_frame["surv"];
          obs_indcovb1 = indcovb1_frame["obs"];
          sizea_indcovb1 = indcovb1_frame["sizea"];
          sizeb_indcovb1 = indcovb1_frame["sizeb"];
          sizec_indcovb1 = indcovb1_frame["sizec"];
          repst_indcovb1 = indcovb1_frame["repst"];
          fec_indcovb1 = indcovb1_frame["fec"];
          
          jsurv_indcovb1 = indcovb1_frame["jsurv"];
          jobs_indcovb1 = indcovb1_frame["jobs"];
          jsizea_indcovb1 = indcovb1_frame["jsizea"];
          jsizeb_indcovb1 = indcovb1_frame["jsizeb"];
          jsizec_indcovb1 = indcovb1_frame["jsizec"];
          jrepst_indcovb1 = indcovb1_frame["jrepst"];
          jmatst_indcovb1 = indcovb1_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovb1_zi = indcovb1_frame["sizea_zi"];
            sizeb_indcovb1_zi = indcovb1_frame["sizeb_zi"];
            sizec_indcovb1_zi = indcovb1_frame["sizec_zi"];
            fec_indcovb1_zi = indcovb1_frame["fec_zi"];
            jsizea_indcovb1_zi = indcovb1_frame["jsizea_zi"];
            jsizeb_indcovb1_zi = indcovb1_frame["jsizeb_zi"];
            jsizec_indcovb1_zi = indcovb1_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc2_frame")) {
          DataFrame indcovc2_frame = as<DataFrame>(msuite["indcovc2_frame"]);
          
          indcovc_names = indcovc2_frame["indcovc"];
          
          surv_indcovc2 = indcovc2_frame["surv"];
          obs_indcovc2 = indcovc2_frame["obs"];
          sizea_indcovc2 = indcovc2_frame["sizea"];
          sizeb_indcovc2 = indcovc2_frame["sizeb"];
          sizec_indcovc2 = indcovc2_frame["sizec"];
          repst_indcovc2 = indcovc2_frame["repst"];
          fec_indcovc2 = indcovc2_frame["fec"];
          
          jsurv_indcovc2 = indcovc2_frame["jsurv"];
          jobs_indcovc2 = indcovc2_frame["jobs"];
          jsizea_indcovc2 = indcovc2_frame["jsizea"];
          jsizeb_indcovc2 = indcovc2_frame["jsizeb"];
          jsizec_indcovc2 = indcovc2_frame["jsizec"];
          jrepst_indcovc2 = indcovc2_frame["jrepst"];
          jmatst_indcovc2 = indcovc2_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovc2_zi = indcovc2_frame["sizea_zi"];
            sizeb_indcovc2_zi = indcovc2_frame["sizeb_zi"];
            sizec_indcovc2_zi = indcovc2_frame["sizec_zi"];
            fec_indcovc2_zi = indcovc2_frame["fec_zi"];
            jsizea_indcovc2_zi = indcovc2_frame["jsizea_zi"];
            jsizeb_indcovc2_zi = indcovc2_frame["jsizeb_zi"];
            jsizec_indcovc2_zi = indcovc2_frame["jsizec_zi"];
          }
        }
        
        if (stringcompare_hard(as<std::string>(modelsuite_names[i]), "indcovc1_frame")) {
          DataFrame indcovc1_frame = as<DataFrame>(msuite["indcovc1_frame"]);
          
          indcovc_names = indcovc1_frame["indcovc"];
          
          surv_indcovc1 = indcovc1_frame["surv"];
          obs_indcovc1 = indcovc1_frame["obs"];
          sizea_indcovc1 = indcovc1_frame["sizea"];
          sizeb_indcovc1 = indcovc1_frame["sizeb"];
          sizec_indcovc1 = indcovc1_frame["sizec"];
          repst_indcovc1 = indcovc1_frame["repst"];
          fec_indcovc1 = indcovc1_frame["fec"];
          
          jsurv_indcovc1 = indcovc1_frame["jsurv"];
          jobs_indcovc1 = indcovc1_frame["jobs"];
          jsizea_indcovc1 = indcovc1_frame["jsizea"];
          jsizeb_indcovc1 = indcovc1_frame["jsizeb"];
          jsizec_indcovc1 = indcovc1_frame["jsizec"];
          jrepst_indcovc1 = indcovc1_frame["jrepst"];
          jmatst_indcovc1 = indcovc1_frame["jmatst"];
          
          if (zi_yn) {
            sizea_indcovc1_zi = indcovc1_frame["sizea_zi"];
            sizeb_indcovc1_zi = indcovc1_frame["sizeb_zi"];
            sizec_indcovc1_zi = indcovc1_frame["sizec_zi"];
            fec_indcovc1_zi = indcovc1_frame["fec_zi"];
            jsizea_indcovc1_zi = indcovc1_frame["jsizea_zi"];
            jsizeb_indcovc1_zi = indcovc1_frame["jsizeb_zi"];
            jsizec_indcovc1_zi = indcovc1_frame["jsizec_zi"];
          }
        }
      }
      
      CharacterVector list_names = {"fixed_slopes", "year_slopes", "patch_slopes",
        "group2_slopes", "dist", "zi", "fixed_zi", "year_zi", "patch_zi",
        "group2_zi", "indcova_names", "indcova2_slopes", "indcova2_zi",
        "indcovb_names", "indcovb2_slopes", "indcovb2_zi", "indcovc_names",
        "indcovc2_slopes", "indcovc2_zi", "year_names", "patch_names",
        "group_names", "main_effect_1", "main_effect_2", "sigma_theta",
        "effects_names", "group1_slopes", "group1_zi", "indcova1_slopes",
        "indcovb1_slopes", "indcovc1_slopes", "indcova1_zi", "indcovb1_zi",
        "indcovc1_zi"};
        
      List surv_list(34);
      surv_list(0) = surv_num;
      surv_list(1) = surv_year;
      surv_list(2) = surv_patch;
      surv_list(3) = surv_group2;
      surv_list(4) = surv_dist;
      surv_list(5) = false;
      surv_list(6) = dud_zi;
      surv_list(7) = dud_zi;
      surv_list(8) = dud_zi;
      surv_list(9) = dud_zi;
      surv_list(10) = indcova_names;
      surv_list(11) = surv_indcova2;
      surv_list(12) = dud_zi;
      surv_list(13) = indcovb_names;
      surv_list(14) = surv_indcovb2;
      surv_list(15) = dud_zi;
      surv_list(16) = indcovc_names;
      surv_list(17) = surv_indcovc2;
      surv_list(18) = dud_zi;
      surv_list(19) = year_names;
      surv_list(20) = patch_names;
      surv_list(21) = group_names;
      surv_list(22) = main_effect_1;
      surv_list(23) = main_effect_2;
      surv_list(24) = 1.0;
      surv_list(25) = effects_names;
      surv_list(26) = surv_group1;
      surv_list(27) = dud_zi;
      surv_list(28) = surv_indcova1;
      surv_list(29) = surv_indcovb1;
      surv_list(30) = surv_indcovc1;
      surv_list(31) = dud_zi;
      surv_list(32) = dud_zi;
      surv_list(33) = dud_zi;
      
      List obs_list(34);
      obs_list(0) = obs_num;
      obs_list(1) = obs_year;
      obs_list(2) = obs_patch;
      obs_list(3) = obs_group2;
      obs_list(4) = obs_dist;
      obs_list(5) = false;
      obs_list(6) = dud_zi;
      obs_list(7) = dud_zi;
      obs_list(8) = dud_zi;
      obs_list(9) = dud_zi;
      obs_list(10) = indcova_names;
      obs_list(11) = obs_indcova2;
      obs_list(12) = dud_zi;
      obs_list(13) = indcovb_names;
      obs_list(14) = obs_indcovb2;
      obs_list(15) = dud_zi;
      obs_list(16) = indcovc_names;
      obs_list(17) = obs_indcovc2;
      obs_list(18) = dud_zi;
      obs_list(19) = year_names;
      obs_list(20) = patch_names;
      obs_list(21) = group_names;
      obs_list(22) = main_effect_1;
      obs_list(23) = main_effect_2;
      obs_list(24) = 1.0;
      obs_list(25) = effects_names;
      obs_list(26) = obs_group1;
      obs_list(27) = dud_zi;
      obs_list(28) = obs_indcova1;
      obs_list(29) = obs_indcovb1;
      obs_list(30) = obs_indcovc1;
      obs_list(31) = dud_zi;
      obs_list(32) = dud_zi;
      obs_list(33) = dud_zi;
      
      List sizea_list(34);
      sizea_list(0) = sizea_num;
      sizea_list(1) = sizea_year;
      sizea_list(2) = sizea_patch;
      sizea_list(3) = sizea_group2;
      sizea_list(4) = sizea_dist;
      sizea_list(5) = zi_yn;
      sizea_list(6) = sizea_zi;
      sizea_list(7) = year_sizea_zi;
      sizea_list(8) = patch_sizea_zi;
      sizea_list(9) = group2_sizea_zi;
      sizea_list(10) = indcova_names;
      sizea_list(11) = sizea_indcova2;
      sizea_list(12) = sizea_indcova2_zi;
      sizea_list(13) = indcovb_names;
      sizea_list(14) = sizea_indcovb2;
      sizea_list(15) = sizea_indcovb2_zi;
      sizea_list(16) = indcovc_names;
      sizea_list(17) = sizea_indcovc2;
      sizea_list(18) = sizea_indcovc2_zi;
      sizea_list(19) = year_names;
      sizea_list(20) = patch_names;
      sizea_list(21) = group_names;
      sizea_list(22) = main_effect_1;
      sizea_list(23) = main_effect_2;
      sizea_list(24) = sizea_st;
      sizea_list(25) = effects_names;
      sizea_list(26) = sizea_group1;
      sizea_list(27) = group1_sizea_zi;
      sizea_list(28) = sizea_indcova1;
      sizea_list(29) = sizea_indcovb1;
      sizea_list(30) = sizea_indcovc1;
      sizea_list(31) = sizea_indcova1_zi;
      sizea_list(32) = sizea_indcovb1_zi;
      sizea_list(33) = sizea_indcovc1_zi;
      
      List sizeb_list(34);
      sizeb_list(0) = sizeb_num;
      sizeb_list(1) = sizeb_year;
      sizeb_list(2) = sizeb_patch;
      sizeb_list(3) = sizeb_group2;
      sizeb_list(4) = sizeb_dist;
      sizeb_list(5) = zi_yn;
      sizeb_list(6) = sizeb_zi;
      sizeb_list(7) = year_sizeb_zi;
      sizeb_list(8) = patch_sizeb_zi;
      sizeb_list(9) = group2_sizeb_zi;
      sizeb_list(10) = indcova_names;
      sizeb_list(11) = sizeb_indcova2;
      sizeb_list(12) = sizeb_indcova2_zi;
      sizeb_list(13) = indcovb_names;
      sizeb_list(14) = sizeb_indcovb2;
      sizeb_list(15) = sizeb_indcovb2_zi;
      sizeb_list(16) = indcovc_names;
      sizeb_list(17) = sizeb_indcovc2;
      sizeb_list(18) = sizeb_indcovc2_zi;
      sizeb_list(19) = year_names;
      sizeb_list(20) = patch_names;
      sizeb_list(21) = group_names;
      sizeb_list(22) = main_effect_1;
      sizeb_list(23) = main_effect_2;
      sizeb_list(24) = sizeb_st;
      sizeb_list(25) = effects_names;
      sizeb_list(26) = sizeb_group1;
      sizeb_list(27) = group1_sizeb_zi;
      sizeb_list(28) = sizeb_indcova1;
      sizeb_list(29) = sizeb_indcovb1;
      sizeb_list(30) = sizeb_indcovc1;
      sizeb_list(31) = sizeb_indcova1_zi;
      sizeb_list(32) = sizeb_indcovb1_zi;
      sizeb_list(33) = sizeb_indcovc1_zi;
      
      List sizec_list(34);
      sizec_list(0) = sizec_num;
      sizec_list(1) = sizec_year;
      sizec_list(2) = sizec_patch;
      sizec_list(3) = sizec_group2;
      sizec_list(4) = sizec_dist;
      sizec_list(5) = zi_yn;
      sizec_list(6) = sizec_zi;
      sizec_list(7) = year_sizec_zi;
      sizec_list(8) = patch_sizec_zi;
      sizec_list(9) = group2_sizec_zi;
      sizec_list(10) = indcova_names;
      sizec_list(11) = sizec_indcova2;
      sizec_list(12) = sizec_indcova2_zi;
      sizec_list(13) = indcovb_names;
      sizec_list(14) = sizec_indcovb2;
      sizec_list(15) = sizec_indcovb2_zi;
      sizec_list(16) = indcovc_names;
      sizec_list(17) = sizec_indcovc2;
      sizec_list(18) = sizec_indcovc2_zi;
      sizec_list(19) = year_names;
      sizec_list(20) = patch_names;
      sizec_list(21) = group_names;
      sizec_list(22) = main_effect_1;
      sizec_list(23) = main_effect_2;
      sizec_list(24) = sizec_st;
      sizec_list(25) = effects_names;
      sizec_list(26) = sizec_group1;
      sizec_list(27) = group1_sizec_zi;
      sizec_list(28) = sizec_indcova1;
      sizec_list(29) = sizec_indcovb1;
      sizec_list(30) = sizec_indcovc1;
      sizec_list(31) = sizec_indcova1_zi;
      sizec_list(32) = sizec_indcovb1_zi;
      sizec_list(33) = sizec_indcovc1_zi;
      
      List repst_list(34);
      repst_list(0) = repst_num;
      repst_list(1) = repst_year;
      repst_list(2) = repst_patch;
      repst_list(3) = repst_group2;
      repst_list(4) = repst_dist;
      repst_list(5) = false;
      repst_list(6) = dud_zi;
      repst_list(7) = dud_zi;
      repst_list(8) = dud_zi;
      repst_list(9) = dud_zi;
      repst_list(10) = indcova_names;
      repst_list(11) = repst_indcova2;
      repst_list(12) = dud_zi;
      repst_list(13) = indcovb_names;
      repst_list(14) = repst_indcovb2;
      repst_list(15) = dud_zi;
      repst_list(16) = indcovc_names;
      repst_list(17) = repst_indcovc2;
      repst_list(18) = dud_zi;
      repst_list(19) = year_names;
      repst_list(20) = patch_names;
      repst_list(21) = group_names;
      repst_list(22) = main_effect_1;
      repst_list(23) = main_effect_2;
      repst_list(24) = 1.0;
      repst_list(25) = effects_names;
      repst_list(26) = repst_group1;
      repst_list(27) = dud_zi;
      repst_list(28) = repst_indcova1;
      repst_list(29) = repst_indcovb1;
      repst_list(30) = repst_indcovc1;
      repst_list(31) = dud_zi;
      repst_list(32) = dud_zi;
      repst_list(33) = dud_zi;
      
      List fec_list(34);
      fec_list(0) = fec_num;
      fec_list(1) = fec_year;
      fec_list(2) = fec_patch;
      fec_list(3) = fec_group2;
      fec_list(4) = fec_dist;
      fec_list(5) = zi_yn;
      fec_list(6) = fec_zi;
      fec_list(7) = year_fec_zi;
      fec_list(8) = patch_fec_zi;
      fec_list(9) = group2_fec_zi;
      fec_list(10) = indcova_names;
      fec_list(11) = fec_indcova2;
      fec_list(12) = fec_indcova2_zi;
      fec_list(13) = indcovb_names;
      fec_list(14) = fec_indcovb2;
      fec_list(15) = fec_indcovb2_zi;
      fec_list(16) = indcovc_names;
      fec_list(17) = fec_indcovc2;
      fec_list(18) = fec_indcovc2_zi;
      fec_list(19) = year_names;
      fec_list(20) = patch_names;
      fec_list(21) = group_names;
      fec_list(22) = main_effect_1;
      fec_list(23) = main_effect_2;
      fec_list(24) = fec_st;
      fec_list(25) = effects_names;
      fec_list(26) = fec_group1;
      fec_list(27) = group1_fec_zi;
      fec_list(28) = fec_indcova1;
      fec_list(29) = fec_indcovb1;
      fec_list(30) = fec_indcovc1;
      fec_list(31) = fec_indcova1_zi;
      fec_list(32) = fec_indcovb1_zi;
      fec_list(33) = fec_indcovc1_zi;
      
      List jsurv_list(34);
      jsurv_list(0) = jsurv_num;
      jsurv_list(1) = jsurv_year;
      jsurv_list(2) = jsurv_patch;
      jsurv_list(3) = jsurv_group2;
      jsurv_list(4) = jsurv_dist;
      jsurv_list(5) = false;
      jsurv_list(6) = dud_zi;
      jsurv_list(7) = dud_zi;
      jsurv_list(8) = dud_zi;
      jsurv_list(9) = dud_zi;
      jsurv_list(10) = indcova_names;
      jsurv_list(11) = jsurv_indcova2;
      jsurv_list(12) = dud_zi;
      jsurv_list(13) = indcovb_names;
      jsurv_list(14) = jsurv_indcovb2;
      jsurv_list(15) = dud_zi;
      jsurv_list(16) = indcovc_names;
      jsurv_list(17) = jsurv_indcovc2;
      jsurv_list(18) = dud_zi;
      jsurv_list(19) = year_names;
      jsurv_list(20) = patch_names;
      jsurv_list(21) = group_names;
      jsurv_list(22) = main_effect_1;
      jsurv_list(23) = main_effect_2;
      jsurv_list(24) = 1.0;
      jsurv_list(25) = effects_names;
      jsurv_list(26) = jsurv_group1;
      jsurv_list(27) = dud_zi;
      jsurv_list(28) = jsurv_indcova1;
      jsurv_list(29) = jsurv_indcovb1;
      jsurv_list(30) = jsurv_indcovc1;
      jsurv_list(31) = dud_zi;
      jsurv_list(32) = dud_zi;
      jsurv_list(33) = dud_zi;
      
      List jobs_list(34);
      jobs_list(0) = jobs_num;
      jobs_list(1) = jobs_year;
      jobs_list(2) = jobs_patch;
      jobs_list(3) = jobs_group2;
      jobs_list(4) = jobs_dist;
      jobs_list(5) = false;
      jobs_list(6) = dud_zi;
      jobs_list(7) = dud_zi;
      jobs_list(8) = dud_zi;
      jobs_list(9) = dud_zi;
      jobs_list(10) = indcova_names;
      jobs_list(11) = jobs_indcova2;
      jobs_list(12) = dud_zi;
      jobs_list(13) = indcovb_names;
      jobs_list(14) = jobs_indcovb2;
      jobs_list(15) = dud_zi;
      jobs_list(16) = indcovc_names;
      jobs_list(17) = jobs_indcovc2;
      jobs_list(18) = dud_zi;
      jobs_list(19) = year_names;
      jobs_list(20) = patch_names;
      jobs_list(21) = group_names;
      jobs_list(22) = main_effect_1;
      jobs_list(23) = main_effect_2;
      jobs_list(24) = 1.0;
      jobs_list(25) = effects_names;
      jobs_list(26) = jobs_group1;
      jobs_list(27) = dud_zi;
      jobs_list(28) = jobs_indcova1;
      jobs_list(29) = jobs_indcovb1;
      jobs_list(30) = jobs_indcovc1;
      jobs_list(31) = dud_zi;
      jobs_list(32) = dud_zi;
      jobs_list(33) = dud_zi;
      
      List jsizea_list(34);
      jsizea_list(0) = jsizea_num;
      jsizea_list(1) = jsizea_year;
      jsizea_list(2) = jsizea_patch;
      jsizea_list(3) = jsizea_group2;
      jsizea_list(4) = jsizea_dist;
      jsizea_list(5) = zi_yn;
      jsizea_list(6) = jsizea_zi;
      jsizea_list(7) = year_jsizea_zi;
      jsizea_list(8) = patch_jsizea_zi;
      jsizea_list(9) = group2_jsizea_zi;
      jsizea_list(10) = indcova_names;
      jsizea_list(11) = jsizea_indcova2;
      jsizea_list(12) = jsizea_indcova2_zi;
      jsizea_list(13) = indcovb_names;
      jsizea_list(14) = jsizea_indcovb2;
      jsizea_list(15) = jsizea_indcovb2_zi;
      jsizea_list(16) = indcovc_names;
      jsizea_list(17) = jsizea_indcovc2;
      jsizea_list(18) = jsizea_indcovc2_zi;
      jsizea_list(19) = year_names;
      jsizea_list(20) = patch_names;
      jsizea_list(21) = group_names;
      jsizea_list(22) = main_effect_1;
      jsizea_list(23) = main_effect_2;
      jsizea_list(24) = jsizea_st;
      jsizea_list(25) = effects_names;
      jsizea_list(26) = jsizea_group1;
      jsizea_list(27) = group1_jsizea_zi;
      jsizea_list(28) = jsizea_indcova1;
      jsizea_list(29) = jsizea_indcovb1;
      jsizea_list(30) = jsizea_indcovc1;
      jsizea_list(31) = jsizea_indcova1_zi;
      jsizea_list(32) = jsizea_indcovb1_zi;
      jsizea_list(33) = jsizea_indcovc1_zi;
      
      List jsizeb_list(34);
      jsizeb_list(0) = jsizeb_num;
      jsizeb_list(1) = jsizeb_year;
      jsizeb_list(2) = jsizeb_patch;
      jsizeb_list(3) = jsizeb_group2;
      jsizeb_list(4) = jsizeb_dist;
      jsizeb_list(5) = zi_yn;
      jsizeb_list(6) = jsizeb_zi;
      jsizeb_list(7) = year_jsizeb_zi;
      jsizeb_list(8) = patch_jsizeb_zi;
      jsizeb_list(9) = group2_jsizeb_zi;
      jsizeb_list(10) = indcova_names;
      jsizeb_list(11) = jsizeb_indcova2;
      jsizeb_list(12) = jsizeb_indcova2_zi;
      jsizeb_list(13) = indcovb_names;
      jsizeb_list(14) = jsizeb_indcovb2;
      jsizeb_list(15) = jsizeb_indcovb2_zi;
      jsizeb_list(16) = indcovc_names;
      jsizeb_list(17) = jsizeb_indcovc2;
      jsizeb_list(18) = jsizeb_indcovc2_zi;
      jsizeb_list(19) = year_names;
      jsizeb_list(20) = patch_names;
      jsizeb_list(21) = group_names;
      jsizeb_list(22) = main_effect_1;
      jsizeb_list(23) = main_effect_2;
      jsizeb_list(24) = jsizeb_st;
      jsizeb_list(25) = effects_names;
      jsizeb_list(26) = jsizeb_group1;
      jsizeb_list(27) = group1_jsizeb_zi;
      jsizeb_list(28) = jsizeb_indcova1;
      jsizeb_list(29) = jsizeb_indcovb1;
      jsizeb_list(30) = jsizeb_indcovc1;
      jsizeb_list(31) = jsizeb_indcova1_zi;
      jsizeb_list(32) = jsizeb_indcovb1_zi;
      jsizeb_list(33) = jsizeb_indcovc1_zi;
      
      List jsizec_list(34);
      jsizec_list(0) = jsizec_num;
      jsizec_list(1) = jsizec_year;
      jsizec_list(2) = jsizec_patch;
      jsizec_list(3) = jsizec_group2;
      jsizec_list(4) = jsizec_dist;
      jsizec_list(5) = zi_yn;
      jsizec_list(6) = jsizec_zi;
      jsizec_list(7) = year_jsizec_zi;
      jsizec_list(8) = patch_jsizec_zi;
      jsizec_list(9) = group2_jsizec_zi;
      jsizec_list(10) = indcova_names;
      jsizec_list(11) = jsizec_indcova2;
      jsizec_list(12) = jsizec_indcova2_zi;
      jsizec_list(13) = indcovb_names;
      jsizec_list(14) = jsizec_indcovb2;
      jsizec_list(15) = jsizec_indcovb2_zi;
      jsizec_list(16) = indcovc_names;
      jsizec_list(17) = jsizec_indcovc2;
      jsizec_list(18) = jsizec_indcovc2_zi;
      jsizec_list(19) = year_names;
      jsizec_list(20) = patch_names;
      jsizec_list(21) = group_names;
      jsizec_list(22) = main_effect_1;
      jsizec_list(23) = main_effect_2;
      jsizec_list(24) = jsizec_st;
      jsizec_list(25) = effects_names;
      jsizec_list(26) = jsizec_group1;
      jsizec_list(27) = group1_jsizec_zi;
      jsizec_list(28) = jsizec_indcova1;
      jsizec_list(29) = jsizec_indcovb1;
      jsizec_list(30) = jsizec_indcovc1;
      jsizec_list(31) = jsizec_indcova1_zi;
      jsizec_list(32) = jsizec_indcovb1_zi;
      jsizec_list(33) = jsizec_indcovc1_zi;
      
      List jrepst_list(34);
      jrepst_list(0) = jrepst_num;
      jrepst_list(1) = jrepst_year;
      jrepst_list(2) = jrepst_patch;
      jrepst_list(3) = jrepst_group2;
      jrepst_list(4) = jrepst_dist;
      jrepst_list(5) = false;
      jrepst_list(6) = dud_zi;
      jrepst_list(7) = dud_zi;
      jrepst_list(8) = dud_zi;
      jrepst_list(9) = dud_zi;
      jrepst_list(10) = indcova_names;
      jrepst_list(11) = jrepst_indcova2;
      jrepst_list(12) = dud_zi;
      jrepst_list(13) = indcovb_names;
      jrepst_list(14) = jrepst_indcovb2;
      jrepst_list(15) = dud_zi;
      jrepst_list(16) = indcovc_names;
      jrepst_list(17) = jrepst_indcovc2;
      jrepst_list(18) = dud_zi;
      jrepst_list(19) = year_names;
      jrepst_list(20) = patch_names;
      jrepst_list(21) = group_names;
      jrepst_list(22) = main_effect_1;
      jrepst_list(23) = main_effect_2;
      jrepst_list(24) = 1.0;
      jrepst_list(25) = effects_names;
      jrepst_list(26) = jrepst_group1;
      jrepst_list(27) = dud_zi;
      jrepst_list(28) = jrepst_indcova1;
      jrepst_list(29) = jrepst_indcovb1;
      jrepst_list(30) = jrepst_indcovc1;
      jrepst_list(31) = dud_zi;
      jrepst_list(32) = dud_zi;
      jrepst_list(33) = dud_zi;
      
      List jmatst_list(34);
      jmatst_list(0) = jmatst_num;
      jmatst_list(1) = jmatst_year;
      jmatst_list(2) = jmatst_patch;
      jmatst_list(3) = jmatst_group2;
      jmatst_list(4) = jmatst_dist;
      jmatst_list(5) = false;
      jmatst_list(6) = dud_zi;
      jmatst_list(7) = dud_zi;
      jmatst_list(8) = dud_zi;
      jmatst_list(9) = dud_zi;
      jmatst_list(10) = indcova_names;
      jmatst_list(11) = jmatst_indcova2;
      jmatst_list(12) = dud_zi;
      jmatst_list(13) = indcovb_names;
      jmatst_list(14) = jmatst_indcovb2;
      jmatst_list(15) = dud_zi;
      jmatst_list(16) = indcovc_names;
      jmatst_list(17) = jmatst_indcovc2;
      jmatst_list(18) = dud_zi;
      jmatst_list(19) = year_names;
      jmatst_list(20) = patch_names;
      jmatst_list(21) = group_names;
      jmatst_list(22) = main_effect_1;
      jmatst_list(23) = main_effect_2;
      jmatst_list(24) = 1.0;
      jmatst_list(25) = effects_names;
      jmatst_list(26) = jmatst_group1;
      jmatst_list(27) = dud_zi;
      jmatst_list(28) = jmatst_indcova1;
      jmatst_list(29) = jmatst_indcovb1;
      jmatst_list(30) = jmatst_indcovc1;
      jmatst_list(31) = dud_zi;
      jmatst_list(32) = dud_zi;
      jmatst_list(33) = dud_zi;
      
      surv_list.attr("names") = list_names;
      obs_list.attr("names") = list_names;
      sizea_list.attr("names") = list_names;
      sizeb_list.attr("names") = list_names;
      sizec_list.attr("names") = list_names;
      repst_list.attr("names") = list_names;
      fec_list.attr("names") = list_names;
      jsurv_list.attr("names") = list_names;
      jobs_list.attr("names") = list_names;
      jsizea_list.attr("names") = list_names;
      jsizeb_list.attr("names") = list_names;
      jsizec_list.attr("names") = list_names;
      jrepst_list.attr("names") = list_names;
      jmatst_list.attr("names") = list_names;
      
      surmodl = surv_list;
      obsmodl = obs_list;
      sizmodl = sizea_list;
      sibmodl = sizeb_list;
      sicmodl = sizec_list;
      repmodl = repst_list;
      fecmodl = fec_list;
      
      jsurmodl = jsurv_list;
      jobsmodl = jobs_list;
      jsizmodl = jsizea_list;
      jsibmodl = jsizeb_list;
      jsicmodl = jsizec_list;
      jrepmodl = jrepst_list;
      jmatmodl = jmatst_list;
      
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
      CharacterVector mainparams = {"year2", "individ", "patch", "surv3",
        "obs3", "size3", "sizeb3", "sizec3", "repst3", "fec3", "fec2", "size2",
        "size1", "sizeb2", "sizeb1", "sizec2", "sizec1", "repst2", "repst1",
        "matst3", "matst2", "age", "density", "indcova2", "indcova1",
        "indcovb2", "indcovb1", "indcovc2", "indcovc1", "group2", "group1",
        "annucova2", "annucova1", "annucovb2", "annucovb1", "annucovc2",
        "annucovc1"};
      CharacterVector modelparams = clone(mainparams);
      
      DataFrame pm_names = DataFrame::create(_["parameter_names"] = parameter_names,
        _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
      pmnames = pm_names;
      pmn_provided = true;
    } else {
      String eat_my_shorts = "Argument modelsuite must be set to a lefkoMod ";
      eat_my_shorts += "or vrm_input object, or vital rate models must be ";
      eat_my_shorts += "supplied with modelsuite left unentered.";
      
      throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
    }
  } else {
    if (surv_model.isNotNull()) {
      RObject sum_intermediate = RObject(surv_model);
      surmodl = sum_intermediate;
      modelcheck++;
    } else {
      surmodl = clone(model1);
    }
    if (obs_model.isNotNull()) {
      RObject obm_intermediate = RObject(obs_model);
      obsmodl = obm_intermediate;
      modelcheck++;
    } else {
      obsmodl = clone(model1);
    }
    if (size_model.isNotNull()) {
      RObject sim_intermediate = RObject(size_model);
      sizmodl = sim_intermediate;
      modelcheck++;
    } else {
      sizmodl = clone(model1);
    }
    if (sizeb_model.isNotNull()) {
      RObject sbm_intermediate = RObject(sizeb_model);
      sibmodl = sbm_intermediate;
      modelcheck++;
    } else {
      sibmodl = clone(model1);
    }
    if (sizec_model.isNotNull()) {
      RObject scm_intermediate = RObject(sizec_model);
      sicmodl = scm_intermediate;
      modelcheck++;
    } else {
      sicmodl = clone(model1);
    }
    if (repst_model.isNotNull()) {
      RObject rpm_intermediate = RObject(repst_model);
      repmodl = rpm_intermediate;
      modelcheck++;
    } else {
      repmodl = clone(model1);
    }
    if (fec_model.isNotNull()) {
      RObject fem_intermediate = RObject(fec_model);
      fecmodl = fem_intermediate;
      modelcheck++;
    } else {
      fecmodl = clone(model1);
    }
    
    if (jsurv_model.isNotNull()) {
      RObject jsum_intermediate = RObject(jsurv_model);
      jsurmodl = jsum_intermediate;
      modelcheck++;
    } else {
      jsurmodl = clone(model1);
    }
    if (jobs_model.isNotNull()) {
      RObject jobm_intermediate = RObject(jobs_model);
      jobsmodl = jobm_intermediate;
      modelcheck++;
    } else {
      jobsmodl = clone(model1);
    }
    if (jsize_model.isNotNull()) {
      RObject jsim_intermediate = RObject(jsize_model);
      jsizmodl = jsim_intermediate;
      modelcheck++;
    } else {
      jsizmodl = clone(model1);
    }
    if (jsizeb_model.isNotNull()) {
      RObject jsbm_intermediate = RObject(jsizeb_model);
      jsibmodl = jsbm_intermediate;
      modelcheck++;
    } else {
      jsibmodl = clone(model1);
    }
    if (jsizec_model.isNotNull()) {
      RObject jscm_intermediate = RObject(jsizec_model);
      jsicmodl = jscm_intermediate;
      modelcheck++;
    } else {
      jsicmodl = clone(model1);
    }
    if (jrepst_model.isNotNull()) {
      RObject jrpm_intermediate = RObject(jrepst_model);
      jrepmodl = jrpm_intermediate;
      modelcheck++;
    } else {
      jrepmodl = clone(model1);
    }
    if (jmatst_model.isNotNull()) {
      RObject jmat_intermediate = RObject(jmatst_model);
      jmatmodl = jmat_intermediate;
      modelcheck++;
    } else {
      jmatmodl = clone(model1);
    }
    
    if (paramnames.isNotNull()) {
      DataFrame pmn_intermediate = DataFrame(paramnames);
      pmnames = pmn_intermediate;
      
      pmn_provided = true;
    }
  }
  
  if (modelcheck == 0 && !nodata) {
    String eat_my_shorts = "Argument modelsuite must be set to a lefkoMod ";
    eat_my_shorts += "or vrm_input object, or vital rate models must be ";
    eat_my_shorts += "supplied with modelsuite left unentered.";
    
    throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
  }
  if (!pmn_provided) pop_error("paramnames", "lefkoMod object", "", 15);
  
  // Stageframe
  DataFrame sframe;
  RObject maingroups;
  
  if (stageframe.isNotNull()) {
    DataFrame sf_intermediate = DataFrame(stageframe);
    sframe = sf_intermediate;
    
    CharacterVector maingroups_int = as<CharacterVector>(sframe["group"]);
    CharacterVector maingroups_int2 = unique(maingroups_int);
    
    CharacterVector maingroups_sorted = stringsort(maingroups_int2);
    maingroups = as<RObject>(maingroups_sorted);
    
  } else if (format < 5) {
    throw Rcpp::exception("A stageframe is required for all MPM formats except Leslie MPMs.", false);
  } else {
    sframe = R_NilValue;
    
    CharacterVector maingroups_sorted = {"0"};
    maingroups = as<RObject>(maingroups_sorted);
  }
  
  // Demographic data
  int used_yearcol {-1};
  int used_patchcol {-1};
  int used_agecol {-1};
  int used_indacol {-1};
  int used_indbcol {-1};
  int used_indccol {-1};
  
  CharacterVector modelparam_names = pmnames["modelparams"];
  std::string year2var = as<std::string>(modelparam_names(0));
  std::string patchvar = as<std::string>(modelparam_names(2));
  std::string agevar = as<std::string>(modelparam_names(21));
  std::string indcova2var = as<std::string>(modelparam_names(23));
  std::string indcova1var = as<std::string>(modelparam_names(24));
  std::string indcovb2var = as<std::string>(modelparam_names(25));
  std::string indcovb1var = as<std::string>(modelparam_names(26));
  std::string indcovc2var = as<std::string>(modelparam_names(27));
  std::string indcovc1var = as<std::string>(modelparam_names(28));
  std::string group2var = as<std::string>(modelparam_names(29));
  std::string group1var = as<std::string>(modelparam_names(30));
  
  IntegerVector mainyears_int;
  CharacterVector mainyears;
  CharacterVector mainpatches;
  DataFrame true_data;
  
  if (data.isNotNull() && !nodata) {
    true_data = DataFrame(data);
    int no_vars = true_data.length();
    CharacterVector data_names = true_data.attr("names");
    
    for (int i = 0; i < no_vars; i++) {
      if (stringcompare_hard(as<std::string>(data_names(i)), year2var)) {
        used_yearcol = i;
      }
      if (stringcompare_hard(as<std::string>(data_names(i)), patchvar)) {
        used_patchcol = i;
      }
      if (stringcompare_hard(as<std::string>(data_names(i)), agevar)) {
        used_agecol = i;
      }
      if (stringcompare_hard(as<std::string>(data_names(i)), indcova2var)) {
        used_indacol = i;
      }
      if (stringcompare_hard(as<std::string>(data_names(i)), indcovb2var)) {
        used_indbcol = i;
      }
      if (stringcompare_hard(as<std::string>(data_names(i)), indcovc2var)) {
        used_indccol = i;
      }
    }
    if (used_yearcol < 0 || used_yearcol > (no_vars - 1)) pop_error("time t", "", "", 16);
    CharacterVector all_years = as<CharacterVector>(true_data[used_yearcol]);
    mainyears = sort_unique(all_years);
    mainyears_int = seq_along(mainyears);
    
    if (used_patchcol > -1) {
      CharacterVector all_patches = as<CharacterVector>(true_data[used_patchcol]);
      CharacterVector mainpatches_int = unique(all_patches);
      
      mainpatches = stringsort(mainpatches_int);
    } else {
      mainpatches = {NA_STRING};
    }
  
  } else if (!data.isNotNull() && !nodata) {
    pop_error("data", "a lefkoMod", "", 11);
      
  } else {
    DataFrame year_frame = as<DataFrame>(msuite["year_frame"]);
    
    CharacterVector year_names = as<CharacterVector>(year_frame["years"]);
    mainyears = year_names;
    mainyears_int = seq_along(mainyears);
    
    DataFrame patch_frame = as<DataFrame>(msuite["patch_frame"]);
    CharacterVector all_patches = as<CharacterVector>(patch_frame["patches"]);
    CharacterVector mainpatches_int = unique(all_patches);
    mainpatches = stringsort(mainpatches_int);
  }
  
  IntegerVector mainages;
  IntegerVector actualages;
  int age_limit {0};
  
  if (format > 3) { // Age and age-by-stage cases
    if (IntegerVector::is_na(start_age)) {
      if (prebreeding) {
        start_age = 1;
      } else {
        start_age = 0;
      }
    }
    
    if (data.isNotNull() && !nodata) {
      if (used_agecol > -1) {
        IntegerVector all_ages = as<IntegerVector>(true_data[used_agecol]);
        IntegerVector mainages_pre = unique(all_ages);
        mainages = int_sort(mainages_pre);
        
      }
    } else {
      if (IntegerVector::is_na(last_age)) pop_error("last_age", "a vrm_input", "", 11);
      
      if (last_age <= start_age) {
        throw Rcpp::exception("Argument last_age must be greater than start_age.", false);
      }
      mainages = seq(start_age, last_age);
    }
    
    age_limit = max(mainages) + 1;
    
    if (IntegerVector::is_na(last_age)) {
      last_age = max(mainages) + 1;
    }
    if (IntegerVector::is_na(fecage_min)) {
      fecage_min = min(mainages);
    }
    if (IntegerVector::is_na(fecage_max)) {
      fecage_max = last_age;
    }
    
    if (start_age > age_limit || last_age > age_limit) {
      if (!quiet) Rf_warningcall(R_NilValue,
        "Entered start_age or last_age is beyond what is found in the dataset.");
    }
    if (fecage_min > age_limit || fecage_max > age_limit) {
      if (!quiet) Rf_warningcall(R_NilValue,
        "Entered fecage_min or fecage_max is beyond what is found in the dataset.");
    }
    
    if (last_age < (start_age + 1)) {
      throw Rcpp::exception("Argument last_age must be greater than start_age.", false);
    }
    if (fecage_max < fecage_min) {
      throw Rcpp::exception("Argument fecage_max must be greater than or equal to fecage_min.", false);
    }
  }
  
  // Ind_names objects
  RObject inda_names;
  RObject indb_names;
  RObject indc_names;
  CharacterVector inda_names_ch;
  CharacterVector indb_names_ch;
  CharacterVector indc_names_ch;
  
  if (!nodata) {
    if (used_indacol > -1) {
      if (random_inda) {
        CharacterVector inda_data = as<CharacterVector>(true_data[used_indacol]);
        CharacterVector unique_a = unique(inda_data);
        CharacterVector sorted_a = stringsort(unique_a);
        inda_names_ch = sorted_a;
        
        inda_names = as<RObject>(sorted_a);
        
      } else {
        IntegerVector notrandom_a = {0};
        inda_names = as<RObject>(notrandom_a);
      }
    } else {
      IntegerVector notrandom_a = {0};
      inda_names = as<RObject>(notrandom_a);
      random_inda = false;
    }
    
    if (used_indbcol > -1) {
      if (random_indb) {
        CharacterVector indb_data = as<CharacterVector>(true_data[used_indbcol]);
        CharacterVector unique_b = unique(indb_data);
        CharacterVector sorted_b = stringsort(unique_b);
        indb_names_ch = sorted_b;
        
        indb_names = as<RObject>(sorted_b);
        
      } else {
        IntegerVector notrandom_b = {0};
        indb_names = as<RObject>(notrandom_b);
      }
    } else {
      IntegerVector notrandom_b = {0};
      indb_names = as<RObject>(notrandom_b);
      random_indb = false;
    }
    
    if (used_indccol > -1) {
      if (random_indc) {
        CharacterVector indc_data = as<CharacterVector>(true_data[used_indccol]);
        CharacterVector unique_c = unique(indc_data);
        CharacterVector sorted_c = stringsort(unique_c);
        indc_names_ch = sorted_c;
        
        indc_names = as<RObject>(sorted_c);
        
      } else {
        IntegerVector notrandom_c = {0};
        indc_names = as<RObject>(notrandom_c);
      }
    } else {
      IntegerVector notrandom_c = {0};
      indc_names = as<RObject>(notrandom_c);
      random_indc = false;
    }
  } else {
    CharacterVector ms_names = msuite.attr("names");
    bool indcova_inc = false;
    bool indcovb_inc = false;
    bool indcovc_inc = false;
    
    for (int i = 0; i < ms_names.length(); i++) {
      if (stringcompare_simple(as<std::string>(ms_names(i)), "indcova2_frame", false)) {
        indcova_inc = true;
      }
      if (stringcompare_simple(as<std::string>(ms_names(i)), "indcovb2_frame", false)) {
        indcovb_inc = true;
      }
      if (stringcompare_simple(as<std::string>(ms_names(i)), "indcovc2_frame", false)) {
        indcovc_inc = true;
      }
    }
    
    if (indcova_inc) {
      DataFrame ica_frame = as<DataFrame>(msuite["indcova2_frame"]);
      CharacterVector sorted_a = as<CharacterVector>(ica_frame["indcova"]);
      
      inda_names_ch = sorted_a;
      
      inda_names = as<RObject>(sorted_a);
    } else {
      IntegerVector notrandom_a = {0};
      inda_names = as<RObject>(notrandom_a);
      random_inda = false;
    }
    
    if (indcovb_inc) {
      DataFrame icb_frame = as<DataFrame>(msuite["indcovb2_frame"]);
      CharacterVector sorted_b = as<CharacterVector>(icb_frame["indcovb"]);
      
      indb_names_ch = sorted_b;
      
      indb_names = as<RObject>(sorted_b);
    } else {
      IntegerVector notrandom_b = {0};
      indb_names = as<RObject>(notrandom_b);
      random_indb = false;
    }
    
    if (indcovc_inc) {
      DataFrame icc_frame = as<DataFrame>(msuite["indcovc2_frame"]);
      CharacterVector sorted_c = as<CharacterVector>(icc_frame["indcovc"]);
      
      indc_names_ch = sorted_c;
      
      indc_names = as<RObject>(sorted_c);
    } else {
      IntegerVector notrandom_c = {0};
      indc_names = as<RObject>(notrandom_c);
      random_indc = false;
    }
  }
  
  // Check years and patches entries
  List all_years_topull(nreps);
  CharacterVector years_topull;
  CharacterMatrix years_projected (times, nreps);
  
  int num_years = static_cast<int>(mainyears.length());
  if (year.isNotNull()) {
    CharacterVector year_char (year);
    
    CharacterVector years_unmatched = setdiff(year_char, mainyears);
    if (static_cast<int>(years_unmatched.length()) > 0) {
      throw Rcpp::exception("Some input year values do not match dataset.", false);
    }
    if (stochastic) {
      String eat_my_shorts = "Projection cannot be stochastic if argument year is set.";
      throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
    }
    years_topull = year_char;
  } else if (!stochastic) {
    CharacterVector year_char = clone(mainyears);
    years_topull = year_char;
    if (!quiet) Rf_warningcall(R_NilValue,
      "Option year not set, so will cycle through existing years.");
  }
  
  // New code, with change to arma
  arma::vec twinput;
  arma::mat twinput_markov;
  if (tweights.isNotNull()) {
    if (Rf_isMatrix(tweights)) {
      twinput_markov = as<arma::mat>(tweights);
      assume_markov = true;
      
      if (static_cast<int>(twinput_markov.n_cols) != num_years) {
        pop_error("tweights", "occasions", "dataset", 20);
      }
      if (twinput_markov.n_cols != twinput_markov.n_rows) {
        pop_error("tweights", "", "", 4);
      }
      
    } else if (is<NumericVector>(tweights)) {
      twinput = as<arma::vec>(tweights);
      
      if (static_cast<int>(twinput.n_elem) != num_years) pop_error("tweights", "occasions", "dataset", 19);
      
    } else pop_error("tweights", "a numeric vector or matrix", "", 1);
    
    if (!stochastic) {
      throw Rcpp::exception("Argument tweights can only be used when stochastic = TRUE.",
        false);
    }
  } else {
    twinput.resize(num_years);
    twinput.ones();
  }
  
  if (stochastic) {
    if (stochastic && !assume_markov) {
      years_topull = Rcpp::RcppArmadillo::sample(mainyears, times * nreps, true, twinput);
    } else if (stochastic && assume_markov) {
      CharacterVector ytp (times * nreps);
      
      twinput = twinput_markov.col(0);
      for (int yr_counter = 0; yr_counter < (times * nreps); yr_counter++) {
        twinput = twinput / sum(twinput);
        CharacterVector ytp_piecemeal = Rcpp::RcppArmadillo::sample(mainyears, 1,
          true, twinput);
        ytp(yr_counter) = ytp_piecemeal(0);
        
        int mnyrs_preassigned {0};
        for (int i = 0; i < num_years; i++) {
          if (ytp_piecemeal(0) == mainyears(i)) mnyrs_preassigned = i;
        }
        twinput = twinput_markov.col(mnyrs_preassigned);
      }
      years_topull = ytp;
    }
  } else {
    CharacterVector true_years_topull (times * nreps);
    int current_year_counter = 0;
    int years_topull_length = static_cast<int>(years_topull.length());
    for (int i = 0; i < (times * nreps); i++) {
      true_years_topull(i) = years_topull(current_year_counter);
      current_year_counter++;
      
      if (current_year_counter == years_topull_length) current_year_counter = 0;
    }
    years_topull = true_years_topull;
  }
  
  CharacterVector patches_topull;
  CharacterVector patches_projected (times);
  int chosenpatch {1};
  if (patch.isNotNull()) {
    CharacterVector patch_int (patch);
    
    CharacterVector patches_unmatched = setdiff(patch_int, mainpatches);
    if (patches_unmatched.length() > 0) pop_error("patch", "dataset", "", 17);
    patches_topull = patch_int;
    int crazy_patch {-1};
    if (patch_int.length() == 1) {
      for (int i = 0; i < mainpatches.length(); i++) {
        if (stringcompare_hard(as<std::string>(patch_int(0)), as<std::string>(mainpatches(i)))) {
          crazy_patch = i;
        }
      }
      if (crazy_patch != -1) chosenpatch = crazy_patch + 1;
    }
  } else {
    CharacterVector patch_one (1);
    patch_one(0) = mainpatches(0);
    patches_topull = patch_one;
    if (!quiet) {
      Rf_warningcall(R_NilValue, "Option patch not set, so will set to first patch/population.");
    }
  }
  
  // Handle spatial density vector
  NumericVector spdensity_topull;
  NumericVector spdensity_projected (times);
  if (sp_density.isNotNull()) {
    NumericVector spdensity_int (sp_density);
    spdensity_topull = spdensity_int;
  } else{
    NumericVector spdensity_int = {0.0};
    spdensity_topull = spdensity_int;
  }
  
  // Ind_terms data frame
  NumericVector f2_inda_values(num_years);
  NumericVector f2_indb_values(num_years);
  NumericVector f2_indc_values(num_years);
  NumericVector f1_inda_values(num_years);
  NumericVector f1_indb_values(num_years);
  NumericVector f1_indc_values(num_years);
  CharacterVector r2_inda_values(num_years);
  CharacterVector r2_indb_values(num_years);
  CharacterVector r2_indc_values(num_years);
  CharacterVector r1_inda_values(num_years);
  CharacterVector r1_indb_values(num_years);
  CharacterVector r1_indc_values(num_years);
  
  NumericVector f_inda_topull;
  NumericVector f_indb_topull;
  NumericVector f_indc_topull;
  CharacterVector r_inda_topull;
  CharacterVector r_indb_topull;
  CharacterVector r_indc_topull;
  
  bool lessthan_warning {false};
  bool greaterthan_warning {false};
  
  if (ind_terms.isNotNull()) {
    DataFrame it_ROint = RObject(ind_terms);
    
    RObject inda_whatever;
    RObject indb_whatever;
    RObject indc_whatever;
    
    if (is<DataFrame>(it_ROint)) {
      DataFrame it_intermediate = DataFrame(it_ROint);
      int it_size = static_cast<int>(it_intermediate.size());
      if (it_size != 3) pop_error("ind_terms", "a data frame or matrix with 3 columns and times rows", "", 1);
      
      inda_whatever = as<RObject>(it_intermediate[0]);
      indb_whatever = as<RObject>(it_intermediate[1]);
      indc_whatever = as<RObject>(it_intermediate[2]);
      
    } else if (is<NumericMatrix>(it_ROint)) {
      NumericMatrix it_intermediate = NumericMatrix(it_ROint);
      int it_rows = it_intermediate.nrow();
      int it_cols = it_intermediate.ncol();
      
      if (it_rows != 3 && it_cols != 3) pop_error("ind_terms",
        "a data frame or matrix with 3 columns and times rows", "", 1);
      
      if (it_cols == 3) {
        NumericVector chuck1 = it_intermediate(_, 0);
        NumericVector chuck2 = it_intermediate(_, 1);
        NumericVector chuck3 = it_intermediate(_, 2);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
      } else if (it_rows == 3) {
        NumericVector chuck1 = it_intermediate(0, _);
        NumericVector chuck2 = it_intermediate(1, _);
        NumericVector chuck3 = it_intermediate(2, _);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
      }
    } else if (is<CharacterMatrix>(it_ROint)) {
      CharacterMatrix it_intermediate = CharacterMatrix(it_ROint);
      int it_rows = it_intermediate.nrow();
      int it_cols = it_intermediate.ncol();
      
      if (it_rows != 3 && it_cols != 3) pop_error("ind_terms",
        "a data frame or matrix with 3 columns and times rows.", "", 1);
      
      if (it_cols == 3) {
        CharacterVector chuck1 = it_intermediate(_, 0);
        CharacterVector chuck2 = it_intermediate(_, 1);
        CharacterVector chuck3 = it_intermediate(_, 2);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
       } else if (it_rows == 3) {
        CharacterVector chuck1 = it_intermediate(0, _);
        CharacterVector chuck2 = it_intermediate(1, _);
        CharacterVector chuck3 = it_intermediate(2, _);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
       }
    }
    
    if (is<NumericVector>(inda_whatever) || is<IntegerVector>(inda_whatever)) {
      if (random_inda) {
        if (!quiet) {
          Rf_warningcall(R_NilValue,
            "Indcov a appears to be numeric. Will assume random_inda = FALSE.");
        }
        random_inda = false;
      }
      NumericVector inda_another_int = as<NumericVector>(inda_whatever);
      
      int check_len = inda_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_inda_topull = inda_another_int;
      r_inda_topull = {"none"};
      
    } else if (is<CharacterVector>(inda_whatever)) {
      if (!random_inda) {
        if (!quiet) {
          Rf_warningcall(R_NilValue,
            "Indcov a appears to be categorical. Will assume random_inda = TRUE.");
        }
        random_inda = true;
      }
      CharacterVector inda_another_int = as<CharacterVector>(inda_whatever);
      
      int inda_inputlength = inda_another_int.length();
      int maininda_length = inda_names_ch.length();
      
      for (int i = 0; i < inda_inputlength; i++) {
        for (int j = 0; j < maininda_length; j++) {
          if (!stringcompare_hard(as<std::string>(inda_another_int(i)), as<std::string>(inda_names_ch(j)))) {
            pop_error("indcov a", "dataset", "", 17);
          }
        }
      }
      int check_len = inda_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      r_inda_topull = inda_another_int;
      f_inda_topull = {0};
    }
    
    if (is<NumericVector>(indb_whatever) || is<IntegerVector>(indb_whatever)) {
      if (random_indb) {
        if (!quiet) {
          Rf_warningcall(R_NilValue,
            "Indcov b appears to be numeric. Will assume random_indb = FALSE.");
        }
        random_indb = false;
      }
      NumericVector indb_another_int = as<NumericVector>(indb_whatever);
      
      int check_len = indb_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_indb_topull = indb_another_int;
      r_indb_topull = {"none"};
      
    } else if (is<CharacterVector>(indb_whatever)) {
      if (!random_indb) {
        if (!quiet) {
          Rf_warningcall(R_NilValue,
            "Indcov b appears to be categorical. Will assume random_indb = TRUE.");
        }
        random_indb = true;
      }
      CharacterVector indb_another_int = as<CharacterVector>(indb_whatever);
      
      int indb_inputlength = indb_another_int.length();
      int mainindb_length = indb_names_ch.length();
      
      for (int i = 0; i < indb_inputlength; i++) {
        for (int j = 0; j < mainindb_length; j++) {
          if (!stringcompare_hard(as<std::string>(indb_another_int(i)), as<std::string>(indb_names_ch(j)))) {
            pop_error("indcov b", "dataset", "", 17);
          }
        }
      }
      int check_len = indb_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      r_indb_topull = indb_another_int;
      f_indb_topull = {0};
    }
    
    if (is<NumericVector>(indc_whatever) || is<IntegerVector>(indc_whatever)) {
      if (random_indc) {
        if (!quiet) {
          Rf_warningcall(R_NilValue,
            "Indcov c appears to be numeric. Will assume random_indc = FALSE.");
        }
        random_indc = false;
      }
      NumericVector indc_another_int = as<NumericVector>(indc_whatever);
      
      int check_len = indc_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_indc_topull = indc_another_int;
      r_indc_topull = {"none"};
      
    } else if (is<CharacterVector>(indc_whatever)) {
      if (!random_indc) {
        if (!quiet) {
          Rf_warningcall(R_NilValue,
            "Indcov c appears to be categorical. Will assume random_indc = TRUE.");
        }
        random_indc = true;
      }
      CharacterVector indc_another_int = as<CharacterVector>(indc_whatever);
      
      int indc_inputlength = indc_another_int.length();
      int mainindc_length = indc_names_ch.length();
      
      for (int i = 0; i < indc_inputlength; i++) {
        for (int j = 0; j < mainindc_length; j++) {
          if (!stringcompare_hard(as<std::string>(indc_another_int(i)), as<std::string>(indc_names_ch(j)))) {
            pop_error("indcov c", "dataset", "", 17);
          }
        }
      }
      int check_len = indc_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      r_indc_topull = indc_another_int;
      f_indc_topull = {0};
    }
    
  } else {
    f_inda_topull = clone(model0);
    r_inda_topull = clone(modelnone);
    f_indb_topull = clone(model0);
    r_indb_topull = clone(modelnone);
    f_indc_topull = clone(model0);
    r_indc_topull = clone(modelnone);
  }
  
  if (!quiet) {
    if (greaterthan_warning) {
      Rf_warningcall(R_NilValue,
        "Too many values of individual covariates have been supplied. Some will be cut.");
    }
    if (lessthan_warning) {
      Rf_warningcall(R_NilValue,
        "Too few values of individual covariates have been supplied. Some will be cycled.");
    }
  }
  
  // Ann_terms data frame
  NumericVector f2_anna_values(num_years);
  NumericVector f2_annb_values(num_years);
  NumericVector f2_annc_values(num_years);
  NumericVector f1_anna_values(num_years);
  NumericVector f1_annb_values(num_years);
  NumericVector f1_annc_values(num_years);
  
  NumericVector f_anna_topull;
  NumericVector f_annb_topull;
  NumericVector f_annc_topull;
  
  lessthan_warning = false;
  greaterthan_warning = false;
  
  if (ann_terms.isNotNull()) {
    DataFrame at_ROint = RObject(ann_terms);
    
    RObject anna_whatever;
    RObject annb_whatever;
    RObject annc_whatever;
    
    if (is<DataFrame>(at_ROint)) {
      DataFrame at_intermediate = DataFrame(at_ROint);
      int at_size = static_cast<int>(at_intermediate.size());
      if (at_size != 3) pop_error("ann_terms",
        "a data frame or matrix with 3 columns and times rows", "", 1);
      
      anna_whatever = as<RObject>(at_intermediate[0]);
      annb_whatever = as<RObject>(at_intermediate[1]);
      annc_whatever = as<RObject>(at_intermediate[2]);
      
    } else if (is<NumericMatrix>(at_ROint)) {
      NumericMatrix at_intermediate = NumericMatrix(at_ROint);
      int at_rows = at_intermediate.nrow();
      int at_cols = at_intermediate.ncol();
      
      if (at_rows != 3 && at_cols != 3) pop_error("ann_terms",
        "a data frame or matrix with 3 columns and times rows", "", 1);
      
      if (at_cols == 3) {
        NumericVector chuck1 = at_intermediate(_, 0);
        NumericVector chuck2 = at_intermediate(_, 1);
        NumericVector chuck3 = at_intermediate(_, 2);
        anna_whatever = as<RObject>(chuck1);
        annb_whatever = as<RObject>(chuck2);
        annc_whatever = as<RObject>(chuck3);
      } else if (at_rows == 3) {
        NumericVector chuck1 = at_intermediate(0, _);
        NumericVector chuck2 = at_intermediate(1, _);
        NumericVector chuck3 = at_intermediate(2, _);
        anna_whatever = as<RObject>(chuck1);
        annb_whatever = as<RObject>(chuck2);
        annc_whatever = as<RObject>(chuck3);
      }
    }
    
    if (is<NumericVector>(anna_whatever)) {
      NumericVector anna_another_int = as<NumericVector>(anna_whatever);
      
      int check_len = anna_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_anna_topull = anna_another_int;
      
    }
    
    if (is<NumericVector>(annb_whatever)) {
      NumericVector annb_another_int = as<NumericVector>(annb_whatever);
      
      int check_len = annb_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_annb_topull = annb_another_int;
      
    }
    
    if (is<NumericVector>(annc_whatever)) {
      NumericVector annc_another_int = as<NumericVector>(annc_whatever);
      
      int check_len = annc_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_annc_topull = annc_another_int;
      
    }
    
  } else {
    f_anna_topull = clone(model0);
    f_annb_topull = clone(model0);
    f_annc_topull = clone(model0);
  }
  
  if (!quiet) {
    if (greaterthan_warning) {
      Rf_warningcall(R_NilValue,
        "Too many values of annual covariates have been supplied. Some will be cut.");
    }
    if (lessthan_warning) {
      Rf_warningcall(R_NilValue,
        "Too few values of annual covariates have been supplied. Some will be cycled.");
    }
  }
  
  // dev_terms data frame or matrix
  NumericVector sur_dev_values(times);
  NumericVector obs_dev_values(times);
  NumericVector siz_dev_values(times);
  NumericVector sib_dev_values(times);
  NumericVector sic_dev_values(times);
  NumericVector rep_dev_values(times);
  NumericVector fec_dev_values(times);
  NumericVector jsur_dev_values(times);
  NumericVector jobs_dev_values(times);
  NumericVector jsiz_dev_values(times);
  NumericVector jsib_dev_values(times);
  NumericVector jsic_dev_values(times);
  NumericVector jrep_dev_values(times);
  NumericVector jmat_dev_values(times);
  
  NumericVector surv_dev_extracted;
  NumericVector obs_dev_extracted;
  NumericVector size_dev_extracted;
  NumericVector sizeb_dev_extracted;
  NumericVector sizec_dev_extracted;
  NumericVector repst_dev_extracted;
  NumericVector fec_dev_extracted;
  NumericVector jsurv_dev_extracted;
  NumericVector jobs_dev_extracted;
  NumericVector jsize_dev_extracted;
  NumericVector jsizeb_dev_extracted;
  NumericVector jsizec_dev_extracted;
  NumericVector jrepst_dev_extracted;
  NumericVector jmatst_dev_extracted;
  int veclimits {0};
  
  bool lessthan_warning_dev = false;
  bool greaterthan_warning_dev = false;
  
  if (dev_terms.isNotNull()) {
    RObject dt_intermediate = RObject(dev_terms);
    
    if (is<NumericMatrix>(dt_intermediate)) {
      NumericMatrix dt_mat = as<NumericMatrix>(dt_intermediate);
      int dt_rows = dt_mat.nrow();
      int dt_cols = dt_mat.ncol();
      
      if (dt_rows != 14 && dt_cols != 14) pop_error("dev_terms",
        "a data frame or matrix with 14 columns", "", 1);
      
      if (dt_rows == 14) {
        surv_dev_extracted = dt_mat(0, _);
        obs_dev_extracted = dt_mat(1, _);
        size_dev_extracted = dt_mat(2, _);
        sizeb_dev_extracted = dt_mat(3, _);
        sizec_dev_extracted = dt_mat(4, _);
        repst_dev_extracted = dt_mat(5, _);
        fec_dev_extracted = dt_mat(6, _);
        jsurv_dev_extracted = dt_mat(7, _);
        jobs_dev_extracted = dt_mat(8, _);
        jsize_dev_extracted = dt_mat(9, _);
        jsizeb_dev_extracted = dt_mat(10, _);
        jsizec_dev_extracted = dt_mat(11, _);
        jrepst_dev_extracted = dt_mat(12, _);
        jmatst_dev_extracted = dt_mat(13, _);
        
        veclimits = dt_cols;
        
        if (dt_cols < times) lessthan_warning_dev = true;
        if (dt_cols > times) greaterthan_warning_dev = true;
        
      } else {
        surv_dev_extracted = dt_mat(_, 0);
        obs_dev_extracted = dt_mat(_, 1);
        size_dev_extracted = dt_mat(_, 2);
        sizeb_dev_extracted = dt_mat(_, 3);
        sizec_dev_extracted = dt_mat(_, 4);
        repst_dev_extracted = dt_mat(_, 5);
        fec_dev_extracted = dt_mat(_, 6);
        jsurv_dev_extracted = dt_mat(_, 7);
        jobs_dev_extracted = dt_mat(_, 8);
        jsize_dev_extracted = dt_mat(_, 9);
        jsizeb_dev_extracted = dt_mat(_, 10);
        jsizec_dev_extracted = dt_mat(_, 11);
        jrepst_dev_extracted = dt_mat(_, 12);
        jmatst_dev_extracted = dt_mat(_, 13);
        
        veclimits = dt_rows;
        
        if (dt_rows < times) lessthan_warning_dev = true;
        if (dt_rows > times) greaterthan_warning_dev = true;
      }
      
    } else if (is<DataFrame>(dt_intermediate)) {
      DataFrame dt_frame = as<DataFrame>(dt_intermediate);
      int dt_vars = static_cast<int>(dt_frame.size());
      
      if (dt_vars != 14) pop_error("dev_terms",
        "a data frame or matrix with 14 columns.", "", 1);
      
      RObject surv_dev_a = as<RObject>(dt_frame[0]);
      RObject obs_dev_a = as<RObject>(dt_frame[1]);
      RObject size_dev_a = as<RObject>(dt_frame[2]);
      RObject sizeb_dev_a = as<RObject>(dt_frame[3]);
      RObject sizec_dev_a = as<RObject>(dt_frame[4]);
      RObject repst_dev_a = as<RObject>(dt_frame[5]);
      RObject fec_dev_a = as<RObject>(dt_frame[6]);
      RObject jsurv_dev_a = as<RObject>(dt_frame[7]);
      RObject jobs_dev_a = as<RObject>(dt_frame[8]);
      RObject jsize_dev_a = as<RObject>(dt_frame[9]);
      RObject jsizeb_dev_a = as<RObject>(dt_frame[10]);
      RObject jsizec_dev_a = as<RObject>(dt_frame[11]);
      RObject jrepst_dev_a = as<RObject>(dt_frame[12]);
      RObject jmatst_dev_a = as<RObject>(dt_frame[13]);
      
      if (is<NumericVector>(surv_dev_a)) {
        surv_dev_extracted = as<NumericVector>(surv_dev_a);
        veclimits = surv_dev_extracted.length();
        
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(obs_dev_a)) {
        obs_dev_extracted = as<NumericVector>(obs_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(size_dev_a)) {
        size_dev_extracted = as<NumericVector>(size_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(sizeb_dev_a)) {
        sizeb_dev_extracted = as<NumericVector>(sizeb_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(sizec_dev_a)) {
        sizec_dev_extracted = as<NumericVector>(sizec_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(repst_dev_a)) {
        repst_dev_extracted = as<NumericVector>(repst_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(fec_dev_a)) {
        fec_dev_extracted = as<NumericVector>(fec_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jsurv_dev_a)) {
        jsurv_dev_extracted = as<NumericVector>(jsurv_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jobs_dev_a)) {
        jobs_dev_extracted = as<NumericVector>(jobs_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jsize_dev_a)) {
        jsize_dev_extracted = as<NumericVector>(jsize_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jsizeb_dev_a)) {
        jsizeb_dev_extracted = as<NumericVector>(jsizeb_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jsizec_dev_a)) {
        jsizec_dev_extracted = as<NumericVector>(jsizec_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jrepst_dev_a)) {
        jrepst_dev_extracted = as<NumericVector>(jrepst_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      if (is<NumericVector>(jmatst_dev_a)) {
        jmatst_dev_extracted = as<NumericVector>(jmatst_dev_a);
      } else throw Rcpp::exception("Some dev_terms columns are not numeric.", false);
      
    } else pop_error("dev_terms", "a data frame with 14 numeric variables", "", 1);
  } else {
    veclimits = 1;
    surv_dev_extracted = clone(model0);
    obs_dev_extracted = clone(model0);
    size_dev_extracted = clone(model0);
    sizeb_dev_extracted = clone(model0);
    sizec_dev_extracted = clone(model0);
    repst_dev_extracted = clone(model0);
    fec_dev_extracted = clone(model0);
    jsurv_dev_extracted = clone(model0);
    jobs_dev_extracted = clone(model0);
    jsize_dev_extracted = clone(model0);
    jsizeb_dev_extracted = clone(model0);
    jsizec_dev_extracted = clone(model0);
    jrepst_dev_extracted = clone(model0);
    jmatst_dev_extracted = clone(model0);
  }
  
  if (!quiet) {
    if (greaterthan_warning_dev) {
      Rf_warningcall(R_NilValue,
        "Too many intercept deviations have been supplied. Some will be cut.");
    }
    if (lessthan_warning_dev) {
      Rf_warningcall(R_NilValue,
        "Too few intercept deviations have been supplied. Some will be cycled.");
    }
  }
  
  // Main for loop adjusting lengths of input vectors
  int year_counter {0};
  int patch_counter {0};
  int spdensity_counter {0};
  int finda_counter {0};
  int findb_counter {0};
  int findc_counter {0};
  int rinda_counter {0};
  int rindb_counter {0};
  int rindc_counter {0};
  
  int fanna_counter {0};
  int fannb_counter {0};
  int fannc_counter {0};
  
  int dev_counter {0};
  
  int year_limit = static_cast<int>(years_topull.length()) / nreps;
  int patch_limit = static_cast<int>(patches_topull.length());
  int spdensity_limit = static_cast<int>(spdensity_topull.length());
  int finda_limit = static_cast<int>(f_inda_topull.length());
  int findb_limit = static_cast<int>(f_indb_topull.length());
  int findc_limit = static_cast<int>(f_indc_topull.length());
  int rinda_limit = static_cast<int>(r_inda_topull.length());
  int rindb_limit = static_cast<int>(r_indb_topull.length());
  int rindc_limit = static_cast<int>(r_indc_topull.length());
  
  int fanna_limit = static_cast<int>(f_anna_topull.length());
  int fannb_limit = static_cast<int>(f_annb_topull.length());
  int fannc_limit = static_cast<int>(f_annc_topull.length());
  
  for (int i = 0; i < num_years; i++) {
    if (finda_counter >= finda_limit) finda_counter = 0;
    if (rinda_counter >= rinda_limit) rinda_counter = 0;
    if (findb_counter >= findb_limit) findb_counter = 0;
    if (rindb_counter >= rindb_limit) rindb_counter = 0;
    if (findc_counter >= findc_limit) findc_counter = 0;
    if (rindc_counter >= rindc_limit) rindc_counter = 0;
    
    if (fanna_counter >= fanna_limit) fanna_counter = 0;
    if (fannb_counter >= fannb_limit) fannb_counter = 0;
    if (fannc_counter >= fannc_limit) fannc_counter = 0;
    
    f2_inda_values(i) = f_inda_topull(finda_counter);
    r2_inda_values(i) = r_inda_topull(rinda_counter);
    f2_indb_values(i) = f_indb_topull(findb_counter);
    r2_indb_values(i) = r_indb_topull(rindb_counter);
    f2_indc_values(i) = f_indc_topull(findc_counter);
    r2_indc_values(i) = r_indc_topull(rindc_counter);
    
    f2_anna_values(i) = f_anna_topull(fanna_counter);
    f2_annb_values(i) = f_annb_topull(fanna_counter);
    f2_annc_values(i) = f_annc_topull(fanna_counter);
    
    if (i > 0) {
      f1_inda_values(i) = f2_inda_values(i-1);
      r1_inda_values(i) = r2_inda_values(i-1);
      f1_indb_values(i) = f2_indb_values(i-1);
      r1_indb_values(i) = r2_indb_values(i-1);
      f1_indc_values(i) = f2_indc_values(i-1);
      r1_indc_values(i) = r2_indc_values(i-1);
      
      f1_anna_values(i) = f2_anna_values(i-1);
      f1_annb_values(i) = f2_annb_values(i-1);
      f1_annc_values(i) = f2_annc_values(i-1);
      
    } else {
      r1_inda_values(i) = modelnone(0);
      r1_indb_values(i) = modelnone(0);
      r1_indc_values(i) = modelnone(0);
    }
    
    finda_counter++;
    rinda_counter++;
    findb_counter++;
    rindb_counter++;
    findc_counter++;
    rindc_counter++;
    
    fanna_counter++;
    fannb_counter++;
    fannc_counter++;
  }
  
  for (int i = 0; i < times; i++) {
    if (year_counter == year_limit) year_counter = 0;
    if (patch_counter == patch_limit) patch_counter = 0;
    if (spdensity_counter == spdensity_limit) spdensity_counter = 0;
    if (dev_counter >= veclimits) dev_counter = 0;
    
    for (int j = 0; j < nreps; j++) {
      years_projected(i, j) = years_topull(year_counter + times * j);
    }
    
    patches_projected(i) = patches_topull(patch_counter);
    spdensity_projected(i) = spdensity_topull(spdensity_counter);
    if (NumericVector::is_na(spdensity_projected(i))) spdensity_projected(i) = 0.0;
    
    sur_dev_values(i) = surv_dev_extracted(dev_counter);
    obs_dev_values(i) = obs_dev_extracted(dev_counter);
    siz_dev_values(i) = size_dev_extracted(dev_counter);
    sib_dev_values(i) = sizeb_dev_extracted(dev_counter);
    sic_dev_values(i) = sizec_dev_extracted(dev_counter);
    rep_dev_values(i) = repst_dev_extracted(dev_counter);
    fec_dev_values(i) = fec_dev_extracted(dev_counter);
    jsur_dev_values(i) = jsurv_dev_extracted(dev_counter);
    jobs_dev_values(i) = jobs_dev_extracted(dev_counter);
    jsiz_dev_values(i) = jsize_dev_extracted(dev_counter);
    jsib_dev_values(i) = jsizeb_dev_extracted(dev_counter);
    jsic_dev_values(i) = jsizec_dev_extracted(dev_counter);
    jrep_dev_values(i) = jrepst_dev_extracted(dev_counter);
    jmat_dev_values(i) = jmatst_dev_extracted(dev_counter);
    
    year_counter++;
    patch_counter++;
    spdensity_counter++;
    dev_counter++;
  }
  
  // Allstages
  DataFrame new_stageframe;
  arma::mat new_repmatrix;
  DataFrame new_ovtable;
  DataFrame allstages;
  
  if (format < 5) {
    bool agemat = false;
    bool historical = false;
    int ehrlen {1};
    int style {0};
    int filter {1};
    
    if (format == 2) ehrlen = 2;
    if (format == 3) style = 1;
    if (format == 4) {
      agemat = true;
      style = 2;
      filter = 2;
    }
    if (format < 3) historical = true;
    
    List melchett = LefkoMats::sf_reassess_internal(sframe, supplement,
      overwrite, repmatrix, agemat, historical, ehrlen);
    new_stageframe = as<DataFrame>(melchett["stageframe"]);
    new_repmatrix = as<arma::mat>(melchett["repmatrix"]);
    if (format < 4) {
      DataFrame new_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
      if (new_ovtable_temp.containsElementNamed("stage3")) {
        new_ovtable = new_ovtable_temp;
      } else {
        StringVector nsst3 = {};
        IntegerVector nsa2 = {};
        NumericVector nsgr = {};
        
        DataFrame intro_ovtable = DataFrame::create(_["stage3"] = nsst3,
          _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
          _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
          _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
          _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
          _["offset"] = clone(nsgr), _["multiplier"] = clone(nsgr),
          _["convtype"] = clone(nsa2), _["convtype_t12"] = clone(nsa2),
          _["pop"] = clone(nsst3), _["patch"] = clone(nsst3),
          _["year2"] = clone(nsst3));
        new_ovtable = intro_ovtable;
      }
    } else {
      DataFrame new_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
      if (new_ovtable_temp.containsElementNamed("stage3")) {
        new_ovtable = LefkoMats::age_expanded(new_ovtable_temp, start_age, last_age);
      } else {
        StringVector nsst3 = {};
        IntegerVector nsa2 = {};
        NumericVector nsgr = {};
        
        DataFrame intro_ovtable = DataFrame::create(_["stage3"] = nsst3,
          _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
          _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
          _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
          _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
          _["offset"] = clone(nsgr), _["multiplier"] = clone(nsgr),
          _["convtype"] = clone(nsa2), _["convtype_t12"] = clone(nsa2),
          _["pop"] = clone(nsst3), _["patch"] = clone(nsst3),
          _["year2"] = clone(nsst3));
        new_ovtable = intro_ovtable;
      }
    }
    
    // the old pizzle needs to be called
    DataFrame allstages_pre = theoldpizzle(new_stageframe, new_ovtable,
      new_repmatrix, start_age, last_age, ehrlen, style, cont, filter);
    allstages = allstages_pre;
    
  } else {
    DataFrame melchett = LefkoMats::sf_leslie(start_age, last_age, fecage_min,
      fecage_max, cont);
    new_stageframe = melchett;
    allstages = melchett;
    
    if (supplement.isNotNull()) {
      new_ovtable = LefkoMats::age_expanded(supplement, start_age, last_age);
    }
    
    CharacterVector maingroups_ch = {"0"};
    maingroups = as<RObject>(maingroups_ch);
    actualages = seq(start_age, last_age);
  }
  
  double maxsize {0.0};
  double maxsizeb {0.0};
  double maxsizec {0.0};
  
  if (format < 5) {
    NumericVector size3 = as<NumericVector>(allstages["size3"]);
    NumericVector size2n = as<NumericVector>(allstages["size2n"]);
    NumericVector size2o = as<NumericVector>(allstages["size2o"]);
    NumericVector sizeb3 = as<NumericVector>(allstages["sizeb3"]);
    NumericVector sizeb2n = as<NumericVector>(allstages["sizeb2n"]);
    NumericVector sizeb2o = as<NumericVector>(allstages["sizeb2o"]);
    NumericVector sizec3 = as<NumericVector>(allstages["sizec3"]);
    NumericVector sizec2n = as<NumericVector>(allstages["sizec2n"]);
    NumericVector sizec2o = as<NumericVector>(allstages["sizec2o"]);
    
    NumericVector maxveca = {max(size3), max(size2n), max(size2o)};
    NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
    NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
    
    maxsize = max(maxveca);
    maxsizeb = max(maxvecb);
    maxsizec = max(maxvecc);
  }
  
  // modelextract proxy lists
  CharacterVector my_char = mainyears;
  List surv_proxy = LefkoUtils::modelextract(surmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List obs_proxy = LefkoUtils::modelextract(obsmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List size_proxy = LefkoUtils::modelextract(sizmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List sizeb_proxy = LefkoUtils::modelextract(sibmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List sizec_proxy = LefkoUtils::modelextract(sicmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List repst_proxy = LefkoUtils::modelextract(repmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List fec_proxy = LefkoUtils::modelextract(fecmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jsurv_proxy = LefkoUtils::modelextract(jsurmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jobs_proxy = LefkoUtils::modelextract(jobsmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jsize_proxy = LefkoUtils::modelextract(jsizmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jsizeb_proxy = LefkoUtils::modelextract(jsibmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jsizec_proxy = LefkoUtils::modelextract(jsicmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jrepst_proxy = LefkoUtils::modelextract(jrepmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  List jmatst_proxy = LefkoUtils::modelextract(jmatmodl, pmnames, my_char,
    mainpatches, maingroups, inda_names, indb_names, indc_names, nodata);
  
  // Main projection set-up
  int yearnumber {0};
  int patchnumber {0};
  NumericVector used_devs;
  
  List madsexmadrigal_oneyear;
  List madsexmadrigal_forward;
  arma::mat Amat;
  arma::mat Umat;
  arma::mat Fmat;
  
  used_devs = {sur_dev_values(0), obs_dev_values(0), siz_dev_values(0), sib_dev_values(0),
    sic_dev_values(0), rep_dev_values(0), fec_dev_values(0), jsur_dev_values(0),
    jobs_dev_values(0), jsiz_dev_values(0), jsib_dev_values(0), jsic_dev_values(0),
    jrep_dev_values(0), jmat_dev_values(0)};
  
  // Vital rate density dependence inputs
  Rcpp::DataFrame dvr_frame;
  LogicalVector dvr_yn = {false, false, false, false, false, false, false, false,
    false, false, false, false, false, false};
  IntegerVector dvr_style = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  IntegerVector dvr_delay  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector dvr_alpha  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  NumericVector dvr_beta  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0};
  bool dens_vr = false;
  
  if (density_vr.isNotNull()) {
    if (!is<DataFrame>(density_vr)) pop_error("density_vr",
      "a data frame of class lefkoDensVR", "", 1);
    Rcpp::DataFrame dens_vr_thru(density_vr);
    dvr_frame = dens_vr_thru;
    
    int dvr_vars = static_cast<int>(dens_vr_thru.size());
    int dvr_rows = dens_vr_thru.nrows();
    
    if (dvr_vars != 6 || dvr_rows != 14) {
      String eat_my_shorts = "Data frame input for option density_vr does not match the ";
      eat_my_shorts += "dimensions of output from function density_vr().";
      
      throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
    }
    
    LogicalVector dvr_yn_ = as<LogicalVector>(dens_vr_thru["density_yn"]);
    IntegerVector dvr_style_ = as<IntegerVector>(dens_vr_thru["style"]);
    IntegerVector dvr_delay_ = as<IntegerVector>(dens_vr_thru["time_delay"]);
    NumericVector dvr_alpha_ = as<NumericVector>(dens_vr_thru["alpha"]);
    NumericVector dvr_beta_ = as<NumericVector>(dens_vr_thru["beta"]);
    
    for (int i = 0; i < dvr_style_.length(); i++) {
      if (dvr_yn_(i) == true) {
        if (dvr_style_(i) < 1 || dvr_style_(i) > 4) pop_error("density inputs", "", "", 21);
        
        if (NumericVector::is_na(dvr_alpha_(i))) pop_error("beta", "argument density", "", 8);
        if (NumericVector::is_na(dvr_beta_(i))) pop_error("beta", "argument density", "", 8);
        
        if (dvr_style_(i) == 1) {
          if (dvr_beta_(i) > exp_tol) {
            Rf_warningcall(R_NilValue, "Ricker beta outside limits. Resetting to exp_tol.");
            
            dvr_beta_(i) = exp_tol;
          } else if (dvr_beta_(i) < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue, "Ricker beta outside limits. Resetting to negative exp_tol.");
            
            dvr_beta_(i) = -1 * exp_tol;
          }
        } else if (dvr_style_(i) == 3) {
          double summed_stuff = dvr_alpha_(i) + dvr_beta_(i);
          
          if (summed_stuff > exp_tol) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high.");
          } else if (summed_stuff < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too low.");
          }
        }
      }
    }
    
    dvr_yn = dvr_yn_;
    dvr_style = dvr_style_;
    dvr_delay = dvr_delay_;
    dvr_alpha = dvr_alpha_;
    dvr_beta = dvr_beta_;
    
    dens_vr = true;
  }
  
  // Initial matrices to develop certain variables - not used in projection
  if (format < 5) {
    NumericVector st_dvr_dens = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0};
    
    madsexmadrigal_oneyear = jerzeibalowski(allstages, new_stageframe,
      format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
      repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
      jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
      f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
      f2_anna_values, f1_anna_values, f2_annb_values, f1_annb_values,
      f2_annc_values, f1_annc_values, r2_inda_values, r1_inda_values,
      r2_indb_values, r1_indb_values, r2_indc_values, r1_indc_values,
      r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
      r2_indc_values, r1_indc_values, used_devs, dens_vr, dvr_yn, dvr_style,
      dvr_alpha, dvr_beta, st_dvr_dens, spdensity_projected(0),
      repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
      yearnumber, patchnumber, exp_tol, theta_tol, ipm_cdf, err_check,
      true, sparse_switch, false);
    
  } else {
    NumericVector st_dvr_dens = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0};
    
    madsexmadrigal_oneyear = motherbalowski(actualages, new_stageframe,
      surv_proxy, fec_proxy, f2_inda_values, f1_inda_values, f2_indb_values,
      f1_indb_values, f2_indc_values, f1_indc_values, f2_anna_values,
      f1_anna_values, f2_annb_values, f1_annb_values, f2_annc_values,
      f1_annc_values, r2_inda_values, r1_inda_values, r2_indb_values,
      r1_indb_values, r2_indc_values, r1_indc_values, r2_inda_values,
      r1_inda_values, r2_indb_values, r1_indb_values, r2_indc_values,
      r1_indc_values, sur_dev_values(0), fec_dev_values(0),
      spdensity_projected(0), repmod, last_age, false, yearnumber, patchnumber,
      dens_vr, dvr_yn, dvr_style, dvr_alpha, dvr_beta, st_dvr_dens, exp_tol,
      theta_tol, true, sparse_switch, new_ovtable);
  }
  
  arma::sp_mat Umat_sp;
  int meanmatrows {0};
  
  if (sparse_switch == 1) {
    Umat_sp = as<arma::sp_mat>(madsexmadrigal_oneyear["U"]);
    meanmatrows = Umat_sp.n_rows;
  } else {
    Umat = as<arma::mat>(madsexmadrigal_oneyear["U"]);
    meanmatrows = Umat.n_rows;
  }
  
  DataFrame ahstages;
  if (format < 5) {
    ahstages = sframe;
  } else {
    ahstages = new_stageframe;
  }
  
  DataFrame hstages;
  if (format < 3) {
    hst_maker(hstages, ahstages, format);
  } else {
    hstages = R_NilValue;
  }
  
  DataFrame agestages;
  if (format == 4) {
    agestages = age_maker(ahstages, start_age, last_age);
  } else {
    agestages = R_NilValue;
  }
  
  // Check if matrix is large and sparse
  if (sparse_auto) {
    int test_elems {0};
    arma::uvec nonzero_elems;
    
    if (sparse_switch == 1) {
      test_elems = static_cast<int>(Umat_sp.n_elem);
      nonzero_elems = find(Umat_sp);
    } else {
      test_elems = static_cast<int>(Umat.n_elem);
      nonzero_elems = find(Umat);
    }
    
    int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
    double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
    if (sparse_check <= 0.5 && meanmatrows > 50) {
      sparse_switch = 1;
    } else { 
      sparse_switch = 0;
    }
  }
  
  // Start data frame and vector
  arma::vec startvec;
  if(start_frame.isNotNull()) {
    if (!is<DataFrame>(start_frame)) pop_error("start_frame",
      "data frame of class lefkoSV", "", 1);
    Rcpp::DataFrame start_thru(start_frame);
    startvec.set_size(meanmatrows);
    startvec.zeros();
    arma::uvec start_elems = as<arma::uvec>(start_thru["row_num"]);
    start_elems = start_elems - 1;
    arma::vec start_values = as<arma::vec>(start_thru["value"]);
    
    if (static_cast<int>(start_elems.max()) > (meanmatrows - 1)) {
      throw Rcpp::exception("Start vector input frame includes element indices too high for this MPM.", false);
    }
    for (int i = 0; i < static_cast<int>(start_elems.n_elem); i++) {
      startvec(start_elems(i)) = start_values(i);
    }
    
  } else if (start_vec.isNotNull()) {
    startvec = as<arma::vec>(start_vec);
    if (static_cast<int>(startvec.n_elem) != meanmatrows) {
      pop_error("start_vec", "rows", "MPM matrices", 19);
    }
    
  } else {
    startvec.set_size(meanmatrows);
    startvec.ones();
  }
  
  // Matrix element density dependence inputs
  Rcpp::DataFrame dens_input;
  List dens_index;
  double pop_size {0};
  double changing_element_U {0.0};
  double changing_element_F {0.0};
  
  int time_delay {1};
  bool dens_elems = false;
  bool warn_trigger_neg = false;
  bool warn_trigger_1 = false;
  
  arma::uvec dyn_index321;
  arma::uvec dyn_index_col;
  arma::uvec dyn_style;
  arma::vec dyn_alpha;
  arma::vec dyn_beta;
  arma::uvec dyn_delay;
  arma::uvec dyn_type;
  int n_dyn_elems {0};
  
  if (density.isNotNull()) {
    if (!is<DataFrame>(density)) pop_error("density",
      "data frame of class lefkoDens", "", 1);
    Rcpp::DataFrame dens_thru(density);
    dens_input = dens_thru;
    
    Rcpp::StringVector di_stage3 = as<StringVector>(dens_input["stage3"]);
    Rcpp::StringVector di_stage2 = as<StringVector>(dens_input["stage2"]);
    Rcpp::StringVector di_stage1 = as<StringVector>(dens_input["stage1"]);
    int di_size = di_stage3.length();
    
    if (format < 3) { // Historical matrices
      StringVector stage3 = as<StringVector>(hstages["stage_2"]);
      StringVector stage2r = as<StringVector>(hstages["stage_1"]);
      StringVector stage2c = as<StringVector>(hstages["stage_2"]);
      StringVector stage1 = as<StringVector>(hstages["stage_1"]);
      int hst_size = stage3.length();
      
      arma::uvec hst_3(hst_size, fill::zeros);
      arma::uvec hst_2r(hst_size, fill::zeros);
      arma::uvec hst_2c(hst_size, fill::zeros);
      arma::uvec hst_1(hst_size, fill::zeros);
      
      arma::uvec di_stage32_id(di_size, fill::zeros);
      arma::uvec di_stage21_id(di_size, fill::zeros);
      arma::uvec di_index(di_size, fill::zeros);
      
      for (int i = 0; i < di_size; i++) {
        for (int j = 0; j < hst_size; j++) {
          if (di_stage3(i) == stage3(j)) {
            hst_3(j) = 1;
          } else {
            hst_3(j) = 0;
          }
        }
        
        for (int j = 0; j < hst_size; j++) {
          if (di_stage2(i) == stage2r(j)) {
            hst_2r(j) = 1;
          } else {
            hst_2r(j) = 0;
          }
        }
        
        for (int j = 0; j < hst_size; j++) {
          if (di_stage2(i) == stage2c(j)) {
            hst_2c(j) = 1;
          } else {
            hst_2c(j) = 0;
          }
        }
        
        for (int j = 0; j < hst_size; j++) {
          if (di_stage1(i) == stage1(j)) {
            hst_1(j) = 1;
          } else {
            hst_1(j) = 0;
          }
        }
        
        arma::uvec find_hst3 = find(hst_3);
        arma::uvec find_hst2r = find(hst_2r);
        arma::uvec find_hst2c = find(hst_2c);
        arma::uvec find_hst1 = find(hst_1);
        
        arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
        arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
        
        di_stage32_id(i) = pop_32(0);
        di_stage21_id(i) = pop_21(0);
        di_index(i) = pop_32(0) + (pop_21(0) * hst_size);
        
        hst_3.zeros();
        hst_2r.zeros();
        hst_2c.zeros();
        hst_1.zeros();
      }
      
      dens_index = Rcpp::List::create(_["index32"] = di_stage32_id,
        _["index21"] = di_stage21_id, _["index321"] = di_index);
      
    } else { // Ahistorical and age-based matrices
      ///// Need an age-by-stage version handling multiple age inputs
      
      
      StringVector stage3 = as<StringVector>(ahstages["stage"]);
      StringVector stage2 = as<StringVector>(ahstages["stage"]);
      int ahst_size = stage3.length();
      
      arma::uvec ahst_3(ahst_size, fill::zeros);
      arma::uvec ahst_2(ahst_size, fill::zeros);
      
      arma::uvec di_stage32_id(di_size, fill::zeros);
      arma::uvec di_stage21_id(di_size, fill::zeros);
      arma::uvec di_index(di_size, fill::zeros);
      
      for (int i = 0; i < di_size; i++) {
        for (int j = 0; j < ahst_size; j++) {
          if (di_stage3(i) == stage3(j)) {
            ahst_3(j) = 1;
          } else {
            ahst_3(j) = 0;
          }
        }
        
        for (int j = 0; j < ahst_size; j++) {
          if (di_stage2(i) == stage2(j)) {
            ahst_2(j) = 1;
          } else {
            ahst_2(j) = 0;
          }
        }
        
        arma::uvec find_ahst3 = find(ahst_3);
        arma::uvec find_ahst2 = find(ahst_2);
        di_stage32_id(i) = find_ahst3(0);
        di_stage21_id(i) = find_ahst2(0);
        di_index(i) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
        
        ahst_3.zeros();
        ahst_2.zeros();
      }
      
      dens_index = Rcpp::List::create(_["index3"] = di_stage32_id,
        _["index2"] = di_stage21_id, _["index321"] = di_index);
    }
    dyn_index321 = as<arma::uvec>(dens_index["index321"]);
    dyn_index_col = as<arma::uvec>(dens_index[1]);
    dyn_style = as<arma::uvec>(dens_input["style"]);
    dyn_alpha = as<arma::vec>(dens_input["alpha"]);
    dyn_beta = as<arma::vec>(dens_input["beta"]);
    dyn_delay = as<arma::uvec>(dens_input["time_delay"]);
    dyn_type = as<arma::uvec>(dens_input["type"]);
    n_dyn_elems = static_cast<int>(dyn_index321.n_elem);
    
    for (int i = 0; i < static_cast<int>(dyn_style.n_elem); i++) {
      if (dyn_style(i) < 1 || dyn_style(i) > 4) pop_error("density inputs", "", "", 21);
      
      if (dyn_style(i) == 1) {
        if (dyn_beta(i) > exp_tol) {
          Rf_warningcall(R_NilValue, "Ricker beta outside limits. Resetting to exp_tol.");
          
          dyn_beta(i) = exp_tol;
        } else if (dyn_beta(i) < (-1.0 * exp_tol)) {
          Rf_warningcall(R_NilValue, "Ricker beta outside limits. Resetting to negative exp_tol.");
          
          dyn_beta(i) = -1 * exp_tol;
        }
      } else if (dyn_style(i) == 3) {
        double summed_stuff = dyn_alpha(i) + dyn_beta(i);
        
        if (summed_stuff > exp_tol) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high.");
        } else if (summed_stuff < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high.");
        }
      }
    }
    dens_elems = true;
  }
  
  if (dens_vr && dens_elems) {
    Rf_warningcall(R_NilValue,
      "Density dependence should be set either to vital rate models or matrix elements, not both.");
  }
  
  // Main projection loop
  List A_all(nreps);
  List F_all(nreps);
  List U_all(nreps);
  List out_all(nreps);
  List A_mats(times);
  List F_mats(times);
  List U_mats(times);
  List out_mats(times);
  
  IntegerMatrix years_int = refsort_str_m(years_projected, mainyears);
  IntegerVector patches_int = refsort_str(patches_projected, mainpatches);
  
  arma::vec theseventhson = startvec;
  arma::rowvec theseventhgrandson = startvec.as_row();
  arma::mat popproj(meanmatrows, (times + 1), fill::zeros); // Population vector
  arma::mat wpopproj(meanmatrows, (times + 1), fill::zeros); // Population w vector
  arma::mat vpopproj(meanmatrows, (times + 1)); // Population v vector
  arma::mat Rvecmat(1, (times + 1), fill::zeros);
  arma::mat thesecondprophecy;
  
  List all_projections (nreps);
  List all_stagedist (nreps);
  List all_repvalues (nreps);
  arma::mat all_R (nreps, times+1);
  
  NumericVector usable_densities = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0};
  
  if (sparse_switch == 0 || format == 5) {
    for (int rep = 0; rep < nreps; rep++) {
      Rvecmat.zeros();
      popproj.zeros();
      wpopproj.zeros();
      vpopproj.zeros();
      
      popproj.col(0) = startvec;
      Rvecmat(0) = sum(startvec);
      if (!growthonly) {
        wpopproj.col(0) = startvec / sum(startvec);
      }
      
      theseventhson = startvec;
      arma::rowvec theseventhgrandson = startvec.as_row();
      
      for (int i = 0; i < times; i++) {
        if (i % 10 == 0) Rcpp::checkUserInterrupt();
        
        yearnumber = years_int(i, rep) - 1;
        patchnumber = patches_int(i) - 1;
        
        used_devs = {sur_dev_values(i), obs_dev_values(i), siz_dev_values(i), 
          sib_dev_values(i), sic_dev_values(i), rep_dev_values(i),
          fec_dev_values(i), jsur_dev_values(i), jobs_dev_values(i),
          jsiz_dev_values(i), jsib_dev_values(i), jsic_dev_values(i),
          jrep_dev_values(i), jmat_dev_values(i)};
        
        if (dens_vr) {
          for (int j = 0; j < 14; j++) {
            if (dvr_delay(j) <= i) {
              usable_densities(j) = Rvecmat(i - dvr_delay(j));
            }
          }
        }
        
        if (format < 5) {
          madsexmadrigal_oneyear = jerzeibalowski(allstages, new_stageframe,
            format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
            repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
            jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
            f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
            f2_anna_values, f1_anna_values, f2_annb_values, f1_annb_values,
            f2_annc_values, f1_annc_values, r2_inda_values, r1_inda_values,
            r2_indb_values, r1_indb_values, r2_indc_values, r1_indc_values,
            r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
            r2_indc_values, r1_indc_values, used_devs, dens_vr, dvr_yn,
            dvr_style, dvr_alpha, dvr_beta, usable_densities,
            spdensity_projected(i), repmod, maxsize, maxsizeb, maxsizec,
            start_age, last_age, false, yearnumber, patchnumber, exp_tol,
            theta_tol, ipm_cdf, err_check, true, sparse_switch, false);
          
        } else {
          madsexmadrigal_oneyear = motherbalowski(actualages, new_stageframe,
            surv_proxy, fec_proxy, f2_inda_values, f1_inda_values, f2_indb_values,
            f1_indb_values, f2_indc_values, f1_indc_values, f2_anna_values,
            f1_anna_values, f2_annb_values, f1_annb_values, f2_annc_values,
            f1_annc_values, r2_inda_values, r1_inda_values, r2_indb_values,
            r1_indb_values, r2_indc_values, r1_indc_values, r2_inda_values,
            r1_inda_values, r2_indb_values, r1_indb_values, r2_indc_values,
            r1_indc_values, sur_dev_values(i), fec_dev_values(i),
            spdensity_projected(i), repmod, last_age, false, yearnumber,
            patchnumber, dens_vr, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
            usable_densities, exp_tol, theta_tol, true, sparse_switch, new_ovtable);
        }
        
        Umat = as<arma::mat>(madsexmadrigal_oneyear["U"]);
        Fmat = as<arma::mat>(madsexmadrigal_oneyear["F"]);
        
        for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
          time_delay = dyn_delay(j);
          if (time_delay > 0) time_delay = time_delay - 1;
          
          if (i >= time_delay) {
            pop_size = sum(popproj.col(i - time_delay));
            
            if (dyn_style(j) == 1) { // Ricker
              changing_element_U = Umat(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              changing_element_F = Fmat(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              
            } else if (dyn_style(j) == 2) { // Beverton-Holt
              changing_element_U = Umat(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              changing_element_F = Fmat(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              
            } else if (dyn_style(j) == 3) { // Usher function
              changing_element_U = Umat(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              changing_element_F = Fmat(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              
            } else if (dyn_style(j) == 4) { // Logistic function
              double used_popsize = pop_size;
              if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
                used_popsize = dyn_alpha(j);
              }
              changing_element_U = Umat(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
              changing_element_F = Fmat(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
            }
            
            if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element_U > 1.0) {
                changing_element_U = 1.0;
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
            } else if (substoch == 2 && dyn_type(j) == 1) {
              arma::vec given_col = Umat.col(dyn_index_col(j));
              arma::uvec gc_negs = find(given_col < 0.0);
              arma::uvec gc_pos = find(given_col > 0.0);
              
              double barnyard_antics = sum(given_col(gc_pos)) - Umat(dyn_index321(j)) +
                changing_element_U;
              
              if (barnyard_antics > 1.0 && changing_element_U > 0.0) {
                double proposed_element_U = changing_element_U - barnyard_antics *
                  (changing_element_U / barnyard_antics);
                
                if (proposed_element_U >= 0.0) {
                  changing_element_U = proposed_element_U;
                } else {
                  changing_element_U = 0.0;
                }
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
            } else if (substoch > 0 && dyn_type(j) == 2) {
              if (changing_element_F < 0.0) {
                changing_element_F = 0.0;
              }
            }
            Umat(dyn_index321(j)) = changing_element_U;
            Fmat(dyn_index321(j)) = changing_element_F;
            
            if (dyn_type(j) == 1 && Umat(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              if (!quiet) Rf_warningcall(R_NilValue,
                "Some probabilities with value > 1.0 produced during density adjustment.");
            } else if ((Umat(dyn_index321(j)) < 0.0 || Fmat(dyn_index321(j)) < 0.0) && 
              !warn_trigger_neg) {
              warn_trigger_neg = true;
              if (!quiet) Rf_warningcall(R_NilValue,
                "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
        
        Amat = Umat + Fmat;
        
        if (err_check) {
          F_mats(i) = Fmat;
          U_mats(i) = Umat;
          A_mats(i) = Amat;
          if (format < 5) out_mats(i) = as<DataFrame>(madsexmadrigal_oneyear["out"]);
        }
        
        theseventhson = Amat * theseventhson;
        if (integeronly) {
          theseventhson = floor(theseventhson);
        }
        popproj.col(i+1) = theseventhson;
        Rvecmat(i+1) = sum(theseventhson);
        
        if (Rvecmat(i+1) <= 0.0) break;
        
        if (standardize) {
          theseventhson = theseventhson / sum(theseventhson);
        }
        if (!growthonly) {
          wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
          
          if (repvalue && !dens_vr) { // Currently does not handle den dep repvalue
            NumericVector second_devs = {sur_dev_values(times - (i+1)),
              obs_dev_values(times - (i+1)), siz_dev_values(times - (i+1)),
              sib_dev_values(times - (i+1)), sic_dev_values(times - (i+1)),
              rep_dev_values(times - (i+1)), fec_dev_values(times - (i+1)),
              jsur_dev_values(times - (i+1)), jobs_dev_values(times - (i+1)),
              jsiz_dev_values(times - (i+1)), jsib_dev_values(times - (i+1)),
              jsic_dev_values(times - (i+1)), jrep_dev_values(times - (i+1)),
              jmat_dev_values(times - (i+1))};
            
            if (format < 5) {
              madsexmadrigal_forward = jerzeibalowski(allstages, new_stageframe,
                format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
                repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
                jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
                f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
                f2_anna_values, f1_anna_values, f2_annb_values, f1_annb_values,
                f2_annc_values, f1_annc_values, r2_inda_values, r1_inda_values,
                r2_indb_values, r1_indb_values, r2_indc_values, r1_indc_values,
                r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
                r2_indc_values, r1_indc_values, used_devs, false, dvr_yn, dvr_style,
                dvr_alpha, dvr_beta, usable_densities, spdensity_projected(times - (i+1)),
                repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
                yearnumber, patchnumber, exp_tol, theta_tol, ipm_cdf, err_check,
                true, sparse_switch, false);
              
            } else {
              madsexmadrigal_forward = motherbalowski(actualages, new_stageframe,
                surv_proxy, fec_proxy, f2_inda_values, f1_inda_values, f2_indb_values,
                f1_indb_values, f2_indc_values, f1_indc_values, f2_anna_values,
                f1_anna_values, f2_annb_values, f1_annb_values, f2_annc_values,
                f1_annc_values, r2_inda_values, r1_inda_values, r2_indb_values,
                r1_indb_values, r2_indc_values, r1_indc_values, r2_inda_values,
                r1_inda_values, r2_indb_values, r1_indb_values, r2_indc_values,
                r1_indc_values, sur_dev_values(times - (i+1)), fec_dev_values(times - (i+1)),
                spdensity_projected(times - (i+1)), repmod, last_age, false, yearnumber,
                patchnumber, false, dvr_yn, dvr_style, dvr_alpha, dvr_beta,
                usable_densities, exp_tol, theta_tol, true, sparse_switch, new_ovtable);
            }
            arma::mat second_U = as<arma::mat>(madsexmadrigal_forward["U"]);
            arma::mat second_F = as<arma::mat>(madsexmadrigal_forward["F"]);
            thesecondprophecy = second_U + second_F;
            theseventhgrandson = theseventhgrandson * thesecondprophecy;
            
            double seventhgrandsum = sum(theseventhgrandson);
            arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
            theseventhgrandson = theseventhgrandson / seventhgrandsum;
            vpopproj.col(times - (i+1)) = midwife;
          } // if repvalue
        } // if growthonly
      } // times for loop
      
      all_projections(rep) = popproj;
      all_R.row(rep) = Rvecmat;
      
      if (!growthonly) {
        all_stagedist(rep) = wpopproj;
        if (repvalue) all_repvalues(rep) = vpopproj;
      }
      
      if (err_check) {
        A_all(rep) = A_mats;
        F_all(rep) = F_mats;
        U_all(rep) = U_mats;
        if (format < 5) out_all(rep) = out_mats;
      }
    } // nreps for loop
  } else {
    arma::sp_mat Fmat_sp;
    arma::sp_mat Amat_sp;
    arma::sp_mat thesecondprophecy_sp;
    arma::sp_mat theseventhson_sp(theseventhson);
    arma::sp_mat theseventhgrandson_sp(theseventhgrandson);
    
    for (int rep = 0; rep < nreps; rep++) {
      Rvecmat.zeros();
      popproj.zeros();
      wpopproj.zeros();
      vpopproj.zeros();
      
      popproj.col(0) = startvec;
      Rvecmat(0) = sum(startvec);
      if (!growthonly) {
        wpopproj.col(0) = startvec / sum(startvec);
      }
      
      for (int i = 0; i < times; i++) {
        if (i % 25 == 0) Rcpp::checkUserInterrupt();
        
        yearnumber = years_int(i, rep) - 1;
        patchnumber = patches_int(i) - 1;
        
        used_devs = {sur_dev_values(i), obs_dev_values(i), siz_dev_values(i),
          sib_dev_values(i), sic_dev_values(i), rep_dev_values(i),
          fec_dev_values(i), jsur_dev_values(i), jobs_dev_values(i),
          jsiz_dev_values(i), jsib_dev_values(i), jsic_dev_values(i),
          jrep_dev_values(i), jmat_dev_values(i)};
        
        if (dens_vr) {
          for (int j = 0; j < 14; j++) {
            if (dvr_delay(j) <= i) {
              usable_densities(j) = Rvecmat(i - dvr_delay(j));
            }
          }
        }
        
        madsexmadrigal_oneyear = jerzeibalowski(allstages, new_stageframe,
          format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
          repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
          jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
          f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
          f2_anna_values, f1_anna_values, f2_annb_values, f1_annb_values,
          f2_annc_values, f1_annc_values, r2_inda_values, r1_inda_values,
          r2_indb_values, r1_indb_values, r2_indc_values, r1_indc_values,
          r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
          r2_indc_values, r1_indc_values, used_devs, dens_vr, dvr_yn, dvr_style,
          dvr_alpha, dvr_beta, usable_densities, spdensity_projected(i), repmod,
          maxsize, maxsizeb, maxsizec, start_age, last_age, false, yearnumber,
          patchnumber, exp_tol, theta_tol, ipm_cdf, err_check, true,
          sparse_switch, false);
        
        Umat_sp = as<arma::sp_mat>(madsexmadrigal_oneyear["U"]);
        Fmat_sp = as<arma::sp_mat>(madsexmadrigal_oneyear["F"]);
        
       for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
          time_delay = dyn_delay(j);
          if (time_delay > 0) time_delay = time_delay - 1;
          
          if (i >= time_delay) {
            pop_size = sum(popproj.col(i - time_delay));
            
            if (dyn_style(j) == 1) { // Ricker
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              
            } else if (dyn_style(j) == 2) { // Beverton-Holt
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              
            } else if (dyn_style(j) == 3) { // Usher function
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              
            } else if (dyn_style(j) == 4) { // Logistic function
              double used_popsize = pop_size;
              if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
                used_popsize = dyn_alpha(j);
              }
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
            }
            
            if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element_U > 1.0) {
                changing_element_U = 1.0;
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
              
            } else if (substoch == 2 && dyn_type(j) == 1) {
              arma::vec given_col = arma::vec(Umat_sp.col(dyn_index_col(j)));
              arma::uvec gc_negs = find(given_col < 0.0);
              arma::uvec gc_pos = find(given_col > 0.0);
              
              double barnyard_antics = sum(given_col(gc_pos)) - Umat_sp(dyn_index321(j)) +
                changing_element_U;
              
              if (barnyard_antics > 1.0 && changing_element_U > 0.0) {
                double proposed_element_U = changing_element_U - barnyard_antics *
                  (changing_element_U / barnyard_antics);
                
                if (proposed_element_U >= 0.0) {
                  changing_element_U = proposed_element_U;
                } else {
                  changing_element_U = 0.0;
                }
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
            } else if (substoch > 0 && dyn_type(j) == 2) {
              if (changing_element_F < 0.0) {
                changing_element_F = 0.0;
              }
            }
            Umat_sp(dyn_index321(j)) = changing_element_U;
            Fmat_sp(dyn_index321(j)) = changing_element_F;
            
            if (dyn_type(j) == 1 && Umat_sp(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              if (!quiet) Rf_warningcall(R_NilValue,
                "Some probabilities with value > 1.0 produced during density adjustment.");
            } else if ((Umat_sp(dyn_index321(j)) < 0.0 || Fmat_sp(dyn_index321(j)) < 0.0) &&
              !warn_trigger_neg) {
              warn_trigger_neg = true;
              if (!quiet) Rf_warningcall(R_NilValue,
                "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
        
        Amat_sp = Umat_sp + Fmat_sp;
        
        if (err_check) {
          F_mats(i) = Fmat_sp;
          U_mats(i) = Umat_sp;
          A_mats(i) = Amat_sp;
          if (format < 5) out_mats(i) = as<DataFrame>(madsexmadrigal_oneyear["out"]);
        }
        
        theseventhson_sp = Amat_sp * theseventhson_sp;
        if (integeronly) {
          theseventhson_sp = floor(theseventhson_sp);
        }
        popproj.col(i+1) = arma::vec(theseventhson_sp);
        Rvecmat(i+1) = sum(popproj.col(i+1));
        
        if (Rvecmat(i+1) <= 0.0) break;
        
        if (standardize) {
          theseventhson_sp = theseventhson_sp / Rvecmat(i+1);
        }
        if (!growthonly) {
          wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
          
          if (repvalue && !dens_vr) { // Currently does not handle dens dep repvalue
            NumericVector second_devs = {sur_dev_values(times - (i+1)),
              obs_dev_values(times - (i+1)), siz_dev_values(times - (i+1)),
              sib_dev_values(times - (i+1)), sic_dev_values(times - (i+1)),
              rep_dev_values(times - (i+1)), fec_dev_values(times - (i+1)),
              jsur_dev_values(times - (i+1)), jobs_dev_values(times - (i+1)),
              jsiz_dev_values(times - (i+1)), jsib_dev_values(times - (i+1)),
              jsic_dev_values(times - (i+1)), jrep_dev_values(times - (i+1)),
              jmat_dev_values(times - (i+1))};
            
            madsexmadrigal_forward = jerzeibalowski(allstages, new_stageframe,
              format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
              repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
              jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
              f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
              f2_anna_values, f1_anna_values, f2_annb_values, f1_annb_values,
              f2_annc_values, f1_annc_values, r2_inda_values, r1_inda_values,
              r2_indb_values, r1_indb_values, r2_indc_values, r1_indc_values,
              r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
              r2_indc_values, r1_indc_values, used_devs, false, dvr_yn, dvr_style,
              dvr_alpha, dvr_beta, usable_densities, spdensity_projected(times - (i+1)),
              repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
              yearnumber, patchnumber, exp_tol, theta_tol, ipm_cdf, err_check,
              true, true, false);
              
            arma::sp_mat second_U = as<arma::sp_mat>(madsexmadrigal_forward["U"]);
            arma::sp_mat second_F = as<arma::sp_mat>(madsexmadrigal_forward["F"]);
            thesecondprophecy_sp = second_U + second_F;
            theseventhgrandson_sp = theseventhgrandson_sp * thesecondprophecy_sp;
            
            arma::rowvec partway(theseventhgrandson_sp);
            double seventhgrandsum = sum(partway);
            arma::sp_mat midwife = theseventhgrandson_sp / seventhgrandsum;
            theseventhgrandson_sp = midwife;
            vpopproj.col(times - (i+1)) = midwife;
          } // if repvalue
        } // if growthonly
      } // times for loop
      
      all_projections(rep) = popproj;
      all_R.row(rep) = Rvecmat;
      
      if (!growthonly) {
        all_stagedist(rep) = wpopproj;
        if (repvalue) all_repvalues(rep) = vpopproj;
      }
      
      if (err_check) {
        A_all(rep) = A_mats;
        F_all(rep) = F_mats;
        U_all(rep) = U_mats;
        if (format < 5) out_all(rep) = out_mats;
      }
    } // nreps for loop
  }
  
  // Final output prep
  IntegerVector control = {nreps, times};
  DataFrame newlabels = DataFrame::create(_["pop"] = 1, _["patch"] = chosenpatch);
  CharacterVector output_class = {"lefkoProj"};
  List output;
  
  List all_R_list(1);
  all_R_list(0) = all_R;
  
  List projections(1);
  projections(0) = all_projections;
  
  List stagedist(1);
  stagedist(0) = all_stagedist;
  
  List repvalues(1);
  repvalues(0) = all_repvalues;
  
  if (err_check) {
    List output_err(38);
    
    output_err(0) = projections;
    output_err(1) = stagedist;
    output_err(2) = repvalues;
    output_err(3) = all_R_list;
    output_err(4) = ahstages;
    output_err(5) = hstages;
    output_err(6) = agestages;
    output_err(7) = newlabels;
    output_err(8) = control;
    output_err(9) = dens_input;
    output_err(10) = dvr_frame;
    output_err(11) = pmnames;
    output_err(12) = mainyears;
    output_err(13) = mainpatches;
    output_err(14) = maingroups;
    output_err(15) = mainages;
    output_err(16) = allstages;
    output_err(17) = surv_proxy;
    output_err(18) = obs_proxy;
    output_err(19) = size_proxy;
    output_err(20) = sizeb_proxy;
    output_err(21) = sizec_proxy;
    output_err(22) = repst_proxy;
    output_err(23) = fec_proxy;
    output_err(24) = jsurv_proxy;
    output_err(25) = jobs_proxy;
    output_err(26) = jsize_proxy;
    output_err(27) = jsizeb_proxy;
    output_err(28) = jsizec_proxy;
    output_err(29) = jrepst_proxy;
    output_err(30) = jmatst_proxy;
    output_err(31) = years_projected;
    output_err(32) = patches_projected;
    output_err(33) = spdensity_projected;
    output_err(34) = A_all;
    output_err(35) = U_all;
    output_err(36) = F_all;
    output_err(37) = out_all;
    
    CharacterVector output_err_names = {"projection", "stage_dist", "rep_value",
      "pop_size", "ahstages", "hstages", "agestages", "labels", "control", "density",
      "density_vr", "paramnames", "mainyears", "mainpatches", "maingroups", "mainages",
      "allstages", "surv_proxy", "obs_proxy", "size_proxy", "sizeb_proxy",
      "sizec_proxy", "repst_proxy", "fec_proxy", "jsurv_proxy", "jobs_proxy",
      "jsize_proxy", "jsizeb_proxy", "jsizec_proxy", "jrepst_proxy", "jmatst_proxy",
      "years_projected", "patches_projected", "spdensity_projected", "A_all",
      "U_all", "F_all", "out_all"};
    output_err.attr("names") = output_err_names;
    output_err.attr("class") = output_class;
    
    output = output_err;
  } else {
    List output_noerr(11);
    
    output_noerr(0) = projections;
    output_noerr(1) = stagedist;
    output_noerr(2) = repvalues;
    output_noerr(3) = all_R_list;
    output_noerr(4) = ahstages;
    output_noerr(5) = hstages;
    output_noerr(6) = agestages;
    output_noerr(7) = newlabels;
    output_noerr(8) = control;
    output_noerr(9) = dens_input;
    output_noerr(10) = dvr_frame;
    
    CharacterVector output_noerr_names = {"projection", "stage_dist", "rep_value",
      "pop_size", "ahstages", "hstages", "agestages", "labels", "control", "density",
      "density_vr"};
    output_noerr.attr("names") = output_noerr_names;
    output_noerr.attr("class") = output_class;
    
    output = output_noerr;
  }
  
  return output;
}

//' General Matrix Projection Model Creation
//' 
//' Function \code{mpm_create()} is the core workhorse function that creates
//' all flavors of MPM in \code{lefko3}. All other MPM creation functions act
//' as wrappers for this function. As such, this function provides the most
//' general and most detailed control over the MPM creation process.
//' 
//' @name mpm_create
//' 
//' @param historical A logical value indicating whether to build a historical
//' MPM. Defaults to \code{FALSE}.
//' @param stage A logical value indicating whether to build a stage-based MPM.
//' If both \code{stage = TRUE} and \code{age = TRUE}, then will proceed to
//' build an age-by-stage MPM. Defaults to \code{TRUE}.
//' @param age A logical value indicating whether to build an age-based MPM. If
//' both \code{stage = TRUE} and \code{age = TRUE}, then will proceed to build
//' an age-by-stage MPM. Defaults to \code{FALSE}.
//' @param devries  A logical value indicating whether to use deVries format
//' for historical MPMs. Defaults to \code{FALSE}, in which case historical MPMs
//' are created in Ehrlen format.
//' @param reduce A logical value denoting whether to remove ages, ahistorical
//' stages, or historical stages associated exclusively with zero transitions.
//' These are removed only if the respective row and column sums in ALL matrices
//' estimated equal 0. Defaults to \code{FALSE}.
//' @param simple A logical value indicating whether to produce \code{A},
//' \code{U}, and \code{F} matrices, or only the latter two. Defaults to
//' \code{FALSE}, in which case all three are output.
//' @param err_check A logical value indicating whether to append extra
//' information used in matrix calculation within the output list. Defaults to
//' \code{FALSE}.
//' @param data A data frame of class \code{hfvdata}. Required for all MPMs,
//' except for function-based MPMs in which \code{modelsuite} is set to a
//' \code{vrm_input} object.
//' @param year A variable corresponding to observation occasion, or a set of
//' such values, given in values associated with the \code{year} term used in
//' vital rate model development. Can also equal \code{"all"}, in which case
//' matrices will be estimated for all occasions. Defaults to \code{"all"}.
//' @param pop A variable designating which populations will have matrices
//' estimated. Should be set to specific population names, or to \code{"all"} if
//' all populations should have matrices estimated. Only used in raw MPMs.
//' @param patch A variable designating which patches or subpopulations will have
//' matrices estimated. Should be set to specific patch names, or to \code{"all"}
//' if matrices should be estimated for all patches. Defaults to \code{NULL}, in
//' which case patch designations are ignored.
//' @param stageframe An object of class \code{stageframe}. These objects are
//' generated by function \code{\link{sf_create}()}, and include information on
//' the size, observation status, propagule status, reproduction status,
//' immaturity status, maturity status, stage group, size bin widths, and other
//' key characteristics of each ahistorical stage. Not needed for purely
//' age-based MPMs.
//' @param supplement An optional data frame of class \code{lefkoSD} that
//' provides supplemental data that should be incorporated into the MPM. Three
//' kinds of data may be integrated this way: transitions to be estimated via the
//' use of proxy transitions, transition overwrites from the literature or
//' supplemental studies, and transition multipliers for survival and fecundity.
//' This data frame should be produced using the \code{\link{supplemental}()}
//' function. Can be used in place of or in addition to an overwrite table (see 
//' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
//' below).
//' @param overwrite An optional data frame developed with the
//' \code{\link{overwrite}()} function describing transitions to be overwritten
//' either with given values or with other estimated transitions. Note that this
//' function supplements overwrite data provided in \code{supplement}.
//' @param repmatrix An optional reproduction matrix. This matrix is composed
//' mostly of \code{0}s, with non-zero entries acting as element identifiers and
//' multipliers for fecundity (with \code{1} equaling full fecundity). If left
//' blank, and no \code{supplement} is provided, then all stages marked as
//' reproductive produce offspring at 1x that of estimated fecundity, and that
//' offspring production will yield the first stage noted as propagule or
//' immature. May be the dimensions of either a historical or an ahistorical
//' matrix. If the latter, then all stages will be used in occasion \emph{t}-1
//' for each suggested ahistorical transition. Not used in purely age-based
//' MPMs.
//' @param alive A vector of names of binomial variables corresponding to status
//' as alive (\code{1}) or dead (\code{0}) in occasions \emph{t}+1, \emph{t},
//' and \emph{t}-1, respectively. Defaults to 
//' \code{c("alive3", "alive2", "alive1")} for historical MPMs, and
//' \code{c("alive3", "alive2")} for ahistorical MPMs. Only needed for raw MPMs.
//' @param obsst A vector of names of binomial variables corresponding to
//' observation status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
//' respectively. Defaults to \code{c("obsstatus3", "obsstatus2", "obsstatus1")}
//' for historical MPMs, and \code{c("obsstatus3", "obsstatus2")} for
//' ahistorical MPMs. Only needed for raw MPMs.
//' @param size A vector of names of variables coding the primary size variable
//' in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
//' \code{c("sizea3", "sizea2", "sizea1")} for historical MPMs, and
//' \code{c("sizea3", "sizea2")} for ahistorical MPMs. Only needed for raw,
//' stage-based MPMs.
//' @param sizeb A vector of names of variables coding the secondary size
//' variable in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
//' Defaults to an empty set, assuming that secondary size is not used. Only
//' needed for raw, stage-based MPMs.
//' @param sizec A vector of names of variables coding the tertiary size
//' variable in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
//' Defaults to an empty set, assuming that tertiary size is not used. Only
//' needed for raw, stage-based MPMs.
//' @param repst A vector of names of binomial variables corresponding to
//' reproductive status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
//' respectively. Defaults to \code{c("repstatus3", "repstatus2", "repstatus1")}
//' for historical MPMs, and \code{c("repstatus3", "repstatus2")} for
//' ahistorical MPMs. Only needed for raw MPMs.
//' @param matst A vector of names of binomial variables corresponding to
//' maturity status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
//' respectively. Defaults to \code{c("matstatus3", "matstatus2", "matstatus1")}
//' for historical MPMs, and \code{c("matstatus3", "matstatus2")} for
//' ahistorical MPMs. Must be provided if building raw MPMs, and \code{stages}
//' is not provided.
//' @param fec A vector of names of variables coding for fecundity in occasions
//' \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to
//' \code{c("feca3", "feca2", "feca1")} for historical MPMs, and
//' \code{c("feca3", "feca2")} for ahistorical MPMs. Only needed for raw,
//' stage-based MPMs.
//' @param stages An optional vector denoting the names of the variables within
//' the main vertical dataset coding for the stages of each individual in
//' occasions \emph{t}+1 and \emph{t}, and \emph{t}-1, if historical. The names
//' of stages in these variables should match those used in the
//' \code{stageframe} exactly. If left blank, then \code{rlefko3()} will attempt
//' to infer stages by matching values of \code{alive}, \code{obsst},
//' \code{size}, \code{sizev}, \code{sizec}, \code{repst}, and \code{matst} to
//' characteristics noted in the associated \code{stageframe}. Only used in raw,
//' stage-based MPMs.
//' @param yearcol The variable name or column number corresponding to occasion
//' \emph{t} in the dataset. Defaults to \code{"year2"}. Only needed for raw
//' MPMs.
//' @param popcol The variable name or column number corresponding to the
//' identity of the population. Defaults to \code{"popid"} if a value is
//' provided for \code{pop}; otherwise empty. Only needed for raw MPMs.
//' @param patchcol The variable name or column number corresponding to patch in 
//' the dataset. Defaults to \code{"patchid"} if a value is provided for
//' \code{patch}; otherwise empty.  Only needed for raw MPMs.
//' @param indivcol The variable name or column number coding individual
//' identity. Only needed for raw MPMs.
//' @param agecol The variable name or column corresponding to age in time
//' \emph{t}. Defaults to \code{"obsage"}. Only used in raw age-based and
//' age-by-stage MPMs.
//' @param censorcol The variable name or column number denoting the censor
//' status. Only needed in raw MPMs, and only if \code{censor = TRUE}.
//' @param modelsuite One of three kinds of lists. The first is a
//' \code{lefkoMod} object holding the vital rate models and associated
//' metadata. Alternatively, an object of class \code{vrm_input} may be
//' provided. Finally, this argument may simply be a list of models used to
//' parameterize the MPM. In the final scenario, \code{data} and
//' \code{paramnames} must also be given, and all variable names must match
//' across all objects. If entered, then a function-based MPM will be developed.
//' Otherwise, a raw MPM will be developed. Only used in function-based MPMs.
//' @param paramnames A data frame with three columns, the first describing all
//' terms used in linear modeling, the second (must be called \code{mainparams})
//' giving the general model terms that will be used in matrix creation, and the
//' third showing the equivalent terms used in modeling (must be named
//' \code{modelparams}). Function \code{\link{create_pm}()} can be used to
//' create a skeleton \code{paramnames} object, which can then be edited. Only
//' required to build function-based MPMs if \code{modelsuite} is neither a
//' \code{lefkoMod} object nor a \code{vrm_input} object.
//' @param inda Can be a single value to use for individual covariate \code{a}
//' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1
//' in historical matrices, or a vector of such values corresponding to each
//' occasion in the dataset. Defaults to \code{NULL}. Only used in
//' function-based MPMs.
//' @param indb Can be a single value to use for individual covariate \code{b}
//' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1
//' in historical matrices, or a vector of such values corresponding to each
//' occasion in the dataset. Defaults to \code{NULL}. Only used in
//' function-based MPMs.
//' @param indc Can be a single value to use for individual covariate \code{c}
//' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1
//' in historical matrices, or a vector of such values corresponding to each
//' occasion in the dataset. Defaults to \code{NULL}. Only used in
//' function-based MPMs.
//' @param annua Can be a single value to use for annual covariate \code{a}
//' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1
//' in historical matrices, or a vector of such values corresponding to each
//' occasion in the dataset. Defaults to \code{NULL}. Only used in
//' function-based MPMs.
//' @param annub Can be a single value to use for annual covariate \code{b}
//' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1
//' in historical matrices, or a vector of such values corresponding to each
//' occasion in the dataset. Defaults to \code{NULL}. Only used in
//' function-based MPMs.
//' @param annuc Can be a single value to use for annual covariate \code{c}
//' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1
//' in historical matrices, or a vector of such values corresponding to each
//' occasion in the dataset. Defaults to \code{NULL}. Only used in
//' function-based MPMs.
//' @param dev_terms A numeric vector of 2 elements in the case of a Leslie MPM,
//' and of 14 elements in all other cases. Consists of scalar additions to the
//' y-intercepts of vital rate linear models used to estimate vital rates in
//' function-based MPMs. Defaults to \code{0} values for all vital rates.
//' @param density A numeric value indicating density value to use to propagate
//' matrices. Only needed if density is an explanatory term used in one or more
//' vital rate models. Defaults to \code{NA}. Only used in function_based MPMs.
//' @param CDF A logical value indicating whether to use the cumulative
//' distribution function to estimate size transition probabilities in
//' function-based MPMs. Defaults to \code{TRUE}, and should only be changed to
//' \code{FALSE} if approximate probabilities calculated via the midpoint method
//' are preferred.
//' @param random_inda A logical value denoting whether to treat individual
//' covariate \code{a} as a random, categorical variable. Otherwise is treated
//' as a fixed, numeric variable. Defaults to \code{FALSE}. Only used in
//' function-based MPMs.
//' @param random_indb A logical value denoting whether to treat individual
//' covariate \code{b} as a random, categorical variable. Otherwise is treated
//' as a fixed, numeric variable. Defaults to \code{FALSE}. Only used in
//' function-based MPMs.
//' @param random_indc A logical value denoting whether to treat individual
//' covariate \code{c} as a random, categorical variable. Otherwise is treated
//' as a fixed, numeric variable. Defaults to \code{FALSE}. Only used in
//' function-based MPMs.
//' @param negfec A logical value denoting whether fecundity values estimated to
//' be negative should be reset to \code{0}. Defaults to \code{FALSE}.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}. Only used in function-based MPMs.
//' @param theta_tol A numeric value used to indicate a maximum value to theta
//' as used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' Only used in function-based MPMs.
//' @param censor If \code{TRUE}, then data will be removed according to the
//' variable set in \code{censorcol}, such that only data with censor values
//' equal to \code{censorkeep} will remain. Defaults to \code{FALSE}. Only
//' used in raw MPMs.
//' @param censorkeep The value of the censor variable denoting data elements to
//' keep. Defaults to \code{0}. Only used in raw MPMs.
//' @param start_age The age from which to start the matrix. Defaults to
//' \code{NULL}, in which case age \code{1} is used if
//' \code{prebreeding = TRUE}, and age \code{0} is used if
//' \code{prebreeding = FALSE}. Only used in age-based MPMs.
//' @param last_age The final age to use in the matrix. Defaults to \code{NULL},
//' in which case the highest age in the dataset is used. Only used in age-based
//' and age-by-stage MPMs.
//' @param fecage_min The minimum age at which reproduction is possible.
//' Defaults to \code{NULL}, which is interpreted to mean that fecundity should
//' be assessed starting in the minimum age observed in the dataset. Only used
//' in age-based MPMs.
//' @param fecage_max The maximum age at which reproduction is possible.
//' Defaults to \code{NULL}, which is interpreted to mean that fecundity should
//' be assessed until the final observed age. Only used in age-based MPMs.
//' @param fectime  An integer indicating whether to estimate fecundity using
//' the variable given for \code{fec} in time \emph{t} (\code{2}) or time
//' \emph{t}+1 (\code{3}). Only used for purely age-based MPMs. Defaults to
//' \code{2}.
//' @param fecmod A scalar multiplier for fecundity. Only used for purely
//' age-based MPMs. Defaults to \code{1.0}.
//' @param cont A logical value designating whether to allow continued survival
//' of individuals past the final age noted in age-based and age-by-stage MPMs,
//' using the demographic characteristics of the final age. Defaults to
//' \code{TRUE}.
//' @param prebreeding A logical value indicating whether the life history model
//' is a pre-breeding model. Defaults to \code{TRUE}.
//' @param stage_NRasRep A logical value indicating whether to treat
//' non-reproductive individuals as reproductive. Used only in raw, stage-based
//' MPMs in cases where stage assignment must still be handled. Not used in
//' function-based MPMs, and in stage-based MPMs in which a valid \code{hfvdata}
//' class data frame with stages already assigned is provided.
//' @param sparse_output A logical value indicating whether to output matrices
//' in sparse format. Defaults to \code{FALSE}, in which case all matrices are
//' output in standard matrix format.
//' 
//' @return An object of class \code{lefkoMat}. This is a list that holds the
//' matrix projection model and all of its metadata. The structure has the
//' following elements:
//' 
//' \item{A}{A list of full projection matrices in order of sorted patches and
//' occasion times. All matrices output in R's \code{matrix} class, or in
//' the \code{dgCMatrix} class from the \code{Matrix} package if sparse.}
//' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
//' matrices output in R's \code{matrix} class, or in the \code{dgCMatrix} class
//' from the \code{Matrix} package if sparse.}
//' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
//' output in R's \code{matrix} class, or in the \code{dgCMatrix} class from the
//' \code{Matrix} package if sparse.}
//' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
//' used to create historical stage pairs. Only used in historical MPMs.}
//' \item{agestages}{A data frame showing age-stage pairs. Only used in
//' age-by-stage MPMs.}
//' \item{ahstages}{A data frame detailing the characteristics of associated
//' ahistorical stages, in the form of a modified stageframe that includes
//' status as an entry stage through reproduction. Used in all stage-based and
//' age-by-stage MPMs.}
//' \item{labels}{A data frame giving the population, patch, and year of each
//' matrix in order.}
//' \item{dataqc}{A vector showing the numbers of individuals and rows in the
//' vertical dataset used as input.}
//' \item{matrixqc}{A short vector describing the number of non-zero elements in
//' \code{U} and \code{F} matrices, and the number of annual matrices.}
//' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite}
//' input.}
//' \item{prob_out}{An optional element only added if \code{err_check = TRUE}.
//' This is a list of vital rate probability matrices, with 7 columns in the
//' order of survival, observation probability, reproduction probability, primary
//' size transition probability, secondary size transition probability, tertiary
//' size transition probability, and probability of juvenile transition to
//' maturity.}
//' \item{allstages}{An optional element only added if \code{err_check = TRUE}.
//' This is a data frame giving the values used to determine each matrix element
//' capable of being estimated.}
//' \item{data}{An optional element only added if \code{err_check = TRUE} and a
//' raw MPM is requested. This consists of the original dataset as edited by
//' this function for indexing purposes.}
//' 
//' @section General Notes:
//' This function automatically determines whether to create a raw or
//' function-based MPM given inputs supplied by the user.
//' 
//' If used, the reproduction matrix (field \code{repmatrix}) may be supplied as
//' either historical or ahistorical. If provided as historical, then
//' a historical MPM must be estimated.
//' 
//' If neither a supplement nor a reproduction matrix are used, and the MPM
//' to create is stage-based, then fecundity will be assumed to occur from all
//' reproductive stages to all propagule and immature stages.
//' 
//' @section Function-based MPM Notes:
//' Users may at times wish to estimate MPMs using a dataset incorporating
//' multiple patches or subpopulations, but without discriminating between those
//' patches or subpopulations. Should the aim of analysis be a general MPM that
//' does not distinguish these patches or subpopulations, the
//' \code{modelsearch()} run should not include patch terms.
//' 
//' Input options including multiple variable names must be entered in the order
//' of variables in occasion \emph{t}+1, \emph{t}, and \emph{t}-1. Rearranging
//' the order will lead to erroneous calculations, and will may lead to fatal
//' errors.
//' 
//' This function provides two different means of estimating the probability of
//' size transition. The midpoint method (\code{CDF = FALSE}) refers to the
//' method in which the probability is estimated by first estimating the
//' probability associated with transition from the exact size at the midpoint
//' of the size class using the corresponding probability density function, and
//' then multiplying that value by the bin width of the size class. Doak et al.
//' 2021 (Ecological Monographs) noted that this method can produce biased
//' results, with total size transitions associated with a specific size not
//' totaling to 1.0 and even specific size transition probabilities capable of
//' being estimated at values greater than 1.0. The alternative and default
//' method (\code{CDF = TRUE}) uses the cumulative density function to estimate
//' the probability of size transition as the cumulative probability of size
//' transition at the greater limit of the size class minus the cumulative
//' probability of size transition at the lower limit of the size class. This
//' latter method avoids this bias. Note, however, that both methods are exact
//' and unbiased for negative binomial and Poisson distributions.
//' 
//' Under the Gaussian and gamma size distributions, the number of estimated
//' parameters may differ between the two \code{ipm_method} settings. Because
//' the midpoint method has a tendency to incorporate upward bias in the
//' estimation of size transition probabilities, it is more likely to yield non-
//' zero values when the true probability is extremely close to 0. This will
//' result in the \code{summary.lefkoMat()} function yielding higher numbers of
//' estimated parameters than the \code{ipm_method = "CDF"} yields in some cases.
//' 
//' @examples
//' \donttest{
//' data(lathyrus)
//' 
//' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
//'   9)
//' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
//'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
//'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
//' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
//'   0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
//'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
//' 
//' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
//'   propstatus = propvector)
//' 
//' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
//'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988", 
//'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
//' 
//' lathvertln$feca2 <- round(lathvertln$feca2)
//' lathvertln$feca1 <- round(lathvertln$feca1)
//' lathvertln$feca3 <- round(lathvertln$feca3)
//' 
//' lathvertln_adults <- subset(lathvertln, stage2index > 2)
//' surv_model <- glm(alive3 ~ sizea2 + sizea1 + as.factor(patchid) +
//'   as.factor(year2), data = lathvertln_adults, family = "binomial")
//' 
//' obs_data <- subset(lathvertln_adults, alive3 == 1)
//' obs_model <- glm(obsstatus3 ~ as.factor(patchid), data = obs_data,
//'   family = "binomial")
//' 
//' size_data <- subset(obs_data, obsstatus3 == 1)
//' siz_model <- lm(sizea3 ~ sizea2 + sizea1 + repstatus1 + as.factor(patchid) +
//'   as.factor(year2), data = size_data)
//' 
//' reps_model <- glm(repstatus3 ~ sizea2 + sizea1 + as.factor(patchid) +
//'   as.factor(year2), data = size_data, family = "binomial")
//' 
//' fec_data <- subset(lathvertln_adults, repstatus2 == 1)
//' fec_model <- glm(feca2 ~ sizea2 + sizea1 + repstatus1 + as.factor(patchid),
//'   data = fec_data, family = "poisson")
//' 
//' lathvertln_juvs <- subset(lathvertln, stage2index < 3)
//' jsurv_model <- glm(alive3 ~ as.factor(patchid), data = lathvertln_juvs,
//'   family = "binomial")
//' 
//' jobs_data <- subset(lathvertln_juvs, alive3 == 1)
//' jobs_model <- glm(obsstatus3 ~ 1, family = "binomial", data = jobs_data)
//' 
//' jsize_data <- subset(jobs_data, obsstatus3 == 1)
//' jsiz_model <- lm(sizea3 ~ as.factor(year2), data = jsize_data)
//' 
//' jrepst_model <- 0
//' jmatst_model <- 1
//' 
//' mod_params <- create_pm(name_terms = TRUE)
//' mod_params$modelparams[3] <- "patchid"
//' mod_params$modelparams[4] <- "alive3"
//' mod_params$modelparams[5] <- "obsstatus3"
//' mod_params$modelparams[6] <- "sizea3"
//' mod_params$modelparams[9] <- "repstatus3"
//' mod_params$modelparams[11] <- "feca2"
//' mod_params$modelparams[12] <- "sizea2"
//' mod_params$modelparams[13] <- "sizea1"
//' mod_params$modelparams[18] <- "repstatus2"
//' mod_params$modelparams[19] <- "repstatus1"
//' 
//' used_models <- list(survival_model = surv_model, observation_model = obs_model,
//'   size_model = siz_model, sizeb_model = 1, sizec_model = 1,
//'   repstatus_model = reps_model, fecundity_model = fec_model,
//'   juv_survival_model = jsurv_model, juv_observation_model = jobs_model,
//'   juv_size_model = jsiz_model, juv_sizeb_model = 1, juv_sizec_model = 1,
//'   juv_reproduction_model = 0, juv_maturity_model = 1, paramnames = mod_params)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "mat", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "Sdl", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "Sd", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, "mat", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, "Sdl", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, "Sdl", NA, NA),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframeln, historical = TRUE)
//' 
//' # While we do not use MPMs to initialize f_projections3(), we do use MPMs to
//' # initialize functions start_input() and density_input().
//' lathmat3ln <- mpm_create(historical = TRUE, year = "all", patch = "all",
//'   data = lathvertln, stageframe = lathframeln, supplement = lathsupp3,
//'   modelsuite = used_models, reduce = FALSE)
//' }
//' 
//' @export mpm_create
// [[Rcpp::export(mpm_create)]]
Rcpp::List mpm_create(bool historical = false, bool stage = true, bool age = false,
  bool devries = false, bool reduce = false, bool simple = false,
  bool err_check = false, Nullable<RObject> data = R_NilValue, 
  Nullable<RObject> year = R_NilValue, Nullable<RObject> pop = R_NilValue,
  Nullable<RObject> patch = R_NilValue, Nullable<RObject> stageframe = R_NilValue,
  Nullable<RObject> supplement = R_NilValue, Nullable<RObject> overwrite = R_NilValue,
  Nullable<RObject> repmatrix = R_NilValue, Nullable<RObject> alive = R_NilValue,
  Nullable<RObject> obsst = R_NilValue, Nullable<RObject> size = R_NilValue,
  Nullable<RObject> sizeb = R_NilValue, Nullable<RObject> sizec = R_NilValue,
  Nullable<RObject> repst = R_NilValue, Nullable<RObject> matst = R_NilValue,
  Nullable<RObject> fec = R_NilValue, Nullable<RObject> stages = R_NilValue,
  Nullable<RObject> yearcol = R_NilValue, Nullable<RObject> popcol = R_NilValue,
  Nullable<RObject> patchcol = R_NilValue, Nullable<RObject> indivcol = R_NilValue,
  Nullable<RObject> agecol = R_NilValue, Nullable<RObject> censorcol = R_NilValue,
  
  Nullable<RObject> modelsuite = R_NilValue, Nullable<RObject> paramnames = R_NilValue,
  Nullable<RObject> inda = R_NilValue, Nullable<RObject> indb = R_NilValue,
  Nullable<RObject> indc = R_NilValue, Nullable<RObject> annua = R_NilValue,
  Nullable<RObject> annub = R_NilValue, Nullable<RObject> annuc = R_NilValue,
  Nullable<RObject> dev_terms = R_NilValue, double density = NA_REAL,
  bool CDF = true, bool random_inda = false, bool random_indb = false,
  bool random_indc = false, bool negfec = false, int exp_tol = 700, int theta_tol = 1e8,
  
  bool censor = false, Nullable<RObject> censorkeep = R_NilValue, int start_age = NA_INTEGER,
  int last_age = NA_INTEGER, int fecage_min = NA_INTEGER, int fecage_max = NA_INTEGER,
  int fectime = 2, double fecmod = 1.0, bool cont = true, bool prebreeding = true,
  bool stage_NRasRep = false, bool sparse_output = false) {
  
  bool raw {true};
  bool nodata {true};
  int data_vars_no {0};
  int data_points {0};
  StringVector data_vars;
  DataFrame data_;
  IntegerVector dataqc_ = {0, 0};
  List err_check_proxies;
  
  if (data.isNotNull()) {
    RObject data_input (data);
    
    if (is<DataFrame>(data_input)) {
      data_ = as<DataFrame>(data_input);
      
      StringVector data_class (as<StringVector>(data_.attr("class")));
      
      int no_classes {static_cast<int>(data_class.length())};
      int matches {0};
      for (int i = 0; i < no_classes; i++) {
        if (stringcompare_simple(String(data_class(i)), "hfv", false)) {
          matches++;
        }
      }
      if (matches == 0) {
        throw Rcpp::exception("This function requires a data frame in hfv format.", false);
      }
      
      data_vars = as<StringVector>(data_.attr("names"));
      data_vars_no = static_cast<int>(data_vars.length());
      data_points = static_cast<int>(data_.nrows());
      dataqc_(1) = data_points;
      nodata = false;
      
    } else {
      throw Rcpp::exception("This function requires a data frame in hfv format.", false);
    }
  }
  
  DataFrame stageframe_;
  if (stageframe.isNotNull()) {
    RObject sf_input (stageframe);
    
    if (is<DataFrame>(sf_input)) {
      stageframe_ = as<DataFrame>(sf_input);
      
      StringVector sf_class (as<StringVector>(stageframe_.attr("class")));
      
      int no_classes {static_cast<int>(sf_class.length())};
      int matches {0};
      for (int i = 0; i < no_classes; i++) {
        if (stringcompare_simple(String(sf_class(i)), "stage", false)) {
          matches++;
        }
      }
      if (matches == 0) throw Rcpp::exception("Argument stageframe is not valid.", false);
      
      StringVector sf_vars (as<StringVector>(stageframe_.attr("names")));
      
      int no_sf_vars {static_cast<int>(sf_vars.length())};
      StringVector sf_var_check {"size", "size_b", "size_c", "repstatus", "obsstatus", "matstatus", "indataset"};
      
      matches = 0;
      for (int i = 0; i < static_cast<int>(sf_var_check.length()); i++) {
        for (int j = 0; j < no_sf_vars; j++) {
          if (stringcompare_hard(String(sf_var_check(i)), String(sf_vars(j)))) {
            matches++;
          }
        }
      }
      if (matches != 7) throw Rcpp::exception("Argument stageframe is not valid.", false);
      
    } else throw Rcpp::exception("Argument stageframe is not valid.", false);
  } else {
    if (stage) throw Rcpp::exception("Valid stageframe is required for all stage-based MPMs.", false);
  }
  
  DataFrame supplement_;
  bool supplement_used {false};
  
  LefkoInputs::RObj_DFr_input_check("supplement", "lefkoSD", supplement_,
    supplement_used, true, true, supplement);
  
  DataFrame overwrite_;
  bool overwrite_used {false};
  
  if (overwrite.isNotNull()) {
    RObject overwrite_input (overwrite);
    if (is<DataFrame>(overwrite_input)) {
      overwrite_ = as<DataFrame>(overwrite_input);
      overwrite_used = true;
      
      int no_ovr_vars {static_cast<int>(overwrite_.length())};
      if (no_ovr_vars != 9) throw Rcpp::exception("Argument overwrite is not valid.", false);
      
    } else throw Rcpp::exception("Argument overwrite is not valid.", false);
  }
  
  NumericMatrix repmatrix_;
  bool repmatrix_used {false};
  
  if (repmatrix.isNotNull()) {
    RObject repmatrix_input (repmatrix);
    if (is<NumericMatrix>(repmatrix_input) || is<IntegerMatrix>(repmatrix_input)) {
      repmatrix_ = as<NumericMatrix>(repmatrix_input);
      
      int repm_cols {static_cast<int>(repmatrix_.ncol())};
      int repm_rows {static_cast<int>(repmatrix_.nrow())};
      
      if (repm_cols != repm_rows) pop_error("repmatrix", "", "", 4);
      repmatrix_used = true;
      
    } else pop_error("repmatrix", "", "", 4);
  }
  
  NumericVector dev_terms_;
  if (dev_terms.isNotNull()) {
    RObject dev_terms_input = as<RObject>(dev_terms);
    
    if (is<NumericVector>(dev_terms_input) || is<IntegerVector>(dev_terms_input)) {
      NumericVector dvi = as<NumericVector>(dev_terms_input);
      int dvi_length = static_cast<int>(dvi.length());
      
      if (age && !stage) {
        if (dvi_length != 2) {
          pop_error("dev_terms", "a numeric vector of 2 elements if developing Leslie MPMs", "", 1);
        }
      } else {
        if (dvi_length != 14) pop_error("dev_terms", "a numeric vector of 14 elements", "", 1);
      }
      dev_terms_ = dvi;
    }
  } else {
    NumericVector basic_devs = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0};
    dev_terms_ = basic_devs;
  }
  
  // Parameter identity through modelsuite, paramnames, raw variable entry
  StringVector alive_;
  IntegerVector alive_int;
  
  StringVector obsst_;
  IntegerVector obsst_int;
  
  StringVector size_;
  IntegerVector size_int;
  StringVector sizeb_;
  IntegerVector sizeb_int;
  StringVector sizec_;
  IntegerVector sizec_int;
  
  StringVector repst_;
  IntegerVector repst_int;
  
  StringVector matst_;
  IntegerVector matst_int;
  
  StringVector fec_;
  IntegerVector fec_int;
  
  bool obsst_used {false};
  bool sizeb_used {false};
  bool sizec_used {false};
  bool repst_used {false};
  bool matst_used {false};
  
  String year_var;
  int year_var_int {-1};
  StringVector mainyears_;
  
  String pop_var;
  int pop_var_int {-1};
  String pop_var_type {"n"};
  StringVector mainpops_;
  
  String patch_var;
  int patch_var_int {-1};
  String patch_var_type {"n"};
  StringVector mainpatches_;
  
  String indiv_var;
  int indiv_var_int {-1};
  
  String age_var;
  int age_var_int {-1};
  
  String censor_var;
  int censorcol_int {-1};
    
  NumericVector cs_keep_n;
  IntegerVector cs_keep_i;
  StringVector cs_keep_s;
  LogicalVector cs_keep_l;
  
  bool cs_n {false};
  bool cs_i {false};
  bool cs_s {false};
  bool cs_l {false};
  bool censorkeep_is_NA {false};
  
  bool modelsuite_provided {false};
  bool paramnames_provided {false};
  bool modelsuite_vrm {false};
  bool modelsuite_lM {false};
  List modelsuite_;
  DataFrame paramnames_;
  DataFrame mod_qc_;
  
  if (modelsuite.isNotNull()) {
    RObject modelsuite_entered = as<RObject>(modelsuite);
    
    if (is<List>(modelsuite_entered)) {
      modelsuite_ = as<List>(modelsuite_entered);
      
    } else throw Rcpp::exception("Argument modelsuite is not valid.", false);
    
    bool ms_has_class = modelsuite_.hasAttribute("class");
    StringVector modelsuite_class;
    if (ms_has_class) {
      modelsuite_class = wrap(modelsuite_.attr("class")); // altered
    }
    for (int i = 0; i < modelsuite_class.length(); i++) {
      if (stringcompare_hard(String(modelsuite_class(i)), "vrm_input")) modelsuite_vrm = true;
      if (stringcompare_hard(String(modelsuite_class(i)), "lefkoMod")) modelsuite_lM = true;
    }
    
    if (modelsuite_lM && !modelsuite_vrm) {
      // lefkoMod
      paramnames_ = modelsuite_["paramnames"];
      
      StringVector modelparams_ = paramnames_["modelparams"];
      String year_var_ = String(modelparams_(0));
      String indiv_var_ = String(modelparams_(1));
      String patch_var_ = String(modelparams_(2));
      String age_var_ = String(modelparams_(21));
      
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_hard(String(data_vars(i)), indiv_var_)) {
          indiv_var = indiv_var_;
          indiv_var_int = i;
        }
        if (stringcompare_hard(String(data_vars(i)), year_var_)) {
          year_var = year_var_;
          year_var_int = i;
        }
        if (stringcompare_hard(String(data_vars(i)), patch_var_)) {
          patch_var = patch_var_;
          patch_var_int = i;
        }
        if (stringcompare_hard(String(data_vars(i)), age_var_)) {
          age_var = age_var_;
          age_var_int = i;
        }
      }
      
      mod_qc_ = as<DataFrame>(modelsuite_["qc"]);
      
    } else if (modelsuite_vrm) {
      // vrm_input
      nodata = true;
      year_var_int = 0;
      patch_var_int = 0;
      
      DataFrame pm_new = paramnames_skeleton(false);
      
      CharacterVector parameter_names = as<CharacterVector>(pm_new["parameter_names"]);
      CharacterVector mainparams = as<CharacterVector>(pm_new["mainparams"]);
      
      CharacterVector modelparams = {"year2", "individ", "patch", "none", "none",
        "none", "none", "none", "none", "none", "none", "none", "none", "none",
        "none", "none", "none", "none", "none", "none", "none", "none", "none",
        "none", "none", "none", "none", "none", "none", "none", "none", "none",
        "none", "none", "none", "none", "none"};
      if (age) modelparams(21) = "age";
      
      DataFrame paramnames_created = DataFrame::create(_["parameter_names"] = parameter_names,
        _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
      paramnames_ = paramnames_created;
      
      modelsuite_["paramnames"] = paramnames_;
    } else {
      // List of models
      // Check modelsuite list structure
      arma::ivec used_name_vector = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1};
      StringVector ms_names_expected = {"survival_model", "observation_model",
        "size_model", "sizeb_model", "sizec_model", "repstatus_model",
        "fecundity_model", "juv_survival_model", "juv_observation_model",
        "juv_size_model", "juv_sizeb_model", "juv_sizec_model",
        "juv_reproduction_model", "juv_maturity_model", "paramnames"};
      int modelsuite_length = modelsuite_.length();
      StringVector modelsuite_names = wrap(modelsuite_.attr("names"));
      
      StringVector ms_names_true (15);
      
      for (int i = 0; i < modelsuite_length; i++) {
        if (stringcompare_simple(String(modelsuite_names(i)), "j", true)) {
          // Juvenile models
          if (stringcompare_simple(String(modelsuite_names(i)), String("sur"), true)) {
            used_name_vector(7) = i;
            ms_names_true(i) = ms_names_expected(7);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("obs"), true)) {
            used_name_vector(8) = i;
            ms_names_true(i) = ms_names_expected(8);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("size_"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("sizea"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("siza"), true)) {
            used_name_vector(9) = i;
            ms_names_true(i) = ms_names_expected(9);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("zeb"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("zb"), true)) {
            used_name_vector(10) = i;
            ms_names_true(i) = ms_names_expected(10);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("zec"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("zc"), true)) {
            used_name_vector(11) = i;
            ms_names_true(i) = ms_names_expected(11);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("rep"), true)) {
            used_name_vector(12) = i;
            ms_names_true(i) = ms_names_expected(12);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("mat"), true)) {
            used_name_vector(13) = i;
            ms_names_true(i) = ms_names_expected(13);
          }
          
        } else {
          // Adult models
          if (stringcompare_simple(String(modelsuite_names(i)), String("sur"), true)) {
            used_name_vector(0) = i;
            ms_names_true(i) = ms_names_expected(0);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("obs"), true)) {
            used_name_vector(1) = i;
            ms_names_true(i) = ms_names_expected(1);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("size_"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("sizea"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("siza"), true)) {
            used_name_vector(2) = i;
            ms_names_true(i) = ms_names_expected(2);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("zeb"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("zb"), true)) {
            used_name_vector(3) = i;
            ms_names_true(i) = ms_names_expected(3);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("zec"), true) ||
            stringcompare_simple(String(modelsuite_names(i)), String("zc"), true)) {
            used_name_vector(4) = i;
            ms_names_true(i) = ms_names_expected(4);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("rep"), true)) {
            used_name_vector(5) = i;
            ms_names_true(i) = ms_names_expected(5);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("fec"), true)) {
            used_name_vector(6) = i;
            ms_names_true(i) = ms_names_expected(6);
          } else if (stringcompare_simple(String(modelsuite_names(i)), String("par"), true)) {
            used_name_vector(14) = i;
            ms_names_true(i) = ms_names_expected(14);
            
            paramnames_ = as<DataFrame>(modelsuite_(i));
            StringVector modelparams_ = paramnames_["modelparams"];
            String indiv_var_ = String(modelparams_(1));
            String year_var_ = String(modelparams_(0));
            String patch_var_ = String(modelparams_(2));
            String age_var_ = String(modelparams_(21));
            
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), indiv_var_)) {
                indiv_var = indiv_var_;
                indiv_var_int = i;
              }
              if (stringcompare_hard(String(data_vars(i)), year_var_)) {
                year_var = year_var_;
                year_var_int = i;
              }
              if (stringcompare_hard(String(data_vars(i)), patch_var_)) {
                patch_var = patch_var_;
                patch_var_int = i;
              }
              if (stringcompare_hard(String(data_vars(i)), age_var_)) {
                age_var = age_var_;
                age_var_int = i;
              }
            }
            
            paramnames_provided = true;
          }
        }
        if (modelsuite_(i) == R_NilValue) modelsuite_(i) = 1;
      }
      
      if (paramnames.isNotNull() && !paramnames_provided) {
        RObject paramnames_entered (paramnames);
        
        if (is<DataFrame>(paramnames_entered)) {
          paramnames_ = as<DataFrame>(paramnames_entered);
        }
        StringVector modelparams_ = paramnames_["modelparams"];
        String indiv_var_ = String(modelparams_(1));
        String year_var_ = String(modelparams_(0));
        String patch_var_ = String(modelparams_(2));
        String age_var_ = String(modelparams_(21));
        
        for (int i = 0; i < data_vars_no; i++) {
          if (stringcompare_hard(String(data_vars(i)), indiv_var_)) {
            indiv_var = indiv_var_;
            indiv_var_int = i;
          }
          if (stringcompare_hard(String(data_vars(i)), year_var_)) {
            year_var = year_var_;
            year_var_int = i;
          }
          if (stringcompare_hard(String(data_vars(i)), patch_var_)) {
            patch_var = patch_var_;
            patch_var_int = i;
          }
          if (stringcompare_hard(String(data_vars(i)), age_var_)) {
            age_var = age_var_;
            age_var_int = i;
          }
        }
        paramnames_provided = true;
      } else if (!paramnames_provided) {
        throw Rcpp::exception("Argument paramnames is required with list of vital rate models.", false);
      }
      
      NumericVector one_vec = {1};
      List new_modelsuite (15);
      
      for (int i = 0; i < 14; i++) {
        if (used_name_vector(i) != -1) {
          new_modelsuite(i) = modelsuite_(used_name_vector(i));
        } else {
          new_modelsuite(i) = one_vec;
        }
      }
      
      new_modelsuite(14) = paramnames_;
      modelsuite_ = new_modelsuite;
      modelsuite_.attr("names") = ms_names_true; // new
    }
    modelsuite_provided = true;
    paramnames_provided = true;
    raw = false;
    
    if (patch_var_int == -1) {
      StringVector mainpatches = {NA_STRING};
      mainpatches_ = mainpatches;
    }
  }
  
  // Processing age
  if (agecol.isNotNull() && !paramnames_provided) {
    RObject age_input (agecol);
    
    if (is<StringVector>(age_input)) {
      StringVector age_ = as<StringVector>(age_input);
      int age_no {static_cast<int>(age_.length())};
      
      if (age_no > 1) {
        pop_error("agecol", "single string", "variable coding for age at time t", 22);
      }
      
      age_var = String(age_(0));
      
      int matches {0};
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_simple(age_var, String(data_vars(i)), false)) {
          age_var_int = i;
          matches++;
        }
      }
      if (matches != age_no) pop_error("age", "", "", 16);
      
    } else if (is<IntegerVector>(age_input) || is<NumericVector>(age_input)) {
      IntegerVector age_int_(age_input);
      int age_int_no {static_cast<int>(age_int_.length())};
      
      if (age_int_no > 1 || age_int_(0) < 1 || age_int_(0) > data_vars_no) {
        throw Rcpp::exception("Invalid entry given for age variable.", false);
      }
      
      age_var_int = age_int_(0) - 1;
      age_var = data_vars(age_var_int);
      
    } else {
      throw Rcpp::exception("Please enter the name of the variable coding age at time t.", false);
    }
  } else if (age && raw && !paramnames_provided) {
    String age_var_ {"obsage"};
    age_var = age_var_;
    
    int matches {0};
    for (int i = 0; i < data_vars_no; i++) {
      if (stringcompare_simple(age_var, String(data_vars(i)), false)) {
        age_var_int = i;
        matches++;
      }
    }
    if (matches != 1) {
      throw Rcpp::exception("Default agecol variable name is not correct. Enter correct variable.", false);
    }
  }
  
  if (age && age_var_int > -1) {
    IntegerVector data_ages = data_[age_var_int];
    int data_age_min = min(data_ages);
    int data_age_max = max(data_ages);
    
    if (IntegerVector::is_na(start_age)) {
      if (prebreeding) {
        start_age = 1;
      } else {
        start_age = 0;
      }
    }
    
    if (prebreeding && start_age == 0) {
      Rf_warningcall(R_NilValue,
        "Switching to post-breeding model to account for start_age = 0.");
      prebreeding = false;
    }
    
    if (IntegerVector::is_na(last_age)) last_age = data_age_max + 1;
    if (IntegerVector::is_na(fecage_min)) fecage_min = data_age_min;
    if (IntegerVector::is_na(fecage_max)) fecage_max = last_age;
  } else if (age) {
    if (IntegerVector::is_na(start_age)) {
      if (prebreeding) {
        start_age = 1;
      } else {
        start_age = 0;
      }
    }
    
    if (IntegerVector::is_na(last_age)) {
      throw Rcpp::exception("Age-based matrices require a valid final age.", false);
    }
    if (IntegerVector::is_na(fecage_min)) fecage_min = start_age;
    if (IntegerVector::is_na(fecage_max)) fecage_max = last_age;
  }
  
  // Edited and reformatted stageframe, supplement, and repmatrix
  DataFrame melchett_stageframe_;
  DataFrame melchett_ovtable_;
  NumericMatrix melchett_repmatrix_;
  arma::vec melchett_stageframe_size_;
  arma::vec melchett_stageframe_sizeb_;
  arma::vec melchett_stageframe_sizec_;
  arma::uvec melchett_stageframe_repst_;
  arma::uvec melchett_stageframe_obsst_;
  arma::uvec melchett_stageframe_matst_;
  arma::uvec melchett_stageframe_indataset_;
  arma::uvec melchett_stageframe_alive_;
  IntegerVector melchett_stageframe_stageid_;
  IntegerVector melchett_stageframe_group_;
  StringVector melchett_stageframe_stage_;
  
  IntegerVector maingroups_;
  
  int melchett_stageframe_length {0};
  int melchett_ovtable_length {0};
  int format_int {1};
  
  if (!historical) {
    // Ahistorical MPMs
    
    if (stage) {
      if (!age) {
        if (repmatrix_used && supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            supplement_, R_NilValue, repmatrix_, false, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (repmatrix_used && !supplement_used && !overwrite_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, R_NilValue, repmatrix_, false, false, 1);
            
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (repmatrix_used && !supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, overwrite_, repmatrix_, false, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (!repmatrix_used && supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            supplement_, R_NilValue, R_NilValue, false, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, R_NilValue, R_NilValue, false, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        }
      } else {
        // Age-by-stage MPMs
        DataFrame melchett_ovtable_temp;
        if (repmatrix_used && supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            supplement_, R_NilValue, repmatrix_, true, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (repmatrix_used && !supplement_used && !overwrite_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, R_NilValue, repmatrix_, true, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (repmatrix_used && !supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, overwrite_, repmatrix_, true, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (!repmatrix_used && supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            supplement_, R_NilValue, R_NilValue, true, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, R_NilValue, R_NilValue, true, false, 1);
          
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_temp = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        }
        
        if (supplement_used) {
          melchett_ovtable_ = LefkoMats::age_expanded(melchett_ovtable_temp,
            start_age, last_age);
        } else {
          StringVector nsst3 = {};
          IntegerVector nsa2 = {};
          NumericVector nsgr = {};
          
          melchett_ovtable_ = DataFrame::create(_["stage3"] = nsst3,
            _["stage2"] = clone(nsst3), _["stage1"] = clone(nsst3),
            _["age2"] = nsa2, _["eststage3"] = clone(nsst3),
            _["eststage2"] = clone(nsst3), _["eststage1"] = clone(nsst3),
            _["estage2"] = clone(nsa2), _["givenrate"] = nsgr,
            _["offset"] = clone(nsgr), _["multiplier"] = clone(nsgr),
            _["convtype"] = clone(nsa2), _["convtype_t12"] = clone(nsa2),
            _["pop"] = clone(nsst3), _["patch"] = clone(nsst3),
            _["year2"] = clone(nsst3));
        }
      }
      
      if (supplement_used) {
        melchett_ovtable_length = static_cast<int>(melchett_ovtable_.nrows());
        
        DataFrame mov_short = DataFrame::create(_["0"] = as<StringVector>(melchett_ovtable_[0]),
          _["1"] = as<StringVector>(melchett_ovtable_[1]),
          _["2"] = as<StringVector>(melchett_ovtable_[2]),
          _["3"] = as<IntegerVector>(melchett_ovtable_[3]),
          _["4"] = as<IntegerVector>(melchett_ovtable_[11]));
        
        if (LefkoUtils::df_duplicates(mov_short)) {
          Rf_warningcall(R_NilValue,
            "Supplement contains multiple entries for the same transitions.");
        }
      }
    } else {
      if (age) {
        // Pure age-based MPMs
        if (supplement_used) {
          melchett_ovtable_ = LefkoMats::age_expanded(supplement_, start_age, last_age);
        }
        melchett_stageframe_ = LefkoMats:: sf_leslie(start_age, last_age,
          fecage_min, fecage_max, cont);
      }
    }
  } else {
    // Historical MPMs
    if (stage) {
      if (!age) {
        if (devries) format_int = 2;
        
        if (repmatrix_used && supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            supplement_, R_NilValue, repmatrix_, false, true, format_int);
            
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (repmatrix_used && !supplement_used && !overwrite_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, R_NilValue, repmatrix_, false, true, format_int);
            
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (repmatrix_used && !supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, overwrite_, repmatrix_, false, true, format_int);
            
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else if (!repmatrix_used && supplement_used) {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            supplement_, R_NilValue, R_NilValue, false, true, format_int);
            
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
          
        } else {
          List melchett = LefkoMats::sf_reassess_internal(stageframe_,
            R_NilValue, R_NilValue, R_NilValue, false, true, format_int);
            
          melchett_stageframe_ = as<DataFrame>(melchett["stageframe"]);
          melchett_ovtable_ = as<DataFrame>(melchett["ovtable"]);
          melchett_repmatrix_ = as<NumericMatrix>(melchett["repmatrix"]);
        }
      } else throw Rcpp::exception("Age-by-stage MPMs cannot be historical.", false);
      
      melchett_ovtable_length = static_cast<int>(melchett_ovtable_.nrows());
      
      DataFrame mov_short = DataFrame::create(_["0"] = as<StringVector>(melchett_ovtable_[0]),
        _["1"] = as<StringVector>(melchett_ovtable_[1]),
        _["2"] = as<StringVector>(melchett_ovtable_[2]),
        _["3"] = as<StringVector>(melchett_ovtable_[10]));
      
      if (LefkoUtils::df_duplicates(mov_short)) {
        Rf_warningcall(R_NilValue, "Supplement contains multiple entries for s0me transitions.");
      }
    } else {
      if (age) throw Rcpp::exception("Age-based MPMs cannot be historical.", false);
    }
  }
  
  melchett_stageframe_size_ = as<arma::vec>(melchett_stageframe_["original_size"]);
  melchett_stageframe_sizeb_ = as<arma::vec>(melchett_stageframe_["original_size_b"]);
  melchett_stageframe_sizec_ = as<arma::vec>(melchett_stageframe_["original_size_c"]);
  melchett_stageframe_repst_ = as<arma::uvec>(melchett_stageframe_["repstatus"]);
  melchett_stageframe_obsst_ = as<arma::uvec>(melchett_stageframe_["obsstatus"]);
  melchett_stageframe_matst_ = as<arma::uvec>(melchett_stageframe_["matstatus"]);
  melchett_stageframe_indataset_ = as<arma::uvec>(melchett_stageframe_["indataset"]);
  melchett_stageframe_alive_ = as<arma::uvec>(melchett_stageframe_["alive"]);
  melchett_stageframe_stage_ = as<StringVector>(melchett_stageframe_["stage"]);
  melchett_stageframe_stageid_ = as<IntegerVector>(melchett_stageframe_["stage_id"]);
  melchett_stageframe_length = static_cast<int>(melchett_stageframe_alive_.n_elem);
  melchett_stageframe_group_ = as<IntegerVector>(melchett_stageframe_["group"]);
  maingroups_ = seq(min(melchett_stageframe_group_), max(melchett_stageframe_group_));
  
  if (!paramnames_provided) {
    // Processing data variables when paramnames is not provided
    if (alive.isNotNull() && raw) {
      RObject alive_input (alive);
      
      if (is<StringVector>(alive_input)) {
        alive_ = as<StringVector>(alive_input);
        int alive_no {static_cast<int>(alive_.length())};
        IntegerVector alive_int_ (alive_no);
        
        if (alive_no < 2 || alive_no > 3) pop_error("alive", "status", "", 10);
        
        int matches {0};
        for (int i = 0; i < alive_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(alive_(i)), String(data_vars(j)), false)) {
              alive_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != alive_no) pop_error("alive status", "hfv data frame", "", 9);
        alive_int = alive_int_;
        
      } else if (is<IntegerVector>(alive_input) || is<NumericVector>(alive_input)) {
        IntegerVector alive_int_(alive_input);
        
        if (min(alive_int_) < 1 || max(alive_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for alive status.", false);
        }
        
        int alive_int_no {static_cast<int>(alive_int_.length())};
        StringVector alive_string (alive_int_no);
        
        for (int i = 0; i < alive_int_no; i++) {
          alive_string(i) = data_vars(alive_int_(i) - 1);
        }
        alive_ = alive_string;
        alive_int = alive_int_;
        
      } else pop_error("alive status", "hfv data frame", "", 9);
    } else if (raw) {
      int alive_no {3};
      if (historical) {
        alive_ = {"alive3", "alive2", "alive1"};
      } else {
        alive_ = {"alive3", "alive2"};
        alive_no = 2;
      }
      
      int matches {0};
      for (int i = 0; i < alive_no; i++) {
        for (int j = 0; j < data_vars_no; j++) {
          if (stringcompare_simple(String(alive_(i)), String(data_vars(j)), false)) {
            matches++;
          }
        }
      }
      if (matches != alive_no) {
          pop_error("alive status", "", "", 16);
      }
    }
    
    // Alive status info used to subset data for raw MPMs
    if (raw) {
      NumericVector living_only = {1.0};
      StringVector alive_var_used = {alive_(1)};
      
      data_ = df_subset(data_, as<RObject>(living_only), false,
        true, false, false, true, as<RObject>(alive_var_used));
      
      data_points = static_cast<int>(data_.nrows());
    }
    
    if (obsst.isNotNull() && raw && stage) {
      RObject obsst_input (obsst);
      obsst_used = true;
      
      if (is<StringVector>(obsst_input)) {
        obsst_ = as<StringVector>(obsst_input);
        int obsst_no {static_cast<int>(obsst_.length())};
        IntegerVector obsst_int_ (obsst_no);
        
        if (obsst_no < 2 || obsst_no > 3) pop_error("obsst", "observation status", "", 10);
        
        int matches {0};
        for (int i = 0; i < obsst_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(obsst_(i)), String(data_vars(j)), false)) {
              obsst_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != obsst_no) pop_error("observation status", "", "", 16);
        obsst_int = obsst_int_;
        
      } else if (is<IntegerVector>(obsst_input) || is<NumericVector>(obsst_input)) {
        IntegerVector obsst_int_(obsst_input);
        
        if (min(obsst_int_) < 1 || max(obsst_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for observation status.", false);
        }
        
        int obsst_int_no {static_cast<int>(obsst_int_.length())};
        StringVector obsst_string (obsst_int_no);
        
        for (int i = 0; i < obsst_int_no; i++) {
          obsst_string(i) = data_vars(obsst_int_(i) - 1);
        }
        obsst_ = obsst_string;
        obsst_int = obsst_int_;
        
      } else pop_error("obsst", "observation status", "", 10);
      
    } else if (raw && stage) {
      arma::uvec mso;
      
      if (melchett_stageframe_alive_(melchett_stageframe_length - 1) == 0 && 
          melchett_stageframe_length > 2) {
        mso = melchett_stageframe_obsst_.subvec(0, (melchett_stageframe_length - 2));
      } else if (melchett_stageframe_alive_(melchett_stageframe_length - 2) == 0 &&
          melchett_stageframe_length > 2) {
        mso = melchett_stageframe_obsst_.subvec(0, (melchett_stageframe_length - 3));
      } else {
        mso = melchett_stageframe_obsst_;
      }
      
      for (int i = 0; i < static_cast<int>(mso.n_elem); i++) {
        if (IntegerVector::is_na(mso(i))) mso(i) = 0.0;
      }
      
      arma::uvec sf_obs_unique = unique(mso);
      int obs_states = static_cast<int>(sf_obs_unique.n_elem);
      
      if (obs_states > 1) {
        obsst_used = true;
        int obsst_no {3};
        
        if (historical) {
          obsst_ = {"obsstatus3", "obsstatus2", "obsstatus1"};
        } else {
          obsst_ = {"obsstatus3", "obsstatus2"};
          obsst_no = 2;
        }
        
        int matches {0};
        for (int i = 0; i < obsst_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(obsst_(i)), String(data_vars(j)), false)) {
              matches++;
            }
          }
        }
        if (matches != obsst_no) {
          pop_error("observation status", "", "", 16);
        }
      }
    }
    
    if (size.isNotNull() && raw && stage) {
      RObject size_input (size);
      
      if (is<StringVector>(size_input)) {
        size_ = as<StringVector>(size_input);
        int size_no {static_cast<int>(size_.length())};
        IntegerVector size_int_ (size_no);
        
        if (size_no < 2 || size_no > 3) pop_error("size", "primary size", "", 10);
        
        int matches {0};
        for (int i = 0; i < size_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(size_(i)), String(data_vars(j)), false)) {
              size_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != size_no) pop_error("size", "", "", 16);
        size_int = size_int_;
        
      } else if (is<IntegerVector>(size_input) || is<NumericVector>(size_input)) {
        IntegerVector size_int_(size_input);
        
        if (min(size_int_) < 1 || max(size_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for size.", false);
        }
        
        int size_int_no {static_cast<int>(size_int_.length())};
        StringVector size_string (size_int_no);
        
        for (int i = 0; i < size_int_no; i++) {
          size_string(i) = data_vars(size_int_(i) - 1);
        }
        size_ = size_string;
        size_int = size_int_;
        
      } else pop_error("size", "primary size", "", 10);
    } else if (raw && stage) {
      int size_no {3};
      if (historical) {
        size_ = {"sizea3", "sizea2", "sizea1"};
      } else {
        size_ = {"sizea3", "sizea2"};
        size_no = 2;
      }
      
      int matches {0};
      for (int i = 0; i < size_no; i++) {
        for (int j = 0; j < data_vars_no; j++) {
          if (stringcompare_simple(String(size_(i)), String(data_vars(j)), false)) {
            matches++;
          }
        }
      }
      if (matches != size_no) {
          pop_error("primary size", "", "", 16);
      }
    }
    
    if (sizeb.isNotNull() && raw && stage) {
      RObject sizeb_input (sizeb);
      sizeb_used = true;
      
      if (is<StringVector>(sizeb_input)) {
        sizeb_ = as<StringVector>(sizeb_input);
        int sizeb_no {static_cast<int>(sizeb_.length())};
        IntegerVector sizeb_int_ (sizeb_no);
        
        if (sizeb_no < 2 || sizeb_no > 3) pop_error("sizeb", "secondary size", "", 10);
        
        int matches {0};
        for (int i = 0; i < sizeb_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(sizeb_(i)), String(data_vars(j)), false)) {
              sizeb_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != sizeb_no) pop_error("sizeb", "", "", 16);
        sizeb_int = sizeb_int_;
        
      } else if (is<IntegerVector>(sizeb_input) || is<NumericVector>(sizeb_input)) {
        IntegerVector sizeb_int_(sizeb_input);
        
        if (min(sizeb_int_) < 1 || max(sizeb_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for sizeb.", false);
        }
        
        int sizeb_int_no {static_cast<int>(sizeb_int_.length())};
        StringVector sizeb_string (sizeb_int_no);
        
        for (int i = 0; i < sizeb_int_no; i++) {
          sizeb_string(i) = data_vars(sizeb_int_(i) - 1);
        }
        sizeb_ = sizeb_string;
        sizeb_int = sizeb_int_;
        
      } else pop_error("sizeb", "secondary size", "", 10);
    } else if (raw && stage) {
      arma::vec mssb;
      
      if (melchett_stageframe_alive_(melchett_stageframe_length - 1) == 0 && 
          melchett_stageframe_length > 2) {
        mssb = melchett_stageframe_sizeb_.subvec(0, (melchett_stageframe_length - 2));
        
      } else if (melchett_stageframe_alive_(melchett_stageframe_length - 2) == 0 && 
          melchett_stageframe_length > 2) {
        mssb = melchett_stageframe_sizeb_.subvec(0, (melchett_stageframe_length - 3));
        
      } else {
        mssb = melchett_stageframe_sizeb_;
      }
      
      for (int i = 0; i < static_cast<int>(mssb.n_elem); i++) {
        if (NumericVector::is_na(mssb(i))) mssb(i) = 0.0;
      }
      
      arma::vec sf_sizeb_unique = unique(mssb);
      int sizeb_states = static_cast<int>(sf_sizeb_unique.n_elem);
      
      if (sizeb_states > 1) {
        sizeb_used = true;
        int size_no {3};
        
        if (historical) {
          sizeb_ = {"sizeb3", "sizeb2", "sizeb1"};
        } else {
          sizeb_ = {"sizeb3", "sizeb2"};
          size_no = 2;
        }
        
        int matches {0};
        for (int i = 0; i < size_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(sizeb_(i)), String(data_vars(j)), false)) {
              matches++;
            }
          }
        }
        if (matches != size_no) {
          pop_error("secondary size", "", "", 16);
        }
      } else {
        if (historical) {
          sizeb_ = StringVector::create(NA_STRING, NA_STRING, NA_STRING);
        } else {
          sizeb_ = StringVector::create(NA_STRING, NA_STRING);
        }
      }
    }
    
    if (sizec.isNotNull() && raw && stage) {
      RObject sizec_input (sizec);
      sizec_used = true;
      
      if (is<StringVector>(sizec_input)) {
        sizec_ = as<StringVector>(sizec_input);
        int sizec_no {static_cast<int>(sizec_.length())};
        IntegerVector sizec_int_ (sizec_no);
        
        if (sizec_no < 2 || sizec_no > 3) pop_error("sizec", "tertiary size", "", 10);
        
        int matches {0};
        for (int i = 0; i < sizec_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(sizec_(i)), String(data_vars(j)), false)) {
              sizec_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != sizec_no) pop_error("sizec", "", "", 16);
        sizec_int = sizec_int_;
        
      } else if (is<IntegerVector>(sizec_input) || is<NumericVector>(sizec_input)) {
        IntegerVector sizec_int_(sizec_input);
        
        if (min(sizec_int_) < 1 || max(sizec_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for tertiary size.", false);
        }
        
        int sizec_int_no {static_cast<int>(sizec_int_.length())};
        StringVector sizec_string (sizec_int_no);
        
        for (int i = 0; i < sizec_int_no; i++) {
          sizec_string(i) = data_vars(sizec_int_(i) - 1);
        }
        sizec_ = sizec_string;
        sizec_int = sizec_int_;
        
      } else pop_error("sizec", "tertiary size", "", 10);
    } else if (raw && stage) {
      arma::vec mssc;
      
      if (melchett_stageframe_alive_(melchett_stageframe_length - 1) == 0 &&
          melchett_stageframe_length > 2) {
        mssc = melchett_stageframe_sizec_.subvec(0, (melchett_stageframe_length - 2));
      } else if (melchett_stageframe_alive_(melchett_stageframe_length - 2) == 0 &&
          melchett_stageframe_length > 2) {
        mssc = melchett_stageframe_sizec_.subvec(0, (melchett_stageframe_length - 3));
      } else {
        mssc = melchett_stageframe_sizec_;
      }
      
      for (int i = 0; i < static_cast<int>(mssc.n_elem); i++) {
        if (NumericVector::is_na(mssc(i))) mssc(i) = 0.0;
      }
      
      arma::vec sf_sizec_unique = unique(mssc);
      int sizec_states = static_cast<int>(sf_sizec_unique.n_elem);
      
      if (sizec_states > 1) {
        sizec_used = true;
        int size_no {3};
        
        if (historical) {
          sizec_ = {"sizec3", "sizec2", "sizec1"};
        } else {
          sizec_ = {"sizec3", "sizec2"};
          size_no = 2;
        }
        
        int matches {0};
        for (int i = 0; i < size_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(sizec_(i)), String(data_vars(j)), false)) {
              matches++;
            }
          }
        }
        if (matches != size_no) {
          pop_error("tertiary size", "", "", 16);
        }
      } else {
        if (historical) {
          sizec_ = StringVector::create(NA_STRING, NA_STRING, NA_STRING);
        } else {
          sizec_ = StringVector::create(NA_STRING, NA_STRING);
        }
      }
    }
    
    if (repst.isNotNull() && raw) {
      RObject repst_input (repst);
      repst_used = true;
      
      if (is<StringVector>(repst_input)) {
        repst_ = as<StringVector>(repst_input);
        int repst_no {static_cast<int>(repst_.length())};
        IntegerVector repst_int_ (repst_no);
        
        if (repst_no < 2 || repst_no > 3) pop_error("repst", "reproductive status", "", 10);
        
        int matches {0};
        for (int i = 0; i < repst_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(repst_(i)), String(data_vars(j)), false)) {
              repst_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != repst_no) pop_error("reproductive status", "", "", 16);
        repst_int = repst_int_;
        
      } else if (is<IntegerVector>(repst_input) || is<NumericVector>(repst_input)) {
        IntegerVector repst_int_(repst_input);
        
        if (min(repst_int_) < 1 || max(repst_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for reproductive status.", false);
        }
        
        int repst_int_no {static_cast<int>(repst_int_.length())};
        StringVector repst_string (repst_int_no);
        
        for (int i = 0; i < repst_int_no; i++) {
          repst_string(i) = data_vars(repst_int_(i) - 1);
        }
        repst_ = repst_string;
        repst_int = repst_int_;
        
      } else pop_error("repst", "reproductive status", "", 10);
    } else if (raw) {
      arma::uvec msr;
      
      if (melchett_stageframe_alive_(melchett_stageframe_length - 1) == 0 &&
          melchett_stageframe_length > 2) {
        msr = melchett_stageframe_repst_.subvec(0, (melchett_stageframe_length - 2));
      } else if (melchett_stageframe_alive_(melchett_stageframe_length - 2) == 0 &&
          melchett_stageframe_length > 2) {
        msr = melchett_stageframe_repst_.subvec(0, (melchett_stageframe_length - 3));
      } else {
        msr = melchett_stageframe_repst_;
      }
      
      for (int i = 0; i < static_cast<int>(msr.n_elem); i++) {
        if (IntegerVector::is_na(msr(i))) msr(i) = 0.0;
      }
      
      arma::uvec sf_repst_unique = unique(msr);
      int repst_states = static_cast<int>(sf_repst_unique.n_elem);
      
      if (repst_states > 1) {
        repst_used = true;
        int repst_no {3};
        
        if (historical) {
          repst_ = {"repstatus3", "repstatus2", "repstatus1"};
        } else {
          repst_ = {"repstatus3", "repstatus2"};
          repst_no = 2;
        }
        
        int matches {0};
        for (int i = 0; i < repst_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(repst_(i)), String(data_vars(j)), false)) {
              matches++;
            }
          }
        }
        if (matches != repst_no) {
          pop_error("reproductive status", "", "", 16);
        }
      }
    }
    
    if (matst.isNotNull() && raw && stage) {
      RObject matst_input (matst);
      matst_used = true;
      
      if (is<StringVector>(matst_input)) {
        matst_ = as<StringVector>(matst_input);
        int matst_no {static_cast<int>(matst_.length())};
        IntegerVector matst_int_ (matst_no);
        
        if (matst_no < 2 || matst_no > 3) pop_error("matst", "maturity status", "", 10);
        
        int matches {0};
        for (int i = 0; i < matst_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(matst_(i)), String(data_vars(j)), false)) {
              matst_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != matst_no) pop_error("maturity status", "", "", 16);
        matst_int = matst_int_;
        
      } else if (is<IntegerVector>(matst_input) || is<NumericVector>(matst_input)) {
        IntegerVector matst_int_(matst_input);
        
        if (min(matst_int_) < 1 || max(matst_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for maturity status.", false);
        }
        
        int matst_int_no {static_cast<int>(matst_int_.length())};
        StringVector matst_string (matst_int_no);
        
        for (int i = 0; i < matst_int_no; i++) {
          matst_string(i) = data_vars(matst_int_(i) - 1);
        }
        matst_ = matst_string;
        matst_int = matst_int_;
        
      } else pop_error("matst", "maturity status", "", 10);
    } else if (raw && stage) {
      arma::uvec msm;
      
      if (melchett_stageframe_alive_(melchett_stageframe_length - 1) == 0 &&
          melchett_stageframe_length > 2) {
        msm = melchett_stageframe_matst_.subvec(0, (melchett_stageframe_length - 2));
      } else if (melchett_stageframe_alive_(melchett_stageframe_length - 2) == 0 &&
          melchett_stageframe_length > 2) {
        msm = melchett_stageframe_matst_.subvec(0, (melchett_stageframe_length - 3));
      } else {
        msm = melchett_stageframe_matst_;
      }
      
      for (int i = 0; i < static_cast<int>(msm.n_elem); i++) {
        if (IntegerVector::is_na(msm(i))) msm(i) = 0.0;
      }
      
      arma::uvec sf_matst_unique = unique(msm);
      int matst_states = static_cast<int>(sf_matst_unique.n_elem);
      
      if (matst_states > 1) {
        matst_used = true;
        int matst_no {3};
        
        if (historical) {
          matst_ = {"matstatus3", "matstatus2", "matstatus1"};
        } else {
          matst_ = {"matstatus3", "matstatus2"};
          matst_no = 2;
        }
        
        int matches {0};
        for (int i = 0; i < matst_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(matst_(i)), String(data_vars(j)), false)) {
              matches++;
            }
          }
        }
        if (matches != matst_no) {
          pop_error("maturity status", "", "", 16);
        }
      }
    }
    
    if (fec.isNotNull() && raw) {
      RObject fec_input (fec);
      
      if (is<StringVector>(fec_input)) {
        fec_ = as<StringVector>(fec_input);
        int fec_no {static_cast<int>(fec_.length())};
        IntegerVector fec_int_ (fec_no);
        
        if (fec_no < 2 || fec_no > 3) pop_error("fec", "fecundity", "", 10);
        
        int matches {0};
        for (int i = 0; i < fec_no; i++) {
          for (int j = 0; j < data_vars_no; j++) {
            if (stringcompare_simple(String(fec_(i)), String(data_vars(j)), false)) {
              fec_int_(i) = j;
              matches++;
            }
          }
        }
        if (matches != fec_no) pop_error("fecundity", "", "", 16);
        fec_int = fec_int_;
        
      } else if (is<IntegerVector>(fec_input) || is<NumericVector>(fec_input)) {
        IntegerVector fec_int_(fec_input);
        
        if (min(fec_int_) < 1 || max(fec_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for fecundity.", false);
        }
        
        int fec_int_no {static_cast<int>(fec_int_.length())};
        StringVector fec_string (fec_int_no);
        
        for (int i = 0; i < fec_int_no; i++) {
          fec_string(i) = data_vars(fec_int_(i) - 1);
        }
        fec_ = fec_string;
        fec_int = fec_int_;
        
      } else pop_error("fec", "fecundity", "", 10);
    } else if (raw) {
      int fec_no {3};
      if (historical) {
        fec_ = {"feca3", "feca2", "feca1"};
      } else {
        fec_ = {"feca3", "feca2"};
        fec_no = 2;
      }
      
      IntegerVector fec_int_(fec_no);
      int matches {0};
      for (int i = 0; i < fec_no; i++) {
        for (int j = 0; j < data_vars_no; j++) {
          if (stringcompare_simple(String(fec_(i)), String(data_vars(j)), false)) {
            fec_int_(i) = j;
            matches++;
          }
        }
      }
      fec_int = fec_int_;
      
      if (matches != fec_no) {
          pop_error("fecundity", "", "", 16);
      }
    }
    
    // Extra control variables for raw MPMs
    if (yearcol.isNotNull()) {
      RObject year_input (yearcol);
      
      if (is<StringVector>(year_input)) {
        StringVector year_ = as<StringVector>(year_input);
        int year_no {static_cast<int>(year_.length())};
        
        if (year_no > 1) {
          pop_error("year", "single string", "variable coding for time t.", 22);
        }
        
        year_var = String(year_(0));
        
        int matches {0};
        for (int i = 0; i < data_vars_no; i++) {
          if (stringcompare_simple(year_var, String(data_vars(i)), false)) {
            year_var_int = i;
            matches++;
          }
        }
        if (matches != year_no) pop_error("time t", "", "", 16);
        
      } else if (is<IntegerVector>(year_input) || is<NumericVector>(year_input)) {
        IntegerVector year_int_(year_input);
        int year_int_no {static_cast<int>(year_int_.length())};
        
        if (year_int_no > 1 || year_int_(0) < 1 || year_int_(0) > data_vars_no) {
          throw Rcpp::exception("Invalid entry given for year at time t.", false);
        }
        
        year_var_int = year_int_(0) - 1;
        year_var = data_vars(year_var_int);
        
      } else {
        throw Rcpp::exception("Please enter string with name of variable coding for time t.",
          false);
      }
      StringVector data_year = as<StringVector>(data_[year_var_int]);
      StringVector mainyears = sort_unique(data_year);
      mainyears_ = mainyears;
      
    } else {
      // Need to adjust this for vrm_input
      String year_default {"year2"};
      
      int matches {0};
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_simple(year_default, String(data_vars(i)), false)) {
          year_var_int = i;
          matches++;
        }
      }
      if (matches != 1) {
        throw Rcpp::exception("Could not locate variable coding for time t.", false);
      }
      
      StringVector data_year = as<StringVector>(data_[year_var_int]);
      StringVector mainyears = sort_unique(data_year);
      mainyears_ = mainyears;
    }
    
    if (popcol.isNotNull()) {
      RObject pop_input (popcol);
      
      if (!raw) {
        throw Rcpp::exception("Function-based MPMs cannot handle population terms.", false);
      }
      
      if (is<StringVector>(pop_input)) {
        StringVector pop_ = as<StringVector>(pop_input);
        int pop_no {static_cast<int>(pop_.length())};
        
        if (pop_no > 1) {
          pop_error("popcol", "single string", "variable coding for population identity", 22);
        }
        
        pop_var = String(pop_(0));
        
        int matches {0};
        for (int i = 0; i < data_vars_no; i++) {
          if (stringcompare_simple(pop_var, String(data_vars(i)), false)) {
            pop_var_int = i;
            matches++;
          }
        }
        if (matches != pop_no) pop_error("population", "", "", 16);
        pop_var_type = "s";
        
      } else if (is<NumericVector>(pop_input)) {
        NumericVector pop_int_(pop_input);
        int pop_int_no {static_cast<int>(pop_int_.length())};
        
        if (pop_int_no > 1 || pop_int_(0) < 1 || pop_int_(0) > data_vars_no) {
          throw Rcpp::exception("Invalid entry given for population identity.", false);
        }
        
        pop_var_int = pop_int_(0) - 1;
        pop_var = data_vars(pop_var_int);
        pop_var_type = "d";
        
      } else if (is<IntegerVector>(pop_input)) {
        IntegerVector pop_int_(pop_input);
        
        int pop_int_no {static_cast<int>(pop_int_.length())};
        
        if (pop_int_no > 1 || pop_int_(0) < 1 || pop_int_(0) > data_vars_no) {
          throw Rcpp::exception("Invalid entry given for population identity.", false);
        }
        
        pop_var_int = pop_int_(0) - 1;
        pop_var = data_vars(pop_var_int);
        
        bool pop_fac = pop_int_.hasAttribute("levels");
        if (pop_fac) {
          pop_var_type = "f";
        } else {
          pop_var_type = "i";
        }
        
      } else {
        throw Rcpp::exception("Please enter string with name of variable coding for population.", false);
      }
      
      StringVector data_pops = as<StringVector>(data_[pop_var_int]);
      StringVector mainpops = sort_unique(data_pops);
      mainpops_ = mainpops;
    }
    
    if (patchcol.isNotNull()) {
      RObject patch_input (patchcol);
      
      if (is<StringVector>(patch_input)) {
        StringVector patch_ = as<StringVector>(patch_input);
        int patch_no {static_cast<int>(patch_.length())};
        
        if (patch_no > 1) {
          pop_error("patchcol", "single string", "variable coding for patch identity.", 22);
        }
        
        patch_var = String(patch_(0));
        
        int matches {0};
        for (int i = 0; i < data_vars_no; i++) {
          if (stringcompare_simple(patch_var, String(data_vars(i)), false)) {
            patch_var_int = i;
            matches++;
          }
        }
        if (matches != patch_no) pop_error("patch", "", "", 16);
          patch_var_type = "s";
        
      } else if (is<NumericVector>(patch_input)) {
        NumericVector patch_int_(patch_input);
        int patch_int_no {static_cast<int>(patch_int_.length())};
        
        if (patch_int_no > 1 || patch_int_(0) < 1 || patch_int_(0) > data_vars_no) {
          throw Rcpp::exception("Invalid entry given for patch identity.", false);
        }
        
        patch_var_int = patch_int_(0) - 1;
        patch_var = data_vars(patch_var_int);
        
        patch_var_type = "d";
        
      } else if (is<IntegerVector>(patch_input)) {
        IntegerVector patch_int_(patch_input);
        int patch_int_no {static_cast<int>(patch_int_.length())};
        
        if (patch_int_no > 1 || patch_int_(0) < 1 || patch_int_(0) > data_vars_no) {
          throw Rcpp::exception("Invalid entry given for patch identity.", false);
        }
        
        patch_var_int = patch_int_(0) - 1;
        patch_var = data_vars(patch_var_int);
        
        bool patch_fac = patch_int_.hasAttribute("levels");
        if (patch_fac) {
          patch_var_type = "f";
        } else {
          patch_var_type = "i";
        }
        
      } else {
        throw Rcpp::exception("Please enter string with name of variable coding for patch .", false);
      }
      
      StringVector data_patch = as<StringVector>(data_[patch_var_int]);
      StringVector mainpatches = sort_unique(data_patch);
      mainpatches_ = mainpatches;
    } else {
      StringVector mainpatches = {NA_STRING};
      mainpatches_ = mainpatches;
    }
    
    if (indivcol.isNotNull()) {
      RObject indiv_input (indivcol);
      
      if (is<StringVector>(indiv_input)) {
        StringVector indiv_ = as<StringVector>(indiv_input);
        int indiv_no {static_cast<int>(indiv_.length())};
        
        if (indiv_no > 1) {
          pop_error("indivcol", "single string", "variable coding for individual identity.", 22);
        }
        
        indiv_var = String(indiv_(0));
        
        int matches {0};
        for (int i = 0; i < data_vars_no; i++) {
          if (stringcompare_simple(indiv_var, String(data_vars(i)), false)) {
            indiv_var_int = i;
            matches++;
          }
        }
        if (matches != indiv_no) pop_error("individual identity", "", "", 16);
        
      } else if (is<IntegerVector>(indiv_input) || is<NumericVector>(indiv_input)) {
        IntegerVector indiv_int_(indiv_input);
        int indiv_int_no {static_cast<int>(indiv_int_.length())};
        
        if (indiv_int_no > 1 || indiv_int_(0) < 1 || indiv_int_(0) > data_vars_no) {
          throw Rcpp::exception("Invalid entry given for individual identity.", false);
        }
        
        indiv_var_int = indiv_int_(0) - 1;
        indiv_var = data_vars(indiv_var_int);
        
      } else {
        throw Rcpp::exception("Please enter string with name of variable coding for individual.", false);
      }
    }
    
    if (censorcol.isNotNull() && censor) {
      RObject censorcol_input (censorcol);
      
      if (is<StringVector>(censorcol_input)) {
        StringVector censorcol_ = as<StringVector>(censorcol_input);
        int censorcol_no {static_cast<int>(censorcol_.length())};
        
        if (censorcol_no > 1) {
          pop_error("censorcol", "a single variable", "", 1);
        }
        
        String censorcol_0 (censorcol_(0));
        censor_var = censorcol_0;
        
        int matches {0};
        for (int j = 0; j < data_vars_no; j++) {
          if (stringcompare_simple(censorcol_0, String(data_vars(j)), false)) {
            censorcol_int = j;
            matches++;
          }
        }
        if (matches != 1) {
          String eat_my_shorts = "Censor variable either does not match entered hfv data ";
          String eat_my_shorts1 = "frame, or matches more than one variable name.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
      } else if (is<IntegerVector>(censorcol_input) || is<NumericVector>(censorcol_input)) {
        IntegerVector censorcol_int_(censorcol_input);
        
        if (min(censorcol_int_) < 1 || max(censorcol_int_) > data_vars_no) {
          throw Rcpp::exception("Invalid entries given for censorcol variable.", false);
        }
        
        int censorcol_int_no {static_cast<int>(censorcol_int_.length())};
        if (censorcol_int_no != 1) {
          throw Rcpp::exception("Please enter a single variable name for option censorcol.", false);
        }
        censorcol_int = censorcol_int_(0) - 1;
        censor_var = data_vars(censorcol_int);
        
      } else {
        throw Rcpp::exception("Please enter name of variable coding for censorcol.", false);
      }
    } else if (censor) {
      throw Rcpp::exception("Please enter name of variable coding for censorcol.", false);
    }
    
    if (censorkeep.isNotNull() && censor) {
      RObject censorkeep_input (censorkeep);
      
      if (is<NumericVector>(censorkeep_input)) {
        cs_keep_n = as<NumericVector>(censorkeep_input);
        cs_n = true;
        
        if (NumericVector::is_na(cs_keep_n(0))) censorkeep_is_NA = true;
        
      } else if (is<IntegerVector>(censorkeep_input)) {
        cs_keep_i = as<IntegerVector>(censorkeep_input);
        cs_i = true;
        
        if (IntegerVector::is_na(cs_keep_i(0))) censorkeep_is_NA = true;
        
      } else if (is<StringVector>(censorkeep_input)) {
        cs_keep_s = as<StringVector>(censorkeep_input);
        cs_s = true;
        
        if (StringVector::is_na(cs_keep_s(0))) censorkeep_is_NA = true;
        
      } else if (is<LogicalVector>(censorkeep_input)) {
        cs_keep_l = as<LogicalVector>(censorkeep_input);
        cs_l = true;
        
        if (LogicalVector::is_na(cs_keep_l(0))) censorkeep_is_NA = true;
        
      } else {
        throw Rcpp::exception("Argument censorkeep is not a recognized input type.", false);
      }
      
      StringVector data_censor (as<StringVector>(data_[censorcol_int]));
    }
    
    if (raw && censor) {
      IntegerVector cscol_vec = {censorcol_int};
      
      if (censorkeep_is_NA) {
        if (cs_n) {
          data_ = df_subset(data_, as<RObject>(cs_keep_n), true,
            true, false, false, true, as<RObject>(cscol_vec));
        } else if (cs_i) {
          data_ = df_subset(data_, as<RObject>(cs_keep_i), true,
            true, false, false, true, as<RObject>(cscol_vec));
        } else if (cs_l) {
          data_ = df_subset(data_, as<RObject>(cs_keep_l), true,
            true, false, false, true, as<RObject>(cscol_vec));
        } else if (cs_s) {
          data_ = df_subset(data_, as<RObject>(cs_keep_s), true,
            true, false, false, true, as<RObject>(cscol_vec));
        }
      } else {
        if (cs_n) {
          data_ = df_subset(data_, as<RObject>(cs_keep_n), false,
            true, false, false, true, as<RObject>(cscol_vec));
        } else if (cs_i) {
          data_ = df_subset(data_, as<RObject>(cs_keep_i), false,
            true, false, false, true, as<RObject>(cscol_vec));
        } else if (cs_l) {
          data_ = df_subset(data_, as<RObject>(cs_keep_l), false,
            true, false, false, true, as<RObject>(cscol_vec));
        } else if (cs_s) {
          data_ = df_subset(data_, as<RObject>(cs_keep_s), false,
            true, false, false, true, as<RObject>(cscol_vec));
        }
      }
      data_points = static_cast<int>(data_.nrows());
      if (data_points == 0) throw Rcpp::exception("Data censoring led to empty data subsets.", false);
    }
  } else if (nodata && modelsuite_provided) {
    // Section for vrm_input
    DataFrame year_frame = as<DataFrame>(modelsuite_["year_frame"]);
    DataFrame patch_frame = as<DataFrame>(modelsuite_["patch_frame"]);
    
    mainyears_ = as<StringVector>(year_frame["years"]);
    mainpatches_ = as<StringVector>(patch_frame["patches"]);
  } else {
    // Guided by paramnames when both data and modelsuites are provided
    if (pop_var_int > -1) {
      StringVector data_pops = as<StringVector>(data_[pop_var_int]);
      StringVector mainpops = sort_unique(data_pops);
      mainpops_ = mainpops;
    }
    
    if (patch_var_int > -1) {
      StringVector data_patch = as<StringVector>(data_[patch_var_int]);
      StringVector mainpatches = sort_unique(data_patch);
      mainpatches_ = mainpatches;
    }
    
    if (year_var_int > -1) { 
      StringVector data_year = as<StringVector>(data_[year_var_int]);
      StringVector mainyears = sort_unique(data_year);
      mainyears_ = mainyears;
    }
  }
  
  if (NumericVector::is_na(density)) {
    density = 0.0;
  }
  
  // Count individuals for dataqc
  if (!nodata) {
    int used_indiv_var = indiv_var_int;
    
    if (indiv_var_int == -1) {
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_simple(String(data_vars(i)), "indiv", false)) {
          used_indiv_var = i;
        }
      }
    }
    
    if (used_indiv_var > -1) {
      StringVector data_indiv = as<StringVector>(data_[used_indiv_var]);
      StringVector mainindivs = unique(data_indiv);
      int no_indivs = static_cast<int>(mainindivs.length());
      
      dataqc_(0) = no_indivs;
    }
  }
  
  // LOY parameters - developed for all MPMs
  StringVector year_;
  
  if (year.isNotNull()) {
    StringVector year_input(year);
    
    if (year_var_int == -1) {
      throw Rcpp::exception("Please input name of variable coding for time t.", false);
    }
    
    if (stringcompare_simple(String(year_input(0)), "all", false)) {
      if (raw) {
        int no_mainyears = static_cast<int>(mainyears_.length());
        int hist_adj {0};
        if (historical) hist_adj = 1;
        
        StringVector chosen_years (no_mainyears - (hist_adj)); // Used to be hist_adj + 1
        for (int i = 0; i < (no_mainyears - (hist_adj)); i++) { // Used to be hist_adj + 1
          chosen_years(i) = mainyears_(i + hist_adj);
        }
        year_ = chosen_years;
        
      } else {
        year_ = mainyears_;
      }
      
    } else {
      for (int i = 0; i < year_input.length(); i++) {
        int matches {0};
        for (int j = 0; j < mainyears_.length(); j++) {
          if (stringcompare_simple(String(year_input(i)), String(mainyears_(j)), false)) matches++;
        }
        if (matches == 0) throw Rcpp::exception("Year value not found in data.", false);
      }
      
      year_ = year_input;
    }
  } else {
    if (year_var_int == -1) {
      throw Rcpp::exception("Please input name of variable coding for time t.", false);
    }
    
    if (raw) {
      int no_mainyears = static_cast<int>(mainyears_.length());
      int hist_adj {0};
      if (historical) hist_adj = 1;
      
      StringVector chosen_years (no_mainyears - hist_adj);
      for (int i = 0; i < (no_mainyears - hist_adj); i++) {
        chosen_years(i) = mainyears_(i + hist_adj);
      }
      year_ = chosen_years;
      
    } else {
      year_ = mainyears_;
    }
  }
  
  StringVector pop_;
  
  if (pop.isNotNull()) {
    StringVector pop_input(pop);
    
    if (pop_var_int == -1 && !paramnames_provided) {
      String pop_var_default {"popid"};
      
      int matches {0};
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_simple(String(data_vars(i)), pop_var_default, false)) {
          pop_var = String(data_vars(i));
          pop_var_int = i;
          matches++;
        }
      }
      if (matches != 1) {
        throw Rcpp::exception("Please input name of variable coding for population.", false);
      }
      
      StringVector data_pop = as<StringVector>(data_[pop_var_int]);
      StringVector mainpops = unique(data_pop);
      mainpops_ = mainpops;
    }
    
    if (StringVector::is_na(pop_input(0))) pop_input(0) = "all";
    
    if (stringcompare_simple(String(pop_input(0)), "all", false)) {
      pop_ = mainpops_;
      
    } else {
      for (int i = 0; i < pop_input.length(); i++) {
        int matches {0};
        for (int j = 0; j < mainpops_.length(); j++) {
          if (stringcompare_simple(String(pop_input(i)), String(mainpops_(j)), false)) matches++;
        }
        if (matches == 0) {
          throw Rcpp::exception("Population identity value not found in data.", false);
        }
      }
      pop_ = pop_input;
    }
  } else if (pop_var_int != -1) {
    pop_ = mainpops_;
  }
  
  StringVector patch_;
  
  if (patch.isNotNull()) {
    StringVector patch_input(patch);
    
    if (patch_var_int == -1 && !paramnames_provided) {
      String patch_var_default {"patchid"};
      
      int matches {0};
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_simple(String(data_vars(i)), patch_var_default, false)) {
          patch_var = String(data_vars(i));
          patch_var_int = i;
          matches++;
        }
      }
      if (matches != 1) {
        throw Rcpp::exception("Please input name of variable coding for patch.", false);
      }
      
      StringVector data_patch = as<StringVector>(data_[patch_var_int]);
      StringVector mainpatches = unique(data_patch);
      mainpatches_ = mainpatches;
    }
    
    if (StringVector::is_na(patch_input(0))) patch_input(0) = "all";
    
    if (stringcompare_simple(String(patch_input(0)), "all", false)) {
      patch_ = mainpatches_;
      
    } else {
      for (int i = 0; i < patch_input.length(); i++) {
        int matches {0};
        for (int j = 0; j < mainpatches_.length(); j++) {
          if (stringcompare_simple(String(patch_input(i)), String(mainpatches_(j)), false)) matches++;
        }
        if (matches == 0) {
          throw Rcpp::exception("Patch identity value not found in data.", false);
        }
      }
      patch_ = patch_input;
      
    }
  } else if (patch_var_int != -1) {
    patch_ = mainpatches_;
  }
  
  // LOY (list of years) calculator
  int no_pops {static_cast<int>(pop_.length())};
  int no_patches {static_cast<int>(patch_.length())};
  int no_years {static_cast<int>(year_.length())};
  
  IntegerVector poporder_;
  IntegerVector patchorder_;
  IntegerVector yearorder_;
  
  bool loy_pop_used {false};
  bool loy_patch_used {false};
  
  DataFrame list_of_years;
  DataFrame labels_;
  StringVector loy_pop_;
  StringVector loy_patch_;
  StringVector loy_year2_;
  
  {
    StringVector pop_new;
    if (no_pops == 0) {
      StringVector pop_new_default {"1"};
      pop_new = pop_new_default;
      no_pops = 1;
      
    } else {
      pop_new = pop_;
      loy_pop_used = true;
    }
    
    StringVector patch_new;
    if (no_patches == 0) {
      StringVector patch_new_default {"1"};
      patch_new = patch_new_default;
      no_patches = 1;
      
    } else {
      patch_new = patch_;
      loy_patch_used = true;
    }
    
    StringVector year_new;
    if (no_years == 0) {
      StringVector year_new_default {"1"};
      year_new = year_new_default;
      no_years = 1;
      
    } else {
      year_new = year_;
    }
    
    int loy_length = no_pops * no_patches * no_years;
    
    StringVector loy_pop (loy_length);
    StringVector loy_patch (loy_length);
    StringVector loy_year2 (loy_length);
    
    IntegerVector poporder (loy_length);
    IntegerVector patchorder (loy_length);
    IntegerVector yearorder (loy_length);
    
    for (int i = 0; i < no_pops; i++) {
      for (int j = 0; j < no_patches; j++) {
        for (int k = 0; k < no_years; k++) {
          loy_pop((i * no_patches * no_years) + (j * no_years) + k) = pop_new(i);
          loy_patch((i * no_patches * no_years) + (j * no_years) + k) = patch_new(j);
          loy_year2((i * no_patches * no_years) + (j * no_years) + k) = year_new(k);
          
          poporder((i * no_patches * no_years) + (j * no_years) + k) = i;
          patchorder((i * no_patches * no_years) + (j * no_years) + k) = j;
          yearorder((i * no_patches * no_years) + (j * no_years) + k) = k;
        }
      }
    }
    
    DataFrame loy = DataFrame::create(_["pop"] = loy_pop, _["patch"] = loy_patch,
      _["year2"] = loy_year2, _["poporder"] = poporder, _["patchorder"] = patchorder,
      _["yearorder"] = yearorder);
    list_of_years = loy;
    DataFrame labels = DataFrame::create(_["pop"] = loy_pop, _["patch"] = loy_patch,
      _["year2"] = loy_year2);
    labels_ = labels;
    
    poporder_ = poporder;
    patchorder_ = patchorder;
    yearorder_ = yearorder;
    
    loy_pop_ = loy_pop;
    loy_patch_ = loy_patch;
    loy_year2_ = loy_year2;
  }
  
  // Stages in raw MPMs
  StringVector stages_;
  IntegerVector stages_int;
  
  IntegerVector new_stageid3 (data_points);
  IntegerVector new_stageid2 (data_points);
  IntegerVector new_stageid1 (data_points);
  StringVector new_stage3 (data_points);
  StringVector new_stage2 (data_points);
  StringVector new_stage1 (data_points);
  
  bool new_stages_needed {false};
  bool new_stage_indices_needed {false}; // For cases with stage calls without indices
  
  if (stages.isNotNull() && raw && stage) {
    RObject stages_input (stages);
    
    if (is<StringVector>(stages_input)) {
      stages_ = as<StringVector>(stages_input);
      int stages_no {static_cast<int>(stages_.length())};
      IntegerVector stages_int_ (stages_no);
      
      if (stages_no < 2 || stages_no > 3) pop_error("stages", "stage identity", "", 10);
      
      int matches {0};
      for (int i = 0; i < stages_no; i++) {
        for (int j = 0; j < data_vars_no; j++) {
          if (stringcompare_simple(String(stages_(i)), String(data_vars(j)), false)) {
            stages_int_(i) = j;
            matches++;
          }
        }
      }
      if (matches != stages_no) {
        throw Rcpp::exception("Names in argument stages do not match data.", false);
      }
      stages_int = stages_int_;
      
    } else if (is<IntegerVector>(stages_input) || is<NumericVector>(stages_input)) {
      IntegerVector stages_int_(stages_input);
      
      if (min(stages_int_) < 1 || max(stages_int_) > data_vars_no) {
        throw Rcpp::exception("Invalid entries given for option stages.", false);
      }
      
      int stages_int_no {static_cast<int>(stages_int_.length())};
      StringVector stages_string (stages_int_no);
      
      for (int i = 0; i < stages_int_no; i++) {
        stages_string(i) = data_vars(stages_int_(i) - 1);
      }
      stages_ = stages_string;
      stages_int = stages_int_;
      
    } else pop_error("stages", "stage identity", "", 10);
    
    new_stage3 = as<StringVector>(data_[stages_int(0)]);
    new_stage2 = as<StringVector>(data_[stages_int(1)]);
    if (historical) new_stage1 = as<StringVector>(data_[stages_int(2)]);
    
    for (int i = 0; i < data_points; i++) {
      for (int j = 0; j < melchett_stageframe_length; j++) {
        if (stringcompare_hard(String(new_stage3(i)), String(melchett_stageframe_stage_(j)))) {
          new_stageid3(i) = j + 1;
        } else if (stringcompare_hard(String(new_stage3(i)), "NotAlive") || 
            stringcompare_hard(String(new_stage3(i)), "Dead")) {
          new_stageid3(i) = melchett_stageframe_length;
        }
        if (stringcompare_hard(String(new_stage2(i)), String(melchett_stageframe_stage_(j)))) {
          new_stageid2(i) = j + 1;
        } else if (stringcompare_hard(String(new_stage2(i)), "NotAlive") || 
            stringcompare_hard(String(new_stage2(i)), "Dead")) {
          new_stageid2(i) = melchett_stageframe_length;
        }
        
        if (historical) {
          if (stringcompare_hard(String(new_stage1(i)), String(melchett_stageframe_stage_(j)))) {
            new_stageid1(i) = j + 1;
          } else if (stringcompare_hard(String(new_stage1(i)), "NotAlive") || 
              stringcompare_hard(String(new_stage1(i)), "Dead")) {
            new_stageid1(i) = melchett_stageframe_length;
          }
        }
      }
    }
    
    // Check stage indices
    int ind_matches {0};
    for (int i = 0; i < data_vars_no; i++) {
      if (stringcompare_simple(String(data_vars(i)), "index"), false) {
        if (is<IntegerVector>(data_[i])) ind_matches++;
      }
    }
    if (ind_matches < 2 || ind_matches > 3) new_stage_indices_needed = true;
    
    StringVector d_rn = data_.attr("row.names");
    
    data_["usedstage3"] = new_stage3;
    data_["usedstage2"] = new_stage2;
    if (historical) data_["usedstage1"] = new_stage1;
    
    data_["index3"] = new_stageid3 - 1;
    data_["index2"] = new_stageid2 - 1;
    if (historical) data_["index1"] = new_stageid1 - 1;
    
    data_.attr("class") = "data.frame";
    data_.attr("row.names") = d_rn;
    
  } else if (raw && stage) {
    int stages_no {3};
    
    // Check for default values, move on to assignment if nothing found
    IntegerVector stages_int;
    if (historical) {
      IntegerVector stages3_int {-1, -1, -1};
      stages_int = stages3_int;
      
      StringVector stages3 {"stage3", "stage2", "stage1"};
      stages_ = stages3;
    } else {
      IntegerVector stages2_int {-1, -1};
      stages_int = stages2_int;
      
      StringVector stages2 {"stage3", "stage2"};
      stages_ = stages2;
      stages_no = 2;
    }
    
    int matches {0};
    for (int i = 0; i < stages_no; i++) {
      for (int j = 0; j < data_vars_no; j++) {
        if (stringcompare_simple(String(stages_(i)), String(data_vars(j)), false)) {
          if ((historical && matches < 4) || (!historical && matches < 3)) {
            stages_int(matches) = j;
          } else new_stages_needed = true;
          
          matches++;
        }
      }
    }
    
    if (!new_stages_needed) {
      IntegerVector unique_elems = unique(stages_int);
      
      for (int i = 0; i < static_cast<int>(unique_elems.length()); i++) {
        if (stages_int(i) == -1) new_stages_needed = true;
      }
    }
    
    if (!new_stages_needed) {
      StringVector stages3_trial = as<StringVector>(data_[stages_int(0)]);
      StringVector stages3_trial_unique = unique(stages3_trial);
      
      int s3_tu_length = stages3_trial_unique.length();
      if (s3_tu_length < 2) new_stages_needed = true;
    }
    
    // Stage calls
    IntegerVector msf_stageid = as<IntegerVector>(melchett_stageframe_["stage_id"]);
    StringVector msf_stage = as<StringVector>(melchett_stageframe_["stage"]);
    
    if (new_stages_needed) {
      // Imports all needed data from hfv data frame
      arma::uvec data_alive3_ = as<arma::uvec>(data_[String(alive_(0))]);
      arma::uvec data_alive2_ = as<arma::uvec>(data_[String(alive_(1))]);
      arma::uvec data_alive1_;
      
      arma::vec data_size3_ = as<arma::vec>(data_[String(size_(0))]);
      arma::vec data_size2_ = as<arma::vec>(data_[String(size_(1))]);
      arma::vec data_size1_;
      
      if (historical) {
        data_alive1_ = as<arma::uvec>(data_[String(alive_(2))]);
        data_size1_ = as<arma::vec>(data_[String(size_(2))]);
      }
      
      arma::uvec data_obsst3_;
      arma::uvec data_obsst2_;
      arma::uvec data_obsst1_;
      
      if (repst_used) {
        data_obsst3_ = as<arma::uvec>(data_[String(obsst_(0))]);
        data_obsst2_ = as<arma::uvec>(data_[String(obsst_(1))]);

        if (historical) {
          data_obsst1_ = as<arma::uvec>(data_[String(obsst_(2))]);
        }
      }
      
      arma::vec data_sizeb3_;
      arma::vec data_sizeb2_;
      arma::vec data_sizeb1_;
      
      if (sizeb_used) {
        data_sizeb3_ = as<arma::vec>(data_[String(sizeb_(0))]);
        data_sizeb2_ = as<arma::vec>(data_[String(sizeb_(1))]);

        if (historical) {
          data_sizeb1_ = as<arma::vec>(data_[String(sizeb_(2))]);
        }
      }
      
      arma::vec data_sizec3_;
      arma::vec data_sizec2_;
      arma::vec data_sizec1_;
      
      if (sizec_used) {
        data_sizec3_ = as<arma::vec>(data_[String(sizec_(0))]);
        data_sizec2_ = as<arma::vec>(data_[String(sizec_(1))]);

        if (historical) {
          data_sizec1_ = as<arma::vec>(data_[String(sizec_(2))]);
        }
      }
      
      arma::uvec data_repst3_;
      arma::uvec data_repst2_;
      arma::uvec data_repst1_;
      
      if (repst_used) {
        data_repst3_ = as<arma::uvec>(data_[String(repst_(0))]);
        data_repst2_ = as<arma::uvec>(data_[String(repst_(1))]);

        if (historical) {
          data_repst1_ = as<arma::uvec>(data_[String(repst_(2))]);
        }
      }
      
      arma::uvec data_matst3_;
      arma::uvec data_matst2_;
      arma::uvec data_matst1_;
      
      if (matst_used) {
        data_matst3_ = as<arma::uvec>(data_[String(matst_(0))]);
        data_matst2_ = as<arma::uvec>(data_[String(matst_(1))]);

        if (historical) {
          data_matst1_ = as<arma::uvec>(data_[String(matst_(2))]);
        }
      }
      
      arma::vec msf_size_min = as<arma::vec>(melchett_stageframe_["sizebin_min"]);
      arma::vec msf_size_max = as<arma::vec>(melchett_stageframe_["sizebin_max"]);
      arma::vec msf_sizeb_min = as<arma::vec>(melchett_stageframe_["sizebinb_min"]);
      arma::vec msf_sizeb_max = as<arma::vec>(melchett_stageframe_["sizebinb_max"]);
      arma::vec msf_sizec_min = as<arma::vec>(melchett_stageframe_["sizebinc_min"]);
      arma::vec msf_sizec_max = as<arma::vec>(melchett_stageframe_["sizebinc_max"]);
      arma::uvec msf_matstatus = as<arma::uvec>(melchett_stageframe_["matstatus"]);
      arma::uvec msf_obsstatus = as<arma::uvec>(melchett_stageframe_["obsstatus"]);
      arma::uvec msf_repstatus = as<arma::uvec>(melchett_stageframe_["repstatus"]);
      arma::uvec msf_indataset = as<arma::uvec>(melchett_stageframe_["indataset"]);
      arma::uvec msf_alive = as<arma::uvec>(melchett_stageframe_["alive"]);
      arma::uvec ind_stages = find(msf_indataset == 1);
      
      // Assign stages across dataset
      for (int i = 0; i < data_points; i++) {
        // Stage 1
        if (historical) {
          
          arma::uvec al_stages1a = find(msf_alive == data_alive1_(i));
          
          if (NumericVector::is_na(data_size1_(i))) data_size1_(i) = 0.0;
          arma::uvec lo_stages1 = find(msf_size_min < data_size1_(i));
          arma::uvec hi_stages1 = find(msf_size_max >= data_size1_(i));
          arma::uvec mainstages1 = intersect(lo_stages1, hi_stages1);
          mainstages1 = intersect(mainstages1, al_stages1a);
          
          if (sizeb_used) {
            arma::uvec lo_stages1b = find(msf_sizeb_min < data_sizeb1_(i));
            arma::uvec hi_stages1b = find(msf_sizeb_max >= data_sizeb1_(i));
            arma::uvec mainstages1b = intersect(lo_stages1b, hi_stages1b);
            
            mainstages1 = intersect(mainstages1, mainstages1b);
          }
          
          if (sizec_used) {
            arma::uvec lo_stages1c = find(msf_sizec_min < data_sizec1_(i));
            arma::uvec hi_stages1c = find(msf_sizec_max >= data_sizec1_(i));
            arma::uvec mainstages1c = intersect(lo_stages1c, hi_stages1c);
            
            mainstages1 = intersect(mainstages1, mainstages1c);
          }
          
          if (matst_used) {
            arma::uvec mat_stages1m = find(msf_matstatus == data_matst1_(i));
            mainstages1 = intersect(mainstages1, mat_stages1m);
          }
          
          if (obsst_used) {
            arma::uvec obs_stages1m = find(msf_obsstatus == data_obsst1_(i));
            mainstages1 = intersect(mainstages1, obs_stages1m);
          }
          
          if (repst_used && !stage_NRasRep) {
            arma::uvec rep_stages1m = find(msf_repstatus == data_repst1_(i));
            mainstages1 = intersect(mainstages1, rep_stages1m);
          }
          
          mainstages1 = intersect(mainstages1, ind_stages);
          
          int no_mainstages1 = static_cast<int>(mainstages1.n_elem);
          if (no_mainstages1 > 0) {
            new_stageid1(i) = msf_stageid(static_cast<int>(mainstages1(0)));
            new_stage1(i) = msf_stage(static_cast<int>(mainstages1(0)));
          } else {
            new_stageid1(i) = msf_stageid(melchett_stageframe_length - 1);
            new_stage1(i) = msf_stage(melchett_stageframe_length - 1);
          }
        }
        
        // Stage 2
        arma::uvec al_stages2a = find(msf_alive == data_alive2_(i));
        
        if (NumericVector::is_na(data_size2_(i))) data_size2_(i) = 0.0;
        arma::uvec lo_stages2 = find(msf_size_min < data_size2_(i));
        arma::uvec hi_stages2 = find(msf_size_max >= data_size2_(i));
        arma::uvec mainstages2 = intersect(lo_stages2, hi_stages2);
        mainstages2 = intersect(mainstages2, al_stages2a);
        
        if (sizeb_used) {
          arma::uvec lo_stages2b = find(msf_sizeb_min < data_sizeb2_(i));
          arma::uvec hi_stages2b = find(msf_sizeb_max >= data_sizeb2_(i));
          arma::uvec mainstages2b = intersect(lo_stages2b, hi_stages2b);
          
          mainstages2 = intersect(mainstages2, mainstages2b);
        }
        
        if (sizec_used) {
          arma::uvec lo_stages2c = find(msf_sizec_min < data_sizec2_(i));
          arma::uvec hi_stages2c = find(msf_sizec_max >= data_sizec2_(i));
          arma::uvec mainstages2c = intersect(lo_stages2c, hi_stages2c);
          
          mainstages2 = intersect(mainstages2, mainstages2c);
        }
        
        if (matst_used) {
          arma::uvec mat_stages2m = find(msf_matstatus == data_matst2_(i));
          mainstages2 = intersect(mainstages2, mat_stages2m);
        }
        
        if (obsst_used) {
          arma::uvec obs_stages2m = find(msf_obsstatus == data_obsst2_(i));
          mainstages2 = intersect(mainstages2, obs_stages2m);
        }
        
        if (repst_used && !stage_NRasRep) {
          arma::uvec rep_stages2m = find(msf_repstatus == data_repst2_(i));
          mainstages2 = intersect(mainstages2, rep_stages2m);
        }
        
        mainstages2 = intersect(mainstages2, ind_stages);
        
        int no_mainstages2 = static_cast<int>(mainstages2.n_elem);
        if (no_mainstages2 > 0) {
          new_stageid2(i) = msf_stageid(static_cast<int>(mainstages2(0)));
          new_stage2(i) = msf_stage(static_cast<int>(mainstages2(0)));
        } else {
          new_stageid2(i) = msf_stageid(melchett_stageframe_length - 1);
          new_stage2(i) = msf_stage(melchett_stageframe_length - 1);
        }
        
        // Stage 3
        arma::uvec al_stages3a = find(msf_alive == data_alive3_(i));
        
        if (NumericVector::is_na(data_size3_(i))) data_size3_(i) = 0.0;
        arma::uvec lo_stages3 = find(msf_size_min < data_size3_(i));
        arma::uvec hi_stages3 = find(msf_size_max >= data_size3_(i));
        arma::uvec mainstages3 = intersect(lo_stages3, hi_stages3);
        mainstages3 = intersect(mainstages3, al_stages3a);
        
        if (sizeb_used) {
          arma::uvec lo_stages3b = find(msf_sizeb_min < data_sizeb3_(i));
          arma::uvec hi_stages3b = find(msf_sizeb_max >= data_sizeb3_(i));
          arma::uvec mainstages3b = intersect(lo_stages3b, hi_stages3b);
          
          mainstages3 = intersect(mainstages3, mainstages3b);
        }
        
        if (sizec_used) {
          arma::uvec lo_stages3c = find(msf_sizec_min < data_sizec3_(i));
          arma::uvec hi_stages3c = find(msf_sizec_max >= data_sizec3_(i));
          arma::uvec mainstages3c = intersect(lo_stages3c, hi_stages3c);
          
          mainstages3 = intersect(mainstages3, mainstages3c);
        }
        
        if (matst_used) {
          arma::uvec mat_stages3m = find(msf_matstatus == data_matst3_(i));
          mainstages3 = intersect(mainstages3, mat_stages3m);
        }
        
        if (obsst_used) {
          arma::uvec obs_stages3m = find(msf_obsstatus == data_obsst3_(i));
          mainstages3 = intersect(mainstages3, obs_stages3m);
        }
        
        if (repst_used && !stage_NRasRep) {
          arma::uvec rep_stages3m = find(msf_repstatus == data_repst3_(i));
          mainstages3 = intersect(mainstages3, rep_stages3m);
        }
        
        mainstages3 = intersect(mainstages3, ind_stages);
        
        int no_mainstages3 = static_cast<int>(mainstages3.n_elem);
        if (no_mainstages3 > 0) {
          new_stageid3(i) = msf_stageid(static_cast<int>(mainstages3(0)));
          new_stage3(i) = msf_stage(static_cast<int>(mainstages3(0)));
        } else {
          new_stageid3(i) = msf_stageid(melchett_stageframe_length - 1);
          new_stage3(i) = msf_stage(melchett_stageframe_length - 1);
        }
      }
      
      // Put data frame back together
      StringVector d_rn = data_.attr("row.names");
      
      data_["usedstage3"] = new_stage3;
      data_["usedstage2"] = new_stage2;
      if (historical) data_["usedstage1"] = new_stage1;
      
      data_["index3"] = new_stageid3 - 1;
      data_["index2"] = new_stageid2 - 1;
      if (historical) data_["index1"] = new_stageid1 - 1;
      
      bool new_s3i_needed {false};
      int s3i_loc {-1};
      StringVector data_var_allnames = data_.attr("names");
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_hard(String(data_var_allnames(i)), "stage3index")) {
          s3i_loc = i;
        }
      }
      if (s3i_loc != -1) {
        IntegerVector s3i_check = as<IntegerVector>(data_[s3i_loc]);
        IntegerVector s3i_unique = unique(s3i_check);
        
        if (s3i_unique.length() < 2) new_s3i_needed = true;
      } else {
        new_s3i_needed = true;
      }
      if (new_s3i_needed) {
        data_["stage3index"] = new_stageid3;
        data_["stage2index"] = new_stageid2;
        if (historical) data_["stage1index"] = new_stageid1;
      }
      
      data_.attr("class") = "data.frame";
      data_.attr("row.names") = d_rn;
      
    } else {
      int ind_matches {0};
      for (int i = 0; i < data_vars_no; i++) {
        if (stringcompare_simple(String(data_vars(i)), "index"), false) {
          if (is<IntegerVector>(data_[i])) ind_matches++;
        }
      }
      if (ind_matches < 2 || ind_matches > 3) new_stage_indices_needed = true;
      
      if (new_stage_indices_needed) {
        IntegerVector new_s3index (data_points);
        IntegerVector new_s2index (data_points);
        IntegerVector new_s1index (data_points);
        
        StringVector stages3_trial = as<StringVector>(data_[stages_int(0)]);
        StringVector stages2_trial = as<StringVector>(data_[stages_int(1)]);
        
        StringVector stages1_trial;
        if (historical) stages1_trial = as<StringVector>(data_[stages_int(2)]);
        
        for (int i = 0; i < data_points; i++) {
          for (int j = 0; j < melchett_stageframe_length; j++) {
            if (stringcompare_hard(String(stages3_trial(i)), String(msf_stage(j)))) new_s3index(i) = (j + 1);
            if (stringcompare_hard(String(stages2_trial(i)), String(msf_stage(j)))) new_s2index(i) = (j + 1);
            
            if (j == (melchett_stageframe_length - 1) && new_s3index(i) == 0) {
              new_s3index(i) = melchett_stageframe_length;
            }
            
            if (j == (melchett_stageframe_length - 1) && new_s2index(i) == 0) {
              new_s2index(i) = melchett_stageframe_length;
            }
            
            if (historical) {
              if (stringcompare_hard(String(stages1_trial(i)), String(msf_stage(j)))) new_s1index(i) = (j + 1);
              
              if (j == (melchett_stageframe_length - 1) && new_s1index(i) == 0) {
                new_s1index(i) = melchett_stageframe_length;
              }
            }
          }
        }
        
        StringVector d_rn = data_.attr("row.names");
        
        data_["stage3index"] = new_s3index;
        data_["stage2index"] = new_s2index;
        if (historical) data_["stage1index"] = new_s1index;
        
        data_["index3"] = new_s3index - 1;
        data_["index2"] = new_s2index - 1;
        if (historical) data_["index1"] = new_s1index - 1;
        
        data_.attr("class") = "data.frame";
        data_.attr("row.names") = d_rn;
      }
    }
  }
  
  // Individual covariate vectors for function-based MPMs
  // Individual covariate a
  StringVector inda_names; // All indcova categories, if factor
  NumericVector f1_inda_num; // Numeric values entered as input
  NumericVector f2_inda_num;
  StringVector f1_inda_cat; // Categorical (factor) values entered as input - fixed
  StringVector f2_inda_cat;
  StringVector r1_inda; // Categorical (factor) values entered as input - random
  StringVector r2_inda;
  
  if (inda.isNotNull() && !raw) {
    RObject inda_input = as<RObject>(inda);
    
    NumericVector inda_num;
    StringVector inda_cat;
    
    int no_mainyears = mainyears_.length();
    
    if (!paramnames_provided) {
      throw Rcpp::exception("Use of individual covariates requires a valid paramnames object", false);
    }
    
    if (!random_inda) {
      // Fixed covariate
      if (is<StringVector>(inda_input)) {
        inda_cat = as<StringVector>(inda_input);
        int inda_cat_length = static_cast<int>(inda_cat.length());
        
        if (inda_cat_length != 1 && inda_cat_length != 2) {
          if (inda_cat_length != no_mainyears) {
            String eat_my_shorts = "Individual covariate vector a must be empty, or include ";
            eat_my_shorts += "1, 2, or as many elements as occasions in the dataset.";
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
        
        if (!nodata) {
          StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
          StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
          
          int no_params = static_cast<int>(mainparams.length());
          
          int indacol_pm {-1};
          int indacol {-1};
          for (int i = 0; i < no_params; i++) {
            if (stringcompare_hard(String(mainparams(i)), "indcova2")) {
              if (!stringcompare_hard(String(modelparams(i)), "none")) {
                indacol_pm = i;
              }
            }
          }
          
          if (indacol_pm != -1) {
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), String(modelparams(indacol_pm)))) {
                indacol = i;
              }
            }
          }
          
          if (indacol == -1) pop_error("Individual covariate a", "paramnames", "modelsuite", 13);
          
          RObject test_inda(data_[indacol]);
          if (is<IntegerVector>(test_inda)) {
            
            IntegerVector iac_int = as<IntegerVector>(test_inda);
            if (iac_int.hasAttribute("levels")) {
              inda_names = as<StringVector>(iac_int.attr("levels"));
            } else {
              StringVector inda_values = as<StringVector>(test_inda);
              inda_names = sort_unique(inda_values);
            }
          } else {
            StringVector inda_values = as<StringVector>(data_[indacol]);
            inda_names = sort_unique(inda_values);
          }
        } else {
          if (!modelsuite_vrm || !modelsuite_provided) {
            throw Rcpp::exception("Individual covariate modeling requires a valid modelsuite.", false);
          }
          
          StringVector ms_names = modelsuite_.attr("names");
          int ms_length = ms_names.length();
          
          int indcova_frame_elem {-1};
          for (int i = 0; i < ms_length; i++) {
            if (stringcompare_hard(String(ms_names(i)), "indcova2_frame")) indcova_frame_elem = i;
          }
          
          if (indcova_frame_elem == -1) {
            String eat_my_shorts = "This function cannot use inda input with a vrm_input ";
            String eat_my_shorts1 = "object that does not include an indcova_frame element.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          
          DataFrame indca2 = as<DataFrame>(modelsuite_["indcova2_frame"]);
          inda_names = as<StringVector>(indca2["indcova"]);
        }
        
        int inda_names_length = static_cast<int>(inda_names.length());
        
        for (int i = 0; i < inda_cat_length; i++) {
          int inda_matches {0};
          
          for (int j = 0; j < inda_names_length; j++) {
            if (stringcompare_hard(String(inda_cat(i)), String(inda_names(j)))) inda_matches++;
          }
          if (inda_matches == 0) pop_error("inda", "dataset", "", 17);
        }
        
        if (inda_cat_length == 1) {
          StringVector sub_f1(no_mainyears);
          StringVector sub_f2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_f1(i) = inda_cat(0);
            sub_f2(i) = inda_cat(0);
          }
          
          f1_inda_cat = sub_f1;
          f2_inda_cat = sub_f2;
          
        } else if (inda_cat_length == 2 && no_years != 2) {
          StringVector sub_f1(no_mainyears);
          StringVector sub_f2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_f1(i) = inda_cat(0);
            sub_f2(i) = inda_cat(1);
          }
          
          f1_inda_cat = sub_f1;
          f2_inda_cat = sub_f2;
          
        } else {
          StringVector sub_f1(inda_cat_length);
          sub_f1(0) = "none";
          
          for (int i = 0; i < (inda_cat_length - 1); i++) {
            sub_f1(i + 1) = inda_cat(i);
          }
          f1_inda_cat = sub_f1;
          
          f2_inda_cat = inda_cat;
        }
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = "none";
          sub_r2(i) = "none";
        }
        
        r1_inda = sub_r1;
        r2_inda = sub_r2;
        
        f2_inda_num = rep(0, no_mainyears);
        f1_inda_num = rep(0, no_mainyears);
        
      } else if (is<IntegerVector>(inda_input) || is<NumericVector>(inda_input)) {
        // Handles possibility of factor variables
        inda_num = as<NumericVector>(inda_input);
        int inda_num_length = static_cast<int>(inda_num.length());
        bool factor_variable {false};
        
        if (inda_num_length != 1 && inda_num_length != 2) {
          if (inda_num_length != no_mainyears) {
            String eat_my_shorts = "Individual covariate vector a must be empty, or include ";
            String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
        
        if (!nodata) {
          StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
          StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
          
          int no_params = static_cast<int>(mainparams.length());
          
          int indacol_pm {-1};
          int indacol {-1};
          for (int i = 0; i < no_params; i++) {
            if (stringcompare_hard(String(mainparams(i)), "indcova2")) {
              if (!stringcompare_hard(String(modelparams(i)), "none")) {
                indacol_pm = i;
              }
            }
          }
          
          if (indacol_pm != -1) {
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), String(modelparams(indacol_pm)))) {
                indacol = i;
              }
            }
          }
          
          if (indacol == -1) pop_error("Individual covariate a", "paramnames", "modelsuite", 13);
          
          RObject test_inda(data_[indacol]);
          if (is<IntegerVector>(test_inda)) {
            
            IntegerVector iac_int = as<IntegerVector>(test_inda);
            if (iac_int.hasAttribute("levels")) {
              inda_names = as<StringVector>(iac_int.attr("levels"));
              inda_cat = as<StringVector>(inda_input);
              factor_variable = true;
              
            } else {
              StringVector inda_values = as<StringVector>(test_inda);
              inda_names = sort_unique(inda_values);
            }
          }
        } // else { /* need vrm input code*/ }
        
        if (!factor_variable) {
          inda_names = {"0"};
          
          if (inda_num_length == 1) {
            f1_inda_num = rep(inda_num(0), no_mainyears); 
            f2_inda_num = rep(inda_num(0), no_mainyears);
            
          } else if (inda_num_length == 2 && no_years != 2) {
            f1_inda_num = rep(inda_num(0), no_mainyears);
            f2_inda_num = rep(inda_num(1), no_mainyears);
            
          } else {
            NumericVector sub_f1(inda_num_length);
            sub_f1(0) = 0;
            
            for (int i = 0; i < (inda_num_length - 1); i++) {
              sub_f1(i + 1) = inda_num(i);
            }
            f1_inda_num = sub_f1;
            
            f2_inda_num = inda_num;
          }
          StringVector sub_r1(no_mainyears);
          StringVector sub_r2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_r1(i) = "none";
            sub_r2(i) = "none";
          }
          
          r1_inda = sub_r1;
          r2_inda = sub_r2;
          
          f1_inda_cat = clone(sub_r1);
          f2_inda_cat = clone(sub_r2);
          
        } else { // If indcova is a factor variable
          int inda_cat_length = inda_cat.length();
          
          if (inda_cat_length == 1) {
            StringVector sub_f1(no_mainyears);
            StringVector sub_f2(no_mainyears);
            
            for (int i = 0; i < no_mainyears; i++) {
              sub_f1(i) = inda_cat(0);
              sub_f2(i) = inda_cat(0);
            }
            
            f1_inda_cat = sub_f1;
            f2_inda_cat = sub_f2;
            
          } else if (inda_cat_length == 2 && no_years != 2) {
            StringVector sub_f1(no_mainyears);
            StringVector sub_f2(no_mainyears);
            
            for (int i = 0; i < no_mainyears; i++) {
              sub_f1(i) = inda_cat(0);
              sub_f2(i) = inda_cat(1);
            }
            
            f1_inda_cat = sub_f1;
            f2_inda_cat = sub_f2;
            
          } else {
            StringVector sub_f1(inda_cat_length);
            sub_f1(0) = "none";
            
            for (int i = 0; i < (inda_cat_length - 1); i++) {
              sub_f1(i + 1) = inda_cat(i);
            }
            f1_inda_cat = sub_f1;
            
            f2_inda_cat = inda_cat;
          }
          StringVector sub_r1(no_mainyears);
          StringVector sub_r2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_r1(i) = "none";
            sub_r2(i) = "none";
          }
          
          r1_inda = sub_r1;
          r2_inda = sub_r2;
          
          f2_inda_num = rep(0, no_mainyears);
          f1_inda_num = rep(0, no_mainyears);
        }
      } else {
        throw Rcpp::exception("Format of indcova not recognized.", false);
      }
      
    } else {
      // Random covariate
      inda_cat = as<StringVector>(inda_input);
      int inda_cat_length = static_cast<int>(inda_cat.length());
      
      if (inda_cat_length != 1 && inda_cat_length != 2) {
        if (inda_cat_length != no_mainyears) {
          String eat_my_shorts = "Individual covariate vector a must be empty, or include ";
          String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      
      if (!nodata) {
        StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
        StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
        
        int no_params = static_cast<int>(mainparams.length());
        
        int indacol_pm {-1};
        int indacol {-1};
        for (int i = 0; i < no_params; i++) {
          if (stringcompare_hard(String(mainparams(i)), "indcova2")) {
            if (!stringcompare_hard(String(modelparams(i)), "none")) {
              indacol_pm = i;
            }
          }
        }
        
        if (indacol_pm != -1) {
          for (int i = 0; i < data_vars_no; i++) {
            if (stringcompare_hard(String(data_vars(i)), String(modelparams(indacol_pm)))) {
              indacol = i;
            }
          }
        }
        
        if (indacol == -1) pop_error("Individual covariate a", "paramnames", "modelsuite", 13);
        
        RObject test_inda(data_[indacol]);
        if (is<IntegerVector>(test_inda)) {
          
          IntegerVector iac_int = as<IntegerVector>(test_inda);
          if (iac_int.hasAttribute("levels")) {
            inda_names = as<StringVector>(iac_int.attr("levels"));
          } else {
            StringVector inda_values = as<StringVector>(test_inda);
            inda_names = sort_unique(inda_values);
          }
        } else {
          StringVector inda_values = as<StringVector>(data_[indacol]);
          inda_names = sort_unique(inda_values);
        }
      } else {
        if (!modelsuite_vrm || !modelsuite_provided) {
          throw Rcpp::exception("Individual covariate modeling requires a valid modelsuite.", false);
        }
        
        StringVector ms_names = modelsuite_.attr("names");
        int ms_length = ms_names.length();
        
        int indcova_frame_elem {-1};
        for (int i = 0; i < ms_length; i++) {
          if (stringcompare_hard(String(ms_names(i)), "indcova2_frame")) indcova_frame_elem = i;
        }
        
        if (indcova_frame_elem == -1) {
          String eat_my_shorts = "This function cannot use inda input with a vrm_input object ";
          String eat_my_shorts1 = "that does not include an indcova_frame element.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
        DataFrame indca2 = as<DataFrame>(modelsuite_["indcova2_frame"]);
        inda_names = as<StringVector>(indca2["indcova"]);
      }
      
      int inda_names_length = static_cast<int>(inda_names.length());
      
      for (int i = 0; i < inda_cat_length; i++) {
        int inda_matches {0};
        
        for (int j = 0; j < inda_names_length; j++) {
          if (stringcompare_hard(String(inda_cat(i)), String(inda_names(j)))) inda_matches++;
        }
        if (inda_matches == 0) pop_error("inda", "dataset", "", 17);
      }
      
      if (inda_cat_length == 1) {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = inda_cat(0);
          sub_r2(i) = inda_cat(0);
        }
        
        r1_inda = sub_r1;
        r2_inda = sub_r2;
        
      } else if (inda_cat_length == 2 && no_years != 2) {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = inda_cat(0);
          sub_r2(i) = inda_cat(1);
        }
        
        r1_inda = sub_r1;
        r2_inda = sub_r2;
        
      } else {
        StringVector sub_r1(inda_cat_length);
        sub_r1(0) = "none";
        
        for (int i = 0; i < (inda_cat_length - 1); i++) {
          sub_r1(i + 1) = inda_cat(i);
        }
        r1_inda = sub_r1;
        
        r2_inda = inda_cat;
      }
      f1_inda_num = rep(0, no_mainyears);
      f2_inda_num = rep(0, no_mainyears);
      
      {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = "none";
          sub_r2(i) = "none";
        }
        
        f1_inda_cat = sub_r1;
        f2_inda_cat = sub_r2;
      }
    }
  } else {
    int no_mainyears = mainyears_.length();
    inda_names = {"0"};
    
    f1_inda_num = rep(0, no_mainyears);
    f2_inda_num = rep(0, no_mainyears);
    
    StringVector sub_r1(no_mainyears);
    StringVector sub_r2(no_mainyears);
    
    for (int i = 0; i < no_mainyears; i++) {
      sub_r1(i) = "none";
      sub_r2(i) = "none";
    }
    
    r1_inda = sub_r1;
    r2_inda = sub_r2;
    
    f1_inda_cat = clone(sub_r1);
    f2_inda_cat = clone(sub_r2);
  }
  
  // Individual covariate b
  StringVector indb_names;
  NumericVector f1_indb_num;
  NumericVector f2_indb_num;
  StringVector f1_indb_cat;
  StringVector f2_indb_cat;
  StringVector r1_indb;
  StringVector r2_indb;
  
  if (indb.isNotNull() && !raw) {
    RObject indb_input = as<RObject>(indb);
    
    NumericVector indb_num;
    StringVector indb_cat;
    
    int no_mainyears = mainyears_.length();
    
    if (!paramnames_provided) {
      throw Rcpp::exception("Use of individual covariates requires a valid paramnames object",
        false);
    }
    
    if (!random_indb) {
      // Fixed covariate
      if (is<StringVector>(indb_input)) { 
        indb_cat = as<StringVector>(indb_input);
        int indb_cat_length = static_cast<int>(indb_cat.length());
        
        if (indb_cat_length != 1 && indb_cat_length != 2) {
          if (indb_cat_length != no_mainyears) {
            String eat_my_shorts = "Individual covariate vector b must be empty, or include ";
            String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
        
        if (!nodata) {
          StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
          StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
          
          int no_params = static_cast<int>(mainparams.length());
          
          int indbcol_pm {-1};
          int indbcol {-1};
          for (int i = 0; i < no_params; i++) {
            if (stringcompare_hard(String(mainparams(i)), "indcovb2")) {
              if (!stringcompare_hard(String(modelparams(i)), "none")) {
                indbcol_pm = i;
              }
            }
          }
          
          if (indbcol_pm != -1) {
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), String(modelparams(indbcol_pm)))) {
                indbcol = i;
              }
            }
          }
          
          if (indbcol == -1) pop_error("Individual covariate b", "paramnames", "modelsuite", 13);
          
          RObject test_indb(data_[indbcol]);
          if (is<IntegerVector>(test_indb)) {
            
            IntegerVector ibc_int = as<IntegerVector>(test_indb);
            if (ibc_int.hasAttribute("levels")) {
              indb_names = as<StringVector>(ibc_int.attr("levels"));
            } else {
              StringVector indb_values = as<StringVector>(test_indb);
              indb_names = sort_unique(indb_values);
            }
          } else {
            StringVector indb_values = as<StringVector>(data_[indbcol]);
            indb_names = sort_unique(indb_values);
          }
        } else {
          if (!modelsuite_vrm || !modelsuite_provided) {
            throw Rcpp::exception("Individual covariate modeling requires a valid modelsuite.",
              false);
          }
          
          StringVector ms_names = modelsuite_.attr("names");
          int ms_length = ms_names.length();
          
          int indcovb_frame_elem {-1};
          for (int i = 0; i < ms_length; i++) {
            if (stringcompare_hard(String(ms_names(i)), "indcovb2_frame")) indcovb_frame_elem = i;
          }
          
          if (indcovb_frame_elem == -1) {
            String eat_my_shorts = "This function cannot use indb input with a vrm_input ";
            String eat_my_shorts1 = "object that does not include an indcovb_frame element.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          
          DataFrame indcb2 = as<DataFrame>(modelsuite_["indcovb2_frame"]);
          indb_names = as<StringVector>(indcb2["indcovb"]);
        }
        
        int indb_names_length = static_cast<int>(indb_names.length());
        
        for (int i = 0; i < indb_cat_length; i++) {
          int indb_matches {0};
          
          for (int j = 0; j < indb_names_length; j++) {
            if (stringcompare_hard(String(indb_cat(i)), String(indb_names(j)))) indb_matches++;
          }
          if (indb_matches == 0) pop_error("indb", "dataset", "", 17);
        }
        
        if (indb_cat_length == 1) {
          StringVector sub_f1(no_mainyears);
          StringVector sub_f2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_f1(i) = indb_cat(0);
            sub_f2(i) = indb_cat(0);
          }
          
          f1_indb_cat = sub_f1;
          f2_indb_cat = sub_f2;
          
        } else if (indb_cat_length == 2 && no_years != 2) {
          StringVector sub_f1(no_mainyears);
          StringVector sub_f2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_f1(i) = indb_cat(0);
            sub_f2(i) = indb_cat(1);
          }
          
          f1_indb_cat = sub_f1;
          f2_indb_cat = sub_f2;
          
        } else {
          StringVector sub_f1(indb_cat_length);
          sub_f1(0) = "none";
          
          for (int i = 0; i < (indb_cat_length - 1); i++) {
            sub_f1(i + 1) = indb_cat(i);
          }
          f1_indb_cat = sub_f1;
          
          f2_indb_cat = indb_cat;
        }
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = "none";
          sub_r2(i) = "none";
        }
        
        r1_indb = sub_r1;
        r2_indb = sub_r2;
        
      } else if (is<NumericVector>(indb_input) || is <IntegerVector>(indb_input)) {
        // Handles possibility of factor variables
        indb_num = as<NumericVector>(indb_input);
        int indb_num_length = static_cast<int>(indb_num.length());
        bool factor_variable {false};
        
        if (indb_num_length != 1 && indb_num_length != 2) {
          if (indb_num_length != no_mainyears) {
            String eat_my_shorts = "Individual covariate vector b must be empty, or include ";
            String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
        
        if (!nodata) {
          StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
          StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
          
          int no_params = static_cast<int>(mainparams.length());
          
          int indbcol_pm {-1};
          int indbcol {-1};
          for (int i = 0; i < no_params; i++) {
            if (stringcompare_hard(String(mainparams(i)), "indcovb2")) {
              if (!stringcompare_hard(String(modelparams(i)), "none")) {
                indbcol_pm = i;
              }
            }
          }
          
          if (indbcol_pm != -1) {
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), String(modelparams(indbcol_pm)))) {
                indbcol = i;
              }
            }
          }
          
          if (indbcol == -1) pop_error("Individual covariate b", "paramnames", "modelsuite", 13);
          
          RObject test_indb(data_[indbcol]);
          if (is<IntegerVector>(test_indb)) {
            
            IntegerVector ibc_int = as<IntegerVector>(test_indb);
            if (ibc_int.hasAttribute("levels")) {
              indb_names = as<StringVector>(ibc_int.attr("levels"));
              indb_cat = as<StringVector>(indb_input);
              factor_variable = true;
              
            } else {
              StringVector indb_values = as<StringVector>(test_indb);
              indb_names = sort_unique(indb_values);
            }
          }
        } // else { /* need vrm input code*/ }
        
        if (!factor_variable) {
          indb_names = {"0"};
          
          if (indb_num_length == 1) {
            f1_indb_num = rep(indb_num(0), no_mainyears);
            f2_indb_num = rep(indb_num(0), no_mainyears);
            
          } else if (indb_num_length == 2 && no_years != 2) {
            f1_indb_num = rep(indb_num(0), no_mainyears);
            f2_indb_num = rep(indb_num(1), no_mainyears);
            
          } else {
            NumericVector sub_f1(indb_num_length);
            sub_f1(0) = 0;
            
            for (int i = 0; i < (indb_num_length - 1); i++) {
              sub_f1(i + 1) = indb_num(i);
            }
            f1_indb_num = sub_f1;
            
            f2_indb_num = indb_num;
          }
          StringVector sub_r1(no_mainyears);
          StringVector sub_r2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_r1(i) = "none";
            sub_r2(i) = "none";
          }
          
          r1_indb = sub_r1;
          r2_indb = sub_r2;
          
          f1_indb_cat = clone(sub_r1);
          f2_indb_cat = clone(sub_r2);
          
        } else { // If indcovb is a factor variable
          int indb_cat_length = indb_cat.length();
          
          if (indb_cat_length == 1) {
            StringVector sub_f1(no_mainyears);
            StringVector sub_f2(no_mainyears);
            
            for (int i = 0; i < no_mainyears; i++) {
              sub_f1(i) = indb_cat(0);
              sub_f2(i) = indb_cat(0);
            }
            
            f1_indb_cat = sub_f1;
            f2_indb_cat = sub_f2;
            
          } else if (indb_cat_length == 2 && no_years != 2) {
            StringVector sub_f1(no_mainyears);
            StringVector sub_f2(no_mainyears);
            
            for (int i = 0; i < no_mainyears; i++) {
              sub_f1(i) = indb_cat(0);
              sub_f2(i) = indb_cat(1);
            }
            
            f1_indb_cat = sub_f1;
            f2_indb_cat = sub_f2;
            
          } else {
            StringVector sub_f1(indb_cat_length);
            sub_f1(0) = "none";
            
            for (int i = 0; i < (indb_cat_length - 1); i++) {
              sub_f1(i + 1) = indb_cat(i);
            }
            f1_indb_cat = sub_f1;
            
            f2_indb_cat = indb_cat;
          }
          StringVector sub_r1(no_mainyears);
          StringVector sub_r2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_r1(i) = "none";
            sub_r2(i) = "none";
          }
          
          r1_indb = sub_r1;
          r2_indb = sub_r2;
          
          f2_indb_num = rep(0, no_mainyears);
          f1_indb_num = rep(0, no_mainyears);
        }
      } else {
        throw Rcpp::exception("Format of indcovb not recognized.", false);
      }
      
    } else {
      // Random covariate
      indb_cat = as<StringVector>(indb_input);
      int indb_cat_length = static_cast<int>(indb_cat.length());
      
      if (indb_cat_length != 1 && indb_cat_length != 2) {
        if (indb_cat_length != no_mainyears) {
          String eat_my_shorts = "Individual covariate vector b must be empty, or include ";
          String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      
      if (!nodata) {
        StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
        StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
        
        int no_params = static_cast<int>(mainparams.length());
        
        int indbcol_pm {-1};
        int indbcol {-1};
        for (int i = 0; i < no_params; i++) {
          if (stringcompare_hard(String(mainparams(i)), "indcovb2")) {
            if (!stringcompare_hard(String(modelparams(i)), "none")) {
              indbcol_pm = i;
            }
          }
        }
        
        if (indbcol_pm != -1) {
          for (int i = 0; i < data_vars_no; i++) {
            if (stringcompare_hard(String(data_vars(i)), String(modelparams(indbcol_pm)))) {
              indbcol = i;
            }
          }
        }
        
        if (indbcol == -1) pop_error("Individual covariate b", "paramnames", "modelsuite", 13);
        
        RObject test_indb(data_[indbcol]);
        if (is<IntegerVector>(test_indb)) {
          
          IntegerVector ibc_int = as<IntegerVector>(test_indb);
          if (ibc_int.hasAttribute("levels")) {
            indb_names = as<StringVector>(ibc_int.attr("levels"));
          } else {
            StringVector indb_values = as<StringVector>(test_indb);
            indb_names = sort_unique(indb_values);
          }
        } else {
          StringVector indb_values = as<StringVector>(data_[indbcol]);
          indb_names = sort_unique(indb_values);
        }
        
      } else {
        if (!modelsuite_vrm || !modelsuite_provided) {
          throw Rcpp::exception("Individual covariate modeling requires a valid modelsuite.",
            false);
        }
        
        StringVector ms_names = modelsuite_.attr("names");
        int ms_length = ms_names.length();
        
        int indcovb_frame_elem {-1};
        for (int i = 0; i < ms_length; i++) {
          if (stringcompare_hard(String(ms_names(i)), "indcovb2_frame")) indcovb_frame_elem = i;
        }
        
        if (indcovb_frame_elem == -1) {
          String eat_my_shorts = "This function cannot use indb input with a vrm_input object ";
          String eat_my_shorts1 = "that does not include an indcovb_frame element.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
        DataFrame indcb2 = as<DataFrame>(modelsuite_["indcovb2_frame"]);
        indb_names = as<StringVector>(indcb2["indcovb"]);
      }
      
      int indb_names_length = static_cast<int>(indb_names.length());
      
      for (int i = 0; i < indb_cat_length; i++) {
        int indb_matches {0};
        
        for (int j = 0; j < indb_names_length; j++) {
          if (stringcompare_hard(String(indb_cat(i)), String(indb_names(j)))) indb_matches++;
        }
        if (indb_matches == 0) pop_error("indb", "dataset", "", 17);
      }
      
      if (indb_cat_length == 1) {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = indb_cat(0);
          sub_r2(i) = indb_cat(0);
        }
        
        r1_indb = sub_r1;
        r2_indb = sub_r2;
        
      } else if (indb_cat_length == 2 && no_years != 2) {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = indb_cat(0);
          sub_r2(i) = indb_cat(1);
        }
        
        r1_indb = sub_r1;
        r2_indb = sub_r2;
        
      } else {
        StringVector sub_r1(indb_cat_length);
        sub_r1(0) = "none";
        
        for (int i = 0; i < (indb_cat_length - 1); i++) {
          sub_r1(i + 1) = indb_cat(i);
        }
        r1_indb = sub_r1;
        
        r2_indb = indb_cat;
      }
      f1_indb_num = rep(0, no_mainyears);
      f2_indb_num = rep(0, no_mainyears);
      
      {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = "none";
          sub_r2(i) = "none";
        }
        
        f1_indb_cat = sub_r1;
        f2_indb_cat = sub_r2;
      }
    }
  } else {
    int no_mainyears = mainyears_.length();
    indb_names = {"0"};
    
    f1_indb_num = rep(0, no_mainyears);
    f2_indb_num = rep(0, no_mainyears);
    
    StringVector sub_r1(no_mainyears);
    StringVector sub_r2(no_mainyears);
    
    for (int i = 0; i < no_mainyears; i++) {
      sub_r1(i) = "none";
      sub_r2(i) = "none";
    }
    
    r1_indb = sub_r1;
    r2_indb = sub_r2;
    
    f1_indb_cat = clone(sub_r1);
    f2_indb_cat = clone(sub_r2);
  }
  
  // Individual covariate c
  StringVector indc_names;
  NumericVector f1_indc_num;
  NumericVector f2_indc_num;
  StringVector f1_indc_cat;
  StringVector f2_indc_cat;
  StringVector r1_indc;
  StringVector r2_indc;
  
  if (indc.isNotNull() && !raw) {
    RObject indc_input = as<RObject>(indc);
    
    NumericVector indc_num;
    StringVector indc_cat;
    
    int no_mainyears = mainyears_.length();
    
    if (!paramnames_provided) {
      throw Rcpp::exception("Use of individual covariates requires a valid paramnames object",
        false);
    }
    
    if (!random_indc) {
      // Fixed covariate
      if (is<StringVector>(indc_input)) { 
        indc_cat = as<StringVector>(indc_input);
        int indc_cat_length = static_cast<int>(indc_cat.length());
        
        if (indc_cat_length != 1 && indc_cat_length != 2) {
          if (indc_cat_length != no_mainyears) {
            String eat_my_shorts = "Individual covariate vector c must be empty, or include 1, 2, ";
            String eat_my_shorts1 = "or as many elements as occasions in the dataset.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
        
        if (!nodata) {
          StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
          StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
          
          int no_params = static_cast<int>(mainparams.length());
          
          int indccol_pm {-1};
          int indccol {-1};
          for (int i = 0; i < no_params; i++) {
            if (stringcompare_hard(String(mainparams(i)), "indcovc2")) {
              if (!stringcompare_hard(String(modelparams(i)), "none")) {
                indccol_pm = i;
              }
            }
          }
          
          if (indccol_pm != -1) {
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), String(modelparams(indccol_pm)))) {
                indccol = i;
              }
            }
          }
          
          if (indccol == -1) pop_error("Individual covariate c", "paramnames", "modelsuite", 13);
          
          RObject test_indc(data_[indccol]);
          if (is<IntegerVector>(test_indc)) {
            
            IntegerVector icc_int = as<IntegerVector>(test_indc);
            if (icc_int.hasAttribute("levels")) {
              indc_names = as<StringVector>(icc_int.attr("levels"));
            } else {
              StringVector indc_values = as<StringVector>(test_indc);
              indc_names = sort_unique(indc_values);
            }
          } else {
            StringVector indc_values = as<StringVector>(data_[indccol]);
            indc_names = sort_unique(indc_values);
          }
          
        } else {
          if (!modelsuite_vrm || !modelsuite_provided) {
            throw Rcpp::exception("Individual covariate modeling requires a valid modelsuite.", false);
          }
          
          StringVector ms_names = modelsuite_.attr("names");
          int ms_length = ms_names.length();
          
          int indcovc_frame_elem {-1};
          for (int i = 0; i < ms_length; i++) {
            if (stringcompare_hard(String(ms_names(i)), "indcovc2_frame")) indcovc_frame_elem = i;
          }
          
          if (indcovc_frame_elem == -1) {
            String eat_my_shorts = "This function cannot use indc input with a vrm_input object ";
            String eat_my_shorts1 = "that does not include an indcovc_frame element.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          
          DataFrame indcc2 = as<DataFrame>(modelsuite_["indcovc2_frame"]);
          indc_names = as<StringVector>(indcc2["indcovc"]);
        }
        
        int indc_names_length = static_cast<int>(indc_names.length());
        
        for (int i = 0; i < indc_cat_length; i++) {
          int indc_matches {0};
          
          for (int j = 0; j < indc_names_length; j++) {
            if (stringcompare_hard(String(indc_cat(i)), String(indc_names(j)))) indc_matches++;
          }
          if (indc_matches == 0) pop_error("indc", "dataset", "", 17);
        }
        
        if (indc_cat_length == 1) {
          StringVector sub_f1(no_mainyears);
          StringVector sub_f2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_f1(i) = indc_cat(0);
            sub_f2(i) = indc_cat(0);
          }
          
          f1_indc_cat = sub_f1;
          f2_indc_cat = sub_f2;
          
        } else if (indc_cat_length == 2 && no_years != 2) {
          StringVector sub_f1(no_mainyears);
          StringVector sub_f2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_f1(i) = indc_cat(0);
            sub_f2(i) = indc_cat(1);
          }
          
          f1_indc_cat = sub_f1;
          f2_indc_cat = sub_f2;
          
        } else {
          StringVector sub_f1(indc_cat_length);
          sub_f1(0) = "none";
          
          for (int i = 0; i < (indc_cat_length - 1); i++) {
            sub_f1(i + 1) = indc_cat(i);
          }
          f1_indc_cat = sub_f1;
          
          f2_indc_cat = indc_cat;
        }
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = "none";
          sub_r2(i) = "none";
        }
        
        r1_indc = sub_r1;
        r2_indc = sub_r2;
        
      } else if (is<NumericVector>(indc_input) || is <IntegerVector>(indc_input)) {
       // Handles the possibility of factor variables
        indc_num = as<NumericVector>(indc_input);
        int indc_num_length = static_cast<int>(indc_num.length());
        bool factor_variable {false};
        
        
        if (indc_num_length != 1 && indc_num_length != 2) {
          if (indc_num_length != no_mainyears) {
            String eat_my_shorts = "Individual covariate vector c must be empty, or include ";
            String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
            eat_my_shorts += eat_my_shorts1;
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
        }
        
        if (!nodata) {
          StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
          StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
          
          int no_params = static_cast<int>(mainparams.length());
          
          int indccol_pm {-1};
          int indccol {-1};
          for (int i = 0; i < no_params; i++) {
            if (stringcompare_hard(String(mainparams(i)), "indcovc2")) {
              if (!stringcompare_hard(String(modelparams(i)), "none")) {
                indccol_pm = i;
              }
            }
          }
          
          if (indccol_pm != -1) {
            for (int i = 0; i < data_vars_no; i++) {
              if (stringcompare_hard(String(data_vars(i)), String(modelparams(indccol_pm)))) {
                indccol = i;
              }
            }
          }
          
          if (indccol == -1) pop_error("Individual covariate c", "paramnames", "modelsuite", 13);
          
          RObject test_indc(data_[indccol]);
          if (is<IntegerVector>(test_indc)) {
            
            IntegerVector icc_int = as<IntegerVector>(test_indc);
            if (icc_int.hasAttribute("levels")) {
              indc_names = as<StringVector>(icc_int.attr("levels"));
              indc_cat = as<StringVector>(indc_input);
              factor_variable = true;
              
            } else {
              StringVector indc_values = as<StringVector>(test_indc);
              indc_names = sort_unique(indc_values);
            }
          }
        } // else { /* need vrm input code*/ }
        
        if (!factor_variable) {
          indc_names = {"0"};
          
          if (indc_num_length == 1) {
            f1_indc_num = rep(indc_num(0), no_mainyears);
            f2_indc_num = rep(indc_num(0), no_mainyears);
            
          } else if (indc_num_length == 2 && no_years != 2) {
            f1_indc_num = rep(indc_num(0), no_mainyears);
            f2_indc_num = rep(indc_num(1), no_mainyears);
            
          } else {
            NumericVector sub_f1(indc_num_length);
            sub_f1(0) = 0;
            
            for (int i = 0; i < (indc_num_length - 1); i++) {
              sub_f1(i + 1) = indc_num(i);
            }
            f1_indc_num = sub_f1;
            
            f2_indc_num = indc_num;
          }
          StringVector sub_r1(no_mainyears);
          StringVector sub_r2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_r1(i) = "none";
            sub_r2(i) = "none";
          }
          
          r1_indc = sub_r1;
          r2_indc = sub_r2;
          
          f1_indc_cat = clone(sub_r1);
          f2_indc_cat = clone(sub_r2);
          
        } else { // If indcovc is a factor variable
          int indc_cat_length = indc_cat.length();
          
          if (indc_cat_length == 1) {
            StringVector sub_f1(no_mainyears);
            StringVector sub_f2(no_mainyears);
            
            for (int i = 0; i < no_mainyears; i++) {
              sub_f1(i) = indc_cat(0);
              sub_f2(i) = indc_cat(0);
            }
            
            f1_indc_cat = sub_f1;
            f2_indc_cat = sub_f2;
            
          } else if (indc_cat_length == 2 && no_years != 2) {
            StringVector sub_f1(no_mainyears);
            StringVector sub_f2(no_mainyears);
            
            for (int i = 0; i < no_mainyears; i++) {
              sub_f1(i) = indc_cat(0);
              sub_f2(i) = indc_cat(1);
            }
            
            f1_indc_cat = sub_f1;
            f2_indc_cat = sub_f2;
            
          } else {
            StringVector sub_f1(indc_cat_length);
            sub_f1(0) = "none";
            
            for (int i = 0; i < (indc_cat_length - 1); i++) {
              sub_f1(i + 1) = indc_cat(i);
            }
            f1_indc_cat = sub_f1;
            
            f2_indc_cat = indc_cat;
          }
          StringVector sub_r1(no_mainyears);
          StringVector sub_r2(no_mainyears);
          
          for (int i = 0; i < no_mainyears; i++) {
            sub_r1(i) = "none";
            sub_r2(i) = "none";
          }
          
          r1_indc = sub_r1;
          r2_indc = sub_r2;
          
          f2_indc_num = rep(0, no_mainyears);
          f1_indc_num = rep(0, no_mainyears);
        }
      } else {
        throw Rcpp::exception("Format of indcovc not recognized.", false);
      }
      
    } else {
      // Random covariate
      indc_cat = as<StringVector>(indc_input);
      int indc_cat_length = static_cast<int>(indc_cat.length());
      
      if (indc_cat_length != 1 && indc_cat_length != 2) {
        if (indc_cat_length != no_mainyears) {
          String eat_my_shorts = "Individual covariate vector c must be empty, or include ";
          String eat_my_shorts1 = "1, 2, or as many elements as occasions in the dataset.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      
      if (!nodata) {
        StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
        StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
        
        int no_params = static_cast<int>(mainparams.length());
        
        int indccol_pm {-1};
        int indccol {-1};
        for (int i = 0; i < no_params; i++) {
          if (stringcompare_hard(String(mainparams(i)), "indcovc2")) {
            if (!stringcompare_hard(String(modelparams(i)), "none")) {
              indccol_pm = i;
            }
          }
        }
        
        if (indccol_pm != -1) {
          for (int i = 0; i < data_vars_no; i++) {
            if (stringcompare_hard(String(data_vars(i)), String(modelparams(indccol_pm)))) {
              indccol = i;
            }
          }
        }
        
        if (indccol == -1) pop_error("Individual covariate c", "paramnames", "modelsuite", 13);
        
        RObject test_indc(data_[indccol]);
        if (is<IntegerVector>(test_indc)) {
          
          IntegerVector icc_int = as<IntegerVector>(test_indc);
          if (icc_int.hasAttribute("levels")) {
            indc_names = as<StringVector>(icc_int.attr("levels"));
          } else {
            StringVector indc_values = as<StringVector>(test_indc);
            indc_names = sort_unique(indc_values);
          }
        } else {
          StringVector indc_values = as<StringVector>(data_[indccol]);
          indc_names = sort_unique(indc_values);
        }
        
      } else {
        if (!modelsuite_vrm || !modelsuite_provided) {
          throw Rcpp::exception("Individual covariate modeling requires a valid modelsuite.", false);
        }
        
        StringVector ms_names = modelsuite_.attr("names");
        int ms_length = ms_names.length();
        
        int indcovc_frame_elem {-1};
        for (int i = 0; i < ms_length; i++) {
          if (stringcompare_hard(String(ms_names(i)), "indcovc2_frame")) indcovc_frame_elem = i;
        }
        
        if (indcovc_frame_elem == -1) {
          String eat_my_shorts = "This function cannot use indc input with a vrm_input object ";
          String eat_my_shorts1 = "that does not include an indcovc_frame element.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
        DataFrame indcc2 = as<DataFrame>(modelsuite_["indcovc2_frame"]);
        indc_names = as<StringVector>(indcc2["indcovc"]);
      }
      
      int indc_names_length = static_cast<int>(indc_names.length());
      
      for (int i = 0; i < indc_cat_length; i++) {
        int indc_matches {0};
        
        for (int j = 0; j < indc_names_length; j++) {
          if (stringcompare_hard(String(indc_cat(i)), String(indc_names(j)))) indc_matches++;
        }
        if (indc_matches == 0) pop_error("indc", "dataset", "", 17);
      }
      
      if (indc_cat_length == 1) {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = indc_cat(0);
          sub_r2(i) = indc_cat(0);
        }
        
        r1_indc = sub_r1;
        r2_indc = sub_r2;
        
      } else if (indc_cat_length == 2 && no_years != 2) {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = indc_cat(0);
          sub_r2(i) = indc_cat(1);
        }
        
        r1_indc = sub_r1;
        r2_indc = sub_r2;
        
      } else {
        StringVector sub_r1(indc_cat_length);
        sub_r1(0) = "none";
        
        for (int i = 0; i < (indc_cat_length - 1); i++) {
          sub_r1(i + 1) = indc_cat(i);
        }
        r1_indc = sub_r1;
        
        r2_indc = indc_cat;
      }
      f1_indc_num = rep(0, no_mainyears);
      f2_indc_num = rep(0, no_mainyears);
      
      {
        StringVector sub_r1(no_mainyears);
        StringVector sub_r2(no_mainyears);
        
        for (int i = 0; i < no_mainyears; i++) {
          sub_r1(i) = "none";
          sub_r2(i) = "none";
        }
        
        f1_indc_cat = sub_r1;
        f2_indc_cat = sub_r2;
      }
    }
  } else {
    int no_mainyears = mainyears_.length();
    indc_names = {"0"};
    
    f1_indc_num = rep(0, no_mainyears);
    f2_indc_num = rep(0, no_mainyears);
    
    StringVector sub_r1(no_mainyears);
    StringVector sub_r2(no_mainyears);
    
    for (int i = 0; i < no_mainyears; i++) {
      sub_r1(i) = "none";
      sub_r2(i) = "none";
    }
    
    r1_indc = sub_r1;
    r2_indc = sub_r2;
    f1_indc_cat = clone(sub_r1);
    f2_indc_cat = clone(sub_r2);
  }
  
  // Annual covariate vectors for function-based MPMs
  // Annual covariate a
  NumericVector f1_annua_num; // Numeric values entered as input
  NumericVector f2_annua_num;
  
  if (annua.isNotNull() && !raw) {
    RObject annua_input = as<RObject>(annua);
    
    NumericVector annua_num;
    
    int no_mainyears = mainyears_.length();
    
    if (!paramnames_provided) {
      throw Rcpp::exception("Use of annual covariates requires a valid paramnames object", false);
    }
    
    // Fixed covariate
    if (is<StringVector>(annua_input)) {
      String eat_my_shorts = "Annual covariate vector a must be empty, or include ";
      eat_my_shorts += "1, 2, or as many numbers as occasions in the dataset.";
      
      throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      
    } else if (is<IntegerVector>(annua_input) || is<NumericVector>(annua_input)) {
      // Handles possibility of factor variables
      annua_num = as<NumericVector>(annua_input);
      int annua_num_length = static_cast<int>(annua_num.length());
      
      if (annua_num_length != 1 && annua_num_length != 2) {
        if (annua_num_length != no_mainyears) {
          String eat_my_shorts = "Annual covariate vector a must be empty, or include ";
          String eat_my_shorts1 = "1, 2, or as many numbers as occasions in the dataset.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      
      StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
      StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
      
      if (annua_num_length == 1) {
        f1_annua_num = rep(annua_num(0), no_mainyears); 
        f2_annua_num = rep(annua_num(0), no_mainyears);
        
      } else if (annua_num_length == 2 && no_years != 2) {
        f1_annua_num = rep(annua_num(0), no_mainyears);
        f2_annua_num = rep(annua_num(1), no_mainyears);
        
      } else {
        NumericVector sub_f1(annua_num_length);
        sub_f1(0) = 0;
        
        for (int i = 0; i < (annua_num_length - 1); i++) {
          sub_f1(i + 1) = annua_num(i);
        }
        f1_annua_num = sub_f1;
        
        f2_annua_num = annua_num;
      }
      
    } else {
      throw Rcpp::exception("Format of annucova not recognized.", false);
    }
  } else {
    int no_mainyears = mainyears_.length();
    
    f1_annua_num = rep(0, no_mainyears);
    f2_annua_num = rep(0, no_mainyears);
  }
  
  // Annual covariate b
  NumericVector f1_annub_num; // Numeric values entered as input
  NumericVector f2_annub_num;
  
  if (annub.isNotNull() && !raw) {
    RObject annub_input = as<RObject>(annub);
    
    NumericVector annub_num;
    
    int no_mainyears = mainyears_.length();
    
    if (!paramnames_provided) {
      throw Rcpp::exception("Use of annual covariates requires a valid paramnames object", false);
    }
    
    // Fixed covariate
    if (is<StringVector>(annub_input)) {
      String eat_my_shorts = "Annual covariate vector b must be empty, or include ";
      eat_my_shorts += "1, 2, or as many numbers as occasions in the dataset.";
      
      throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      
    } else if (is<IntegerVector>(annub_input) || is<NumericVector>(annub_input)) {
      // Handles possibility of factor variables
      annub_num = as<NumericVector>(annub_input);
      int annub_num_length = static_cast<int>(annub_num.length());
      
      if (annub_num_length != 1 && annub_num_length != 2) {
        if (annub_num_length != no_mainyears) {
          String eat_my_shorts = "Annual covariate vector b must be empty, or include ";
          String eat_my_shorts1 = "1, 2, or as many numbers as occasions in the dataset.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      
      StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
      StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
      
      if (annub_num_length == 1) {
        f1_annub_num = rep(annub_num(0), no_mainyears); 
        f2_annub_num = rep(annub_num(0), no_mainyears);
        
      } else if (annub_num_length == 2 && no_years != 2) {
        f1_annub_num = rep(annub_num(0), no_mainyears);
        f2_annub_num = rep(annub_num(1), no_mainyears);
        
      } else {
        NumericVector sub_f1(annub_num_length);
        sub_f1(0) = 0;
        
        for (int i = 0; i < (annub_num_length - 1); i++) {
          sub_f1(i + 1) = annub_num(i);
        }
        f1_annub_num = sub_f1;
        
        f2_annub_num = annub_num;
      }
      
    } else {
      throw Rcpp::exception("Format of annucovb not recognized.", false);
    }
  } else {
    int no_mainyears = mainyears_.length();
    
    f1_annub_num = rep(0, no_mainyears);
    f2_annub_num = rep(0, no_mainyears);
  }
  
  // Annual covariate c
  NumericVector f1_annuc_num; // Numeric values entered as input
  NumericVector f2_annuc_num;
  
  if (annuc.isNotNull() && !raw) {
    RObject annuc_input = as<RObject>(annuc);
    
    NumericVector annuc_num;
    
    int no_mainyears = mainyears_.length();
    
    if (!paramnames_provided) {
      throw Rcpp::exception("Use of annual covariates requires a valid paramnames object", false);
    }
    
    // Fixed covariate
    if (is<StringVector>(annuc_input)) {
      String eat_my_shorts = "Annual covariate vector c must be empty, or include ";
      eat_my_shorts += "1, 2, or as many numbers as occasions in the dataset.";
      
      throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      
    } else if (is<IntegerVector>(annuc_input) || is<NumericVector>(annuc_input)) {
      // Handles possibility of factor variables
      annuc_num = as<NumericVector>(annuc_input);
      int annuc_num_length = static_cast<int>(annuc_num.length());
      
      if (annuc_num_length != 1 && annuc_num_length != 2) {
        if (annuc_num_length != no_mainyears) {
          String eat_my_shorts = "Annual covariate vector c must be empty, or include ";
          String eat_my_shorts1 = "1, 2, or as many numbers as occasions in the dataset.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
      
      StringVector mainparams = as<StringVector>(paramnames_["mainparams"]);
      StringVector modelparams = as<StringVector>(paramnames_["modelparams"]);
      
      if (annuc_num_length == 1) {
        f1_annuc_num = rep(annuc_num(0), no_mainyears); 
        f2_annuc_num = rep(annuc_num(0), no_mainyears);
        
      } else if (annuc_num_length == 2 && no_years != 2) {
        f1_annuc_num = rep(annuc_num(0), no_mainyears);
        f2_annuc_num = rep(annuc_num(1), no_mainyears);
        
      } else {
        NumericVector sub_f1(annuc_num_length);
        sub_f1(0) = 0;
        
        for (int i = 0; i < (annuc_num_length - 1); i++) {
          sub_f1(i + 1) = annuc_num(i);
        }
        f1_annuc_num = sub_f1;
        
        f2_annuc_num = annuc_num;
      }
      
    } else {
      throw Rcpp::exception("Format of annucovc not recognized.", false);
    }
  } else {
    int no_mainyears = mainyears_.length();
    
    f1_annuc_num = rep(0, no_mainyears);
    f2_annuc_num = rep(0, no_mainyears);
  }
  
  // Stageexpansions, data frame prep, MPM estimation
  DataFrame agestages;
  DataFrame ahstages;
  DataFrame hstages;
  List output_draft;
  
  NumericVector NA_empty = {NA_REAL};
  DataFrame NA_empty_df = DataFrame::create(_["X1"] = NA_empty);
  
  if (raw) {
    if (!historical) {
      if (stage && !age) {
        // Stage-only ahistorical raw
        IntegerVector removal_row = {melchett_stageframe_length};
        StringVector removal_var = {"stage_id"};
        DataFrame ahstages_now = LefkoUtils::df_remove(melchett_stageframe_,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        DataFrame stageexpansion3 = theoldpizzle(melchett_stageframe_, melchett_ovtable_,
          as<arma::mat>(melchett_repmatrix_), 0, 0, 1, 1, 0, 0);
        
        IntegerVector mel_sid = as<IntegerVector>(melchett_stageframe_["stage_id"]);
        
        DataFrame stageexpansion2;
        
        if (err_check) {
          NumericVector mel_sz2 = as<NumericVector>(melchett_stageframe_["sizebin_center"]);
          NumericVector mel_szb2 = as<NumericVector>(melchett_stageframe_["sizebinb_center"]);
          NumericVector mel_szc2 = as<NumericVector>(melchett_stageframe_["sizebinc_center"]);
          NumericVector mel_rep = as<NumericVector>(melchett_stageframe_["repstatus"]);
          NumericVector mel_ind = as<NumericVector>(melchett_stageframe_["indataset"]);
          IntegerVector mel_ind2 = mel_sid - 1;
          
          int repmat_rows = static_cast<int>(melchett_repmatrix_.nrow());
          IntegerVector mel_fec3 (melchett_stageframe_length);
          for (int i = 0; i < repmat_rows; i++) {
            if (sum(melchett_repmatrix_.row(i)) > 0.0) mel_fec3 = 1;
          }
          
          stageexpansion2 = DataFrame::create(_["stage2"] = mel_sid,
            _["size2"] = mel_sz2, _["sizeb2"] = mel_szb2, _["sizec2"] = mel_szc2,
            _["rep2"] = mel_rep, _["indata2"] = mel_ind, _["index2"] = mel_ind2,
            _["fec3"] = mel_fec3);
        }
        
        StringVector d_rn = data_.attr("row.names");
        
        IntegerVector index3;
        IntegerVector index2;
        if (!new_stages_needed) {
          index3 = as<IntegerVector>(data_["stage3index"]) - 1;
          index2 = as<IntegerVector>(data_["stage2index"]) - 1;
          
          data_["index3"] = index3;
          data_["index2"] = index2;
        } else {
          index3 = as<IntegerVector>(data_["index3"]);
          index2 = as<IntegerVector>(data_["index2"]);
        }
        
        IntegerVector index32 (data_points);
        for (int i = 0; i < data_points; i++) {
          if (index3(i) < 0) index3(i) = melchett_stageframe_length - 1;
          if (index2(i) < 0) index2(i) = melchett_stageframe_length - 1;
          index32(i) = index3(i) + (index2(i) * melchett_stageframe_length);
        }
        
        if (fectime == 3) {
          data_["usedfec"] = data_[fec_int(0)];
        } else {
          data_["usedfec"] = data_[fec_int(1)];
        }
        
        data_["index32"] = index32;
        data_.attr("class") = "data.frame";
        data_.attr("row.names") = d_rn;
        
        // Matrix estimation
        List madsexmadrigal = normalpatrolgroup(stageexpansion3,
          as<arma::ivec>(mel_sid), data_, melchett_stageframe_, err_check,
          loy_pop_, loy_patch_, loy_year2_, yearorder_, pop_var_int,
          patch_var_int, year_var_int, loy_pop_used, loy_patch_used, simple,
          sparse_output);
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(madsexmadrigal, mat_qc, ahstages_now, NA_empty_df,
          NA_empty_df, false, true, false, reduce, simple, sparse_output);
        
        madsexmadrigal["ahstages"] = ahstages_now;
        madsexmadrigal["agestages"] = NA_empty_df;
        madsexmadrigal["hstages"] = NA_empty_df;
        madsexmadrigal["labels"] = labels_;
        madsexmadrigal["matrixqc"] = mat_qc;
        madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          madsexmadrigal["sge3"] = stageexpansion3;
          madsexmadrigal["sge2"] = stageexpansion2;
          madsexmadrigal["supplement"] = melchett_ovtable_;
        }
        output_draft = madsexmadrigal;
        
      } else if (stage && age) {
        // Age-stage ahistorical raw
        IntegerVector removal_row = {melchett_stageframe_length};
        StringVector removal_var = {"stage_id"};
        DataFrame ahstages_now = LefkoUtils::df_remove(melchett_stageframe_,
          removal_row, false, true, false, false, true,
          as<RObject>(removal_var));
        
        DataFrame stageexpansion3 = theoldpizzle(melchett_stageframe_, melchett_ovtable_,
          as<arma::mat>(melchett_repmatrix_), start_age, last_age, 1, 2, cont, 0);
        
        IntegerVector agevec = seq(start_age, last_age);
        int totalages = static_cast<int>(agevec.length());
        
        IntegerVector mel_sid = rep(as<IntegerVector>(melchett_stageframe_["stage_id"]), totalages);
        
        IntegerVector mel_ages (totalages * melchett_stageframe_length);
        IntegerVector ahage_stage_id (totalages * (melchett_stageframe_length - 1));
        StringVector ahage_stage (totalages * (melchett_stageframe_length - 1));
        IntegerVector ahage_age (totalages * (melchett_stageframe_length - 1));
        
        for (int i = 0; i < totalages; i++) {
          for (int j = 0; j < melchett_stageframe_length; j++) {
            mel_ages(j + (i * melchett_stageframe_length)) = agevec(i);
            
            if (j != (melchett_stageframe_length - 1)) {
              ahage_stage_id(j + (i * (melchett_stageframe_length - 1))) = mel_sid(j);
              ahage_stage(j + (i * (melchett_stageframe_length - 1))) = melchett_stageframe_stage_(j);
              ahage_age(j + (i * (melchett_stageframe_length - 1))) = agevec(i);
            }
          }
        }
        DataFrame agestages_now = DataFrame::create(_["stage_id"] = ahage_stage_id,
          _["stage"] = ahage_stage, _["age"] = ahage_age);
        
        IntegerVector mel_idx2 = mel_sid - 1;
        
        IntegerVector mel_idx21;
        {
          IntegerVector mel_idx21_1 = mel_ages - start_age;
          IntegerVector mel_idx21_2 = mel_idx21_1 * melchett_stageframe_length;
          mel_idx21 = mel_idx2 + mel_idx21_2;
        }
        
        if (err_check) {
          NumericVector mel_sz2 = rep(as<NumericVector>(melchett_stageframe_["sizebin_center"]), totalages);
          NumericVector mel_szb2 = rep(as<NumericVector>(melchett_stageframe_["sizebinb_center"]), totalages);
          NumericVector mel_szc2 = rep(as<NumericVector>(melchett_stageframe_["sizebinc_center"]), totalages);
          NumericVector mel_rep = rep(as<NumericVector>(melchett_stageframe_["repstatus"]), totalages);
          NumericVector mel_ind = rep(as<NumericVector>(melchett_stageframe_["indataset"]), totalages);
          
          NumericVector fec_sums;
          {
            arma::vec rep_sums = arma::sum(as<arma::mat>(melchett_repmatrix_), 1);
            int mle_rep_dim = as<arma::mat>(melchett_repmatrix_).n_rows;
            
            rep_sums.resize(mle_rep_dim + 1);
            arma::uvec fec3_nonzeros = find(rep_sums);
            rep_sums.elem(fec3_nonzeros).ones();
            
            fec_sums = as<NumericVector>(wrap(rep_sums));
          }
          
          DataFrame stageexpansion2 = DataFrame::create(_["stage2"] = mel_sid,
            _["size2"] = mel_sz2, _["sizeb2"] = mel_szb2, _["sizec2"] = mel_szc2,
            _["rep2"] = mel_rep, _["indata2"] = mel_ind, _["index2"] = mel_idx2,
            _["fec3"] = fec_sums, _["age2"] = mel_ages, _["index21"] = mel_idx21);
        }
        
        IntegerVector usedindex3 = data_["index3"];
        IntegerVector usedindex2 = data_["index2"];
        IntegerVector usedage2 = data_[age_var_int];
        
        IntegerVector index321 (data_points);
        IntegerVector index21 (data_points);
        
        for (int i = 0; i < data_points; i++) {
          if (usedindex3(i) < 0) usedindex3(i) = melchett_stageframe_length - 1;
          if (usedindex2(i) < 0) usedindex2(i) = melchett_stageframe_length - 1;
          
          index321(i) = usedindex3(i) + (((usedage2(i) + 1) - start_age) * 
              melchett_stageframe_length) +
            (usedindex2(i) * melchett_stageframe_length * totalages) +
            ((usedage2(i) - start_age) * melchett_stageframe_length *
              melchett_stageframe_length * totalages);
          index21(i) = usedindex2(i) + ((usedage2(i) - start_age) *
              melchett_stageframe_length);
        }
        
        StringVector d_rn = data_.attr("row.names");
        
        data_["usedage"] = usedage2;
        data_["index321"] = index321;
        data_["index21"] = index21;
        
        if (fectime == 3) {
          data_["usedfec"] = data_[fec_int(0)];
        } else {
          data_["usedfec"] = data_[fec_int(1)];
        }
        
        data_.attr("class") = "data.frame";
        data_.attr("row.names") = d_rn;
        
        List madsexmadrigal = subvertedpatrolgroup(stageexpansion3,
          as<arma::ivec>(mel_idx21), data_, melchett_stageframe_, start_age,
          last_age, cont, err_check, loy_pop_, loy_patch_, loy_year2_,
          yearorder_, pop_var_int, patch_var_int, year_var_int, loy_pop_used,
          loy_patch_used, simple, sparse_output);
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(madsexmadrigal, mat_qc, ahstages_now, NA_empty_df,
          agestages_now, true, true, false, reduce, simple, sparse_output);
        
        madsexmadrigal["ahstages"] = ahstages_now;
        madsexmadrigal["agestages"] = agestages_now;
        madsexmadrigal["hstages"] = NA_empty_df;
        madsexmadrigal["labels"] = labels_;
        madsexmadrigal["matrixqc"] = mat_qc;
        madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          madsexmadrigal["sge3"] = stageexpansion3;
          madsexmadrigal["supplement"] = melchett_ovtable_;
        }
        output_draft = madsexmadrigal;
        
      } else if (!stage && age) {
        // Age-only raw
        ahstages = melchett_stageframe_;
        StringVector d_rn = data_.attr("row.names");
        
        IntegerVector usedage2 = data_[age_var_int];
        data_["usedage"] = usedage2;
        
        if (fectime == 3) {
          data_["usedfec"] = data_[fec_int(0)];
        } else {
          data_["usedfec"] = data_[fec_int(1)];
        }
        
        data_.attr("class") = "data.frame";
        data_.attr("row.names") = d_rn;
        
        // Matrix estimation
        List madsexmadrigal = minorpatrolgroup(data_, melchett_stageframe_,
          melchett_ovtable_, start_age, last_age, cont, fecmod, err_check,
          loy_pop_, loy_patch_, loy_year2_, yearorder_, pop_var_int,
          patch_var_int, year_var_int, loy_pop_used, loy_patch_used, simple,
          sparse_output);
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(madsexmadrigal, mat_qc, ahstages, NA_empty_df,
          NA_empty_df, true, false, false, reduce, simple, sparse_output);
        
        madsexmadrigal["ahstages"] = ahstages;
        madsexmadrigal["agestages"] = NA_empty_df;
        madsexmadrigal["hstages"] = NA_empty_df;
        madsexmadrigal["labels"] = labels_;
        madsexmadrigal["matrixqc"] = mat_qc;
        madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          //madsexmadrigal["sge3"] = stageexpansion3;
          //madsexmadrigal["sge2"] = stageexpansion2;
          madsexmadrigal["supplement"] = melchett_ovtable_;
        }
        output_draft = madsexmadrigal;
        
      }
    } else {
      // Historical stage-only raw
      if (stage && !age) {
        List sge3;
        
        IntegerVector removal_row = {melchett_stageframe_length};
        StringVector removal_var = {"stage_id"};
        DataFrame ahstages_now = LefkoUtils::df_remove(melchett_stageframe_,
          removal_row, false, true, false, false, true,
          as<RObject>(removal_var));
        
        StringVector mel_ovt_eststage1 = as<StringVector>(melchett_ovtable_["eststage1"]);
        StringVector mel_ovt_eststage2 = as<StringVector>(melchett_ovtable_["eststage2"]);
        StringVector mel_ovt_eststage3 = as<StringVector>(melchett_ovtable_["eststage3"]);
        
        StringVector mel_ovt_stage1 = as<StringVector>(melchett_ovtable_["stage1"]);
        StringVector mel_ovt_stage2 = as<StringVector>(melchett_ovtable_["stage2"]);
        StringVector mel_ovt_stage3 = as<StringVector>(melchett_ovtable_["stage3"]);
        
        IntegerVector usedstage3index = data_["stage3index"];
        IntegerVector usedstage2index = data_["stage2index"];
        IntegerVector usedstage1index = data_["stage1index"];
        
        IntegerVector usedindex3 = data_["index3"];
        IntegerVector usedindex2 = data_["index2"];
        IntegerVector usedindex1 = data_["index1"];
        
        {
          StringVector data_stage1 = as<StringVector>(data_["stage1"]);
          StringVector data_stage2 = as<StringVector>(data_["stage2"]);
          StringVector data_stage3 = as<StringVector>(data_["stage3"]);
          
          for (int j = 0; j < usedindex3.length(); j++) {
            for (int k = 0; k < melchett_stageframe_length; k++) {
              if (stringcompare_hard(String(data_stage3(j)), 
                  String(melchett_stageframe_stage_(k)))) {
                usedstage3index(j) = melchett_stageframe_stageid_(k);
                usedindex3(j) = melchett_stageframe_stageid_(k) - 1;
              }
              
              if (stringcompare_hard(String(data_stage2(j)), 
                  String(melchett_stageframe_stage_(k)))) {
                usedstage2index(j) = melchett_stageframe_stageid_(k);
                usedindex2(j) = melchett_stageframe_stageid_(k) - 1;
              }
              
              if (stringcompare_hard(String(data_stage1(j)), 
                      String(melchett_stageframe_stage_(k)))) {
                usedstage1index(j) = melchett_stageframe_stageid_(k);
                usedindex1(j) = melchett_stageframe_stageid_(k) - 1;
              }
            }
          }
          
          data_["index3"] = usedindex3;
          data_["index2"] = usedindex2;
          data_["index1"] = usedindex1;
          
          data_["stage1"] = data_stage1;
          data_["stage2"] = data_stage2;
          data_["stage3"] = data_stage3;
        }
        
        IntegerVector mel_ovt_es1_notalive (melchett_ovtable_length);
        int mov_notalive_count {0};
        
        // Replace transitions involving NotAlive in t-1 with replacements in ovtable
        for (int i = 0; i < melchett_ovtable_length; i++) {
          String nolovelost = String(mel_ovt_eststage1(i));
          if (LefkoUtils::stringcompare_simple(nolovelost, "notalive", true)) {
            mel_ovt_es1_notalive(i) = 1;
            mov_notalive_count++;
          }
        }
        
        if (mov_notalive_count > 0) {
          arma::uvec flubble_indices = find(as<arma::ivec>(mel_ovt_es1_notalive));
          
          StringVector data_stage1 = as<StringVector>(data_["stage1"]);
          StringVector data_stage2 = as<StringVector>(data_["stage2"]);
          StringVector data_stage3 = as<StringVector>(data_["stage3"]);
          
          IntegerVector data_stage1_notalive = LefkoUtils::index_l3(data_stage1, "NotAlive");
          int all_stage1_notalives = static_cast<int>(data_stage1_notalive.length());
          
          for (int i = 0; i < mov_notalive_count; i++) {
            for (int j = 0; j < all_stage1_notalives; j++) {
              
              if (LefkoUtils::stringcompare_hard(String(data_stage2(data_stage1_notalive(j))),   
                  String(mel_ovt_eststage2(flubble_indices(i))))) {
                if (LefkoUtils::stringcompare_hard(String(data_stage3(data_stage1_notalive(j))), 
                    String(mel_ovt_eststage3(flubble_indices(i))))) {
                  
                  data_stage3(data_stage1_notalive(j)) = mel_ovt_stage3(flubble_indices(i));
                  data_stage2(data_stage1_notalive(j)) = mel_ovt_stage2(flubble_indices(i));
                  data_stage1(data_stage1_notalive(j)) = mel_ovt_stage1(flubble_indices(i));
                  
                  for (int k = 0; k < melchett_stageframe_length; k++) {
                    if (stringcompare_hard(String(data_stage3(data_stage1_notalive(j))), 
                        String(melchett_stageframe_stage_(k)))) {
                      usedstage3index(data_stage1_notalive(j)) = melchett_stageframe_stageid_(k);
                      usedindex3(data_stage1_notalive(j)) = melchett_stageframe_stageid_(k) - 1;
                    }
                    
                    if (stringcompare_hard(String(data_stage2(data_stage1_notalive(j))), 
                        String(melchett_stageframe_stage_(k)))) {
                      usedstage2index(data_stage1_notalive(j)) = melchett_stageframe_stageid_(k);
                      usedindex2(data_stage1_notalive(j)) = melchett_stageframe_stageid_(k) - 1;
                    }
                    
                    if (stringcompare_hard(String(data_stage1(data_stage1_notalive(j))), 
                        String(melchett_stageframe_stage_(k)))) {
                      usedstage1index(data_stage1_notalive(j)) = melchett_stageframe_stageid_(k);
                      usedindex1(data_stage1_notalive(j)) = melchett_stageframe_stageid_(k) - 1;
                    }
                  }
                }
                
                data_["stage3index"] = usedstage3index;
                data_["stage2index"] = usedstage2index;
                data_["stage1index"] = usedstage1index;
                
                data_["index3"] = usedindex3;
                data_["index2"] = usedindex2;
                data_["index1"] = usedindex1;
                
                data_["stage1"] = data_stage1;
                data_["stage2"] = data_stage2;
                data_["stage3"] = data_stage3;
              }
            }
          }
          // Remove ovtable transitions involving movement from NotAlive in time t-1
          DataFrame movt_new = clone(melchett_ovtable_);
          movt_new["to_remove"] = mel_ovt_es1_notalive;
          melchett_ovtable_ = movt_new;
          
          StringVector mel_check_name = {"to_remove", "to_remove"};
          IntegerVector mel_check_int = {1, 1};
          DataFrame mel_ov_new = LefkoUtils::df_remove(melchett_ovtable_,
            as<RObject>(mel_check_int), false, true, false, false, true,
            as<RObject>(mel_check_name));
          melchett_ovtable_ = mel_ov_new;
        }
        
        // Large and small matrix element indices
        DataFrame stageexpansion9 = theoldpizzle(melchett_stageframe_, melchett_ovtable_,
          as<arma::mat>(melchett_repmatrix_), 0, 0, format_int, 0, 0, 0);
        
        IntegerVector mel_sid = as<IntegerVector>(melchett_stageframe_["stage_id"]);
        List se36 = LefkoUtils::exp_grd(as<RObject>(mel_sid), as<RObject>(mel_sid));
        IntegerVector se3_index21 = (as<IntegerVector>(se36(0)) - 1) + 
          ((as<IntegerVector>(se36(1))) - 1) * melchett_stageframe_length;
        
        if (err_check) {
          NumericVector mel_sz2 = as<NumericVector>(melchett_stageframe_["sizebin_center"]);
          NumericVector mel_szb2 = as<NumericVector>(melchett_stageframe_["sizebinb_center"]);
          NumericVector mel_szc2 = as<NumericVector>(melchett_stageframe_["sizebinc_center"]);
          NumericVector mel_rep = as<NumericVector>(melchett_stageframe_["repstatus"]);
          NumericVector mel_ind = as<NumericVector>(melchett_stageframe_["indataset"]);
          List se31 = LefkoUtils::exp_grd(as<RObject>(mel_sz2), as<RObject>(mel_sz2));
          List se32 = LefkoUtils::exp_grd(as<RObject>(mel_szb2), as<RObject>(mel_szb2));
          List se33 = LefkoUtils::exp_grd(as<RObject>(mel_szc2), as<RObject>(mel_szc2));
          List se34 = LefkoUtils::exp_grd(as<RObject>(mel_rep), as<RObject>(mel_rep));
          List se35 = LefkoUtils::exp_grd(as<RObject>(mel_ind), as<RObject>(mel_ind));
          
          List stageexpansion3 (15);
          stageexpansion3(0) = as<NumericVector>(se31(0));
          stageexpansion3(1) = as<NumericVector>(se31(1));
          stageexpansion3(2) = as<NumericVector>(se32(0));
          stageexpansion3(3) = as<NumericVector>(se32(1));
          stageexpansion3(4) = as<NumericVector>(se33(0));
          stageexpansion3(5) = as<NumericVector>(se33(1));
          stageexpansion3(6) = as<NumericVector>(se34(0));
          stageexpansion3(7) = as<NumericVector>(se34(1));
          stageexpansion3(8) = as<NumericVector>(se35(0));
          stageexpansion3(9) = as<NumericVector>(se35(1));
          stageexpansion3(10) = as<IntegerVector>(se36(0));
          stageexpansion3(11) = as<IntegerVector>(se36(1));
          
          NumericVector harpoon;
          {
            arma::mat melchett_repmatrix_arma = as<arma::mat>(melchett_repmatrix_);
            if (devries) {
              arma::mat mra_colSums = arma::sum(melchett_repmatrix_arma, 0);
              
              arma::uvec mra_cS_ones = find(mra_colSums);
              mra_colSums.elem(mra_cS_ones).ones();
              int mra_dim = static_cast<int>(mra_colSums.n_elem);
              
              arma::mat mra_zerorow = zeros(1, (mra_dim));
              arma::mat mra_zero2col = zeros((mra_dim + 2), 2);
              
              melchett_repmatrix_arma.insert_rows(mra_dim, mra_colSums);
              melchett_repmatrix_arma.insert_rows((mra_dim + 1), mra_zerorow);
              melchett_repmatrix_arma.insert_cols((mra_dim), mra_zero2col);
              
              arma::vec mra_vec = arma::vectorise(melchett_repmatrix_arma);
              harpoon = as<NumericVector>(wrap(mra_vec));
            } else {
              int mra_dim = static_cast<int>(melchett_repmatrix_arma.n_rows);
              
              arma::mat mra_zerorow = zeros(1, mra_dim);
              arma::mat mra_zerocol = zeros((mra_dim + 1), 1);
              
              melchett_repmatrix_arma.insert_rows((mra_dim), mra_zerorow);
              melchett_repmatrix_arma.insert_cols((mra_dim), mra_zerocol);
              
              arma::vec mra_vec = arma::vectorise(melchett_repmatrix_arma);
              harpoon = as<NumericVector>(wrap(mra_vec));
            }
          }
          stageexpansion3(12) = harpoon;
          
          NumericVector se3_indata32n = as<NumericVector>(se35(0)) * as<NumericVector>(se35(1));
          stageexpansion3(13) = se3_indata32n;
          stageexpansion3(14) = se3_index21;
          
          StringVector se3_names = {"size3", "size2n", "sizeb3", "sizeb2n",
            "sizec3", "sizec2n", "rep3", "rep2n", "indata3", "indata2n",
            "stage3", "stage2n", "fec32n", "indata32n", "index21"};
          stageexpansion3.attr("names") = se3_names;
          stageexpansion3.attr("class") = "data.frame";
          
          stageexpansion3.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, 
              (melchett_stageframe_length * melchett_stageframe_length));
          
          if (err_check) sge3 = stageexpansion3;
        }
        
        
        
        
        int format_alteration {0};
        if (!devries) format_alteration = 1;
        
        IntegerVector hst_sid_2 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        IntegerVector hst_sid_1 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        StringVector hst_stage_2 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        StringVector hst_stage_1 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        
        int hst_n_counter {0};
        for (int i = 0; i < (melchett_stageframe_length - 1); i++) {
          for (int j = 0; j < (melchett_stageframe_length - 2 + format_alteration); j++) {
            hst_sid_2(hst_n_counter) = mel_sid(j);
            hst_sid_1(hst_n_counter) = mel_sid(i);
            
            hst_stage_2(hst_n_counter) = melchett_stageframe_stage_(j);
            hst_stage_1(hst_n_counter) = melchett_stageframe_stage_(i);
            
            hst_n_counter++;
          }
        }
        
        DataFrame hstages_now = DataFrame::create(_["stage_id_2"] = hst_sid_2,
          _["stage_id_1"] = hst_sid_1, _["stage_2"] = hst_stage_2,
          _["stage_1"] = hst_stage_1);
        
        IntegerVector index321 (data_points);
        IntegerVector index21 (data_points);
        
        if (format_int == 2) {
          for (int i = 0; i < data_points; i++) {
            if (usedindex3(i) < 0) usedindex3(i) = melchett_stageframe_length - 1;
            if (usedindex2(i) < 0) usedindex2(i) = melchett_stageframe_length - 1;
            if (usedindex1(i) < 0) usedindex1(i) = melchett_stageframe_length - 1;
            
            index321(i) = usedindex3(i) + (usedindex2(i) * melchett_stageframe_length) +
              (usedindex2(i) * melchett_stageframe_length * melchett_stageframe_length) +
              (usedindex1(i) * melchett_stageframe_length * melchett_stageframe_length *
              melchett_stageframe_length);
            index21(i) = usedindex2(i) + (usedindex1(i) * melchett_stageframe_length);
          }
        } else {
          int sl_small = melchett_stageframe_length - 1;
          for (int i = 0; i < data_points; i++) {
            if (usedindex3(i) < 0) usedindex3(i) = melchett_stageframe_length - 1;
            if (usedindex2(i) < 0) usedindex2(i) = melchett_stageframe_length - 1;
            if (usedindex1(i) < 0) usedindex1(i) = melchett_stageframe_length - 1;
            
            index321(i) = usedindex3(i) + (usedindex2(i) * sl_small) +
              (usedindex2(i) * sl_small * sl_small) +
              (usedindex1(i) * sl_small * sl_small * sl_small);
            index21(i) = usedindex2(i) + (usedindex1(i) * melchett_stageframe_length);
          }
        }
        
        StringVector d_rn = data_.attr("row.names");
        
        data_["index321"] = index321;
        data_["pairindex21"] = index21;
        
        if (fectime == 3) {
          data_["usedfec"] = data_[fec_int(0)];
        } else {
          data_["usedfec"] = data_[fec_int(1)];
        }
        
        data_.attr("class") = "data.frame";
        data_.attr("row.names") = d_rn;
        
        // Matrix estimation
        List madsexmadrigal = specialpatrolgroup(stageexpansion9,
          as<arma::ivec>(se3_index21), data_, melchett_stageframe_, format_int,
          err_check, loy_pop_, loy_patch_, loy_year2_, yearorder_, pop_var_int,
          patch_var_int, year_var_int, loy_pop_used, loy_patch_used, simple,
          sparse_output);
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(madsexmadrigal, mat_qc, ahstages_now, hstages_now,
          NA_empty_df, false, true, true, reduce, simple, sparse_output);
        
        madsexmadrigal["ahstages"] = ahstages_now;
        madsexmadrigal["agestages"] = NA_empty_df;
        madsexmadrigal["hstages"] = hstages_now;
        madsexmadrigal["labels"] = labels_;
        madsexmadrigal["matrixqc"] = mat_qc;
        madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          madsexmadrigal["sge9"] = stageexpansion9;
          madsexmadrigal["sge3"] = sge3;
          madsexmadrigal["supplement"] = melchett_ovtable_;
        }
        output_draft = madsexmadrigal;
      }
    }
  } else {
    // Function-based MPMs
    
    if (!historical) {
      if (stage && !age) {
        // Stage-only ahistorical function-based
        IntegerVector removal_row = {melchett_stageframe_length};
        StringVector removal_var = {"stage_id"};
        DataFrame ahstages_now = LefkoUtils::df_remove(melchett_stageframe_,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        List new_madsexmadrigal = raymccooney(list_of_years, modelsuite_,
          mainyears_, mainpatches_, as<RObject>(maingroups_),
          as<RObject>(inda_names), as<RObject>(indb_names),
          as<RObject>(indc_names), melchett_stageframe_, melchett_ovtable_,
          as<arma::mat>(melchett_repmatrix_), f2_inda_num, f1_inda_num,
          f2_indb_num, f1_indb_num, f2_indc_num, f1_indc_num, f2_annua_num,
          f1_annua_num, f2_annub_num, f1_annub_num, f2_annuc_num, f1_annuc_num,
          f2_inda_cat, f1_inda_cat, f2_indb_cat, f1_indb_cat, f2_indc_cat,
          f1_indc_cat, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc,
          dev_terms_, density, fecmod, 0, 0, 1, 1, 0, 1, negfec, nodata,
          exp_tol, theta_tol, CDF, err_check, simple, sparse_output);
        if (err_check) err_check_proxies = new_madsexmadrigal["proxies"];
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(new_madsexmadrigal, mat_qc, ahstages_now, NA_empty_df,
          NA_empty_df, false, true, false, reduce, simple, sparse_output);
        
        new_madsexmadrigal["ahstages"] = ahstages_now;
        new_madsexmadrigal["agestages"] = NA_empty_df;
        new_madsexmadrigal["hstages"] = NA_empty_df;
        new_madsexmadrigal["labels"] = labels_;
        new_madsexmadrigal["matrixqc"] = mat_qc;
        new_madsexmadrigal["modelqc"] = mod_qc_;
        new_madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          new_madsexmadrigal["proxies"] = err_check_proxies;
        }
        output_draft = new_madsexmadrigal;
        
      } else if (stage && age) {
        // Age-stage function-based
        IntegerVector removal_row = {melchett_stageframe_length};
        StringVector removal_var = {"stage_id"};
        DataFrame ahstages_now = LefkoUtils::df_remove(melchett_stageframe_,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        IntegerVector agevec = seq(start_age, last_age);
        
        int totalages = static_cast<int>(agevec.length());
        
        IntegerVector mel_sid = rep(as<IntegerVector>(melchett_stageframe_["stage_id"]), totalages);
        
        IntegerVector mel_ages (totalages * melchett_stageframe_length);
        IntegerVector ahage_stage_id (totalages * (melchett_stageframe_length - 1));
        StringVector ahage_stage (totalages * (melchett_stageframe_length - 1));
        IntegerVector ahage_age (totalages * (melchett_stageframe_length - 1));
        
        for (int i = 0; i < totalages; i++) {
          for (int j = 0; j < melchett_stageframe_length; j++) {
            mel_ages(j + (i * melchett_stageframe_length)) = agevec(i);
            
            if (j != (melchett_stageframe_length - 1)) {
              ahage_stage_id(j + (i * (melchett_stageframe_length - 1))) = mel_sid(j);
              ahage_stage(j + (i * (melchett_stageframe_length - 1))) = melchett_stageframe_stage_(j);
              ahage_age(j + (i * (melchett_stageframe_length - 1))) = agevec(i);
            }
          }
        }
        DataFrame agestages_now = DataFrame::create(_["stage_id"] = ahage_stage_id,
          _["stage"] = ahage_stage, _["age"] = ahage_age);
        
        List new_madsexmadrigal = raymccooney(list_of_years, modelsuite_,
          mainyears_, mainpatches_, as<RObject>(maingroups_),
          as<RObject>(inda_names), as<RObject>(indb_names), as<RObject>(indc_names),
          melchett_stageframe_, melchett_ovtable_, as<arma::mat>(melchett_repmatrix_),
          f2_inda_num, f1_inda_num, f2_indb_num, f1_indb_num, f2_indc_num,
          f1_indc_num, f2_annua_num, f1_annua_num, f2_annub_num, f1_annub_num,
          f2_annuc_num, f1_annuc_num, f2_inda_cat, f1_inda_cat, f2_indb_cat, f1_indb_cat,
          f2_indc_cat, f1_indc_cat, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
          r1_indc, dev_terms_, density, fecmod, start_age, last_age, 1, 2, cont,
          2, negfec, nodata, exp_tol, theta_tol, CDF, err_check, simple,
          sparse_output);
        if (err_check) err_check_proxies = new_madsexmadrigal["proxies"];
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(new_madsexmadrigal, mat_qc, ahstages_now, NA_empty_df,
          agestages_now, true, true, false, reduce, simple, sparse_output);
        
        new_madsexmadrigal["ahstages"] = ahstages_now;
        new_madsexmadrigal["agestages"] = agestages_now;
        new_madsexmadrigal["hstages"] = NA_empty_df;
        new_madsexmadrigal["labels"] = labels_;
        new_madsexmadrigal["matrixqc"] = mat_qc;
        new_madsexmadrigal["modelqc"] = mod_qc_;
        new_madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          new_madsexmadrigal["proxies"] = err_check_proxies;
          new_madsexmadrigal["supplement"] = melchett_ovtable_;
        }
        output_draft = new_madsexmadrigal;
        
      } else if (!stage && age) {
        // Age-only function-based
        ahstages = melchett_stageframe_;
        IntegerVector actualages = seq(start_age, last_age);
        
        List new_madsexmadrigal = mothermccooney(list_of_years, modelsuite_,
          actualages, mainyears_, mainpatches_, as<RObject>(maingroups_),
          as<RObject>(inda_names), as<RObject>(indb_names),
          as<RObject>(indc_names), melchett_stageframe_, melchett_ovtable_,
          f2_inda_num, f1_inda_num, f2_indb_num, f1_indb_num, f2_indc_num,
          f1_indc_num, f2_annua_num, f1_annua_num, f2_annub_num, f1_annub_num,
          f2_annuc_num, f1_annuc_num, f2_inda_cat, f1_inda_cat, f2_indb_cat, f1_indb_cat,
          f2_indc_cat, f1_indc_cat, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
          r1_indc, dev_terms_, density, fecmod, last_age, cont, negfec, nodata,
          exp_tol, theta_tol, err_check, simple, sparse_output);
        if (err_check) err_check_proxies = new_madsexmadrigal["proxies"];
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(new_madsexmadrigal, mat_qc, ahstages, NA_empty_df,
          NA_empty_df, true, false, false, reduce, simple, sparse_output);
        
        new_madsexmadrigal["ahstages"] = ahstages;
        new_madsexmadrigal["agestages"] = NA_empty_df;
        new_madsexmadrigal["hstages"] = NA_empty_df;
        new_madsexmadrigal["labels"] = labels_;
        new_madsexmadrigal["matrixqc"] = mat_qc;
        new_madsexmadrigal["modelqc"] = mod_qc_;
        new_madsexmadrigal["dataqc"] = dataqc_;
        output_draft = new_madsexmadrigal;
        
      }
    } else {
      if (stage && !age) {
        // Historical stage-only function-based
        IntegerVector removal_row = {melchett_stageframe_length};
        StringVector removal_var = {"stage_id"};
        DataFrame ahstages_now = LefkoUtils::df_remove(melchett_stageframe_,
          removal_row, false, true, false, false, true, as<RObject>(removal_var));
        
        int format_alteration {0};
        if (!devries) format_alteration = 1;
        
        IntegerVector mel_sid = as<IntegerVector>(melchett_stageframe_["stage_id"]);
        IntegerVector hst_sid_2 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        IntegerVector hst_sid_1 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        StringVector hst_stage_2 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        StringVector hst_stage_1 ((melchett_stageframe_length - 2 + format_alteration) * (melchett_stageframe_length - 1));
        
        int hst_n_counter {0};
        for (int i = 0; i < (melchett_stageframe_length - 1); i++) {
          for (int j = 0; j < (melchett_stageframe_length - 2 + format_alteration); j++) {
            hst_sid_2(hst_n_counter) = mel_sid(j);
            hst_sid_1(hst_n_counter) = mel_sid(i);
            
            hst_stage_2(hst_n_counter) = melchett_stageframe_stage_(j);
            hst_stage_1(hst_n_counter) = melchett_stageframe_stage_(i);
            
            hst_n_counter++;
          }
        }
        
        DataFrame hstages_now = DataFrame::create(_["stage_id_2"] = hst_sid_2,
          _["stage_id_1"] = hst_sid_1, _["stage_2"] = hst_stage_2,
          _["stage_1"] = hst_stage_1);
        
        List new_madsexmadrigal = raymccooney(list_of_years, modelsuite_,
          mainyears_, mainpatches_, as<RObject>(maingroups_),
          as<RObject>(inda_names), as<RObject>(indb_names), as<RObject>(indc_names),
          melchett_stageframe_, melchett_ovtable_, as<arma::mat>(melchett_repmatrix_),
          f2_inda_num, f1_inda_num, f2_indb_num, f1_indb_num, f2_indc_num,
          f1_indc_num, f2_annua_num, f1_annua_num, f2_annub_num, f1_annub_num,
          f2_annuc_num, f1_annuc_num, f2_inda_cat, f1_inda_cat, f2_indb_cat, f1_indb_cat,
          f2_indc_cat, f1_indc_cat, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc,
          r1_indc, dev_terms_, density, fecmod, 0, 0, format_int, 0, 0, 1, negfec,
          nodata, exp_tol, theta_tol, CDF, err_check, simple, sparse_output);
        if (err_check) err_check_proxies = new_madsexmadrigal["proxies"];
        
        IntegerVector mat_qc = {0, 0, 0};
        LefkoUtils::matrix_reducer(new_madsexmadrigal, mat_qc, ahstages_now, hstages_now,
          NA_empty_df, false, true, true, reduce, simple, sparse_output);
        
        new_madsexmadrigal["ahstages"] = ahstages_now;
        new_madsexmadrigal["agestages"] = NA_empty_df;
        new_madsexmadrigal["hstages"] = hstages_now;
        new_madsexmadrigal["labels"] = labels_;
        new_madsexmadrigal["matrixqc"] = mat_qc;
        new_madsexmadrigal["modelqc"] = mod_qc_;
        new_madsexmadrigal["dataqc"] = dataqc_;
        if (err_check) {
          new_madsexmadrigal["proxies"] = err_check_proxies;
        }
        output_draft = new_madsexmadrigal;
      }
    }
  }
  
  StringVector out_classes = {"lefkoMat"};
  output_draft.attr("class") = out_classes;
  
  return output_draft;
}

// Pop Dynamics

//' Estimate Stable Stage Distribution of Any Population Matrix
//' 
//' \code{ss3matrix()} returns the stable stage distribution for a 
//' dense or sparse population matrix.
//' 
//' @name .ss3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix.
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ss3matrix)]]
arma::vec ss3matrix(const arma::mat& Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = LefkoMats::decomp3sp(Amat);
  } else {
    eigenstuff = LefkoMats::decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.0000000001);
  
  double rvsum = sum(realrightvec);
  realrightvec = realrightvec / rvsum;
  
  return realrightvec;
}

  //' Estimate Stable Stage Distribution of Any Population Matrix
//' 
//' \code{ss3matrix_sp()} returns the stable stage distribution for a sparse
//' population matrix.
//' 
//' @name ss3matrix_sp
//' 
//' @param Amat A population projection matrix of class \code{dgCMatrix}.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix.
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ss3matrix_sp)]]
arma::vec ss3matrix_sp(const arma::sp_mat& Amat) {
  
  List eigenstuff = LefkoMats::decomp3sp_inp(Amat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.0000000001);
  
  double rvsum = sum(realrightvec);
  realrightvec = realrightvec / rvsum;
  
  return realrightvec;
}

//' Estimate Reproductive Value of Any Population Matrix
//' 
//' \code{rv3matrix()} returns the reproductive values for stages in a
//' dense or sparse population matrix (both provided in dense matrix format).
//' The function provides standard reproductive values, meaning that the overall
//' reproductive values of basic life history stages in a historical matrix are
//' not provided (the \code{\link{repvalue3.lefkoMat}()} function estimates
//' these on the basis of stage description information provided in the
//' \code{lefkoMat} object used as input in that function).
//' 
//' @name .rv3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a vector characterizing the reproductive
//' values for stages of a population projection matrix.
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.rv3matrix)]]
arma::vec rv3matrix(const arma::mat& Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = LefkoMats::decomp3sp(Amat);
  } else {
    eigenstuff = LefkoMats::decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.0000000001);
  
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  
  realleftvec = realleftvec / rlvmin;
  
  return realleftvec;
}

//' Estimate Reproductive Value of Any Population Matrix
//' 
//' \code{rv3matrix_sp()} returns the reproductive values for stages in a
//' sparse population matrix (both provided in dense matrix format).
//' The function provides standard reproductive values, meaning that the overall
//' reproductive values of basic life history stages in a historical matrix are
//' not provided (the \code{\link{repvalue3.lefkoMat}()} function estimates
//' these on the basis of stage description information provided in the
//' \code{lefkoMat} object used as input in that function).
//' 
//' @name .rv3matrix_sp
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' 
//' @return This function returns a vector characterizing the reproductive
//' values for stages of a population projection matrix.
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.rv3matrix_sp)]]
arma::vec rv3matrix_sp(const arma::sp_mat& Amat) {
  List eigenstuff;
  
  eigenstuff = LefkoMats::decomp3sp_inp(Amat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.0000000001);
  
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  
  realleftvec = realleftvec / rlvmin;
  
  return realleftvec;
}

//' Estimate Deterministic Sensitivities of Any Population Matrix
//' 
//' \code{sens3matrix()} returns the sensitivity of lambda with respect
//' to each element in a dense or sparse matrix (provided in dense matrix
//' format). This is accomplished via the \code{eig_gen}() and \code{eigs_gen}()
//' functions in the C++ Armadillo library.
//' 
//' @name .sens3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a matrix of deterministic sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3matrix)]]
arma::mat sens3matrix(const arma::mat& Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = LefkoMats::decomp3sp(Amat);
  } else {
    eigenstuff = LefkoMats::decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;
  
  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel, fill::zeros);
  arma::mat smat (rvel, rvel, fill::zeros);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Populate sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of A Spars Matrixe
//' 
//' \code{sens3sp_matrix()} returns the sensitivity of lambda with respect
//' to each element in a sparse matrix, provided in sparse matrix format. This
//' is accomplished via the \code{eigs_gen}() function in the C++ Armadillo
//' library.
//' 
//' @name .sens3sp_matrix
//' 
//' @param Aspmat A population projection matrix in sparse matrix format.
//' @param refmat A sparse matrix used for reference to create associated 0s in
//' the sensitivity matrix.
//' 
//' @return This function returns a sparse matrix of deterministic
//' sensitivities. Zeroes are derived from the reference matrix, and replace
//' non-zero entries that will be zeroed out in the following math. Currently
//' used in LTRE estimation.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3sp_matrix)]]
arma::sp_mat sens3sp_matrix(const arma::sp_mat& Aspmat, const arma::sp_mat& refmat) {
  
  List eigenstuff = LefkoMats::decomp3sp_inp(Aspmat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;
  
  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel, fill::zeros);
  arma::sp_mat smat (rvel, rvel);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Populate sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      if (refmat(i, j) != 0.0) {
        smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
      }
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of Any Population Matrix
//' 
//' \code{sens3matrix_spinp()} returns the sensitivity of lambda with respect
//' to each element in a sparse matrix. This is accomplished via the
//' \code{eigs_gen()} function in the C++ Armadillo library. Output is in
//' standard matrix format.
//' 
//' @name .sens3matrix_spinp
//' 
//' @param Amat A population projection matrix of class \code{dgCMatrix}.
//' 
//' @return This function returns a standard matrix of deterministic
//' sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3matrix_spinp)]]
arma::mat sens3matrix_spinp(const arma::sp_mat& Amat) {
  List eigenstuff;
  
  eigenstuff = LefkoMats::decomp3sp_inp(Amat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;
  
  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel, fill::zeros);
  arma::mat smat (rvel, rvel, fill::zeros);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Populate sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of a Historical LefkoMat Object
//' 
//' \code{sens3hlefko()} returns the sensitivity of lambda with respect
//' to each historical stage-pair in the matrix, and the associated
//' sensitivity for each life history stage. This is accomplished via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @name .sens3hlefko
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic sensitivity
//' matrices:
//' \item{h_smat}{Matrix of sensitivities corresponding to the historical
//' matrix.}
//' \item{ah_smat}{Matrix of sensitivities corresponding to the ahistorical
//' matrix.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3hlefko)]]
List sens3hlefko(const arma::mat& Amat, const DataFrame& ahstages,
  const DataFrame& hstages) {
  
  arma::uvec stage_id = as<arma::uvec>(ahstages["stage_id"]);
  arma::uvec h_stage_2 = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec h_stage_1 = as<arma::uvec>(hstages["stage_id_1"]);
  
  List eigenstuff = LefkoMats::decomp3sp(Amat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;
  
  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel, fill::zeros);
  arma::mat smat (rvel, rvel, fill::zeros);
  int ahstagelength = static_cast<int>(stage_id.n_elem);
  
  arma::vec wcorrah (ahstagelength, fill::zeros);
  arma::vec vcorrah (ahstagelength, fill::zeros);
  arma::vec vwprodah (ahstagelength, fill::zeros);
  arma::mat ahsens(ahstagelength, ahstagelength, fill::zeros);
  
  // Create calar product vw and ahistorical stable stage dist w
  int ahrows {0};
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
    ahrows = h_stage_2(i) - 1;
    wcorrah(ahrows) = wcorrah(ahrows) + realrightvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Create corrected rep value vector
  for (int i = 0; i < rvel; i++) {
    ahrows = h_stage_2(i) - 1;
    
    if (wcorrah(ahrows) != 0) {
      vcorrah(ahrows) = vwprod(i) / wcorrah(ahrows) + vcorrah(ahrows);
    } else {
      // Some stages have expected corrected stable stage proportions of 0
      vcorrah(ahrows) = 0.0 + vcorrah(ahrows);
    }
  }
  
  // Populate sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  // Create ahistorical sensitivity matrix
  for (int i = 0; i < ahstagelength; i++) {
    vwprodah(i) = wcorrah(i) * vcorrah(i);
  }
  double vwscalarah = sum(vwprodah);
  
  for (int i = 0; i < ahstagelength; i++) {
    for (int j = 0; j < ahstagelength; j++) {
      ahsens(i, j) = vcorrah(i) * wcorrah(j) / vwscalarah;
    }
  }
  
  List output = List::create(Named("h_smat") = smat, _["ah_smat"] = ahsens);
  
  return output;
}

//' Estimate Deterministic Sensitivities of a Historical LefkoMat Object
//' 
//' \code{sens3hlefko()} returns the sensitivity of lambda with respect
//' to each historical stage-pair in the matrix, and the associated
//' sensitivity for each life history stage. This is accomplished via the 
//' \code{eigs_gen}() function in the C++ Armadillo library. Input is in sparse
//' matrix format.
//' 
//' @name .sens3hlefko_sp
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic sensitivity
//' matrices:
//' \item{h_smat}{Matrix of sensitivities corresponding to the historical
//' matrix.}
//' \item{ah_smat}{Matrix of sensitivities corresponding to the ahistorical
//' matrix.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3hlefko_sp)]]
List sens3hlefko_sp(const arma::sp_mat& Amat, const DataFrame& ahstages,
  const DataFrame& hstages) {
  
  arma::uvec stage_id = as<arma::uvec>(ahstages["stage_id"]);
  arma::uvec h_stage_2 = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec h_stage_1 = as<arma::uvec>(hstages["stage_id_1"]);
  
  List eigenstuff = LefkoMats::decomp3sp_inp(Amat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;
  
  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel, fill::zeros);
  arma::mat smat (rvel, rvel, fill::zeros);
  int ahstagelength = static_cast<int>(stage_id.n_elem);
  
  arma::vec wcorrah (ahstagelength, fill::zeros);
  arma::vec vcorrah (ahstagelength, fill::zeros);
  arma::vec vwprodah (ahstagelength, fill::zeros);
  arma::mat ahsens(ahstagelength, ahstagelength, fill::zeros);
  
  // Create calar product vw and ahistorical stable stage dist w
  int ahrows {0};
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
    ahrows = h_stage_2(i) - 1;
    wcorrah(ahrows) = wcorrah(ahrows) + realrightvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Create corrected rep value vector
  for (int i = 0; i < rvel; i++) {
    ahrows = h_stage_2(i) - 1;
    
    if (wcorrah(ahrows) != 0) {
      vcorrah(ahrows) = vwprod(i) / wcorrah(ahrows) + vcorrah(ahrows);
    } else {
      // Some stages have expected corrected stable stage proportions of 0
      vcorrah(ahrows) = 0.0 + vcorrah(ahrows);
    }
  }
  
  // Populate sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  // Create ahistorical sensitivity matrix
  for (int i = 0; i < ahstagelength; i++) {
    vwprodah(i) = wcorrah(i) * vcorrah(i);
  }
  double vwscalarah = sum(vwprodah);
  
  for (int i = 0; i < ahstagelength; i++) {
    for (int j = 0; j < ahstagelength; j++) {
      ahsens(i, j) = vcorrah(i) * wcorrah(j) / vwscalarah;
    }
  }
  
  List output = List::create(Named("h_smat") = smat, _["ah_smat"] = ahsens);
  
  return output;
}

//' Estimate Deterministic Elasticities of Any Population Matrix
//' 
//' \code{elas3matrix()} returns the elasticity of lambda with respect
//' to each element in a dense or sparse matrix, both provided in dense matrix
//' format. This is accomplished via the \code{eig_gen}() and \code{eigs_gen}()
//' functions in the C++ Armadillo library.
//' 
//' @name .elas3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a matrix of deterministic elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3matrix)]]
arma::mat elas3matrix(const arma::mat& Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = LefkoMats::decomp3sp(Amat);
  } else {
    eigenstuff = LefkoMats::decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);

  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // Lower threshold than used in w and v
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;

  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // Lower threshold than used in w and v

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel, fill::zeros);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Populate elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Deterministic Elasticities of Any Population Matrix (Sparse Output)
//' 
//' This function conducts a deterministic elasticity analysis but returns the
//' output in sparse matrix format.
//' 
//' @name .elas3sp_matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' 
//' @return This function returns a matrix of deterministic elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3sp_matrix)]]
arma::sp_mat elas3sp_matrix(const arma::sp_mat& Amat) {
  
  List eigenstuff = LefkoMats::decomp3sp_inp(Amat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // Lower threshold than used in w and v
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;

  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // Lower threshold than used in w and v

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::sp_mat emat (rvel, rvel);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // Populate elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Deterministic Elasticities of a Historical LefkoMat Object
//' 
//' \code{elas3hlefko()} returns the elasticity of lambda with respect
//' to each historical stage-pair in the matrix, and the summed elasticities
//' for each life history stage. This is accomplished via the \code{eigs_gen}()
//' function in the C++ Armadillo library.
//' 
//' @name .elas3hlefko
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic elasticity
//' matrices:
//' \item{h_emat}{Matrix of elasticities corresponding to the historical matrix.}
//' \item{ah_emat}{Matrix of elasticities corresponding to the ahistorical
//' matrix, but using summed historical elasticities as the basis of estimation.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3hlefko)]]
List elas3hlefko(const arma::mat& Amat, const DataFrame& ahstages, const
  DataFrame& hstages) {
  arma::uvec stage_id = as<arma::uvec>(ahstages["stage_id"]);
  arma::uvec h_stage_2 = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec h_stage_1 = as<arma::uvec>(hstages["stage_id_1"]);
  
  List eigenstuff = LefkoMats::decomp3sp(Amat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;

  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel, fill::zeros);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  int ahstagelength = static_cast<int>(stage_id.n_elem);
  arma::mat ahelas(ahstagelength, ahstagelength, fill::zeros);
  
  // Populate elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
      ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) =
        ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) + emat(i, j);
    }
  }
  
  List output = List::create(Named("h_emat") = emat, _["ah_emat"] = ahelas);
  
  return output;
}

//' Estimate Deterministic Elasticities of a Historical LefkoMat Object
//' 
//' \code{elas3hlefko()} returns the elasticity of lambda with respect
//' to each historical stage-pair in the matrix, and the summed elasticities
//' for each life history stage. This is accomplished via the \code{eigs_gen}()
//' function in the C++ Armadillo library.
//' 
//' @name .elas3sp_hlefko
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic elasticity
//' matrices:
//' \item{h_emat}{Matrix of elasticities corresponding to the historical matrix.}
//' \item{ah_emat}{Matrix of elasticities corresponding to the ahistorical
//' matrix, but using summed historical elasticities as the basis of estimation.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3sp_hlefko)]]
List elas3sp_hlefko(const arma::sp_mat& Amat, const DataFrame& ahstages, const
  DataFrame& hstages) {
  arma::uvec stage_id = as<arma::uvec>(ahstages["stage_id"]);
  arma::uvec h_stage_2 = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec h_stage_1 = as<arma::uvec>(hstages["stage_id_1"]);
  
  List eigenstuff = LefkoMats::decomp3sp_inp(Amat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);
  
  // w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = static_cast<int>(realrightvec.n_elem);
  realrightvec = realrightvec / rvsum;

  // v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // Identify first non-zero element of rep value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::sp_mat emat (rvel, rvel);
  
  // Create scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  int ahstagelength = static_cast<int>(stage_id.n_elem);
  arma::sp_mat ahelas(ahstagelength, ahstagelength);
  
  // Populate elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
      ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) =
        ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) + emat(i, j);
    }
  }
  
  List output = List::create(Named("h_emat") = emat, _["ah_emat"] = ahelas);
  
  return output;
}

//' Core Time-based Population Matrix Projection Function
//' 
//' Function \code{proj3()} runs the matrix projections used in other functions
//' in package \code{lefko3}.
//' 
//' @name proj3
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{$A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' @param sparse_auto A logical value indicating whether to determine whether
//' to use sparse matrix encoding automatically.
//' @param sparse A logical value indicating whether to use sparse matrix
//' encoding if \code{sparse_auto = FALSE}.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3)]]
arma::mat proj3(const arma::vec& start_vec, const List& core_list,
  const arma::uvec& mat_order, bool standardize, bool growthonly,
  bool integeronly, bool sparse_auto, bool sparse) {
  
  int sparse_switch {0};
  int nostages = static_cast<int>(start_vec.n_elem);
  int theclairvoyant = static_cast<int>(mat_order.n_elem);
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  arma::mat popproj(nostages, (theclairvoyant + 1), fill::zeros); // Population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1), fill::zeros); // Population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1), fill::zeros); // Population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1), fill::zeros);
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  arma::mat finaloutput;
  
  // Check if matrix is large and sparse
  if (sparse_auto) {
    int test_elems = static_cast<int>(as<arma::mat>(core_list(0)).n_elem);
    arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
    int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
    double sparse_check = static_cast<double>(all_nonzeros) /
      static_cast<double>(test_elems);
    if (sparse_check <= 0.5 && start_vec.n_elem > 50) {
      sparse_switch = 1;
    } else sparse_switch = 0;
  } else sparse_switch = sparse;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  if (sparse_switch == 0) {
    // Dense matrix projection
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
      
      theseventhson = theprophecy * theseventhson;
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      popproj.col(i+1) = theseventhson;
      Rvecmat(i+1) = sum(theseventhson);
      
      if (Rvecmat(i+1) <= 0.0) break;
      
      if (standardize) {
        theseventhson = theseventhson / sum(theseventhson);
      }
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        thesecondprophecy = as<arma::mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * thesecondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
    
  } else {
    // Sparse matrix projection
    arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    
    int matlist_length = static_cast<int>(core_list.size());
    arma::mat first_mat;
    arma::sp_mat new_sparse;
    Rcpp::List sparse_list (matlist_length);
    
    for (int i = 0; i < matlist_length; i++) {
      first_mat = as<arma::mat>(core_list(i));
      new_sparse = arma::sp_mat(first_mat);
      sparse_list(i) = new_sparse;
    }
    arma::sp_mat sparse_prophecy;
    arma::sp_mat sparse_secondprophecy;
    
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      sparse_prophecy = as<arma::sp_mat>(sparse_list[(mat_order(i))]);
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
      Rvecmat(i+1) = sum(popproj.col(i+1));
      
      if (Rvecmat(i+1) <= 0.0) break;
      
      if (standardize) {
        sparse_seventhson = sparse_seventhson / sum(popproj.col(i+1));
      }
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        sparse_secondprophecy = as<arma::sp_mat>(sparse_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Slimmed-down Time-based Population Sparse Matrix Projection Function
//' 
//' Function \code{proj3sp()} runs the matrix projections used in some other
//' functions in package \code{lefko3}, but only when the input is sparse. This
//' is a slimmed down version of function \code{proj3()}
//' 
//' @name proj3sp
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to
//' the \code{$A} list within a \code{lefkoMat} object. Matrices must be in
//' \code{arma::sp_mat} format.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3sp)]]
arma::mat proj3sp(const arma::vec& start_vec, const List& core_list,
  const arma::uvec& mat_order, bool standardize, bool growthonly,
  bool integeronly) {
  
  int nostages = static_cast<int>(start_vec.n_elem);
  int theclairvoyant = static_cast<int>(mat_order.n_elem);
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  arma::mat popproj(nostages, (theclairvoyant + 1), fill::zeros); // Population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1), fill::zeros); // Population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1), fill::zeros); // Population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1), fill::zeros);
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
  arma::mat finaloutput;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  // Sparse matrix projection
  arma::sp_mat sparse_prophecy;
  arma::sp_mat sparse_secondprophecy;
    
  for (int i = 0; i < theclairvoyant; i++) {
    if (i % 50 == 0) Rcpp::checkUserInterrupt();
    
    sparse_prophecy = as<arma::sp_mat>(core_list[(mat_order(i))]);
    sparse_seventhson = sparse_prophecy * sparse_seventhson;
    if (integeronly) {
      sparse_seventhson = floor(sparse_seventhson);
    }
    popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
    Rvecmat(i+1) = sum(popproj.col(i+1));
    
    if (Rvecmat(i+1) <= 0.0) break;
    
    if (standardize) {
      sparse_seventhson = sparse_seventhson / sum(popproj.col(i+1));
    }
    
    if (!growthonly) {
      wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
      sparse_secondprophecy = as<arma::sp_mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
      theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
      
      double seventhgrandsum = sum(theseventhgrandson);
      arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
      theseventhgrandson = theseventhgrandson / seventhgrandsum;
      vpopproj.col(theclairvoyant - (i+1)) = midwife;
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Core Time-based Density-Dependent Population Matrix Projection Function
//' 
//' Function \code{proj3dens()} runs density-dependent matrix projections.
//' 
//' @name proj3dens
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0
//' and fecundity to be non-negative, and \code{2} forces all column rows to
//' total no more than 1.0.
//' @param dens_input The original \code{lefkoDens} data frame supplied through
//' the \code{\link{density_input}()} function.
//' @param dens_index A list giving the indices of elements in object
//' \code{dens_input}.
//' @param sparse_auto A logical value indicating whether to determine whether
//' to use sparse matrix encoding automatically.
//' @param sparse A logical value indicating whether to use sparse matrix
//' encoding if \code{sparse_auto = FALSE}.
//' @param sparse_input A logical value indicating whether matrices in the input
//' MPM are in sparse format (class \code{dgCMatrix}). If so, then all
//' projection will be handled in sparse format. Defaults to \code{FALSE}.
//' @param allow_warnings A logical value indicating whether the function should
//' send warnings if estimated values fall outside of the realm of possibility.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @section Notes:
//' There is no option to standardize population vectors here, because density
//' dependence requires the full population size to be tracked.
//' 
//' @keywords internal
//' @noRd
arma::mat proj3dens(const arma::vec& start_vec, const List& core_list,
  const arma::uvec& mat_order, bool growthonly, bool integeronly, int substoch,
  const Rcpp::DataFrame& dens_input, const Rcpp::List& dens_index,
  bool sparse_auto, bool sparse, bool sparse_input = false,
  bool allow_warnings = false) {
  
  int sparse_switch {0};
  int time_delay {1};
  double pop_size {0};
  bool warn_trigger_neg = false;
  bool warn_trigger_1 = false;
  
  int nostages = static_cast<int>(start_vec.n_elem);
  int theclairvoyant = static_cast<int>(mat_order.n_elem);
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  arma::sp_mat theprophecy_sp;
  arma::sp_mat thesecondprophecy_sp;
  
  // Density dependence
  arma::uvec dyn_index321 = as<arma::uvec>(dens_index["index321"]);
  arma::uvec dyn_index_col = as<arma::uvec>(dens_index[1]);
  arma::uvec dyn_style = as<arma::uvec>(dens_input["style"]);
  arma::vec dyn_alpha = as<arma::vec>(dens_input["alpha"]);
  arma::vec dyn_beta = as<arma::vec>(dens_input["beta"]);
  arma::uvec dyn_delay = as<arma::uvec>(dens_input["time_delay"]);
  arma::uvec dyn_type = as<arma::uvec>(dens_input["type"]);
  int n_dyn_elems = static_cast<int>(dyn_index321.n_elem);
  
  // Matrices and vectors for projection results
  arma::mat popproj(nostages, (theclairvoyant + 1), fill::zeros); // Population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1), fill::zeros); // Population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1), fill::zeros); // Population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1), fill::zeros);
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  arma::mat finaloutput;
  
  // Check if matrix is large and sparse
  if (sparse_auto && !sparse_input) {
    int test_elems = static_cast<int>(as<arma::mat>(core_list(0)).n_elem);
    arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
    int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
    double sparse_check = static_cast<double>(all_nonzeros) /
      static_cast<double>(test_elems);
    if (sparse_check <= 0.5 && theseventhson.n_elem > 50) {
      sparse_switch = 1;
    } else sparse_switch = 0;
  } else sparse_switch = sparse;
  
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  double changing_element {0.0};
  
  if (sparse_switch == 0 && !sparse_input) {
    // Dense matrix projection
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        if (i >= time_delay) {
          pop_size = sum(popproj.col(i - time_delay));
          
          if (dyn_style(j) == 1) { // Ricker
            changing_element = theprophecy(dyn_index321(j)) * 
              dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
            
          } else if (dyn_style(j) == 2) { // Beverton-Holt
            changing_element = theprophecy(dyn_index321(j)) * 
              dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
            
          } else if (dyn_style(j) == 3) { // Usher function
            changing_element = theprophecy(dyn_index321(j)) * 
              (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
            
          } else if (dyn_style(j) == 4) { // Logistic function
            double used_popsize = pop_size;
            if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
              used_popsize = dyn_alpha(j);
            }
            changing_element = theprophecy(dyn_index321(j)) * 
              (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
          }
          
          if (substoch == 1 && dyn_type(j) == 1) {
            if (changing_element > 1.0) {
              changing_element = 1.0;
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch == 2 && dyn_type(j) == 1) {
            arma::vec given_col = theprophecy.col(dyn_index_col(j));
            arma::uvec gc_negs = find(given_col < 0.0);
            arma::uvec gc_pos = find(given_col > 0.0);
            
            double barnyard_antics = sum(given_col(gc_pos)) -
              theprophecy(dyn_index321(j)) + changing_element;
            
            if (barnyard_antics > 1.0 && changing_element > 0.0) {
              double proposed_element = changing_element - barnyard_antics *
                (changing_element / barnyard_antics);
              
              if (proposed_element >= 0.0) {
                changing_element = proposed_element;
              } else {
                changing_element = 0.0;
              }
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch > 0 && dyn_type(j) == 2) {
            if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          }
          theprophecy(dyn_index321(j)) = changing_element;
          
          if (allow_warnings) {
            if (dyn_type(j) == 1 && theprophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              Rf_warningcall(R_NilValue,
                "Some probabilities with value > 1.0 produced during density adjustment.");
              
            } else if (theprophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
              warn_trigger_neg = true;
              Rf_warningcall(R_NilValue,
                "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
      }
      theseventhson = theprophecy * theseventhson;
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      popproj.col(i+1) = theseventhson;
      Rvecmat(i+1) = sum(theseventhson);
      
      if (Rvecmat(i+1) <= 0.0) break;
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        thesecondprophecy = as<arma::mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * thesecondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
    
  } else {
    // Sparse matrix projection
    arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    int matlist_length = static_cast<int>(core_list.size());
    Rcpp::List sparse_list(matlist_length);
    
    if (!sparse_input) {
      arma::mat first_mat = core_list(0);
      arma::sp_mat new_sparse = arma::sp_mat(first_mat);
      sparse_list(0) = new_sparse;
      
      if(matlist_length > 1) {
        for (int i = 1; i < matlist_length; i++) {
          first_mat = as<arma::mat>(core_list(i));
          new_sparse = arma::sp_mat(first_mat);
          sparse_list(i) = new_sparse;
        }
      }
    } else sparse_list = core_list;
    
    arma::sp_mat sparse_prophecy;
    arma::sp_mat sparse_secondprophecy;
    
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      sparse_prophecy = as<arma::sp_mat>(sparse_list[(mat_order(i))]);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        if (i >= time_delay) {
          pop_size = sum(popproj.col(i - time_delay));
          
          if (dyn_style(j) == 1) { // Ricker
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
            
          } else if (dyn_style(j) == 2) { // Beverton-Holt
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
            
          } else if (dyn_style(j) == 3) { // Usher function
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
            
          } else if (dyn_style(j) == 4) { // Logistic function
            double used_popsize = pop_size;
            if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
              used_popsize = dyn_alpha(j);
            }
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
          }
          
          if (substoch == 1 && dyn_type(j) == 1) {
            if (changing_element > 1.0) {
              changing_element = 1.0;
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch == 2 && dyn_type(j) == 1) {
            arma::vec given_col = arma::vec(sparse_prophecy.col(dyn_index_col(j)));
            arma::uvec gc_negs = find(given_col < 0.0);
            arma::uvec gc_pos = find(given_col > 0.0);
            
            double barnyard_antics = sum(given_col(gc_pos)) -
              sparse_prophecy(dyn_index321(j)) + changing_element;
            
            if (barnyard_antics > 1.0 && changing_element > 0.0) {
              double proposed_element = changing_element - barnyard_antics *
                (changing_element / barnyard_antics);
              
              if (proposed_element >= 0.0) {
                changing_element = proposed_element;
              } else {
                changing_element = 0.0;
              }
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch > 0 && dyn_type(j) == 2) {
            if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          }
          sparse_prophecy(dyn_index321(j)) = changing_element;
          
          if (allow_warnings) {
            if (dyn_type(j) == 1 && sparse_prophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              Rf_warningcall(R_NilValue,
                "Some probabilities with value > 1.0 produced during density adjustment.");
                
            } else if (sparse_prophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
              warn_trigger_neg = true;
              Rf_warningcall(R_NilValue,
                "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
      }
      
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
      Rvecmat(i+1) = sum(popproj.col(i+1));
      
      if (Rvecmat(i+1) <= 0.0) break;
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        sparse_secondprophecy = as<arma::sp_mat>(sparse_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Conduct Population Projection Simulations
//' 
//' Function \code{projection3()} runs projection simulations. It projects the
//' population and patches forward in time by a user-defined number of
//' occasions. A given set of matrices is utilized and not recreated, although
//' elements may be altered if density dependence is set. Projections may be
//' deterministic or stochastic, and may be density dependent in either case. If
//' deterministic, then projections will be cyclical if matrices exist covering
//' multiple occasions for each population or patch. If stochastic, then annual
//' matrices will be shuffled within patches and populations. Also produces
//' replicates if set.
//' 
//' @name projection3
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param nreps The number of replicate projections.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' 10,000.
//' @param historical An optional logical value only used if object \code{mpm}
//' is a list of matrices, rather than a \code{lefkoMat} object. Defaults to
//' \code{FALSE} for the former case, and overridden by information supplied in
//' the \code{lefkoMat} object for the latter case.
//' @param stochastic A logical value denoting whether to conduct a stochastic
//' projection or a deterministic / cyclical projection.
//' @param standardize A logical value denoting whether to re-standardize the
//' population size to 1.0 at each occasion. Defaults to \code{FALSE}.
//' @param growthonly A logical value indicating whether to produce only the
//' projected population size at each occasion, or a vector showing the stage
//' distribution followed by the reproductive value vector followed by the full
//' population size at each occasion. Defaults to \code{TRUE}.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion to the nearest
//' integer. Defaults to \code{TRUE}.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0,
//' and forces fecundity to be non-negative; and \code{2} forces all column rows
//' in the survival-transition matrices to total no more than 1.0, in addition
//' to the actions outlined for option \code{1}.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param sub_warnings A logical value indicating whether to warn the user if
//' density dependence yields matrix values outside of the realm of possibility.
//' Generally, this means that survival-transition elements altered to values
//' outside of the interval [0, 1], and negative fecundity values, will both
//' yield warnings. Defaults to \code{TRUE}, but becomes \code{FALSE} if
//' \code{quiet = TRUE}.
//' @param quiet A logical value indicating whether to suppress warnings.
//' Defaults to \code{FALSE}.
//' @param year Either a single integer value corresponding to the year to
//' project, or a vector of \code{times} elements with the year to use at each
//' time step. If a vector shorter than \code{times} is supplied, then this
//' vector will be cycled. If not provided, then all annual matrices will be
//' cycled within patches or populations.
//' @param start_vec An optional numeric vector denoting the starting stage
//' distribution for the projection. Defaults to a single individual of each
//' stage.
//' @param start_frame An optional data frame characterizing stages, age-stages,
//' or stage-pairs that should be set to non-zero values in the starting vector,
//' and what those values should be. Can only be used with \code{lefkoMat}
//' objects.
//' @param tweights An optional numeric vector or matrix denoting the
//' probabilities of choosing each matrix in a stochastic projection. If a
//' matrix is input, then a first-order Markovian environment is assumed, in
//' which the probability of choosing a specific annual matrix depends on which
//' annual matrix is currently chosen. If a vector is input, then the choice of
//' annual matrix is assumed to be independent of the current matrix. Defaults
//' to equal weighting among matrices.
//' @param density An optional data frame describing the matrix elements that
//' will be subject to density dependence, and the exact kind of density
//' dependence that they will be subject to. The data frame used should be an
//' object of class \code{lefkoDens}, which is the output from function
//' \code{\link{density_input}()}.
//' @param sparse A text string indicating whether to use sparse matrix encoding
//' (\code{"yes"}) or dense matrix encoding (\code{"no"}), if the
//' \code{lefkoMat} object input as \code{mpm} is composed of standard matrices.
//' Defaults to \code{"auto"}, in which case sparse matrix encoding is used with
//' standard, square matrices with at least 50 rows and no more than 50\% of
//' elements with values greater than zero, or when input \code{lefkoMat}
//' objects include matrices of class \code{dgCMatrix}.
//' 
//' @return A list of class \code{lefkoProj}, which always includes the first
//' three elements of the following, and also includes the remaining elements
//' below when a \code{lefkoMat} object is used as input:
//' \item{projection}{A list of lists of matrices showing the total number of
//' individuals per stage per occasion. The first list corresponds to each
//' pop-patch followed by each population. The inner list corresponds to
//' replicates within each pop-patch or population.}
//' \item{stage_dist}{A list of lists of the actual stage distribution in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{rep_value}{A list of lists of the actual reproductive value in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{pop_size}{A list of matrices showing the total population size in
//' each occasion per replicate (row within matrix) per pop-patch or
//' population (list element).}
//' \item{labels}{A data frame showing the order of populations and patches in
//' item \code{projection}.}
//' \item{ahstages}{The original stageframe used in the study.}
//' \item{hstages}{A data frame showing the order of historical stage pairs.}
//' \item{agestages}{A data frame showing the order of age-stage pairs.}
//' \item{control}{A short vector indicating the number of replicates and the
//' number of occasions projected per replicate.}
//' \item{density}{The data frame input under the density option. Only provided
//' if input by the user.}
//' 
//' @section Notes:
//' Projections are run both at the patch level and at the population level.
//' Population level estimates will be noted at the end of the data frame with
//' \code{0} entries for patch designation.
//' 
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//' 
//' Starting vectors can be input in one of two ways: 1) as \code{start_vec}
//' input, which is a vector of numbers of the numbers of individuals in each
//' stage, stage pair, or age-stage, with the length of the vector necessarily
//' as long as there are rows in the matrices of the MPM; or 2) as
//' \code{start_frame} input, which is a data frame showing only those stages,
//' stage pairs, or age-stages that should begin with more than 0 individuals,
//' and the numbers of individuals that those stages should start with (this
//' object is created using the \code{\link{start_input}()} function). If both
//' are provided, then \code{start_frame} takes precedence and \code{start_vec}
//' is ignored. If neither is provided, then \code{projection3()} automatically
//' assumes that each stage, stage pair, or age-stage begins with a single
//' individual. Importantly, if a \code{lefkoMat} object is not used, and a list
//' of matrices is provided instead, then \code{start_frame} cannot be utilized
//' and a full \code{start_vec} must be provided to conduct a simulation with
//' starting numbers of individuals other than 1 per stage.
//' 
//' The resulting data frames in element \code{projection} are separated by
//' pop-patch according to the order provided in element \code{labels}, but the
//' matrices for each element of \code{projection} have the result of each
//' replicate stacked in order on top of one another without any break or
//' indication. Results for each replicate must be separated using the
//' information provided in elements \code{control} and the 3 stage
//' descriptor elements.
//' 
//' Density dependent projections are automatically set up if object
//' \code{density} is input. If this object is not included, then density
//' independent projections will be set up. Note that currently, density
//' dependent projections can only be performed with \code{lefkoMat} objects.
//' 
//' When running density dependent simulations involving user-set exponents,
//' such as the beta term in the Ricker function and both the alpha and beta
//' terms in the Usher function, values above or below the computer limits may
//' cause unpredictable behavior. Noted odd behavior includes sudden shifts in
//' population size to negative values. This function produces warnings when
//' such values are used, and the values used for warnings may be reset with the
//' \code{exp_tol} term.
//' 
//' The stage distributions and reproductive values produced are not the
//' asymptotic values as would be given by the standardized right and left
//' eigenvectors associated with the dominant eigenvalue of a matrix, but are
//' vectors describing these values at the specific points in time projected.
//' See equations 14.86 and 14.88 and section 14.4 on Sensitivity and Elasticity
//' Analysis under Environmental Stochasticity in Caswell (2001, Matrix
//' Population Models, Sinauer Associates) for more details.
//' 
//' Consistently positive population growth can quickly lead to population size
//' numbers larger than can be handled computationally. In that circumstance, a
//' continuously rising population size will suddenly become \code{NaN} for the
//' remainder of the projection.
//' 
//' Users wishing to run a projection of a single patch in a \code{lefkoMat}
//' object with multiple patches should subset the MPM first to contain only
//' the patch needed. This can be accomplished with the
//' \code{\link{subset_lM}()} function.
//' 
//' Speed can sometimes be increased by shifting from automatic sparse matrix
//' determination to forced dense or sparse matrix projection. This will most
//' likely occur when matrices have between 30 and 300 rows and columns.
//' Defaults work best when matrices are very small and dense, or very large and
//' sparse. Speed can also be maximized by keeping the default setting,
//' \code{integeronly = TRUE}, since the default behavior is to run each
//' projection (replicate) until either the end, or the population size drops
//' to 0. Setting \code{integeronly = FALSE} may increase runtime dramatically,
//' since the population size can reach extremely small levels without dropping
//' to 0.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{density_input}()}
//' @seealso \code{\link{f_projection3}()}
//' @seealso \code{\link{append_lP}()}
//' @seealso \code{\link{summary.lefkoProj}()}
//' @seealso \code{\link{plot.lefkoProj}()}
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
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
//'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
//'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
//' 
//' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
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
//'   supplement = cypsupp3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- projection3(cypmatrix3r, nreps = 5, stochastic = TRUE)
//' 
//' @export projection3
// [[Rcpp::export(projection3)]]
Rcpp::List projection3(const List& mpm, int nreps = 1, int times = 10000,
  bool historical = false, bool stochastic = false, bool standardize = false,
  bool growthonly = true, bool integeronly = true, int substoch = 0,
  double exp_tol = 700.0, bool sub_warnings = true, bool quiet = false,
  Nullable<IntegerVector> year = R_NilValue,
  Nullable<NumericVector> start_vec = R_NilValue,
  Nullable<DataFrame> start_frame = R_NilValue,
  Nullable<RObject> tweights = R_NilValue,
  Nullable<DataFrame> density = R_NilValue,
  Nullable<RObject> sparse = R_NilValue) {
  
  Rcpp::List dens_index;
  Rcpp::DataFrame dens_input;
  
  int theclairvoyant = times;
  int dens_switch {0};
  int sparse_switch {0};
  bool sparse_bool {false};
  bool sparse_auto {true};
  int used_matsize {0};
  int total_projrows {0};
  bool year_override = false;
  bool sparse_input {false};
  bool matrix_input {true};
  bool assume_markov {false};
  bool agestage_format {false};
  
  if (theclairvoyant < 1) pop_error("times", "a positive integer", "", 1);
  if (nreps < 1) pop_error("nreps", "a positive integer", "", 1);
  if (substoch < 0 || substoch > 2) pop_error("substoch", "integer 0, 1, or 2", "", 1);
  
  if (quiet) sub_warnings = false;
  
  arma::uvec theprophecy(theclairvoyant, fill::zeros);
  
  arma::vec startvec;
  arma::mat projection;
  
  if (sparse.isNotNull()) {
    LefkoInputs::yesnoauto_to_logic(as<RObject>(sparse), "sparse", sparse_bool, sparse_auto);
    if (sparse_bool) sparse_switch = 1;
  }
  
  if (!is<NumericMatrix>(mpm(0)) && !is<S4>(mpm(0))) {
    List amats = as<List>(mpm["A"]);
    List umats = as<List>(mpm["U"]);
    List fmats = as<List>(mpm["F"]);
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    if (is<S4>(amats(0))) {
      sparse_input = true;
      matrix_input = false;
      
    } else if (is<NumericMatrix>(amats(0))) {
      sparse_input = false;
      matrix_input = true;
      
    } else { 
      throw Rcpp::exception("Object mpm does not appear to include matrices.", false);
    }
    
    historical = false;
    
    if (hstages.length() > 1) {
      historical = true;
    }
    
    if (agestages.length() > 1) {
      agestage_format = true;
    } 
    
    if (density.isNotNull()) {
      Rcpp::DataFrame dens_thru(density);
      dens_input = dens_thru;
      dens_switch = 1;
      
      Rcpp::StringVector di_stage3 = as<StringVector>(dens_input["stage3"]);
      Rcpp::StringVector di_stage2 = as<StringVector>(dens_input["stage2"]);
      Rcpp::StringVector di_stage1 = as<StringVector>(dens_input["stage1"]);
      int di_size = di_stage3.length();
      
      if (historical) {
        StringVector stage3 = as<StringVector>(hstages["stage_2"]);
        StringVector stage2r = as<StringVector>(hstages["stage_1"]);
        StringVector stage2c = as<StringVector>(hstages["stage_2"]);
        StringVector stage1 = as<StringVector>(hstages["stage_1"]);
        int hst_size = stage3.length();
        
        arma::uvec hst_3(hst_size, fill::zeros);
        arma::uvec hst_2r(hst_size, fill::zeros);
        arma::uvec hst_2c(hst_size, fill::zeros);
        arma::uvec hst_1(hst_size, fill::zeros);
        
        arma::uvec di_stage32_id(di_size, fill::zeros);
        arma::uvec di_stage21_id(di_size, fill::zeros);
        arma::uvec di_index(di_size, fill::zeros);
        
        for (int i = 0; i < di_size; i++) { // Loop through each density_input line
          for (int j = 0; j < hst_size; j++) {
            if (di_stage3(i) == stage3(j)) {
              hst_3(j) = 1;
            } else {
              hst_3(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage2(i) == stage2r(j)) {
              hst_2r(j) = 1;
            } else {
              hst_2r(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage2(i) == stage2c(j)) {
              hst_2c(j) = 1;
            } else {
              hst_2c(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage1(i) == stage1(j)) {
              hst_1(j) = 1;
            } else {
              hst_1(j) = 0;
            }
          }
          
          arma::uvec find_hst3 = find(hst_3);
          arma::uvec find_hst2r = find(hst_2r);
          arma::uvec find_hst2c = find(hst_2c);
          arma::uvec find_hst1 = find(hst_1);
          
          arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
          arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
          
          di_stage32_id(i) = pop_32(0);
          di_stage21_id(i) = pop_21(0);
          di_index(i) = pop_32(0) + (pop_21(0) * hst_size);
          
          hst_3.zeros();
          hst_2r.zeros();
          hst_2c.zeros();
          hst_1.zeros();
        }
        
        dens_index = Rcpp::List::create(_["index32"] = di_stage32_id,
          _["index21"] = di_stage21_id, _["index321"] = di_index);
        
      } else {
        
        StringVector stage3 = as<StringVector>(stageframe["stage"]);
        StringVector stage2 = as<StringVector>(stageframe["stage"]);
        int ahst_size = stage3.length();
        
        arma::uvec ahst_3(ahst_size, fill::zeros);
        arma::uvec ahst_2(ahst_size, fill::zeros);
        
        arma::uvec di_stage32_id(di_size, fill::zeros);
        arma::uvec di_stage21_id(di_size, fill::zeros);
        arma::uvec di_index(di_size, fill::zeros);
        
        for (int i = 0; i < di_size; i++) { // Loop through each density_input
          for (int j = 0; j < ahst_size; j++) {
            if (di_stage3(i) == stage3(j)) {
              ahst_3(j) = 1;
            } else {
              ahst_3(j) = 0;
            }
          }
          
          for (int j = 0; j < ahst_size; j++) {
            if (di_stage2(i) == stage2(j)) {
              ahst_2(j) = 1;
            } else {
              ahst_2(j) = 0;
            }
          }
          
          arma::uvec find_ahst3 = find(ahst_3);
          arma::uvec find_ahst2 = find(ahst_2);
          di_stage32_id(i) = find_ahst3(0);
          di_stage21_id(i) = find_ahst2(0);
          di_index(i) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
          
          ahst_3.zeros();
          ahst_2.zeros();
        }
        
        dens_index = Rcpp::List::create(_["index3"] = di_stage32_id,
          _["index2"] = di_stage21_id, _["index321"] = di_index);
      }
      
      arma::uvec dyn_style = as<arma::uvec>(dens_input["style"]);
      arma::vec dyn_alpha = as<arma::vec>(dens_input["alpha"]);
      arma::vec dyn_beta = as<arma::vec>(dens_input["beta"]);
      
      for (int i = 0; i < static_cast<int>(dyn_style.n_elem); i++) {
        if (dyn_style(i) < 1 || dyn_style(i) > 4) pop_error("density inputs", "", "", 21);
        
        if (dyn_style(i) == 1) {
          if (dyn_beta(i) > exp_tol) {
            Rf_warningcall(R_NilValue,
              "Beta used in Ricker function may be too high. Results may be unpredictable.");
            
          } else if (dyn_beta(i) < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue,
              "Beta used in Ricker function may be too high. Results may be unpredictable.");
            
          }
          
        } else if (dyn_style(i) == 3) {
          double summed_stuff = dyn_alpha(i) + dyn_beta(i);
          
          if (summed_stuff > exp_tol) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high. Results may be unpredictable.");
            
          } else if (summed_stuff < (-1.0 * exp_tol)) {
            Rf_warningcall(R_NilValue,
              "Alpha and beta used in Usher function may be too high. Results may be unpredictable.");
          }
        }
      }
    }
    
    StringVector yearorder;
    StringVector patchorder;
    if (labels.length() < 3) {
      StringVector label_elements = labels.attr("names");
      std::string patch_named = "patch";
      
      for (int i = 0; i < label_elements.length(); i++) {
        if (stringcompare_hard(as<std::string>(label_elements(i)), "patch")) {
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "This lefkoMat object appears to be a set of mean matrices. Will project only the mean.");
          }
        }
      }
      
      StringVector patch_projected = as<StringVector>(labels["patch"]);
      StringVector years_projected(patch_projected.length());
      for (int i = 0; i < patch_projected.length(); i++) {
        years_projected(i) = "1";
      }
      
      patchorder = patch_projected;
      yearorder = years_projected;
    } else {
      patchorder = as<StringVector>(labels["patch"]);
      yearorder = as<StringVector>(labels["year2"]);
    }
    StringVector poporder = as<StringVector>(labels["pop"]);
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    StringVector uniqueyears = sort_unique(yearorder);
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    StringVector years_forward;
    if (year.isNotNull()) {
      if (stochastic) {
        throw Rcpp::exception("Argument year cannot be used when stochastic = TRUE.",
          false);
      }
      
      StringVector years_ = as<StringVector>(year);
      
      int member_sum {0};
      for (int i = 0; i < years_.length(); i++) {
        for (int j = 0; j < yl; j++) {
          if (years_[i] == uniqueyears[j]) member_sum++;
        }
        if (member_sum == 0) {
          pop_error("year", "input lefkoMat object", "", 17);
        }
        member_sum = 0;
      }
      
      StringVector years_pre (times);
      
      int rampant_exigence {0};
      for (int i = 0; i < times; i++) {
        years_pre(i) = years_(rampant_exigence);
        rampant_exigence++;
        
        if (rampant_exigence >= years_.length()) {
          rampant_exigence = 0;
        }
      }
      years_forward = years_pre; // Order of matrices for all times, if years is input
      year_override = true; // Use years or default matrix vectors
    }
    
    arma::vec twinput;
    arma::mat twinput_markov;
    if (tweights.isNotNull()) {
      if (Rf_isMatrix(tweights)) {
        twinput_markov = as<arma::mat>(tweights);
        assume_markov = true;
        
        if (static_cast<int>(twinput_markov.n_cols) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
        }
        if (twinput_markov.n_cols != twinput_markov.n_rows) {
          throw Rcpp::exception("Time weight matrix must be square.", false);
        }
        
      } else if (is<NumericVector>(tweights)) {
        twinput = as<arma::vec>(tweights);
        
        if (static_cast<int>(twinput.n_elem) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
        }
        
      } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
      
      if (!stochastic) {
        throw Rcpp::exception("Argument tweights can only be used when stochastic = TRUE.",
          false);
      }
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize);
    arma::uvec yearsinpatch(loysize);
    patchesinpop.zeros();
    yearsinpatch.zeros();
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = static_cast<int>(summervacation.n_elem);
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = static_cast<int>(motorhead.n_elem);
      
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    List mean_lefkomat;
    
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1, matrix_input, sparse_switch);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1, matrix_input, sparse_switch);
    }
    
    // Run simulation on all patch matrices, estimate all descriptive metrics
    List meanamats = as<List>(mean_lefkomat["A"]);
    List mmlabels = as<List>(mean_lefkomat["labels"]);
    StringVector mmpops_1 = as<StringVector>(mmlabels["pop"]);
    StringVector mmpatches_1 = as<StringVector>(mmlabels["patch"]);
    
    int meanmatsize {0};
    int meanmatrows {0};
    
    if (!sparse_input) {
      arma::mat thechosenone = as<arma::mat>(meanamats(0));
      meanmatsize = static_cast<int>(thechosenone.n_elem);
      meanmatrows = static_cast<int>(thechosenone.n_rows);
    } else {
      arma::sp_mat thechosenone = as<arma::sp_mat>(meanamats(0));
      meanmatsize = static_cast<int>(thechosenone.n_elem);
      meanmatrows = static_cast<int>(thechosenone.n_rows);
    }
    
    arma::vec startvec;
    int trials = meanamats.length();
    used_matsize = meanmatrows;
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = static_cast<int>(allppcs.n_elem);
    
    StringVector mmpops;
    StringVector mmpatches;
    
    bool add_mean {true};
    int mmpatches_1_length = mmpatches_1.length();
    if (mmpatches_1(mmpatches_1_length - 1) == "0" && mmpatches_1_length > 1) {
      if (mmpatches_1(mmpatches_1_length - 2) == "0") {
        add_mean = false;
      }
    }
    
    if (!add_mean) {
      trials = allppcsnem;
      
      StringVector new_mmpops (trials);
      StringVector new_mmpatches (trials);
      
      for (int i = 0; i < trials; i++) {
        new_mmpops(i) = mmpops_1(i);
        new_mmpatches(i) = mmpatches_1(i);
      }
      
      mmpops = new_mmpops;
      mmpatches = new_mmpatches;
    } else {
      mmpops = mmpops_1;
      mmpatches = mmpatches_1;
    }
    
    List plist_hold(allppcsnem);
    int pop_est {1};
    if (allppcsnem > 1) {
      pop_est = trials - allppcsnem;    
    }
    
    List projection_list(trials);
    
    if(start_frame.isNotNull()) {
      Rcpp::DataFrame start_thru(start_frame);
      startvec.set_size(meanmatrows);
      startvec.zeros();
      
      arma::uvec start_elems = as<arma::uvec>(start_thru["row_num"]);
      start_elems = start_elems - 1;
      arma::vec start_values = as<arma::vec>(start_thru["value"]);
      
      // Here we will check to make sure the row designations are consistent
      // with the stageframe, and if not, fix them
      StringVector startframe_stage2 = as<StringVector>(start_thru["stage2"]);
      StringVector startframe_stage1 = as<StringVector>(start_thru["stage1"]);
      IntegerVector startframe_age2 = as<IntegerVector>(start_thru["age2"]);
      int startframe_rows = startframe_stage2.length();
      
      if (agestage_format) {
        StringVector agestages_stage = as<StringVector>(agestages["stage"]);
        IntegerVector agestages_age = as<IntegerVector>(agestages["age"]);
        int agestages_rows = agestages_stage.length();
        
        for (int i = 0; i < startframe_rows; i++) { 
          for (int j = 0; j < agestages_rows; j++) {
            for (int k = 0; k < agestages_rows; k++) {
              if (LefkoUtils::stringcompare_simple(String(agestages_stage(j)),
                    String(startframe_stage2(i)))) {
                if (agestages_age(k) == startframe_age2(i)) {
                  start_elems(i) = k;
                }
              }
            }
          }
        }
      } else if (historical) {
        StringVector hstages_stage2 = as<StringVector>(hstages["stage_2"]);
        StringVector hstages_stage1 = as<StringVector>(hstages["stage_1"]);
        int hstages_rows = hstages_stage2.length();
        
        for (int i = 0; i < startframe_rows; i++) { 
          for (int j = 0; j < hstages_rows; j++) {
            for (int k = 0; k < hstages_rows; k++) {
              if (LefkoUtils::stringcompare_simple(String(hstages_stage2(j)),
                    String(startframe_stage2(i)))) {
                if (LefkoUtils::stringcompare_simple(String(hstages_stage1(k)),
                    String(startframe_stage1(i)))) {
                  start_elems(i) = k;
                }
              }
            }
          }
        }
      } else {
        StringVector stageframe_stage2 = as<StringVector>(stageframe["stage"]);
        int stageframe_rows = stageframe_stage2.length();
        
        for (int i = 0; i < startframe_rows; i++) { 
          for (int j = 0; j < stageframe_rows; j++) {
            if (LefkoUtils::stringcompare_simple(String(stageframe_stage2(j)),
                  String(startframe_stage2(i)))) {
              start_elems(i) = j;
            }
          }
        }
      }
      
      if (static_cast<int>(start_elems.max()) > (meanmatrows - 1)) {
        throw Rcpp::exception("Start vector input frame includes element indices too high for this MPM.",
          false);
      }
      for (int i = 0; i < static_cast<int>(start_elems.n_elem); i++) {
        startvec(start_elems(i)) = start_values(i);
      }
      
    } else if (start_vec.isNotNull()) {
      startvec = as<arma::vec>(start_vec);
      if (static_cast<int>(startvec.n_elem) != meanmatrows) {
        pop_error("start_vec", "rows", "MPM matrices", 19);
      }
      
    } else {
      startvec.set_size(meanmatrows);
      startvec.ones();
    }
    
    if (!assume_markov) twinput = twinput / sum(twinput);
    
    for (int i= 0; i < allppcsnem; i++) {
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      int chosen_yl = static_cast<int>(thenumbersofthebeast.n_elem);
      
      arma::uvec pre_prophecy (theclairvoyant, fill::zeros);
      if (year_override) {
        for (int j = 0; j < theclairvoyant; j++) {
          IntegerVector tnb_year_indices_IV = (match(as<StringVector>(years_forward(j)),
            yearorder) - 1) + (i * yl);
          arma::uvec tnb_year_indices = as<arma::uvec>(tnb_year_indices_IV);
          arma::uvec year_patch_intersect = intersect(thenumbersofthebeast,
            tnb_year_indices);
          
          pre_prophecy(j) = year_patch_intersect(0);
        }
      }
      
      // Replicate loop, creating final data frame of results for each pop-patch
      for (int rep = 0; rep < nreps; rep++) {
        if (stochastic && !assume_markov) {
          theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
            theclairvoyant, true, twinput);
          
        } else if (stochastic && assume_markov) {
          for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
            if (yr_counter == 0) {
              twinput = twinput_markov.col(0);
            }
            twinput = twinput / sum(twinput);
            
            arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
              1, true, twinput);
            theprophecy(yr_counter) = theprophecy_piecemeal(0);
              
            arma::uvec tnotb_preassigned = find(thenumbersofthebeast == theprophecy_piecemeal(0));
            twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
          }
        } else if (year_override) {
          theprophecy = pre_prophecy;
          
        } else {
          theprophecy.set_size(theclairvoyant);
          theprophecy.zeros();
          
          for (int j = 0; j < theclairvoyant; j++) {
            theprophecy(j) = thenumbersofthebeast(j % chosen_yl);
          }
        }
        
        if (dens_switch) {
          if (rep == 0) {
            projection = proj3dens(startvec, amats, theprophecy, growthonly,
              integeronly, substoch, dens_input, dens_index, sparse_auto,
              sparse_switch, sparse_input, sub_warnings);
          } else {
            arma::mat nextproj = proj3dens(startvec, amats, theprophecy,
              growthonly, integeronly, substoch, dens_input, dens_index,
              sparse_auto, sparse_switch, sparse_input, sub_warnings);
            projection = arma::join_cols(projection, nextproj);
          }
        } else {
          if (rep == 0) {
            if (!sparse_input) {
              projection = proj3(startvec, amats, theprophecy, standardize,
                growthonly, integeronly, sparse_auto, sparse_switch);
            } else {
              projection = proj3sp(startvec, amats, theprophecy, standardize,
                growthonly, integeronly);
            }
          } else {
            if (!sparse_input) {
              arma::mat nextproj = proj3(startvec, amats, theprophecy,
                standardize, growthonly, integeronly, sparse_auto, sparse_switch);
              projection = arma::join_cols(projection, nextproj);
            } else {
              arma::mat nextproj = proj3sp(startvec, amats, theprophecy,
                standardize, growthonly, integeronly);
              projection = arma::join_cols(projection, nextproj);
            }
          }
        }
      }
      
      projection_list(i) = projection;
    }
    
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize, fill::zeros);
    arma::uvec yearmatch(loysize, fill::zeros);
    List meanmatyearlist(uniqueyears.length());
    
    if (allppcsnem > 1) { // Checks pop-mean matrices separate from patch means
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // population loop
        for (int j = 0; j < loysize; j++) { // Checks which A matrices match current pop
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        arma::uvec neededmatspop = find(popmatch);
        
        for (int j = 0; j < yl; j++) { // Develops matrix mean across patches
          for (int k = 0; k < loysize; k++) { // Develops vector to find all mats for current year
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          
          // Catches matrix indices matching current year and pop
          int crankybankynem = static_cast<int>(crankybanky.n_elem);
          
          if (!sparse_input) {
            arma::mat crossmat(meanmatsize, crankybankynem, fill::zeros);
            for (int j = 0; j < crankybankynem; j++) {
              crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
            }
            
            arma::mat happymedium(meanmatsize, 1, fill::zeros);
            for (int j = 0; j < meanmatsize; j++) {
              for (int k = 0; k < crankybankynem; k++) {
                happymedium(j, 0) = happymedium(j, 0) + crossmat(j, k) / crankybankynem;
              }
            }
            happymedium.reshape(meanmatrows, meanmatrows);
            meanmatyearlist(j) = happymedium;
            
          } else {
            arma::sp_mat crossmat(meanmatsize, crankybankynem);
            for (int j = 0; j < crankybankynem; j++) {
              crossmat.col(j) = vectorise(as<arma::sp_mat>(amats(crankybanky(j))));
            }
            
            arma::sp_mat happymedium(meanmatsize, 1);
            for (int j = 0; j < meanmatsize; j++) {
              for (int k = 0; k < crankybankynem; k++) {
                happymedium(j, 0) = happymedium(j, 0) + crossmat(j, k) / crankybankynem;
              }
            }
            happymedium.reshape(meanmatrows, meanmatrows);
            meanmatyearlist(j) = happymedium;
          }
        }
        
        int numyearsused = meanmatyearlist.length();
        arma::uvec choicevec = linspace<arma::uvec>(0, (numyearsused - 1),
          numyearsused);
        int chosen_yl = static_cast<int>(choicevec.n_elem);
        
        // Replicate loop, creating final data frame of results for pop means
        for (int rep = 0; rep < nreps; rep++) {
          if (stochastic && !assume_markov) {
            theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant,
              true, twinput);
          } else if (stochastic && assume_markov) {
            for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
              if (yr_counter == 0) {
                twinput = twinput_markov.col(0);
              }
              twinput = twinput / sum(twinput);
              
              arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(choicevec,
                1, true, twinput);
              theprophecy(yr_counter) = theprophecy_piecemeal(0);
                
              arma::uvec tnotb_preassigned = find(choicevec == theprophecy_piecemeal(0));
              twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
            }
          } else {
            theprophecy.zeros();
            for (int j = 0; j < theclairvoyant; j++) {
              theprophecy(j) = choicevec(j % chosen_yl);
            }
          }
          
          if (dens_switch) {
            if (rep == 0) {
              projection = proj3dens(startvec, meanmatyearlist, theprophecy,
                growthonly, integeronly, substoch, dens_input, dens_index,
                sparse_auto, sparse_switch, sparse_input, sub_warnings);
            } else {
              arma::mat nextproj = proj3dens(startvec, meanmatyearlist,
                theprophecy, growthonly, integeronly, substoch, dens_input,
                dens_index, sparse_auto, sparse_switch, sparse_input,
                sub_warnings);
              projection = arma::join_cols(projection, nextproj);
            }
          } else {
            if (rep == 0) {
              if (!sparse_input) {
                projection = proj3(startvec, meanmatyearlist, theprophecy,
                  standardize, growthonly, integeronly, sparse_auto, sparse_switch);
              } else {
                projection = proj3sp(startvec, meanmatyearlist, theprophecy,
                  standardize, growthonly, integeronly);
              }
            } else {
              if (!sparse_input) {
                arma::mat nextproj = proj3(startvec, meanmatyearlist, theprophecy,
                  standardize, growthonly, integeronly, sparse_auto, sparse_switch);
                projection = arma::join_cols(projection, nextproj);
              } else {
                arma::mat nextproj = proj3sp(startvec, meanmatyearlist, theprophecy,
                  standardize, growthonly, integeronly);
                projection = arma::join_cols(projection, nextproj);
              }
            }
          }
        }
        projection_list(allppcsnem + i) = projection;
      }
    }
    
    // Output proj list w/ #elem = nreps nested within list w/ #elem = #poppatches
    List projection_set(nreps);
    List ss_set(nreps);
    List rv_set(nreps);
    arma::mat total_sizes_set(nreps, (times+1), fill::zeros);
    
    int length_ppy = projection_list.length();
    List final_projection(length_ppy);
    List final_ss(length_ppy);
    List final_rv(length_ppy);
    List final_ns(length_ppy);
    
    int out_elements {9};
    if (dens_switch) out_elements++;
    
    List output (out_elements);
    
    if (!growthonly) {
      arma::mat list_proj(total_projrows, (times+1), fill::zeros);
      arma::mat extracted_proj(used_matsize, used_matsize, fill::zeros);
      int diversion = used_matsize * 3 + 1;
      
      for (int j = 0; j < length_ppy; j++) {
        list_proj = as<arma::mat>(projection_list[j]);
        
        for (int i = 0; i < nreps; i++) {
          extracted_proj = list_proj.rows((diversion * i),
            (diversion * i + used_matsize - 1));
          projection_set(i) = extracted_proj;
          
          extracted_proj = list_proj.rows((diversion * i + used_matsize),
            (diversion * i + (2 * used_matsize) - 1));
          ss_set(i) = extracted_proj;
          
          extracted_proj = list_proj.rows((diversion * i + (2 * used_matsize)),
            (diversion * i + (3 * used_matsize) - 1));
          rv_set(i) = extracted_proj;
          
          total_sizes_set.row(i) = list_proj.row(diversion * (i+1) - 1);
        }
        final_projection(j) = clone(projection_set);
        final_ss(j) = clone(ss_set);
        final_rv(j) = clone(rv_set);
        final_ns(j) = total_sizes_set;
      }
      
    } else {
      arma::mat list_proj(total_projrows, (times+1), fill::zeros);
      arma::mat extracted_proj(used_matsize, used_matsize, fill::zeros);
      int diversion = used_matsize;
      
      for (int j = 0; j < length_ppy; j++) {
        list_proj = as<arma::mat>(projection_list[j]);
        
        for (int i = 0; i < nreps; i++) {
          extracted_proj = list_proj.rows((diversion * i),
            (diversion * i + used_matsize - 1));
          projection_set(i) = extracted_proj;
          
          total_sizes_set.row(i) = sum(extracted_proj, 0);
        }
        final_projection(j) = clone(projection_set);
        final_ss(j) = NULL;
        final_rv(j) = NULL;
        final_ns(j) = total_sizes_set;
      }
    }
    
    DataFrame newlabels = DataFrame::create(_["pop"] = mmpops,
      _["patch"] = mmpatches);
    Rcpp::IntegerVector control = {nreps, times};
    
    output(0) = final_projection;
    output(1) = final_ss;
    output(2) = final_rv;
    output(3) = final_ns;
    output(4) = newlabels;
    output(5) = stageframe;
    output(6) = hstages;
    output(7) = agestages;
    output(8) = control;
    
    if (dens_switch) {
      output(9) = dens_input;
      
      CharacterVector namevec = {"projection", "stage_dist", "rep_value", "pop_size",
        "labels", "ahstages", "hstages", "agestages", "control", "density"};
      output.attr("names") = namevec;
    } else {
      CharacterVector namevec = {"projection", "stage_dist", "rep_value", "pop_size",
        "labels", "ahstages", "hstages", "agestages", "control"};
      output.attr("names") = namevec;
    }
    output.attr("class") = "lefkoProj";
    
    return output;
    
  } else {
    // Matrix list input
    
    List projection_list (1);
    List amats = mpm;
    int yl = amats.length();
    
    int matrows {0};
    int matcols {0};
    
    if (is<S4>(amats(0))) {
      sparse_input = true;
      
      arma::sp_mat firstmat = as<arma::sp_mat>(amats(0));
      matrows = firstmat.n_rows;
      matcols = firstmat.n_cols;
      
    } else if (is<NumericMatrix>(amats(0))) {
      arma::mat firstmat = as<arma::mat>(amats(0));
      matrows = firstmat.n_rows;
      matcols = firstmat.n_cols;
      
    } else throw Rcpp::exception("Object mpm does not appear to include matrices.", false);
    
    used_matsize = matrows;
    
    arma::uvec uniqueyears(yl);
    for (int i = 0; i < yl; i++) {
      uniqueyears(i) = i;
    }
    
    if (matrows != matcols) throw Rcpp::exception("Supplied matrices must be square.", false);
    
    arma::vec twinput;
    arma::mat twinput_markov;
    if (tweights.isNotNull()) {
      if (Rf_isMatrix(tweights)) {
        twinput_markov = as<arma::mat>(tweights);
        assume_markov = true;
        
        if (static_cast<int>(twinput_markov.n_cols) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
        }
        if (twinput_markov.n_cols != twinput_markov.n_rows) {
          throw Rcpp::exception("Time weight matrix must be square.", false);
        }
        
      } else if (is<NumericVector>(tweights)) {
        twinput = as<arma::vec>(tweights);
        
        if (static_cast<int>(twinput.n_elem) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
        }
        
      } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
      
      if (!stochastic) {
        throw Rcpp::exception("Argument tweights can only be used when stochastic = TRUE.", false);
      }
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    if (start_vec.isNotNull()) {
      startvec = as<arma::vec>(start_vec);
      if (static_cast<int>(startvec.n_elem) != matrows) {
        pop_error("start_vec", "rows", "MPM matrices", 19);
      }
      
    } else {
      startvec.set_size(matrows);
      startvec.ones();
    }
    
    IntegerVector years_forward;
    if (year.isNotNull()) {
      if (stochastic) {
        throw Rcpp::exception("Argument year cannot be used when stochastic = TRUE.",
          false);
      }
      
      IntegerVector years_ = clone(as<IntegerVector>(year));
      years_ = years_ - 1;
      
      int member_sum {0};
      for (int i = 0; i < years_.length(); i++) {
        for (int j = 0; j < yl; j++) {
          if (years_[i] == static_cast<int>(uniqueyears[j])) member_sum++;
        }
        if (member_sum == 0) pop_error("year", "lefkoMat object", "", 17);
        member_sum = 0;
      }
      
      IntegerVector years_pre (times);
      
      int rampant_exigence {0};
      for (int i = 0; i < times; i++) {
        years_pre(i) = years_(rampant_exigence);
        rampant_exigence++;
        
        if (rampant_exigence >= years_.length()) {
          rampant_exigence = 0;
        }
      }
      years_forward = years_pre; // Order of matrices for all times, if years input
      year_override = true; // Decides whether to use years or default matrix vectors
    }
    
    // Run simulation, estimate descriptive metrics
    if (!assume_markov) twinput = twinput / sum(twinput);
    arma::uvec thenumbersofthebeast = uniqueyears;
    
    // Replicate loop, creating a data frame of results
    for (int rep = 0; rep < nreps; rep++) {
      if (stochastic && !assume_markov) {
        theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
          theclairvoyant, true, twinput);
        
      } else if (stochastic && assume_markov) {
        for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
          if (yr_counter == 0) {
            twinput = twinput_markov.col(0);
          }
          twinput = twinput / sum(twinput);
          
          arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
            1, true, twinput);
          theprophecy(yr_counter) = theprophecy_piecemeal(0);
            
          arma::uvec tnotb_preassigned = find(thenumbersofthebeast == theprophecy_piecemeal(0));
          twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
        }
      } else if (year_override) {
        theprophecy = as<arma::uvec>(years_forward);
        
      } else {
        theprophecy.zeros();
        for (int i = 0; i < theclairvoyant; i++) {
          theprophecy(i) = thenumbersofthebeast(i % yl);
        }
      }
      
      if (!sparse_input) {
        if (rep == 0) {
          projection = proj3(startvec, amats, theprophecy, standardize,
            growthonly, integeronly, sparse_auto, sparse_switch);
          
        } else {
          arma::mat nextproj = proj3(startvec, amats, theprophecy, standardize,
            growthonly, integeronly, sparse_auto, sparse_switch);
          projection = arma::join_cols(projection, nextproj);
        }
      } else {
        if (rep == 0) {
          projection = proj3sp(startvec, amats, theprophecy, standardize,
            growthonly, integeronly);
          
        } else {
          arma::mat nextproj = proj3sp(startvec, amats, theprophecy,
            standardize, growthonly, integeronly);
          projection = arma::join_cols(projection, nextproj);
        }
      }
    }
    projection_list(0) = projection;
    DataFrame newlabels = DataFrame::create(_["pop"] = 1, _["patch"] = 1);
    Rcpp::IntegerVector control = {nreps, times};
    
    Rcpp::List output = List::create(_["projection"] = projection_list,
      _["labels"] = newlabels, _["control"] = control);
    output.attr("class") = "lefkoProj";
    
    return output;
  }
}

//' Estimate Stochastic Population Growth Rate
//' 
//' Function \code{slambda3()} estimates the stochastic population growth rate,
//' \eqn{a}, defined as the long-term arithmetic mean of the log population 
//' growth rate estimated per simulated occasion. This function can handle both
//' lefkoMat objects and lists of full A matrices as input. 
//' 
//' @name slambda3
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param times Number of occasions to iterate. Defaults to \code{10000}.
//' @param historical An optional logical value only used if object \code{mpm}
//' is a list of matrices, rather than a \code{lefkoMat} object. Defaults to
//' \code{FALSE} for the former case, and overridden by information supplied in
//' the \code{lefkoMat} object for the latter case.
//' @param tweights An optional numeric vector or matrix denoting the
//' probabilities of choosing each matrix in a stochastic projection. If a
//' matrix is input, then a first-order Markovian environment is assumed, in
//' which the probability of choosing a specific annual matrix depends on which
//' annual matrix is currently chosen. If a vector is input, then the choice of
//' annual matrix is assumed to be independent of the current matrix. Defaults
//' to equal weighting among matrices.
//' @param force_sparse A text string indicating whether to force sparse matrix 
//' encoding (\code{"yes"}) or not (\code{"no"}) if the MPM is composed of
//' simple matrices. Defaults to \code{"auto"}, in which case sparse matrix
//' encoding is used with simple square matrices with at least 50 rows and no
//' more than 50\% of elements with values greater than zero.
//' 
//' @return A data frame with the following variables:
//' \item{pop}{The identity of the population.}
//' \item{patch}{The identity of the patch.}
//' \item{a}{Estimate of stochastic growth rate, estimated as the arithmetic
//' mean of the log population growth rate across simulated occasions.}
//' \item{var}{The estimated variance of a.}
//' \item{sd}{The standard deviation of a.}
//' \item{se}{The standard error of a.}
//' 
//' @section Notes:
//' The log stochastic population growth rate, \eqn{a}, is as given in equation
//' 2 of Tuljapurkar, Horvitz, and Pascarella 2003. This term is estimated via
//' projection of randomly sampled matrices, similarly to the procedure outlined
//' in Box 7.4 of Morris and Doak (2002).
//' 
//' Stochastic growth rate is estimated both at the patch level and at the
//' population level. Population level estimates will be noted at the end of the
//' data frame with 0 entries for patch designation.
//' 
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//' 
//' Speed can sometimes be increased by shifting from automatic sparse matrix
//' determination to forced dense or sparse matrix projection. This will most
//' likely occur when matrices have between 30 and 300 rows and columns.
//' Defaults work best when matrices are very small and dense, or very large and
//' sparse.
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
//'   supplement = cypsupp3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- slambda3(cypmatrix3r)
//' 
//' @export slambda3
// [[Rcpp::export(slambda3)]]
DataFrame slambda3(const List& mpm, int times = 10000, bool historical = false,
  Nullable<RObject> tweights = R_NilValue,
  Nullable<RObject> force_sparse = R_NilValue) {
  
  int theclairvoyant {0};
  int sparse_switch {0};
  bool sparse_bool {false};
  bool sparse_auto {true};
  bool assume_markov {false};
  
  theclairvoyant = times;
  
  if (theclairvoyant < 1) pop_error("times", "a positive integer", "", 1);
  
  if (force_sparse.isNotNull()) {
    LefkoInputs::yesnoauto_to_logic(as<RObject>(force_sparse), "force_sparse", sparse_bool,
      sparse_auto);
    if (sparse_bool) sparse_switch = 1;
  }
  
  bool matrix_class_input {false};
  
  if (is<List>(mpm(0))) {
    List amats = as<List>(mpm["A"]);
    List umats = as<List>(mpm["U"]);
    List fmats = as<List>(mpm["F"]);
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    if (hstages.length() > 1) {
      historical = true;
    } else {
      historical = false;
    }
    
    if (is<NumericMatrix>(amats(0))) matrix_class_input = true;
    
    // Check if matrix is large and sparse
    if (sparse_auto && matrix_class_input) {
      int test_elems = static_cast<int>(as<arma::mat>(amats(0)).n_elem);
      int Amatrows = as<arma::mat>(amats(0)).n_rows;
      arma::uvec nonzero_elems = find(as<arma::mat>(amats(0)));
      int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
      double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
      if (sparse_check <= 0.5 && Amatrows > 50) {
        sparse_switch = 1;
      } else {
        sparse_switch = 0;
      }
    }
    
    if (labels.length() < 3) {
      throw Rcpp::exception("This function requires annual matrices, and cannot take mean matrices.", false);
    }
    
    StringVector poporder = as<StringVector>(labels["pop"]);
    StringVector patchorder = as<StringVector>(labels["patch"]);
    StringVector yearorder = as<StringVector>(labels["year2"]);
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    StringVector uniqueyears = sort_unique(yearorder);
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    arma::vec twinput;
    arma::mat twinput_markov;
    if (tweights.isNotNull()) {
      if (Rf_isMatrix(tweights)) {
        twinput_markov = as<arma::mat>(tweights);
        assume_markov = true;
        
        if (static_cast<int>(twinput_markov.n_cols) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
        }
        if (twinput_markov.n_cols != twinput_markov.n_rows) {
          throw Rcpp::exception("Time weight matrix must be square.", false);
        }
        
      } else if (is<NumericVector>(tweights)) {
        twinput = as<arma::vec>(tweights);
        
        if (static_cast<int>(twinput.n_elem) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
        }
        
      } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize, fill::zeros);
    arma::uvec yearsinpatch(loysize, fill::zeros);
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = static_cast<int>(summervacation.n_elem);
      
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = static_cast<int>(motorhead.n_elem);
      
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    List mean_lefkomat;
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1, matrix_class_input, sparse_switch);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1, matrix_class_input, sparse_switch);
    }
    
    // Run simulation for each patch, estimate metrics
    List meanamats = as<List>(mean_lefkomat["A"]);
    List mmlabels = as<List>(mean_lefkomat["labels"]);
    StringVector mmpops = as<StringVector>(mmlabels["pop"]);
    StringVector mmpatches = as<StringVector>(mmlabels["patch"]);
    
    arma::vec startvec;
    int trials = meanamats.length();
    int meanmatsize {0};
    int meanmatrows {0};
    
    if (matrix_class_input && sparse_switch == 0) {
      meanmatsize = static_cast<int>(as<arma::mat>(meanamats[0]).n_elem); // thechosenone
      meanmatrows = static_cast<int>(as<arma::mat>(meanamats[0]).n_rows);
    } else {
      meanmatsize = static_cast<int>(as<arma::sp_mat>(meanamats[0]).n_elem);
      meanmatrows = static_cast<int>(as<arma::sp_mat>(meanamats[0]).n_rows);
    }
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = static_cast<int>(allppcs.n_elem);
    
    arma::mat slmat(theclairvoyant, trials, fill::zeros);
    arma::vec sl_mean(trials, fill::zeros);
    arma::vec sl_var(trials, fill::zeros);
    arma::vec sl_sd(trials, fill::zeros);
    arma::vec sl_se(trials, fill::zeros);
    
    arma::uvec theprophecy (theclairvoyant);
    for (int i= 0; i < allppcsnem; i++) {
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      if (!assume_markov) {
        twinput = twinput / sum(twinput);
        
        theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
          theclairvoyant, true, twinput);
          
      } else {
        for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
          if (yr_counter == 0) {
            twinput = twinput_markov.col(0);
          }
          twinput = twinput / sum(twinput);
          
          arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
            1, true, twinput);
          theprophecy(yr_counter) = theprophecy_piecemeal(0);
          
          arma::uvec tnotb_preassigned = find(thenumbersofthebeast == theprophecy_piecemeal(0));
          twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
        }
      }
      
      arma::mat projection;
      if (!matrix_class_input || sparse_switch == 1) {
        startvec = ss3matrix_sp(as<arma::sp_mat>(meanamats[i]));
        
        if (!matrix_class_input) {
          projection = proj3sp(startvec, amats, theprophecy, 1, 1, 0);
        } else {
          projection = proj3(startvec, amats, theprophecy, 1, 1, 0, sparse_auto,
            sparse_switch);
        }
      } else {
        startvec = ss3matrix(as<arma::mat>(meanamats[i]), sparse_switch);
        projection = proj3(startvec, amats, theprophecy, 1, 1, 0, sparse_auto,
          sparse_switch);
      }
      
      for (int j = 0; j < theclairvoyant; j++) {
        double madness = sum(projection.col(j+1));
        slmat(j,i) = log(madness);
      }
      
      sl_mean(i) = mean(slmat.col(i));
      sl_var(i) = var(slmat.col(i));
      sl_sd(i) = stddev(slmat.col(i));
      sl_se(i) = sl_sd(i) / sqrt(static_cast<double>(theclairvoyant));
    }
    
    int pop_est {1};
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize, fill::zeros);
    arma::uvec yearmatch(loysize, fill::zeros);
    List meanmatyearlist(uniqueyears.length());
    
    if (allppcsnem > 1) { // Check pop-mean matrices separate from patch means
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // population loop
        if (!matrix_class_input || sparse_switch == 1) {
          startvec = ss3matrix_sp(as<arma::sp_mat>(meanamats[allppcsnem + i]));
        } else {
          startvec = ss3matrix(as<arma::mat>(meanamats[allppcsnem + i]), sparse_switch);
        }
        
        // Checks which A matrices match current pop
        for (int j = 0; j < loysize; j++) { 
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        arma::uvec neededmatspop = find(popmatch);
        
        // Develop mean annual matrices across patches
        for (int j = 0; j < yl; j++) { 
          for (int k = 0; k < loysize; k++) {
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          
          // Catch matrix indices matching current year and pop
          int crankybankynem = static_cast<int>(crankybanky.n_elem);
          
          if (!matrix_class_input) {
            arma::sp_mat crossmat(meanmatsize, crankybankynem);
            for (int j = 0; j < crankybankynem; j++) {
              crossmat.col(j) = arma::vectorise(as<arma::sp_mat>(amats(crankybanky(j))));
            }
            
            arma::sp_mat happymedium(meanmatsize, 1);
            for (int j = 0; j < meanmatsize; j++) {
              for (int k = 0; k < crankybankynem; k++) {
                happymedium(j, 0) = happymedium(j, 0) + crossmat(j, k) / (crankybankynem);
              }
            }
            
            arma::sp_mat finalyearmat = happymedium;
            finalyearmat.reshape(meanmatrows, meanmatrows);
            meanmatyearlist(j) = finalyearmat;
            
          } else {
            arma::mat crossmat(meanmatsize, crankybankynem, fill::zeros);
            for (int j = 0; j < crankybankynem; j++) {
              crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
            }
            
            arma::vec happymedium(meanmatsize, fill::zeros);
            for (int j = 0; j < meanmatsize; j++) {
              for (int k = 0; k < crankybankynem; k++) {
                happymedium(j) = happymedium(j) + crossmat(j, k) / (crankybankynem);
              }
            }
            
            arma::mat finalyearmat = happymedium;
            finalyearmat.reshape(meanmatrows, meanmatrows);
            meanmatyearlist(j) = finalyearmat;
          }
        }
        
        int numyearsused = meanmatyearlist.length();
        arma::uvec choicevec = linspace<arma::uvec>(0, (numyearsused - 1), numyearsused);
        
        if (!assume_markov) {
          twinput = twinput / sum(twinput);
          
          theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant,
            true, twinput);
            
        } else {
          for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
            if (yr_counter == 0) {
              twinput = twinput_markov.col(0);
            }
            twinput = twinput / sum(twinput);
            
            arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(choicevec,
              1, true, twinput);
            theprophecy(yr_counter) = theprophecy_piecemeal(0);
            
            arma::uvec tnotb_preassigned = find(choicevec == theprophecy_piecemeal(0));
            twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
          }
        }
        
        arma::mat projection;
        if (is<S4>(meanmatyearlist(0))) { 
          projection = proj3sp(startvec, meanmatyearlist, theprophecy, 1, 1, 0);
        } else {
          projection = proj3(startvec, meanmatyearlist, theprophecy, 1, 1, 0,
            sparse_auto, sparse_switch);
        }
        
        for (int j = 0; j < theclairvoyant; j++) {
          double madness = sum(projection.col(j+1));
          slmat(j,(allppcsnem +i)) = log(madness);
        }
        
        sl_mean((allppcsnem + i)) = mean(slmat.col((allppcsnem +i)));
        sl_var((allppcsnem + i)) = var(slmat.col((allppcsnem +i)));
        sl_sd((allppcsnem + i)) = stddev(slmat.col((allppcsnem +i)));
        sl_se((allppcsnem + i)) = sl_sd((allppcsnem +i)) / sqrt(static_cast<double>(theclairvoyant));
      }
    }
    return DataFrame::create(_["pop"] = mmpops, _["patch"] = mmpatches,
      _["a"] = sl_mean, _["var"] = sl_var, _["sd"] = sl_sd, _["se"] = sl_se);
    
  } else {
    List amats = mpm;
    
    if (is<NumericMatrix>(amats(0))) matrix_class_input = true;
    
    // Check if matrix is large and sparse
    if (sparse_auto && matrix_class_input) {
      int test_elems = static_cast<int>(as<arma::mat>(amats(0)).n_elem);
      int Amatrows = static_cast<int>(as<arma::mat>(amats(0)).n_rows);
      arma::uvec nonzero_elems = find(as<arma::mat>(amats(0)));
      int all_nonzeros = static_cast<int>(nonzero_elems.n_elem);
      double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
      if (sparse_check <= 0.5 && Amatrows > 50) {
        sparse_switch = 1;
      } else sparse_switch = 0;
    }
    
    int yl = amats.length();
    arma::mat firstmat = as<arma::mat>(amats[0]);
    int matrows = firstmat.n_rows;
    int matcols = firstmat.n_cols;
    
    arma::uvec uniqueyears(yl);
    for (int i = 0; i < yl; i++) {
      uniqueyears(i) = i;
    }
    
    if (matrows != matcols) throw Rcpp::exception("Supplied matrices must be square.", false);
    
    arma::vec twinput;
    arma::mat twinput_markov;
    if (tweights.isNotNull()) {
      if (Rf_isMatrix(tweights)) {
        twinput_markov = as<arma::mat>(tweights);
        assume_markov = true;
        
        if (static_cast<int>(twinput_markov.n_cols) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
        }
        if (twinput_markov.n_cols != twinput_markov.n_rows) {
          throw Rcpp::exception("Time weight matrix must be square.", false);
        }
        
      } else if (is<NumericVector>(tweights)) {
        twinput = as<arma::vec>(tweights);
        
        if (static_cast<int>(twinput.n_elem) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
        }
        
      } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    // Create mean matrix
    arma::mat thechosenone;
    arma::sp_mat thechosenone_sp;
    
    if (!matrix_class_input) {
      arma::sp_mat thechosenone_(matrows, matcols);
      thechosenone_sp = thechosenone_;
      
      for (int i = 0; i < yl; i++) {
        thechosenone_sp = thechosenone_sp + (as<arma::sp_mat>(amats[i]) / yl);
      }
    } else {
      arma::mat thechosenone_(matrows, matcols, fill::zeros);
      thechosenone = thechosenone_;
      
      for (int i = 0; i < yl; i++) {
        thechosenone = thechosenone + (as<arma::mat>(amats[i]) / yl);
      }
    }
    
    // Run simulation, estimate all metrics
    arma::vec startvec;
    int trials {1};
    
    arma::mat slmat(theclairvoyant, trials);
    slmat.zeros();
    arma::vec sl_mean(trials);
    arma::vec sl_var(trials);
    arma::vec sl_sd(trials);
    arma::vec sl_se(trials);
    sl_mean.zeros();
    sl_var.zeros();
    sl_sd.zeros();
    sl_se.zeros();
    
    arma::uvec thenumbersofthebeast = uniqueyears;
    arma::uvec theprophecy (theclairvoyant);
    if (!assume_markov) {
      twinput = twinput / sum(twinput);
      theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
        theclairvoyant, true, twinput);
      
    } else if (assume_markov) {
      for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
        if (yr_counter == 0) {
          twinput = twinput_markov.col(0);
        }
        twinput = twinput / sum(twinput);
        
        arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(thenumbersofthebeast,
          1, true, twinput);
        theprophecy(yr_counter) = theprophecy_piecemeal(0);
          
        arma::uvec tnotb_preassigned = find(thenumbersofthebeast == theprophecy_piecemeal(0));
        twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
      }
    }
    
    arma::mat projection;
    if (!matrix_class_input) {
      startvec = ss3matrix_sp(thechosenone_sp);
      projection = proj3sp(startvec, amats, theprophecy, 1, 1, 0);
    } else {
      startvec = ss3matrix(thechosenone, sparse_switch);
      projection = proj3(startvec, amats, theprophecy, 1, 1, 0, sparse_auto, sparse_switch);
    }
    
    for (int j = 0; j < theclairvoyant; j++) {
      slmat(j,0) = sum(projection.col(j+1));
    }
    
    sl_mean(0) = mean(slmat.col(0));
    sl_var(0) = var(slmat.col(0));
    sl_sd(0) = stddev(slmat.col(0));
    sl_se(0) = sl_sd(0) / sqrt(static_cast<double>(theclairvoyant));
    
    CharacterVector mmpops(1);
    CharacterVector mmpatches(1);
    mmpops(0) = "1";
    mmpatches(0) = "0";
    
    return DataFrame::create(_["pop"] = mmpops, _["patch"] = mmpatches,
      _["a"] = sl_mean, _["var"] = sl_var, _["sd"] = sl_sd, _["se"] = sl_se);
  }
}

//' Estimate Stochastic Sensitivity or Elasticity of Matrix Set
//' 
//' Function \code{stoch_senselas()} estimates the sensitivity and elasticity to
//' matrix elements of \eqn{a}, defined as the long-term arithmetic mean of the
//' log population growth estimated per simulated occasion (as given in equation 2
//' in Tuljapurkar, Horvitz, and Pascarella 2003). 
//' 
//' @name .stoch_senselas
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' projection matrices.
//' @param times Number of occasions to iterate. Defaults to 10,000.
//' @param historical An optional logical value only used if object \code{mpm}
//' is a list of matrices, rather than a \code{lefkoMat} object. Defaults to
//' \code{FALSE} for the former case, and overridden by information supplied in
//' the \code{lefkoMat} object for the latter case.
//' @param style An integer designating whether to estimate sensitivity matrices
//' (\code{1}) or elasticity matrices (\code{2}). Defaults to \code{1}.
//' @param sparse An integer denoting whether to run the projection in sparse
//' (\code{1}) format or standard (\code{0}) format.
//' @param lefkoProj A logical value indicating whether the input MPM is a
//' \code{lefkoProj} object. Defaults to \code{TRUE}.
//' @param tweights An optional numeric vector or matrix denoting the
//' probabilities of choosing each matrix in a stochastic projection. If a
//' matrix is input, then a first-order Markovian environment is assumed, in
//' which the probability of choosing a specific annual matrix depends on which
//' annual matrix is currently chosen. If a vector is input, then the choice of
//' annual matrix is assumed to be independent of the current matrix. Defaults
//' to equal weighting among matrices.
//' 
//' @return A list of one or two cubes (3d array) where each slice corresponds
//' to a sensitivity or elasticity matrix for a specific pop-patch, followed by
//' the sensitivity or elasticity matrices of all populations (only if multiple
//' pop-patches occur in the input). Two such cubes are only provided when a
//' historical lefkoMat object is used as input, in which case the first
//' element is the historical sensitivity/elasticity matrix, and the second is
//' the ahistorical sensitivity/elasticity matrix.
//' 
//' @section Notes:
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector. Alternatively, a matrix may be
//' supplied, in which case it is assumed that the environment is first-order
//' Markovian.
//' 
//' This function currently requires all patches to have the same occasions, if
//' a \code{lefkoMat} object is used as input. Asymmetry in the number of
//' occasions across patches and/or populations will likely cause errors.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.stoch_senselas)]]
Rcpp::List stoch_senselas(const List& mpm, int times = 10000,
  bool historical = false, int style = 1, int sparse = 0, bool lefkoProj = true,
  Nullable<RObject> tweights = R_NilValue) {
  
  int theclairvoyant = times;
  if (theclairvoyant < 1) pop_error("times", "a positive integer", "", 1);
  
  bool lMat_matrix_class_input {false};
  bool lMat_sparse_class_input {false};
  bool list_input {false};
  bool assume_markov {false};
  
  if (lefkoProj) { 
    List A_list = as<List>(mpm["A"]);
    
    if (is<NumericMatrix>(A_list(0))) { 
      lMat_matrix_class_input = true;
    } else if (is<S4>(A_list(0))) {
      lMat_sparse_class_input = true;
    }
  } else {
    if (is<NumericMatrix>(mpm(0))) {
      list_input = true;
      lMat_matrix_class_input = true;
    } else if (is<S4>(mpm(0))) {
      list_input = true;
      lMat_sparse_class_input = true;
    }
  }
  
  if (!list_input && (lMat_matrix_class_input || lMat_sparse_class_input)) {
    List amats = as<List>(mpm["A"]);
    List umats = as<List>(mpm["U"]);
    List fmats = as<List>(mpm["F"]);
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    if (labels.length() < 3) {
      Rf_warningcall(R_NilValue,
        "This lefkoMat object appears to be a set of mean matrices. Will project only the mean.");
    }
    
    // Assess ahistorical versions of historical sensitivities / elasticities
    arma::uvec ahstages_id = as<arma::uvec>(stageframe["stage_id"]);
    StringVector ahstages_name = as<StringVector>(stageframe["stage"]);
    int ahstages_num = static_cast<int>(ahstages_id.n_elem);
    arma::uvec hstages_id2;
    int hstages_num {0};
    
    if (hstages.length() > 1) {
      arma::uvec hstages_id2_(ahstages_num * ahstages_num, fill::zeros);
      historical = true;
      arma::uvec hstages_id = as<arma::uvec>(hstages["stage_id_2"]);
      hstages_num = static_cast<int>(hstages_id.n_elem);
      
      for (int i = 0; i < hstages_num; i++) {
        hstages_id2_(i) = hstages_id(i);
      }
      
      hstages_id2 = hstages_id2_;
    } else {
      historical = false;
    }
    
    StringVector poporder = as<StringVector>(labels["pop"]);
    StringVector patchorder = as<StringVector>(labels["patch"]);
    StringVector yearorder = as<StringVector>(labels["year2"]);
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    StringVector uniqueyears = sort_unique(yearorder);
    
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    arma::vec twinput;
    arma::mat twinput_markov;
    if (tweights.isNotNull()) {
      if (Rf_isMatrix(tweights)) {
        twinput_markov = as<arma::mat>(tweights);
        assume_markov = true;
        
        if (static_cast<int>(twinput_markov.n_cols) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
        }
        if (twinput_markov.n_cols != twinput_markov.n_rows) {
          throw Rcpp::exception("Time weight matrix must be square.", false);
        }
        
      } else if (is<NumericVector>(tweights)) {
        twinput = as<arma::vec>(tweights);
        
        if (static_cast<int>(twinput.n_elem) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
        }
        
      } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    StringVector theprophecy_allyears (theclairvoyant);
    
    if (!assume_markov) {
      arma::vec twinput_corr = twinput / sum(twinput);
      theprophecy_allyears = Rcpp::RcppArmadillo::sample(uniqueyears,
        theclairvoyant, true, twinput_corr);
      
    } else if (assume_markov) {
      arma::vec twinput_corr;
      for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
        if (yr_counter == 0) {
          twinput = twinput_markov.col(0);
        }
        twinput_corr = twinput / sum(twinput);
        
        StringVector theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(uniqueyears,
          1, true, twinput_corr);
        theprophecy_allyears(yr_counter) = theprophecy_piecemeal(0);
        
        for (int yr_finder = 0; yr_finder < static_cast<int>(uniqueyears.length()); yr_finder++) {
          if (uniqueyears(yr_finder) == theprophecy_piecemeal(0)) {
            twinput = twinput_markov.col(yr_finder);
            break;
          }
        }
      }
    }
    
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize, fill::zeros);
    arma::uvec yearsinpatch(loysize, fill::zeros);
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = static_cast<int>(summervacation.n_elem);
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = static_cast<int>(motorhead.n_elem);
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    // Patch and pop means
    List mean_lefkomat;
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1, lMat_matrix_class_input, sparse);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1, lMat_matrix_class_input, sparse);
    }
    
    // Preliminaries for stochastic simulations
    List meanamats = as<List>(mean_lefkomat["A"]);
    List mmlabels = as<List>(mean_lefkomat["labels"]);
    StringVector mmpops = as<StringVector>(mmlabels["pop"]);
    StringVector mmpatches = as<StringVector>(mmlabels["patch"]);
    
    int meanmatsize {0};
    int meanmatrows {0};
    if (!is<S4>(meanamats(0))) {
      meanmatsize = static_cast<int>(as<arma::mat>(meanamats(0)).n_elem);
      meanmatrows = static_cast<int>(as<arma::mat>(meanamats(0)).n_rows);
    } else { 
      meanmatsize = static_cast<int>(as<arma::sp_mat>(meanamats(0)).n_elem);
      meanmatrows = static_cast<int>(as<arma::sp_mat>(meanamats(0)).n_rows);
    }
    
    arma::vec startvec(meanmatrows, fill::ones);
    startvec = startvec / meanmatrows; // Start vector for w and v calc
    int trials = meanamats.length();
    
    // Two lists for sensitivity/elasticity matrices
    // First general use, second for ahistorical versions of historical matrices
    List senscube (trials);
    List senscube_ah (trials);
    
    arma::mat sens_base;
    arma::mat sens_base_ah;
    arma::sp_mat sens_base_sp;
    arma::sp_mat sens_base_ah_sp;
    
    if (lMat_matrix_class_input && sparse == 0) {
      arma::mat sens_base_ (meanmatrows, meanmatrows, fill::zeros);
      arma::mat sens_base_ah_ (ahstages_num, ahstages_num, fill::zeros);
      sens_base = sens_base_;
      sens_base_ah = sens_base_ah_;
    } else {
      arma::sp_mat sens_base_ (meanmatrows, meanmatrows);
      arma::sp_mat sens_base_ah_ (ahstages_num, ahstages_num);
      sens_base_sp = sens_base_;
      sens_base_ah_sp = sens_base_ah_;
    }
    
    // Year value matrix for each run
    arma::umat yearspulled(trials, theclairvoyant, fill::zeros);
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = static_cast<int>(allppcs.n_elem);
    
    // Matrices and vectors for R values
    arma::mat Rvecmat(trials, theclairvoyant, fill::zeros);
    
    for (int i= 0; i < allppcsnem; i++) { //  pop-patch loop
      arma::uvec theprophecy (theprophecy_allyears.length(), fill::zeros);
      arma::uvec tnotb_patch = find(ppcindex == allppcs(i));
      
      for (int j = 0; j < yl; j++) { // Main index for used matrices
        // Modify in case patches do not have same years
        IntegerVector tnotb_years_IV = match(as<StringVector>(uniqueyears(j)), yearorder) - 1;
        arma::uvec tnotb_years = as<arma::uvec>(wrap(tnotb_years_IV));
        arma::uvec thenumbersofthebeast = intersect(tnotb_patch, tnotb_years);
        
        if (thenumbersofthebeast.n_elem > 0) {
          IntegerVector prophetic_yearindices_IV = match(as<StringVector>(uniqueyears(j)), 
              theprophecy_allyears) - 1;
          arma::uvec prophetic_yearindices = as<arma::uvec>(wrap(prophetic_yearindices_IV));
          
          if (prophetic_yearindices.n_elem > 0) {
            int replacement = thenumbersofthebeast(0);
            theprophecy.elem(prophetic_yearindices).fill(replacement);
          }
        }
      }
      yearspulled.row(i) = theprophecy.t();
      
      // Stable stage and rep value vectors, ahistorical versions of hMPMs
      arma::vec wprojection_ah(ahstages_num, fill::zeros);
      arma::vec vprojection_ah(ahstages_num, fill::zeros);
      
      // Control loop to develop w and v vectors
      arma::mat crazy_prophet;
      
      if (lMat_matrix_class_input) {
        crazy_prophet = proj3(startvec, amats, theprophecy, 1, 0, 0, false, sparse);
      } else {
        crazy_prophet = proj3sp(startvec, amats, theprophecy, 1, 0, 0);
      }
      arma::mat wprojection = crazy_prophet.submat(static_cast<int>(startvec.n_elem), 0,
          ((static_cast<int>(startvec.n_elem) * 2) - 1), theclairvoyant);
      arma::mat vprojection = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 2), 0,
        ((static_cast<int>(startvec.n_elem) * 3) - 1), theclairvoyant);
      Rvecmat.row(i) = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 3), 1,
        (static_cast<int>(startvec.n_elem) * 3), theclairvoyant); // Rvec
      
      // Main sensitivity matrix loop
      for (int j = 0; j < theclairvoyant; j++) {
        // Adds each occasion to matrix for each respective pop-patch
        if (j % 50 == 0) Rcpp::checkUserInterrupt();
        
        arma::vec vtplus1 = vprojection.col(j+1);
        arma::vec wtplus1 = wprojection.col(j+1);
        arma::vec wt = wprojection.col(j);
        
        arma::mat currentsens_num;
        arma::sp_mat currentsens_num_sp;
        arma::mat currentsens_den;
        arma::sp_mat currentsens_den_sp;
        arma::mat currentsens;
        arma::sp_mat currentsens_sp;
        
        if (lMat_matrix_class_input && sparse == 0) {
          currentsens_num = vtplus1 * wt.as_row(); // Key equation numerator
          currentsens_den = (Rvecmat(i, j) * vtplus1.as_row() * wtplus1); // Denominator
          double cd_double = static_cast<double>(currentsens_den(0,0));
          double downward_spiral = (cd_double * static_cast<double>(theclairvoyant));
          
          if (downward_spiral != 0.0) {
            currentsens = currentsens_num / downward_spiral;
          } else {
            arma::mat zero_mat (currentsens_num.n_rows, currentsens_num.n_cols, fill::zeros);
            currentsens = zero_mat;
          }
        } else {
          currentsens_num_sp = vtplus1 * wt.as_row(); // Key equation numerator
          currentsens_den_sp = (Rvecmat(i, j) * vtplus1.as_row() * wtplus1); // Denominator
          double cd_double = static_cast<double>(currentsens_den_sp(0,0));
          double downward_spiral = (cd_double * static_cast<double>(theclairvoyant));
          
          if (downward_spiral != 0.0) {
            currentsens_sp = currentsens_num_sp / downward_spiral;
          } else {
            arma::sp_mat zero_mat_sp (currentsens_num_sp.n_rows, currentsens_num_sp.n_cols);
            currentsens_sp = zero_mat_sp;
          }
        }
        
        // Creates sensitivity matrices
        if (style == 1) {
          if (lMat_matrix_class_input && sparse == 0) {
            sens_base += currentsens;
            senscube(i) = sens_base;
          } else { 
            sens_base_sp += currentsens_sp;
            
            if (sparse == 1) {
              senscube(i) = sens_base_sp;
            } else {
              senscube(i) = arma::mat(sens_base_sp);
            }
          }
          
          if (historical) {
            wprojection_ah.zeros();
            vprojection_ah.zeros();
            
            // Ahistorical stable stage dist for occasion j+1
            for (int k1 = 0; k1 < hstages_num; k1++) {
              int current_stage2 = hstages_id2(k1);
              wprojection_ah(current_stage2 - 1) = wprojection_ah(current_stage2 - 1)  +
                wtplus1(k1);
            } // k1 loop
            
           // Ahistorical rep value vector for occasion j+1
            for (int k2 = 0; k2 < hstages_num; k2++) {
              int current_stage2 = hstages_id2(k2);
              
              if (wprojection_ah(current_stage2 - 1) > 0) {
                vprojection_ah(current_stage2 - 1) = vprojection_ah(current_stage2 - 1) +
                  (vtplus1(k2) * wtplus1(k2) / wprojection_ah(current_stage2 - 1));
              }
            } // k2 loop
            
            // Propagate proj sens matrix, add to main sens matrix
            arma::rowvec wtah_tpose = wprojection_ah.as_row();
            arma::rowvec vtah_tpose = vprojection_ah.as_row();
            
            arma::mat csah_num;
            arma::sp_mat csah_num_sp;
            arma::mat csah_den;
            arma::sp_mat csah_den_sp;
            arma::mat csah;
            arma::sp_mat csah_sp;
            
            if (lMat_matrix_class_input && sparse == 0) {
              csah_num = vprojection_ah * wtah_tpose;
              csah_den = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
              double cdah_double = csah_den(0,0);
              double dismal_showing = cdah_double * static_cast<double>(theclairvoyant);
              
              if (dismal_showing != 0.0) {
                csah = csah_num / dismal_showing;
              } else {
                arma::mat zero_csah (csah_num.n_rows, csah_num.n_cols, fill::zeros);
                csah = zero_csah;
              }
              sens_base_ah += csah;
              senscube_ah(i) = sens_base_ah;
              
            } else {
              csah_num_sp = vprojection_ah * wtah_tpose;
              csah_den_sp = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
              double cdah_double = csah_den_sp(0,0);
              double dismal_showing = cdah_double * static_cast<double>(theclairvoyant);
              
              if (dismal_showing != 0.0) {
                csah_sp = csah_num_sp / dismal_showing;
              } else {
                arma::mat zero_csah_sp (csah_num_sp.n_rows, csah_num_sp.n_cols);
                csah_sp = zero_csah_sp;
              }
              sens_base_ah_sp += csah_sp;
              
              if (sparse == 1) {
                senscube_ah(i) = sens_base_ah_sp;
              } else {
                senscube_ah(i) = arma::mat(sens_base_ah_sp);
              }
            }
          } // if historical statement
        } else {
          // Elasticity matrices
          if (lMat_matrix_class_input && sparse == 0) {
            sens_base += currentsens % as<arma::mat>(amats[theprophecy(j)]);
            senscube(i) = sens_base;
          } else if (lMat_matrix_class_input) {
            arma::mat temp_mat = as<arma::mat>(amats[theprophecy(j)]);
            arma::sp_mat cs_sp = arma::sp_mat(temp_mat);
            
            cs_sp = currentsens_sp % cs_sp;
            sens_base_sp += cs_sp;
            
            if (sparse == 1) {
              senscube(i) = sens_base_sp;
            } else {
              senscube(i) = arma::mat(sens_base_sp);
            }
          } else {
            arma::sp_mat temp_spmat = as<arma::sp_mat>(amats[theprophecy(j)]);
            arma::sp_mat cs_sp = currentsens_sp % temp_spmat;
            sens_base_sp += cs_sp;
            
            if (sparse == 1) {
              senscube(i) = sens_base_sp;
            } else {
              senscube(i) = arma::mat(sens_base_sp);
            }
          }
        }
      }
    }
    
    // Pop means
    arma::mat pop_sens_base;
    arma::mat pop_sens_base_ah;
    arma::sp_mat pop_sens_base_sp;
    arma::sp_mat pop_sens_base_ah_sp;
    
    if (lMat_matrix_class_input && sparse == 0) {
      arma::mat sens_base_ (meanmatrows, meanmatrows, fill::zeros);
      arma::mat sens_base_ah_ (ahstages_num, ahstages_num, fill::zeros);
      pop_sens_base = sens_base_;
      pop_sens_base_ah = sens_base_ah_;
    } else {
      arma::sp_mat sens_base_ (meanmatrows, meanmatrows);
      arma::sp_mat sens_base_ah_ (ahstages_num, ahstages_num);
      pop_sens_base_sp = sens_base_;
      pop_sens_base_ah_sp = sens_base_ah_;
    }
    
    int pop_est {1};
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize, fill::zeros);
    arma::uvec yearmatch(loysize, fill::zeros);
    List meanmatyearlist(yl);
    
    IntegerVector tnotb_all = seq(0, (yl - 1));
    arma::uvec theprophecy (theprophecy_allyears.length(), fill::zeros);
    
    for (int j = 0; j < yl; j++) { // Main index marking matrices to use
      IntegerVector prophetic_yearindices_IV = match(as<StringVector>(uniqueyears(j)), 
          theprophecy_allyears) - 1;
      arma::uvec prophetic_yearindices = as<arma::uvec>(wrap(prophetic_yearindices_IV));
      if (prophetic_yearindices.n_elem > 0) {
        theprophecy.elem(prophetic_yearindices).fill(j);
      }
    }
    
    if (allppcsnem > 1) { // Checks for pop-mean matrices, only if >1 pop-patches
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // Pop loop
        for (int j = 0; j < loysize; j++) { // Finds A matrices matching pop
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        arma::uvec neededmatspop = find(popmatch == 1);
        
        for (int j = 0; j < yl; j++) { // Checks each year, develops patch mean matrix
          for (int k = 0; k < loysize; k++) { // Vector of matrices for current year
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          
          // Catch matrix indices matching current year and pop
          int crankybankynem = static_cast<int>(crankybanky.n_elem);
          arma::mat crossmat;
          arma::sp_mat crossmat_sp;
          
          if (is<NumericMatrix>(amats(0))) {
            arma::mat xmat(meanmatsize, crankybankynem, fill::zeros);
            crossmat = xmat;
            
            for (int j = 0; j < crankybankynem; j++) {
              crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
            }
            
            arma::mat happymedium(meanmatsize, 1);
            happymedium.zeros();
            for (int j = 0; j < meanmatsize; j++) {
              for (int k = 0; k < crankybankynem; k++) {
                happymedium(j, 0) = happymedium(j, 0) + crossmat(j, k) / crankybankynem;
              }
            }
            
            happymedium.reshape(meanmatrows, meanmatrows);
            meanmatyearlist(j) = happymedium;
          } else {
            arma::sp_mat xmat_sp(meanmatsize, crankybankynem);
            crossmat_sp = xmat_sp;
            
            for (int j = 0; j < crankybankynem; j++) {
              crossmat_sp.col(j) = arma::vectorise(as<arma::sp_mat>(amats(crankybanky(j))));
            }
            
            arma::sp_mat happymedium(meanmatsize, 1);
            for (int j = 0; j < meanmatsize; j++) {
              for (int k = 0; k < crankybankynem; k++) {
                happymedium(j, 0) = happymedium(j, 0) + crossmat_sp(j, k) / crankybankynem;
              }
            }
            
            happymedium.reshape(meanmatrows, meanmatrows);
            meanmatyearlist(j) = happymedium;
          }
        }
        yearspulled.row(allppcsnem + i) = theprophecy.t();
        
        arma::vec wprojection_ah(ahstages_num, fill::zeros);
        arma::vec vprojection_ah(ahstages_num, fill::zeros);
        
        // Control loop for w and v
        arma::mat crazy_prophet;
        if (is<NumericMatrix>(meanmatyearlist(0))) {
          crazy_prophet = proj3(startvec, meanmatyearlist, theprophecy, 1, 0, 0, false, sparse);
        } else { 
          crazy_prophet = proj3sp(startvec, meanmatyearlist, theprophecy, 1, 0, 0);
        }
        arma::mat wprojection = crazy_prophet.submat(static_cast<int>(startvec.n_elem),
            0, ((static_cast<int>(startvec.n_elem) * 2) - 1), theclairvoyant);
        arma::mat vprojection = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 2),
            0, ((static_cast<int>(startvec.n_elem) * 3) - 1), theclairvoyant);
        Rvecmat.row(allppcsnem + i) = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 3),
            1, (static_cast<int>(startvec.n_elem) * 3), theclairvoyant);
        
        // Sens matrix loop, adding each occasion to each pop-patch matrix
        for (int j = 0; j < theclairvoyant; j++) {
          arma::vec vtplus1 = vprojection.col(j+1);
          arma::vec wtplus1 = wprojection.col(j+1);
          arma::vec wt = wprojection.col(j);
          
          arma::mat currentsens_num;
          arma::sp_mat currentsens_num_sp;
          arma::mat currentsens_den;
          arma::sp_mat currentsens_den_sp;
          arma::mat currentsens;
          arma::sp_mat currentsens_sp;
          
          if (lMat_matrix_class_input && sparse == 0) {
            currentsens_num = vtplus1 * wt.as_row(); // Key equation numerator
            currentsens_den = (Rvecmat((allppcsnem + i), j) *
            vtplus1.as_row() * wtplus1); // Denominator
            double cd_double = currentsens_den(0,0);
            double downward_spiral = (cd_double * static_cast<double>(theclairvoyant));
            
            if (downward_spiral != 0.0) {
              currentsens = currentsens_num / downward_spiral;
            } else {
              arma::mat zero_mat (currentsens_num.n_rows, currentsens_num.n_cols, fill::zeros);
              currentsens = zero_mat;
            }
          } else {
            currentsens_num_sp = vtplus1 * wt.as_row(); // Key equation numerator
            currentsens_den_sp = (Rvecmat((allppcsnem + i), j) * vtplus1.as_row() *
              wtplus1); // Denominator
            double cd_double = currentsens_den_sp(0,0);
            double downward_spiral = (cd_double * static_cast<double>(theclairvoyant));
            
            if (downward_spiral != 0.0) {
              currentsens_sp = currentsens_num_sp / downward_spiral;
            } else {
              arma::sp_mat zero_mat_sp (currentsens_num.n_rows, currentsens_num.n_cols);
              currentsens_sp = zero_mat_sp;
            }
          }
          
          if (style == 1) {
            // Sensitivity matrix
            if (lMat_matrix_class_input && sparse == 0) {
              pop_sens_base += currentsens;
              senscube(allppcsnem + i) = pop_sens_base; 
            } else {
              pop_sens_base_sp += currentsens_sp;
              
              if (sparse == 1) {
                senscube(allppcsnem + i) = pop_sens_base_sp;
              } else {
                senscube(allppcsnem + i) = arma::mat(pop_sens_base_sp);
              }
            }
            
            if (historical) {
              wprojection_ah.zeros();
              vprojection_ah.zeros();
              
              // Ahistorical stable stage vector for occasion j+1
              for (int k1 = 0; k1 < hstages_num; k1++) {
                int current_stage2 = hstages_id2(k1);
                wprojection_ah(current_stage2 - 1) = wprojection_ah(current_stage2 - 1)  +
                  wtplus1(k1);
              } // k1 loop
              
              // Ahistorical rep value vector for occasion j+1
              for (int k2 = 0; k2 < hstages_num; k2++) {
                int current_stage2 = hstages_id2(k2);
                
                if (wprojection_ah(current_stage2 - 1) > 0) {
                  vprojection_ah(current_stage2 - 1) = vprojection_ah(current_stage2 - 1) +
                    (vtplus1(k2) * wtplus1(k2) / wprojection_ah(current_stage2 - 1));
                }
              } // k2 loop
              
              // Proj sens matrix, added to main sens matrix
              arma::rowvec wtah_tpose = wprojection_ah.as_row();
              arma::rowvec vtah_tpose = vprojection_ah.as_row();
              
              arma::mat csah_num;
              arma::sp_mat csah_num_sp;
              arma::mat csah_den;
              arma::sp_mat csah_den_sp;
              arma::mat csah;
              arma::sp_mat csah_sp;
              
              if (lMat_matrix_class_input && sparse == 0) {
                csah_num = vprojection_ah * wtah_tpose;
                csah_den = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
                
                double cdah_double = csah_den(0,0);
                double downward_spiral = (cdah_double * static_cast<double>(theclairvoyant));
                
                if (downward_spiral != 0.0) {
                  csah = csah_num / downward_spiral;
                } else {
                  arma::mat zero_csah (csah_num.n_rows, csah_num.n_cols, fill::zeros);
                  csah = zero_csah;
                }
                pop_sens_base_ah += csah;
                senscube_ah(allppcsnem + i) = pop_sens_base_ah;
                
              } else {
                csah_num_sp = vprojection_ah * wtah_tpose;
                csah_den_sp = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
                
                double cdah_double = csah_den_sp(0,0);
                double downward_spiral = (cdah_double * static_cast<double>(theclairvoyant));
                
                if (downward_spiral > 0.0) {
                  csah_sp = csah_num_sp / downward_spiral;
                } else {
                  arma::sp_mat zero_csah_sp (csah_num_sp.n_rows, csah_num_sp.n_cols);
                  csah_sp = zero_csah_sp;
                }
                pop_sens_base_ah_sp += csah_sp;
                
                if (sparse == 1) {
                  senscube_ah(allppcsnem + i) = pop_sens_base_ah_sp;
                } else {
                  senscube_ah(allppcsnem + i) = arma::mat(pop_sens_base_ah_sp);
                }
              }
            } // if historical statement
          } else {
            // Elasticity matrix
            if (lMat_matrix_class_input && sparse == 0) {
              pop_sens_base += currentsens % as<arma::mat>(meanmatyearlist[(theprophecy(j))]);
              senscube(allppcsnem + i) = pop_sens_base;
            } else if (is <S4>(meanmatyearlist[(theprophecy(j))])) {
              arma::sp_mat cs_sp = currentsens_sp % as<arma::sp_mat>(meanmatyearlist[(theprophecy(j))]);
              pop_sens_base_sp += cs_sp;
              
              if (sparse == 1) {
                senscube(allppcsnem + i) = pop_sens_base_sp;
              } else {
                senscube(allppcsnem + i) = arma::mat(pop_sens_base_sp);
              }
            } else {
              arma::sp_mat cs_sp = arma::sp_mat(as<arma::mat>(meanmatyearlist[(theprophecy(j))]));
              cs_sp = currentsens_sp % cs_sp;
              pop_sens_base_sp += cs_sp;
              
              if (sparse == 1) {
                senscube(allppcsnem + i) = pop_sens_base_sp;
              } else {
                senscube(allppcsnem + i) = arma::mat(pop_sens_base_sp);
              }
            }
          }
        }
      } // for loop i, for populations
    } // if statement, checks if >1 patch, determines if pop means to be estimated
    
    if (historical && style == 2) {
      for (int k = 0; k < trials; k++) {
        arma::uvec hstages_id1 = as<arma::uvec>(hstages["stage_id_1"]);
        arma::uvec hstages_id2 = as<arma::uvec>(hstages["stage_id_2"]);
        
        if (sparse == 0) {
          arma::mat elasah(ahstages_num, ahstages_num, fill::zeros);
          
          arma::mat hslice = as<arma::mat>(senscube(k));
          for (int i = 0; i < hstages_num; i++) {
            for (int j = 0; j < hstages_num; j++) {
              elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) =
                elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) + hslice(j, i);
            }
          }
          senscube_ah(k) = elasah;
        } else {
          arma::sp_mat elasah(ahstages_num, ahstages_num);
          
          arma::sp_mat hslice = as<arma::sp_mat>(senscube(k));
          for (int i = 0; i < hstages_num; i++) {
            for (int j = 0; j < hstages_num; j++) {
              elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) =
                elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) + hslice(j, i);
            }
          }
          
          senscube_ah(k) = elasah;
        }
      }
    }
    
    return Rcpp::List::create(_["maincube"] = senscube, _["ahcube"] = senscube_ah);
    
  } else {
    // List of A matrices input
    List amats = mpm;
    
    int matrows {0};
    int matcols {0};
    
    if (lMat_matrix_class_input) {
      arma::mat firstmat = as<arma::mat>(amats(0));
      matrows = firstmat.n_rows;
      matcols = firstmat.n_cols;
      
    } else {
      arma::sp_mat firstmat = as<arma::sp_mat>(amats(0));
      matrows = firstmat.n_rows;
      matcols = firstmat.n_cols;
    }
    
    int yl = amats.length();
    
    IntegerVector uniqueyears = seq(0, (yl - 1));
    arma::uvec uniqueyears_arma = as<arma::uvec>(uniqueyears);
    
    if (matrows != matcols) throw Rcpp::exception("Input matrices must be square.", false);
    
    arma::vec twinput;
    arma::mat twinput_markov;
    if (tweights.isNotNull()) {
      if (Rf_isMatrix(tweights)) {
        twinput_markov = as<arma::mat>(tweights);
        assume_markov = true;
        
        if (static_cast<int>(twinput_markov.n_cols) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
        }
        if (twinput_markov.n_cols != twinput_markov.n_rows) {
          throw Rcpp::exception("Time weight matrix must be square.", false);
        }
        
      } else if (is<NumericVector>(tweights)) {
        twinput = as<arma::vec>(tweights);
        
        if (static_cast<int>(twinput.n_elem) != yl) {
          pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
        }
        
      } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    arma::uvec theprophecy (theclairvoyant);
    arma::vec twinput_corr;
    if (!assume_markov) {
      twinput_corr = twinput / sum(twinput);
      theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears_arma, theclairvoyant,
        true, twinput_corr);
      
    } else if (assume_markov) {
      for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
        if (yr_counter == 0) {
          twinput = twinput_markov.col(0);
        }
        twinput_corr = twinput / sum(twinput);
        
        arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(uniqueyears_arma,
          1, true, twinput_corr);
        theprophecy(yr_counter) = theprophecy_piecemeal(0);
        
        arma::uvec tnotb_preassigned = find(uniqueyears_arma == theprophecy_piecemeal(0));
        twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
      }
    } 
    
    // Initialize empty matrix, start vector for w and v. Matrix updated at each occasion
    arma::vec startvec(matrows);
    startvec.ones();
    startvec = startvec / matrows; // Start vector for w and v calculations
    int trials {1};
    
    // Flat cube to hold sensitivity or elasticity matrix
    List senscube (trials);
    arma::mat sens_base;
    arma::sp_mat sens_base_sp;
    
    if (lMat_matrix_class_input && sparse == 0) {
      arma::mat sens_base_ (matrows, matrows, fill::zeros);
      sens_base = sens_base_;
    } else {
      arma::sp_mat sens_base_ (matrows, matrows);
      sens_base_sp = sens_base_;
    }
    
    
    // Matrices to hold R values
    arma::mat Rvecmat(trials, theclairvoyant, fill::zeros);
    
    // Loop to develop w and v values
    arma::mat crazy_prophet;
    if (lMat_matrix_class_input) {
      crazy_prophet = proj3(startvec, amats, theprophecy, 1, 0, 0, false, sparse);
    } else {
      crazy_prophet = proj3sp(startvec, amats, theprophecy, 1, 0, 0);
    }
    arma::mat wprojection = crazy_prophet.submat(static_cast<int>(startvec.n_elem), 0,
        ((static_cast<int>(startvec.n_elem) * 2) - 1), theclairvoyant);
    arma::mat vprojection = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 2), 0,
      ((static_cast<int>(startvec.n_elem) * 3) - 1), theclairvoyant);
    Rvecmat.row(0) = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 3), 1,
      (static_cast<int>(startvec.n_elem) * 3), theclairvoyant); // Rvec
    
    // Main loop for sensitivity matrices
    for (int j = 0; j < theclairvoyant; j++) {
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num;
      arma::sp_mat currentsens_num_sp;
      arma::mat currentsens_den;
      arma::sp_mat currentsens_den_sp;
      arma::mat currentsens;
      arma::sp_mat currentsens_sp;
      
      if (lMat_matrix_class_input && sparse == 0) {
        currentsens_num = vtplus1 * wt.as_row(); // Key matrix equation numerator
        currentsens_den = (Rvecmat(0, j) * vtplus1.as_row() * 
          wprojection.col(j+1)); // Denominator
        double cd_double = currentsens_den(0,0);
        double downward_spiral = (cd_double * static_cast<double>(theclairvoyant));
        
        if (downward_spiral != 0.0) {
          currentsens = currentsens_num / downward_spiral;
        } else {
          arma::mat zero_cursen (currentsens_num.n_rows, currentsens_num.n_cols, fill::zeros);
          currentsens = zero_cursen;
        }
          
        if (style == 1) {
          sens_base += currentsens;
          senscube(0) = sens_base; // Sensitivity matrix
        } else {
          sens_base += currentsens % as<arma::mat>(amats[(theprophecy(j))]);
          senscube(0) = sens_base; // Elasticity matrix
        }
      } else {
        currentsens_num_sp = vtplus1 * wt.as_row(); // Key matrix equation numerator
        currentsens_den_sp = (Rvecmat(0, j) * vtplus1.as_row() * 
          wprojection.col(j+1)); // Denominator
        double cd_double = currentsens_den_sp(0,0);
        double downward_spiral = (cd_double * static_cast<double>(theclairvoyant));

        if (downward_spiral > 0.0) {
          currentsens_sp = currentsens_num_sp / downward_spiral;
        } else {
          arma::sp_mat zero_cursen_sp (currentsens_num_sp.n_rows, currentsens_num_sp.n_cols);
          currentsens_sp = zero_cursen_sp;
        }
        
        if (style == 1) {
          sens_base_sp += currentsens_sp;
          
          if (sparse == 1) {
            senscube(0) = sens_base_sp; // Sensitivity matrix
          } else {
            senscube(0) = arma::mat(sens_base_sp);
          }
        } else {
          if (lMat_sparse_class_input) {
            arma::sp_mat cs_sp = currentsens_sp % as<arma::sp_mat>(amats[(theprophecy(j))]);
            sens_base_sp += cs_sp;
            
            if (sparse == 1) {
              senscube(0) = sens_base_sp; // Elasticity matrix
            } else {
              senscube(0) = arma::mat(sens_base_sp); // Elasticity matrix
            }
          } else {
            arma::sp_mat cs_sp = currentsens_sp % arma::sp_mat(as<arma::mat>(amats[(theprophecy(j))]));
            sens_base_sp += cs_sp;
            
            if (sparse == 1) {
              senscube(0) = sens_base_sp; // Elasticity matrix
            } else {
              senscube(0) = arma::mat(sens_base_sp); // Elasticity matrix
            }
          }
        }
      }
    }
    
  return Rcpp::List::create(_["maincube"] = senscube);
  }
}

//' Estimate LTRE of Any Population Matrix
//' 
//' \code{ltre3matrix()} returns the one-way fixed deterministic LTRE matrix of
//' a dense or sparse set of input matrices.
//' 
//' @name .ltre3matrix
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from_ \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param mean A logical value indicating whether to use the element-wise mean
//' matrix as the reference.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a one-element list with a list of LTRE
//' contributions, each element a matrix of contributions corresponding to each
//' input matrix in \code{Amats}. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ltre3matrix)]]
Rcpp::List ltre3matrix(const List& Amats, Rcpp::IntegerVector refnum,
  Nullable<Rcpp::List> refmats_ = R_NilValue, bool mean = true,
  bool sparse = false) {
  
  bool sparse_input {false};
  if (is<S4>(Amats(0))) sparse_input = true;
  
  int Amatnum = Amats.length();
  int matdim {0};
  if (sparse_input) {
    matdim = as<arma::sp_mat>(Amats(0)).n_rows;
  } else {
    matdim = as<arma::mat>(Amats(0)).n_rows;
  }
  
  // Check refnum vector
  int refmatnum = refnum.length();
  for (int i = 0; i < refmatnum; i++) {
    refnum(i) = refnum(i) - 1;
    if (refnum(i) < 0) {
      throw Rcpp::exception("Use R indexing for reference matrices (min: 1).", false);
    }
  }
  
  List eigenstuff;
  List mean_set;
  List cont_list(Amatnum);
  
  if (!sparse && !sparse_input) {
    // Dense matrix analysis
    arma::mat finalrefmat(matdim, matdim, fill::zeros);
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::mat>(refmats(refnum(i))) /
            static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = as<arma::mat>(refmats(refnum(0)));
      }
      
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::mat>(Amats(refnum(i))) /
            static_cast<double>(Amatnum));
        }
      } else {
        finalrefmat = as<arma::mat>(Amats(refnum(0)));
      }
    }
    
    // Create halfway matrices and run sensitivities
    arma::mat halfmat(matdim, matdim, fill::zeros);
    arma::mat diffmat(matdim, matdim, fill::zeros);
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (as<arma::mat>(Amats(i)))) / static_cast<double>(2.0);
      diffmat = (as<arma::mat>(Amats(i))) - finalrefmat;
      
      cont_list(i) = diffmat % sens3matrix(halfmat, 0);
    }
  } else if (!sparse_input) {
    // Sparse matrix analysis with dense input
    arma::sp_mat finalrefmat(matdim, matdim);
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) / 
            static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = (arma::sp_mat(as<arma::mat>(refmats(refnum(0)))));
      }
      
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (arma::sp_mat(as<arma::mat>(Amats(refnum(i)))) / 
            static_cast<double>(Amatnum));
        }
        
      } else {
        finalrefmat = (arma::sp_mat(as<arma::mat>(Amats(refnum(0)))));
      }
    }
    
    // Create halfway matrices and run sensitivities
    arma::sp_mat halfmat(matdim, matdim);
    arma::sp_mat diffmat(matdim, matdim);
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (arma::sp_mat(as<arma::mat>(Amats(i))))) /
        static_cast<double>(2.0);
      diffmat = (arma::sp_mat(as<arma::mat>(Amats(i)))) - finalrefmat;
      arma::sp_mat set_up = diffmat % sens3sp_matrix(halfmat, diffmat);
      cont_list(i) = set_up;
    }
  } else {
    // Sparse matrix analysis with sparse input
    arma::sp_mat finalrefmat(matdim, matdim);
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::sp_mat>(refmats(refnum(i))) / 
            static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = (as<arma::sp_mat>(refmats(refnum(0))));
      }
      
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::sp_mat>(Amats(refnum(i))) / 
            static_cast<double>(Amatnum));
        }
        
      } else {
        finalrefmat = (as<arma::sp_mat>(Amats(refnum(0))));
      }
    }
    
    // Create halfway matrices and run sensitivities
    arma::sp_mat halfmat(matdim, matdim);
    arma::sp_mat diffmat(matdim, matdim);
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (as<arma::sp_mat>(Amats(i)))) /
        static_cast<double>(2.0);
      diffmat = (as<arma::sp_mat>(Amats(i))) - finalrefmat;
      
      arma::sp_mat outmat = diffmat % sens3sp_matrix(halfmat, diffmat);
      cont_list(i) = outmat;
    }
  }
  
  List output = List::create(_["cont_mean"] = cont_list);
  CharacterVector out_class = {"lefkoLTRE"};
  output.attr("class") = out_class;
  
  return output;
}

//' Estimate sLTRE of Any Population Matrix
//' 
//' \code{sltre3matrix()} returns the one-way stochastic LTRE of a dense or
//' sparse set of input matrices.
//' 
//' @name .sltre3matrix
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object).
//' @param labels The data frame included in the input \code{lefkoMat} object.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param tweights_ Numeric vector or matrix denoting the probabilistic
//' weightings of annual matrices. Defaults to equal weighting among occasions.
//' @param steps The number of occasions to project the stochastic simulation
//' forward, if performing an sLTRE. Defaults to \code{10000}. Note that the
//' total number of occasions projected equals this number plus the number given
//' in object \code{burnin}.
//' @param burnin The number of initial occasions to project the population
//' without calculating population metrics. Defaults to \code{3000}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' @param tol_used A double precision numeric value indicating a lower positive
//' limit to matrix element values used in calculations. Matrix elements lower
//' than this value will be treated as \code{0.0} values.
//' 
//' @return This function returns a list of two lists of matrices. The first,
//' \code{cont_mean}, holds the sLTRE contributions of shifts in mean elements.
//' The second, \code{cont_sd}, holds the sLTRE contributions of shifts in
//' temporal standard deviations of matrix elements.
//' 
//' @section Notes:
//' This function uses the simulation approach developed in Davison et al.
//' (2010), which is a good approximation though not an analytical solution.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sltre3matrix)]]
Rcpp::List sltre3matrix(const List& Amats, const DataFrame& labels,
  Rcpp::IntegerVector refnum, Nullable<Rcpp::List> refmats_ = R_NilValue,
  Nullable<RObject> tweights_ = R_NilValue, int steps = 10000,
  int burnin = 3000, bool sparse = false, double tol_used = 1e-30) {
  
  bool sparse_input {false};
  bool assume_markov {false};
  if (is<S4>(Amats(0))) sparse_input = true;
  
  int refmatnum = refnum.length();
  for (int i = 0; i < refmatnum; i++) {
    refnum(i) = refnum(i) - 1;
    if (refnum(i) < 0) {
      throw Rcpp::exception("Use R indexing for reference matrices (min: 1).", false);
    }
  }
  
  int theclairvoyant = steps + burnin;
  
  int Amatnum = Amats.length();
  int matdim {0};
  int matlength {0};
  if (sparse_input) {
    matdim = static_cast<int>(as<arma::sp_mat>(Amats(0)).n_rows);
    matlength = static_cast<int>(as<arma::sp_mat>(Amats(0)).n_elem);
  } else {
    matdim = static_cast<int>(as<arma::mat>(Amats(0)).n_rows);
    matlength = static_cast<int>(as<arma::mat>(Amats(0)).n_elem);
  }
  
  List eigenstuff;
  List mean_set;
  List poppatch_meanmat;
  List poppatch_sdmat;
  List ref_byyear;
  List cont_meanmat;
  List cont_sdmat;
  
  // Create order vectors for pop-patch-year
  DataFrame loy = loy_inator(labels, true);
  
  StringVector pops = as<StringVector>(loy["pop"]);
  arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
  StringVector patches = as<StringVector>(loy["patch"]);
  arma::uvec year2c = as<arma::uvec>(loy["year2c"]);
  arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
  arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
  arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  arma::uvec uniqueyears = unique(year2c);
  int numpoppatches = static_cast<int>(uniquepoppatches.n_elem);
  int numyears = static_cast<int>(uniqueyears.n_elem);
  
  List poppatch_meanmat_temp (numpoppatches);
  List poppatch_sdmat_temp (numpoppatches);
  
  arma::mat ref_matmean;
  arma::mat ref_matsd;
  arma::sp_mat ref_sp_matmean;
  arma::sp_mat ref_sp_matsd;
  
  // First the pop/patch means and sds
  for (int i = 0; i < numpoppatches; i++) {
    arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
    int numpoppatch_chosen = static_cast<int>(poppatch_chosen.n_elem);
    arma::vec mat_elems(numpoppatch_chosen);
    
    if (!sparse && !sparse_input) {
      arma::mat mat_mean(matdim, matdim, fill::zeros);
      arma::mat mat_sd(matdim, matdim, fill::zeros);
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (as<arma::mat>(Amats(poppatch_chosen(j))) / 
          static_cast<double>(numpoppatch_chosen));
      }
      poppatch_meanmat_temp(i) = mat_mean;
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        bool found_greater {false};
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          double elem_check_mat_elems = as<arma::mat>(Amats(poppatch_chosen(k)))(j);
          if (elem_check_mat_elems > tol_used) {
            mat_elems(k) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(j) = arma::stddev(mat_elems, 0);
      }
      poppatch_sdmat_temp(i) = mat_sd;
      
    } else if (!sparse_input) {
      arma::sp_mat mat_mean(matdim, matdim);
      arma::sp_mat mat_sd(matdim, matdim);
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(poppatch_chosen(j)))) / 
          static_cast<double>(numpoppatch_chosen));
      }
      poppatch_meanmat_temp(i) = mat_mean;
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        bool found_greater {false};
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          double elem_check_mat_elems = arma::sp_mat(as<arma::mat>(Amats(poppatch_chosen(k))))(j);
          if (elem_check_mat_elems > tol_used) {
            mat_elems(k) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(j) = arma::stddev(mat_elems, 0);
      }
      poppatch_sdmat_temp(i) = mat_sd;
      
    } else {
      arma::sp_mat mat_mean(matdim, matdim);
      arma::sp_mat mat_sd(matdim, matdim);
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (as<arma::sp_mat>(Amats(poppatch_chosen(j))) / 
          static_cast<double>(numpoppatch_chosen));
      }
      poppatch_meanmat_temp(i) = mat_mean;
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        bool found_greater {false};
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          double elem_check_mat_elems = as<arma::sp_mat>(Amats(poppatch_chosen(k)))(j);
          if (elem_check_mat_elems > tol_used) {
            mat_elems(k) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(j) = arma::stddev(mat_elems, 0);
      }
      poppatch_sdmat_temp(i) = mat_sd;
    }
  }
  poppatch_meanmat = poppatch_meanmat_temp;
  poppatch_sdmat = poppatch_sdmat_temp;
  
  // Create reference matrix sets
  if (refmats_.isNotNull()) {
    Rcpp::List refmats(refmats_);
    ref_byyear = refmats;
    
    if (!sparse && !sparse_input) {
      arma::mat mat_mean(matdim, matdim, fill::zeros);
      arma::mat mat_sd(matdim, matdim, fill::zeros);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (as<arma::mat>(refmats(refnum(i))) / static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum, fill::zeros);
        
        bool found_greater {false};
        
        for (int j = 0; j < refmatnum; j++) {
          double elem_check_mat_elems = as<arma::mat>(refmats(refnum(j)))(i);
          if (elem_check_mat_elems > tol_used) {
            mat_elems(j) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
      }
      ref_matmean = mat_mean;
      ref_matsd = mat_sd;
      
    } else if (!sparse_input) {
      arma::sp_mat mat_mean(matdim, matdim);
      arma::sp_mat mat_sd(matdim, matdim);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) /
          static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum, fill::zeros);
        
        bool found_greater {false};
        
        for (int j = 0; j < refmatnum; j++) {
          double elem_check_mat_elems = arma::sp_mat(as<arma::mat>(refmats(refnum(j))))(i);
          if (elem_check_mat_elems > tol_used) {
            mat_elems(j) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
      }
      ref_sp_matmean = mat_mean;
      ref_sp_matsd = mat_sd;
      
    } else {
      arma::sp_mat mat_mean(matdim, matdim);
      arma::sp_mat mat_sd(matdim, matdim);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (as<arma::sp_mat>(refmats(refnum(i))) /
          static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum, fill::zeros);
        
        bool found_greater {false};
        
        for (int j = 0; j < refmatnum; j++) {
          double elem_check_mat_elems = as<arma::sp_mat>(refmats(refnum(j)))(i);
          if (elem_check_mat_elems > tol_used) {
            mat_elems(j) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
      }
      ref_sp_matmean = mat_mean;
      ref_sp_matsd = mat_sd;
    }
  } else {
    if (refmatnum == Amatnum) {
      
      // Reference by year
      List ref_byyear_temp (numyears);
      
      for (int i = 0; i < numyears; i++) {
        arma::uvec year_chosen = find(year2c == uniqueyears(i));
        int numyear_chosen = static_cast<int>(year_chosen.n_elem);
        
        if (!sparse && !sparse_input) {
          arma::mat mat_mean(matdim, matdim, fill::zeros);
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (as<arma::mat>(Amats(year_chosen(j))) /
              static_cast<double>(numyear_chosen));
          }
          ref_byyear_temp(i) = mat_mean;
        } else if (!sparse_input) {
          arma::sp_mat mat_mean(matdim, matdim);
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(year_chosen(j)))) /
              static_cast<double>(numyear_chosen));
          }
          ref_byyear_temp(i) = mat_mean;
        } else {
          arma::sp_mat mat_mean(matdim, matdim);
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (as<arma::sp_mat>(Amats(year_chosen(j))) /
              static_cast<double>(numyear_chosen));
          }
          ref_byyear_temp(i) = mat_mean;
        }
      }
      ref_byyear = ref_byyear_temp;
      
      // Reference mean and sd
      if (!sparse && !sparse_input) {
        arma::mat mat_mean(matdim, matdim, fill::zeros);
        arma::mat mat_sd(matdim, matdim, fill::zeros);
        
        for (int i = 0; i < numyears; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) / static_cast<double>(numyears));
        }
        
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < numyears; j++) {
            double elem_check_mat_elems = as<arma::mat>(ref_byyear(j))(i);
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        
      } else {
        arma::sp_mat mat_mean(matdim, matdim);
        arma::sp_mat mat_sd(matdim, matdim);
        
        for (int i = 0; i < numyears; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) / static_cast<double>(numyears));
        }
        
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < numyears; j++) {
            double elem_check_mat_elems = as<arma::sp_mat>(ref_byyear(j))(i);
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
        }
        ref_sp_matmean = mat_mean;
        ref_sp_matsd = mat_sd;
      }
    } else {
      // Reference by year
      List ref_byyear_temp (refmatnum);
      
      for (int i = 0; i < refmatnum; i++) {
        if (!sparse && !sparse_input) {
          ref_byyear_temp(i) = as<arma::mat>(Amats(i));
        } else {
          ref_byyear_temp(i) = arma::sp_mat(as<arma::mat>(Amats(i)));
        }
      }
      
      // Reference mean and sd
      if (!sparse && !sparse_input) {
        arma::mat mat_mean(matdim, matdim, fill::zeros);
        arma::mat mat_sd(matdim, matdim, fill::zeros);
        
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) /
            static_cast<double>(refmatnum));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < refmatnum; j++) {
            double elem_check_mat_elems = as<arma::mat>(ref_byyear(j))(i);
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        
      } else {
        arma::sp_mat mat_mean(matdim, matdim);
        arma::sp_mat mat_sd(matdim, matdim);
        
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) /
            static_cast<double>(refmatnum));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < refmatnum; j++) {
            double elem_check_mat_elems = as<arma::sp_mat>(ref_byyear(j))(i);
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(i) = arma::stddev(mat_elems, 0);
        }
        ref_sp_matmean = mat_mean;
        ref_sp_matsd = mat_sd;
      }
    }
  }
  
  // Time weights
  arma::vec twinput;
  arma::mat twinput_markov;
  if (tweights_.isNotNull()) {
    if (Rf_isMatrix(tweights_)) {
      twinput_markov = as<arma::mat>(tweights_);
      assume_markov = true;
      
      if (static_cast<int>(twinput_markov.n_cols) != numyears) {
        pop_error("tweights", "occasions", "lefkoMat object used as input", 20);
      }
      if (twinput_markov.n_cols != twinput_markov.n_rows) {
        throw Rcpp::exception("Time weight matrix must be square.", false);
      }
      
    } else if (is<NumericVector>(tweights_)) {
      twinput = as<arma::vec>(tweights_);
      
      if (static_cast<int>(twinput.n_elem) != numyears) {
        pop_error("tweights", "occasions", "lefkoMat object used as input", 19);
      }
      
    } else pop_error("tweights", "a non-negative numeric vector or matrix", "", 1);
  } else {
    twinput.resize(numyears);
    twinput.ones();
  }
  
  // Vector of chosen occasions, sampled from all possible occasions
  arma::uvec theprophecy (theclairvoyant);
  if (!assume_markov) {
    arma::vec twinput_corr = twinput / sum(twinput);
    theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears, theclairvoyant, true,
      twinput_corr);
    
  } else if (assume_markov) {
    arma::vec twinput_corr;
    for (int yr_counter = 0; yr_counter < theclairvoyant; yr_counter++) {
      if (yr_counter == 0) {
        twinput = twinput_markov.col(0);
      }
      twinput_corr = twinput / sum(twinput);
      
      arma::uvec theprophecy_piecemeal = Rcpp::RcppArmadillo::sample(uniqueyears,
        1, true, twinput_corr);
      theprophecy(yr_counter) = theprophecy_piecemeal(0);
      
      arma::uvec tnotb_preassigned = find(uniqueyears == theprophecy_piecemeal(0));
      twinput = twinput_markov.col(static_cast<int>(tnotb_preassigned(0)));
    }
  }
  
  arma::vec startvec(matdim, fill::ones);
  startvec = startvec / matdim; // This is the start vector for w and v calculations
  
  // Stochastic elasticities
  // Control loop creating w and v values from reference annual matrices
  arma::mat crazy_prophet;
  if (!sparse && !sparse_input) {
    crazy_prophet = proj3(startvec, ref_byyear, theprophecy, 1, 0, 0, false, sparse);
  } else {
    crazy_prophet = proj3sp(startvec, ref_byyear, theprophecy, 1, 0, 0);
  }
  arma::mat wprojection = crazy_prophet.submat(static_cast<int>(startvec.n_elem), 0,
      ((static_cast<int>(startvec.n_elem) * 2) - 1), theclairvoyant);
  arma::mat vprojection = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 2), 0,
    ((static_cast<int>(startvec.n_elem) * 3) - 1), theclairvoyant);
  arma::rowvec Rvecmat = crazy_prophet.submat((static_cast<int>(startvec.n_elem) * 3), 1,
      (static_cast<int>(startvec.n_elem) * 3), theclairvoyant); // Rvec
  
  if (!sparse && !sparse_input) {
    arma::mat sensmat(matdim, matdim, fill::zeros);
    arma::mat elasmean(matdim, matdim, fill::zeros);
    arma::mat elassd(matdim, matdim, fill::zeros);
    
    // Time loop for sensitivity matrices
    for (int j = burnin; j < theclairvoyant; j++) { 
      if (j % 10 == 0) Rcpp::checkUserInterrupt();
      
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wtplus1 = wprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // Key matrix equation numerator
      arma::mat currentsens_den = (Rvecmat(j) * vtplus1.as_row() * wtplus1); // Denominator
      double cd_double = currentsens_den(0,0);
      sensmat = currentsens_num / (cd_double);
      
      elasmean = elasmean + ((sensmat % ref_matmean) /
        static_cast<double>(theclairvoyant - burnin));
      elassd = elassd + ((sensmat % ((as<arma::mat>(ref_byyear(theprophecy(j)))) - ref_matmean)) / 
        static_cast<double>(theclairvoyant - burnin));
    }
    // Difference matrix estimation
    arma::mat diffmean1(matdim, matdim, fill::zeros);
    arma::mat diffsd1(matdim, matdim, fill::zeros);
    
    // Creates indices of 0 elements for next control loop
    diffmean1 = log(ref_matmean);
    diffsd1 = log(ref_matsd);
    diffmean1.elem(find_nonfinite(diffmean1)).zeros();
    diffsd1.elem(find_nonfinite(diffsd1)).zeros();
    ref_matmean = diffmean1;
    ref_matsd = diffsd1;
    
    List cont_meanmat_temp (numpoppatches);
    List cont_sdmat_temp (numpoppatches);
    
    for (int i = 0; i < numpoppatches; i++) {
      arma::mat diffmean2(matdim, matdim, fill::zeros);
      arma::mat diffsd2(matdim, matdim, fill::zeros);
      
      diffmean2 = log(as<arma::mat>(poppatch_meanmat(i)));
      diffsd2 = log(as<arma::mat>(poppatch_sdmat(i)));
      diffmean2.elem(find_nonfinite(diffmean2)).zeros();
      diffsd2.elem(find_nonfinite(diffsd2)).zeros();
      
      diffmean2 = diffmean2 - ref_matmean;
      diffsd2 = diffsd2 - ref_matsd;
      
      // Contributions
      diffmean2 = diffmean2 % elasmean;
      diffsd2 = diffsd2 % elassd;
      
      cont_meanmat_temp(i) = diffmean2;
      cont_sdmat_temp(i) = diffsd2;
    }
    cont_meanmat = cont_meanmat_temp;
    cont_sdmat = cont_sdmat_temp;
    
  } else {
    arma::sp_mat sensmat(matdim, matdim);
    arma::sp_mat elasmean(matdim, matdim);
    arma::sp_mat elassd(matdim, matdim);
    
    // Time loop for sensitivity matrices
    for (int j = burnin; j < theclairvoyant; j++) {
      if (j % 10 == 0) Rcpp::checkUserInterrupt();
      
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wtplus1 = wprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::sp_mat currentsens_num = arma::sp_mat(vtplus1) * arma::sp_mat(wt.as_row()); // Key matrix equation numerator 
      arma::sp_mat currentsens_den = arma::sp_mat(Rvecmat(j) * vtplus1.as_row() * wtplus1); // Denominator
      double cd_double = currentsens_den(0,0);
      sensmat = currentsens_num / (cd_double);
      
      elasmean = elasmean + ((sensmat % ref_sp_matmean) /
        static_cast<double>(theclairvoyant - burnin));
      elassd = elassd + ((sensmat % ((as<arma::sp_mat>(ref_byyear(theprophecy(j)))) - ref_sp_matmean)) / 
        static_cast<double>(theclairvoyant - burnin));
    }
    // Difference matrix estimation
    ref_sp_matmean = LefkoUtils::spmat_log(ref_sp_matmean);
    ref_sp_matsd = LefkoUtils::spmat_log(ref_sp_matsd);
    
    List cont_meanmat_temp (numpoppatches);
    List cont_sdmat_temp (numpoppatches);
    
    for (int i = 0; i < numpoppatches; i++) {
      arma::sp_mat diffmean1 = LefkoUtils::spmat_log(as<arma::sp_mat>(poppatch_meanmat(i))); // (matdim, matdim)
      arma::sp_mat diffsd1 = LefkoUtils::spmat_log(as<arma::sp_mat>(poppatch_sdmat(i))); // (matdim, matdim)
      
      diffmean1 = diffmean1 - ref_sp_matmean;
      diffsd1 = diffsd1 - ref_sp_matsd;
      
      // Contributions
      diffmean1 = diffmean1 % elasmean;
      diffsd1 = diffsd1 % elassd;
      
      cont_meanmat_temp(i) = diffmean1;
      cont_sdmat_temp(i) = diffsd1;
    }
    cont_meanmat = cont_meanmat_temp;
    cont_sdmat = cont_sdmat_temp;
  }
  
  List cont_list = List::create(Named("cont_mean") = cont_meanmat,
    _["cont_sd"] = cont_sdmat);
  StringVector output_class = {"lefkoLTRE"};
  cont_list.attr("class") = output_class;
  
  return cont_list;
}

//' Estimate SNA-LTRE of Any Population Matrix
//' 
//' \code{snaltre3matrix()} returns the one-way small noise approximation LTRE
//' of a dense or sparse set of input matrices.
//' 
//' @name .snaltre3matrix
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object).
//' @param labels The data frame included in the input \code{lefkoMat} object.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param tweights_ Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' @param tol_used A double precision numeric value indicating a lower positive
//' limit to matrix element values used in calculations. Matrix elements lower
//' than this value will be treated as \code{0.0} values.
//' 
//' @return This function returns a list of four lists of matrices. The first,
//' \code{cont_mean}, holds the sLTRE contributions of shifts in mean elements.
//' The second, \code{cont_elas}, holds the sLTRE contributions of shifts in
//' deterministic elasticity across matrices. The third, \code{cont_cv}, holds
//' the sLTRE contributions of shifts in temporal coefficients of variation of
//' matrix elements. The fourth, \code{cont_corr}, holds the contributions of
//' shifts in temporal correlations across matrices.
//' 
//' @section Notes:
//' This function uses the simulation approach developed in Davison et al.
//' (2019), which provides an analytical solution to the stochastic LTRE.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.snaltre3matrix)]]
Rcpp::List snaltre3matrix(const List& Amats, const DataFrame& labels,
  Rcpp::IntegerVector refnum, Nullable<Rcpp::List> refmats_ = R_NilValue,
  Nullable<arma::vec> tweights_ = R_NilValue, bool sparse = false,
  double tol_used = 1e-30) {
  
  bool sparse_input {false};
  if (is<S4>(Amats(0))) sparse_input = true;
  
  int refmatnum = refnum.length();
  for (int i = 0; i < refmatnum; i++) {
    refnum(i) = refnum(i) - 1;
    if (refnum(i) < 0) {
      throw Rcpp::exception("Use R indexing for reference matrices (min: 1).", false);
    }
  }
  
  arma::vec tweights;
  arma::ivec tweights_int;
  
  int Amatnum = Amats.length();
  int matdim {0};
  if (sparse_input) {
    matdim = as<arma::sp_mat>(Amats(0)).n_rows;
  } else {
    matdim = as<arma::mat>(Amats(0)).n_rows;
  }

  List eigenstuff;
  List mean_set;
  List ref_byyear;
  
  // Order vectors for pop-patch-year
  DataFrame loy = loy_inator(labels, true);
  
  StringVector pops = as<StringVector>(loy["pop"]);
  arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
  StringVector patches = as<StringVector>(loy["patch"]);
  arma::uvec year2c = as<arma::uvec>(loy["year2c"]);
  arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
  arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
  arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  arma::uvec uniqueyears = unique(year2c);
  int numpoppatches = static_cast<int>(uniquepoppatches.n_elem);
  int numyears = static_cast<int>(uniqueyears.n_elem);
  
  arma::vec twinput_corr;
  arma::uvec theprophecy;
  
  List poppatch_meanmat (numpoppatches);
  List poppatch_elasmat (numpoppatches);
  List poppatch_cvmat (numpoppatches);
  List poppatch_corrmat (numpoppatches);
  
  arma::mat ref_matmean;
  arma::mat ref_matelas;
  arma::mat ref_matcv;
  arma::mat ref_matcorr;
  arma::sp_mat ref_sp_matmean;
  arma::sp_mat ref_sp_matelas;
  arma::sp_mat ref_sp_matcv;
  arma::sp_mat ref_sp_matcorr;
  
  arma::uvec mat_index_main = LefkoMats::general_index(Amats, tol_used, true);
  int mat_index_main_nonzeros = static_cast<int>(mat_index_main.n_elem);
  
  if (tweights_.isNotNull()) {
    tweights = as<arma::vec>(tweights_);
    
    if (static_cast<int>(tweights.n_elem) != numyears) {
      pop_error("tweights", "occasions", "reference matrix set", 19);
    }
    
    double tw_multiplier = 1.0;
    bool all_ints {true};
    for (int i = 0; i < numyears; i++) {
      if ((static_cast<int>(tweights(i)) % 1) != 0) all_ints = false;
    }
    if (all_ints) {
      arma::uvec nonzero_tweights_index_pre = find(tweights);
      arma::vec nonzero_tweights_pre = tweights.elem(nonzero_tweights_index_pre);
      tw_multiplier = min(nonzero_tweights_pre);
    }
    
    double tw_sum = sum(tweights);
    tweights = tweights / tw_sum;
    
    arma::uvec unseemly_neg = find(tweights < 0.0);
    if (unseemly_neg.n_elem > 0) pop_error("tweights", "non-negative numeric vector or matrix", "", 1);
    
    arma::uvec nonzero_tweights_index = find(tweights);
    arma::vec nonzero_tweights = tweights.elem(nonzero_tweights_index);
    
    double tweights_min = min(nonzero_tweights);
    tweights = tweights / tweights_min;
    tweights = tweights * tw_multiplier;
    tweights_int = arma::conv_to<arma::ivec>::from(tweights);
    
  } else {
    arma::ivec tweights_int_temp (numyears, fill::ones);
    tweights_int = tweights_int_temp;
  }
  int tweights_total_mats = sum(tweights_int);
  
  NumericVector rvals_poppatch (numpoppatches);
  double rvals_ref {0.0};
  
  // First pop/patch means and sds
  for (int i = 0; i < numpoppatches; i++) {
    arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
    
    arma::uvec poppatch_extended (tweights_total_mats);
    int tw_counter {0};
    for (int j = 0; j < numyears; j++) {
      for (int k = 0; k < tweights_int(j); k++) { 
        poppatch_extended(tw_counter) = poppatch_chosen(j);
        tw_counter++;
      }
    }
    
    int numpoppatch_extended = static_cast<int>(poppatch_extended.n_elem);
    arma::vec mat_elems (numpoppatch_extended, fill::zeros);
    arma::vec mat_elems_forcorr (numpoppatch_extended, fill::zeros);
    
    if (!sparse && !sparse_input) {
      arma::mat mat_mean (matdim, matdim, fill::zeros);
      arma::mat mat_sd (matdim, matdim, fill::zeros);
      arma::mat mat_cv (matdim, matdim, fill::zeros);
      arma::mat mat_corr (matdim * matdim, matdim * matdim, fill::zeros);
      
      for (int j = 0; j < numpoppatch_extended; j++) {
        mat_mean = mat_mean + (as<arma::mat>(Amats(poppatch_extended(j))) / 
          static_cast<double>(numpoppatch_extended));
      }
      poppatch_meanmat(i) = mat_mean;
      
      List eigenstuff = LefkoMats::decomp3(mat_mean);
      arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
      double lambda = max(realeigenvals);
      rvals_poppatch(i) = log(lambda);
      
      arma::mat mat_elas = elas3matrix(mat_mean, false);
      poppatch_elasmat(i) = mat_elas;
      
      for (int j = 0; j < mat_index_main_nonzeros; j++) {
        mat_elems.zeros();
        bool found_greater {false};
        
        for (int k = 0; k < numpoppatch_extended; k++) {
          double elem_check_mat_elems = as<arma::mat>(Amats(poppatch_extended(k)))(mat_index_main(j));
          if (elem_check_mat_elems > tol_used) {
            mat_elems(k) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(mat_index_main(j)) = arma::stddev(mat_elems, 0);
        
        if (mat_sd(mat_index_main(j)) > tol_used) {
          if (mat_mean(mat_index_main(j)) > tol_used) {
            mat_cv(mat_index_main(j)) = mat_sd(mat_index_main(j)) / mat_mean(mat_index_main(j));
          } else {
            mat_cv(mat_index_main(j)) = 0.0;
          }
          
          for (int k = 0; k < mat_index_main_nonzeros; k++) {
            if (k % 10 == 0) Rcpp::checkUserInterrupt();
            
            mat_elems_forcorr.zeros();
            bool found_greater_forcorr {false};
            
            for (int l = 0; l < numpoppatch_extended; l++) {
              double elem_check_mef = as<arma::mat>(Amats(poppatch_extended(l)))(mat_index_main(k));
              if (elem_check_mef > tol_used) {
                mat_elems_forcorr(l) = elem_check_mef;
                found_greater_forcorr = true;
              }
            }
            
            double mefc_sd {0.0};
            if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
            
            if (mefc_sd > 0.0) {
              arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
              double corr_bit = new_cor(0);
              if (corr_bit > tol_used) mat_corr(mat_index_main(j), mat_index_main(k)) = corr_bit;
            }
          }
        }
      }
      poppatch_cvmat(i) = mat_cv;
      poppatch_corrmat(i) = mat_corr;
      
    } else if (!sparse_input) {
      arma::sp_mat mat_mean(matdim, matdim);
      arma::sp_mat mat_sd(matdim, matdim);
      arma::sp_mat mat_cv (matdim, matdim);
      arma::sp_mat mat_corr (matdim * matdim, matdim * matdim);
      
      for (int j = 0; j < numpoppatch_extended; j++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(poppatch_extended(j)))) / 
          static_cast<double>(numpoppatch_extended));
      }
      poppatch_meanmat(i) = mat_mean;
      
      List eigenstuff = LefkoMats::decomp3sp_inp(mat_mean);
      arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
      double lambda = max(realeigenvals);
      rvals_poppatch(i) = log(lambda);
      
      arma::sp_mat mat_elas = elas3sp_matrix(mat_mean);
      poppatch_elasmat(i) = mat_elas;
      
      for (int j = 0; j < mat_index_main_nonzeros; j++) {
        mat_elems.zeros();
        bool found_greater {false};
        
        for (int k = 0; k < numpoppatch_extended; k++) {
          double elem_check_mat_elems = arma::sp_mat(as<arma::mat>(Amats(poppatch_extended(k))))(mat_index_main(j));
          if (elem_check_mat_elems > tol_used) {
            mat_elems(k) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(mat_index_main(j)) = arma::stddev(mat_elems, 0);
        
        if (mat_sd(mat_index_main(j)) > tol_used) {
          if (mat_mean(mat_index_main(j)) > tol_used) {
            mat_cv(mat_index_main(j)) = mat_sd(mat_index_main(j)) / mat_mean(mat_index_main(j));
          } else {
            mat_cv(mat_index_main(j)) = 0.0;
          }
          
          for (int k = 0; k < mat_index_main_nonzeros; k++) {
            if (k % 10 == 0) Rcpp::checkUserInterrupt();
            
            mat_elems_forcorr.zeros();
            bool found_greater_forcorr {false};
            
            for (int l = 0; l < numpoppatch_extended; l++) {
              double elem_check_mef = as<arma::mat>(Amats(poppatch_extended(l)))(mat_index_main(k));
              if (elem_check_mef > tol_used) {
                mat_elems_forcorr(l) = elem_check_mef;
                found_greater_forcorr = true;
              }
            }
            
            double mefc_sd {0.0};
            if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
            
            if (mefc_sd > 0.0) {
              arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
              double corr_bit = new_cor(0);
              if (corr_bit > tol_used) mat_corr(mat_index_main(j), mat_index_main(k)) = corr_bit;
            }
          }
        }
      }
      poppatch_cvmat(i) = mat_cv;
      poppatch_corrmat(i) = mat_corr;
      
    } else {
      arma::sp_mat mat_mean(matdim, matdim);
      arma::sp_mat mat_sd(matdim, matdim);
      arma::sp_mat mat_cv (matdim, matdim);
      arma::sp_mat mat_corr (matdim * matdim, matdim * matdim);
      
      for (int j = 0; j < numpoppatch_extended; j++) {
        mat_mean = mat_mean + (as<arma::sp_mat>(Amats(poppatch_extended(j))) / 
          static_cast<double>(numpoppatch_extended));
      }
      poppatch_meanmat(i) = mat_mean;
      
      List eigenstuff = LefkoMats::decomp3sp_inp(mat_mean);
      arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
      double lambda = max(realeigenvals);
      rvals_poppatch(i) = log(lambda);
      
      arma::sp_mat mat_elas = elas3sp_matrix(mat_mean);
      poppatch_elasmat(i) = mat_elas;
      
      for (int j = 0; j < mat_index_main_nonzeros; j++) {
        mat_elems.zeros();
        bool found_greater {false};
        
        for (int k = 0; k < numpoppatch_extended; k++) {
          double elem_check_mat_elems = as<arma::sp_mat>(Amats(poppatch_extended(k)))(mat_index_main(j));
          if (elem_check_mat_elems > tol_used) {
            mat_elems(k) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(mat_index_main(j)) = arma::stddev(mat_elems, 0);
        
        if (mat_sd(mat_index_main(j)) > tol_used) {
          if (mat_mean(mat_index_main(j)) > tol_used) {
            mat_cv(mat_index_main(j)) = mat_sd(mat_index_main(j)) / mat_mean(mat_index_main(j));
          } else {
            mat_cv(mat_index_main(j)) = 0.0;
          }
          
          for (int k = 0; k < mat_index_main_nonzeros; k++) {
            if (k % 10 == 0) Rcpp::checkUserInterrupt();
            
            mat_elems_forcorr.zeros();
            bool found_greater_forcorr {false};
            
            for (int l = 0; l < numpoppatch_extended; l++) {
              double elem_check_mef = as<arma::sp_mat>(Amats(poppatch_extended(l)))(mat_index_main(k));
              if (elem_check_mef > tol_used) {
                mat_elems_forcorr(l) = elem_check_mef;
                found_greater_forcorr = true;
              }
            }
            
            double mefc_sd {0.0};
            if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
            
            if (mefc_sd > 0.0) {
              arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
              double corr_bit = new_cor(0);
              if (corr_bit > tol_used) mat_corr(mat_index_main(j), mat_index_main(k)) = corr_bit;
            }
          }
        }
      }
      poppatch_cvmat(i) = mat_cv;
      poppatch_corrmat(i) = mat_corr;
    }
  }
  
  // Create reference matrix sets
  if (refmats_.isNotNull()) {
    Rcpp::List refmats(refmats_);
    ref_byyear = refmats;
    
    arma::uvec ref_index_main = LefkoMats::general_index(refmats, tol_used, true);
    int ref_index_main_nonzeros = static_cast<int>(ref_index_main.n_elem);
    
    if (!sparse && !sparse_input) {
      arma::mat mat_mean (matdim, matdim, fill::zeros);
      arma::mat mat_sd (matdim, matdim, fill::zeros);
      arma::mat mat_cv (matdim, matdim, fill::zeros);
      arma::mat mat_corr (matdim * matdim, matdim * matdim, fill::zeros);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (as<arma::mat>(refmats(refnum(i))) / static_cast<double>(refmatnum));
      }
      
      List eigenstuff = LefkoMats::decomp3(mat_mean);
      arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
      double lambda = max(realeigenvals);
      rvals_ref = log(lambda);
      
      for (int i = 0; i < ref_index_main_nonzeros; i++) {
        arma::vec mat_elems(refmatnum, fill::zeros);
        arma::vec mat_elems_forcorr (refmatnum, fill::zeros);
        
        bool found_greater {false};
        
        for (int j = 0; j < refmatnum; j++) {
          double elem_check_mat_elems = as<arma::mat>(refmats(refnum(j)))(ref_index_main(i));
          if (elem_check_mat_elems > tol_used) {
            mat_elems(j) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(ref_index_main(i)) = arma::stddev(mat_elems, 0);
        
        if (mat_sd(ref_index_main(i)) > tol_used) {
          if (mat_mean(ref_index_main(i)) > tol_used) {
            mat_cv(ref_index_main(i)) = mat_sd(ref_index_main(i)) / mat_mean(ref_index_main(i));
          } else {
            mat_cv(ref_index_main(i)) = 0.0;
          }
          
          for (int k = 0; k < ref_index_main_nonzeros; k++) {
            if (k % 10 == 0) Rcpp::checkUserInterrupt();
            
            mat_elems_forcorr.zeros();
            bool found_greater_forcorr {false};
            
            for (int l = 0; l < refmatnum; l++) {
              double elem_check_mef = as<arma::mat>(refmats(refnum(l)))(ref_index_main(k));
              if (elem_check_mef > tol_used) {
                mat_elems_forcorr(l) = elem_check_mef;
                found_greater_forcorr = true;
              }
            }
            
            double mefc_sd {0.0};
            if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
            
            if (mefc_sd > 0.0) {
              arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
              double corr_bit = new_cor(0);
              if (corr_bit > tol_used) mat_corr(ref_index_main(i), ref_index_main(k)) = corr_bit;
            }
          }
        }
      }
      
      ref_matmean = mat_mean;
      ref_matelas = elas3matrix(mat_mean, false);
      ref_matcv = mat_cv;
      ref_matcorr = mat_corr;
      
    } else if (!sparse_input) {
      arma::sp_mat mat_mean (matdim, matdim);
      arma::sp_mat mat_sd (matdim, matdim);
      arma::sp_mat mat_cv (matdim, matdim);
      arma::sp_mat mat_corr (matdim * matdim, matdim * matdim);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) /
          static_cast<double>(refmatnum));
      }
      
      List eigenstuff = LefkoMats::decomp3sp_inp(mat_mean);
      arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
      double lambda = max(realeigenvals);
      rvals_ref = log(lambda);
      
      for (int i = 0; i < ref_index_main_nonzeros; i++) {
        arma::vec mat_elems(refmatnum, fill::zeros);
        arma::vec mat_elems_forcorr (refmatnum, fill::zeros);
        
        bool found_greater {false};
        
        for (int j = 0; j < refmatnum; j++) {
          double elem_check_mat_elems = arma::sp_mat(as<arma::mat>(refmats(refnum(j))))(ref_index_main(i));
          if (elem_check_mat_elems > tol_used) {
            mat_elems(j) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(ref_index_main(i)) = arma::stddev(mat_elems, 0);
        
        if (mat_sd(ref_index_main(i)) > tol_used) {
          if (mat_mean(ref_index_main(i)) > tol_used) {
            mat_cv(ref_index_main(i)) = mat_sd(ref_index_main(i)) / mat_mean(ref_index_main(i));
          } else {
            mat_cv(ref_index_main(i)) = 0.0;
          }
          
          for (int k = 0; k < ref_index_main_nonzeros; k++) {
            if (k % 10 == 0) Rcpp::checkUserInterrupt();
            
            mat_elems_forcorr.zeros();
            bool found_greater_forcorr {false};
            
            for (int l = 0; l < refmatnum; l++) {
              double elem_check_mef = as<arma::mat>(refmats(refnum(l)))(ref_index_main(k));
              if (elem_check_mef > tol_used) {
                mat_elems_forcorr(l) = elem_check_mef;
                found_greater_forcorr = true;
              }
            }
            
            double mefc_sd {0.0};
            if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
            
            if (mefc_sd > 0.0) {
              arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
              double corr_bit = new_cor(0);
              if (corr_bit > tol_used) mat_corr(ref_index_main(i), ref_index_main(k)) = corr_bit;
            }
          }
        }
      }
      
      ref_sp_matmean = mat_mean;
      ref_sp_matelas = elas3sp_matrix(mat_mean);
      ref_sp_matcv = mat_cv;
      ref_sp_matcorr = mat_corr;
      
    } else {
      arma::sp_mat mat_mean (matdim, matdim);
      arma::sp_mat mat_sd (matdim, matdim);
      arma::sp_mat mat_cv (matdim, matdim);
      arma::sp_mat mat_corr (matdim * matdim, matdim * matdim);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (as<arma::sp_mat>(refmats(refnum(i))) /
          static_cast<double>(refmatnum));
      }
      
      List eigenstuff = LefkoMats::decomp3sp_inp(mat_mean);
      arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
      double lambda = max(realeigenvals);
      rvals_ref = log(lambda);
      
      for (int i = 0; i < ref_index_main_nonzeros; i++) {
        arma::vec mat_elems(refmatnum, fill::zeros);
        arma::vec mat_elems_forcorr (refmatnum, fill::zeros);
        
        bool found_greater {false};
        
        for (int j = 0; j < refmatnum; j++) {
          double elem_check_mat_elems = as<arma::sp_mat>(refmats(refnum(j)))(ref_index_main(i));
          if (elem_check_mat_elems > tol_used) {
            mat_elems(j) = elem_check_mat_elems;
            found_greater = true;
          }
        }
        
        if (found_greater) mat_sd(ref_index_main(i)) = arma::stddev(mat_elems, 0);
        
        if (mat_sd(ref_index_main(i)) > tol_used) {
          if (mat_mean(ref_index_main(i)) > tol_used) {
            mat_cv(ref_index_main(i)) = mat_sd(ref_index_main(i)) / mat_mean(ref_index_main(i));
          } else {
            mat_cv(ref_index_main(i)) = 0.0;
          }
          
          for (int k = 0; k < ref_index_main_nonzeros; k++) {
            if (k % 10 == 0) Rcpp::checkUserInterrupt();
            
            mat_elems_forcorr.zeros();
            bool found_greater_forcorr {false};
            
            for (int l = 0; l < refmatnum; l++) {
              double elem_check_mef = as<arma::sp_mat>(refmats(refnum(l)))(ref_index_main(k));
              if (elem_check_mef > tol_used) {
                mat_elems_forcorr(l) = elem_check_mef;
                found_greater_forcorr = true;
              }
            }
            
            double mefc_sd {0.0};
            if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
            
            if (mefc_sd > 0.0) {
              arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
              double corr_bit = new_cor(0);
              if (corr_bit > tol_used) mat_corr(ref_index_main(i), ref_index_main(k)) = corr_bit;
            }
          }
        }
      }
      
      ref_sp_matmean = mat_mean;
      ref_sp_matelas = elas3sp_matrix(mat_mean);
      ref_sp_matcv = mat_cv;
      ref_sp_matcorr = mat_corr;
    }
  } else {
    if (refmatnum == Amatnum) {
      
      // Reference by year
      List ref_byyear_temp (numyears);
      
      for (int i = 0; i < numyears; i++) {
        arma::uvec year_chosen = find(year2c == uniqueyears(i));
        int numyear_chosen = static_cast<int>(year_chosen.n_elem);
        
        if (!sparse && !sparse_input) {
          arma::mat mat_mean(matdim, matdim, fill::zeros);
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (as<arma::mat>(Amats(year_chosen(j))) /
              static_cast<double>(numyear_chosen));
          }
          ref_byyear_temp(i) = mat_mean;
          
        } else if (!sparse_input) {
          arma::sp_mat mat_mean(matdim, matdim);
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(year_chosen(j)))) /
              static_cast<double>(numyear_chosen));
          }
          ref_byyear_temp(i) = mat_mean;
          
        } else {
          arma::sp_mat mat_mean(matdim, matdim);
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (as<arma::sp_mat>(Amats(year_chosen(j))) /
              static_cast<double>(numyear_chosen));
          }
          ref_byyear_temp(i) = mat_mean;
        }
      }
      ref_byyear = ref_byyear_temp;
      
      arma::uvec ref_extended (tweights_total_mats);
      int tw_counter {0};
      for (int j = 0; j < numyears; j++) {
        for (int k = 0; k < tweights_int(j); k++) { 
          ref_extended(tw_counter) = j;
          tw_counter++;
        }
      }
      int ref_extended_length = static_cast<int>(ref_extended.n_elem);
      
      // Reference mean and sd
      if (!sparse && !sparse_input) {
        arma::mat mat_mean(matdim, matdim, fill::zeros);
        arma::mat mat_sd(matdim, matdim, fill::zeros);
        arma::mat mat_cv (matdim, matdim, fill::zeros);
        arma::mat mat_corr (matdim * matdim, matdim * matdim, fill::zeros);
        
        for (int i = 0; i < ref_extended_length; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(ref_extended(i))) /
            static_cast<double>(ref_extended_length));
        }
        
        List eigenstuff = LefkoMats::decomp3(mat_mean);
        arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
        double lambda = max(realeigenvals);
        rvals_ref = log(lambda);
        
        for (int i = 0; i < mat_index_main_nonzeros; i++) {
          arma::vec mat_elems(ref_extended_length, fill::zeros);
          arma::vec mat_elems_forcorr (ref_extended_length, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < ref_extended_length; j++) {
            double elem_check_mat_elems = as<arma::mat>(ref_byyear(ref_extended(j)))(mat_index_main(i));
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(mat_index_main(i)) = arma::stddev(mat_elems, 0);
          
          if (mat_sd(mat_index_main(i)) > tol_used) {
            if (mat_mean(mat_index_main(i)) > tol_used) {
              mat_cv(mat_index_main(i)) = mat_sd(mat_index_main(i)) / mat_mean(mat_index_main(i));
            } else {
              mat_cv(mat_index_main(i)) = 0.0;
            }
            
            for (int k = 0; k < mat_index_main_nonzeros; k++) {
              if (k % 10 == 0) Rcpp::checkUserInterrupt();
              
              mat_elems_forcorr.zeros();
              bool found_greater_forcorr {false};
              
              for (int l = 0; l < ref_extended_length; l++) {
                double elem_check_mef = as<arma::mat>(ref_byyear(ref_extended(l)))(mat_index_main(k));
                if (elem_check_mef > tol_used) {
                  mat_elems_forcorr(l) = elem_check_mef;
                  found_greater_forcorr = true;
                }
              }
              
              double mefc_sd {0.0};
              if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
              
              if (mefc_sd > 0.0) {
                arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
                double corr_bit = new_cor(0);
                if (corr_bit > tol_used) mat_corr(mat_index_main(i), mat_index_main(k)) = corr_bit;
              }
            }
          }
        }
        
        ref_matmean = mat_mean;
        ref_matelas = elas3matrix(mat_mean, false);
        ref_matcv = mat_cv;
        ref_matcorr = mat_corr;
        
      } else {
        arma::sp_mat mat_mean(matdim, matdim);
        arma::sp_mat mat_sd(matdim, matdim);
        arma::sp_mat mat_cv (matdim, matdim);
        arma::sp_mat mat_corr (matdim * matdim, matdim * matdim);
        
        for (int i = 0; i < ref_extended_length; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(ref_extended(i))) /
            static_cast<double>(ref_extended_length));
        }
        
        List eigenstuff = LefkoMats::decomp3sp_inp(mat_mean);
        arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
        double lambda = max(realeigenvals);
        rvals_ref = log(lambda);
        
        for (int i = 0; i < mat_index_main_nonzeros; i++) {
          arma::vec mat_elems(ref_extended_length, fill::zeros);
          arma::vec mat_elems_forcorr (ref_extended_length, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < ref_extended_length; j++) {
            double elem_check_mat_elems = as<arma::sp_mat>(ref_byyear(ref_extended(j)))(mat_index_main(i));
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(mat_index_main(i)) = arma::stddev(mat_elems, 0);
          
          if (mat_sd(mat_index_main(i)) > tol_used) {
            if (mat_mean(mat_index_main(i)) > tol_used) {
              mat_cv(mat_index_main(i)) = mat_sd(mat_index_main(i)) / mat_mean(mat_index_main(i));
            } else {
              mat_cv(mat_index_main(i)) = 0.0;
            }
            
            // This is crazy - correlation matrix is ultra-high dimension...
            for (int k = 0; k < mat_index_main_nonzeros; k++) {
              if (k % 10 == 0) Rcpp::checkUserInterrupt();
              
              mat_elems_forcorr.zeros();
              bool found_greater_forcorr {false};
              
              for (int l = 0; l < ref_extended_length; l++) {
                double elem_check_mef = as<arma::sp_mat>(ref_byyear(ref_extended(l)))(mat_index_main(k));
                if (elem_check_mef > tol_used) {
                  mat_elems_forcorr(l) = elem_check_mef;
                  found_greater_forcorr = true;
                }
              }
              
              double mefc_sd {0.0};
              if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
              
              if (mefc_sd > 0.0) {
                arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
                double corr_bit = new_cor(0);
                if (corr_bit > tol_used) mat_corr(mat_index_main(i), mat_index_main(k)) = corr_bit;
              }
            }
          }
        }
        
        ref_sp_matmean = mat_mean;
        ref_sp_matelas = elas3sp_matrix(mat_mean);
        ref_sp_matcv = mat_cv;
        ref_sp_matcorr = mat_corr;
      }
      
    } else {
      // Reference by year
      List ref_byyear_temp (refmatnum);
      
      for (int i = 0; i < refmatnum; i++) {
        if (!sparse && !sparse_input) {
          ref_byyear_temp(i) = as<arma::mat>(Amats(i));
        } else if (!sparse_input) {
          ref_byyear_temp(i) = arma::sp_mat(as<arma::mat>(Amats(i)));
        } else {
          ref_byyear_temp(i) = as<arma::sp_mat>(Amats(i));
        }
      }
      ref_byyear = ref_byyear_temp;
      
      // Reference mean and sd
      if (!sparse) {
        arma::mat mat_mean (matdim, matdim, fill::zeros);
        arma::mat mat_sd (matdim, matdim, fill::zeros);
        arma::mat mat_cv (matdim, matdim, fill::zeros);
        arma::mat mat_corr (matdim * matdim, matdim * matdim, fill::zeros);
        
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) /
            static_cast<double>(refmatnum));
        }
        
        List eigenstuff = LefkoMats::decomp3(mat_mean);
        arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
        double lambda = max(realeigenvals);
        rvals_ref = log(lambda);
        
        for (int i = 0; i < mat_index_main_nonzeros; i++) {
          arma::vec mat_elems(refmatnum, fill::zeros);
          arma::vec mat_elems_forcorr (refmatnum, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < refmatnum; j++) {
            double elem_check_mat_elems = as<arma::mat>(ref_byyear(j))(mat_index_main(i));
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(mat_index_main(i)) = arma::stddev(mat_elems, 0);
          
          if (mat_sd(mat_index_main(i)) > tol_used) {
            if (mat_mean(mat_index_main(i)) > tol_used) {
              mat_cv(mat_index_main(i)) = mat_sd(mat_index_main(i)) / mat_mean(mat_index_main(i));
            } else {
              mat_cv(mat_index_main(i)) = 0.0;
            }
            
            for (int k = 0; k < mat_index_main_nonzeros; k++) {
              if (k % 10 == 0) Rcpp::checkUserInterrupt();
              
              mat_elems_forcorr.zeros();
              bool found_greater_forcorr {false};
              
              for (int l = 0; l < refmatnum; l++) {
                double elem_check_mef = as<arma::mat>(ref_byyear(l))(mat_index_main(k));
                if (elem_check_mef > tol_used) {
                  mat_elems_forcorr(l) = elem_check_mef;
                  found_greater_forcorr = true;
                }
              }
              
              double mefc_sd {0.0};
              if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
              
              if (mefc_sd > 0.0) {
                arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
                double corr_bit = new_cor(0);
                if (corr_bit > tol_used) mat_corr(mat_index_main(i), mat_index_main(k)) = corr_bit;
              }
            }
          }
        }
        
        ref_matmean = mat_mean;
        ref_matelas = elas3matrix(mat_mean, false);
        ref_matcv = mat_cv;
        ref_matcorr = mat_corr;
        
      } else {
        arma::sp_mat mat_mean(matdim, matdim);
        arma::sp_mat mat_sd(matdim, matdim);
        arma::sp_mat mat_cv (matdim, matdim);
        arma::sp_mat mat_corr (matdim * matdim, matdim * matdim);
        
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) /
            static_cast<double>(refmatnum));
        }
        
        List eigenstuff = LefkoMats::decomp3sp_inp(mat_mean);
        arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
        double lambda = max(realeigenvals);
        rvals_ref = log(lambda);
        
        for (int i = 0; i < mat_index_main_nonzeros; i++) {
          arma::vec mat_elems(refmatnum, fill::zeros);
          arma::vec mat_elems_forcorr (refmatnum, fill::zeros);
          
          bool found_greater {false};
          
          for (int j = 0; j < refmatnum; j++) {
            double elem_check_mat_elems = as<arma::sp_mat>(ref_byyear(j))(mat_index_main(i));
            if (elem_check_mat_elems > tol_used) {
              mat_elems(j) = elem_check_mat_elems;
              found_greater = true;
            }
          }
          
          if (found_greater) mat_sd(mat_index_main(i)) = arma::stddev(mat_elems, 0);
          
          if (mat_sd(mat_index_main(i)) > tol_used) {
            if (mat_mean(mat_index_main(i)) > tol_used) {
              mat_cv(mat_index_main(i)) = mat_sd(mat_index_main(i)) / mat_mean(mat_index_main(i));
            } else {
              mat_cv(mat_index_main(i)) = 0.0;
            }
            
            // This is crazy - correlation matrix is ultra-high dimension...
            for (int k = 0; k < mat_index_main_nonzeros; k++) {
              if (k % 10 == 0) Rcpp::checkUserInterrupt();
              
              mat_elems_forcorr.zeros();
              bool found_greater_forcorr {false};
              
              for (int l = 0; l < refmatnum; l++) {
                double elem_check_mef = as<arma::sp_mat>(ref_byyear(l))(mat_index_main(k));
                if (elem_check_mef > tol_used) {
                  mat_elems_forcorr(l) = elem_check_mef;
                  found_greater_forcorr = true;
                }
              }
              
              double mefc_sd {0.0};
              if (found_greater_forcorr) mefc_sd = stddev(mat_elems_forcorr);
              
              if (mefc_sd > 0.0) {
                arma::mat new_cor = arma::cor(mat_elems, mat_elems_forcorr);
                double corr_bit = new_cor(0);
                if (corr_bit > tol_used) mat_corr(mat_index_main(i), mat_index_main(k)) = corr_bit;
              }
            }
          }
        }
        
        ref_sp_matmean = mat_mean;
        ref_sp_matelas = elas3sp_matrix(mat_mean);
        ref_sp_matcv = mat_cv;
        ref_sp_matcorr = mat_corr;
      }
    }
  }
  
  // Contributions
  List poppatch_contmeans (numpoppatches);
  List poppatch_diffelas (numpoppatches);
  List poppatch_diffcv (numpoppatches);
  List poppatch_diffcorr (numpoppatches);
  
  if (!sparse && !sparse_input) {
    arma::mat log_ref_matmean = arma::log(ref_matmean);
    arma::mat ref_cv_col = vectorise(ref_matcv);
    arma::mat big_ref_cv_col = ref_cv_col * ref_cv_col.t();
    arma::mat ref_elas_col = vectorise(ref_matelas);
    arma::mat big_ref_elas_col = ref_elas_col * ref_elas_col.t();
    
    for (int i = 0; i < numpoppatches; i++) {
      Rcpp::checkUserInterrupt();
      
      // Contributions of differences in element means
      arma::mat current_contmean = ref_matelas % // Originally (0.5 * (ref_matelas + as<arma::mat>(poppatch_elasmat(i))))
        (arma::log(as<arma::mat>(poppatch_meanmat(i))) - log_ref_matmean); // (log_ref_matmean - arma::log(as<arma::mat>(poppatch_meanmat(i)))) 
      current_contmean.elem(find_nonfinite(current_contmean)).zeros();
      
      poppatch_contmeans(i) = current_contmean;
      
      // Contributions of differences in matrix element elasticities
      arma::mat cv_col = vectorise(as<arma::mat>(poppatch_cvmat(i)));
      arma::mat big_cv_col = cv_col * cv_col.t();
      arma::mat elas_col = vectorise(as<arma::mat>(poppatch_elasmat(i)));
      arma::mat big_elas_col = elas_col * elas_col.t();
      arma::mat current_contelas = -0.5 * 0.5 * ((big_ref_cv_col % ref_matcorr) +
        (big_cv_col % as<arma::mat>(poppatch_corrmat(i)))) %
        (big_ref_elas_col - big_elas_col);
      
      poppatch_diffelas(i) = current_contelas;
      
      // Contributions of differences in CV
      arma::mat current_meanelas = (0.5 * (big_ref_elas_col + big_elas_col));
      arma::mat current_diffcv = -0.5 * current_meanelas % (big_ref_cv_col - big_cv_col) %
        (0.5 * (ref_matcorr + as<arma::mat>(poppatch_corrmat(i))));
      poppatch_diffcv(i) = current_diffcv;
      
      // Contributions of differences in correlations
      arma::mat current_diffcorr = -0.5 * current_meanelas % (0.5 * (big_ref_cv_col + big_cv_col)) %
        (ref_matcorr - as<arma::mat>(poppatch_corrmat(i)));
      poppatch_diffcorr(i) = current_diffcorr;
    }
  } else {
    arma::sp_mat log_ref_matmean = LefkoUtils::spmat_log(ref_sp_matmean);
    
    arma::sp_mat ref_cv_col = vectorise(ref_sp_matcv);
    arma::sp_mat big_ref_cv_col = ref_cv_col * ref_cv_col.t();
    arma::sp_mat ref_elas_col = vectorise(ref_sp_matelas);
    arma::sp_mat big_ref_elas_col = ref_elas_col * ref_elas_col.t();
    
    for (int i = 0; i < numpoppatches; i++) {
      Rcpp::checkUserInterrupt();
      
      // Contributions of differences in element means
      arma::sp_mat log_spizzle= LefkoUtils::spmat_log(as<arma::sp_mat>(poppatch_meanmat(i)));
      
      arma::sp_mat current_contmean = ref_sp_matelas % (log_spizzle - log_ref_matmean);  // Originally (0.5 * (ref_sp_matelas + as<arma::sp_mat>(poppatch_elasmat(i)))) % (log_ref_matmean - log_spizzle)
      //current_contmean.elem(find_nonfinite(current_contmean)).zeros();
      
      poppatch_contmeans(i) = current_contmean;
      
      // Contributions of differences in matrix element elasticities
      arma::sp_mat cv_col = vectorise(as<arma::sp_mat>(poppatch_cvmat(i)));
      arma::sp_mat big_cv_col = cv_col * cv_col.t();
      arma::sp_mat elas_col = vectorise(as<arma::sp_mat>(poppatch_elasmat(i)));
      arma::sp_mat big_elas_col = elas_col * elas_col.t();
      arma::sp_mat current_contelas = -0.5 * 0.5 * ((big_ref_cv_col % ref_sp_matcorr) +
        (big_cv_col % as<arma::sp_mat>(poppatch_corrmat(i)))) %
        (big_ref_elas_col - big_elas_col);
      
      poppatch_diffelas(i) = current_contelas;
      
      // Contributions of differences in CV
      arma::sp_mat current_meanelas = (0.5 * (big_ref_elas_col + big_elas_col));
      arma::sp_mat current_diffcv = -0.5 * current_meanelas % (big_ref_cv_col - big_cv_col) %
        (0.5 * (ref_sp_matcorr + as<arma::sp_mat>(poppatch_corrmat(i))));
      
      poppatch_diffcv(i) = current_diffcv;
      
      // Contributions of differences in correlations
      arma::sp_mat current_diffcorr = -0.5 * current_meanelas % (0.5 * (big_ref_cv_col + big_cv_col)) %
        (ref_sp_matcorr - as<arma::sp_mat>(poppatch_corrmat(i)));
      
      poppatch_diffcorr(i) = current_diffcorr;
    }
  }
  
  List output_now(6);
  output_now(0) = poppatch_contmeans;
  output_now(1) = poppatch_diffelas;
  output_now(2) = poppatch_diffcv;
  output_now(3) = poppatch_diffcorr;
  output_now(4) = rvals_poppatch;
  output_now(5) = rvals_ref;

  StringVector output_names = {"cont_mean", "cont_elas", "cont_cv", "cont_corr",
    "r_values_m", "r_value_ref"};
  output_now.attr("names") = output_names;
  StringVector output_class = {"lefkoLTRE"};
  output_now.attr("class") = output_class;
  
  return output_now;
}

//' Creates Vector of Times Based on First-Order Markov Transition Matrix
//' 
//' Creates a vector of randomly sampled years / times to be used in projection.
//' Random sampling requires a 1st order Markovian transition matrix, showing
//' the probability of transitioning to each time from each time. Note that this
//' function is not required if the probability of transitioning to a particular
//' time does not vary with time.
//' 
//' @name markov_run
//' 
//' @param main_times An integer vector giving the years / times to use.
//' @param mat A matrix giving the transition probabilities from each time to
//' each time. Must have the same number of columns and rows as there are
//' elements in vector \code{times}.
//' @param times The number of times to project forward. Defaults to 10000.
//' @param start The start time to use. Defaults to the first time in vector
//' \code{main_times}.
//' 
//' @return An integer vector giving the order of times / years to use in
//' projection. This can be used as input in the \code{year} option in
//' functions \code{\link{projection3}()} and \code{\link{f_projection3}()}.
//' 
//' @seealso \code{\link{projection3}()}
//' @seealso \code{\link{f_projection3}()}
//' 
//' @examples
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
//'   year = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"), 
//'   supplement = cypsupp3r, yearcol = "year2", indivcol = "individ")
//' 
//' used_years <-c(2005, 2006, 2007, 2008)
//' 
//' yr_tx_vec <- c(0.4, 0.2, 0.2, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.4, 0.2,
//'   0.2, 0.2, 0.2, 0.4)
//' yr_tx_mat <- matrix(yr_tx_vec, 4, 4)
//' 
//' set.seed(1)
//' cyp_markov_vec_1 <- markov_run(main_times = used_years, mat = yr_tx_mat,
//'   times = 100)
//' 
//' set.seed(2)
//' cyp_markov_vec_2 <- markov_run(main_times = used_years, mat = yr_tx_mat,
//'   times = 100)
//' 
//' set.seed(3)
//' cyp_markov_vec_3 <- markov_run(main_times = used_years, mat = yr_tx_mat,
//'   times = 100)
//' 
//' cypstoch_1 <- projection3(cypmatrix3r, nreps = 1, times = 100,
//'   year = cyp_markov_vec_1)
//' cypstoch_2 <- projection3(cypmatrix3r, nreps = 1, times = 100,
//'   year = cyp_markov_vec_2)
//' cypstoch_3 <- projection3(cypmatrix3r, nreps = 1, times = 100,
//'   year = cyp_markov_vec_3)
//' 
//' @export markov_run
// [[Rcpp::export(markov_run)]]
Rcpp::IntegerVector markov_run(Rcpp::IntegerVector main_times,
  Rcpp::NumericMatrix mat, int times = 10000,
  Nullable<IntegerVector> start = R_NilValue) {
  
  int start_time {0};
  bool start_time_found {false};
  int start_position {0};
  
  int mat_rows = mat.nrow();
  int mat_cols = mat.ncol();
  int main_times_length = main_times.length();
  
  if (mat_rows != mat_cols) throw Rcpp::exception("Input matrix must be square.", false);
  if (mat_rows != main_times_length) {
    throw Rcpp::exception("Input matrix must have same dimensions as length of vector main_times.", 
      false);
  }
  
  if (start.isNotNull()) {
    IntegerVector start_entered(start);
    
    if (start_entered.length() != 1) pop_error("start", "a single integer", "", 1);
    
    for (int i = 0; i < main_times_length; i++) {
      if (main_times(i) == start_entered(0)) {
        start_time_found = true;
        start_position = i;
        start_time = main_times(i);
      }
    }
    
    if (!start_time_found && start_entered(0) == 0) {
      start_time_found = true;
      start_position = 0;
      start_time = main_times(0);
    }
  } else {
    start_time_found = true;
    start_position = 0;
    start_time = main_times(0);
  }
  
  if (!start_time_found) {
    throw Rcpp::exception("Argument main_times does not include input start_time values",
      false);
  }
  
  // Matrix standardization
  NumericVector mat_colsums (mat_rows);
  for (int i = 0; i < mat_cols; i++) {
    for (int j = 0; j < mat_rows; j++) {
      if (j == 0) mat_colsums(i) = 0.0;
      mat_colsums(i) += mat(j, i);
    }
  }
  
  for (int i = 0; i < mat_cols; i++) {
    for (int j = 0; j < mat_rows; j++) {
      mat(j, i) = mat(j, i) / mat_colsums(i);
    }
  }
  
  // Creation of output vector
  IntegerVector out_vector (times);
  IntegerVector time_pulled_vec;
  int time_pulled;
  IntegerVector choices = Rcpp::seq(0, (main_times_length - 1));
  
  out_vector(0) = start_time;
  for (int i = 1; i < times; i++) {
    NumericVector mat_col_used = mat(_, start_position);
    arma::vec arma_mat_col_used = arma::mat(mat_col_used);
    time_pulled_vec = Rcpp::RcppArmadillo::sample(choices, 1, true, mat_col_used);
    time_pulled = time_pulled_vec(0);
    
    out_vector(i) = main_times(time_pulled);
    start_position = time_pulled;
  }
  
  return (out_vector);
}

