#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;



// Index of functions
// 
// 1. List sf_reassess  Standardize Stageframe For MPM Analysis
// 2. DataFrame sf_skeleton  Create Skeleton Stageframe
// 3. List thefifthhousemate  Create Historically Structured Version of ahMPM
// 4. hist_null  Create Historical MPMs Assuming No Influence of Individual History
// 5. List lmean  Estimate Mean Projection Matrices
// 6. List add_stage  Add a New Stage to an Existing LefkoMat Object
// 7. List cycle_check  Check Continuity of Life Cycle through Matrices in lefkoMat Objects


//' Standardize Stageframe For MPM Analysis
//' 
//' Function \code{sf_reassess()} takes a stageframe as input, and uses
//' information supplied there and through the supplement, reproduction and
//' overwrite tables to rearrange this into a format usable by the matrix
//' creation functions, \code{mpm_create()}, \code{flefko3()},
//' \code{flefko2()}, \code{aflefko2()}, \code{rlefko3()}, and \code{rlefko2()}.
//' This is performed through a call to \code{sf_reassess_internal()}.
//' 
//' @name .sf_reassess
//' 
//' @param stageframe The original stageframe.
//' @param supplement The original supplemental data input (class
//' \code{lefkoSD}). Can also equal NA.
//' @param overwrite An overwrite table.
//' @param repmatrix The original reproduction matrix. Can also equal \code{NA},
//' \code{0}, or \code{NULL} (the last value by default).
//' @param agemat A logical value indicating whether MPM is age-by-stage.
//' @param historical A logical value indicating whether MPM is historical.
//' @param format An integer indicating whether matrices will be in Ehrlen
//' format (if set to 1), or deVries format (if set to 2). Setting to deVries
//' format adds one extra stage to account for the prior status of newborns.
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
// [[Rcpp::export(.sf_reassess)]]
Rcpp::List sf_reassess(const DataFrame& stageframe,
  Nullable<DataFrame> supplement = R_NilValue,
  Nullable<DataFrame> overwrite = R_NilValue,
  Nullable<NumericMatrix> repmatrix = R_NilValue,
  bool agemat = false, bool historical = false, int format = 1) {
  
  List sf_output = LefkoMats::sf_reassess_internal(stageframe, supplement,
    overwrite, repmatrix, agemat, historical, format);
  
  return sf_output;
}

//' Create Skeleton Stageframe
//' 
//' Function \code{sf_skeleton()} creates a skeleton \code{stageframe} object.
//' 
//' @name sf_skeleton
//' 
//' @param stages The number of stages, as an integer.
//' @param standard A logical value indicating whether to create a standard
//' \code{stageframe} object (\code{TRUE}, the default), or a reassessed
//' \code{stageframe} object as created by function \code{mpm_create()}
//' (\code{FALSE}).
//' 
//' @return A data frame of class \code{stageframe}.
//' 
//' @export
// [[Rcpp::export(sf_skeleton)]]
Rcpp::DataFrame sf_skeleton(int stages, bool standard = true) { 
  
  bool reassessed = false;
  if (!standard) reassessed = true;
  
  DataFrame new_sf = LefkoMats::sf_core(stages, reassessed, false);
  
  return new_sf;
}

//' Create Element Index for Matrix Estimation
//' 
//' Function \code{simplepizzle()} creates a data frame object used by function
//' \code{\link{hist_null}()} to provide an index for estimation of null
//' historical matrices from ahistorical MPM inputs.
//' 
//' @name simplepizzle
//' 
//' @param StageFrame The stageframe object identifying the life history model
//' being operationalized.
//' @param format Integer indicating whether historical matrices should be in
//' (1) Ehrlen or (2) deVries format.
//' 
//' @return The output is composed of three elements:
//' \item{ahstages}{A new stageframe, which only differs from the input
//' stageframe in deVries format.}
//' \item{hstages}{A new historical stage-pair index for the new historical
//' matrices.}
//' \item{allstages}{A large data frame describing every element to be estimated
//' in the new historical matrices.}
//' 
//' @keywords internal
//' @noRd
Rcpp::List simplepizzle(DataFrame StageFrame, int format) {
  
  arma::vec newstageid = as<arma::vec>(StageFrame["stage_id"]);
  StringVector origstageid = as<StringVector>(StageFrame["stage"]);
  arma::vec binsizectr = as<arma::vec>(StageFrame["sizebin_center"]);
  arma::vec repstatus = as<arma::vec>(StageFrame["repstatus"]);
  arma::vec obsstatus = as<arma::vec>(StageFrame["obsstatus"]);
  arma::vec immstatus = as<arma::vec>(StageFrame["immstatus"]);
  arma::vec matstatus = as<arma::vec>(StageFrame["matstatus"]);
  arma::vec indata = as<arma::vec>(StageFrame["indataset"]);
  arma::vec binsizewidth = as<arma::vec>(StageFrame["sizebin_width"]);
  arma::vec minage = as<arma::vec>(StageFrame["min_age"]);
  arma::vec maxage = as<arma::vec>(StageFrame["max_age"]);
  arma::vec group = as<arma::vec>(StageFrame["group"]);
  arma::vec entrystage = as<arma::vec>(StageFrame["entrystage"]);
  
  arma::vec binsizebctr = as<arma::vec>(StageFrame["sizebinb_center"]);
  arma::vec binsizecctr = as<arma::vec>(StageFrame["sizebinc_center"]);
  arma::vec binsizebwidth = as<arma::vec>(StageFrame["sizebinb_width"]);
  arma::vec binsizecwidth = as<arma::vec>(StageFrame["sizebinc_width"]);
  
  // This section determines the length of the matrix map data frame
  int nostages = static_cast<int>(newstageid.n_elem);
  int nostages_nounborn = nostages;
  int totallength {0};
  if (format == 2)  nostages = nostages + 1;
  
  arma::vec almostborn (nostages, fill::zeros);
  
  if (format == 2) {
    arma::vec oldorigsize = as<arma::vec>(StageFrame["original_size"]);
    arma::vec oldorigbsize = as<arma::vec>(StageFrame["original_size_b"]);
    arma::vec oldorigcsize = as<arma::vec>(StageFrame["original_size_c"]);
    arma::vec oldpropstatus = as<arma::vec>(StageFrame["propstatus"]);
    arma::vec oldbinhalfwidthraw = as<arma::vec>(StageFrame["binhalfwidth_raw"]);
    arma::vec oldsizebinmin = as<arma::vec>(StageFrame["sizebin_min"]);
    arma::vec oldsizebinmax = as<arma::vec>(StageFrame["sizebin_max"]);
    arma::vec oldbinhalfwidthbraw = as<arma::vec>(StageFrame["binhalfwidthb_raw"]);
    arma::vec oldsizebinbmin = as<arma::vec>(StageFrame["sizebinb_min"]);
    arma::vec oldsizebinbmax = as<arma::vec>(StageFrame["sizebinb_max"]);
    arma::vec oldbinhalfwidthcraw = as<arma::vec>(StageFrame["binhalfwidthc_raw"]);
    arma::vec oldsizebincmin = as<arma::vec>(StageFrame["sizebinc_min"]);
    arma::vec oldsizebincmax = as<arma::vec>(StageFrame["sizebinc_max"]);
    Rcpp::StringVector oldcomments = as<StringVector>(StageFrame["comments"]);
    arma::vec oldentrystage = as<arma::vec>(StageFrame["entrystage"]);
    
    nostages_nounborn = nostages - 1;
    totallength = (2 * nostages_nounborn * nostages_nounborn *
      nostages);
    
    arma::vec newstageidvec(nostages, fill::zeros);
    Rcpp::StringVector newstagevec(nostages);
    arma::vec neworigsizevec(nostages, fill::zeros);
    arma::vec neworigsizebvec(nostages, fill::zeros);
    arma::vec neworigsizecvec(nostages, fill::zeros);
    arma::vec newminagevec(nostages, fill::zeros);
    arma::vec newmaxagevec(nostages, fill::zeros);
    arma::vec newrepstatusvec(nostages, fill::zeros);
    arma::vec newobsstatusvec(nostages);
    arma::vec newpropstatusvec(nostages, fill::zeros);
    arma::vec newimmstatusvec(nostages, fill::zeros);
    arma::vec newmatstatusvec(nostages, fill::zeros);
    arma::vec newindatasetvec(nostages, fill::zeros);
    arma::vec newbinhalfwidthrawvec(nostages, fill::zeros);
    arma::vec newsizebinminvec(nostages, fill::zeros);
    arma::vec newsizebinmaxvec(nostages, fill::zeros);
    arma::vec newsizebincentervec(nostages, fill::zeros);
    arma::vec newsizebinwidthvec(nostages, fill::zeros);
    
    arma::vec newbinhalfwidthbrawvec(nostages, fill::zeros);
    arma::vec newsizebinbminvec(nostages, fill::zeros);
    arma::vec newsizebinbmaxvec(nostages, fill::zeros);
    arma::vec newsizebinbcentervec(nostages, fill::zeros);
    arma::vec newsizebinbwidthvec(nostages, fill::zeros);
    
    arma::vec newbinhalfwidthcrawvec(nostages, fill::zeros);
    arma::vec newsizebincminvec(nostages, fill::zeros);
    arma::vec newsizebincmaxvec(nostages, fill::zeros);
    arma::vec newsizebinccentervec(nostages, fill::zeros);
    arma::vec newsizebincwidthvec(nostages, fill::zeros);
    
    arma::vec newgroupvec(nostages, fill::zeros);
    Rcpp::StringVector newcomments(nostages);
    arma::vec newentrystage(nostages, fill::zeros);
    arma::vec newalmostborn(nostages, fill::zeros);
    
    for (int i = 0; i < nostages_nounborn; i++) {
      newstageidvec(i) = newstageid(i);
      newstagevec(i) = origstageid(i);
      neworigsizevec(i) = oldorigsize(i);
      neworigsizebvec(i) = oldorigbsize(i);
      neworigsizecvec(i) = oldorigcsize(i);
      newminagevec(i) = minage(i);
      newmaxagevec(i) = maxage(i);
      newrepstatusvec(i) = repstatus(i);
      newobsstatusvec(i) = obsstatus(i);
      newpropstatusvec(i) = oldpropstatus(i);
      newimmstatusvec(i) = immstatus(i);
      newmatstatusvec(i) = matstatus(i);
      newindatasetvec(i) = indata(i);
      
      newbinhalfwidthrawvec(i) = oldbinhalfwidthraw(i);
      newsizebinminvec(i) = oldsizebinmin(i);
      newsizebinmaxvec(i) = oldsizebinmax(i);
      newsizebincentervec(i) = binsizectr(i);
      newsizebinwidthvec(i) = binsizewidth(i);
      
      newbinhalfwidthbrawvec(i) = oldbinhalfwidthbraw(i);
      newsizebinbminvec(i) = oldsizebinbmin(i);
      newsizebinbmaxvec(i) = oldsizebinbmax(i);
      newsizebinbcentervec(i) = binsizebctr(i);
      newsizebinbwidthvec(i) = binsizebwidth(i);
      
      newbinhalfwidthcrawvec(i) = oldbinhalfwidthcraw(i);
      newsizebincminvec(i) = oldsizebincmin(i);
      newsizebincmaxvec(i) = oldsizebincmax(i);
      newsizebinccentervec(i) = binsizecctr(i);
      newsizebincwidthvec(i) = binsizecwidth(i);
      
      newgroupvec(i) = group(i);
      newcomments(i) = oldcomments(i);
      newentrystage(i) = oldentrystage(i);
    }
    
    newstageidvec(nostages_nounborn) = newstageid(nostages_nounborn - 1.0) + 1.0;
    newstagevec(nostages_nounborn) = "AlmostBorn";
    neworigsizevec(nostages_nounborn) = 0.0;
    neworigsizebvec(nostages_nounborn) = 0.0;
    neworigsizecvec(nostages_nounborn) = 0.0;
    newminagevec(nostages_nounborn) = NA_REAL;
    newmaxagevec(nostages_nounborn) = NA_REAL;
    newrepstatusvec(nostages_nounborn) = 0.0;
    newobsstatusvec(nostages_nounborn) = 1.0;
    newpropstatusvec(nostages_nounborn) = 0.0;
    newimmstatusvec(nostages_nounborn) = 1.0;
    newmatstatusvec(nostages_nounborn) = 0.0;
    newindatasetvec(nostages_nounborn) = 1.0;
    
    newbinhalfwidthrawvec(nostages_nounborn) = 0.0;
    newsizebinminvec(nostages_nounborn) = 0.0;
    newsizebinmaxvec(nostages_nounborn) = 0.0;
    newsizebincentervec(nostages_nounborn) = 0.0;
    newsizebinwidthvec(nostages_nounborn) = 0.0;
    
    newbinhalfwidthbrawvec(nostages_nounborn) = 0.0;
    newsizebinbminvec(nostages_nounborn) = 0.0;
    newsizebinbmaxvec(nostages_nounborn) = 0.0;
    newsizebinbcentervec(nostages_nounborn) = 0.0;
    newsizebinbwidthvec(nostages_nounborn) = 0.0;
    
    newbinhalfwidthcrawvec(nostages_nounborn) = 0.0;
    newsizebincminvec(nostages_nounborn) = 0.0;
    newsizebincmaxvec(nostages_nounborn) = 0.0;
    newsizebinccentervec(nostages_nounborn) = 0.0;
    newsizebincwidthvec(nostages_nounborn) = 0.0;
    
    newgroupvec(nostages_nounborn) = 0.0;
    newcomments(nostages_nounborn) = "Almost Born";
    newentrystage(nostages_nounborn) = 0.0;
    newalmostborn(nostages_nounborn) = 1.0;
    
    Rcpp::List new_stageframe(32);
    
    new_stageframe(0) = newstageidvec;
    new_stageframe(1) = newstagevec;
    new_stageframe(2) = neworigsizevec;
    new_stageframe(3) = neworigsizebvec;
    new_stageframe(4) = neworigsizecvec;
    new_stageframe(5) = newminagevec;
    new_stageframe(6) = newmaxagevec;
    new_stageframe(7) = newrepstatusvec;
    new_stageframe(8) = newobsstatusvec;
    new_stageframe(9) = newpropstatusvec;
    new_stageframe(10) = newimmstatusvec;
    new_stageframe(11) = newmatstatusvec;
    new_stageframe(12) = newindatasetvec;
    
    new_stageframe(13) = newbinhalfwidthrawvec;
    new_stageframe(14) = newsizebinminvec;
    new_stageframe(15) = newsizebinmaxvec;
    new_stageframe(16) = newsizebincentervec;
    new_stageframe(17) = newsizebinwidthvec;
    
    new_stageframe(18) = newbinhalfwidthbrawvec;
    new_stageframe(19) = newsizebinbminvec;
    new_stageframe(20) = newsizebinbmaxvec;
    new_stageframe(21) = newsizebinbcentervec;
    new_stageframe(22) = newsizebinbwidthvec;
    
    new_stageframe(23) = newbinhalfwidthcrawvec;
    new_stageframe(24) = newsizebincminvec;
    new_stageframe(25) = newsizebincmaxvec;
    new_stageframe(26) = newsizebinccentervec;
    new_stageframe(27) = newsizebincwidthvec;
    
    new_stageframe(28) = newgroupvec;
    new_stageframe(29) = newcomments;
    new_stageframe(30) = newentrystage;
    new_stageframe(31) = newalmostborn;
    
    CharacterVector sfnamevec = {"stage_id", "stage", "original_size", "size_b",
      "size_c", "min_age", "max_age", "repstatus", "obsstatus", "propstatus",
      "immstatus", "matstatus", "indataset", "binhalfwidth_raw", "sizebin_min",
      "sizebin_max", "sizebin_center", "sizebin_width", "binhalfwidthb_raw",
      "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
      "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max", "sizebinc_center",
      "sizebinc_width", "group", "comments", "entrystage", "almostborn"};
    
    new_stageframe.attr("names") = sfnamevec;
    new_stageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      static_cast<int>(newstageidvec.n_elem));
    new_stageframe.attr("class") = "data.frame";
    
    StageFrame = new_stageframe;
    
    newstageid = as<arma::vec>(new_stageframe["stage_id"]);
    origstageid = as<StringVector>(new_stageframe["stage"]);
    binsizectr = as<arma::vec>(new_stageframe["sizebin_center"]);
    repstatus = as<arma::vec>(new_stageframe["repstatus"]);
    obsstatus = as<arma::vec>(new_stageframe["obsstatus"]);
    immstatus = as<arma::vec>(new_stageframe["immstatus"]);
    matstatus = as<arma::vec>(new_stageframe["matstatus"]);
    indata = as<arma::vec>(new_stageframe["indataset"]);
    binsizewidth = as<arma::vec>(new_stageframe["sizebin_width"]);
    minage = as<arma::vec>(new_stageframe["min_age"]);
    maxage = as<arma::vec>(new_stageframe["max_age"]);
    group = as<arma::vec>(new_stageframe["group"]);
    
    binsizebctr = as<arma::vec>(new_stageframe["sizebinb_center"]);
    binsizecctr = as<arma::vec>(new_stageframe["sizebinc_center"]);
    binsizebwidth = as<arma::vec>(new_stageframe["sizebinb_width"]);
    binsizecwidth = as<arma::vec>(new_stageframe["sizebinc_width"]);
    almostborn = as<arma::vec>(new_stageframe["almostborn"]);
    entrystage = as<arma::vec>(new_stageframe["entrystage"]);
    
  } else {
    totallength = (nostages * nostages * nostages);
  }
  
  // New let's create the new hstages
  Rcpp::List hstages(4);
  
  int hstages_length = nostages * nostages_nounborn;
  
  Rcpp::IntegerVector stid2(hstages_length);
  Rcpp::IntegerVector stid1(hstages_length);
  Rcpp::StringVector st2(hstages_length);
  Rcpp::StringVector st1(hstages_length);
  
  for (int j = 0; j < nostages; j++) {
    for (int i = 0; i < nostages_nounborn; i++) {
      stid2((j * (nostages_nounborn)) + i) = i + 1;
      stid1((j * (nostages_nounborn)) + i) = j + 1;
      st2((j * (nostages_nounborn)) + i) = origstageid(i);
      st1((j * (nostages_nounborn)) + i) = origstageid(j);
    }
  }
  
  hstages(0) = stid2;
  hstages(1) = stid1;
  hstages(2) = st2;
  hstages(3) = st1;
  
  CharacterVector hsnamevec = {"stage_id_2", "stage_id_1", "stage_2", "stage_1"};
  hstages.attr("names") = hsnamevec;
  hstages.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
    stid2.length());
  hstages.attr("class") = "data.frame";
  
  // Set up vectors to put together into matrix map data frame
  arma::vec stage3(totallength, fill::zeros);
  arma::vec stage2n(totallength, fill::zeros);
  arma::vec stage2o(totallength, fill::zeros);
  arma::vec stage1(totallength, fill::zeros);
  
  arma::vec size3(totallength, fill::zeros);
  arma::vec size2n(totallength, fill::zeros);
  arma::vec size2o(totallength, fill::zeros);
  arma::vec size1(totallength, fill::zeros);
  
  arma::vec sizeb3(totallength, fill::zeros);
  arma::vec sizeb2n(totallength, fill::zeros);
  arma::vec sizeb2o(totallength, fill::zeros);
  arma::vec sizeb1(totallength, fill::zeros);
  
  arma::vec sizec3(totallength, fill::zeros);
  arma::vec sizec2n(totallength, fill::zeros);
  arma::vec sizec2o(totallength, fill::zeros);
  arma::vec sizec1(totallength, fill::zeros);
  
  arma::vec obs3(totallength, fill::zeros);
  arma::vec obs2n(totallength, fill::zeros);
  arma::vec obs2o(totallength, fill::zeros);
  arma::vec obs1(totallength, fill::zeros);
  
  arma::vec rep3(totallength, fill::zeros);
  arma::vec rep2n(totallength, fill::zeros);
  arma::vec rep2o(totallength, fill::zeros);
  arma::vec rep1(totallength, fill::zeros);
  
  arma::vec mat3(totallength, fill::zeros);
  arma::vec mat2n(totallength, fill::zeros);
  arma::vec mat2o(totallength, fill::zeros);
  arma::vec mat1(totallength, fill::zeros);
  
  arma::vec imm3(totallength, fill::zeros);
  arma::vec imm2n(totallength, fill::zeros);
  arma::vec imm2o(totallength, fill::zeros);
  arma::vec imm1(totallength, fill::zeros);
  
  arma::vec repentry3(totallength, fill::zeros);
  arma::vec repentry2o(totallength, fill::zeros);
  arma::vec almostborn2n(totallength, fill::zeros);
  arma::vec almostborn1(totallength, fill::zeros);
  
  arma::vec binwidth(totallength, fill::zeros);
  arma::vec binbwidth(totallength, fill::zeros);
  arma::vec bincwidth(totallength, fill::zeros);
  
  arma::vec indata3(totallength, fill::zeros);
  arma::vec indata2n(totallength, fill::zeros);
  arma::vec indata2o(totallength, fill::zeros);
  arma::vec indata1(totallength, fill::zeros);
  
  arma::vec minage3(totallength, fill::zeros);
  arma::vec minage2(totallength, fill::zeros);
  arma::vec maxage3(totallength, fill::zeros);
  arma::vec maxage2(totallength, fill::zeros);
  
  arma::vec grp3(totallength, fill::zeros);
  arma::vec grp2n(totallength, fill::zeros);
  arma::vec grp2o(totallength, fill::zeros);
  arma::vec grp1(totallength, fill::zeros);
  
  arma::vec actualage(totallength, fill::zeros);
  arma::vec index321u(totallength);
  arma::vec index321f(totallength);
  arma::vec index21(totallength);
  arma::vec indatalong(totallength, fill::zeros);
  arma::vec aliveequal(totallength);
  arma::vec included(totallength, fill::zeros);
  index321u.fill(-1);
  index321f.fill(-1);
  index21.fill(-1);
  aliveequal.fill(-1);
  
  long long int currentindex {0};
  
  // Main data frame creation loops
  // If style = 0, this will create AllStages for the historical case
  if (format == 2) {
    for (int time1 = 0; time1 < nostages; time1++) {
      for (int time2o = 0; time2o < nostages_nounborn; time2o++) {
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            if (time3 != nostages_nounborn) {
              if (time2n == time2o || time2n == nostages_nounborn) {
                
                included(currentindex) = 1;
                
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
                
                if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0;
                if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0;
                if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0;
                if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0;
                
                sizec3(currentindex) = binsizecctr(time3);
                sizec2n(currentindex) = binsizecctr(time2n);
                sizec2o(currentindex) = binsizecctr(time2o);
                sizec1(currentindex) = binsizecctr(time1);
                
                if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0;
                if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0;
                if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0;
                if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0;
                
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
                
                repentry3(currentindex) = entrystage(time3);
                repentry2o(currentindex) = entrystage(time2o);
                almostborn2n(currentindex) = almostborn(time2n);
                almostborn1(currentindex) = almostborn(time1);
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2o);
                indata1(currentindex) = indata(time1);
                
                binwidth(currentindex) = binsizewidth(time3);
                binbwidth(currentindex) = binsizebwidth(time3);
                bincwidth(currentindex) = binsizecwidth(time3);
                
                if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0;
                if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0;
                
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2o);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2o);
                actualage(currentindex) = 0;
                
                grp3(currentindex) = group(time3);
                grp2n(currentindex) = group(time2n);
                grp2o(currentindex) = group(time2o);
                grp1(currentindex) = group(time1);
                
                aliveequal(currentindex) = (stage3(currentindex) - 1) + 
                  ((stage2n(currentindex) - 1) * nostages_nounborn) + 
                  ((stage2o(currentindex) - 1) * nostages * nostages_nounborn) + 
                  ((stage1(currentindex) - 1) * nostages_nounborn * nostages *
                    nostages_nounborn);
                
                index321u(currentindex) = (stage3(currentindex) - 1) + 
                  ((stage2n(currentindex) - 1) * nostages_nounborn) + 
                  ((stage2o(currentindex) - 1) * nostages * nostages_nounborn) + 
                  ((stage1(currentindex) - 1) * nostages_nounborn * nostages *
                    nostages_nounborn);
                index321f(currentindex) = index321u(currentindex);
                  
                index21(currentindex) = (stage3(currentindex) - 1) + 
                  ((stage2o(currentindex) - 1) * nostages_nounborn); // Used to be 2o and 1
                
                indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
                  indata2o(currentindex) * indata1(currentindex);
                
                // Now a few corrections
                if ((almostborn1(currentindex) > 0.0 && repentry2o(currentindex) < 1.0)) {
                    index321u(currentindex) = -1.0;
                    index321f(currentindex) = -1.0;
                }
                if (almostborn2n(currentindex) > 0.0) {
                  if (repentry3(currentindex) < 1.0 || rep2o(currentindex) < 1.0) {
                    index321u(currentindex) = -1.0;
                    index321f(currentindex) = -1.0;
                  }
                }
                if (time2n == time2o) {
                  index321f(currentindex) = -1.0;
                }
                currentindex += 1;
              } // if (time2n == tim2o || time2n == nostages_nounborn) statement
            } // if (time3n != dead_stage) statement
          } // time3 loop
        } // time2n loop
      } // time2o loop
    } // time1 loop 
    
  } else if (format == 1) { // Historical MPM in Ehrlen format
    for (int time1 = 0; time1 < nostages; time1++) {
      for (int time2o = 0; time2o < nostages; time2o++) {
        for (int time3 = 0; time3 < nostages; time3++) {
          
          included(currentindex) = 1;
          
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
          
          if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0;
          if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0;
          if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0;
          if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0;
                
          sizec3(currentindex) = binsizecctr(time3);
          sizec2n(currentindex) = binsizecctr(time2o);
          sizec2o(currentindex) = binsizecctr(time2o);
          sizec1(currentindex) = binsizecctr(time1);
          
          if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0;
          if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0;
          if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0;
          if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0;
          
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
          
          repentry3(currentindex) = 0;
          repentry2o(currentindex) = entrystage(time2o);
          almostborn1(currentindex) = almostborn(time1);
          
          indata3(currentindex) = indata(time3);
          indata2n(currentindex) = indata(time2o);
          indata2o(currentindex) = indata(time2o);
          indata1(currentindex) = indata(time1);
          
          binwidth(currentindex) = binsizewidth(time3);
          binbwidth(currentindex) = binsizebwidth(time3);
          bincwidth(currentindex) = binsizecwidth(time3);
          
          if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0;
          if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0;
          
          minage3(currentindex) = minage(time3);
          minage2(currentindex) = minage(time2o);
          maxage3(currentindex) = maxage(time3);
          maxage2(currentindex) = maxage(time2o);
          actualage(currentindex) = 0;
          
          grp3(currentindex) = group(time3);
          grp2n(currentindex) = group(time2o);
          grp2o(currentindex) = group(time2o);
          grp1(currentindex) = group(time1);
          
          aliveequal(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * 
              (nostages - 1)) + ((stage2o(currentindex) - 1) * (nostages - 1) * (nostages - 1)) + 
            ((stage1(currentindex) - 1) * (nostages - 1) * (nostages - 1) * (nostages - 1));
            
          index321u(currentindex) = (stage3(currentindex) - 1) + 
            ((stage2n(currentindex) - 1) * nostages) + 
            ((stage2n(currentindex) - 1) * nostages * nostages) + 
            ((stage1(currentindex) - 1) * nostages * nostages * nostages);
          index321f(currentindex) = index321u(currentindex);
          index21(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * nostages);
          
          indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
            indata2o(currentindex) * indata1(currentindex);
          
          currentindex += 1;
        } // time3 loop
      } // time2o loop
    } // time1 loop 
  }
  
  int stage3_length = static_cast<int>(stage3.n_elem);
  
  Rcpp::List output_longlist(63);
  
  output_longlist(0) = Rcpp::NumericVector(stage3.begin(), stage3.end());
  output_longlist(1) = Rcpp::NumericVector(stage2n.begin(), stage2n.end());
  output_longlist(2) = Rcpp::NumericVector(stage2o.begin(), stage2o.end());
  output_longlist(3) = Rcpp::NumericVector(stage1.begin(), stage1.end());
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
  output_longlist(50) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(51) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(52) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(53) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(54) = Rcpp::NumericVector(stage3_length, 1.0);
  output_longlist(55) = Rcpp::NumericVector(stage3_length, 1.0);
  
  output_longlist(56) = Rcpp::NumericVector(aliveequal.begin(), aliveequal.end());
  output_longlist(57) = Rcpp::NumericVector(index321u.begin(), index321u.end());
  output_longlist(58) = Rcpp::NumericVector(index21.begin(), index21.end());
  output_longlist(59) = Rcpp::NumericVector(repentry2o.begin(), repentry2o.end());
  output_longlist(60) = Rcpp::NumericVector(almostborn2n.begin(), almostborn2n.end());
  output_longlist(61) = Rcpp::NumericVector(almostborn1.begin(), almostborn1.end());
  output_longlist(62) = Rcpp::NumericVector(index321f.begin(), index321f.end());
  
  CharacterVector namevec = {"stage3", "stage2n", "stage2o", "stage1", "size3",
    "size2n", "size2o", "size1", "sizeb3", "sizeb2n", "sizeb2o", "sizeb1", 
    "sizec3", "sizec2n", "sizec2o", "sizec1", "obs3", "obs2n", "obs2o", "obs1",
    "rep3", "rep2n", "rep2o", "rep1", "mat3", "mat2n", "mat2o", "mat1", "imm3",
    "imm2n", "imm2o", "imm1", "repentry3", "indata3", "indata2n", "indata2o",
    "indata1", "binwidth", "binbwidth", "bincwidth", "minage3", "minage2",
    "maxage3", "maxage2", "actualage", "group3", "group2n", "group2o", "group1",
    "indata", "ovgiven_t", "ovest_t", "ovgiven_f", "ovest_f", "ovsurvmult",
    "ovfecmult", "aliveandequal", "index321u", "index21", "repentry2o",
    "almostborn2n", "almostborn1", "index321f"};
  output_longlist.attr("names") = namevec;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stage3_length);
  output_longlist.attr("class") = "data.frame";
  
  Rcpp::List output = Rcpp::List::create(Named("ahstages") = StageFrame,
    _["hstages"] = hstages, _["allstages"] = output_longlist);
  return output;
}

//' Create Historically Structured Version of ahMPM
//' 
//' Function \code{thefifthhousemate()} takes an ahistorical MPM as input, and
//' uses the \code{allstages} index to create a historically structured version
//' of it.
//' 
//' @name thefifthhousemate
//' 
//' @param mpm The original ahMPM, supplied as a \code{lefkoMat} object.
//' @param allstages The index dataframe named \code{allstages}, in the third
//' element of output developed by \code{simplepizzle()}.
//' @param hstages The index dataframe named \code{hstages}, in the second
//' element of output developed by \code{simplepizzle()}.
//' @param stageframe The ahistorical stageframe supplied by
//' \code{simplepizzle()}.
//' @param format Integer indicating whether historical matrices should be in
//' (1) Ehrlen or (2) deVries format.
//' 
//' @return This will return a list of lists. The first list is composed of all
//' new \code{A} matrices. The second list is composed of all new \code{U}
//' matrices. The third list is composed of all new \code{F} matrices.
//' 
//' @keywords internal
//' @noRd
Rcpp::List thefifthhousemate(List mpm, DataFrame allstages, DataFrame hstages,
  DataFrame stageframe, int format) {
  Rcpp::List old_Umats = as<List>(mpm["U"]);
  Rcpp::List old_Fmats = as<List>(mpm["F"]);
  
  Rcpp::IntegerVector stageid = stageframe["stage_id"];
  int nostages = stageid.length();
  int nocols = nostages * nostages;
  
  if (format == 2) nocols = (nostages -1) * nostages;
  
  Rcpp::IntegerVector old_index = allstages["index21"];
  Rcpp::IntegerVector new_indexu = allstages["index321u"];
  Rcpp::IntegerVector new_indexf = allstages["index321f"];
  arma::ivec stage2o_all = as<arma::ivec>(allstages["stage2o"]);
  arma::ivec stage1_all = as<arma::ivec>(allstages["stage1"]);
  arma::ivec hst_s_id_1 = as<arma::ivec>(hstages["stage_id_1"]);
  arma::ivec hst_s_id_2 = as<arma::ivec>(hstages["stage_id_2"]);
  
  int num_mats = old_Umats.length();
  int index_elems = new_indexu.length();
  
  Rcpp::List new_Umats(num_mats);
  Rcpp::List new_Fmats(num_mats);
  Rcpp::List new_Amats(num_mats);
  
  arma::mat new_U(nocols, nocols, fill::zeros);
  arma::mat new_F(nocols, nocols, fill::zeros);
  arma::mat new_A(nocols, nocols, fill::zeros);
  arma::mat old_U(nostages, nostages, fill::zeros);
  arma::mat old_F(nostages, nostages, fill::zeros);
  
  for (int i = 0; i < num_mats; i++) {
    new_U.zeros();
    new_F.zeros();
    new_A.zeros();
    old_U.zeros();
    old_F.zeros();
    
    old_U = as<arma::mat>(old_Umats(i));
    old_F = as<arma::mat>(old_Fmats(i));
    
    for (int j = 0; j < index_elems; j++) {
      if (new_indexu(j) > -1.0) {
        new_U(new_indexu(j)) = old_U(old_index(j));
      }
      if (new_indexf(j) > -1.0) {
        new_F(new_indexf(j)) = old_F(old_index(j));
      }
    }
    
    new_A = new_U + new_F;
    new_Umats(i) = new_U;
    new_Fmats(i) = new_F;
    new_Amats(i) = new_A;
  }
  
  Rcpp::List output = List::create(Named("A") = new_Amats, _["U"] = new_Umats,
    _["F"] = new_Fmats);
  return output;
}

//' Create Historical MPMs Assuming No Influence of Individual History
//' 
//' Function \code{hist_null()} uses ahistorical MPMs to create the equivalent
//' MPMs in the structure of historical MPMs. These MPMs have the same
//' dimensions and stage structure of hMPMs but assume no influence of
//' individual history, and so can be compared to actual hMPMs.
//' 
//' @name hist_null
//' 
//' @param mpm An ahistorical MPM of class \code{lefkoMat}.
//' @param format An integer stipulating whether historical matrices should be
//' produced in Ehrlen format (\code{1}) or deVries format (\code{2}).
//' @param err_check A logical value indicating whether to output the main index
//' data frames used to sort elements in the matrices.
//' 
//' @return An object of class \code{lefkoMat}, with the same list structure as
//' the input object, but with \code{A}, \code{U}, and \code{F} elements
//' replaced with lists of historically-structured matrices, and with element
//' \code{hstages} changed from \code{NA} to an index of stage pairs
//' corresponding to the rows and columns of the new matrices. If
//' \code{err_check = TRUE}, then a list of three data frames showing the values
//' used to determine matrix element index values is also exported.
//' 
//' @section Notes:
//' This function does not currently identify biologically impossible
//' transitions. Ahistorical transition values are placed in all theoretically
//' possible historical transitions.
//' 
//' @examples
//' sizevector <- c(1, 1, 2, 3)
//' stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
//' repvector <- c(0, 0, 1, 1)
//' obsvector <- c(1, 1, 1, 1)
//' matvector <- c(0, 1, 1, 1)
//' immvector <- c(1, 0, 0, 0)
//' propvector <- c(0, 0, 0, 0)
//' indataset <- c(1, 1, 1, 1)
//' binvec <- c(0.5, 0.5, 0.5, 0.5)
//' 
//' anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' # POPN C 2003-2004
//' XC3 <- matrix(c(0, 0, 1.74, 1.74,
//' 0.208333333, 0, 0, 0.057142857,
//' 0.041666667, 0.076923077, 0, 0,
//' 0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
//' 
//' # 2004-2005
//' XC4 <- matrix(c(0, 0, 0.3, 0.6,
//' 0.32183908, 0.142857143, 0, 0,
//' 0.16091954, 0.285714286, 0, 0,
//' 0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
//' 
//' mats_list <- list(XC3, XC4)
//' yr_ord <- c(1, 2)
//' pch_ord <- c(1, 1)
//' 
//' anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA, historical = FALSE,
//'   poporder = 1, patchorder = pch_ord, yearorder = yr_ord)
//'   
//' nullmodel1 <- hist_null(anth_lefkoMat, 1) # Ehrlen format
//' nullmodel2 <- hist_null(anth_lefkoMat, 2) # deVries format
//' 
//' @export hist_null
// [[Rcpp::export(hist_null)]]
Rcpp::List hist_null (RObject mpm, int format = 1, bool err_check = false) {
  
  List mpm_;
  DataFrame mpm_ahstages;
  DataFrame mpm_agestages;
  DataFrame mpm_labels;
  IntegerVector mpm_data_qc;
  DataFrame mpm_model_qc;
  
  if (is<List>(mpm)) {
    List mpm_thru(mpm);
    mpm_ = mpm_thru;
    
    StringVector mpm_class_vec = mpm_.attr("class");
    std::string mpm_class = as<std::string>(mpm_class_vec[0]);
    if (!LefkoUtils::stringcompare_hard(mpm_class, "lefkoMat")) {
      throw Rcpp::exception("Object mpm is of unrecognized class.", false);
    }
    
    DataFrame mpm_ahstages_ = as<DataFrame>(mpm_["ahstages"]);
    DataFrame mpm_hstages_ = as<DataFrame>(mpm_["hstages"]);
    DataFrame mpm_agestages_ = as<DataFrame>(mpm_["agestages"]);
    DataFrame mpm_labels_ = as<DataFrame>(mpm_["labels"]);
    
    if (mpm_ahstages.length() == 1) {
      throw Rcpp::exception("Object mpm is of unrecognized class.", false);
    }
    mpm_ahstages = mpm_ahstages_;
    
    if (mpm_hstages_.length() != 1) {
      throw Rcpp::exception("Input MPM must be ahistorical.", false);
    }

    if (mpm_agestages_.length() != 1) {
      throw Rcpp::exception("Input MPM must be ahistorical, and cannot be age-by-stage.",
        false);
    }
    mpm_agestages = mpm_agestages_;
    
    if (mpm_labels_.length() == 1) {
      throw Rcpp::exception("Object mpm is of unrecognized class.", false);
    }
    mpm_labels = mpm_labels_;
    
  } else {
    throw Rcpp::exception("Object mpm is of unrecognized class.", false);
  }
  
  int tot_dims = mpm_.length();
  int datqc_position = -1;
  int modqc_position = -1;
  
  StringVector mpm_components = mpm_.names();
  for (int i = 0; i < tot_dims; i++) {
    if (stringcompare_hard(as<std::string>(mpm_components(i)), "modelqc")) {
      modqc_position = i;
    }
    if (stringcompare_hard(as<std::string>(mpm_components(i)), "dataqc")) {
      datqc_position = i;
    }
  }
  
  if (datqc_position != -1) {
    IntegerVector mpm_data_qc_ = as<IntegerVector>(mpm_["dataqc"]);
    mpm_data_qc = mpm_data_qc_;
  }
  if (modqc_position != -1) {
    DataFrame mpm_model_qc_ = as<DataFrame>(mpm_["modelqc"]);
    mpm_model_qc = mpm_model_qc_;
  }
  
  List allstages = simplepizzle(mpm_ahstages, format);
  
  DataFrame allindices = as<DataFrame>(allstages["allstages"]);
  DataFrame all_hstages = as<DataFrame>(allstages["hstages"]);
  DataFrame all_ahstages = as<DataFrame>(allstages["ahstages"]);
  
  List redone_mpms = thefifthhousemate(mpm_, allindices, all_hstages,
    all_ahstages, format);
  
  List A_mpms = as<List>(redone_mpms["A"]);
  List U_mpms = as<List>(redone_mpms["U"]);
  List F_mpms = as<List>(redone_mpms["F"]);
  
  int tot_U_elems {0};
  int tot_F_elems {0};
  int tot_mats = A_mpms.length();
  
  for (int i = 0; i < tot_mats; i++) {
    arma::mat current_Umat = as<arma::mat>(U_mpms(i));
    arma::mat current_Fmat = as<arma::mat>(F_mpms(i));
    
    arma::uvec found_U_uvec = find(current_Umat);
    arma::uvec found_F_uvec = find(current_Fmat);
    int found_U = static_cast<int>(found_U_uvec.n_elem);
    int found_F = static_cast<int>(found_F_uvec.n_elem);
    
    tot_U_elems += found_U;
    tot_F_elems += found_F;
  }
  
  IntegerVector new_mqc = {tot_U_elems, tot_F_elems, tot_mats};
  
  if (err_check) tot_dims++;
  
  List final_output(tot_dims);
  CharacterVector listnames (tot_dims);
  
  final_output(0) = A_mpms;
  listnames(0) = "A";
  final_output(1) = U_mpms;
  listnames(1) = "U";
  final_output(2) = F_mpms;
  listnames(2) = "F";
  final_output(3) = all_hstages;
  listnames(3) = "hstages";
  final_output(4) = mpm_agestages;
  listnames(4) = "agestages";
  final_output(5) = all_ahstages;
  listnames(5) = "ahstages";
  final_output(6) = mpm_labels;
  listnames(6) = "labels";
  final_output(7) = new_mqc;
  listnames(7) = "matrixqc";
  
  int next_dude = 8;
  if (modqc_position > -1) {
    final_output(next_dude) = mpm_model_qc;
    listnames(next_dude) = "modelqc";
    next_dude++;
  }
  if (datqc_position > -1) {
    final_output(next_dude) = mpm_data_qc;
    listnames(next_dude) = "dataqc";
  }
  
  if (err_check) {
    final_output(tot_dims - 1) = allstages;
    listnames(tot_dims - 1) = "err_check";
  }
  
  final_output.attr("names") = listnames;
  StringVector needed_classes {"lefkoMat"};
  final_output.attr("class") = needed_classes;
  
  return final_output;
}

//' Estimate Mean Projection Matrices
//' 
//' Function \code{lmean()} estimates mean projection matrices as element-wise
//' arithmetic means. It produces \code{lefkoMat} objects if provided with them,
//' or single matrices in a simple one-element list if provided a list of
//' matrices.
//' 
//' @name lmean
//' 
//' @param mats A \code{lefkoMat} object, or a list of square matrices of equal
//' dimension.
//' @param matsout A string identifying which means to estimate. Option
//' \code{"pop"} indicates population-level only, \code{"patch"} indicates
//' patch-level only, and \code{"all"} indicates that both patch- and
//' population-level means should be estimated. Defaults to \code{"all"}.
//' @param force_sparse A logical value identifying whether to output the mean
//' matrices in sparse format, if input as standard matrices.
//' 
//' @return Yields a \code{lefkoMat} object with the following characteristics:
//' 
//' \item{A}{A list of full mean projection matrices in order of sorted
//' populations, patches, and years. These are typically estimated as the sums
//' of the associated mean \code{U} and \code{F} matrices. All matrices output
//' in either the \code{matrix} class, or the \code{dgCMatrix} class.}
//' \item{U}{A list of mean survival-transition matrices sorted as in \code{A}.
//' All matrices output in the \code{matrix} class.}
//' \item{F}{A list of mean fecundity matrices sorted as in \code{A}. All
//' matrices output in the \code{matrix} class.}
//' \item{hstages}{A data frame showing the pairing of ahistorical stages used
//' to create historical stage pairs. Given if the MPM is historical.}
//' \item{ahstages}{A data frame detailing the characteristics of associated
//' ahistorical stages.}
//' \item{labels}{A data frame detailing the order of population, patch, and
//' year of each mean matrix. If \code{pop}, \code{patch}, or \code{year2} are
//' \code{NA} in the original \code{labels} set, then these will be re-labeled
//' as \code{A}, \code{1}, or \code{1}, respectively.}
//' \item{matrixqc}{A short vector describing the number of non-zero elements in
//' \code{U} and \code{F} mean matrices, and the number of annual matrices.}
//' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite} input.
//' Only output from \code{lefkoMat} objects resulting from function-based
//' estimation.}
//' \item{dataqc}{A vector showing the numbers of individuals and rows in the
//' vertical dataset used as input. Only output from \code{lefkoMat} objects
//' resulting from raw matrix estimation.}
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
//' cyp2mean <- lmean(cypmatrix2r)
//' 
//' @export
// [[Rcpp::export(lmean)]]
Rcpp::List lmean(RObject mats, Nullable<String> matsout = R_NilValue,
  bool force_sparse = false) {
  
  StringVector mats_class;
  StringVector mats_elements;
  String matsout_;
  
  List u_mats;
  List f_mats;
  List a_mats;
  
  DataFrame listofyears;
  DataFrame ahstages;
  DataFrame agestages;
  DataFrame hstages;
  DataFrame modelqc;
  IntegerVector dataqc;
  
  bool historical {false};
  bool poponly {false};
  bool patchonly {false};
  bool dataqc_found {false};
  bool modelqc_found {false};
  bool mat_input {true};
  
  List gd_output;
  
  if (matsout.isNotNull()) {
    String matsout__(matsout);
    matsout_ = matsout__;
    StringVector matsout_choices = {"pop", "patch", "all"};
    int found_one {0};
    for (int i = 0; i < 3; i++) {
      if (stringcompare_simple(matsout_, String(matsout_choices(i)))) {
        found_one++;
      }
    }
    if (found_one == 0 && matsout_ == "") {
      matsout_ = "all";
    } else if (found_one != 1) { 
      throw Rcpp::exception("Argument matsout must equal 'all', 'pop', or 'patch'.",
        false);
    }
  } else {
    matsout_ = "all";
  }
  
  if (stringcompare_simple(matsout_, "all")) { 
    poponly = true;
    patchonly = true;
    
  } else if (stringcompare_simple(matsout_, "pat")) {
    patchonly = true;
    
  } else if (stringcompare_simple(matsout_, "pop")) { 
    poponly = true;
    
  }
  
  if (is<NumericMatrix>(mats)) {
    throw Rcpp::exception("Matrix averaging cannot be performed on a single matrix.",
      false);
    
  } else if (is<List>(mats)) {
    List mats_ = as<List>(mats);
    if (mats_.hasAttribute("class")) {
      mats_class = as<StringVector>(mats_.attr("class"));
    }
    if (mats_.hasAttribute("names")) {
      mats_elements = as<StringVector>(mats_.attr("names"));
    }
    
    int mats_class_length = mats_class.length();
    bool found_it {false};
    for (int i = 0; i < mats_class_length; i++) {
      if (stringcompare_hard(String(mats_class(i)), "lefkoMat")) {
        found_it = true;
      }
    }
    
    if (found_it) {
      // lefkoMat input
      StringVector lefkoMat_required = {"A", "U", "F", "ahstages", "labels"};
      int found_element {0};
      
      for (int i = 0; i < lefkoMat_required.length(); i++) {
        for (int j = 0; j < mats_elements.length(); j++) {
          if (stringcompare_hard(String(lefkoMat_required(i)), String(mats_elements(j)))) {
            found_element++;
          }
        }
      }
      if (found_element < 5) {
        throw Rcpp::exception("Argument mats does not appear to be a lefkoMat object.",
          false);
      }
      
      u_mats = as<List>(mats_["U"]);
      f_mats = as<List>(mats_["F"]);
      ahstages = as<DataFrame>(mats_["ahstages"]);
      agestages = as<DataFrame>(mats_["agestages"]);
      hstages = as<DataFrame>(mats_["hstages"]);
      
      CharacterVector hstages_names = hstages.names();
      
      if (hstages_names.length() > 1) historical = true;
      
      for (int i = 0; i < mats_elements.length(); i++) {
        if (stringcompare_hard(String(mats_elements(i)), "dataqc")) {
          dataqc = as<IntegerVector>(mats_["dataqc"]);
          dataqc_found = true;
        }
        
        if (stringcompare_hard(String(mats_elements(i)), "modelqc")) {
          modelqc = as<DataFrame>(mats_["modelqc"]);
          modelqc_found = true;
        }
      }
      
      // Quality control for listofyears
      DataFrame labels = as<DataFrame>(mats_["labels"]);
      listofyears = loy_inator(labels, true);
      
      if (is<S4>(u_mats(0))) mat_input = false;
      
      // Matrix averaging
      if (historical) {
        gd_output = turbogeodiesel(listofyears, u_mats, f_mats, hstages, 
          agestages, ahstages, patchonly, poponly, mat_input, 0);
        
      } else {
        gd_output = geodiesel(listofyears, u_mats, f_mats, agestages, ahstages,
          patchonly, poponly, mat_input, 0);
      }
    
      gd_output["ahstages"] = ahstages;
      gd_output["hstages"] = hstages;
      gd_output["agestages"] = agestages;
      if (dataqc_found) gd_output["dataqc"] = dataqc;
      if (modelqc_found) gd_output["modelqc"] = modelqc;
      
      StringVector gd_class = {"lefkoMat"};
      gd_output.attr("class") = gd_class;
      
    } else {
      // List of NumericMatrix input
      int mats_length = mats_.length();
      
      if (is<S4>(mats_(0))) mat_input = false;
      
      arma::mat core_mat;
      arma::sp_mat core_mat_sp;
      int mat_rows {0};
      int mat_cols {0};
      
      for (int i = 0; i < mats_length; i++) {
        RObject test_bit = as<RObject>(mats_[i]);
        
        if (!is<NumericMatrix>(test_bit) && !is<S4>(test_bit)) {
          String eat_my_shorts = "Argument mats must be either a lefkoMat object ";
          String eat_my_shorts1 = "or a list of numeric matrices.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
        if (i == 0) {
          if (mat_input) {
            core_mat = as<arma::mat>(mats_[0]);
            
            mat_rows = core_mat.n_rows;
            mat_cols = core_mat.n_cols;
          } else {
            core_mat_sp = as<arma::sp_mat>(mats_[0]);
            
            mat_rows = core_mat_sp.n_rows;
            mat_cols = core_mat_sp.n_cols;
          }
          
          if (mat_rows != mat_cols) {
            throw Rcpp::exception("Input matrices must be square.", false);
          }
          
          if (mat_input) {
            core_mat = core_mat / mats_length;
          } else {
            core_mat_sp = core_mat_sp / mats_length;
          }
        } else {
          if (mat_input) {
            arma::mat next_mat = as<arma::mat>(mats_[i]);
            
            int next_rows = next_mat.n_rows;
            int next_cols = next_mat.n_cols;
            
            if (next_rows != mat_rows || next_cols != mat_cols) {
              throw Rcpp::exception("All input matrices must have the same dimensions.",
                false);
            }
            
            core_mat = core_mat + (next_mat / mats_length);
          } else {
            arma::sp_mat next_mat_sp = as<arma::sp_mat>(mats_[i]);
            
            int next_rows = next_mat_sp.n_rows;
            int next_cols = next_mat_sp.n_cols;
            
            if (next_rows != mat_rows || next_cols != mat_cols) {
              throw Rcpp::exception("All input matrices must have the same dimensions.",
                false);
            }
            
            core_mat_sp = core_mat_sp + (next_mat_sp / mats_length);
          }
        }
        
        if (mat_input && force_sparse) {
          arma::sp_mat core_mat_sp_(core_mat);
          core_mat_sp = core_mat_sp_;
        }
        
        if (mat_input && !force_sparse) {
          List A = List::create(_["A"] = core_mat);
          gd_output = A;
        } else {
          List A = List::create(_["A"] = core_mat_sp);
          gd_output = A;
        }
      }
    }
  }
  
  return gd_output;
}

//' Add a New Stage to an Existing LefkoMat Object
//' 
//' Function \code{add_stage()} adds a new stage to an existing \code{lefkoMat}
//' object. In addition to altering the \code{ahstages} object within the MPM,
//' it alters the \code{hstages} and \code{agestages} objects and adds the
//' appropriate number of new rows and columns depending on the kind of MPM
//' input.
//' 
//' @name add_stage
//' 
//' @param mpm The \code{lefkoMat} object to add a stage to.
//' @param add_before The index of the stage to insert a new stage before. This
//' index should be derived from the \code{ahstages} of the input \code{mpm}.
//' Cannot be set if \code{add_after} is to be used.
//' @param add_after The index of the stage to insert a new stage after. This
//' index should be derived from the \code{ahstages} of the input \code{mpm}.
//' Cannot be set if \code{add_before} is to be used.
//' @param stage_name The name of the new stage to add. Defaults to
//' \code{new_stage}. 
//' 
//' @return A new copy of the original MPM edited to include new rows and
//' columns in the associated matrices, and with \code{ahstages},
//' \code{agestages}, and \code{hstages} objects edited to include the new
//' stage.
//' 
//' @seealso \code{\link{edit_lM}()}
//' 
//' @examples
//' data(cypdata)
//' 
//' cyp_lesl_data <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4, 
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04", 
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04", 
//'   stagesize = "sizeadded", NAas0 = TRUE, age_offset = 2)
//' 
//' cyp_survival <- glm(alive3 ~ obsage + as.factor(year2), data = cyp_lesl_data,
//'   family = "binomial")
//' cyp_fecundity <- glm(feca2 ~ 1 + obsage + as.factor(year2),
//'   data = cyp_lesl_data, family = "poisson")
//' 
//' mod_params <- create_pm(name_terms = TRUE)
//' mod_params$modelparams[22] <- "obsage"
//' 
//' germination <- 0.08
//' protocorm_to_seedling <- 0.10
//' seeding_to_adult <- 0.20
//' seeds_per_fruit <- 8000
//' 
//' cyp_lesl_supp <- supplemental(historical = FALSE, stagebased = FALSE,
//'   agebased = TRUE, age2 = c(1, 2), type = c(1, 1),
//'   givenrate = c(protocorm_to_seedling, seeding_to_adult))
//' 
//' cyp_lesl_fb_mpm <- fleslie(data = cyp_lesl_data, surv_model = cyp_survival,
//'   fec_model = cyp_fecundity, paramnames = mod_params, last_age = 7,
//'   fecage_min = 3, fecmod = (germination * seeds_per_fruit),
//'   supplement = cyp_lesl_supp)
//' 
//' altered1 <- add_stage(cyp_lesl_fb_mpm, add_before = 1, stage_name = "DS")
//' 
//' @export add_stage
// [[Rcpp::export(add_stage)]]
Rcpp::List add_stage (const RObject mpm, int add_before = 0,
  int add_after = 0, Nullable<CharacterVector> stage_name  = R_NilValue) {
  
  if (add_before > 0 && add_after > 0) {
    throw Rcpp::exception("Please set either add_before or add_after, but not both.",
      false);
  }
  
  String stage_name_to_add;
  if (stage_name.isNotNull()) {
    CharacterVector stage_name_input = as<CharacterVector>(stage_name);
    
    if (stage_name_input.length() > 1) {
      throw Rcpp::exception("Please enter only one new stage name.", false);
    } else if (stage_name_input.length() == 0) {
      stage_name_to_add = "new_stage";
    } else {
      stage_name_to_add = stage_name_input(0);
    }
  } else {
    stage_name_to_add = "new_stage";
  }
  
  List mpm_list;
  
  if (is<List>(mpm)) mpm_list = mpm;
  StringVector mpm_class = mpm_list.attr("class");
  
  String mpm_error = "Please enter a lefkoMat object as input.";
  bool mpm_yes {false};
  
  if (!mpm_list.containsElementNamed("ahstages")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  if (!mpm_list.containsElementNamed("hstages")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  if (!mpm_list.containsElementNamed("agestages")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  if (!mpm_list.containsElementNamed("labels")) {
    throw Rcpp::exception(mpm_error.get_cstring(), false);
  }
  for (int i = 0; i < static_cast<int>(mpm_class.length()); i++) {
    if (mpm_class(i) == "lefkoMat") mpm_yes = true;
  }
  if (!mpm_yes) throw Rcpp::exception(mpm_error.get_cstring(), false);
  
  Rcpp::DataFrame ahstages = as<DataFrame>(mpm_list["ahstages"]);
  Rcpp::DataFrame hstages = as<DataFrame>(mpm_list["hstages"]);
  Rcpp::DataFrame agestages = as<DataFrame>(mpm_list["agestages"]);
  Rcpp::DataFrame labels = as<DataFrame>(mpm_list["labels"]);
  
  int wtf = LefkoUtils::whichbrew(ahstages, hstages, agestages);
  // wtf possible results: \code{0}: historical MPM, \code{1}:
  // ahistorical MPM, \code{2}: age-by-stage MPM, and \code{3}: age-based MPM
  
  StringVector stagevec = as<StringVector>(ahstages["stage"]);
  int num_stages = stagevec.length();
  
  for (int i = 0; i < num_stages; i++) {
    if (LefkoUtils::stringcompare_hard(stage_name_to_add, as<std::string>(stagevec(i)))) {
      throw Rcpp::exception("Entered stage_name cannot be the same as an existing stage.",
        false);
    }
  }
  
  if (add_before > num_stages) {
    throw Rcpp::exception("Fewer stages exist than suggested by number input in option add_before.",
      false);
  } else if (add_after > (num_stages + 1)) {
    throw Rcpp::exception("Fewer stages exist than suggested by number input in option add_after.",
      false);
  }
  
  // New stageframe
  IntegerVector stageidvec = as<IntegerVector>(ahstages["stage_id"]);
  NumericVector origsizevec = as<NumericVector>(ahstages["original_size"]);
  NumericVector origsizebvec = as<NumericVector>(ahstages["original_size_b"]);
  NumericVector origsizecvec = as<NumericVector>(ahstages["original_size_c"]);
  NumericVector minagevec = as<NumericVector>(ahstages["min_age"]);
  NumericVector maxagevec = as<NumericVector>(ahstages["max_age"]);
  IntegerVector repvec = as<IntegerVector>(ahstages["repstatus"]);
  IntegerVector obsvec = as<IntegerVector>(ahstages["obsstatus"]);
  IntegerVector propvec = as<IntegerVector>(ahstages["propstatus"]);
  IntegerVector immvec = as<IntegerVector>(ahstages["immstatus"]);
  IntegerVector matvec = as<IntegerVector>(ahstages["matstatus"]);
  IntegerVector repentryvec = as<IntegerVector>(ahstages["entrystage"]);
  IntegerVector indvec = as<IntegerVector>(ahstages["indataset"]);
  NumericVector binvec = as<NumericVector>(ahstages["binhalfwidth_raw"]);
  NumericVector binbvec = as<NumericVector>(ahstages["binhalfwidthb_raw"]);
  NumericVector bincvec = as<NumericVector>(ahstages["binhalfwidthc_raw"]);
  NumericVector sizeminvec = as<NumericVector>(ahstages["sizebin_min"]);
  NumericVector sizemaxvec = as<NumericVector>(ahstages["sizebin_max"]);
  NumericVector sizectrvec = as<NumericVector>(ahstages["sizebin_center"]);
  NumericVector sizewidthvec = as<NumericVector>(ahstages["sizebin_width"]);
  NumericVector sizebminvec = as<NumericVector>(ahstages["sizebinb_min"]);
  NumericVector sizebmaxvec = as<NumericVector>(ahstages["sizebinb_max"]);
  NumericVector sizebctrvec = as<NumericVector>(ahstages["sizebinb_center"]);
  NumericVector sizebwidthvec = as<NumericVector>(ahstages["sizebinb_width"]);
  NumericVector sizecminvec = as<NumericVector>(ahstages["sizebinc_min"]);
  NumericVector sizecmaxvec = as<NumericVector>(ahstages["sizebinc_max"]);
  NumericVector sizecctrvec = as<NumericVector>(ahstages["sizebinc_center"]);
  NumericVector sizecwidthvec = as<NumericVector>(ahstages["sizebinc_width"]);
  IntegerVector groupvec = as<IntegerVector>(ahstages["group"]);
  StringVector commentsvec = as<StringVector>(ahstages["comments"]);
  IntegerVector alivevec = as<IntegerVector>(ahstages["alive"]);
  IntegerVector almostbornvec = as<IntegerVector>(ahstages["almostborn"]);
  
  IntegerVector new_stageidvec (num_stages + 1);
  int ns_counter {0};
  int ns_place {0};
  for (int i = 0; i < (num_stages + 1); i++) {
    if (i == add_before - 1 || i == add_after) {
      new_stageidvec(i) = num_stages + 1;
      ns_place = i;
    } else {
      new_stageidvec(i) = stageidvec(ns_counter);
      ns_counter++;
    }
  }
  
  StringVector new_stagevec(num_stages + 1);
  NumericVector new_origsizevec(num_stages + 1);
  NumericVector new_origsizebvec(num_stages + 1);
  NumericVector new_origsizecvec(num_stages + 1);
  NumericVector new_minagevec(num_stages + 1);
  NumericVector new_maxagevec(num_stages + 1);
  IntegerVector new_repvec(num_stages + 1);
  IntegerVector new_obsvec(num_stages + 1);
  IntegerVector new_propvec(num_stages + 1);
  IntegerVector new_immvec(num_stages + 1);
  IntegerVector new_matvec(num_stages + 1);
  IntegerVector new_repentryvec(num_stages + 1);
  IntegerVector new_indvec(num_stages + 1);
  NumericVector new_binvec(num_stages + 1);
  NumericVector new_binbvec(num_stages + 1);
  NumericVector new_bincvec(num_stages + 1);
  NumericVector new_sizeminvec(num_stages + 1);
  NumericVector new_sizemaxvec(num_stages + 1);
  NumericVector new_sizectrvec(num_stages + 1);
  NumericVector new_sizewidthvec(num_stages + 1);
  NumericVector new_sizebminvec(num_stages + 1);
  NumericVector new_sizebmaxvec(num_stages + 1);
  NumericVector new_sizebctrvec(num_stages + 1);
  NumericVector new_sizebwidthvec(num_stages + 1);
  NumericVector new_sizecminvec(num_stages + 1);
  NumericVector new_sizecmaxvec(num_stages + 1);
  NumericVector new_sizecctrvec(num_stages + 1);
  NumericVector new_sizecwidthvec(num_stages + 1);
  IntegerVector new_groupvec(num_stages + 1);
  StringVector new_commentsvec(num_stages + 1);
  IntegerVector new_alivevec(num_stages + 1);
  IntegerVector new_almostbornvec(num_stages + 1);
  
  int place_correction {0};
  for (int i = 0; i < (num_stages + 1); i++) {
    if ((i == (add_before - 1) && add_before > 0) || (i == add_after && add_after > 0)) {
      new_stagevec(i) = stage_name_to_add;
      new_origsizevec(i) = 0.0;
      new_origsizebvec(i) = 0.0;
      new_origsizecvec(i) = 0.0;
      new_minagevec(i) = 0;
      new_maxagevec(i) = 0;
      new_obsvec(i) = 1;
      new_propvec(i) = 0;
      
      if (add_before == 1) {
        new_immvec(i) = 1;
        new_matvec(i) = 0;
        new_repvec(i) = 0;
        new_repentryvec(i) = 1;
      } else {
        new_immvec(i) = 0;
        new_matvec(i) = 1;
        new_repvec(i) = 1;
        new_repentryvec(i) = 0;
      }
      
      new_indvec(i) = 0;
      new_binvec(i) = 1.0;
      new_binbvec(i) = 1.0;
      new_bincvec(i) = 1.0;
      new_sizeminvec(i) = 0.5;
      new_sizemaxvec(i) = 1.5;
      new_sizectrvec(i) = 1.0;
      new_sizewidthvec(i) = 0.5;
      new_sizebminvec(i) = 0.5;
      new_sizebmaxvec(i) = 1.5;
      new_sizebctrvec(i) = 1.0;
      new_sizebwidthvec(i) = 0.5;
      new_sizecminvec(i) = 0.5;
      new_sizecmaxvec(i) = 1.5;
      new_sizecctrvec(i) = 1.0;
      new_sizecwidthvec(i) = 0.5;
      new_groupvec(i) = 0;
      new_commentsvec(i) = "new stage";
      new_alivevec(i) = 1.0;
      new_almostbornvec(i) = 0.0;
      
      place_correction++;
    } else {
      new_stagevec(i) = stagevec(i - place_correction);
      new_origsizevec(i) = origsizevec(i - place_correction);
      new_origsizebvec(i) = origsizebvec(i - place_correction);
      new_origsizecvec(i) = origsizecvec(i - place_correction);
      new_minagevec(i) = minagevec(i - place_correction);
      new_maxagevec(i) = maxagevec(i - place_correction);
      new_repvec(i) = repvec(i - place_correction);
      new_obsvec(i) = obsvec(i - place_correction);
      new_propvec(i) = propvec(i - place_correction);
      new_immvec(i) = immvec(i - place_correction);
      new_matvec(i) = matvec(i - place_correction);
      new_repentryvec(i) = repentryvec(i - place_correction);
      new_indvec(i) = indvec(i - place_correction);
      new_binvec(i) = binvec(i - place_correction);
      new_binbvec(i) = binbvec(i - place_correction);
      new_bincvec(i) = bincvec(i - place_correction);
      new_sizeminvec(i) = sizeminvec(i - place_correction);
      new_sizemaxvec(i) = sizemaxvec(i - place_correction);
      new_sizectrvec(i) = sizectrvec(i - place_correction);
      new_sizewidthvec(i) = sizewidthvec(i - place_correction);
      new_sizebminvec(i) = sizebminvec(i - place_correction);
      new_sizebmaxvec(i) = sizebmaxvec(i - place_correction);
      new_sizebctrvec(i) = sizebctrvec(i - place_correction);
      new_sizebwidthvec(i) = sizebwidthvec(i - place_correction);
      new_sizecminvec(i) = sizecminvec(i - place_correction);
      new_sizecmaxvec(i) = sizecmaxvec(i - place_correction);
      new_sizecctrvec(i) = sizecctrvec(i - place_correction);
      new_sizecwidthvec(i) = sizecwidthvec(i - place_correction);
      new_groupvec(i) = groupvec(i - place_correction);
      new_commentsvec(i) = commentsvec(i - place_correction);
      new_alivevec(i) = alivevec(i - place_correction);
      new_almostbornvec(i) = almostbornvec(i - place_correction);
    }
  }
  
  Rcpp::List new_stageframe(33);
  
  new_stageframe(0) = new_stageidvec;
  new_stageframe(1) = new_stagevec;
  new_stageframe(2) = new_origsizevec;
  new_stageframe(3) = new_origsizebvec;
  new_stageframe(4) = new_origsizecvec;
  new_stageframe(5) = new_minagevec;
  new_stageframe(6) = new_maxagevec;
  new_stageframe(7) = new_repvec;
  new_stageframe(8) = new_obsvec;
  new_stageframe(9) = new_propvec;
  new_stageframe(10) = new_immvec;
  new_stageframe(11) = new_matvec;
  new_stageframe(12) = new_repentryvec;
  new_stageframe(13) = new_indvec;
  new_stageframe(14) = new_binvec;
  new_stageframe(15) = new_sizeminvec;
  new_stageframe(16) = new_sizemaxvec;
  new_stageframe(17) = new_sizectrvec;
  new_stageframe(18) = new_sizewidthvec;
  new_stageframe(19) = new_binbvec;
  new_stageframe(20) = new_sizebminvec;
  new_stageframe(21) = new_sizebmaxvec;
  new_stageframe(22) = new_sizebctrvec;
  new_stageframe(23) = new_sizebwidthvec;
  new_stageframe(24) = new_bincvec;
  new_stageframe(25) = new_sizecminvec;
  new_stageframe(26) = new_sizecmaxvec;
  new_stageframe(27) = new_sizecctrvec;
  new_stageframe(28) = new_sizecwidthvec;
  new_stageframe(29) = new_groupvec;
  new_stageframe(30) = new_commentsvec;
  new_stageframe(31) = new_alivevec;
  new_stageframe(32) = new_almostbornvec;
  
  CharacterVector namevec = {"stage_id", "stage", "original_size",
    "original_size_b", "original_size_c", "min_age", "max_age", "repstatus",
    "obsstatus", "propstatus", "immstatus", "matstatus", "entrystage",
    "indataset", "binhalfwidth_raw", "sizebin_min", "sizebin_max",
    "sizebin_center", "sizebin_width", "binhalfwidthb_raw", "sizebinb_min",
    "sizebinb_max", "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw",
    "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
    "group", "comments", "alive", "almostborn"};
  CharacterVector new_classes = {"data.frame", "stageframe"};
  new_stageframe.attr("names") = namevec;
  new_stageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      (num_stages + 1));
  new_stageframe.attr("class") = new_classes;
  
  DataFrame new_hstages;
  DataFrame new_agestages;
  
  // Vector of insertion rows/columns
  arma::uvec rows_cols_to_add_vec;
  unsigned int num_rows_to_add {0};
  unsigned int num_stages_base {0};
  
  if (wtf == 1 || wtf == 3) {
    unsigned int new_col_row {0};
    if (add_before > 0) {
      new_col_row = add_before - 1;
    } else if (add_after > 0) {
      new_col_row = add_after;
    }
    
    arma::uvec new_cols_rows(1, fill::zeros);
    new_cols_rows(0) = new_col_row;
    rows_cols_to_add_vec = new_cols_rows;
    num_rows_to_add = 1;
    num_stages_base = num_stages;
    
    DataFrame old_hstages = as<DataFrame>(mpm_list["hstages"]);
    DataFrame old_agestages = as<DataFrame>(mpm_list["agestages"]);
    new_hstages = old_hstages;
    new_agestages = old_agestages;
    
  } else if (wtf == 0) {
    Rcpp::DataFrame hstages = as<DataFrame>(mpm_list["hstages"]);
    arma::ivec sid2 = as<arma::ivec>(hstages["stage_id_2"]);
    arma::ivec sid1 = as<arma::ivec>(hstages["stage_id_1"]);
    StringVector h_stage2 = as<StringVector>(hstages["stage_2"]);
    StringVector h_stage1 = as<StringVector>(hstages["stage_1"]);
    
    arma::uvec first_stretch; // For stage at time t
    arma::uvec second_stretch(num_stages + 1, fill::zeros); // For stage at time t-1
    
    unsigned int first_new_col_row {0};
    if (add_before > 0) {
      arma::uvec found_at_t2 = find(sid2 == add_before);
      first_stretch = found_at_t2;
      
      first_new_col_row = add_before - 1;
      
    } else if (add_after > 0) {
      arma::uvec found_at_t2 = find(sid2 == add_after);
      first_stretch = found_at_t2 + 1;
      
      first_new_col_row = add_after;
    }
    
    for (int i = 0; i < (num_stages + 1); i++) {
      second_stretch(i) = first_new_col_row * (num_stages + 1) + 1 + i;
    }
    
    
    int num_first_stretch = first_stretch.n_elem;
    int num_second_stretch = second_stretch.n_elem;
    num_rows_to_add = num_first_stretch + num_second_stretch;
    
    int hstages_dim = sid2.n_elem;
    num_stages_base = hstages_dim;
    
    DataFrame old_agestages = as<DataFrame>(mpm_list["agestages"]);
    new_agestages = old_agestages;
    
    int new_vec_length = num_rows_to_add + num_stages_base;
    
    IntegerVector new_sid2 (new_vec_length);
    IntegerVector new_sid1 (new_vec_length);
    StringVector new_h_stage2 (new_vec_length);
    StringVector new_h_stage1 (new_vec_length);
    
    int vec_length_counter = num_stages_base - 1;
    int addition_adj = num_rows_to_add;
    bool ratcheted {false};
    
    arma::uvec rctav1 (num_first_stretch);
    int rctav1_adj {1};
    for (int i = new_vec_length - 1; i >= 0; i--) {
      int check_index = i - addition_adj + 1;
      arma::uvec found_elems_first = find(first_stretch == check_index);
      arma::uvec found_elems_second = find(second_stretch == check_index);
      
      arma::uvec found_greater_second = find(second_stretch > i);
      int i_corr = 0;
      if (found_greater_second.n_elem > 0) i_corr = second_stretch.n_elem;
      
      if ((i - i_corr) > -1) {
        if (vec_length_counter >= 0 || (add_before == 1 && vec_length_counter >= -1)) {
          if (found_elems_first.n_elem > 0 && !ratcheted) {
            new_sid2(i - i_corr) = num_stages + 1;
            new_h_stage2(i - i_corr) = stage_name_to_add;
            
            if (!(add_before == 1)) {
              new_sid1(i - i_corr) = sid1(vec_length_counter);
              new_h_stage1(i - i_corr) = h_stage1(vec_length_counter);
            } else {
              new_sid1(i - i_corr) = sid1(vec_length_counter + 1);
              new_h_stage1(i - i_corr) = h_stage1(vec_length_counter + 1);
            }
            ratcheted = true;
            addition_adj--;
            
            rctav1(num_first_stretch - rctav1_adj) = i - i_corr;
            rctav1_adj++;
            
          } else if (vec_length_counter >= 0) {
            new_sid2(i - i_corr) = sid2(vec_length_counter);
            new_sid1(i - i_corr) = sid1(vec_length_counter);
            new_h_stage2(i - i_corr) = h_stage2(vec_length_counter);
            new_h_stage1(i - i_corr) = h_stage1(vec_length_counter);
            
            vec_length_counter--;
            ratcheted = false;
          }
        }
      }
    }
    
    arma::uvec rctav2 (num_second_stretch);
    for (int i = 0; i < num_second_stretch; i++) {
      new_sid2(second_stretch(i) - 1) = new_stageidvec(i);
      new_h_stage2(second_stretch(i) - 1) = new_stagevec(i);
      new_sid1(second_stretch(i) - 1) = num_stages + 1;
      new_h_stage1(second_stretch(i) - 1) = new_stagevec(ns_place);
      rctav2(i) = second_stretch(i) - 1;
    }
    
    // Create first_stretch and second_stretch, find unique elements, & sort
    arma::uvec all_stretch = arma::join_cols(rctav1, rctav2);
    arma::uvec all_sorted = arma::sort(all_stretch);
    rows_cols_to_add_vec = all_sorted;
    
    List new_hstages_pre (4);
    new_hstages_pre(0) = new_sid2;
    new_hstages_pre(1) = new_sid1;
    new_hstages_pre(2) = new_h_stage2;
    new_hstages_pre(3) = new_h_stage1;
    
    CharacterVector hstage_namevec = {"stage_id_2", "stage_id_1", "stage_2",
      "stage_1"};
    CharacterVector hstage_classes = {"data.frame"};
    new_hstages_pre.attr("names") = hstage_namevec;
    new_hstages_pre.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      new_vec_length);
    new_hstages_pre.attr("class") = hstage_classes;
    new_hstages = new_hstages_pre;
    
  } else if (wtf == 2) {
    Rcpp::DataFrame agestages = as<DataFrame>(mpm_list["agestages"]);
    arma::ivec sid = as<arma::ivec>(agestages["stage_id"]);
    StringVector age_stage = as<StringVector>(agestages["stage"]);
    IntegerVector age_age = as<IntegerVector>(agestages["age"]);
    
    arma::uvec first_stretch; // For stage at time t
    
    if (add_before > 0) {
      arma::uvec found_at_t2 = find(sid == add_before);
      first_stretch = found_at_t2; // Might need to add -1
      
    } else if (add_after > 0) {
      arma::uvec found_at_t2 = find(sid == add_after);
      first_stretch = found_at_t2 + 1; // Might need to remove +1
    }
    
    rows_cols_to_add_vec = first_stretch;
    num_rows_to_add = first_stretch.n_elem;
    
    int agestages_dim = sid.n_elem;
    num_stages_base = agestages_dim;
    
    DataFrame old_hstages = as<DataFrame>(mpm_list["hstages"]);
    new_hstages = old_hstages;
    
    int new_vec_length = num_rows_to_add + num_stages_base;
    
    IntegerVector new_sid (new_vec_length);
    StringVector new_age_stage (new_vec_length);
    IntegerVector new_age_age (new_vec_length);
    
    int vec_length_counter = num_stages_base - 1;
    int addition_adj = num_rows_to_add;
    bool ratcheted {false};
    
    for (int i = new_vec_length - 1; i >= 0; i--) {
      int check_index = i - addition_adj + 1;
      arma::uvec found_elems = find(rows_cols_to_add_vec == check_index);
      
      if (found_elems.n_elem > 0 && !ratcheted) {
        new_sid(i) = num_stages + 1;
        new_age_stage(i) = stage_name_to_add;
        
        if (!(add_before == 1)) {
          new_age_age(i) = age_age(vec_length_counter - 1);
        } else {
          new_age_age(i) = 0;
        }
        ratcheted = true;
        addition_adj--;
        
      } else {
        new_sid(i) = sid(vec_length_counter);
        new_age_stage(i) = age_stage(vec_length_counter);
        new_age_age(i) = age_age(vec_length_counter);
        
        vec_length_counter--;
        ratcheted = false;
      }
    }
    
    List new_agestages_pre (3);
    new_agestages_pre(0) = new_sid;
    new_agestages_pre(1) = new_age_stage;
    new_agestages_pre(2) = new_age_age;
    
    CharacterVector agestage_namevec = {"stage_id", "stage", "age"};
    CharacterVector agestage_classes = {"data.frame"};
    new_agestages_pre.attr("names") = agestage_namevec;
    new_agestages_pre.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      new_vec_length);
    new_agestages_pre.attr("class") = agestage_classes;
    new_agestages = new_agestages_pre;
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
  
  if (mpm_list.containsElementNamed("A")) {
    if (is<List>(mpm_list["A"])) {
      List new_A_mats = as<List>(mpm_list["A"]);
      A_mats = clone(new_A_mats);
      A_used = true;
      mat_num = static_cast<int>(A_mats.length());
      
      if (is<S4>(A_mats(0))) mat_sparse = true;
    }
    if (is<List>(mpm_list["U"])) {
      List new_U_mats = as<List>(mpm_list["U"]);
      U_mats = clone(new_U_mats);
      U_used = true;
      mat_num = static_cast<int>(U_mats.length());
      
      if (is<S4>(U_mats(0))) mat_sparse = true;
    }
    if (is<List>(mpm_list["F"])) {
      List new_F_mats = as<List>(mpm_list["F"]);
      F_mats = clone(new_F_mats);
      F_used = true;
      mat_num = static_cast<int>(F_mats.length());
      
      if (is<S4>(F_mats(0))) mat_sparse = true;
    }
  }
  
  if (!A_used && !U_used && !F_used) {
    throw Rcpp::exception("The input object does not appear to hold matrices.", false);
  }
  
  if (U_used && F_used) {
    for (int i = 0; i < mat_num; i++) {
      if (!mat_sparse) {
        arma::mat current_U = as<arma::mat>(U_mats(i));
        arma::mat current_F = as<arma::mat>(F_mats(i));
        
        unsigned int old_size = current_U.n_cols;
        unsigned int total_size = num_stages_base + num_rows_to_add;
        unsigned int found_new_row {0};
        
        arma::mat new_U (total_size, total_size);
        arma::mat new_F (total_size, total_size);
        
        for (unsigned int j = 0; j < total_size; j++) {
          arma::uvec current_row_found = find(rows_cols_to_add_vec == j);
          
          if (current_row_found.n_elem == 0) {
            unsigned int found_new_col {0};
            for (unsigned int k = 0; k < total_size; k++) {
              arma::uvec current_col_found = find(rows_cols_to_add_vec == k);
              
              if (current_col_found.n_elem == 0) {
                if ((j - found_new_row) < old_size && (k - found_new_col) < old_size) {
                  new_U (j, k) = current_U(j - found_new_row, k - found_new_col);
                  new_F (j, k) = current_F(j - found_new_row, k - found_new_col);
                }
              } else {
                found_new_col++;
              }
            }
          } else {
            found_new_row++;
          }
        }
        
        arma::mat new_A = new_U + new_F;
        A_mats(i) = new_A;
        U_mats(i) = new_U;
        F_mats(i) = new_F;
        
      } else {
        arma::sp_mat current_U = as<arma::sp_mat>(U_mats(i));
        arma::sp_mat current_F = as<arma::sp_mat>(F_mats(i));
        
        unsigned int old_size = current_U.n_cols;
        unsigned int total_size = num_stages_base + num_rows_to_add;
        unsigned int found_new_row {0};
        
        arma::sp_mat new_U (total_size, total_size);
        arma::sp_mat new_F (total_size, total_size);
        
        for (unsigned int j = 0; j < total_size; j++) {
          arma::uvec current_row_found = find(rows_cols_to_add_vec == j);
          
          if (current_row_found.n_elem == 0) {
            unsigned int found_new_col {0};
            for (unsigned int k = 0; k < total_size; k++) {
              arma::uvec current_col_found = find(rows_cols_to_add_vec == k);
              
              if (current_col_found.n_elem == 0) {
                if ((j - found_new_row) < old_size && (k - found_new_col) < old_size) {
                  new_U (j, k) = current_U(j - found_new_row, k - found_new_col);
                  new_F (j, k) = current_F(j - found_new_row, k - found_new_col);
                }
              } else {
                found_new_col++;
              }
            }
          } else {
            found_new_row++;
          }
        }
        
        arma::sp_mat new_A = new_U + new_F;
        A_mats(i) = new_A;
        U_mats(i) = new_U;
        F_mats(i) = new_F;
      }
    }
  } else {
    for (int i = 0; i < mat_num; i++) {
      if (!mat_sparse) {
        arma::mat current_A = as<arma::mat>(A_mats(i));
        
        unsigned int old_size = current_A.n_cols;
        unsigned int total_size = num_stages_base + num_rows_to_add;
        unsigned int found_new_row {0};
        
        arma::mat new_A (total_size, total_size);
        
        for (unsigned int j = 0; j < total_size; j++) {
          arma::uvec current_row_found = find(rows_cols_to_add_vec == j);
          
          if (current_row_found.n_elem == 0) {
            unsigned int found_new_col {0};
            for (unsigned int k = 0; k < total_size; k++) {
              arma::uvec current_col_found = find(rows_cols_to_add_vec == k);
              
              if (current_col_found.n_elem == 0) {
                if ((j - found_new_row) < old_size && (k - found_new_col) < old_size) {
                  new_A (j, k) = current_A(j - found_new_row, k - found_new_col);
                }
              } else {
                found_new_col++;
              }
            }
          } else {
            found_new_row++;
          }
        }
        
        A_mats(i) = new_A;
        
      } else {
        arma::sp_mat current_A = as<arma::sp_mat>(A_mats(i));
        
        unsigned int old_size = current_A.n_cols;
        unsigned int total_size = num_stages_base + num_rows_to_add;
        unsigned int found_new_row {0};
        
        arma::sp_mat new_A (total_size, total_size);
        
        for (unsigned int j = 0; j < total_size; j++) {
          arma::uvec current_row_found = find(rows_cols_to_add_vec == j);
          
          if (current_row_found.n_elem == 0) {
            unsigned int found_new_col {0};
            for (unsigned int k = 0; k < total_size; k++) {
              arma::uvec current_col_found = find(rows_cols_to_add_vec == k);
              
              if (current_col_found.n_elem == 0) {
                if ((j - found_new_row) < old_size && (k - found_new_col) < old_size) {
                  new_A (j, k) = current_A(j - found_new_row, k - found_new_col);
                }
              } else {
                found_new_col++;
              }
            }
          } else {
            found_new_row++;
          }
        }
        
        A_mats(i) = new_A;
      }
    }
  }
  
  // MPM final processing
  bool found_modelqc = mpm_list.containsElementNamed("modelqc");
  bool found_dataqc = mpm_list.containsElementNamed("dataqc");
  
  List true_output;
  if (!found_modelqc && !found_dataqc) {
    DataFrame new_labels = as<DataFrame>(mpm_list["labels"]);
    IntegerVector new_matrixqc = as<IntegerVector>(mpm_list["matrixqc"]);
    
    List new_mpm (8);
    new_mpm(0) = A_mats;
    new_mpm(1) = U_mats;
    new_mpm(2) = F_mats;
    new_mpm(3) = new_stageframe;
    new_mpm(4) = new_agestages;
    new_mpm(5) = new_hstages;
    new_mpm(6) = new_labels;
    new_mpm(7) = new_matrixqc;
    
    CharacterVector new_mpm_namevec = {"A", "U", "F", "ahstages", "agestages",
      "hstages", "labels", "matrixqc"};
    CharacterVector new_mpm_classes = {"lefkoMat"};
    new_mpm.attr("names") = new_mpm_namevec;
    new_mpm.attr("class") = new_mpm_classes;
    
    true_output = new_mpm;
    
  } else if (!found_modelqc) {
    DataFrame new_labels = as<DataFrame>(mpm_list["labels"]);
    IntegerVector new_matrixqc = as<IntegerVector>(mpm_list["matrixqc"]);
    IntegerVector new_dataqc = as<IntegerVector>(mpm_list["dataqc"]);
    
    List new_mpm (9);
    new_mpm(0) = A_mats;
    new_mpm(1) = U_mats;
    new_mpm(2) = F_mats;
    new_mpm(3) = new_stageframe;
    new_mpm(4) = new_agestages;
    new_mpm(5) = new_hstages;
    new_mpm(6) = new_labels;
    new_mpm(7) = new_matrixqc;
    new_mpm(8) = new_dataqc;
    
    CharacterVector new_mpm_namevec = {"A", "U", "F", "ahstages", "agestages",
      "hstages", "labels", "matrixqc", "dataqc"};
    CharacterVector new_mpm_classes = {"lefkoMat"};
    new_mpm.attr("names") = new_mpm_namevec;
    new_mpm.attr("class") = new_mpm_classes;
    
    true_output = new_mpm;
    
  } else if (!found_dataqc) {
    DataFrame new_labels = as<DataFrame>(mpm_list["labels"]);
    IntegerVector new_matrixqc = as<IntegerVector>(mpm_list["matrixqc"]);
    DataFrame new_modelqc = as<DataFrame>(mpm_list["modelqc"]);
    
    List new_mpm (9);
    new_mpm(0) = A_mats;
    new_mpm(1) = U_mats;
    new_mpm(2) = F_mats;
    new_mpm(3) = new_stageframe;
    new_mpm(4) = new_agestages;
    new_mpm(5) = new_hstages;
    new_mpm(6) = new_labels;
    new_mpm(7) = new_matrixqc;
    new_mpm(8) = new_modelqc;
    
    CharacterVector new_mpm_namevec = {"A", "U", "F", "ahstages", "agestages",
      "hstages", "labels", "matrixqc", "modelqc"};
    CharacterVector new_mpm_classes = {"lefkoMat"};
    new_mpm.attr("names") = new_mpm_namevec;
    new_mpm.attr("class") = new_mpm_classes;
    
    true_output = new_mpm;
    
  } else {
    DataFrame new_labels = as<DataFrame>(mpm_list["labels"]);
    IntegerVector new_matrixqc = as<IntegerVector>(mpm_list["matrixqc"]);
    DataFrame new_modelqc = as<DataFrame>(mpm_list["modelqc"]);
    IntegerVector new_dataqc = as<IntegerVector>(mpm_list["dataqc"]);
    
    List new_mpm (10);
    new_mpm(0) = A_mats;
    new_mpm(1) = U_mats;
    new_mpm(2) = F_mats;
    new_mpm(3) = new_stageframe;
    new_mpm(4) = new_agestages;
    new_mpm(5) = new_hstages;
    new_mpm(6) = new_labels;
    new_mpm(7) = new_matrixqc;
    new_mpm(8) = new_modelqc;
    new_mpm(9) = new_dataqc;
    
    CharacterVector new_mpm_namevec = {"A", "U", "F", "ahstages", "agestages",
      "hstages", "labels", "matrixqc", "modelqc", "dataqc"};
    CharacterVector new_mpm_classes = {"lefkoMat"};
    new_mpm.attr("names") = new_mpm_namevec;
    new_mpm.attr("class") = new_mpm_classes;
    
    true_output = new_mpm;
  }
  
  return true_output;;
}

//' Check Continuity of Life Cycle through Matrices in lefkoMat Objects
//' 
//' Function \code{cycle_check()} tests whether stages, stage-pairs, or
//' age-stages connect in matrices within \code{lefkoMat} objects.
//' 
//' @name cycle_check
//' 
//' @param mpm An object of class lefkoMat, a matrix, or a list of matrices.
//' @param quiet A logical variable indicating whether to suppress diagnostic
//' messages. Defaults to \code{FALSE}.
//' 
//' @return Returns a list with two elements, both of which are also lists.
//' The first list, \code{no_in}, contains as many elements as matrices, with
//' each element containing an integer vector showing the identification numbers
//' of stages, stage-pairs, or age-stages, in each matrix that do not show any
//' transitions leading to them. The second list, \code{no_out}, is structured
//' similarly to the first, but shows stages, stage-pairs, or age-stages from
//' which there are no transitions leading out.
//' 
//' @section Notes:
//' This function tests whether stages, stage-pairs, and age-stages are
//' connected to others in matrices used for projection. Whether stages,
//' stage-pairs, or age-stages are shown depends on whether the MPM is
//' ahistorical / age-based, historical stage-based, or age-by-stage,
//' respectively. Checks are performed by testing whether each column in a
//' matrix includes non-zero transitions to other columns, and by testing
//' whether any columns have no transitions to them from other columns. If any
//' such columns are found, then function \code{cycle_check}
//' will export an integer vector giving the column numbers with problems.
//' These column numbers may then be checked against the \code{stage_id} column
//' of the associated stageframe in the case of a ahistorical or age-based MPM,
//' against the row number of the associated \code{hstages} data frame in the
//' case of a historical MPM, or against the row number of the associated
//' \code{agestages} data frame in the case of an age-by-stage MPM.
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
//' cycle_check(cypmatrix2r)
//' 
//' @export cycle_check
// [[Rcpp::export(cycle_check)]]
List cycle_check(RObject mpm, Nullable<RObject> quiet = R_NilValue) {
  
  List lopv (2);
  bool quiet_bool = false;
  
  if (quiet.isNotNull()) {
    if (is<RObject>(quiet)) {
      quiet_bool = LefkoInputs::yesno_to_logic(as<RObject>(quiet), "quiet");
    }
  }
  
  if (is<List>(mpm)) {
    List mpm_list = as<List>(mpm);
    int mat_type {0}; // 0: ahist/age; 1: hist; 2:agestage
    // Need to include a test to determine if lefkoMat object or list of matrices
    
    CharacterVector stage_terms_caps = {"Stages", "Stage-pairs", "Age-stages"};
    CharacterVector stage_terms_small = {"stage", "stage-pair", "age-stage"};
    
    CharacterVector list_names = mpm_list.attr("names");
    int all_found {0};
    for (int i = 0; i < list_names.length(); i++) {
      if (list_names(i) == "A") all_found++;
      if (list_names(i) == "ahstages") all_found++;
      if (list_names(i) == "hstages") all_found++;
      if (list_names(i) == "agestages") all_found++;
      if (list_names(i) == "labels") all_found++;
    }
    if (all_found < 5) pop_error("mpm", "a lefkoMat object, matrix, or list of matrices", "", 1);
    
    List amats = as<List>(mpm_list["A"]);
    DataFrame stageframe = as<DataFrame>(mpm_list["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm_list["hstages"]);
    DataFrame labels = as<DataFrame>(mpm_list["labels"]);
    DataFrame agestages = as<DataFrame>(mpm_list["agestages"]);
    
    StringVector poporder = as<StringVector>(labels["pop"]);
    int loysize = poporder.length();
    
    StringVector stages = as<StringVector>(stageframe["stage"]);
    IntegerVector stage_id = as<IntegerVector>(stageframe["stage_id"]);
    int num_stages = stages.length();
    
    if (hstages.length() > 1) {
      mat_type = 1;
    } else if (agestages.length() > 1) mat_type = 2;
    
    if (mat_type == 1) {
      IntegerVector sid2_hstages = as<IntegerVector>(hstages["stage_id_2"]);
      num_stages = sid2_hstages.length();
    } else if (mat_type == 2) {
      IntegerVector sid_agestages = as<IntegerVector>(agestages["stage_id"]);
      num_stages = sid_agestages.length();
    }
    
    List list_of_no_in (loysize);
    List list_of_no_out (loysize);
    
    arma::uvec mats_probs_a_arma (loysize, fill::zeros);
    arma::uvec mats_probs_b_arma (loysize, fill::zeros);
    
    for (int mat = 0; mat < loysize; mat++) {
      arma::uvec stage_check (num_stages, fill::zeros); // Stages without inward transitions (as 0)
      arma::uvec problem_vector (num_stages, fill::zeros); // Stages without outward transitions (as 1)
      
      if (is<NumericMatrix>(amats(mat))) {
        arma::mat current_mat = as<arma::mat>(amats(mat));
        int num_cols = current_mat.n_cols;
        
        for (int current_col = 0; current_col < num_cols; current_col++) {
          arma::vec check_col = current_mat.col(current_col);
          arma::uvec nonzeros_vec = find(check_col);
          
          int nonzeros_num = nonzeros_vec.n_elem;
          if (check_col(current_col) != 0.) nonzeros_num -= 1;
          
          if (nonzeros_num == 0) {
            problem_vector(current_col) = 1; // No out}
            mats_probs_a_arma(mat) = 1;
          }
          
          arma::mat cmat_nocurrent = current_mat;
          cmat_nocurrent.shed_col(current_col);
          arma::vec cmat_nc_sums = sum(cmat_nocurrent, 1);
          
          if (cmat_nc_sums(current_col) > 0.) stage_check(current_col) += 1;
          
        }
        
      } else if (is<S4>(amats(mat))) {
        arma::sp_mat current_mat = as<arma::sp_mat>(amats(mat));
        int num_cols = current_mat.n_cols;
        
        if (num_cols != num_stages) {
          throw Rcpp::exception("This function currently handles only ahistorical and age-based MPMs", false);
        }
        
        for (int current_col = 0; current_col < num_cols; current_col++) {
          arma::vec check_col = arma::vec(current_mat.col(current_col));
          arma::uvec nonzeros_vec = find(check_col);
          
          int nonzeros_num = nonzeros_vec.n_elem;
          if (check_col(current_col) != 0.) nonzeros_num -= 1;
          
          if (nonzeros_num == 0) problem_vector(current_col) = 1; // No out
          
          arma::sp_mat cmat_nocurrent = current_mat;
          cmat_nocurrent.shed_col(current_col);
          arma::vec cmat_nc_sums = arma::vec(sum(cmat_nocurrent, 1));
          
          if (cmat_nc_sums(current_col) > 0.) {
            stage_check(current_col) += 1;
            mats_probs_b_arma(mat) = 1;
          }
          
        }
        
      }
      
      arma::uvec probs_a_arma = find(problem_vector); // No out
      arma::uvec probs_b_arma = find(stage_check == 0); // No in
      
      probs_a_arma += 1;
      probs_b_arma += 1;
      
      IntegerVector unique_probs_a = as<IntegerVector>(wrap(probs_a_arma));
      IntegerVector unique_probs_b = as<IntegerVector>(wrap(probs_b_arma));
      
      list_of_no_out(mat) = unique_probs_a;
      list_of_no_in(mat) = unique_probs_b;
    }
    
    if (!quiet_bool) {
      arma::uvec find_em_a = find(mats_probs_a_arma);
      arma::uvec find_em_b = find(mats_probs_b_arma);
      
      find_em_a += 1;
      find_em_b += 1;
      
      if (find_em_a.n_elem > 0) {
        Rcout << stage_terms_caps(mat_type)
          << " without connections leading to the rest of the life cycle "
          << "found in matrices: " << as<IntegerVector>(wrap(find_em_a)) << endl;
      }
      if (find_em_b.n_elem > 0) {
        Rcout << stage_terms_caps(mat_type)
          << " without connections leading to them found in matrices: "  
          << as<IntegerVector>(wrap(find_em_b)) << endl;
      }
      
      if (find_em_a.n_elem == 0 && find_em_b.n_elem == 0) {
        Rcout << "All matrices have no " << stage_terms_small(mat_type)
          << " discontinuities." << endl;
      }
    }
    
    lopv(0) = list_of_no_in;
    lopv(1) = list_of_no_out;
    
  } else if (is<NumericMatrix>(mpm)) {
    NumericMatrix mpm_matrix = as<NumericMatrix>(mpm);
    int num_stages = mpm_matrix.ncol();
    List list_of_no_in (1);
    List list_of_no_out (1);
    
    arma::uvec stage_check (num_stages, fill::zeros); // Stages without inward transitions (as 0)
    arma::uvec problem_vector (num_stages, fill::zeros); // Stages without outward transitions (as 1)
    
    arma::mat current_mat = as<arma::mat>(mpm_matrix);
    int num_cols = current_mat.n_cols;
    
    for (int current_col = 0; current_col < num_cols; current_col++) {
      arma::vec check_col = current_mat.col(current_col);
      arma::uvec nonzeros_vec = find(check_col);
      
      int nonzeros_num = nonzeros_vec.n_elem;
      if (check_col(current_col) != 0.) nonzeros_num -= 1;
      
      if (nonzeros_num == 0) problem_vector(current_col) = 1; // No out
      
      arma::mat cmat_nocurrent = current_mat;
      cmat_nocurrent.shed_col(current_col);
      arma::vec cmat_nc_sums = sum(cmat_nocurrent, 1);
      
      if (cmat_nc_sums(current_col) > 0.) stage_check(current_col) += 1;
      
    }
    
    arma::uvec probs_a_arma = find(problem_vector); // No out
    arma::uvec probs_b_arma = find(stage_check == 0); // No in
    
    probs_a_arma += 1;
    probs_b_arma += 1;
    
    IntegerVector unique_probs_a = as<IntegerVector>(wrap(probs_a_arma));
    IntegerVector unique_probs_b = as<IntegerVector>(wrap(probs_b_arma));
    
    list_of_no_out(0) = unique_probs_a;
    list_of_no_in(0) = unique_probs_b;
    
    if (!quiet_bool) {
      if (probs_a_arma.n_elem > 0) {
        Rcout << "Stages without connections leading to the rest of the "  
          << "life cycle found in matrix" << endl;
      }
      if (probs_b_arma.n_elem > 0) {
        Rcout << "Stages without connections leading to them found in matrix"
          << endl;
      }
      
      if (probs_a_arma.n_elem == 0 && probs_b_arma.n_elem == 0) {
        Rcout << "Matrix has no stage discontinuities." << endl;
      }
    }
    
    lopv(0) = list_of_no_in;
    lopv(1) = list_of_no_out;
    
  } else if (is<S4>(mpm)) {
    arma::sp_mat current_mat = as<arma::sp_mat>(mpm);
    int num_stages = current_mat.n_cols;
    int num_cols = num_stages;
    List list_of_no_in (1);
    List list_of_no_out (1);
    
    arma::uvec stage_check (num_stages, fill::zeros); // Stages without inward transitions (as 0)
    arma::uvec problem_vector (num_stages, fill::zeros); // Stages without outward transitions (as 1)
    
    for (int current_col = 0; current_col < num_cols; current_col++) {
      arma::vec check_col = arma::vec(current_mat.col(current_col));
      arma::uvec nonzeros_vec = find(check_col);
      
      int nonzeros_num = nonzeros_vec.n_elem;
      if (check_col(current_col) != 0.) nonzeros_num -= 1;
      
      if (nonzeros_num == 0) problem_vector(current_col) = 1; // No out
      
      arma::sp_mat cmat_nocurrent = current_mat;
      cmat_nocurrent.shed_col(current_col);
      arma::vec cmat_nc_sums = arma::vec(sum(cmat_nocurrent, 1));
      
      if (cmat_nc_sums(current_col) > 0.) stage_check(current_col) += 1;
      
    }
    
    arma::uvec probs_a_arma = find(problem_vector); // No out
    arma::uvec probs_b_arma = find(stage_check == 0); // No in
    
    probs_a_arma += 1;
    probs_b_arma += 1;
    
    IntegerVector unique_probs_a = as<IntegerVector>(wrap(probs_a_arma));
    IntegerVector unique_probs_b = as<IntegerVector>(wrap(probs_b_arma));
    
    list_of_no_out(0) = unique_probs_a;
    list_of_no_in(0) = unique_probs_b;
    
    if (!quiet_bool) {
      if (probs_a_arma.n_elem > 0) {
        Rcout << "Stages without connections leading to the rest of the "  
          << "life cycle found in matrix" << endl;
      }
      if (probs_b_arma.n_elem > 0) {
        Rcout << "Stages without connections leading to them found in matrix" << endl;
      }
      
      if (probs_a_arma.n_elem == 0 && probs_b_arma.n_elem == 0) {
        Rcout << "Matrix has no stage discontinuities." << endl;
      }
    }
    
    lopv(0) = list_of_no_in;
    lopv(1) = list_of_no_out;
    
  } else {
    pop_error("mpm", "a lefkoMat object, matrix, or list of matrices", "", 1);
  }
  
  CharacterVector lopv_names {"no_in", "no_out"};
  lopv.attr("names") = lopv_names;
  
  return lopv;
}

