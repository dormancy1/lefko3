// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;


// Index of functions
// 
// 1. List pfj  Create Vertical Structure for Horizontal Data Frame Input
// 2. List jpf  Create Historical Vertical Structure for Ahistorical Vertical Data Frame
// 3. NumericVector density3  Estimate Radial Density in Cartesian Space
// 4. List bootstrap3  Bootstrap Standardized hfv_data Datasets




//' Create Vertical Structure for Horizontal Data Frame Input
//' 
//' Function \code{pfj()} powers the R function \code{\link{verticalize3}()},
//' creating the vertical structure and rearranging the data in that shape.
//' 
//' @name .pfj
//' 
//' @param data The horizontal data file.
//' @param stageframe The stageframe object describing the life history model.
//' This should be the full stageframe.
//' @param noyears The number of years or observation periods in the dataset.
//' @param firstyear The first year or time of observation.
//' @param popidcol Column number corresponding to the identity of the
//' population for each individual.
//' @param patchidcol Column number corresponding to the identity of the patch
//' for each individual.
//' @param individcol Column number corresponding to the identity of each 
//' individual.
//' @param blocksize The number of variables corresponding to each time step in 
//' the input dataset designated in \code{data}.
//' @param xcol Vector of column numbers corresponding to the x coordinate of
//' each individual in Cartesian space.
//' @param ycol Vector of column numbers corresponding to the y coordinate of
//' each individual in Cartesian space.
//' @param juvcol Vector of column numbers that marks individuals in immature
//' stages within the dataset.
//' @param sizeacol Vector of column numbers corresponding to the first or main
//' size variable associated with the first year or observation time in the
//' dataset.
//' @param sizebcol Vector of column numbers corresponding to the second size
//' variable associated with the first year or observation time in the dataset.
//' @param sizeccol Vector of column numbers corresponding to the third size
//' variable associated with the first year or observation time in the dataset.
//' @param repstracol Vector of column numbers corresponding to the main 
//' variable coding the production of reproductive structures associated with
//' the first year or observation period in the input dataset.
//' @param repstrbcol Vector of column numbers corresponding to a second
//' variable coding the production of reproductive structures associated with
//' the first year or observation period in the input dataset.
//' @param fecacol Vector of column numbers corresponding to the main variable
//' coding for fecundity associated with the first year or observation period in
//' the dataset.
//' @param fecbcol Vector of column numbers corresponding to a second variable
//' coding for fecundity associated with the first year or observation period in
//' the dataset.
//' @param indcovacol Vector of column numbers corresponding to an individual
//' covariate.
//' @param indcovbcol Vector of column numbers corresponding to an individual
//' covariate.
//' @param indcovccol Vector of column numbers corresponding to an individual
//' covariate.
//' @param aliveacol Vector of column numbers that details whether an individual
//' is alive at a given time.
//' @param deadacol Vector of column numbers that details whether an individual
//' is dead at a given time.
//' @param obsacol Vector of column numbers that details whether an individual
//' is in an observable stage at a given time.
//' @param nonobsacol Vector of column numbers that details whether an
//' individual is in an unobservable stage at a given time.
//' @param censorcol Vector of column numbers corresponding to the first entry
//' of a censor variable.
//' @param stagecol Vector of column numbers corresponding to the first entry of
//' a column designating stages.
//' @param repstrrel This is a scalar modifier for that makes the variable in
//' \code{repstrbcol} equivalent to \code{repstracol}.
//' @param fecrel This is a scalar modifier for that makes the variable in
//' \code{fecbcol} equivalent to \code{fecacol}.
//' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
//' will be set to 0.
//' @param NRasRep If TRUE, then will treat non-reproductive but mature
//' individuals as reproductive during stage assignment.
//' @param RepasObs If TRUE, then will treat individuals with size 0 as observed
//' if and only if they are reproductive. Otherwise, all individuals with size 0
//' are treated as not observed.
//' @param NOasObs If TRUE, then will treat unobserved individuals as observed
//' during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Integer describing which size variable or combination of size
//' variables to use in stage estimation.
//' @param censorkeep The value of the censoring variable identifying data
//' that should be included in analysis. Defaults to 0, but may take any numeric
//' value including NA.
//' @param censbool A logical variable determining whether NA denotes the value
//' of the censoring variable identifying data to keep. If used, then will set
//' all NAs to 0 and all other values to 1, treating 0 as the value to keep.
//' @param censrepeat A logical value indicating whether censor variable is a
//' single static column, or whether censor variables repeat across blocks.
//' @param coordsrepeat A logical value indicating whether coordinate variables
//' are single static columns, or whether they repeat across blocks.
//' @param retain_alive0 A logical variable indicating whether to keep or remove
//' data rows for individuals not alive in time \emph{t}.
//' @param reduce A logical variable determining whether unused variables and
//' some invariant state variables should be removed from the output dataset.
//' Defaults to \code{TRUE}.
//' @param quiet A logical value indicating whether to silense warnings.
//' 
//' @return The output is currently a 7 element list, where each element is a
//' data frame with the same number of rows.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.pfj)]]
Rcpp::List pfj(const DataFrame& data, const DataFrame& stageframe,
  const int noyears, int firstyear, const int popidcol, const int patchidcol,
  const int individcol, const int blocksize, arma::ivec xcol,
  arma::ivec ycol, arma::ivec juvcol, arma::ivec sizeacol, arma::ivec sizebcol,
  arma::ivec sizeccol, arma::ivec repstracol, arma::ivec repstrbcol,
  arma::ivec fecacol, arma::ivec fecbcol, arma::ivec indcovacol,
  arma::ivec indcovbcol, arma::ivec indcovccol, arma::ivec aliveacol,
  arma::ivec deadacol, arma::ivec obsacol, arma::ivec nonobsacol,
  arma::ivec censorcol, arma::ivec stagecol, double repstrrel, double fecrel,
  bool NAas0, bool NRasRep, bool RepasObs, bool NOasObs, bool stassign,
  int stszcol, double censorkeep, bool censbool, bool censrepeat,
  bool coordsrepeat, bool retain_alive0, bool reduce, bool quiet) {
  
  Rcpp::List output_longlist (81);
  int noindivs = data.nrows();
  int novars = data.length();
  
  // Vector length checks and standardizing vectors of data frame columns
  if (xcol.n_elem == 1) {
    xcol.resize(noyears);
    
    if (!coordsrepeat || blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        xcol(i) = xcol(0);
      }
    } else if (xcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        xcol(i) = xcol(i-1) + blocksize;
        if (xcol(i) >= novars) {
          String eat_my_shorts = "Vector xcol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      xcol.fill(-1);
    }
  }
  
  if (ycol.n_elem == 1) {
    ycol.resize(noyears);
    
    if (!coordsrepeat || blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        ycol(i) = ycol(0);
      }
    } else if (ycol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        ycol(i) = ycol(i-1) + blocksize;
        if (ycol(i) >= novars) {
          String eat_my_shorts = "Vector ycol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      ycol.fill(-1);
    }
  }
  
  if (juvcol.n_elem == 1) {
    juvcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        juvcol(i) = juvcol(0);
      }
    } else if (juvcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        juvcol(i) = juvcol(i-1) + blocksize;
        if (juvcol(i) >= novars) {
          String eat_my_shorts = "Vector juvcol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      juvcol.fill(-1);
    }
  }
  
  if (sizeacol.n_elem == 1) {
    sizeacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizeacol(i) = sizeacol(0);
      }
    } else if (sizeacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        sizeacol(i) = sizeacol(i-1) + blocksize;
        if (sizeacol(i) >= novars) {
          String eat_my_shorts = "Vector sizeacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      sizeacol.fill(-1);
    }
  }
  
  if (sizebcol.n_elem == 1) {
    sizebcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizebcol(i) = sizebcol(0);
      }
    } else if (sizebcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        sizebcol(i) = sizebcol(i-1) + blocksize;
        if (sizebcol(i) >= novars) {
          String eat_my_shorts = "Vector sizebcol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      sizebcol.fill(-1);
    }
  }
  
  if (sizeccol.n_elem == 1) {
    sizeccol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizeccol(i) = sizeccol(0);
      }
    } else if (sizeccol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        sizeccol(i) = sizeccol(i-1) + blocksize;
        if (sizeccol(i) >= novars) {
          String eat_my_shorts = "Vector sizeccol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      sizeccol.fill(-1);
    }
  }
  
  if (repstracol.n_elem == 1) {
    repstracol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        repstracol(i) = repstracol(0);
      }
    } else if (repstracol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        repstracol(i) = repstracol(i-1) + blocksize;
        if (repstracol(i) >= novars) {
          String eat_my_shorts = "Vector repstracol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      repstracol.fill(-1);
    }
  }
  
  if (repstrbcol.n_elem == 1) {
    repstrbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        repstrbcol(i) = repstrbcol(0);
      }
    } else if (repstrbcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        repstrbcol(i) = repstrbcol(i-1) + blocksize;
        if (repstrbcol(i) >= novars) {
          String eat_my_shorts = "Vector repstrbcol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      repstrbcol.fill(-1);
    }
  }
  
  if (fecacol.n_elem == 1) {
    fecacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        fecacol(i) = fecacol(0);
      }
    } else if (fecacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        fecacol(i) = fecacol(i-1) + blocksize;
        if (fecacol(i) >= novars) {
          String eat_my_shorts = "Vector fecacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      fecacol.fill(-1);
    }
  }
  
  if (fecbcol.n_elem == 1) {
    fecbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        fecbcol(i) = fecbcol(0);
      }
    } else if (fecbcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        fecbcol(i) = fecbcol(i-1) + blocksize;
        if (fecbcol(i) >= novars) {
          String eat_my_shorts = "Vector fecbcol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      fecbcol.fill(-1);
    }
  }
  
  if (indcovacol.n_elem == 1) {
    indcovacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovacol(i) = indcovacol(0);
      }
    } else if (indcovacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        indcovacol(i) = indcovacol(i-1) + blocksize;
        if (indcovacol(i) >= novars) {
          String eat_my_shorts = "Vector indcovacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      indcovacol.fill(-1);
    }
  }
  
  if (indcovbcol.n_elem == 1) {
    indcovbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovbcol(i) = indcovbcol(0);
      }
    } else if (indcovbcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        indcovbcol(i) = indcovbcol(i-1) + blocksize;
        if (indcovbcol(i) >= novars) {
          String eat_my_shorts = "Vector indcovbcol is too small. ither noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      indcovbcol.fill(-1);
    }
  }
  
  if (indcovccol.n_elem == 1) {
    indcovccol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovccol(i) = indcovccol(0);
      }
    } else if (indcovccol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        indcovccol(i) = indcovccol(i-1) + blocksize;
        if (indcovccol(i) >= novars) {
          String eat_my_shorts = "Vector indcovccol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      indcovccol.fill(-1);
    }
  }
  
  if (aliveacol.n_elem == 1) {
    aliveacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        aliveacol(i) = aliveacol(0);
      }
    } else if (aliveacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        aliveacol(i) = aliveacol(i-1) + blocksize;
        if (aliveacol(i) >= novars) {
          String eat_my_shorts = "Vector aliveacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      aliveacol.fill(-1);
    }
  }
  
  if (deadacol.n_elem == 1) {
    deadacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        deadacol(i) = deadacol(0);
      }
    } else if (deadacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        deadacol(i) = deadacol(i-1) + blocksize;
        if (deadacol(i) >= novars) {
          String eat_my_shorts = "Vector deadacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      deadacol.fill(-1);
    }
  }
  
  if (obsacol.n_elem == 1) {
    obsacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        obsacol(i) = obsacol(0);
      }
    } else if (obsacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        obsacol(i) = obsacol(i-1) + blocksize;
        if (obsacol(i) >= novars) {
          String eat_my_shorts = "Vector obsacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      obsacol.fill(-1);
    }
  }
  
  if (nonobsacol.n_elem == 1) {
    nonobsacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        nonobsacol(i) = nonobsacol(0);
      }
    } else if (nonobsacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        nonobsacol(i) = nonobsacol(i-1) + blocksize;
        if (nonobsacol(i) >= novars) {
          String eat_my_shorts = "Vector nonobsacol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      nonobsacol.fill(-1);
    }
  }
  
  if (censorcol.n_elem == 1) {
    censorcol.resize(noyears);
    
    if (blocksize == 0 || !censrepeat) {
      for (int i = 1; i < noyears; i++) {
        censorcol(i) = censorcol(0);
      }
    } else if (blocksize > 0 && censorcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        censorcol(i) = censorcol(i-1) + blocksize;
        if (censorcol(i) >= novars) {
          String eat_my_shorts = "Vector censorcol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      censorcol.fill(-1);
    }
  }
  
  if (stagecol.n_elem == 1) {
    stagecol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        stagecol(i) = stagecol(0);
      }
    } else if (stagecol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        stagecol(i) = stagecol(i-1) + blocksize;
        if (stagecol(i) >= novars) {
          String eat_my_shorts = "Vector stagecol is too small. Either noyears is too small or not ";
          String eat_my_shorts1 = "enough data blocks have been input.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      }
    } else {
      stagecol.fill(-1);
    }
  }
  
  // Variable initialization for most vectors
  double livcheck1 {0.0};
  double livcheck2 {0.0};
  double livcheck3 {0.0};
  
  double nonobssub1 {0.0};
  double nonobssub2 {0.0};
  double nonobssub3 {0.0};
  
  double stagesize1 {0.0};
  double stagesize2 {0.0};
  double stagesize3 {0.0};
  double stagesize1b {0.0};
  double stagesize2b {0.0};
  double stagesize3b {0.0};
  double stagesize1c {0.0};
  double stagesize2c {0.0};
  double stagesize3c {0.0};
  
  bool indcova_as_int {false};
  bool indcovb_as_int {false};
  bool indcovc_as_int {false};
  
  arma::uvec stagemini1;
  arma::uvec stagemaxi1;
  arma::uvec stagemini2;
  arma::uvec stagemaxi2;
  arma::uvec stagemini3;
  arma::uvec stagemaxi3;
  arma::uvec stageobs;
  arma::uvec stagerep;
  arma::uvec stagemat;
  arma::uvec cs1;
  arma::uvec cs2;
  arma::uvec cs3;
  arma::uvec cs4;
  int choicestage;
  
  Rcpp::StringVector sfname = as<StringVector>(stageframe["stage"]);
  Rcpp::NumericVector repstat = as<NumericVector>(stageframe["repstatus"]);
  Rcpp::NumericVector obsstat = as<NumericVector>(stageframe["obsstatus"]);
  Rcpp::NumericVector matstat = as<NumericVector>(stageframe["matstatus"]);
  arma::vec indataset = as<arma::vec>(stageframe["indataset"]);
  arma::vec sfszmin = as<arma::vec>(stageframe["sizebin_min"]);
  arma::vec sfszmax = as<arma::vec>(stageframe["sizebin_max"]);
  arma::vec sfszminb = as<arma::vec>(stageframe["sizebinb_min"]);
  arma::vec sfszmaxb = as<arma::vec>(stageframe["sizebinb_max"]);
  arma::vec sfszminc = as<arma::vec>(stageframe["sizebinc_min"]);
  arma::vec sfszmaxc = as<arma::vec>(stageframe["sizebinc_max"]);
  
  arma::vec repstatarma = as<arma::vec>(repstat);
  arma::vec obsstatarma = as<arma::vec>(obsstat);
  arma::vec matstatarma = as<arma::vec>(matstat);
  arma::vec sfszminarma = sfszmin;
  arma::vec sfszmaxarma = sfszmax;
  arma::vec sfszminarmab = sfszminb;
  arma::vec sfszmaxarmab = sfszmaxb;
  arma::vec sfszminarmac = sfszminc;
  arma::vec sfszmaxarmac = sfszmaxc;
  int stagenum = static_cast<int>(sfszmaxarma.n_elem);
  
  arma::uvec instages = find(indataset == 1);
  int instagenum = static_cast<int>(instages.n_elem);
  
  arma::uvec stageid (stagenum);
  arma::uvec instageid (instagenum);
  arma::vec insfszminarma (instagenum);
  arma::vec insfszmaxarma (instagenum);
  arma::vec inrepstatarma (instagenum);
  arma::vec inobsstatarma (instagenum);
  arma::vec inmatstatarma (instagenum);
  arma::vec insfszminarmab (instagenum);
  arma::vec insfszmaxarmab (instagenum);
  arma::vec insfszminarmac (instagenum);
  arma::vec insfszmaxarmac (instagenum);

  int inplace {0};
  for (int i = 0; i < stagenum; i++) {
    stageid(i) = i+1;
    
    if (indataset(i) == 1) {
      instageid(inplace) = i + 1;
      insfszminarma(inplace) = sfszminarma(i);
      insfszmaxarma(inplace) = sfszmaxarma(i);
      insfszminarmab(inplace) = sfszminarmab(i);
      insfszmaxarmab(inplace) = sfszmaxarmab(i);
      insfszminarmac(inplace) = sfszminarmac(i);
      insfszmaxarmac(inplace) = sfszmaxarmac(i);
      inrepstatarma(inplace) = repstatarma(i);
      inobsstatarma(inplace) = obsstatarma(i);
      inmatstatarma(inplace) = matstatarma(i);
      
      inplace++;
    }
  }
  
  Rcpp::StringVector popidx (noindivs);
  Rcpp::StringVector patchidx (noindivs);
  Rcpp::StringVector individx (noindivs);
  
  Rcpp::IntegerVector popidx_int (noindivs);
  Rcpp::IntegerVector patchidx_int (noindivs);
  Rcpp::IntegerVector individx_int (noindivs);
  
  StringVector popid_class;
  StringVector patchid_class;
  StringVector individ_class;
  
  StringVector popid_levels;
  StringVector patchid_levels;
  StringVector individ_levels;
  
  int popid_type {0}; // 0: string; 1: integer, 2: factor
  int patchid_type {0}; // 0: string; 1: integer, 2: factor
  int individ_type {0}; // 0: string; 1: integer, 2: factor
  
  Rcpp::StringVector stage1x (noindivs);
  Rcpp::StringVector stage2x (noindivs);
  Rcpp::StringVector stage3x (noindivs);
  
  Rcpp::NumericVector xpos1x (noindivs);
  Rcpp::NumericVector ypos1x (noindivs);
  Rcpp::NumericVector xpos2x (noindivs);
  Rcpp::NumericVector ypos2x (noindivs);
  Rcpp::NumericVector xpos3x (noindivs);
  Rcpp::NumericVector ypos3x (noindivs);
  
  Rcpp::NumericVector sizea1x (noindivs);
  Rcpp::NumericVector sizea2x (noindivs);
  Rcpp::NumericVector sizea3x (noindivs);
  Rcpp::NumericVector repstra1x (noindivs);
  Rcpp::NumericVector repstra2x (noindivs);
  Rcpp::NumericVector repstra3x (noindivs);
  Rcpp::NumericVector feca1x (noindivs);
  Rcpp::NumericVector feca2x (noindivs);
  Rcpp::NumericVector feca3x (noindivs);
  Rcpp::NumericVector sizeb1x (noindivs);
  Rcpp::NumericVector sizeb2x (noindivs);
  Rcpp::NumericVector sizeb3x (noindivs);
  Rcpp::NumericVector repstrb1x (noindivs);
  Rcpp::NumericVector repstrb2x (noindivs);
  Rcpp::NumericVector repstrb3x (noindivs);
  Rcpp::NumericVector fecb1x (noindivs);
  Rcpp::NumericVector fecb2x (noindivs);
  Rcpp::NumericVector fecb3x (noindivs);
  Rcpp::NumericVector sizec1x (noindivs);
  Rcpp::NumericVector sizec2x (noindivs);
  Rcpp::NumericVector sizec3x (noindivs);
  
  Rcpp::NumericVector indcova1x;
  Rcpp::NumericVector indcova2x;
  Rcpp::NumericVector indcova3x;
  Rcpp::NumericVector indcovb1x;
  Rcpp::NumericVector indcovb2x;
  Rcpp::NumericVector indcovb3x;
  Rcpp::NumericVector indcovc1x;
  Rcpp::NumericVector indcovc2x;
  Rcpp::NumericVector indcovc3x;
  
  Rcpp::IntegerVector indcova1x_int;
  Rcpp::IntegerVector indcova2x_int;
  Rcpp::IntegerVector indcova3x_int;
  Rcpp::IntegerVector indcovb1x_int;
  Rcpp::IntegerVector indcovb2x_int;
  Rcpp::IntegerVector indcovb3x_int;
  Rcpp::IntegerVector indcovc1x_int;
  Rcpp::IntegerVector indcovc2x_int;
  Rcpp::IntegerVector indcovc3x_int;
  
  Rcpp::NumericVector censor1x (noindivs);
  Rcpp::NumericVector censor2x (noindivs);
  Rcpp::NumericVector censor3x (noindivs);
  Rcpp::NumericVector alivegiven1x (noindivs);
  Rcpp::NumericVector alivegiven2x (noindivs);
  Rcpp::NumericVector alivegiven3x (noindivs);
  Rcpp::NumericVector deadgiven1x (noindivs);
  Rcpp::NumericVector deadgiven2x (noindivs);
  Rcpp::NumericVector deadgiven3x (noindivs);
  Rcpp::NumericVector obsgiven1x (noindivs);
  Rcpp::NumericVector obsgiven2x (noindivs);
  Rcpp::NumericVector obsgiven3x (noindivs);
  Rcpp::NumericVector nonobsgiven1x (noindivs, -1.0);
  Rcpp::NumericVector nonobsgiven2x (noindivs, -1.0); 
  Rcpp::NumericVector nonobsgiven3x (noindivs, -1.0);
  Rcpp::NumericVector juvgiven1x (noindivs);
  Rcpp::NumericVector juvgiven2x (noindivs);
  Rcpp::NumericVector juvgiven3x (noindivs);
  
  Rcpp::NumericVector zerovec (noindivs, 0.0);
  Rcpp::IntegerVector NAvec_int (noindivs, NA_INTEGER);
  Rcpp::NumericVector onevec (noindivs, 1.0);
  Rcpp::NumericVector negonevec (noindivs, -1.0);
  
  int ndflength = noindivs * (noyears - 1);
  
  Rcpp::IntegerVector rowid (ndflength);
  Rcpp::StringVector popid (ndflength);
  Rcpp::IntegerVector popid_int (ndflength);
  Rcpp::StringVector patchid (ndflength);
  Rcpp::IntegerVector patchid_int (ndflength);
  Rcpp::StringVector individ (ndflength);
  Rcpp::IntegerVector individ_int (ndflength);
  Rcpp::NumericVector year2 (ndflength);
  Rcpp::NumericVector xpos1 (ndflength);
  Rcpp::NumericVector ypos1 (ndflength);
  Rcpp::NumericVector xpos2 (ndflength); 
  Rcpp::NumericVector ypos2 (ndflength);
  Rcpp::NumericVector xpos3 (ndflength);
  Rcpp::NumericVector ypos3 (ndflength);
  
  Rcpp::NumericVector sizea1 (ndflength);
  Rcpp::NumericVector sizea2 (ndflength);
  Rcpp::NumericVector sizea3 (ndflength);
  Rcpp::NumericVector sizea10 (ndflength);
  Rcpp::NumericVector sizea20 (ndflength);
  Rcpp::NumericVector sizea30 (ndflength);
  Rcpp::NumericVector sizeb1 (ndflength);
  Rcpp::NumericVector sizeb2 (ndflength); 
  Rcpp::NumericVector sizeb3 (ndflength);
  Rcpp::NumericVector sizeb10 (ndflength);
  Rcpp::NumericVector sizeb20 (ndflength); 
  Rcpp::NumericVector sizeb30 (ndflength);
  Rcpp::NumericVector sizec1 (ndflength);
  Rcpp::NumericVector sizec2 (ndflength);
  Rcpp::NumericVector sizec3 (ndflength);
  Rcpp::NumericVector sizec10 (ndflength);
  Rcpp::NumericVector sizec20 (ndflength);
  Rcpp::NumericVector sizec30 (ndflength);
  
  Rcpp::NumericVector repstra1 (ndflength);
  Rcpp::NumericVector repstra2 (ndflength);
  Rcpp::NumericVector repstra3 (ndflength);
  Rcpp::NumericVector repstra10 (ndflength);
  Rcpp::NumericVector repstra20 (ndflength);
  Rcpp::NumericVector repstra30 (ndflength);
  Rcpp::NumericVector repstrb1 (ndflength);
  Rcpp::NumericVector repstrb2 (ndflength);
  Rcpp::NumericVector repstrb3 (ndflength);
  Rcpp::NumericVector repstrb10 (ndflength);
  Rcpp::NumericVector repstrb20 (ndflength);
  Rcpp::NumericVector repstrb30 (ndflength);
  
  Rcpp::NumericVector feca1 (ndflength);
  Rcpp::NumericVector feca2 (ndflength);
  Rcpp::NumericVector feca3 (ndflength);
  Rcpp::NumericVector feca10 (ndflength);
  Rcpp::NumericVector feca20 (ndflength);
  Rcpp::NumericVector feca30 (ndflength);
  Rcpp::NumericVector fecb1 (ndflength);
  Rcpp::NumericVector fecb2 (ndflength);
  Rcpp::NumericVector fecb3 (ndflength);
  Rcpp::NumericVector fecb10 (ndflength);
  Rcpp::NumericVector fecb20 (ndflength);
  Rcpp::NumericVector fecb30 (ndflength);
  
  Rcpp::NumericVector indcova1 (ndflength, 0.0);
  Rcpp::NumericVector indcova2 (ndflength, 0.0);
  Rcpp::NumericVector indcova3 (ndflength, 0.0);
  Rcpp::NumericVector indcovb1 (ndflength, 0.0);
  Rcpp::NumericVector indcovb2 (ndflength, 0.0);
  Rcpp::NumericVector indcovb3 (ndflength, 0.0);
  Rcpp::NumericVector indcovc1 (ndflength, 0.0);
  Rcpp::NumericVector indcovc2 (ndflength, 0.0);
  Rcpp::NumericVector indcovc3 (ndflength, 0.0);
  
  Rcpp::IntegerVector indcova1_int (ndflength, 0);
  Rcpp::IntegerVector indcova2_int (ndflength, 0);
  Rcpp::IntegerVector indcova3_int (ndflength, 0);
  Rcpp::IntegerVector indcovb1_int (ndflength, 0);
  Rcpp::IntegerVector indcovb2_int (ndflength, 0);
  Rcpp::IntegerVector indcovb3_int (ndflength, 0);
  Rcpp::IntegerVector indcovc1_int (ndflength, 0);
  Rcpp::IntegerVector indcovc2_int (ndflength, 0);
  Rcpp::IntegerVector indcovc3_int (ndflength, 0);
  
  StringVector indcova_class;
  StringVector indcovb_class;
  StringVector indcovc_class;
  StringVector indcova_levels;
  StringVector indcovb_levels;
  StringVector indcovc_levels;
  
  Rcpp::NumericVector censor1 (ndflength);
  Rcpp::NumericVector censor2 (ndflength);
  Rcpp::NumericVector censor3 (ndflength);
  Rcpp::NumericVector alivegiven1 (ndflength);
  Rcpp::NumericVector alivegiven2 (ndflength);
  Rcpp::NumericVector alivegiven3 (ndflength);
  Rcpp::NumericVector deadgiven1 (ndflength);
  Rcpp::NumericVector deadgiven2 (ndflength);
  Rcpp::NumericVector deadgiven3 (ndflength);
  Rcpp::NumericVector obsgiven1 (ndflength);
  Rcpp::NumericVector obsgiven2 (ndflength);
  Rcpp::NumericVector obsgiven3 (ndflength);
  Rcpp::NumericVector nonobsgiven1 (ndflength);
  Rcpp::NumericVector nonobsgiven2 (ndflength); 
  Rcpp::NumericVector nonobsgiven3 (ndflength);
  Rcpp::NumericVector juvgiven1 (ndflength);
  Rcpp::NumericVector juvgiven2 (ndflength);
  Rcpp::NumericVector juvgiven3 (ndflength);
  
  Rcpp::NumericVector addedsize1 (ndflength);
  Rcpp::NumericVector addedsize2 (ndflength);
  Rcpp::NumericVector addedsize3 (ndflength);
  Rcpp::NumericVector addedrepstr1 (ndflength);
  Rcpp::NumericVector addedrepstr2 (ndflength);
  Rcpp::NumericVector addedrepstr3 (ndflength);
  Rcpp::NumericVector addedfec1 (ndflength);
  Rcpp::NumericVector addedfec2 (ndflength);
  Rcpp::NumericVector addedfec3 (ndflength);
  
  Rcpp::NumericVector spryn1 (ndflength);
  Rcpp::NumericVector spryn2 (ndflength);
  Rcpp::NumericVector spryn3 (ndflength);
  Rcpp::NumericVector repyn1 (ndflength);
  Rcpp::NumericVector repyn2 (ndflength);
  Rcpp::NumericVector repyn3 (ndflength);
  Rcpp::NumericVector fecyn1 (ndflength);
  Rcpp::NumericVector fecyn2 (ndflength);
  Rcpp::NumericVector fecyn3 (ndflength);
  Rcpp::NumericVector matstat1 (ndflength);
  Rcpp::NumericVector matstat2 (ndflength);
  Rcpp::NumericVector matstat3 (ndflength);
  
  Rcpp::NumericVector alive1 (ndflength);
  Rcpp::NumericVector alive2 (ndflength);
  Rcpp::NumericVector alive3 (ndflength);
  
  Rcpp::StringVector stage1 (ndflength);
  Rcpp::StringVector stage2 (ndflength);
  Rcpp::StringVector stage3 (ndflength);
  Rcpp::NumericVector stage1num (ndflength);
  Rcpp::NumericVector stage2num (ndflength);
  Rcpp::NumericVector stage3num (ndflength);
  
  Rcpp::NumericVector firstseenx (noindivs, 0.0);
  Rcpp::NumericVector lastseenx (noindivs, 0.0);
  Rcpp::NumericVector firstseen (ndflength);
  Rcpp::NumericVector lastseen (ndflength);
  Rcpp::NumericVector obsage (ndflength);
  Rcpp::NumericVector obslifespan (ndflength);
  
  for (int i = 0; i < noindivs; i++) {
    Rcpp::NumericVector temp1 (noindivs, 0.0);
    Rcpp::NumericVector temp2 (noindivs, 0.0);
    Rcpp::NumericVector temp3 (noindivs, 0.0);
    Rcpp::NumericVector temp4 (noindivs, 0.0);
    Rcpp::NumericVector temp5 (noindivs, 0.0);
    
    for (int k = 0; k < noyears; k++) {
      if (sizeacol(k) != -1 ) {
        temp1 = (as<NumericVector>(data[sizeacol(k)]));
      } else {temp1 = clone(zerovec);}
      if (sizebcol(k) != -1 ) {
        temp2 = (as<NumericVector>(data[sizebcol(k)]));
      } else {temp2 = clone(zerovec);}
      if (sizeccol(k) != -1 ) {
        temp3 = (as<NumericVector>(data[sizeccol(k)]));
      } else {temp3 = clone(zerovec);}
      if (repstracol(k) != -1 ) {
        temp4 = (as<NumericVector>(data[repstracol(k)]));
      } else {temp4 = clone(zerovec);}
      if (repstrbcol(k) != -1 ) {
        temp5 = (as<NumericVector>(data[repstrbcol(k)]));
      } else {temp5 = clone(zerovec);}
      
      double temp1_i {0.0};
      double temp2_i {0.0};
      double temp3_i {0.0};
      double temp4_i {0.0};
      double temp5_i {0.0};
      
      if (!NumericVector::is_na(temp1(i))) {temp1_i = temp1(i);}
      if (!NumericVector::is_na(temp2(i))) {temp2_i = temp2(i);}
      if (!NumericVector::is_na(temp3(i))) {temp3_i = temp3(i);}
      if (!NumericVector::is_na(temp4(i))) {temp4_i = temp4(i);}
      if (!NumericVector::is_na(temp5(i))) {temp5_i = temp5(i);}
      
      if (temp1_i > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp2_i > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp3_i > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp4_i > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp5_i > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      }
      
      if (temp1_i > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp2_i > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp3_i > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp4_i > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp5_i > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      }
    }
  }
  
  // Main loop
  for (int j = 0; j < (noyears - 1); j++) {
    Rcpp::checkUserInterrupt();
    
    if (popidcol > -1) {
      RObject popid_test = as<RObject>(data[popidcol]);
      
      if (is<IntegerVector>(popid_test)) { 
        popid_type = 1;
        popidx_int = as<IntegerVector>(data[popidcol]);
        
        if (popidx_int.hasAttribute("class")) {
          popid_class = popidx_int.attr("class");
        }
        
        if (popidx_int.hasAttribute("levels")) { 
          popid_type = 2;
          popid_levels = popidx_int.attr("levels");
        }
      } else {
        popidx = as<StringVector>(data[popidcol]);
      }
    }
    
    if (patchidcol > -1) {
      RObject patchid_test = as<RObject>(data[patchidcol]);
      
      if (is<IntegerVector>(patchid_test)) { 
        patchid_type = 1;
        patchidx_int = as<IntegerVector>(data[patchidcol]);
        
        if (patchidx_int.hasAttribute("class")){
          patchid_class = patchidx_int.attr("class");
        }
        
        if (patchidx_int.hasAttribute("levels")) { 
          patchid_type = 2;
          patchid_levels = patchidx_int.attr("levels");
        }
      } else {
        patchidx = as<StringVector>(data[patchidcol]);
      }
    }
    
    if (individcol > -1) {
      RObject individ_test = as<RObject>(data[individcol]);
      
      if (is<IntegerVector>(individ_test)) { 
        individ_type = 1;
        individx_int = as<IntegerVector>(data[individcol]);
        
        if (individx_int.hasAttribute("class")) {
          individ_class = individx_int.attr("class");
        }
        
        if (individx_int.hasAttribute("levels")) { 
          individ_type = 2;
          individ_levels = individx_int.attr("levels");
        }
      } else {
        individx = as<StringVector>(data[individcol]);
      }
    }
    
    if (j == 0) { // Implemented in the first year
      
      sizea1x = clone(zerovec);
      sizea2x = data[sizeacol(0)];
      sizea3x = data[sizeacol(1)];
      
      if (stagecol(0) > -1) {
        stage1x.fill("NotAlive");
        stage2x = data[stagecol(0)];
        stage3x = data[stagecol(1)];
      }
      
      if (xcol(0) > -1) {
        xpos1x = clone(zerovec);
        xpos2x = data[xcol(0)];
        xpos3x = data[xcol(1)];
      }
      
      if (ycol(0) > -1) {
        ypos1x = clone(zerovec);
        ypos2x = data[ycol(0)];
        ypos3x = data[ycol(1)];
      }
      
      if (repstracol(0) > -1) {
        repstra1x = clone(zerovec);
        repstra2x = data[repstracol(0)];
        repstra3x = data[repstracol(1)];
      }
      
      if (fecacol(0) > -1) {
        feca1x = clone(zerovec);
        feca2x = data[fecacol(0)];
        feca3x = data[fecacol(1)];
      }
      
      if (sizebcol(0) > -1) {
        sizeb1x = clone(zerovec);
        sizeb2x = data[sizebcol(0)];
        sizeb3x = data[sizebcol(1)];
      }
      
      if (repstrbcol(0) > -1) {
        repstrb1x = clone(zerovec);
        repstrb2x = data[repstrbcol(0)];
        repstrb3x = data[repstrbcol(1)];
      }
      
      if (fecbcol(0) > -1) {
        fecb1x = clone(zerovec);
        fecb2x = data[fecbcol(0)];
        fecb3x = data[fecbcol(1)];
      }
      
      if (sizeccol(0) > -1) {
        sizec1x = clone(zerovec);
        sizec2x = data[sizeccol(0)];
        sizec3x = data[sizeccol(1)];
      }
      
      if (indcovacol(0) > -1) {
        RObject indcova_test = as<RObject>(data[indcovacol(0)]);
        
        if (is<IntegerVector>(indcova_test)) {
          indcova3x_int = as<IntegerVector>(data[indcovacol(1)]);
          indcova2x_int = as<IntegerVector>(data[indcovacol(0)]);
          indcova1x_int = clone(NAvec_int);
          indcova_as_int = true;
          
          if (indcova2x_int.hasAttribute("class")) {
            indcova_class = indcova2x_int.attr("class");
            
            indcova1_int.attr("class") = indcova_class;
            indcova2_int.attr("class") = indcova_class;
            indcova3_int.attr("class") = indcova_class;
          }
          
          if (indcova2x_int.hasAttribute("levels")) {
            indcova_levels = indcova2x_int.attr("levels");
            
            indcova1_int.attr("levels") = indcova_levels;
            indcova2_int.attr("levels") = indcova_levels;
            indcova3_int.attr("levels") = indcova_levels;
          }
          
        } else if (is<NumericVector>(indcova_test)) {
          indcova3x = data[indcovacol(1)];
          indcova2x = data[indcovacol(0)];
          indcova1x = clone(zerovec);
          
        } else {
          throw Rcpp::exception("Individual covariate a must be a numeric, integer, or factor variable.", 
            false);
        }
      }
      
      if (indcovbcol(0) > -1) {
        RObject indcovb_test = as<RObject>(data[indcovbcol(0)]);
        
        if (is<IntegerVector>(indcovb_test)) {
          indcovb3x_int = as<IntegerVector>(data[indcovbcol(1)]);
          indcovb2x_int = as<IntegerVector>(data[indcovbcol(0)]);
          indcovb1x_int = clone(NAvec_int);
          indcovb_as_int = true;
          
          if (indcovb2x_int.hasAttribute("class")) {
            indcovb_class = indcovb2x_int.attr("class");
            
            indcovb1_int.attr("class") = indcovb_class;
            indcovb2_int.attr("class") = indcovb_class;
            indcovb3_int.attr("class") = indcovb_class;
          }
          
          if (indcovb2x_int.hasAttribute("levels")) {
            indcovb_levels = indcovb2x_int.attr("levels");
            
            indcovb1_int.attr("levels") = indcovb_levels;
            indcovb2_int.attr("levels") = indcovb_levels;
            indcovb3_int.attr("levels") = indcovb_levels;
          }
          
        } else if (is<NumericVector>(indcovb_test)) {
          indcovb3x = data[indcovbcol(1)];
          indcovb2x = data[indcovbcol(0)];
          indcovb1x = clone(zerovec);
          
        } else {
          throw Rcpp::exception("Individual covariate b must be a numeric, integer, or factor variable.", 
            false);
        }
      }
      
      if (indcovccol(0) > -1) {
        RObject indcovc_test = as<RObject>(data[indcovccol(0)]);
        
        if (is<IntegerVector>(indcovc_test)) {
          indcovc3x_int = as<IntegerVector>(data[indcovccol(1)]);
          indcovc2x_int = as<IntegerVector>(data[indcovccol(0)]);
          indcovc1x_int = clone(NAvec_int);
          indcovc_as_int = true;
          
          if (indcovc2x_int.hasAttribute("class")) {
            indcovc_class = indcovc2x_int.attr("class");
            
            indcovc1_int.attr("class") = indcovc_class;
            indcovc2_int.attr("class") = indcovc_class;
            indcovc3_int.attr("class") = indcovc_class;
          }
          
          if (indcovc2x_int.hasAttribute("levels")) {
            indcovc_levels = indcovc2x_int.attr("levels");
            
            indcovc1_int.attr("levels") = indcovc_levels;
            indcovc2_int.attr("levels") = indcovc_levels;
            indcovc3_int.attr("levels") = indcovc_levels;
          }
          
        } else if (is<NumericVector>(indcovc_test)) {
          indcovc3x = data[indcovccol(1)];
          indcovc2x = data[indcovccol(0)];
          indcovc1x = clone(zerovec);
          
        } else {
          throw Rcpp::exception("Individual covariate c must be a numeric, integer, or factor variable.", 
            false);
        }
      }
      
      // This section decides how to deal with censor variables
      if (censorcol(0) > -1) {
        
        if (censrepeat) {
        
          if (censorkeep == 0) {
            censor1x = clone(onevec);
          } else {
            censor1x = clone(zerovec);
          }
          
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(1)];
          
        } else {
        
          if (censorkeep == 0) {
            censor1x = clone(onevec);
          } else {
            censor1x = clone(zerovec);
          }
          
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(0)];
        }
      }
      
      if (aliveacol(0) > -1) {
        alivegiven1x = clone(zerovec);
        alivegiven2x = data[aliveacol(0)];
        alivegiven3x = data[aliveacol(1)];
      }
      
      if (deadacol(0) > -1) {
        deadgiven1x = clone(zerovec);
        deadgiven2x = data[deadacol(0)];
        deadgiven3x = data[deadacol(1)];
      }
      
      if (obsacol(0) > -1) {
        obsgiven1x = clone(zerovec);
        obsgiven2x = data[obsacol(0)];
        obsgiven3x = data[obsacol(1)];
      }
      
      if (nonobsacol(0) > -1) {
        nonobsgiven1x = negonevec;
        nonobsgiven2x = data[nonobsacol(0)];
        nonobsgiven3x = data[nonobsacol(1)];
      }
      
      if (juvcol(0) > -1) {
        juvgiven1x = clone(zerovec);
        juvgiven2x = data[juvcol(0)];
        juvgiven3x = data[juvcol(1)];
      }
      
    } else { // Establishes what gets done after the 1st year
      sizea1x = data[sizeacol(j - 1)];
      sizea2x = data[sizeacol(j)];
      sizea3x = data[sizeacol(j + 1)];
      
      if (stagecol(0) > -1) {
        stage1x = data[stagecol(j - 1)];
        stage2x = data[stagecol(j)];
        stage3x = data[stagecol(j + 1)];
      }
      
      if (xcol(0) > -1) {
        xpos1x = data[xcol(j - 1)];
        xpos2x = data[xcol(j)];
        xpos3x = data[xcol(j + 1)];
      }
      
      if (ycol(0) > -1) {
        ypos1x = data[ycol(j - 1)];
        ypos2x = data[ycol(j)];
        ypos3x = data[ycol(j + 1)];
      }
      
      if (repstracol(0) > -1) {
        repstra1x = data[repstracol(j - 1)];
        repstra2x = data[repstracol(j)];
        repstra3x = data[repstracol(j + 1)];
      }
      
      if (fecacol(0) > -1) {
        feca1x = data[fecacol(j - 1)];
        feca2x = data[fecacol(j)];
        feca3x = data[fecacol(j + 1)];
      }
      
      if (sizebcol(0) > -1) {
        sizeb1x = data[sizebcol(j - 1)];
        sizeb2x = data[sizebcol(j)];
        sizeb3x = data[sizebcol(j + 1)];
      }
      
      if (repstrbcol(0) > -1) {
        repstrb1x = data[repstrbcol(j - 1)];
        repstrb2x = data[repstrbcol(j)];
        repstrb3x = data[repstrbcol(j + 1)];
      }
      
      if (fecbcol(0) > -1) {
        fecb1x = data[fecbcol(j - 1)];
        fecb2x = data[fecbcol(j)];
        fecb3x = data[fecbcol(j + 1)];
      }
      
      if (sizeccol(0) > -1) {
        sizec1x = data[sizeccol(j - 1)];
        sizec2x = data[sizeccol(j)];
        sizec3x = data[sizeccol(j + 1)];
      }
      
      if (indcovacol(0) > -1) {
        if (indcova_as_int) {
          indcova3x_int = as<IntegerVector>(data[indcovacol(j + 1)]);
          indcova2x_int = as<IntegerVector>(data[indcovacol(j)]);
          indcova1x_int = as<IntegerVector>(data[indcovacol(j - 1)]);

        } else {
          indcova3x = as<NumericVector>(data[indcovacol(j + 1)]);
          indcova2x = as<NumericVector>(data[indcovacol(j)]);
          indcova1x = as<NumericVector>(data[indcovacol(j - 1)]);
        }
      }
      
      if (indcovbcol(0) > -1) {
        if (indcovb_as_int) {
          indcovb3x_int = as<IntegerVector>(data[indcovbcol(j + 1)]);
          indcovb2x_int = as<IntegerVector>(data[indcovbcol(j)]);
          indcovb1x_int = as<IntegerVector>(data[indcovbcol(j - 1)]);

        } else {
          indcovb3x = as<NumericVector>(data[indcovbcol(j + 1)]);
          indcovb2x = as<NumericVector>(data[indcovbcol(j)]);
          indcovb1x = as<NumericVector>(data[indcovbcol(j - 1)]);
        }
      }
      
      if (indcovccol(0) > -1) {
        if (indcovc_as_int) {
          indcovc3x_int = as<IntegerVector>(data[indcovccol(j + 1)]);
          indcovc2x_int = as<IntegerVector>(data[indcovccol(j)]);
          indcovc1x_int = as<IntegerVector>(data[indcovccol(j - 1)]);

        } else {
          indcovc3x = as<NumericVector>(data[indcovccol(j + 1)]);
          indcovc2x = as<NumericVector>(data[indcovccol(j)]);
          indcovc1x = as<NumericVector>(data[indcovccol(j - 1)]);
        }
      }
      
      // Censor variable handling
      if (censorcol(0) > -1) {
        if (censrepeat) {
          censor1x = data[censorcol(j - 1)];
          censor2x = data[censorcol(j)];
          censor3x = data[censorcol(j + 1)];
        } else {
          censor1x = data[censorcol(0)];
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(0)];
        }
      }
      
      if (aliveacol(0) > -1) {
        alivegiven1x = data[aliveacol(j - 1)];
        alivegiven2x = data[aliveacol(j)];
        alivegiven3x = data[aliveacol(j + 1)];
      }
      
      if (deadacol(0) > -1) {
        deadgiven1x = data[deadacol(j - 1)];
        deadgiven2x = data[deadacol(j)];
        deadgiven3x = data[deadacol(j + 1)];
      }
      
      if (obsacol(0) > -1) {
        obsgiven1x = data[obsacol(j - 1)];
        obsgiven2x = data[obsacol(j)];
        obsgiven3x = data[obsacol(j + 1)];
      }
      
      if (nonobsacol(0) > -1) {
        nonobsgiven1x = data[nonobsacol(j - 1)];
        nonobsgiven2x = data[nonobsacol(j)];
        nonobsgiven3x = data[nonobsacol(j + 1)];
      }
      
      if (juvcol(0) > -1) {
        juvgiven1x = data[juvcol(j - 1)];
        juvgiven2x = data[juvcol(j)];
        juvgiven3x = data[juvcol(j + 1)];
      }
    }
    
    for (int i = 0; i < noindivs; i++) {
      livcheck1 = 0.0;
      livcheck2 = 0.0;
      livcheck3 = 0.0;
      
      nonobssub1 = -1.0;
      nonobssub2 = -1.0;
      nonobssub3 = -1.0;
      
      rowid[(i + (j * noindivs))] = i + 1;
      
      if (popid_type > 0) { 
        popid_int[(i + (j * noindivs))] = popidx_int[i];
      } else { 
        popid[(i + (j * noindivs))] = popidx[i];
      }
      if (patchid_type > 0) { 
        patchid_int[(i + (j * noindivs))] = patchidx_int[i];
      } else { 
        patchid[(i + (j * noindivs))] = patchidx[i];
      }
      if (individ_type > 0) { 
        individ_int[(i + (j * noindivs))] = individx_int[i];
      } else { 
        individ[(i + (j * noindivs))] = individx[i];
      }
      
      year2[(i + (j * noindivs))] = static_cast<double>(firstyear) + static_cast<double>(j);
      xpos1[(i + (j * noindivs))] = xpos1x[i];
      ypos1[(i + (j * noindivs))] = ypos1x[i];
      xpos2[(i + (j * noindivs))] = xpos2x[i];
      ypos2[(i + (j * noindivs))] = ypos2x[i];
      xpos3[(i + (j * noindivs))] = xpos3x[i];
      ypos3[(i + (j * noindivs))] = ypos3x[i];
      
      sizea1[(i + (j * noindivs))] = sizea1x[i];
      sizea2[(i + (j * noindivs))] = sizea2x[i];
      sizea3[(i + (j * noindivs))] = sizea3x[i];
      sizea10[(i + (j * noindivs))] = sizea1x[i];
      sizea20[(i + (j * noindivs))] = sizea2x[i];
      sizea30[(i + (j * noindivs))] = sizea3x[i];
      sizeb1[(i + (j * noindivs))] = sizeb1x[i];
      sizeb2[(i + (j * noindivs))] = sizeb2x[i];
      sizeb3[(i + (j * noindivs))] = sizeb3x[i];
      sizeb10[(i + (j * noindivs))] = sizeb1x[i];
      sizeb20[(i + (j * noindivs))] = sizeb2x[i];
      sizeb30[(i + (j * noindivs))] = sizeb3x[i];
      sizec1[(i + (j * noindivs))] = sizec1x[i];
      sizec2[(i + (j * noindivs))] = sizec2x[i];
      sizec3[(i + (j * noindivs))] = sizec3x[i];
      sizec10[(i + (j * noindivs))] = sizec1x[i];
      sizec20[(i + (j * noindivs))] = sizec2x[i];
      sizec30[(i + (j * noindivs))] = sizec3x[i];
      
      repstra1[(i + (j * noindivs))] = repstra1x[i];
      repstra2[(i + (j * noindivs))] = repstra2x[i];
      repstra3[(i + (j * noindivs))] = repstra3x[i];
      repstra10[(i + (j * noindivs))] = repstra1x[i];
      repstra20[(i + (j * noindivs))] = repstra2x[i];
      repstra30[(i + (j * noindivs))] = repstra3x[i];
      repstrb1[(i + (j * noindivs))] = repstrb1x[i];
      repstrb2[(i + (j * noindivs))] = repstrb2x[i];
      repstrb3[(i + (j * noindivs))] = repstrb3x[i];
      repstrb10[(i + (j * noindivs))] = repstrb1x[i];
      repstrb20[(i + (j * noindivs))] = repstrb2x[i];
      repstrb30[(i + (j * noindivs))] = repstrb3x[i];
      
      feca1[(i + (j * noindivs))] = feca1x[i];
      feca2[(i + (j * noindivs))] = feca2x[i];
      feca3[(i + (j * noindivs))] = feca3x[i];
      feca10[(i + (j * noindivs))] = feca1x[i];
      feca20[(i + (j * noindivs))] = feca2x[i];
      feca30[(i + (j * noindivs))] = feca3x[i];
      fecb1[(i + (j * noindivs))] = fecb1x[i];
      fecb2[(i + (j * noindivs))] = fecb2x[i];
      fecb3[(i + (j * noindivs))] = fecb3x[i];
      fecb10[(i + (j * noindivs))] = fecb1x[i];
      fecb20[(i + (j * noindivs))] = fecb2x[i];
      fecb30[(i + (j * noindivs))] = fecb3x[i];
      
      if (!indcova_as_int && indcovacol(0) > -1) {
        indcova1[(i + (j * noindivs))] = indcova1x[i];
        indcova2[(i + (j * noindivs))] = indcova2x[i];
        indcova3[(i + (j * noindivs))] = indcova3x[i];
      } else if (indcovacol(0) > -1) {
        indcova1_int[(i + (j * noindivs))] = indcova1x_int[i];
        indcova2_int[(i + (j * noindivs))] = indcova2x_int[i];
        indcova3_int[(i + (j * noindivs))] = indcova3x_int[i];
      }
      
      if (!indcovb_as_int && indcovbcol(0) > -1) {
        indcovb1[(i + (j * noindivs))] = indcovb1x[i];
        indcovb2[(i + (j * noindivs))] = indcovb2x[i];
        indcovb3[(i + (j * noindivs))] = indcovb3x[i];
      } else if (indcovbcol(0) > -1) {
        indcovb1_int[(i + (j * noindivs))] = indcovb1x_int[i];
        indcovb2_int[(i + (j * noindivs))] = indcovb2x_int[i];
        indcovb3_int[(i + (j * noindivs))] = indcovb3x_int[i];
      }
      
      if (!indcovc_as_int && indcovccol(0) > -1) {
        indcovc1[(i + (j * noindivs))] = indcovc1x[i];
        indcovc2[(i + (j * noindivs))] = indcovc2x[i];
        indcovc3[(i + (j * noindivs))] = indcovc3x[i];
      } else if (indcovccol(0) > -1) {
        indcovc1_int[(i + (j * noindivs))] = indcovc1x_int[i];
        indcovc2_int[(i + (j * noindivs))] = indcovc2x_int[i];
        indcovc3_int[(i + (j * noindivs))] = indcovc3x_int[i];
      }
      
      if (stagecol(0) > -1) {
        stage1[(i + (j * noindivs))] = stage1x[i];
        stage2[(i + (j * noindivs))] = stage2x[i];
        stage3[(i + (j * noindivs))] = stage3x[i];
        
        stage1num[(i + (j * noindivs))] = 0.0;
        stage2num[(i + (j * noindivs))] = 0.0;
        stage3num[(i + (j * noindivs))] = 0.0;
      }
      
      if (censbool) { // Develops the censor variable
        if (NumericVector::is_na(censor1x[i])) {
          censor1[(i + (j * noindivs))] = 0.0;
        } else if (j != 0) {
          censor1[(i + (j * noindivs))] = 1.0;
        }
        
        if (NumericVector::is_na(censor2x[i])) {
          censor2[(i + (j * noindivs))] = 0.0;
        } else {
          censor2[(i + (j * noindivs))] = 1.0;
        }
        
        if (NumericVector::is_na(censor3x[i])) {
          censor3[(i + (j * noindivs))] = 0.0;
        } else {
          censor3[(i + (j * noindivs))] = 1.0;
        }
      } else {
        censor1[(i + (j * noindivs))] = censor1x[i];
        censor2[(i + (j * noindivs))] = censor2x[i];
        censor3[(i + (j * noindivs))] = censor3x[i];
      }
      
      alivegiven1[(i + (j * noindivs))] = alivegiven1x[i];
      alivegiven2[(i + (j * noindivs))] = alivegiven2x[i];
      alivegiven3[(i + (j * noindivs))] = alivegiven3x[i];
      deadgiven1[(i + (j * noindivs))] = deadgiven1x[i];
      deadgiven2[(i + (j * noindivs))] = deadgiven2x[i];
      deadgiven3[(i + (j * noindivs))] = deadgiven3x[i];
      obsgiven1[(i + (j * noindivs))] = obsgiven1x[i];
      obsgiven2[(i + (j * noindivs))] = obsgiven2x[i];
      obsgiven3[(i + (j * noindivs))] = obsgiven3x[i];
      nonobsgiven1[(i + (j * noindivs))] = nonobsgiven1x[i];
      nonobsgiven2[(i + (j * noindivs))] = nonobsgiven2x[i];
      nonobsgiven3[(i + (j * noindivs))] = nonobsgiven3x[i];
      juvgiven1[(i + (j * noindivs))] = juvgiven1x[i];
      juvgiven2[(i + (j * noindivs))] = juvgiven2x[i];
      juvgiven3[(i + (j * noindivs))] = juvgiven3x[i];
      
      if (NumericVector::is_na(sizea10[(i + (j * noindivs))])) sizea10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizea20[(i + (j * noindivs))])) sizea20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizea30[(i + (j * noindivs))])) sizea30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizeb10[(i + (j * noindivs))])) sizeb10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizeb20[(i + (j * noindivs))])) sizeb20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizeb30[(i + (j * noindivs))])) sizeb30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizec10[(i + (j * noindivs))])) sizec10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizec20[(i + (j * noindivs))])) sizec20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizec30[(i + (j * noindivs))])) sizec30[(i + (j * noindivs))] = 0.0;
      
      if (NumericVector::is_na(repstra10[(i + (j * noindivs))])) repstra10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstra20[(i + (j * noindivs))])) repstra20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstra30[(i + (j * noindivs))])) repstra30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstrb10[(i + (j * noindivs))])) repstrb10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstrb20[(i + (j * noindivs))])) repstrb20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstrb30[(i + (j * noindivs))])) repstrb30[(i + (j * noindivs))] = 0.0;
      
      if (NumericVector::is_na(feca10[(i + (j * noindivs))])) feca10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(feca20[(i + (j * noindivs))])) feca20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(feca30[(i + (j * noindivs))])) feca30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(fecb10[(i + (j * noindivs))])) fecb10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(fecb20[(i + (j * noindivs))])) fecb20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(fecb30[(i + (j * noindivs))])) fecb30[(i + (j * noindivs))] = 0.0;
      
      if (!indcova_as_int) {
        if (NumericVector::is_na(indcova1[(i + (j * noindivs))])) indcova1[(i + (j * noindivs))] = 0.0;
        if (NumericVector::is_na(indcova2[(i + (j * noindivs))])) indcova2[(i + (j * noindivs))] = 0.0;
        if (NumericVector::is_na(indcova3[(i + (j * noindivs))])) indcova3[(i + (j * noindivs))] = 0.0;
      }
      
      if (!indcovb_as_int) {
        if (NumericVector::is_na(indcovb1[(i + (j * noindivs))])) indcovb1[(i + (j * noindivs))] = 0.0;
        if (NumericVector::is_na(indcovb2[(i + (j * noindivs))])) indcovb2[(i + (j * noindivs))] = 0.0;
        if (NumericVector::is_na(indcovb3[(i + (j * noindivs))])) indcovb3[(i + (j * noindivs))] = 0.0;
      }
      
      if (!indcovc_as_int) {
        if (NumericVector::is_na(indcovc1[(i + (j * noindivs))])) indcovc1[(i + (j * noindivs))] = 0.0;
        if (NumericVector::is_na(indcovc2[(i + (j * noindivs))])) indcovc2[(i + (j * noindivs))] = 0.0;
        if (NumericVector::is_na(indcovc3[(i + (j * noindivs))])) indcovc3[(i + (j * noindivs))] = 0.0;
      }
      
      if (NumericVector::is_na(juvgiven1[(i + (j * noindivs))])) {
        juvgiven1[(i + (j * noindivs))] = 0.0;
      } else if (juvgiven1[(i + (j * noindivs))] != 0.0) {
        juvgiven1[(i + (j * noindivs))] = 1.0;
      }
      
      if (NumericVector::is_na(juvgiven2[(i + (j * noindivs))])) {
        juvgiven2[(i + (j * noindivs))] = 0.0;
      } else if (juvgiven2[(i + (j * noindivs))] != 0.0) {
        juvgiven2[(i + (j * noindivs))] = 1.0;
      }
      
      if (NumericVector::is_na(juvgiven3[(i + (j * noindivs))])) {
        juvgiven3[(i + (j * noindivs))] = 0.0;
      } else if (juvgiven3[(i + (j * noindivs))] != 0.0) {
        juvgiven3[(i + (j * noindivs))] = 1.0;
      }
      
      matstat1[(i + (j * noindivs))] = 1.0 - juvgiven1[(i + (j * noindivs))];
      matstat2[(i + (j * noindivs))] = 1.0 - juvgiven2[(i + (j * noindivs))];
      matstat3[(i + (j * noindivs))] = 1.0 - juvgiven3[(i + (j * noindivs))];
      
      addedsize1[(i + (j * noindivs))] = sizea10[(i + (j * noindivs))] + 
        sizeb10[(i + (j * noindivs))] + sizec10[(i + (j * noindivs))];
      addedsize2[(i + (j * noindivs))] = sizea20[(i + (j * noindivs))] + 
        sizeb20[(i + (j * noindivs))] + sizec20[(i + (j * noindivs))];
      addedsize3[(i + (j * noindivs))] = sizea30[(i + (j * noindivs))] + 
        sizeb30[(i + (j * noindivs))] + sizec30[(i + (j * noindivs))];
      
      addedrepstr1[(i + (j * noindivs))] = repstra10[(i + (j * noindivs))] + 
        (repstrb10[(i + (j * noindivs))] * repstrrel);
      addedrepstr2[(i + (j * noindivs))] = repstra20[(i + (j * noindivs))] + 
        (repstrb20[(i + (j * noindivs))] * repstrrel);
      addedrepstr3[(i + (j * noindivs))] = repstra30[(i + (j * noindivs))] + 
        (repstrb30[(i + (j * noindivs))] * repstrrel);
      
      addedfec1[(i + (j * noindivs))] = feca10[(i + (j * noindivs))] + 
        (fecb10[(i + (j * noindivs))] * fecrel);
      addedfec2[(i + (j * noindivs))] = feca20[(i + (j * noindivs))] + 
        (fecb20[(i + (j * noindivs))] * fecrel);
      addedfec3[(i + (j * noindivs))] = feca30[(i + (j * noindivs))] + 
        (fecb30[(i + (j * noindivs))] * fecrel);
      
      if (NumericVector::is_na(nonobsgiven1[(i + (j * noindivs))])) {
        nonobssub1 = -1.0;
      } else if (nonobsgiven1[(i + (j * noindivs))] >= 0) {
        nonobsgiven1[(i + (j * noindivs))] = 1.0;
        nonobssub1 = 1.0;
      }
      
      if (NumericVector::is_na(nonobsgiven2[(i + (j * noindivs))])) {
        nonobssub2 = -1.0;
      } else if (nonobsgiven2[(i + (j * noindivs))] >= 0) {
        nonobsgiven2[(i + (j * noindivs))] = 1.0;
        nonobssub2 = 1.0;
      }
      
      if (NumericVector::is_na(nonobsgiven3[(i + (j * noindivs))])) {
        nonobssub3 = -1.0;
      } else if (nonobsgiven3[(i + (j * noindivs))] >= 0) {
        nonobsgiven3[(i + (j * noindivs))] = 1.0;
        nonobssub3 = 1.0;
      }
      
      // Observation status
      if (addedsize1[(i + (j * noindivs))] > 0 || obsgiven1[(i + (j * noindivs))] > 0 ||
          nonobssub1 == 0.0) {
        spryn1[(i + (j * noindivs))] = 1.0;
      }
      if (addedsize2[(i + (j * noindivs))] > 0 || obsgiven2[(i + (j * noindivs))] > 0 ||
          nonobssub2 == 0.0) {
        spryn2[(i + (j * noindivs))] = 1.0;
      } 
      if (addedsize3[(i + (j * noindivs))] > 0 || obsgiven3[(i + (j * noindivs))] > 0 ||
          nonobssub3 == 0.0) {
        spryn3[(i + (j * noindivs))] = 1.0;
      }
      
      // Correction if size 0 individuals that reproduce are observed
      if (RepasObs == 1 && addedsize1[(i + (j * noindivs))] == 0 &&
          addedrepstr1[(i + (j * noindivs))] != 0.0) {
        spryn1[(i + (j * noindivs))] = 1.0;
      }
      if (RepasObs == 1 && addedsize2[(i + (j * noindivs))] == 0 &&
          addedrepstr2[(i + (j * noindivs))] != 0.0) {
        spryn2[(i + (j * noindivs))] = 1.0;
      }
      if (RepasObs == 1 && addedsize3[(i + (j * noindivs))] == 0 &&
          addedrepstr3[(i + (j * noindivs))] != 0.0) {
        spryn3[(i + (j * noindivs))] = 1.0;
      }
      
      if (addedrepstr1[(i + (j * noindivs))] > 0) repyn1[(i + (j * noindivs))] = 1.0;
      if (addedrepstr2[(i + (j * noindivs))] > 0) repyn2[(i + (j * noindivs))] = 1.0;
      if (addedrepstr3[(i + (j * noindivs))] > 0) repyn3[(i + (j * noindivs))] = 1.0;
      
      if (addedfec1[(i + (j * noindivs))] > 0) fecyn1[(i + (j * noindivs))] = 1.0;
      if (addedfec2[(i + (j * noindivs))] > 0) fecyn2[(i + (j * noindivs))] = 1.0;
      if (addedfec3[(i + (j * noindivs))] > 0) fecyn3[(i + (j * noindivs))] = 1.0;
      
      if (nonobssub1 >= 0) {
        livcheck1 = addedsize1[(i + (j * noindivs))] + addedrepstr1[(i + (j * noindivs))] + 
            spryn1[(i + (j * noindivs))] + nonobssub1;
      } else {
        livcheck1 = addedsize1[(i + (j * noindivs))] + addedrepstr1[(i + (j * noindivs))] + 
            spryn1[(i + (j * noindivs))];
      }
      if (nonobssub2 >= 0) {
        livcheck2 = addedsize2[(i + (j * noindivs))] + addedrepstr2[(i + (j * noindivs))] + 
            spryn2[(i + (j * noindivs))] + nonobssub2;
      } else {
        livcheck2 = addedsize2[(i + (j * noindivs))] + addedrepstr2[(i + (j * noindivs))] + 
            spryn2[(i + (j * noindivs))];
      }
      if (nonobssub3 >= 0) {
        livcheck3 = addedsize3[(i + (j * noindivs))] + addedrepstr3[(i + (j * noindivs))] + 
            spryn3[(i + (j * noindivs))] + nonobssub3;
      } else {
        livcheck3 = addedsize3[(i + (j * noindivs))] + addedrepstr3[(i + (j * noindivs))] + 
            spryn3[(i + (j * noindivs))];
      }
      
      if (livcheck1 > 0) alive1[(i + (j * noindivs))] = 1.0;
      if (livcheck2 > 0) alive2[(i + (j * noindivs))] = 1.0;
      if (livcheck3 > 0) alive3[(i + (j * noindivs))] = 1.0;
      
      if (alivegiven1[(i + (j * noindivs))] > 0) alive1[(i + (j * noindivs))] = 1.0;
      if (alivegiven2[(i + (j * noindivs))] > 0) alive2[(i + (j * noindivs))] = 1.0;
      if (alivegiven3[(i + (j * noindivs))] > 0) alive3[(i + (j * noindivs))] = 1.0;
      
      if (deadgiven1[(i + (j * noindivs))] > 0) alive1[(i + (j * noindivs))] = 0.0;
      if (deadgiven2[(i + (j * noindivs))] > 0) alive2[(i + (j * noindivs))] = 0.0;
      if (deadgiven3[(i + (j * noindivs))] > 0) alive3[(i + (j * noindivs))] = 0.0;
      
      // Corrections to living status, lifespan, and age based on later sightings
      if (firstseenx[i] <= (j + firstyear) && lastseenx[i] >= (j + firstyear)) {
        firstseen[(i + (j * noindivs))] = firstseenx[i];
        lastseen[(i + (j * noindivs))] = lastseenx[i];
        
        obsage[(i + (j * noindivs))] = year2[(i + (j * noindivs))] - firstseen[(i + (j * noindivs))];
        obslifespan[(i + (j * noindivs))] = lastseen[(i + (j * noindivs))] - firstseen[(i + (j * noindivs))];
        
        alive2[(i + (j * noindivs))] = 1;
        
        if (lastseenx[i] >= (j + firstyear + 1)) {
          alive3[(i + (j * noindivs))] = 1.0;
        }
        if (firstseenx[i] <= (j + firstyear - 1)) {
          alive1[(i + (j * noindivs))] = 1.0;
        }
      }
      
      // Stage assignments
      if (stassign && stagecol(0) == -1) {
        if (stszcol == 8) {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          stagesize1b = sizeb10[(i + (j * noindivs))];
          stagesize2b = sizeb20[(i + (j * noindivs))];
          stagesize3b = sizeb30[(i + (j * noindivs))];
          stagesize1c = sizec10[(i + (j * noindivs))];
          stagesize2c = sizec20[(i + (j * noindivs))];
          stagesize3c = sizec30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini1c = find(insfszminarmac < stagesize1c);
          arma::uvec stagemaxi1c = find(insfszmaxarmac >= stagesize1c);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini2c = find(insfszminarmac < stagesize2c);
          arma::uvec stagemaxi2c = find(insfszmaxarmac >= stagesize2c);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          arma::uvec stagemini3c = find(insfszminarmac < stagesize3c);
          arma::uvec stagemaxi3c = find(insfszmaxarmac >= stagesize3c);
          
          arma::uvec stagemini1d = intersect(stagemini1a, stagemini1b);
          stagemini1 = intersect(stagemini1c, stagemini1d);
          arma::uvec stagemaxi1d = intersect(stagemaxi1a, stagemaxi1b);
          stagemaxi1 = intersect(stagemaxi1c, stagemaxi1d);
          arma::uvec stagemini2d = intersect(stagemini2a, stagemini2b);
          stagemini2 = intersect(stagemini2c, stagemini2d);
          arma::uvec stagemaxi2d = intersect(stagemaxi2a, stagemaxi2b);
          stagemaxi2 = intersect(stagemaxi2c, stagemaxi2d);
          arma::uvec stagemini3d = intersect(stagemini3a, stagemini3b);
          stagemini3 = intersect(stagemini3c, stagemini3d);
          arma::uvec stagemaxi3d = intersect(stagemaxi3a, stagemaxi3b);
          stagemaxi3 = intersect(stagemaxi3c, stagemaxi3d);
          
        } else if (stszcol == 7) {
          stagesize1 = sizeb10[(i + (j * noindivs))];
          stagesize2 = sizeb20[(i + (j * noindivs))];
          stagesize3 = sizeb30[(i + (j * noindivs))];
          stagesize1b = sizec10[(i + (j * noindivs))];
          stagesize2b = sizec20[(i + (j * noindivs))];
          stagesize3b = sizec30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 6) {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          stagesize1b = sizec10[(i + (j * noindivs))];
          stagesize2b = sizec20[(i + (j * noindivs))];
          stagesize3b = sizec30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 5) {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          stagesize1b = sizeb10[(i + (j * noindivs))];
          stagesize2b = sizeb20[(i + (j * noindivs))];
          stagesize3b = sizeb30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 4) {
          stagesize1 = addedsize1[(i + (j * noindivs))];
          stagesize2 = addedsize2[(i + (j * noindivs))];
          stagesize3 = addedsize3[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else if (stszcol == 3) {
          stagesize1 = sizec10[(i + (j * noindivs))];
          stagesize2 = sizec20[(i + (j * noindivs))];
          stagesize3 = sizec30[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else if (stszcol == 2) {
          stagesize1 = sizeb10[(i + (j * noindivs))];
          stagesize2 = sizeb20[(i + (j * noindivs))];
          stagesize3 = sizeb30[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
        }
        
        // Stage 2
        stagerep = find(inrepstatarma == repyn2[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat2[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn2[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini2, stagemaxi2);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive2[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage2num[(i + (j * noindivs))] = static_cast<double>(choicestage) + 1.0;
          
          stage2[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage2[(i + (j * noindivs))] = "NotAlive";
        }
        
        // Stage 1
        stagerep = find(inrepstatarma == repyn1[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat1[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn1[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini1, stagemaxi1);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive1[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage1num[(i + (j * noindivs))] = static_cast<double>(choicestage) + 1.0;
          
          stage1[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage1[(i + (j * noindivs))] = "NotAlive";
          matstat1[(i + (j * noindivs))] = 0.0;
        }
        
        // Stage 3
        stagerep = find(inrepstatarma == repyn3[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat3[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn3[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini3, stagemaxi3);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        // Exceptions based on stage assignment problems in time t+1
        if (cs4.n_elem == 1 && alive3[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage3num[(i + (j * noindivs))] = static_cast<double>(choicestage) + 1.0;
          
          stage3[(i + (j * noindivs))] = sfname[choicestage];
          
        } else if (alive3[(i + (j * noindivs))] != 1) {
          stage3[(i + (j * noindivs))] = "NotAlive";
          
        } else if (cs4.n_elem == 0) {
          stage3[(i + (j * noindivs))] = "NoMatch";
          
          if (!quiet) Rf_warningcall(R_NilValue, 
              "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
          
        } else if (cs4.n_elem > 1) {
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages in stageframe have same description. All stages should be unique.");
          
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
      } // stassign if statement
    } // i loop
  } // j loop
  
  if (NAas0) {
    sizea1 = sizea10;
    sizea2 = sizea20;
    sizea3 = sizea30;
    
    sizeb1 = sizeb10;
    sizeb2 = sizeb20;
    sizeb3 = sizeb30;
    
    sizec1 = sizec10;
    sizec2 = sizec20;
    sizec3 = sizec30;
    
    repstra1 = repstra10;
    repstra2 = repstra20;
    repstra3 = repstra30;
    
    repstrb1 = repstrb10;
    repstrb2 = repstrb20;
    repstrb3 = repstrb30;
    
    feca1 = feca10;
    feca2 = feca20;
    feca3 = feca30;
    
    fecb1 = fecb10;
    fecb2 = fecb20;
    fecb3 = fecb30;
  }
  
  // Final output
  List final_output;
  
  if (reduce) { 
    bool xpos1_used {false}, ypos1_used {false}, xpos2_used {false};
    bool ypos2_used {false}, xpos3_used {false}, ypos3_used {false};
    bool censor1_used {false}, censor2_used {false}, censor3_used {false};
    bool sizea1_used {false}, sizea2_used {false}, sizea3_used {false};
    bool sizeb1_used {false}, sizeb2_used {false}, sizeb3_used {false};
    bool sizec1_used {false}, sizec2_used {false}, sizec3_used {false};
    bool size1added_used {false}, size2added_used {false}, size3added_used {false};
    bool repstra1_used {false}, repstra2_used {false}, repstra3_used {false};
    bool repstrb1_used {false}, repstrb2_used {false}, repstrb3_used {false};
    bool repstr1added_used {false}, repstr2added_used {false}, repstr3added_used {false};
    bool feca1_used {false}, feca2_used {false}, feca3_used {false};
    bool fecb1_used {false}, fecb2_used {false}, fecb3_used {false};
    bool fec1added_used {false}, fec2added_used {false}, fec3added_used {false};
    bool indcova1_used {false}, indcova2_used {false}, indcova3_used {false};
    bool indcovb1_used {false}, indcovb2_used {false}, indcovb3_used {false};
    bool indcovc1_used {false}, indcovc2_used {false}, indcovc3_used {false};
    bool juvgiven1_used {false}, juvgiven2_used {false}, juvgiven3_used {false};
    
    NumericVector xpos1_u = unique(xpos1);
    if (xpos1_u.length() > 1) xpos1_used = true;
    NumericVector ypos1_u = unique(ypos1);
    if (ypos1_u.length() > 1) ypos1_used = true;
    
    NumericVector xpos2_u = unique(xpos2);
    if (xpos2_u.length() > 1) xpos2_used = true;
    NumericVector ypos2_u = unique(ypos2);
    if (ypos2_u.length() > 1) ypos2_used = true;
    
    NumericVector xpos3_u = unique(xpos3);
    if (xpos3_u.length() > 1) xpos3_used = true;
    NumericVector ypos3_u = unique(ypos3);
    if (ypos3_u.length() > 1) ypos3_used = true;
    
    NumericVector censor1_u = unique(censor1);
    if (censor1_u.length() > 1) censor1_used = true;
    NumericVector censor2_u = unique(censor2);
    if (censor2_u.length() > 1) censor2_used = true;
    NumericVector censor3_u = unique(censor3);
    if (censor3_u.length() > 1) censor3_used = true;
    
    NumericVector sizea1_u = unique(sizea1);
    if (sizea1_u.length() > 1) sizea1_used = true;
    NumericVector sizea2_u = unique(sizea2);
    if (sizea2_u.length() > 1) sizea2_used = true;
    NumericVector sizea3_u = unique(sizea3);
    if (sizea3_u.length() > 1) sizea3_used = true;
    
    NumericVector sizeb1_u = unique(sizeb1);
    if (sizeb1_u.length() > 1) sizeb1_used = true;
    NumericVector sizeb2_u = unique(sizeb2);
    if (sizeb2_u.length() > 1) sizeb2_used = true;
    NumericVector sizeb3_u = unique(sizeb3);
    if (sizeb3_u.length() > 1) sizeb3_used = true;
    
    NumericVector sizec1_u = unique(sizec1);
    if (sizec1_u.length() > 1) sizec1_used = true;
    NumericVector sizec2_u = unique(sizec2);
    if (sizec2_u.length() > 1) sizec2_used = true;
    NumericVector sizec3_u = unique(sizec3);
    if (sizec3_u.length() > 1) sizec3_used = true;
    
    NumericVector size1added_u = unique(addedsize1);
    if (size1added_u.length() > 1) size1added_used = true;
    NumericVector size2added_u = unique(addedsize2);
    if (size2added_u.length() > 1) size2added_used = true;
    NumericVector size3added_u = unique(addedsize3);
    if (size3added_u.length() > 1) size3added_used = true;
    
    NumericVector repstra1_u = unique(repstra1);
    if (repstra1_u.length() > 1) repstra1_used = true;
    NumericVector repstra2_u = unique(repstra2);
    if (repstra2_u.length() > 1) repstra2_used = true;
    NumericVector repstra3_u = unique(repstra3);
    if (repstra3_u.length() > 1) repstra3_used = true;
    
    NumericVector repstrb1_u = unique(repstrb1);
    if (repstrb1_u.length() > 1) repstrb1_used = true;
    NumericVector repstrb2_u = unique(repstrb2);
    if (repstrb2_u.length() > 1) repstrb2_used = true;
    NumericVector repstrb3_u = unique(repstrb3);
    if (repstrb3_u.length() > 1) repstrb3_used = true;
    
    NumericVector repstr1added_u = unique(addedrepstr1);
    if (repstr1added_u.length() > 1) repstr1added_used = true;
    NumericVector repstr2added_u = unique(addedrepstr2);
    if (repstr2added_u.length() > 1) repstr2added_used = true;
    NumericVector repstr3added_u = unique(addedrepstr3);
    if (repstr3added_u.length() > 1) repstr3added_used = true;
    
    NumericVector feca1_u = unique(feca1);
    if (feca1_u.length() > 1) feca1_used = true;
    NumericVector feca2_u = unique(feca2);
    if (feca2_u.length() > 1) feca2_used = true;
    NumericVector feca3_u = unique(feca3);
    if (feca3_u.length() > 1) feca3_used = true;
    
    NumericVector fecb1_u = unique(fecb1);
    if (fecb1_u.length() > 1) fecb1_used = true;
    NumericVector fecb2_u = unique(fecb2);
    if (fecb2_u.length() > 1) fecb2_used = true;
    NumericVector fecb3_u = unique(fecb3);
    if (fecb3_u.length() > 1) fecb3_used = true;
    
    NumericVector fec1added_u = unique(addedfec1);
    if (fec1added_u.length() > 1) fec1added_used = true;
    NumericVector fec2added_u = unique(addedfec2);
    if (fec2added_u.length() > 1) fec2added_used = true;
    NumericVector fec3added_u = unique(addedfec3);
    if (fec3added_u.length() > 1) fec3added_used = true;
    
    NumericVector indcova1_u = unique(indcova1);
    IntegerVector indcova1_int_u = unique(indcova1_int);
    if (indcova1_u.length() > 1 || indcova1_int_u.length() > 1) indcova1_used = true;
    NumericVector indcova2_u = unique(indcova2);
    IntegerVector indcova2_int_u = unique(indcova2_int);
    if (indcova2_u.length() > 1 || indcova2_int_u.length() > 1) indcova2_used = true;
    NumericVector indcova3_u = unique(indcova3);
    IntegerVector indcova3_int_u = unique(indcova3_int);
    if (indcova3_u.length() > 1 || indcova3_int_u.length() > 1) indcova3_used = true;
    
    NumericVector indcovb1_u = unique(indcovb1);
    IntegerVector indcovb1_int_u = unique(indcovb1_int);
    if (indcovb1_u.length() > 1 || indcovb1_int_u.length() > 1) indcovb1_used = true;
    NumericVector indcovb2_u = unique(indcovb2);
    IntegerVector indcovb2_int_u = unique(indcovb2_int);
    if (indcovb2_u.length() > 1 || indcovb2_int_u.length() > 1) indcovb2_used = true;
    NumericVector indcovb3_u = unique(indcovb3);
    IntegerVector indcovb3_int_u = unique(indcovb3_int);
    if (indcovb3_u.length() > 1 || indcovb3_int_u.length() > 1) indcovb3_used = true;
    
    NumericVector indcovc1_u = unique(indcovc1);
    IntegerVector indcovc1_int_u = unique(indcovc1_int);
    if (indcovc1_u.length() > 1 || indcovc1_int_u.length() > 1) indcovc1_used = true;
    NumericVector indcovc2_u = unique(indcovc2);
    IntegerVector indcovc2_int_u = unique(indcovc2_int);
    if (indcovc2_u.length() > 1 || indcovc2_int_u.length() > 1) indcovc2_used = true;
    NumericVector indcovc3_u = unique(indcovc3);
    IntegerVector indcovc3_int_u = unique(indcovc3_int);
    if (indcovc3_u.length() > 1 || indcovc3_int_u.length() > 1) indcovc3_used = true;
    
    NumericVector juvgiven1_u = unique(juvgiven1);
    if (juvgiven1_u.length() > 1) juvgiven1_used = true;
    NumericVector juvgiven2_u = unique(juvgiven2);
    if (juvgiven2_u.length() > 1) juvgiven2_used = true;
    NumericVector juvgiven3_u = unique(juvgiven3);
    if (juvgiven3_u.length() > 1) juvgiven3_used = true;
    
    int total_added = xpos1_used + ypos1_used + xpos2_used + ypos2_used +
      xpos3_used + ypos3_used + censor1_used + censor2_used + censor3_used +
      sizea1_used + sizea2_used + sizea3_used + sizeb1_used + sizeb2_used +
      sizeb3_used + sizec1_used + sizec2_used + sizec3_used + size1added_used +
      size2added_used + size3added_used + repstra1_used + repstra2_used +
      repstra3_used + repstrb1_used + repstrb2_used + repstrb3_used +
      repstr1added_used + repstr2added_used + repstr3added_used + feca1_used +
      feca2_used + feca3_used + fecb1_used + fecb2_used + fecb3_used +
      fec1added_used + fec2added_used + fec3added_used + indcova1_used +
      indcova2_used + indcova3_used + indcovb1_used + indcovb2_used +
      indcovb3_used + indcovc1_used + indcovc2_used + indcovc3_used +
      juvgiven1_used + juvgiven2_used + juvgiven3_used;
    
    int list_length = 30 + total_added;
    List reduced (list_length);
    Rcpp::CharacterVector varnames (list_length);
    
    reduced(0) = rowid;
    varnames(0) = "rowid";
    
    if (popid_type > 0) {
      popid_int.attr("class") = popid_class;
      if (popid_type > 1) {
        popid_int.attr("levels") = popid_levels;
      }
      reduced(1) = popid_int;
      
    } else { 
      reduced(1) = popid;
    }
    varnames(1) = "popid";
    
    if (patchid_type > 0) {
      patchid_int.attr("class") = patchid_class;
      if (patchid_type > 1) {
        patchid_int.attr("levels") = patchid_levels;
      }
      reduced(2) = patchid_int;
      
    } else { 
      reduced(2) = patchid;
    }
    varnames(2) = "patchid";
    
    if (individ_type > 0) {
      individ_int.attr("class") = individ_class;
      if (individ_type > 1) {
        individ_int.attr("levels") = individ_levels;
      }
      reduced(3) = individ_int;
      
    } else { 
      reduced(3) = individ;
    }
    varnames(3) = "individ";
    
    reduced(4) = year2;
    varnames(4) = "year2";
    reduced(5) = firstseen;
    varnames(5) = "firstseen";
    reduced(6) = lastseen;
    varnames(6) = "lastseen";
    reduced(7) = obsage;
    varnames(7) = "obsage";
    reduced(8) = obslifespan;
    varnames(8) = "obslifespan";
    
    int red_tracker {9};
    
    if (xpos1_used) {
      reduced(red_tracker) = xpos1;
      varnames(red_tracker) = "xpos1";
      red_tracker++;
    }
    if (ypos1_used) {
      reduced(red_tracker) = ypos1;
      varnames(red_tracker) = "ypos1";
      red_tracker++;
    }
    if (sizea1_used) {
      reduced(red_tracker) = sizea1;
      varnames(red_tracker) = "sizea1";
      red_tracker++;
    }
    if (sizeb1_used) {
      reduced(red_tracker) = sizeb1;
      varnames(red_tracker) = "sizeb1";
      red_tracker++;
    }
    if (sizec1_used) {
      reduced(red_tracker) = sizec1;
      varnames(red_tracker) = "sizec1";
      red_tracker++;
    }
    if (size1added_used) {
      reduced(red_tracker) = addedsize1;
      varnames(red_tracker) = "size1added";
      red_tracker++;
    }
    if (repstra1_used) {
      reduced(red_tracker) = repstra1;
      varnames(red_tracker) = "repstra1";
      red_tracker++;
    }
    if (repstrb1_used) {
      reduced(red_tracker) = repstrb1;
      varnames(red_tracker) = "repstrb1";
      red_tracker++;
    }
    if (repstr1added_used) {
      reduced(red_tracker) = addedrepstr1;
      varnames(red_tracker) = "repstr1added";
      red_tracker++;
    }
    if (feca1_used) {
      reduced(red_tracker) = feca1;
      varnames(red_tracker) = "feca1";
      red_tracker++;
    }
    if (fecb1_used) {
      reduced(red_tracker) = fecb1;
      varnames(red_tracker) = "fecb1";
      red_tracker++;
    }
    if (fec1added_used) {
      reduced(red_tracker) = addedfec1;
      varnames(red_tracker) = "fec1added";
      red_tracker++;
    }
    
    if (indcova1_used) {
      if (!indcova_as_int) {
        reduced(red_tracker) = indcova1;
      } else {
        reduced(red_tracker) = indcova1_int;
      }
      varnames(red_tracker) = "indcova1";
      red_tracker++;
    }
    
    if (indcovb1_used) {
      if (!indcovb_as_int) {
        reduced(red_tracker) = indcovb1;
      } else {
        reduced(red_tracker) = indcovb1_int;
      }
      varnames(red_tracker) = "indcovb1";
      red_tracker++;
    }
    
    if (indcovc1_used) {
      if (!indcovc_as_int) {
        reduced(red_tracker) = indcovc1;
      } else {
        reduced(red_tracker) = indcovc1_int;
      }
      varnames(red_tracker) = "indcovc1";
      red_tracker++;
    }
    
    if (censor1_used) {
      reduced(red_tracker) = censor1;
      varnames(red_tracker) = "censor1";
      red_tracker++;
    }
    if (juvgiven1_used) {
      reduced(red_tracker) = juvgiven1;
      varnames(red_tracker) = "juvgiven1";
      red_tracker++;
    }
    
    reduced(red_tracker) = spryn1;
    varnames(red_tracker) = "obsstatus1";
    red_tracker++;
    reduced(red_tracker) = repyn1;
    varnames(red_tracker) = "repstatus1";
      
    red_tracker++;
    reduced(red_tracker) = fecyn1;
    varnames(red_tracker) = "fecstatus1";
    red_tracker++;
    reduced(red_tracker) = matstat1;
    varnames(red_tracker) = "matstatus1";
    red_tracker++;
    reduced(red_tracker) = alive1;
    varnames(red_tracker) = "alive1";
    red_tracker++;
    reduced(red_tracker) = stage1;
    varnames(red_tracker) = "stage1";
    red_tracker++;
    reduced(red_tracker) = stage1num;
    varnames(red_tracker) = "stage1index";
    red_tracker++;
    
    if (xpos2_used) {
      reduced(red_tracker) = xpos2;
      varnames(red_tracker) = "xpos2";
      red_tracker++;
    }
    if (ypos2_used) {
      reduced(red_tracker) = ypos2;
      varnames(red_tracker) = "ypos2";
      red_tracker++;
    }
    if (sizea2_used) {
      reduced(red_tracker) = sizea2;
      varnames(red_tracker) = "sizea2";
      red_tracker++;
    }
    if (sizeb2_used) {
      reduced(red_tracker) = sizeb2;
      varnames(red_tracker) = "sizeb2";
      red_tracker++;
    }
    if (sizec2_used) {
      reduced(red_tracker) = sizec2;
      varnames(red_tracker) = "sizec2";
      red_tracker++;
    }
    if (size2added_used) {
      reduced(red_tracker) = addedsize2;
      varnames(red_tracker) = "size2added";
      red_tracker++;
    }
    if (repstra2_used) {
      reduced(red_tracker) = repstra2;
      varnames(red_tracker) = "repstra2";
      red_tracker++;
    }
    if (repstrb2_used) {
      reduced(red_tracker) = repstrb2;
      varnames(red_tracker) = "repstrb2";
      red_tracker++;
    }
    if (repstr2added_used) {
      reduced(red_tracker) = addedrepstr2;
      varnames(red_tracker) = "repstr2added";
      red_tracker++;
    }
    if (feca2_used) {
      reduced(red_tracker) = feca2;
      varnames(red_tracker) = "feca2";
      red_tracker++;
    }
    if (fecb2_used) {
      reduced(red_tracker) = fecb2;
      varnames(red_tracker) = "fecb2";
      red_tracker++;
    }
    if (fec2added_used) {
      reduced(red_tracker) = addedfec2;
      varnames(red_tracker) = "fec2added";
      red_tracker++;
    }
    
    if (indcova2_used) {
      if (!indcova_as_int) {
        reduced(red_tracker) = indcova2;
      } else {
        reduced(red_tracker) = indcova2_int;
      }
      varnames(red_tracker) = "indcova2";
      red_tracker++;
    }
    
    if (indcovb2_used) {
      if (!indcovb_as_int) {
        reduced(red_tracker) = indcovb2;
      } else {
        reduced(red_tracker) = indcovb2_int;
      }
      varnames(red_tracker) = "indcovb2";
      red_tracker++;
    }
    
    if (indcovc2_used) {
      if (!indcovc_as_int) {
        reduced(red_tracker) = indcovc2;
      } else {
        reduced(red_tracker) = indcovc2_int;
      }
      varnames(red_tracker) = "indcovc2";
      red_tracker++;
    }
    
    if (censor2_used) {
      reduced(red_tracker) = censor2;
      varnames(red_tracker) = "censor2";
      red_tracker++;
    }
    if (juvgiven2_used) {
      reduced(red_tracker) = juvgiven2;
      varnames(red_tracker) = "juvgiven2";
      red_tracker++;
    }
    
    reduced(red_tracker) = spryn2;
    varnames(red_tracker) = "obsstatus2";
    red_tracker++;
    reduced(red_tracker) = repyn2;
    varnames(red_tracker) = "repstatus2";
      
    red_tracker++;
    reduced(red_tracker) = fecyn2;
    varnames(red_tracker) = "fecstatus2";
    red_tracker++;
    reduced(red_tracker) = matstat2;
    varnames(red_tracker) = "matstatus2";
    red_tracker++;
    reduced(red_tracker) = alive2;
    varnames(red_tracker) = "alive2";
    red_tracker++;
    reduced(red_tracker) = stage2;
    varnames(red_tracker) = "stage2";
    red_tracker++;
    reduced(red_tracker) = stage2num;
    varnames(red_tracker) = "stage2index";
    red_tracker++;
    
    if (xpos3_used) {
      reduced(red_tracker) = xpos3;
      varnames(red_tracker) = "xpos3";
      red_tracker++;
    }
    if (ypos3_used) {
      reduced(red_tracker) = ypos3;
      varnames(red_tracker) = "ypos3";
      red_tracker++;
    }
    if (sizea3_used) {
      reduced(red_tracker) = sizea3;
      varnames(red_tracker) = "sizea3";
      red_tracker++;
    }
    if (sizeb3_used) {
      reduced(red_tracker) = sizeb3;
      varnames(red_tracker) = "sizeb3";
      red_tracker++;
    }
    if (sizec3_used) {
      reduced(red_tracker) = sizec3;
      varnames(red_tracker) = "sizec3";
      red_tracker++;
    }
    if (size3added_used) {
      reduced(red_tracker) = addedsize3;
      varnames(red_tracker) = "size3added";
      red_tracker++;
    }
    if (repstra3_used) {
      reduced(red_tracker) = repstra3;
      varnames(red_tracker) = "repstra3";
      red_tracker++;
    }
    if (repstrb3_used) {
      reduced(red_tracker) = repstrb3;
      varnames(red_tracker) = "repstrb3";
      red_tracker++;
    }
    if (repstr3added_used) {
      reduced(red_tracker) = addedrepstr3;
      varnames(red_tracker) = "repstr3added";
      red_tracker++;
    }
    if (feca3_used) {
      reduced(red_tracker) = feca3;
      varnames(red_tracker) = "feca3";
      red_tracker++;
    }
    if (fecb3_used) {
      reduced(red_tracker) = fecb3;
      varnames(red_tracker) = "fecb3";
      red_tracker++;
    }
    if (fec3added_used) {
      reduced(red_tracker) = addedfec3;
      varnames(red_tracker) = "fec3added";
      red_tracker++;
    }
    
    if (indcova3_used) {
      if (!indcova_as_int) {
        reduced(red_tracker) = indcova3;
      } else {
        reduced(red_tracker) = indcova3_int;
      }
      varnames(red_tracker) = "indcova3";
      red_tracker++;
    }
    
    if (indcovb3_used) {
      if (!indcovb_as_int) {
        reduced(red_tracker) = indcovb3;
      } else {
        reduced(red_tracker) = indcovb3_int;
      }
      varnames(red_tracker) = "indcovb3";
      red_tracker++;
    }
    
    if (indcovc3_used) {
      if (!indcovc_as_int) {
        reduced(red_tracker) = indcovc3;
      } else {
        reduced(red_tracker) = indcovc3_int;
      }
      varnames(red_tracker) = "indcovc3";
      red_tracker++;
    }
    
    if (censor3_used) {
      reduced(red_tracker) = censor3;
      varnames(red_tracker) = "censor3";
      red_tracker++;
    }
    if (juvgiven3_used) {
      reduced(red_tracker) = juvgiven3;
      varnames(red_tracker) = "juvgiven3";
      red_tracker++;
    }
    
    reduced(red_tracker) = spryn3;
    varnames(red_tracker) = "obsstatus3";
    red_tracker++;
    reduced(red_tracker) = repyn3;
    varnames(red_tracker) = "repstatus3";
      
    red_tracker++;
    reduced(red_tracker) = fecyn3;
    varnames(red_tracker) = "fecstatus3";
    red_tracker++;
    reduced(red_tracker) = matstat3;
    varnames(red_tracker) = "matstatus3";
    red_tracker++;
    reduced(red_tracker) = alive3;
    varnames(red_tracker) = "alive3";
    red_tracker++;
    reduced(red_tracker) = stage3;
    varnames(red_tracker) = "stage3";
    red_tracker++;
    reduced(red_tracker) = stage3num;
    varnames(red_tracker) = "stage3index";
    
    reduced.attr("names") = varnames;
    reduced.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, alive3.length());
    StringVector needed_classes {"data.frame", "hfvdata"};
    reduced.attr("class") = needed_classes;
    
    final_output = reduced;
    
  } else {
    List output_longlist (81);
    
    output_longlist(0) = rowid;
    output_longlist(4) = year2;
    output_longlist(5) = firstseen;
    output_longlist(6) = lastseen;
    output_longlist(7) = obsage;
    output_longlist(8) = obslifespan;
    
    output_longlist(9) = xpos1;
    output_longlist(10) = ypos1;
    output_longlist(11) = sizea1;
    output_longlist(12) = sizeb1;
    output_longlist(13) = sizec1;
    output_longlist(14) = addedsize1;
    output_longlist(15) = repstra1;
    output_longlist(16) = repstrb1;
    output_longlist(17) = addedrepstr1;
    output_longlist(18) = feca1;
    output_longlist(19) = fecb1;
    output_longlist(20) = addedfec1;
    
    output_longlist(24) = censor1;
    output_longlist(25) = juvgiven1;
    
    output_longlist(26) = spryn1;
    output_longlist(27) = repyn1;
    output_longlist(28) = fecyn1;
    output_longlist(29) = matstat1;
    output_longlist(30) = alive1;
    output_longlist(31) = stage1;
    output_longlist(32) = stage1num;
    
    output_longlist(33) = xpos2;
    output_longlist(34) = ypos2;
    output_longlist(35) = sizea2;
    output_longlist(36) = sizeb2;
    output_longlist(37) = sizec2;
    output_longlist(38) = addedsize2;
    output_longlist(39) = repstra2;
    output_longlist(40) = repstrb2;
    output_longlist(41) = addedrepstr2;
    output_longlist(42) = feca2;
    output_longlist(43) = fecb2;
    output_longlist(44) = addedfec2;
    
    output_longlist(48) = censor2;
    output_longlist(49) = juvgiven2;
    
    output_longlist(50) = spryn2;
    output_longlist(51) = repyn2;
    output_longlist(52) = fecyn2;
    output_longlist(53) = matstat2;
    output_longlist(54) = alive2;
    output_longlist(55) = stage2;
    output_longlist(56) = stage2num;
  
    output_longlist(57) = xpos3;
    output_longlist(58) = ypos3;
    output_longlist(59) = sizea3;
    output_longlist(60) = sizeb3;
    output_longlist(61) = sizec3;
    output_longlist(62) = addedsize3;
    output_longlist(63) = repstra3;
    output_longlist(64) = repstrb3;
    output_longlist(65) = addedrepstr3;
    output_longlist(66) = feca3;
    output_longlist(67) = fecb3;
    output_longlist(68) = addedfec3;
    
    output_longlist(72) = censor3;
    output_longlist(73) = juvgiven3;
    
    output_longlist(74) = spryn3;
    output_longlist(75) = repyn3;
    output_longlist(76) = fecyn3;
    output_longlist(77) = matstat3;
    output_longlist(78) = alive3;
    output_longlist(79) = stage3;
    output_longlist(80) = stage3num;
    
    if (popid_type > 0) {
      if (popid_type > 1) {
        popid_int.attr("class") = popid_class;
        popid_int.attr("levels") = popid_levels;
      }
      output_longlist(1) = popid_int;
      
    } else { 
      output_longlist(1) = popid;
    }
    
    if (patchid_type > 0) {
      if (patchid_type > 1) {
        patchid_int.attr("class") = patchid_class;
        patchid_int.attr("levels") = patchid_levels;
      }
      output_longlist(2) = patchid_int;
      
    } else { 
      output_longlist(2) = patchid;
    }

    if (individ_type > 0) {
      if (individ_type > 1) {
        individ_int.attr("class") = individ_class;
        individ_int.attr("levels") = individ_levels;
      }
      output_longlist(3) = individ_int;
      
    } else { 
      output_longlist(3) = individ;
    }
    
    if (!indcova_as_int) {
      output_longlist(21) = indcova1;
      output_longlist(45) = indcova2;
      output_longlist(69) = indcova3;
    } else {
      output_longlist(21) = indcova1_int;
      output_longlist(45) = indcova2_int;
      output_longlist(69) = indcova3_int;
    }
    
    if (!indcovb_as_int) {
      output_longlist(22) = indcovb1;
      output_longlist(46) = indcovb2;
      output_longlist(70) = indcovb3;
    } else {
      output_longlist(22) = indcovb1_int;
      output_longlist(46) = indcovb2_int;
      output_longlist(70) = indcovb3_int;
    }
    
    if (!indcovc_as_int) {
      output_longlist(23) = indcovc1;
      output_longlist(47) = indcovc2;
      output_longlist(71) = indcovc3;
    } else {
      output_longlist(23) = indcovc1_int;
      output_longlist(47) = indcovc2_int;
      output_longlist(71) = indcovc3_int;
    }
    
    Rcpp::CharacterVector varnames {"rowid", "popid", "patchid", "individ",
      "year2", "firstseen", "lastseen", "obsage", "obslifespan","xpos1", "ypos1",
      "sizea1", "sizeb1", "sizec1", "size1added", "repstra1", "repstrb1",
      "repstr1added", "feca1", "fecb1", "fec1added", "indcova1", "indcovb1",
      "indcovc1", "censor1", "juvgiven1", "obsstatus1", "repstatus1",
      "fecstatus1", "matstatus1", "alive1", "stage1", "stage1index",
      "xpos2", "ypos2", "sizea2", "sizeb2", "sizec2", "size2added", "repstra2",
      "repstrb2", "repstr2added", "feca2", "fecb2", "fec2added", "indcova2",
      "indcovb2", "indcovc2", "censor2", "juvgiven2", "obsstatus2", "repstatus2",
      "fecstatus2", "matstatus2", "alive2", "stage2", "stage2index",
      "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "size3added", "repstra3",
      "repstrb3", "repstr3added", "feca3", "fecb3", "fec3added", "indcova3",
      "indcovb3", "indcovc3", "censor3", "juvgiven3", "obsstatus3", "repstatus3",
      "fecstatus3", "matstatus3", "alive3", "stage3", "stage3index"};
    
    output_longlist.attr("names") = varnames;
    output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, alive3.length());
    StringVector needed_classes {"data.frame", "hfvdata"};
    output_longlist.attr("class") = needed_classes;
    
    final_output = output_longlist;
  }
  
  if (!retain_alive0) { 
    IntegerVector retained_bit {1};
    StringVector alive_var_used = {"alive2"};
    List reduced_output = LefkoUtils::df_subset(final_output, retained_bit, false, true,
      false, false, true, as<RObject>(alive_var_used));
      
    return reduced_output;
    
  } else { 
    return final_output;
  }
}

//' Create Historical Vertical Structure for Ahistorical Vertical Data Frame
//' 
//' Function \code{jpf()} is the core kernel for function
//' \code{\link{historicalize3}()}, creating the historical, vertical structure
//' and rearranging the data in that shape.
//' 
//' @name .jpf
//' 
//' @param data The horizontal data file.
//' @param stageframe The stageframe object identifying the life history model
//' being operationalized. This should be the full stageframe.
//' @param popidcol Column number corresponding to the identity of the
//' population for each individual.
//' @param patchidcol Column number corresponding to the identity of the patch
//' for each individual.
//' @param individcol Column number corresponding to the identity of each 
//' individual.
//' @param year2col Column number of year or occasion in time \emph{t}.
//' @param year3col Column number of year or occasion in time \emph{t}+1.
//' @param xcol Column number corresponding to the x coordinate of each
//' individual in Cartesian space.
//' @param ycol Column number corresponding to the y coordinate of each
//' individual in Cartesian space.
//' @param juv2col Column number coding for status as a juvenile in time
//' \emph{t}.
//' @param juv3col Column number coding for status as a juvenile in time
//' \emph{t}+1.
//' @param sizea2col Column number corresponding to the primary size variable in
//' time \emph{t}.
//' @param sizea3col Column number corresponding to the primary size variable in
//' time \emph{t}+1.
//' @param sizeb2col Column number corresponding to the secondary size variable
//' in time \emph{t}.
//' @param sizeb3col Column number corresponding to the secondary size variable
//' in time \emph{t}+1.
//' @param sizec2col Column number corresponding to the tertiary size variable
//' in time \emph{t}.
//' @param sizec3col Column number corresponding to the tertiary size variable
//' in time \emph{t}+1.
//' @param repstra2col Column number corresponding to the main variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}.
//' @param repstra3col Column number corresponding to the main variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}+1.
//' @param repstrb2col Column number corresponding to a second variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}.
//' @param repstrb3col Column number corresponding to a second variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}+1.
//' @param feca2col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}.
//' @param feca3col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}+1.
//' @param fecb2col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}.
//' @param fecb3col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}+1.
//' @param indcova2col Column number corresponding to an individual covariate in
//' time \emph{t}.
//' @param indcova3col Column number corresponding to an individual covariate in
//' time \emph{t}+1.
//' @param indcovb2col Column number corresponding to an individual covariate in
//' time \emph{t}.
//' @param indcovb3col Column number corresponding to an individual covariate in
//' time \emph{t}+1.
//' @param indcovc2col Column number corresponding to an individual covariate in
//' time \emph{t}.
//' @param indcovc3col Column number corresponding to an individual covariate in
//' time \emph{t}+1.
//' @param alive2col Column number detailing whether an individual is alive in 
//' time \emph{t}.
//' @param alive3col Column number detailing whether an individual is alive in 
//' time \emph{t}+1.
//' @param dead2col Column number detailing whether an individual is dead in 
//' time \emph{t}.
//' @param dead3col Column number detailing whether an individual is dead in 
//' time \emph{t}+1.
//' @param obs2col Column number detailing whether an individual is in an
//' observable stage in time \emph{t}.
//' @param obs3col Column number detailing whether an individual is in an
//' observable stage in time \emph{t}+1.
//' @param nonobs2col Column number detailing whether an individual is in an
//' unobservable stage in time \emph{t}.
//' @param nonobs3col Column number detailing whether an individual is in an
//' unobservable stage in time \emph{t}+1.
//' @param repstrrel This is a scalar multiplier for that makes the variable in
//' \code{repstrb2col} equivalent to \code{repstra2col}.
//' @param fecrel This is a scalar multiplier for that makes the variable in
//' \code{fecb2col} equivalent to \code{feca2col}.
//' @param stage2col Column number corresponding to life history stage in time
//' \emph{t}.
//' @param stage3col Column number corresponding to life history stage in time
//' \emph{t}+1.
//' @param censorcol Column number corresponding to a censor variable within the
//' dataset.
//' @param NAas0 If \code{TRUE}, then all \code{NA} entries for size and
//' fecundity variables will be set to 0.
//' @param NRasRep If \code{TRUE}, then non-reproductive but mature individuals
//' will be treated as reproductive during stage assignment.
//' @param NOasObs If TRUE, then will treat unobserved individuals as observed
//' during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Integer describing which size variable to use in stage 
//' estimation. Numbers 1 through 8 are possible.
//' @param censorkeep Numeric value of censor variable, denoting elements to
//' keep. If \code{NA} is to be used, then set this variable to \code{0} and set
//' \code{censbool = TRUE}.
//' @param censbool A logical variable determining whether \code{NA} denotes the
//' value of the censoring variable identifying data to keep.
//' @param retain_alive0 A logical variable indicating whether to keep or remove
//' data rows for individuals not alive in time \emph{t}.
//' @param reduce A logical variable determining whether unused variables and
//' some invariant state variables should be removed from the output dataset.
//' Defaults to \code{TRUE}.
//' @param quiet A logical value indicating whether to silense warnings.
//' 
//' @return The output is currently a list coerced into the data frame class,
//' and is read as a data frame by R. It is secondarily in class \code{hfvdata}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.jpf)]]
Rcpp::List jpf(const DataFrame& data, const DataFrame& stageframe, int popidcol,
  int patchidcol, int individcol, int year2col, int year3col, int xcol,
  int ycol, int juv2col, int juv3col, int sizea2col, int sizea3col,
  int sizeb2col, int sizeb3col, int sizec2col, int sizec3col, int repstra2col,
  int repstra3col, int repstrb2col, int repstrb3col, int feca2col, int feca3col,
  int fecb2col, int fecb3col, int indcova2col, int indcova3col, int indcovb2col,
  int indcovb3col, int indcovc2col, int indcovc3col, int alive2col,
  int alive3col, int dead2col, int dead3col, int obs2col, int obs3col,
  int nonobs2col, int nonobs3col, double repstrrel, double fecrel,
  int stage2col, int stage3col, int censorcol, bool NAas0, bool NRasRep,
  bool NOasObs, bool stassign, int stszcol, double censorkeep, bool censbool,
  bool retain_alive0, const bool reduce, bool quiet) {
  
  //Rcout << "jpf A" << endl;
  
  int norows = data.nrows();
  int noindivs {0};
  
  Rcpp::NumericVector ckcheck (1);
  ckcheck(0) = censorkeep;
  double crazycensor;
  Rcpp::NumericVector censfillvec (norows);
  if (NumericVector::is_na(ckcheck(0))) {
    crazycensor = 0;
  } else{
    crazycensor = censorkeep;
  }
  
  Rcpp::NumericVector zerovec (norows, 0.0);
  Rcpp::NumericVector negonevec (norows, -1.0);
  
  Rcpp::IntegerVector individx_int;
  Rcpp::StringVector individx = as<StringVector>(data[individcol]);
  int individ_type {0}; // 0: string (or numeric, if default), 1: integer, 2: factor, 3: numeric
  
  StringVector individ_class;
  StringVector individ_levels;
  
  RObject indiv_test = as<RObject>(data[individcol]);
  if (is<IntegerVector>(indiv_test)) {
    individ_type = 1;
    individx_int = as<IntegerVector>(data[individcol]);
    
    if (individx_int.hasAttribute("class")) {
      individ_class = as<StringVector>(individx_int.attr("class"));
    }
    
    if (individx_int.hasAttribute("levels")) {
      individ_type = 2;
      individ_levels = as<StringVector>(individx_int.attr("levels"));
    }
  }
  
  Rcpp::StringVector allindivs = unique(individx);
  noindivs = static_cast<int>(allindivs.size()); // Total no individuals
  
  Rcpp::IntegerVector year3x;
  Rcpp::IntegerVector year2x = as<IntegerVector>(data[year2col]);
  if (year3col != -1) {
    year3x = as<IntegerVector>(data[year3col]);
  } else {year3x = zerovec;}
  
  Rcpp::IntegerVector yearall2x = sort_unique(year2x);
  int firstyear = min(yearall2x);
  
  const int noyears = static_cast<int>(yearall2x.size()); // Total obs periods, minus last period in year3
  
  int ndflength = noyears * noindivs; // Initial length of final hfv dataset
  int currentyear {0};
  int currentindiv {-1};
  int ndfindex {0};
  int prevyrindex {0};
  int nextyrindex {0};
  double livcheck2 {0.0};
  double livcheck3 {0.0};
  
  Rcpp::StringVector sfname = as<StringVector>(stageframe["stage"]);
  Rcpp::NumericVector repstat = as<NumericVector>(stageframe["repstatus"]);
  Rcpp::NumericVector obsstat = as<NumericVector>(stageframe["obsstatus"]);
  Rcpp::NumericVector matstat = as<NumericVector>(stageframe["matstatus"]);
  arma::vec indataset = as<arma::vec>(stageframe["indataset"]);
  arma::vec sfszmin = as<arma::vec>(stageframe["sizebin_min"]);
  arma::vec sfszmax = as<arma::vec>(stageframe["sizebin_max"]);
  arma::vec sfszminb = as<arma::vec>(stageframe["sizebinb_min"]);
  arma::vec sfszmaxb = as<arma::vec>(stageframe["sizebinb_max"]);
  arma::vec sfszminc = as<arma::vec>(stageframe["sizebinc_min"]);
  arma::vec sfszmaxc = as<arma::vec>(stageframe["sizebinc_max"]);
  
  arma::vec repstatarma = as<arma::vec>(repstat);
  arma::vec obsstatarma = as<arma::vec>(obsstat);
  arma::vec matstatarma = as<arma::vec>(matstat);
  arma::vec sfszminarma = sfszmin;
  arma::vec sfszmaxarma = sfszmax;
  arma::vec sfszminarmab = sfszminb;
  arma::vec sfszmaxarmab = sfszmaxb;
  arma::vec sfszminarmac = sfszminc;
  arma::vec sfszmaxarmac = sfszmaxc;
  int stagenum = static_cast<int>(sfszmaxarma.n_elem); // Total no stages in stageframe
  
  arma::uvec instages = find(indataset == 1); 
  int instagenum = static_cast<int>(instages.n_elem); // Total no stages in dataset
  
  arma::uvec stageid(stagenum);
  arma::uvec instageid(instagenum);
  arma::vec insfszminarma(instagenum);
  arma::vec insfszmaxarma(instagenum);
  arma::vec inrepstatarma(instagenum);
  arma::vec inobsstatarma(instagenum);
  arma::vec inmatstatarma(instagenum);
  arma::vec insfszminarmab (instagenum);
  arma::vec insfszmaxarmab (instagenum);
  arma::vec insfszminarmac (instagenum);
  arma::vec insfszmaxarmac (instagenum);
  
  // Create vectors describing stages in dataset
  int inplace {0};
  for (int i = 0; i < stagenum; i++) {
    stageid(i) = i+1;
    
    if (indataset(i) == 1) {
      instageid(inplace) = i + 1;
      insfszminarma(inplace) = sfszminarma(i);
      insfszmaxarma(inplace) = sfszmaxarma(i);
      insfszminarmab(inplace) = sfszminarmab(i);
      insfszmaxarmab(inplace) = sfszmaxarmab(i);
      insfszminarmac(inplace) = sfszminarmac(i);
      insfszmaxarmac(inplace) = sfszmaxarmac(i);
      inrepstatarma(inplace) = repstatarma(i);
      inobsstatarma(inplace) = obsstatarma(i);
      inmatstatarma(inplace) = matstatarma(i);
      
      inplace++;
    }
  }
  
  // Variables defined by row structure of data
  Rcpp::StringVector popidx;
  Rcpp::StringVector patchidx;
  
  Rcpp::IntegerVector popidx_int;
  Rcpp::IntegerVector patchidx_int;
  
  int popid_type {0};
  int patchid_type {0};
  
  StringVector popid_class;
  StringVector patchid_class;
  
  StringVector popid_levels;
  StringVector patchid_levels;
  
  Rcpp::NumericVector xpos2x;
  Rcpp::NumericVector ypos2x;
  Rcpp::NumericVector sizea2x;
  Rcpp::NumericVector sizea3x;
  Rcpp::NumericVector repstra2x;
  Rcpp::NumericVector repstra3x;
  Rcpp::NumericVector feca2x;
  Rcpp::NumericVector feca3x;
  Rcpp::NumericVector sizeb2x;
  Rcpp::NumericVector sizeb3x;
  Rcpp::NumericVector repstrb2x;
  Rcpp::NumericVector repstrb3x;
  Rcpp::NumericVector fecb2x;
  Rcpp::NumericVector fecb3x;
  Rcpp::NumericVector sizec2x;
  Rcpp::NumericVector sizec3x;
  
  Rcpp::NumericVector indcova2x;
  Rcpp::NumericVector indcova3x;
  Rcpp::NumericVector indcovb2x;
  Rcpp::NumericVector indcovb3x;
  Rcpp::NumericVector indcovc2x;
  Rcpp::NumericVector indcovc3x;
  
  Rcpp::IntegerVector indcova2x_int;
  Rcpp::IntegerVector indcova3x_int;
  Rcpp::IntegerVector indcovb2x_int;
  Rcpp::IntegerVector indcovb3x_int;
  Rcpp::IntegerVector indcovc2x_int;
  Rcpp::IntegerVector indcovc3x_int;
  
  Rcpp::StringVector indcova2x_str;
  Rcpp::StringVector indcova3x_str;
  Rcpp::StringVector indcovb2x_str;
  Rcpp::StringVector indcovb3x_str;
  Rcpp::StringVector indcovc2x_str;
  Rcpp::StringVector indcovc3x_str;
  
  int indcova_type {3}; // 0: string, 1: integer, 2: factor, 3: numeric
  int indcovb_type {3}; // 0: string, 1: integer, 2: factor, 3: numeric
  int indcovc_type {3}; // 0: string, 1: integer, 2: factor, 3: numeric
  
  StringVector indcova_class;
  StringVector indcovb_class;
  StringVector indcovc_class;
  StringVector indcova_levels;
  StringVector indcovb_levels;
  StringVector indcovc_levels;
  
  Rcpp::NumericVector censor2x;
  Rcpp::NumericVector alivegiven2x;
  Rcpp::NumericVector alivegiven3x;
  Rcpp::NumericVector deadgiven2x;
  Rcpp::NumericVector deadgiven3x;
  Rcpp::NumericVector obsgiven2x;
  Rcpp::NumericVector obsgiven3x;
  Rcpp::NumericVector nonobsgiven2x; 
  Rcpp::NumericVector nonobsgiven3x;
  Rcpp::NumericVector juvgiven2x;
  Rcpp::NumericVector juvgiven3x;
  
  Rcpp::NumericVector firstseenx (noindivs, -1.0);
  Rcpp::NumericVector lastseenx (noindivs, -1.0);
  
  Rcpp::NumericVector sizea20x (norows);
  Rcpp::NumericVector sizeb20x (norows);
  Rcpp::NumericVector sizec20x (norows);
  Rcpp::NumericVector repstra20x (norows);
  Rcpp::NumericVector repstrb20x (norows);
  Rcpp::NumericVector feca20x (norows);
  Rcpp::NumericVector feca30x (norows);
  Rcpp::NumericVector juvgiven20x (norows);
  
  Rcpp::NumericVector sizea30x (norows);
  Rcpp::NumericVector sizeb30x (norows);
  Rcpp::NumericVector sizec30x (norows);
  Rcpp::NumericVector repstra30x (norows);
  Rcpp::NumericVector repstrb30x (norows);
  Rcpp::NumericVector fecb20x (norows);
  Rcpp::NumericVector fecb30x (norows);
  Rcpp::NumericVector juvgiven30x (norows);
  
  Rcpp::StringVector stage2x;
  Rcpp::StringVector stage3x;
  
  // Assign values from the dataset
  if (popidcol != -1) {
    RObject popid_test = as<RObject>(data[popidcol]);
    
    if (is<IntegerVector>(popid_test)) { 
      popid_type = 1;
      popidx_int = as<IntegerVector>(data[popidcol]);
      
      if (popidx_int.hasAttribute("class")) {
        popid_class = as<StringVector>(popidx_int.attr("class"));
      }
      
      if (popidx_int.hasAttribute("levels")) { 
        popid_type = 2;
        popid_levels = as<StringVector>(popidx_int.attr("levels"));
      }
    } else {
      popidx = as<StringVector>(data[popidcol]);
    }
  } else {popidx = as<StringVector>(clone(zerovec));}
  
  if (patchidcol != -1) {
    RObject patchid_test = as<RObject>(data[patchidcol]);
    
    if (is<IntegerVector>(patchid_test)) { 
      patchid_type = 1;
      patchidx_int = as<IntegerVector>(data[patchidcol]);
      
      if (patchidx_int.hasAttribute("class")) {
        patchid_class = as<StringVector>(patchidx_int.attr("class"));
      }
      
      if (patchidx_int.hasAttribute("levels")) {
        patchid_type = 2;
        patchid_levels = as<StringVector>(patchidx_int.attr("levels"));
      }
    } else {
      patchidx = as<StringVector>(data[patchidcol]);
    }
  } else {patchidx = as<StringVector>(clone(zerovec));}
  
  if (censorcol != -1) {
    censor2x = as<NumericVector>(data[censorcol]);
  } else {censor2x = clone(zerovec);}
  
  if (sizea2col != -1) {
    sizea2x = as<NumericVector>(data[sizea2col]);
  } else {sizea2x = clone(zerovec);}
  if (sizeb2col != -1) {
    sizeb2x = as<NumericVector>(data[sizeb2col]);
  } else {sizeb2x = clone(zerovec);}
  if (sizec2col != -1) {
    sizec2x = as<NumericVector>(data[sizec2col]);
  } else {sizec2x = clone(zerovec);}
  if (repstra2col != -1) {
    repstra2x = as<NumericVector>(data[repstra2col]);
  } else {repstra2x = clone(zerovec);}
  if (repstrb2col != -1) {
    repstrb2x = as<NumericVector>(data[repstrb2col]);
  } else {repstrb2x = clone(zerovec);}
  if (feca2col != -1) {
    feca2x = as<NumericVector>(data[feca2col]);
  } else {feca2x = clone(zerovec);}
  if (fecb2col != -1) {
    fecb2x = as<NumericVector>(data[fecb2col]);
  } else {fecb2x = clone(zerovec);}
  
  //Rcout << "jpf B" << endl;
  
  if (indcova2col != -1) {
    RObject indcova2_test = as<RObject>(data[indcova2col]);
    
    if (is<NumericVector>(indcova2_test)) {
      indcova_type = 3;
      indcova2x = as<NumericVector>(data[indcova2col]);
    } else if (is<IntegerVector>(indcova2_test)) {
      indcova_type = 1;
      indcova2x_int = as<IntegerVector>(data[indcova2col]);
      
      if (indcova2x_int.hasAttribute("class")) {
        indcova_class = indcova2x_int.attr("class");
      }
      
      if (indcova2x_int.hasAttribute("levels")) {
        indcova_type = 2;
        indcova_levels = indcova2x_int.attr("levels");
      }
    } else {
      indcova_type = 0;
      indcova2x_str = as<StringVector>(data[indcova2col]);
    }
  } else {
    indcova2x = clone(zerovec);
  }
  
  if (indcova3col != -1) {
    RObject indcova3_test = as<RObject>(data[indcova3col]);
    
    if (is<NumericVector>(indcova3_test)) {
      indcova_type = 3;
      indcova3x = as<NumericVector>(data[indcova3col]);
    } else if (is<IntegerVector>(indcova3_test)) {
      indcova_type = 1;
      indcova3x_int = as<IntegerVector>(data[indcova3col]);
      
      if (indcova3x_int.hasAttribute("class")) {
        indcova_class = indcova3x_int.attr("class");
      }
      
      if (indcova3x_int.hasAttribute("levels")) {
        indcova_type = 2;
        indcova_levels = indcova3x_int.attr("levels");
      }
    } else {
      indcova_type = 0;
      indcova3x_str = as<StringVector>(data[indcova3col]);
    }
  } else {
    indcova3x = clone(zerovec);
  }
  
  if (indcovb2col != -1) {
    RObject indcovb2_test = as<RObject>(data[indcovb2col]);
    
    if (is<NumericVector>(indcovb2_test)) {
      indcovb_type = 3;
      indcovb2x = as<NumericVector>(data[indcovb2col]);
    } else if (is<IntegerVector>(indcovb2_test)) {
      indcovb_type = 1;
      indcovb2x_int = as<IntegerVector>(data[indcovb2col]);
      
      if (indcovb2x_int.hasAttribute("class")) {
        indcovb_class = indcovb2x_int.attr("class");
      }
      
      if (indcovb2x_int.hasAttribute("levels")) {
        indcovb_type = 2;
        indcovb_levels = indcovb2x_int.attr("levels");
      }
    } else {
      indcovb_type = 0;
      indcovb2x_str = as<StringVector>(data[indcovb2col]);
    }
  } else {
    indcovb2x = clone(zerovec);
  }
  
  if (indcovb3col != -1) {
    RObject indcovb3_test = as<RObject>(data[indcovb3col]);
    
    if (is<NumericVector>(indcovb3_test)) {
      indcovb_type = 3;
      indcovb3x = as<NumericVector>(data[indcovb3col]);
    } else if (is<IntegerVector>(indcovb3_test)) {
      indcovb_type = 1;
      indcovb3x_int = as<IntegerVector>(data[indcovb3col]);
      
      if (indcovb3x_int.hasAttribute("class")) {
        indcovb_class = indcovb3x_int.attr("class");
      }
      
      if (indcovb3x_int.hasAttribute("levels")) {
        indcovb_type = 2;
        indcovb_levels = indcovb3x_int.attr("levels");
      }
    } else {
      indcovb_type = 0;
      indcovb3x_str = as<StringVector>(data[indcovb3col]);
    }
  } else {
    indcovb3x = clone(zerovec);
  }
  
  if (indcovc2col != -1) {
    RObject indcovc2_test = as<RObject>(data[indcovc2col]);
    
    if (is<NumericVector>(indcovc2_test)) {
      indcovc_type = 3;
      indcovc2x = as<NumericVector>(data[indcovc2col]);
    } else if (is<IntegerVector>(indcovc2_test)) {
      indcovc_type = 1;
      indcovc2x_int = as<IntegerVector>(data[indcovc2col]);
      
      if (indcovc2x_int.hasAttribute("class")) {
        indcovc_class = indcovc2x_int.attr("class");
      }
      
      if (indcovc2x_int.hasAttribute("levels")) {
        indcovc_type = 2;
        indcovc_levels = indcovc2x_int.attr("levels");
      }
    } else {
      indcovc_type = 0;
      indcovc2x_str = as<StringVector>(data[indcovc2col]);
    }
  } else {
    indcovc2x = clone(zerovec);
  }
  
  if (indcovc3col != -1) {
    RObject indcovc3_test = as<RObject>(data[indcovc3col]);
    
    if (is<NumericVector>(indcovc3_test)) {
      indcovc_type = 3;
      indcovc3x = as<NumericVector>(data[indcovc3col]);
    } else if (is<IntegerVector>(indcovc3_test)) {
      indcovc_type = 1;
      indcovc3x_int = as<IntegerVector>(data[indcovc3col]);
      
      if (indcovc3x_int.hasAttribute("class")) {
        indcovc_class = indcovc3x_int.attr("class");
      }
      
      if (indcovc3x_int.hasAttribute("levels")) {
        indcovc_type = 1;
        indcovc_levels = indcovc3x_int.attr("levels");
      }
    } else {
      indcovc_type = 0;
      indcovc3x_str = as<StringVector>(data[indcovc3col]);
    }
  } else {
    indcovc3x = clone(zerovec);
  }
  
  //Rcout << "jpf C" << endl;
  
  if (juv2col != -1) {
    juvgiven2x = as<NumericVector>(data[juv2col]);
  } else {juvgiven2x = clone(zerovec);}
  if (obs2col != -1) {
    obsgiven2x = as<NumericVector>(data[obs2col]);
  } else {obsgiven2x = clone(negonevec);}
  if (nonobs2col != -1) {
    nonobsgiven2x = as<NumericVector>(data[nonobs2col]);
  } else {nonobsgiven2x = clone(negonevec);}
  if (alive2col != -1) {
    alivegiven2x = as<NumericVector>(data[alive2col]);
  } else {alivegiven2x = clone(negonevec);}
  if (dead2col != -1) {
    deadgiven2x = as<NumericVector>(data[dead2col]);
  } else {deadgiven2x = clone(negonevec);}
  if (xcol != -1) {
    xpos2x = as<NumericVector>(data[xcol]);
  } else {xpos2x = clone(negonevec);}
  if (ycol != -1) {
    ypos2x = as<NumericVector>(data[ycol]);
  } else {ypos2x = clone(negonevec);}
  
  if (sizea3col != -1) {
    sizea3x = as<NumericVector>(data[sizea3col]);
  } else {sizea3x = clone(zerovec);}
  if (sizeb3col != -1) {
    sizeb3x = as<NumericVector>(data[sizeb3col]);
  } else {sizeb3x = clone(zerovec);}
  if (sizec3col != -1) {
    sizec3x = as<NumericVector>(data[sizec3col]);
  } else {sizec3x = clone(zerovec);}
  if (repstra3col != -1) {
    repstra3x = as<NumericVector>(data[repstra3col]);
  } else {repstra3x = clone(zerovec);}
  if (repstrb3col != -1) {
    repstrb3x = as<NumericVector>(data[repstrb3col]);
  } else {repstrb3x = clone(zerovec);}
  if (feca3col != -1) {
    feca3x = as<NumericVector>(data[feca3col]);
  } else {feca3x = clone(zerovec);}
  if (fecb3col != -1) {
    fecb3x = as<NumericVector>(data[fecb3col]);
  } else {fecb3x = clone(zerovec);}
  if (juv3col != -1) {
    juvgiven3x = as<NumericVector>(data[juv3col]);
  } else {juvgiven3x = clone(zerovec);}
  if (obs3col != -1) {
    obsgiven3x = as<NumericVector>(data[obs3col]);
  } else {obsgiven3x = clone(negonevec);}
  if (nonobs3col != -1) {
    nonobsgiven3x = as<NumericVector>(data[nonobs3col]);
  } else {nonobsgiven3x = clone(negonevec);}
  if (alive3col != -1) {
    alivegiven3x = as<NumericVector>(data[alive3col]);
  } else {alivegiven3x = clone(negonevec);}
  if (dead3col != -1) {
    deadgiven3x = as<NumericVector>(data[dead3col]);
  } else {deadgiven3x = clone(negonevec);}
  
  if (stage2col != -1) {
    stage2x = as<StringVector>(data[stage2col]);
  } else {stage2x = as<StringVector>(clone(zerovec));}
  if (stage3col != -1) {
    stage3x = as<StringVector>(data[stage3col]);
  } else {stage3x = as<StringVector>(clone(zerovec));}
  
  Rcpp::IntegerVector rowid (ndflength);
  
  Rcpp::StringVector popid (ndflength);
  Rcpp::StringVector patchid (ndflength);
  Rcpp::StringVector individ (ndflength);
  
  Rcpp::IntegerVector popid_int (ndflength);
  Rcpp::IntegerVector patchid_int (ndflength);
  Rcpp::IntegerVector individ_int (ndflength);
  
  Rcpp::NumericVector censor2 (ndflength);
  Rcpp::NumericVector xpos2 (ndflength, 0.0);
  Rcpp::NumericVector ypos2 (ndflength, 0.0);
  Rcpp::NumericVector sizea2 (ndflength, 0.0);
  Rcpp::NumericVector sizeb2 (ndflength, 0.0);
  Rcpp::NumericVector sizec2 (ndflength, 0.0);
  Rcpp::NumericVector sizeadded2 (ndflength, 0.0);
  Rcpp::NumericVector repstra2 (ndflength, 0.0);
  Rcpp::NumericVector repstrb2 (ndflength, 0.0);
  Rcpp::NumericVector repstradded2 (ndflength, 0.0);
  Rcpp::NumericVector feca2 (ndflength, 0.0);
  Rcpp::NumericVector fecb2 (ndflength, 0.0);
  Rcpp::NumericVector fecadded2 (ndflength, 0.0);
  
  arma::ivec year2 (ndflength, fill::zeros);
  
  Rcpp::NumericVector indcova2 (ndflength, 0.0);
  Rcpp::NumericVector indcovb2 (ndflength, 0.0);
  Rcpp::NumericVector indcovc2 (ndflength, 0.0);
  
  Rcpp::IntegerVector indcova2_int (ndflength, NA_INTEGER);
  Rcpp::IntegerVector indcovb2_int (ndflength, NA_INTEGER);
  Rcpp::IntegerVector indcovc2_int (ndflength, NA_INTEGER);
  
  Rcpp::StringVector indcova2_str (ndflength, NA_STRING);
  Rcpp::StringVector indcovb2_str (ndflength, NA_STRING);
  Rcpp::StringVector indcovc2_str (ndflength, NA_STRING);
  
  Rcpp::IntegerVector repstatus2 (ndflength, 0);
  Rcpp::IntegerVector fecstatus2 (ndflength, 0);
  Rcpp::IntegerVector obsstatus2 (ndflength, 0);
  Rcpp::NumericVector juvgiven2 (ndflength, 0.0);
  Rcpp::IntegerVector matstat2 (ndflength, 0);
  
  censor2.fill(crazycensor);
  
  Rcpp::NumericVector censor3 (ndflength);
  Rcpp::NumericVector xpos3 (ndflength, 0.0);
  Rcpp::NumericVector ypos3 (ndflength, 0.0);
  Rcpp::NumericVector sizea3 (ndflength, 0.0);
  Rcpp::NumericVector sizeb3 (ndflength, 0.0);
  Rcpp::NumericVector sizec3 (ndflength, 0.0);
  Rcpp::NumericVector sizeadded3 (ndflength, 0.0);
  Rcpp::NumericVector repstra3 (ndflength, 0.0);
  Rcpp::NumericVector repstrb3 (ndflength, 0.0);
  Rcpp::NumericVector repstradded3 (ndflength, 0.0);
  Rcpp::NumericVector feca3 (ndflength, 0.0);
  Rcpp::NumericVector fecb3 (ndflength, 0.0);
  Rcpp::NumericVector fecadded3 (ndflength, 0.0);
  
  Rcpp::NumericVector indcova3 (ndflength, 0.0);
  Rcpp::NumericVector indcovb3 (ndflength, 0.0);
  Rcpp::NumericVector indcovc3 (ndflength, 0.0);
  
  Rcpp::IntegerVector indcova3_int (ndflength, NA_INTEGER);
  Rcpp::IntegerVector indcovb3_int (ndflength, NA_INTEGER);
  Rcpp::IntegerVector indcovc3_int (ndflength, NA_INTEGER);
  
  Rcpp::StringVector indcova3_str (ndflength, NA_STRING);
  Rcpp::StringVector indcovb3_str (ndflength, NA_STRING);
  Rcpp::StringVector indcovc3_str (ndflength, NA_STRING);
  
  Rcpp::IntegerVector repstatus3 (ndflength, 0);
  Rcpp::IntegerVector fecstatus3 (ndflength, 0);
  Rcpp::IntegerVector obsstatus3 (ndflength, 0);
  Rcpp::NumericVector juvgiven3 (ndflength, 0.0);
  Rcpp::IntegerVector matstat3 (ndflength, 0);
  
  censor3.fill(crazycensor);
  
  Rcpp::NumericVector censor1 (ndflength);
  Rcpp::NumericVector xpos1 (ndflength, 0.0);
  Rcpp::NumericVector ypos1 (ndflength, 0.0);
  Rcpp::NumericVector sizea1 (ndflength, 0.0);
  Rcpp::NumericVector sizeb1 (ndflength, 0.0);
  Rcpp::NumericVector sizec1 (ndflength, 0.0);
  Rcpp::NumericVector sizeadded1 (ndflength, 0.0);
  Rcpp::NumericVector repstra1 (ndflength, 0.0);
  Rcpp::NumericVector repstrb1 (ndflength, 0.0);
  Rcpp::NumericVector repstradded1 (ndflength, 0.0);
  Rcpp::NumericVector feca1 (ndflength, 0.0);
  Rcpp::NumericVector fecb1 (ndflength, 0.0);
  Rcpp::NumericVector fecadded1 (ndflength, 0.0);
  
  Rcpp::NumericVector indcova1 (ndflength, 0.0);
  Rcpp::NumericVector indcovb1 (ndflength, 0.0);
  Rcpp::NumericVector indcovc1 (ndflength, 0.0);
  
  Rcpp::IntegerVector indcova1_int (ndflength, NA_INTEGER);
  Rcpp::IntegerVector indcovb1_int (ndflength, NA_INTEGER);
  Rcpp::IntegerVector indcovc1_int (ndflength, NA_INTEGER);
  
  Rcpp::StringVector indcova1_str (ndflength, NA_STRING);
  Rcpp::StringVector indcovb1_str (ndflength, NA_STRING);
  Rcpp::StringVector indcovc1_str (ndflength, NA_STRING);
  
  Rcpp::IntegerVector repstatus1 (ndflength, 0);
  Rcpp::IntegerVector fecstatus1 (ndflength, 0);
  Rcpp::IntegerVector obsstatus1 (ndflength, 0);
  Rcpp::NumericVector juvgiven1 (ndflength, 0.0);
  Rcpp::IntegerVector matstat1 (ndflength, 0);
  
  censor1.fill(crazycensor);
  
  // Variables to check whether censor has been checked and set at each step
  arma::uvec indivnum (ndflength, fill::zeros);
  arma::uvec censor2check (ndflength, fill::zeros);
  
  // Derived variables requiring extra looping or other control parameters
  arma::ivec firstseen (ndflength);
  arma::ivec lastseen (ndflength);
  arma::ivec obsage (ndflength);
  arma::ivec obslifespan (ndflength);
  Rcpp::NumericVector alive1 (ndflength, 0.0);
  Rcpp::NumericVector alive2 (ndflength, 0.0);
  Rcpp::NumericVector alive3 (ndflength, 0.0);
  firstseen.fill(-1);
  lastseen.fill(-1);
  
  Rcpp::StringVector stage1 (ndflength);
  Rcpp::StringVector stage2 (ndflength);
  Rcpp::StringVector stage3 (ndflength);
  Rcpp::NumericVector stage1num (ndflength, 0.0);
  Rcpp::NumericVector stage2num (ndflength, 0.0);
  Rcpp::NumericVector stage3num (ndflength, 0.0);
  
  double stagesize1 {0.0};
  double stagesize2 {0.0};
  double stagesize3 {0.0};
  double stagesize1b {0.0};
  double stagesize2b {0.0};
  double stagesize3b {0.0};
  double stagesize1c {0.0};
  double stagesize2c {0.0};
  double stagesize3c {0.0};
  
  arma::uvec stagemini1;
  arma::uvec stagemaxi1;
  arma::uvec stagemini2;
  arma::uvec stagemaxi2;
  arma::uvec stagemini3;
  arma::uvec stagemaxi3;
  arma::uvec stageobs;
  arma::uvec stagerep;
  arma::uvec stagemat;
  arma::uvec cs1;
  arma::uvec cs2;
  arma::uvec cs3;
  arma::uvec cs4;
  int choicestage {0};
  int fsyear_check {0};
  
  //Rcout << "C1" << endl;
  
  // Main loop creating new dataset rows
  // Establishes state in time t for all cases in which individual is observed
  for (int i = 0; i < norows; i++) { // i is row in old dataset
    if (i % 5 == 0) Rcpp::checkUserInterrupt();
    
    for (int j = 0; j < noyears; j++) {
      // Establishes place marker for vectors corresponding to current year
      if (year2x[i] == yearall2x[j]) currentyear = j;
    }
    
    // Establishes place marker corresponding to current individual
    currentindiv = -1;
    for (int k = 0; k < noindivs; k++) {
      if (individx[i] == allindivs[k]) {
        currentindiv = k;
        indivnum[i] = k;
      }
    }
    
    // Establishes row in new dataset
    ndfindex = (noyears * currentindiv) + currentyear;
    
    if (NumericVector::is_na(sizea2x[i])) {
      sizea20x[i] = 0.0;
      if (NAas0) {sizea2x[i] = 0.0;}
    } else {sizea20x[i] = sizea2x[i];}
    if (NumericVector::is_na(sizeb2x[i])) {
      sizeb20x[i] = 0.0;
      if (NAas0) {sizeb2x[i] = 0.0;}
    } else {sizeb20x[i] = sizeb2x[i];}
    if (NumericVector::is_na(sizec2x[i])) {
      sizec20x[i] = 0.0;
      if (NAas0) {sizec2x[i] = 0.0;}
    } else {sizec20x[i] = sizec2x[i];}
    if (NumericVector::is_na(repstra2x[i])) {
      repstra20x[i] = 0.0;
      if (NAas0) {repstra2x[i] = 0.0;}
    } else {repstra20x[i] = repstra2x[i];}
    if (NumericVector::is_na(repstrb2x[i])) {
      repstrb20x[i] = 0.0;
      if (NAas0) {repstrb2x[i] = 0.0;}
    } else {repstrb20x[i] = repstrb2x[i];}
    if (NumericVector::is_na(feca2x[i])) {
      feca20x[i] = 0.0;
      if (NAas0) {feca2x[i] = 0.0;}
    } else {feca20x[i] = feca2x[i];}
    if (NumericVector::is_na(fecb2x[i])) {
      fecb20x[i] = 0.0;
      if (NAas0) {fecb2x[i] = 0.0;}
    } else {fecb20x[i] = fecb2x[i];}
    
    if (NumericVector::is_na(juvgiven2x[i])) {
      juvgiven20x[i] = 0.0;
    } else if (juvgiven2x[i] != 0) {
      juvgiven20x[i] = 1.0;
    } else {juvgiven20x[i] = 0.0;}
    
    // Develops censoring variable
    if (censbool && censorcol != -1) {
      // Provides replacements in cases where NA designates data to keep
      if (NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0.0;
        censor2check[ndfindex] = 1;
      } else {
        censor2[ndfindex] = 1.0;
        censor2check[ndfindex] = 1;
      }
    } else if (censorcol != -1) {
      if (censorkeep == 0 && NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 1.0;
        censor2check[ndfindex] = 1;
      } else if (censorkeep == 1 && NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0.0;
        censor2check[ndfindex] = 1;
      } else {
        censor2[ndfindex] = censor2x[i];
        censor2check[ndfindex] = 1;
      }
    }
    
    rowid[ndfindex] = i;
    
    if (popid_type > 0) { 
      popid_int[ndfindex] = popidx_int[i];
    } else {
      popid[ndfindex] = popidx[i];
    }
    
    if (patchid_type > 0) { 
      patchid_int[ndfindex] = patchidx_int[i];
    } else {
      patchid[ndfindex] = patchidx[i];
    }
    
    if (individ_type > 0) { 
      individ_int[ndfindex] = individx_int[i];
    }
    individ[ndfindex] = allindivs[currentindiv];
    
    year2[ndfindex] = yearall2x[currentyear];
    xpos2[ndfindex] = xpos2x[i];
    ypos2[ndfindex] = ypos2x[i];
    sizea2[ndfindex] = sizea2x[i];
    sizeb2[ndfindex] = sizeb2x[i];
    sizec2[ndfindex] = sizec2x[i];
    sizeadded2[ndfindex] = sizea20x[i] + sizeb20x[i] + sizec20x[i];
    
    repstra2[ndfindex] = repstra2x[i];
    repstrb2[ndfindex] = repstrb2x[i];
    repstradded2[ndfindex] = repstra20x[i] + (repstrb20x[i] * repstrrel);
    
    feca2[ndfindex] = feca2x[i];
    fecb2[ndfindex] = fecb2x[i];
    fecadded2[ndfindex] = feca20x[i] + (fecb20x[i] * fecrel);
    
    //Rcout << "jpf D" << endl;
  
    if (indcova_type == 3) {
      indcova2[ndfindex] = indcova2x[i];
    } else if (indcova_type == 0) {
      indcova2_str[ndfindex] = indcova2x_str[i];
    } else {
      indcova2_int[ndfindex] = indcova2x_int[i];
    }
    
    if (indcovb_type == 3) {
      indcovb2[ndfindex] = indcovb2x[i];
    } else if (indcovb_type == 0) {
      indcovb2_str[ndfindex] = indcovb2x_str[i];
    } else {
      indcovb2_int[ndfindex] = indcovb2x_int[i];
    }
    
    if (indcovc_type == 3) {
      indcovc2[ndfindex] = indcovc2x[i];
    } else if (indcovc_type == 0) {
      indcovc2_str[ndfindex] = indcovc2x_str[i];
    } else {
      indcovc2_int[ndfindex] = indcovc2x_int[i];
    }
    
    //Rcout << "jpf E" << endl;
  
    if (repstradded2[ndfindex] > 0) {
      repstatus2[ndfindex] = 1;} else {repstatus2[ndfindex] = 0;
    }
    if (NumericVector::is_na(obsgiven2x[i])) {
      obsstatus2[ndfindex] = 0;
    } else if (obsgiven2x[i] > 0) {
      obsstatus2[ndfindex] = 1;
    } else if (obsgiven2x[i] == -1 && (sizeadded2[ndfindex] + repstradded2[ndfindex]) > 0) {
      obsstatus2[ndfindex] = 1;
    } else {obsstatus2[ndfindex] = 0;}
    
    if (nonobsgiven2x[i] >= 0) {
      obsstatus2[ndfindex] = 0;
    } else if (nonobsgiven2x[i] == 0) {
      obsstatus2[ndfindex] = 1;
    }
    
    juvgiven2[ndfindex] = juvgiven20x[i];
    matstat2[ndfindex] = 1 - juvgiven2[ndfindex];
    
    if (alivegiven2x[i] > 0) {
      alive2[ndfindex] = 1.0;
    } else if (alivegiven2x[i] == 0) {
      alive2[ndfindex] = 0.0;
    }
    
    if (deadgiven2x[i] > 0) {
      alive2[ndfindex] = 0.0;
    } else if (deadgiven2x[i] == 0) {
      alive2[ndfindex] = 1.0;
    }
    
    livcheck2 = sizeadded2[ndfindex] + repstradded2[ndfindex] + obsstatus2[ndfindex];
    fsyear_check = currentyear + firstyear;
    
    if (livcheck2 > 0) {
      if (firstseenx[currentindiv] == -1) {
        firstseenx[currentindiv] = fsyear_check;
        lastseenx[currentindiv] = fsyear_check;
        alive2[ndfindex] = 1.0;
      } else if (firstseenx[currentindiv] > fsyear_check) {
        firstseenx[currentindiv] = fsyear_check;
        alive2[ndfindex] = 1.0;
      } else if (lastseenx[currentindiv] < fsyear_check) {
        lastseenx[currentindiv] = fsyear_check;
        alive2[ndfindex] = 1.0;
      }
    }
    
    if (alive2[ndfindex] == 1.0 && matstat2[ndfindex] == 1) {
      matstat3[ndfindex] = 1;
    }
    
    if (stage2col != -1 && alive2[ndfindex] == 1.0) {
      stage2[ndfindex] = stage2x[i];
    } else if (stassign) {stage2[ndfindex] = "NotAlive";}
    
    
    // Time t+1 for last time t (2nd to last time) in cases with t+1 columns provided
    if (currentyear == (noyears - 1)) {
      if (censbool && censorcol != -1) { // Censoring variable for final time
        if (NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0.0;
        } else {
          censor3[ndfindex] = 1.0;
        }
      } else if (censorcol != -1) {
        if (censorkeep == 0 && NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 1.0;
        } else if (censorkeep == 1 && NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0.0;
        } else {
          censor3[ndfindex] = censor2x[i];
        }
      }
      
      if (NumericVector::is_na(juvgiven3x[i])) {
        juvgiven30x[i] = 0.0;
      } else if (juvgiven3x[i] != 0) {
        juvgiven30x[i] = 1.0;
      } else {juvgiven30x[i] = 0.0;}
      
      if (sizea3col != -1) {
        if (NumericVector::is_na(sizea3x[i])) {
          sizea30x[i] = 0.0;
          if (NAas0) {sizea3x[i] = 0.0;}
        } else {sizea30x[i] = sizea3x[i];}
        sizea3[ndfindex] = sizea3x[i];
      }
      if (sizeb3col != -1) {
        if (NumericVector::is_na(sizeb3x[i])) {
          sizeb30x[i] = 0.0;
          if (NAas0) {sizeb3x[i] = 0.0;}
        } else {sizeb30x[i] = sizeb3x[i];}
        sizeb3[ndfindex] = sizeb3x[i];
      }
      if (sizec3col != -1) {
        if (NumericVector::is_na(sizec3x[i])) {
          sizec30x[i] = 0.0;
          if (NAas0) {sizec3x[i] = 0.0;}
        } else {sizec30x[i] = sizec3x[i];}
        sizec3[ndfindex] = sizec3x[i];
      }
      sizeadded3[ndfindex] = sizea30x[i] + sizeb30x[i] + sizec30x[i];
      
      if (repstra3col != -1) {
        if (NumericVector::is_na(repstra3x[i])) {
          repstra30x[i] = 0.0;
          if (NAas0) {repstra3x[i] = 0.0;}
        } else {repstra30x[i] = repstra3x[i];}
        repstra3[ndfindex] = repstra3x[i];
      }
      if (repstrb3col != -1) {
        if (NumericVector::is_na(repstrb3x[i])) {
          repstrb30x[i] = 0.0;
          if (NAas0) {repstrb3x[i] = 0.0;}
        } else {repstrb30x[i] = repstrb3x[i];}
        repstrb3[ndfindex] = repstrb3x[i];
      }
      repstradded3[ndfindex] = repstra30x[i] + (repstrb30x[i] * repstrrel);
      
      if (feca3col != -1) {
        if (NumericVector::is_na(feca3x[i])) {
          feca30x[i] = 0.0;
          if (NAas0) {feca3x[i] = 0.0;}
        } else {feca30x[i] = feca3x[i];}
        feca3[ndfindex] = feca3x[i];
      }
      if (fecb3col != -1) {
        if (NumericVector::is_na(fecb3x[i])) {
          fecb30x[i] = 0.0;
          if (NAas0) {feca3x[i] = 0.0;}
        } else {fecb30x[i] = fecb3x[i];}
        fecb3[ndfindex] = fecb3x[i];
      }
      fecadded3[ndfindex] = feca30x[i] + (fecb30x[i] * fecrel);
      if (fecadded3[ndfindex] > 0) {fecstatus3[ndfindex] = 1;}
      
      if (repstradded3[ndfindex] > 0) {
        repstatus3[ndfindex] = 1;
      } else {
        repstatus3[ndfindex] = 0;
      }
      
      //Rcout << "jpf F" << endl;
      
      if (indcova3col != -1) {
        if (indcova_type == 3) {
          indcova3[ndfindex] = indcova3x[i];
        } else if (indcova_type == 0) {
          indcova3_str[ndfindex] = indcova3x_str[i];
        } else {
          indcova3_int[ndfindex] = indcova3x_int[i];
        }
      }
      
      if (indcovb3col != -1) {
        if (indcovb_type == 3) {
          indcovb3[ndfindex] = indcovb3x[i];
        } else if (indcovb_type == 0) {
          indcovb3_str[ndfindex] = indcovb3x_str[i];
        } else {
          indcovb3_int[ndfindex] = indcovb3x_int[i];
        }
      }
      
      if (indcovc3col != -1) {
        if (indcovc_type == 3) {
          indcovc3[ndfindex] = indcovc3x[i];
        } else if (indcovc_type == 0) {
          indcovc3_str[ndfindex] = indcovc3x_str[i];
        } else {
          indcovc3_int[ndfindex] = indcovc3x_int[i];
        }
      }
      
      if (NumericVector::is_na(obsgiven3x[i])) {
        obsstatus3[ndfindex] = 0;
      } else if (obsgiven3x[i] > 0) {
        obsstatus3[ndfindex] = 1;
      } else if (obsgiven3x[i] == -1 && (sizeadded3[ndfindex] + repstradded3[ndfindex]) > 0) {
        obsstatus3[ndfindex] = 1;
      } else {obsstatus3[ndfindex] = 0;}
      
      if (nonobsgiven3x[i] >= 0) {
        obsstatus3[ndfindex] = 0;
      } else if (nonobsgiven3x[i] == 0) {
        obsstatus3[ndfindex] = 1;
      }
      
      juvgiven3[ndfindex] = juvgiven30x[i];
      matstat3[ndfindex] = 1 - juvgiven3[ndfindex];
      
      if (alivegiven3x[i] > 0) {
        alive3[ndfindex] = 1.0;
      } else if (alivegiven3x[i] == 0) {
        alive3[ndfindex] = 0.0;
      }
      if (deadgiven3x[i] > 0) {
        alive3[ndfindex] = 0.0;
      } else if (deadgiven3x[i] == 0) {
        alive3[ndfindex] = 1.0;
      }
      
      livcheck3 = sizeadded3[ndfindex] + repstradded3[ndfindex] + obsstatus3[ndfindex];
      fsyear_check = currentyear + firstyear + 1;
      
      if (livcheck3 > 0) {
        if (firstseenx[currentindiv] == -1) {
          firstseenx[currentindiv] = fsyear_check;
          lastseenx[currentindiv] = fsyear_check;
          alive3[ndfindex] = 1.0;
        } else if (firstseenx[currentindiv] > fsyear_check && alive2[ndfindex] < 1.0) {
          firstseenx[currentindiv] = fsyear_check;
          alive3[ndfindex] = 1.0;
        } else if (lastseenx[currentindiv] < fsyear_check) {
          lastseenx[currentindiv] = fsyear_check;
          alive3[ndfindex] = 1.0;
        }
        
        if (juv3col == -1 && repstradded2[ndfindex] > 0) {
          matstat3[ndfindex] = 1;
        }
      }
      
      if (stage3col != -1 && alive3[ndfindex] == 1.0) {
        stage3[ndfindex] = stage2x[i];
      } else if (stassign) {stage3[ndfindex] = "NotAlive";}
    } // End of currentyear if statement
  } // End of i loop
  
  //Rcout << "jpf G" << endl;
  
  // Loop determining most states in time t+1 and t-1, and stages in all times
  for (int i = 0; i < ndflength; i++) { // i refers to rows in final dataset
    if (i % 5 == 0) Rcpp::checkUserInterrupt();
    
    // Corrects info for individuals unobserved for long periods
    if (i > 0 && rowid[i] == 0) {
      if (year2[i-1] < lastseen[i-1] && (year2[i-1] + 1) < (firstyear + noyears)) {
        
        if (individ_type > 0) { 
          individ_int[i] = individ_int[i-1];
        }
        individ[i] = individ[i-1];
        
        if (popid_type > 0) { 
          popid_int[i] = popid_int[i-1];
        } else {
          popid[i] = popid[i-1];
        }
        
        if (patchid_type > 0) { 
          patchid_int[i] = patchid_int[i-1];
        } else {
          patchid[i] = patchid[i-1];
        }
        
        year2[i] = year2[i-1] + 1;
        xpos2[i] = xpos2[i-1];
        ypos2[i] = ypos2[i-1];
        
        if (matstat2[i-1] == 1) {
          matstat2[i] = 1;
          
          if (year2[i+1] <= lastseen[i]) {
            matstat3[i] = 1;
          }
        }
      }
    }
    
    currentindiv = -1;
    for (int k = 0; k < noindivs; k++) {
      if (individ[i] == allindivs[k]) currentindiv = k;
    }
    
    if (currentindiv != -1) { // Limits to only real individuals in the dataset
      if (year2[i] <= lastseenx[currentindiv] && year2[i] >= firstseenx[currentindiv] && 
          year2[i] < (firstyear + noyears)) {
        firstseen[i] = firstseenx[currentindiv];
        lastseen[i] = lastseenx[currentindiv];
        
        obsage[i] = year2[i] - firstseen[i];
        obslifespan[i] = lastseen[i] - firstseen[i];
      }
      
      if (year2[i] >= firstseen[i] && year2[i] <= lastseen[i] && alive2[i] == 0.0) {
        alive2[i] = 1.0;
      } else if (alive2[i] == 0.0) {
        alive2[i] = 0.0;
      }
      
      if ((year2[i] + 1) >= firstseen[i] && (year2[i] + 1) <= lastseen[i] && alive3[i] == 0.0) {
        alive3[i] = 1.0;
      } else if (alive3[i] == 0.0) {
        alive3[i] = 0.0;
      }
      
      currentyear = year2[i] - firstyear;
      
      if (fecadded2[i] > 0) {
        fecstatus2[i] = 1;
      } else {
        fecstatus2[i] = 0;
      }
      
      if (stage2col != -1 && alive2[i] == 1.0) {
        if (obsstatus2[i] == 0) {
          stage2[i] = "NotObserved";
        }
      } else if (stage2col != -1 && alive2[i] == 0.0 && stassign) {
        stage2[i] = "NotAlive";
      }
      
      if (currentyear < (noyears - 1)) { 
        nextyrindex = (noyears * currentindiv) + (currentyear + 1);
        
        if (censor2[i] == censorkeep && alive2[i] == 1.0) {
          censor3[i] = censorkeep;
        } else {
          censor3[i] = censor2[nextyrindex];
        }
        
        xpos3[i] = xpos2[nextyrindex];
        ypos3[i] = ypos2[nextyrindex];
        
        sizea3[i] = sizea2[nextyrindex];
        sizeb3[i] = sizeb2[nextyrindex];
        sizec3[i] = sizec2[nextyrindex];
        sizeadded3[i] = sizeadded2[nextyrindex];
        
        repstra3[i] = repstra2[nextyrindex];
        repstrb3[i] = repstrb2[nextyrindex];
        repstradded3[i] = repstradded2[nextyrindex];
        
        feca3[i] = feca2[nextyrindex];
        fecb3[i] = fecb2[nextyrindex];
        fecadded3[i] = fecadded2[nextyrindex];
        
        //Rcout << "jpf H" << endl;
        
        if (indcova_type == 3) {
          indcova3[i] = indcova2[nextyrindex];
        } else if (indcova_type == 0) {
          indcova3_str[i] = indcova2_str[nextyrindex];
        } else {
          indcova3_int[i] = indcova2_int[nextyrindex];
        }
        
        if (indcovb_type == 3) {
          indcovb3[i] = indcovb2[nextyrindex];
        } else if (indcovb_type == 0) {
          indcovb3_str[i] = indcovb2_str[nextyrindex];
        } else {
          indcovb3_int[i] = indcovb2_int[nextyrindex];
        }
        
        if (indcovc_type == 3) {
          indcovc3[i] = indcovc2[nextyrindex];
        } else if (indcova_type == 0) {
          indcovc3_str[i] = indcovc2_str[nextyrindex];
        } else {
          indcovc3_int[i] = indcovc2_int[nextyrindex];
        }
        
        if (fecadded3[i] > 0) {
          fecstatus3[i] = 1;
        } else {
          fecstatus3[i] = 0;
        }
        
        repstatus3[i] = repstatus2[nextyrindex];
        obsstatus3[i] = obsstatus2[nextyrindex];
        juvgiven3[i] = juvgiven2[nextyrindex];
        matstat3[i] = matstat2[nextyrindex];
        
        if (matstat2[i] == 1 && stage3[i] != "NotAlive") {
          matstat3[i] = 1;
        }
        
        if (stage2col != -1 && alive3[i] == 1.0 && stassign) {
          if (obsstatus3[i] == 0) {
            stage3[i] = "NotObserved";
          } else {stage3[i] = stage2[nextyrindex];}
        } else if (stage2col != -1 && alive3[i] == 0.0 && stassign) {
          stage3[i] = "NotAlive";
        }
      }
      
      if (currentyear > 0  && year2[i] < (firstyear + noyears)) {
        prevyrindex = (noyears * currentindiv) + (currentyear - 1);
        
        alive1[i] = alive2[prevyrindex];
        
        if (censor2(i) == censorkeep && alive1(i) == 1.0) {
          if (censbool) {
            censor1(i) = 0;
          } else {
            censor1[i] = censorkeep;
          }
        } else {
          censor1[i] = censor2[prevyrindex];
        }
        
        xpos1[i] = xpos2[prevyrindex];
        ypos1[i] = ypos2[prevyrindex];
        
        sizea1[i] = sizea2[prevyrindex];
        sizeb1[i] = sizeb2[prevyrindex];
        sizec1[i] = sizec2[prevyrindex];
        sizeadded1[i] = sizeadded2[prevyrindex];
        
        repstra1[i] = repstra2[prevyrindex];
        repstrb1[i] = repstrb2[prevyrindex];
        repstradded1[i] = repstradded2[prevyrindex];
        
        feca1[i] = feca2[prevyrindex];
        fecb1[i] = fecb2[prevyrindex];
        fecadded1[i] = fecadded2[prevyrindex];
        
        //Rcout << "jpf I" << endl;
        
        if (indcova_type == 3) {
          indcova1[i] = indcova2[prevyrindex];
        } else if (indcova_type == 0) {
          indcova1_str[i] = indcova2_str[prevyrindex];
        } else {
          indcova1_int[i] = indcova2_int[prevyrindex];
        }
        
        if (indcovb_type == 3) {
          indcovb1[i] = indcovb2[prevyrindex];
        } else if (indcovb_type == 0) {
          indcovb1_str[i] = indcovb2_str[prevyrindex];
        } else {
          indcovb1_int[i] = indcovb2_int[prevyrindex];
        }
        
        if (indcovc_type == 3) {
          indcovc1[i] = indcovc2[prevyrindex];
        } else if (indcova_type == 0) {
          indcovc1_str[i] = indcovc2_str[prevyrindex];
        } else {
          indcovc1_int[i] = indcovc2_int[prevyrindex];
        }
        
        if (fecadded1[i] > 0) {
          fecstatus1[i] = 1;
        } else {
          fecstatus1[i] = 0;
        }
        
        repstatus1[i] = repstatus2[prevyrindex];
        obsstatus1[i] = obsstatus2[prevyrindex];
        juvgiven1[i] = juvgiven2[prevyrindex];
        matstat1[i] = matstat2[prevyrindex];
        
        if (stage2col != -1 && alive1[i] == 1.0 && stassign) {
          if (obsstatus1[i] == 0) {
            stage1[i] = "NotObserved";
          } else {
            stage1[i] = stage2[prevyrindex];
          }
        } else if (stage2col != -1 && alive1[i] == 0.0 && stassign) {
          stage1[i] = "NotAlive";
        }
      }
      
      // Coordinate corrections
      if (xpos1[i] != xpos2[i] && xpos1[i] == 0.0) {xpos1[i] = xpos2[i];}
      if (xpos3[i] != xpos2[i] && xpos3[i] == 0.0) {xpos3[i] = xpos2[i];}
      if (ypos1[i] != ypos2[i] && ypos1[i] == 0.0) {ypos1[i] = ypos2[i];}
      if (ypos3[i] != ypos2[i] && ypos3[i] == 0.0) {ypos3[i] = ypos2[i];}
      
      //Rcout << "jpf J" << endl;
      
      // Stage assignments
      if (stassign && stage2col == -1) {
        if (stszcol == 8) {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          stagesize1b = sizeb1[i];
          stagesize2b = sizeb2[i];
          stagesize3b = sizeb3[i];
          stagesize1c = sizec1[i];
          stagesize2c = sizec2[i];
          stagesize3c = sizec3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini1c = find(insfszminarmac < stagesize1c);
          arma::uvec stagemaxi1c = find(insfszmaxarmac >= stagesize1c);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini2c = find(insfszminarmac < stagesize2c);
          arma::uvec stagemaxi2c = find(insfszmaxarmac >= stagesize2c);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          arma::uvec stagemini3c = find(insfszminarmac < stagesize3c);
          arma::uvec stagemaxi3c = find(insfszmaxarmac >= stagesize3c);
          
          arma::uvec stagemini1d = intersect(stagemini1a, stagemini1b);
          stagemini1 = intersect(stagemini1c, stagemini1d);
          arma::uvec stagemaxi1d = intersect(stagemaxi1a, stagemaxi1b);
          stagemaxi1 = intersect(stagemaxi1c, stagemaxi1d);
          arma::uvec stagemini2d = intersect(stagemini2a, stagemini2b);
          stagemini2 = intersect(stagemini2c, stagemini2d);
          arma::uvec stagemaxi2d = intersect(stagemaxi2a, stagemaxi2b);
          stagemaxi2 = intersect(stagemaxi2c, stagemaxi2d);
          arma::uvec stagemini3d = intersect(stagemini3a, stagemini3b);
          stagemini3 = intersect(stagemini3c, stagemini3d);
          arma::uvec stagemaxi3d = intersect(stagemaxi3a, stagemaxi3b);
          stagemaxi3 = intersect(stagemaxi3c, stagemaxi3d);
          
        } else if (stszcol == 7) {
          stagesize1 = sizeb1[i];
          stagesize2 = sizeb2[i];
          stagesize3 = sizeb3[i];
          stagesize1b = sizec1[i];
          stagesize2b = sizec2[i];
          stagesize3b = sizec3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 6) {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          stagesize1b = sizec1[i];
          stagesize2b = sizec2[i];
          stagesize3b = sizec3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 5) {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          stagesize1b = sizeb1[i];
          stagesize2b = sizeb2[i];
          stagesize3b = sizeb3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 4) {
          stagesize1 = sizeadded1[i];
          stagesize2 = sizeadded2[i];
          stagesize3 = sizeadded3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else if (stszcol == 3) {
          stagesize1 = sizec1[i];
          stagesize2 = sizec2[i];
          stagesize3 = sizec3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
        } else if (stszcol == 2) {
          stagesize1 = sizeb1[i];
          stagesize2 = sizeb2[i];
          stagesize3 = sizeb3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
        }
        
        // Stage 2
        stagerep = find(inrepstatarma == repstatus2[i]);
        stagemat = find(inmatstatarma == matstat2[i]);
        stageobs = find(inobsstatarma == obsstatus2[i]);
        
        cs1 = intersect(stagemini2, stagemaxi2);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem == 1 && alive2[i] == 1.0) {
          choicestage = instageid(cs4[0]) - 1;
          stage2num[i] = choicestage + 1;
          
          stage2[i] = sfname[choicestage];
        } else if (cs4.n_elem == 0 && alive2[i] == 1.0) {
          stage2[i] = "NoMatch";
          
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "Some stages in the dataset do not match stage descriptions in the stageframe.");
          }
        } else if (alive2[i] != 1.0)  {
          stage2[i] = "NotAlive";
        } else if (cs4.n_elem > 1) {
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "Some stages in stageframe have same description. All stages should be unique.");
          }
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
        
        // Stage 1
        stagerep = find(inrepstatarma == repstatus1[i]);
        stagemat = find(inmatstatarma == matstat1[i]);
        stageobs = find(inobsstatarma == obsstatus1[i]);
        
        cs1 = intersect(stagemini1, stagemaxi1);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem == 1 && alive1[i] == 1.0) {
          choicestage = instageid(cs4[0]) - 1;
          stage1num[i] = choicestage + 1;
          
          stage1[i] = sfname[choicestage];
        } else if (cs4.n_elem == 0 && alive1[i] == 1.0) {
          stage1[i] = "NoMatch";
          
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "Some stages in the dataset do not match stage descriptions in the stageframe.");
          }
        } else if (alive1[i] != 1.0) {
          stage1[i] = "NotAlive";
          matstat1[i] = 0;
        } else if (cs4.n_elem > 1) {
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "Some stages in stageframe have same description. All stages should be unique.");
          }
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
        
        // Stage 3
        stagerep = find(inrepstatarma == repstatus3[i]);
        stagemat = find(inmatstatarma == matstat3[i]);
        stageobs = find(inobsstatarma == obsstatus3[i]);
        
        cs1 = intersect(stagemini3, stagemaxi3);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        // Create exceptions based on stage assignment problems in time t+1
        if (cs4.n_elem == 1 && alive3[i] == 1.0) {
          choicestage = instageid(cs4[0]) - 1;
          stage3num[i] = choicestage + 1;
          
          stage3[i] = sfname[choicestage];
        } else if (alive3[i] != 1.0) {
          stage3[i] = "NotAlive";
        } else if (cs4.n_elem == 0) {
          stage3[i] = "NoMatch";
          
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "Some stages in the dataset do not match stage descriptions in the stageframe.");
          }
          if (i > 836) throw Rcpp::exception("made it this far 7f", false); 
        } else if (cs4.n_elem > 1) {
          if (!quiet) {
            Rf_warningcall(R_NilValue,
              "Some stages in stageframe have same description. All stages should be unique.");
          }
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
      } // stassign if statement
    } // currentindiv if statement
  } // i loop
  
  // Checks if censor variables have been fully and properly assigned
  if (censorcol != -1) {
    arma::uvec censorzeros = find(censor2check == 0);
    
    for (int i = 0; i < noindivs; i++) {
      arma::uvec indivindices = sort(find(indivnum == i));
      
      arma::uvec indivzeros = intersect(indivindices, censorzeros);
      arma::ivec years_utilized = year2.elem(indivindices);
      
      int newcount = static_cast<int>(indivzeros.n_elem);
      
      if (newcount > 0) {
        for (int j = 0; j < newcount; j++) {
          if (censbool) {
            censor2(indivzeros(j)) = 0.0;
          } else {
            censor2(indivzeros(j)) = censorkeep;
          }
          
          arma::uvec currenttracking = find(indivindices == indivzeros(j));
          
          if (currenttracking.n_elem > 0) {
            
            int yearnow = year2(indivindices(currenttracking(0)));
            int yearprior = yearnow - 1;
            int yearnext = yearnow + 1;
            
            arma::uvec priorvec = find(years_utilized == yearprior);
            arma::uvec nextvec = find(years_utilized == yearnext);
            
            if (nextvec.n_elem > 0) {
              if (censbool) {
                censor1(nextvec(0)) = 0.0;
              } else {
                censor1(nextvec(0)) = censorkeep;
              }
            }
            
            if (priorvec.n_elem > 0) {
              if (censbool) {
                censor3(priorvec(0)) = 0.0;
              } else {
                censor3(priorvec(0)) = censorkeep;
              }
            } // priorvec if statement
          } // currenttracking if statement
        } // newcount for loop
      } // newcount if statement
    } // noindivs for loop
  } // end censor correction section
  
  // Output list creation
  List final_output;
  
  //Rcout << "jpf K" << endl;
  
  if (reduce) { 
    bool xpos1_used {false}, ypos1_used {false}, xpos2_used {false};
    bool ypos2_used {false}, xpos3_used {false}, ypos3_used {false};
    bool censor1_used {false}, censor2_used {false}, censor3_used {false};
    bool sizea1_used {false}, sizea2_used {false}, sizea3_used {false};
    bool sizeb1_used {false}, sizeb2_used {false}, sizeb3_used {false};
    bool sizec1_used {false}, sizec2_used {false}, sizec3_used {false};
    bool size1added_used {false}, size2added_used {false}, size3added_used {false};
    bool repstra1_used {false}, repstra2_used {false}, repstra3_used {false};
    bool repstrb1_used {false}, repstrb2_used {false}, repstrb3_used {false};
    bool repstr1added_used {false}, repstr2added_used {false}, repstr3added_used {false};
    bool feca1_used {false}, feca2_used {false}, feca3_used {false};
    bool fecb1_used {false}, fecb2_used {false}, fecb3_used {false};
    bool fec1added_used {false}, fec2added_used {false}, fec3added_used {false};
    bool indcova1_used {false}, indcova2_used {false}, indcova3_used {false};
    bool indcovb1_used {false}, indcovb2_used {false}, indcovb3_used {false};
    bool indcovc1_used {false}, indcovc2_used {false}, indcovc3_used {false};
    bool juvgiven1_used {false}, juvgiven2_used {false}, juvgiven3_used {false};
    
    NumericVector xpos1_u = unique(xpos1);
    if (xpos1_u.length() > 1) {
      if (xpos1_u.length() > 2) {
        xpos1_used = true;
      } else {
        arma::uvec xplu0_arma = find(as<arma::vec>(xpos1_u) == 0.0);
        arma::uvec xplu1_arma = find(as<arma::vec>(xpos1_u) == -1.0);
        
        if (xplu0_arma.n_elem + xplu1_arma.n_elem != 2) xpos1_used = true;
      }
    }
    NumericVector ypos1_u = unique(ypos1);
    if (ypos1_u.length() > 1) {
      if (ypos1_u.length() > 2) {
        ypos1_used = true;
      } else {
        arma::uvec yplu0_arma = find(as<arma::vec>(ypos1_u) == 0.0);
        arma::uvec yplu1_arma = find(as<arma::vec>(ypos1_u) == -1.0);
        
        if (yplu0_arma.n_elem + yplu1_arma.n_elem != 2) ypos1_used = true;
      }
    }
    
    NumericVector xpos2_u = unique(xpos2);
    if (xpos2_u.length() > 1) {
      if (xpos2_u.length() > 2) {
        xpos2_used = true;
      } else {
        arma::uvec xplu0_arma = find(as<arma::vec>(xpos2_u) == 0.0);
        arma::uvec xplu1_arma = find(as<arma::vec>(xpos2_u) == -1.0);
        
        if (xplu0_arma.n_elem + xplu1_arma.n_elem != 2) xpos2_used = true;
      }
    }
    NumericVector ypos2_u = unique(ypos2);
    if (ypos2_u.length() > 1) {
      if (ypos2_u.length() > 2) {
        ypos2_used = true;
      } else {
        arma::uvec yplu0_arma = find(as<arma::vec>(ypos2_u) == 0.0);
        arma::uvec yplu1_arma = find(as<arma::vec>(ypos2_u) == -1.0);
        
        if (yplu0_arma.n_elem + yplu1_arma.n_elem != 2) ypos2_used = true;
      }
    }
    
    NumericVector xpos3_u = unique(xpos3);
    if (xpos3_u.length() > 1) {
      if (xpos3_u.length() > 2) {
        xpos3_used = true;
      } else {
        arma::uvec xplu0_arma = find(as<arma::vec>(xpos3_u) == 0.0);
        arma::uvec xplu1_arma = find(as<arma::vec>(xpos3_u) == -1.0);
        
        if (xplu0_arma.n_elem + xplu1_arma.n_elem != 2) xpos3_used = true;
      }
    }
    NumericVector ypos3_u = unique(ypos3);
    if (ypos3_u.length() > 1) {
      if (ypos3_u.length() > 2) {
        ypos3_used = true;
      } else {
        arma::uvec yplu0_arma = find(as<arma::vec>(ypos3_u) == 0.0);
        arma::uvec yplu1_arma = find(as<arma::vec>(ypos3_u) == -1.0);
        
        if (yplu0_arma.n_elem + yplu1_arma.n_elem != 2) ypos3_used = true;
      }
    }
    
    //Rcout << "jpf L" << endl;
    
    NumericVector censor1_u = unique(censor1);
    if (censor1_u.length() > 1) censor1_used = true;
    NumericVector censor2_u = unique(censor2);
    if (censor2_u.length() > 1) censor2_used = true;
    NumericVector censor3_u = unique(censor3);
    if (censor3_u.length() > 1) censor3_used = true;
    
    NumericVector sizea1_u = unique(sizea1);
    if (sizea1_u.length() > 1) sizea1_used = true;
    NumericVector sizea2_u = unique(sizea2);
    if (sizea2_u.length() > 1) sizea2_used = true;
    NumericVector sizea3_u = unique(sizea3);
    if (sizea3_u.length() > 1) sizea3_used = true;
    
    NumericVector sizeb1_u = unique(sizeb1);
    if (sizeb1_u.length() > 1) sizeb1_used = true;
    NumericVector sizeb2_u = unique(sizeb2);
    if (sizeb2_u.length() > 1) sizeb2_used = true;
    NumericVector sizeb3_u = unique(sizeb3);
    if (sizeb3_u.length() > 1) sizeb3_used = true;
    
    NumericVector sizec1_u = unique(sizec1);
    if (sizec1_u.length() > 1) sizec1_used = true;
    NumericVector sizec2_u = unique(sizec2);
    if (sizec2_u.length() > 1) sizec2_used = true;
    NumericVector sizec3_u = unique(sizec3);
    if (sizec3_u.length() > 1) sizec3_used = true;
    
    NumericVector size1added_u = unique(sizeadded1);
    if (size1added_u.length() > 1) size1added_used = true;
    NumericVector size2added_u = unique(sizeadded2);
    if (size2added_u.length() > 1) size2added_used = true;
    NumericVector size3added_u = unique(sizeadded3);
    if (size3added_u.length() > 1) size3added_used = true;
    
    NumericVector repstra1_u = unique(repstra1);
    if (repstra1_u.length() > 1) repstra1_used = true;
    NumericVector repstra2_u = unique(repstra2);
    if (repstra2_u.length() > 1) repstra2_used = true;
    NumericVector repstra3_u = unique(repstra3);
    if (repstra3_u.length() > 1) repstra3_used = true;
    
    NumericVector repstrb1_u = unique(repstrb1);
    if (repstrb1_u.length() > 1) repstrb1_used = true;
    NumericVector repstrb2_u = unique(repstrb2);
    if (repstrb2_u.length() > 1) repstrb2_used = true;
    NumericVector repstrb3_u = unique(repstrb3);
    if (repstrb3_u.length() > 1) repstrb3_used = true;
    
    NumericVector repstr1added_u = unique(repstradded1);
    if (repstr1added_u.length() > 1) repstr1added_used = true;
    NumericVector repstr2added_u = unique(repstradded2);
    if (repstr2added_u.length() > 1) repstr2added_used = true;
    NumericVector repstr3added_u = unique(repstradded3);
    if (repstr3added_u.length() > 1) repstr3added_used = true;
    
    NumericVector feca1_u = unique(feca1);
    if (feca1_u.length() > 1) feca1_used = true;
    NumericVector feca2_u = unique(feca2);
    if (feca2_u.length() > 1) feca2_used = true;
    NumericVector feca3_u = unique(feca3);
    if (feca3_u.length() > 1) feca3_used = true;
    
    NumericVector fecb1_u = unique(fecb1);
    if (fecb1_u.length() > 1) fecb1_used = true;
    NumericVector fecb2_u = unique(fecb2);
    if (fecb2_u.length() > 1) fecb2_used = true;
    NumericVector fecb3_u = unique(fecb3);
    if (fecb3_u.length() > 1) fecb3_used = true;
    
    NumericVector fec1added_u = unique(fecadded1);
    if (fec1added_u.length() > 1) fec1added_used = true;
    NumericVector fec2added_u = unique(fecadded2);
    if (fec2added_u.length() > 1) fec2added_used = true;
    NumericVector fec3added_u = unique(fecadded3);
    if (fec3added_u.length() > 1) fec3added_used = true;
    
    //Rcout << "jpf M" << endl;
    
    NumericVector indcova1_u = unique(indcova1);
    IntegerVector indcova1_int_u = unique(indcova1_int);
    StringVector indcova1_str_u = unique(indcova1_str);
    if (indcova1_u.length() > 1 || indcova1_int_u.length() > 1 ||
        indcova1_str_u.length() > 1) indcova1_used = true;
    NumericVector indcova2_u = unique(indcova2);
    IntegerVector indcova2_int_u = unique(indcova2_int);
    StringVector indcova2_str_u = unique(indcova2_str);
    if (indcova2_u.length() > 1 || indcova2_int_u.length() > 1 ||
        indcova2_str_u.length() > 1) indcova2_used = true;
    NumericVector indcova3_u = unique(indcova3);
    IntegerVector indcova3_int_u = unique(indcova3_int);
    StringVector indcova3_str_u = unique(indcova3_str);
    if (indcova3_u.length() > 1 || indcova3_int_u.length() > 1 ||
        indcova3_str_u.length() > 1) indcova3_used = true;
    
    NumericVector indcovb1_u = unique(indcovb1);
    IntegerVector indcovb1_int_u = unique(indcovb1_int);
    StringVector indcovb1_str_u = unique(indcovb1_str);
    if (indcovb1_u.length() > 1 || indcovb1_int_u.length() > 1 ||
        indcovb1_str_u.length() > 1) indcovb1_used = true;
    NumericVector indcovb2_u = unique(indcovb2);
    IntegerVector indcovb2_int_u = unique(indcovb2_int);
    StringVector indcovb2_str_u = unique(indcovb2_str);
    if (indcovb2_u.length() > 1 || indcovb2_int_u.length() > 1 ||
        indcovb2_str_u.length() > 1) indcovb2_used = true;
    NumericVector indcovb3_u = unique(indcovb3);
    IntegerVector indcovb3_int_u = unique(indcovb3_int);
    StringVector indcovb3_str_u = unique(indcovb3_str);
    if (indcovb3_u.length() > 1 || indcovb3_int_u.length() > 1 ||
        indcovb3_str_u.length() > 1) indcovb3_used = true;
    
    NumericVector indcovc1_u = unique(indcovc1);
    IntegerVector indcovc1_int_u = unique(indcovc1_int);
    StringVector indcovc1_str_u = unique(indcovc1_str);
    if (indcovc1_u.length() > 1 || indcovc1_int_u.length() > 1 ||
        indcovc1_str_u.length() > 1) indcovc1_used = true;
    NumericVector indcovc2_u = unique(indcovc2);
    IntegerVector indcovc2_int_u = unique(indcovc2_int);
    StringVector indcovc2_str_u = unique(indcovc2_str);
    if (indcovc2_u.length() > 1 || indcovc2_int_u.length() > 1 ||
        indcovc2_str_u.length() > 1) indcovc2_used = true;
    NumericVector indcovc3_u = unique(indcovc3);
    IntegerVector indcovc3_int_u = unique(indcovc3_int);
    StringVector indcovc3_str_u = unique(indcovc3_str);
    if (indcovc3_u.length() > 1 || indcovc3_int_u.length() > 1 ||
        indcovc3_str_u.length() > 1) indcovc3_used = true;
    
    //Rcout << "jpf N" << endl;
    
    NumericVector juvgiven1_u = unique(juvgiven1);
    if (juvgiven1_u.length() > 1) juvgiven1_used = true;
    NumericVector juvgiven2_u = unique(juvgiven2);
    if (juvgiven2_u.length() > 1) juvgiven2_used = true;
    NumericVector juvgiven3_u = unique(juvgiven3);
    if (juvgiven3_u.length() > 1) juvgiven3_used = true;
    
    int total_added = xpos1_used + ypos1_used + xpos2_used + ypos2_used +
      xpos3_used + ypos3_used + censor1_used + censor2_used + censor3_used +
      sizea1_used + sizea2_used + sizea3_used + sizeb1_used + sizeb2_used +
      sizeb3_used + sizec1_used + sizec2_used + sizec3_used + size1added_used +
      size2added_used + size3added_used + repstra1_used + repstra2_used +
      repstra3_used + repstrb1_used + repstrb2_used + repstrb3_used +
      repstr1added_used + repstr2added_used + repstr3added_used + feca1_used +
      feca2_used + feca3_used + fecb1_used + fecb2_used + fecb3_used +
      fec1added_used + fec2added_used + fec3added_used + indcova1_used +
      indcova2_used + indcova3_used + indcovb1_used + indcovb2_used +
      indcovb3_used + indcovc1_used + indcovc2_used + indcovc3_used +
      juvgiven1_used + juvgiven2_used + juvgiven3_used;
    
    int list_length = 30 + total_added;
    
    List reduced (list_length);
    Rcpp::CharacterVector varnames (list_length);
    
    reduced(0) = rowid;
    varnames(0) = "rowid";
    
    if (popid_type > 0) {
      popid_int.attr("class") = popid_class;
      if (popid_type > 1) {
        popid_int.attr("levels") = popid_levels;
      }
      reduced(1) = popid_int;
      
    } else { 
      reduced(1) = popid;
    }
    varnames(1) = "popid";
    
    if (patchid_type > 0) {
      patchid_int.attr("class") = patchid_class;
      if (patchid_type > 1) {
        patchid_int.attr("levels") = patchid_levels;
      }
      reduced(2) = patchid_int;
      
    } else { 
      reduced(2) = patchid;
    }
    varnames(2) = "patchid";
    
    if (individ_type > 0) {
      individ_int.attr("class") = individ_class;
      if (individ_type > 1) {
        individ_int.attr("levels") = individ_levels;
      }
      reduced(3) = individ_int;
      
    } else { 
      reduced(3) = individ;
    }
    varnames(3) = "individ";
    
    reduced(4) = year2;
    varnames(4) = "year2";
    reduced(5) = firstseen;
    varnames(5) = "firstseen";
    reduced(6) = lastseen;
    varnames(6) = "lastseen";
    reduced(7) = obsage;
    varnames(7) = "obsage";
    reduced(8) = obslifespan;
    varnames(8) = "obslifespan";
    
    int red_tracker {9};
    
    if (xpos1_used) {
      reduced(red_tracker) = xpos1;
      varnames(red_tracker) = "xpos1";
      red_tracker++;
    }
    if (ypos1_used) {
      reduced(red_tracker) = ypos1;
      varnames(red_tracker) = "ypos1";
      red_tracker++;
    }
    if (sizea1_used) {
      reduced(red_tracker) = sizea1;
      varnames(red_tracker) = "sizea1";
      red_tracker++;
    }
    if (sizeb1_used) {
      reduced(red_tracker) = sizeb1;
      varnames(red_tracker) = "sizeb1";
      red_tracker++;
    }
    if (sizec1_used) {
      reduced(red_tracker) = sizec1;
      varnames(red_tracker) = "sizec1";
      red_tracker++;
    }
    if (size1added_used) {
      reduced(red_tracker) = sizeadded1;
      varnames(red_tracker) = "size1added";
      red_tracker++;
    }
    if (repstra1_used) {
      reduced(red_tracker) = repstra1;
      varnames(red_tracker) = "repstra1";
      red_tracker++;
    }
    if (repstrb1_used) {
      reduced(red_tracker) = repstrb1;
      varnames(red_tracker) = "repstrb1";
      red_tracker++;
    }
    if (repstr1added_used) {
      reduced(red_tracker) = repstradded1;
      varnames(red_tracker) = "repstr1added";
      red_tracker++;
    }
    if (feca1_used) {
      reduced(red_tracker) = feca1;
      varnames(red_tracker) = "feca1";
      red_tracker++;
    }
    if (fecb1_used) {
      reduced(red_tracker) = fecb1;
      varnames(red_tracker) = "fecb1";
      red_tracker++;
    }
    if (fec1added_used) {
      reduced(red_tracker) = fecadded1;
      varnames(red_tracker) = "fec1added";
      red_tracker++;
    }
    
    //Rcout << "jpf O" << endl;
    
    if (indcova1_used) {
      if (indcova_type == 3) {
        reduced(red_tracker) = indcova1;
      } else if (indcova_type > 0) {
        indcova1_int.attr("class") = indcova_class;
        if (indcova_type == 2) indcova1_int.attr("levels") = indcova_levels;
        reduced(red_tracker) = indcova1_int;
      } else {
        reduced(red_tracker) = indcova1_str;
      }
      varnames(red_tracker) = "indcova1";
      red_tracker++;
    }
    
    if (indcovb1_used) {
      if (indcovb_type == 3) {
        reduced(red_tracker) = indcovb1;
      } else if (indcovb_type > 0) {
        indcovb1_int.attr("class") = indcovb_class;
        if (indcovb_type == 2) indcovb1_int.attr("levels") = indcovb_levels;
        reduced(red_tracker) = indcovb1_int;
      } else {
        reduced(red_tracker) = indcovb1_str;
      }
      varnames(red_tracker) = "indcovb1";
      red_tracker++;
    }
    
    if (indcovc1_used) {
      if (indcovc_type == 3) {
        reduced(red_tracker) = indcovc1;
      } else if (indcovc_type > 0) {
        indcovc1_int.attr("class") = indcovc_class;
        if (indcovc_type == 2) indcovc1_int.attr("levels") = indcovc_levels;
        reduced(red_tracker) = indcovc1_int;
      } else {
        reduced(red_tracker) = indcovc1_str;
      }
      varnames(red_tracker) = "indcovc1";
      red_tracker++;
    }
    
    //Rcout << "jpf P" << endl;
    
    if (censor1_used) {
      reduced(red_tracker) = censor1;
      varnames(red_tracker) = "censor1";
      red_tracker++;
    }
    if (juvgiven1_used) {
      reduced(red_tracker) = juvgiven1;
      varnames(red_tracker) = "juvgiven1";
      red_tracker++;
    }
    
    reduced(red_tracker) = obsstatus1;
    varnames(red_tracker) = "obsstatus1";
    red_tracker++;
    reduced(red_tracker) = repstatus1;
    varnames(red_tracker) = "repstatus1";
    red_tracker++;
    
    reduced(red_tracker) = fecstatus1;
    varnames(red_tracker) = "fecstatus1";
    red_tracker++;
    reduced(red_tracker) = matstat1;
    varnames(red_tracker) = "matstatus1";
    red_tracker++;
    
    reduced(red_tracker) = alive1;
    varnames(red_tracker) = "alive1";
    red_tracker++;
    reduced(red_tracker) = stage1;
    varnames(red_tracker) = "stage1";
    red_tracker++;
    reduced(red_tracker) = stage1num;
    varnames(red_tracker) = "stage1index";
    red_tracker++;
    
    if (xpos2_used) {
      reduced(red_tracker) = xpos2;
      varnames(red_tracker) = "xpos2";
      red_tracker++;
    }
    if (ypos2_used) {
      reduced(red_tracker) = ypos2;
      varnames(red_tracker) = "ypos2";
      red_tracker++;
    }
    if (sizea2_used) {
      reduced(red_tracker) = sizea2;
      varnames(red_tracker) = "sizea2";
      red_tracker++;
    }
    if (sizeb2_used) {
      reduced(red_tracker) = sizeb2;
      varnames(red_tracker) = "sizeb2";
      red_tracker++;
    }
    if (sizec2_used) {
      reduced(red_tracker) = sizec2;
      varnames(red_tracker) = "sizec2";
      red_tracker++;
    }
    if (size2added_used) {
      reduced(red_tracker) = sizeadded2;
      varnames(red_tracker) = "size2added";
      red_tracker++;
    }
    if (repstra2_used) {
      reduced(red_tracker) = repstra2;
      varnames(red_tracker) = "repstra2";
      red_tracker++;
    }
    if (repstrb2_used) {
      reduced(red_tracker) = repstrb2;
      varnames(red_tracker) = "repstrb2";
      red_tracker++;
    }
    if (repstr2added_used) {
      reduced(red_tracker) = repstradded2;
      varnames(red_tracker) = "repstr2added";
      red_tracker++;
    }
    if (feca2_used) {
      reduced(red_tracker) = feca2;
      varnames(red_tracker) = "feca2";
      red_tracker++;
    }
    if (fecb2_used) {
      reduced(red_tracker) = fecb2;
      varnames(red_tracker) = "fecb2";
      red_tracker++;
    }
    if (fec2added_used) {
      reduced(red_tracker) = fecadded2;
      varnames(red_tracker) = "fec2added";
      red_tracker++;
    }
    
    //Rcout << "jpf Q" << endl;
    
    if (indcova2_used) {
      if (indcova_type == 3) {
        reduced(red_tracker) = indcova2;
      } else if (indcova_type > 0) {
        indcova2_int.attr("class") = indcova_class;
        if (indcova_type == 2) indcova2_int.attr("levels") = indcova_levels;
        reduced(red_tracker) = indcova2_int;
      } else {
        reduced(red_tracker) = indcova2_str;
      }
      varnames(red_tracker) = "indcova2";
      red_tracker++;
    }
    
    if (indcovb2_used) {
      if (indcovb_type == 3) {
        reduced(red_tracker) = indcovb2;
      } else if (indcovb_type > 0) {
        indcovb2_int.attr("class") = indcovb_class;
        if (indcovb_type == 2) indcovb2_int.attr("levels") = indcovb_levels;
        reduced(red_tracker) = indcovb2_int;
      } else {
        reduced(red_tracker) = indcovb2_str;
      }
      varnames(red_tracker) = "indcovb2";
      red_tracker++;
    }
    
    if (indcovc2_used) {
      if (indcovc_type == 3) {
        reduced(red_tracker) = indcovc2;
      } else if (indcovc_type > 0) {
        indcovc2_int.attr("class") = indcovc_class;
        if (indcovc_type == 2) indcovc2_int.attr("levels") = indcovc_levels;
        reduced(red_tracker) = indcovc2_int;
      } else {
        reduced(red_tracker) = indcovc2_str;
      }
      varnames(red_tracker) = "indcovc2";
      red_tracker++;
    }
    
    //Rcout << "jpf R" << endl;
    
    if (censor2_used) {
      reduced(red_tracker) = censor2;
      varnames(red_tracker) = "censor2";
      red_tracker++;
    }
    if (juvgiven2_used) {
      reduced(red_tracker) = juvgiven2;
      varnames(red_tracker) = "juvgiven2";
      red_tracker++;
    }
    
    reduced(red_tracker) = obsstatus2;
    varnames(red_tracker) = "obsstatus2";
    red_tracker++;
    reduced(red_tracker) = repstatus2;
    varnames(red_tracker) = "repstatus2";
    red_tracker++;
    
    reduced(red_tracker) = fecstatus2;
    varnames(red_tracker) = "fecstatus2";
    red_tracker++;
    reduced(red_tracker) = matstat2;
    varnames(red_tracker) = "matstatus2";
    red_tracker++;
    
    reduced(red_tracker) = alive2;
    varnames(red_tracker) = "alive2";
    red_tracker++;
    reduced(red_tracker) = stage2;
    varnames(red_tracker) = "stage2";
    red_tracker++;
    reduced(red_tracker) = stage2num;
    varnames(red_tracker) = "stage2index";
    red_tracker++;
    
    if (xpos3_used) {
      reduced(red_tracker) = xpos3;
      varnames(red_tracker) = "xpos3";
      red_tracker++;
    }
    if (ypos3_used) {
      reduced(red_tracker) = ypos3;
      varnames(red_tracker) = "ypos3";
      red_tracker++;
    }
    if (sizea3_used) {
      reduced(red_tracker) = sizea3;
      varnames(red_tracker) = "sizea3";
      red_tracker++;
    }
    if (sizeb3_used) {
      reduced(red_tracker) = sizeb3;
      varnames(red_tracker) = "sizeb3";
      red_tracker++;
    }
    if (sizec3_used) {
      reduced(red_tracker) = sizec3;
      varnames(red_tracker) = "sizec3";
      red_tracker++;
    }
    if (size3added_used) {
      reduced(red_tracker) = sizeadded3;
      varnames(red_tracker) = "size3added";
      red_tracker++;
    }
    if (repstra3_used) {
      reduced(red_tracker) = repstra3;
      varnames(red_tracker) = "repstra3";
      red_tracker++;
    }
    if (repstrb3_used) {
      reduced(red_tracker) = repstrb3;
      varnames(red_tracker) = "repstrb3";
      red_tracker++;
    }
    if (repstr3added_used) {
      reduced(red_tracker) = repstradded3;
      varnames(red_tracker) = "repstr3added";
      red_tracker++;
    }
    if (feca3_used) {
      reduced(red_tracker) = feca3;
      varnames(red_tracker) = "feca3";
      red_tracker++;
    }
    if (fecb3_used) {
      reduced(red_tracker) = fecb3;
      varnames(red_tracker) = "fecb3";
      red_tracker++;
    }
    if (fec3added_used) {
      reduced(red_tracker) = fecadded3;
      varnames(red_tracker) = "fec3added";
      red_tracker++;
    }
    
    //Rcout << "jpf S" << endl;
    
    if (indcova3_used) {
      if (indcova_type == 3) {
        reduced(red_tracker) = indcova3;
      } else if (indcova_type > 0) {
        indcova3_int.attr("class") = indcova_class;
        if (indcova_type == 2) indcova3_int.attr("levels") = indcova_levels;
        reduced(red_tracker) = indcova3_int;
      } else {
        reduced(red_tracker) = indcova3_str;
      }
      varnames(red_tracker) = "indcova3";
      red_tracker++;
    }
    
    if (indcovb3_used) {
      if (indcovb_type == 3) {
        reduced(red_tracker) = indcovb3;
      } else if (indcovb_type > 0) {
        indcovb3_int.attr("class") = indcovb_class;
        if (indcovb_type == 2) indcovb3_int.attr("levels") = indcovb_levels;
        reduced(red_tracker) = indcovb3_int;
      } else {
        reduced(red_tracker) = indcovb3_str;
      }
      varnames(red_tracker) = "indcovb3";
      red_tracker++;
    }
    
    if (indcovc3_used) {
      if (indcovc_type == 3) {
        reduced(red_tracker) = indcovc3;
      } else if (indcovc_type > 0) {
        indcovc3_int.attr("class") = indcovc_class;
        if (indcovc_type == 2) indcovc3_int.attr("levels") = indcovc_levels;
        reduced(red_tracker) = indcovc3_int;
      } else {
        reduced(red_tracker) = indcovc3_str;
      }
      varnames(red_tracker) = "indcovc3";
      red_tracker++;
    }
    
    //Rcout << "jpf T" << endl;
    
    if (censor3_used) {
      reduced(red_tracker) = censor3;
      varnames(red_tracker) = "censor3";
      red_tracker++;
    }
    if (juvgiven3_used) {
      reduced(red_tracker) = juvgiven3;
      varnames(red_tracker) = "juvgiven3";
      red_tracker++;
    }
    
    reduced(red_tracker) = obsstatus3;
    varnames(red_tracker) = "obsstatus3";
    red_tracker++;
    reduced(red_tracker) = repstatus3;
    varnames(red_tracker) = "repstatus3";
      
    red_tracker++;
    reduced(red_tracker) = fecstatus3;
    varnames(red_tracker) = "fecstatus3";
    red_tracker++;
    reduced(red_tracker) = matstat3;
    varnames(red_tracker) = "matstatus3";
    red_tracker++;
    reduced(red_tracker) = alive3;
    varnames(red_tracker) = "alive3";
    red_tracker++;
    reduced(red_tracker) = stage3;
    varnames(red_tracker) = "stage3";
    red_tracker++;
    reduced(red_tracker) = stage3num;
    varnames(red_tracker) = "stage3index";
    
    reduced.attr("names") = varnames;
    reduced.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, alive3.length());
    StringVector needed_classes {"data.frame", "hfvdata"};
    reduced.attr("class") = needed_classes;
    
    final_output = reduced;
    
  } else {
    Rcpp::List output_longlist(81);
    
    output_longlist(0) = rowid;
    
    output_longlist(4) = Rcpp::IntegerVector(year2.begin(), year2.end());
    output_longlist(5) = Rcpp::IntegerVector(firstseen.begin(), firstseen.end());
    output_longlist(6) = Rcpp::IntegerVector(lastseen.begin(), lastseen.end());
    output_longlist(7) = Rcpp::IntegerVector(obsage.begin(), obsage.end());
    output_longlist(8) = Rcpp::IntegerVector(obslifespan.begin(), obslifespan.end());
    
    output_longlist(9) = xpos1;
    output_longlist(10) = ypos1;
    output_longlist(11) = sizea1;
    output_longlist(12) = sizeb1;
    output_longlist(13) = sizec1;
    output_longlist(14) = sizeadded1;
    output_longlist(15) = repstra1;
    output_longlist(16) = repstrb1;
    output_longlist(17) = repstradded1;
    output_longlist(18) = feca1;
    output_longlist(19) = fecb1;
    output_longlist(20) = fecadded1;
    
    output_longlist(24) = censor1;
    output_longlist(25) = juvgiven1;
    
    output_longlist(26) = obsstatus1;
    output_longlist(27) = repstatus1;
    output_longlist(28) = fecstatus1;
    output_longlist(29) = matstat1;
    output_longlist(30) = alive1;
    output_longlist(31) = stage1;
    output_longlist(32) = stage1num;
    
    output_longlist(33) = xpos2;
    output_longlist(34) = ypos2;
    output_longlist(35) = sizea2;
    output_longlist(36) = sizeb2;
    output_longlist(37) = sizec2;
    output_longlist(38) = sizeadded2;
    output_longlist(39) = repstra2;
    output_longlist(40) = repstrb2;
    output_longlist(41) = repstradded2;
    output_longlist(42) = feca2;
    output_longlist(43) = fecb2;
    output_longlist(44) = fecadded2;
    output_longlist(48) = censor2;
    output_longlist(49) = juvgiven2;
    
    output_longlist(50) = obsstatus2;
    output_longlist(51) = repstatus2;
    output_longlist(52) = fecstatus2;
    output_longlist(53) = matstat2;
    output_longlist(54) = alive2;
    output_longlist(55) = stage2;
    output_longlist(56) = stage2num;
  
    output_longlist(57) = xpos3;
    output_longlist(58) = ypos3;
    output_longlist(59) = sizea3;
    output_longlist(60) = sizeb3;
    output_longlist(61) = sizec3;
    output_longlist(62) = sizeadded3;
    output_longlist(63) = repstra3;
    output_longlist(64) = repstrb3;
    output_longlist(65) = repstradded3;
    output_longlist(66) = feca3;
    output_longlist(67) = fecb3;
    output_longlist(68) = fecadded3;
    output_longlist(72) = censor3;
    output_longlist(73) = juvgiven3;
    
    output_longlist(74) = obsstatus3;
    output_longlist(75) = repstatus3;
    output_longlist(76) = fecstatus3;
    output_longlist(77) = matstat3;
    output_longlist(78) = alive3;
    output_longlist(79) = stage3;
    output_longlist(80) = stage3num;
    
    if (popid_type > 0) {
      if (popid_type > 1) popid_int.attr("levels") = popid_levels;
      output_longlist(1) = popid_int;
      
    } else {
      output_longlist(1) = popid;
    }
    
    if (patchid_type > 0) { 
      if (patchid_type > 1) patchid_int.attr("levels") = patchid_levels;
      output_longlist(2) = patchid_int;
      
    } else {
      output_longlist(2) = patchid;
    }
    
    if (individ_type > 0) { 
      if (individ_type > 1) individ_int.attr("levels") = individ_levels;
      output_longlist(3) = individ_int;
      
    } else {
      output_longlist(3) = individ;
    }
    
    //Rcout << "jpf U" << endl;
    
    if (indcova_type == 1 || indcova_type == 2) {
      output_longlist(21) = indcova1_int;
      output_longlist(45) = indcova2_int;
      output_longlist(69) = indcova3_int;
    } else if (indcova_type == 3) { 
      output_longlist(21) = indcova1;
      output_longlist(45) = indcova2;
      output_longlist(69) = indcova3;
    } else {
      output_longlist(21) = indcova1_str;
      output_longlist(45) = indcova2_str;
      output_longlist(69) = indcova3_str;
    }
    
    if (indcovb_type > 0 || indcovb_type == 2) { 
      output_longlist(22) = indcovb1_int;
      output_longlist(46) = indcovb2_int;
      output_longlist(70) = indcovb3_int;
    } else if (indcovb_type == 3)  { 
      output_longlist(22) = indcovb1;
      output_longlist(46) = indcovb2;
      output_longlist(70) = indcovb3;
    } else {
      output_longlist(22) = indcovb1_str;
      output_longlist(46) = indcovb2_str;
      output_longlist(70) = indcovb3_str;
    }
    
    if (indcovc_type > 0 || indcovc_type == 2) { 
      output_longlist(23) = indcovc1_int;
      output_longlist(47) = indcovc2_int;
      output_longlist(71) = indcovc3_int;
    } else if (indcovc_type == 3)  { 
      output_longlist(23) = indcovc1;
      output_longlist(47) = indcovc2;
      output_longlist(71) = indcovc3;
    } else {
      output_longlist(23) = indcovc1_str;
      output_longlist(47) = indcovc2_str;
      output_longlist(71) = indcovc3_str;
    }
    
    //Rcout << "jpf V" << endl;
    
    Rcpp::CharacterVector varnames {"rowid", "popid", "patchid", "individ",
      "year2", "firstseen", "lastseen", "obsage", "obslifespan","xpos1", "ypos1",
      "sizea1", "sizeb1", "sizec1", "size1added", "repstra1", "repstrb1",
      "repstr1added", "feca1", "fecb1", "fec1added", "indcova1", "indcovb1",
      "indcovc1", "censor1", "juvgiven1", "obsstatus1", "repstatus1",
      "fecstatus1", "matstatus1", "alive1", "stage1", "stage1index",
      "xpos2", "ypos2", "sizea2", "sizeb2", "sizec2", "size2added", "repstra2",
      "repstrb2", "repstr2added", "feca2", "fecb2", "fec2added", "indcova2",
      "indcovb2", "indcovc2", "censor2", "juvgiven2", "obsstatus2", "repstatus2",
      "fecstatus2", "matstatus2", "alive2", "stage2", "stage2index",
      "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "size3added", "repstra3",
      "repstrb3", "repstr3added", "feca3", "fecb3", "fec3added", "indcova3",
      "indcovb3", "indcovc3", "censor3", "juvgiven3", "obsstatus3", "repstatus3",
      "fecstatus3", "matstatus3", "alive3", "stage3", "stage3index"};
    
    output_longlist.attr("names") = varnames;
    output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, alive3.length());
    StringVector needed_classes {"data.frame", "hfvdata"};
    output_longlist.attr("class") = needed_classes;
    
    final_output = output_longlist;
  }
  
  if (!retain_alive0) { 
    NumericVector retained_bit {1.0};
    StringVector alive_var_used = {"alive2"};
    List reduced_output = LefkoUtils::df_subset(final_output, retained_bit, false, true,
      false, false, true, as<RObject>(alive_var_used));
      
    return reduced_output;
    
  } else { 
    return final_output;
  }
}

//' Estimate Radial Density in Cartesian Space
//' 
//' Function \code{.density3()} estimates radial density on the basis of
//' Cartesian coordinates and spacing information supplied as input. It is used
//' internally by \code{\link{historicalize3}()} and
//' \code{\link{verticalize3}()}.
//' 
//' @name .density3
//' 
//' @param data Demographic dataset in historical vertical format.
//' @param xcol Number of column in \code{data} corresponding to x position.
//' @param ycol Number of column in \code{data} corresponding to y position.
//' @param yearcol Number of column in \code{data} corresponding to occasion
//' \emph{t}.
//' @param spacing Resolution of density estimation, as a scalar numeric.
//' 
//' @return This function returns a vector counting the number of individuals
//' within the specified distance of each individual in the historically
//' formatted vertical dataset.
//' 
//' @section Notes:
//' The process used to estimate density is one in which the distances between
//' all pairs of individuals are calculated via the Pythagorean theorem, and
//' then individual density equals the number of these individuals with
//' distances within the number input as \code{spacing}, respectively for each
//' individual.
//' 
//' This function assumes that all individuals are alive in time \emph{t}, and
//' so data should be filtered appropriately beforehand. Any rows with NA in X
//' or Y will not be counted, and density is estimated specific to time \emph{t}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.density3)]]
Rcpp::NumericVector density3(Rcpp::DataFrame data, int xcol, int ycol,
  int yearcol, double spacing) {
  
  int data_size = data.length();
  int data_n = data.nrows();
  
  if (xcol < 0 || ycol < 0 || xcol > data_size || ycol > data_size) {
    throw Rcpp::exception("Input column numbers for X and/or Y are outside the range of the dataset", 
      false);
  }
  
  int xcol_true = xcol - 1;
  int ycol_true = ycol - 1;
  
  NumericVector Xdata = as<NumericVector>(data[xcol_true]);
  NumericVector Ydata = as<NumericVector>(data[ycol_true]);
  NumericVector yeardata = as<NumericVector>(data[(yearcol - 1)]);
  
  double ref_x {0.0};
  double ref_y {0.0};
  double est_a {0.0};
  double est_b {0.0};
  double est_c {0.0};
  int counted_n {0};
  
  arma::vec density(data_n);
  density.zeros();
  
  for (int i = 0; i < data_n; i++) {
    ref_x = Xdata(i);
    ref_y = Ydata(i);
    
    if (!NumericVector::is_na(Xdata(i)) && !NumericVector::is_na(Ydata(i))) {
      for (int j = 0; j < data_n; j++) {
        if (!NumericVector::is_na(Xdata(j)) && !NumericVector::is_na(Ydata(j)) && 
            yeardata(i) == yeardata(j)) {
          
          est_a = Xdata(j) - ref_x;
          est_b = Ydata(j) - ref_y;
          
          est_c = (est_a * est_a) + (est_b * est_b);
          est_c = sqrt(est_c);
          
          if (est_c < spacing) counted_n++;
        }
      }
    } else {
      counted_n = 1;
    }
    
    density(i) = counted_n - 1;
    counted_n = 0;
  }
  
  return Rcpp::NumericVector(density.begin(), density.end());
}

//' Bootstrap Standardized hfv_data Datasets
//' 
//' Function \code{bootstrap3()} takes already standardized \code{hfvdata}
//' datasets and bootstraps them by individual identity, or by row.
//' 
//' @name bootstrap3
//' 
//' @param data A data frame of class \code{hfvdata}.
//' @param by_pop A logical value indicating whether to sample the data frame
//' by population. If \code{TRUE}, then the number of individuals sampled for
//' each population will be set to the respective population's actual number of
//' individuals; otherwise, population identity is ignored. Defaults to
//' \code{TRUE}.
//' @param by_patch A logical value indicating whether to sample the data frame
//' by patch. If \code{TRUE}, then the number of individuals sampled for each
//' patch will be set to the respective patch's actual number of individuals;
//' otherwise, patch identity is ignored. Defaults to \code{TRUE}.
//' @param by_indiv A logical value indicating whether to sample the data frame
//' by individual identity, or by row. If \code{TRUE}, then samples by
//' individual identity. Defaults to \code{TRUE}.
//' @param prop_size A logical value indicating whether to keep the proportions
//' of individuals (if \code{by_indiv = TRUE}) or of rows (if \code{by_indiv =
//' FALSE}) in each bootstrapped dataset to the same proportions across
//' populations (if \code{by_pop = TRUE}, and patches (if
//' \code{by_patch = TRUE}, as in the original dataset. If \code{FALSE}, then
//' allows the specific proportions to be set by argument \code{max_limit}.
//' Defaults to \code{TRUE}.
//' @param max_limit Sets the sample size to pull from the original data frame,
//' if \code{prop_size = FALSE}. Defaults to the size of the original dataset if
//' \code{prop_size = TRUE}, and to 100 if if \code{prop_size = FALSE}. Can also
//' be input as an integer vector giving the number of samples to take by
//' population (if \code{by_pop = TRUE}), patch (if \code{by_patch = TRUE}), or
//' population-patch (if \code{by_pop = TRUE} and \code{by_patch = TRUE}).
//' @param reps The number of bootstrap replicates to produce. Defaults to
//' \code{100}.
//' @param popcol A string denoting the variable name coding for population
//' identity in the data frame. Defaults to \code{"popid"}.
//' @param patchcol A string denoting the variable name coding for patch
//' identity in the data frame. Defaults to \code{"patchid"}.
//' @param indivcol A string denoting the variable name coding for individual
//' identity in the data frame. Defaults to \code{"individ"}.
//' 
//' @return A list of class \code{hfvlist}, which is composed of data frames of
//' class \code{hfvdata}.
//' 
//' @examples
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
//' lathboot <- bootstrap3(lathvert)
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
//' lathproj <- projection3(ehrlen3, nreps = 5, integeronly = TRUE,
//'   stochastic = TRUE)
//' 
//' @export bootstrap3
// [[Rcpp::export(bootstrap3)]]
Rcpp::List bootstrap3(RObject data, Nullable<RObject> by_pop = R_NilValue,
  Nullable<RObject> by_patch = R_NilValue, Nullable<RObject> by_indiv = R_NilValue,
  Nullable<RObject> prop_size = R_NilValue, Nullable<RObject> max_limit = R_NilValue,
  Nullable<RObject> reps = R_NilValue, Nullable<RObject> popcol = R_NilValue,
  Nullable<RObject> patchcol = R_NilValue, Nullable<RObject> indivcol = R_NilValue ) {
  
  DataFrame true_data;
  bool by_pop_bool {true};
  bool by_patch_bool {true};
  bool by_indiv_bool {true};
  bool prop_size_bool {true};
  bool max_limit_bool {false};
  int reps_true {100};
  int default_sample_set {100};
  
  StringVector popcol_str;
  StringVector patchcol_str;
  StringVector indivcol_str;
  int popcol_int {-1};
  int patchcol_int {-1};
  int indivcol_int {-1};
  bool found_pop_col {false};
  bool found_patch_col {false};
  bool found_indiv_col {false};
  
  int total_pops {1};
  int total_patches {1};
  int total_poppatches {1};
  int total_indivs {0};
  int total_rows {0};
  
  StringVector pops_allrows_str;
  StringVector patches_allrows_str;
  StringVector poppatches_allrows_str;
  StringVector individuals_allrows_str;
  StringVector unique_pops_str;
  StringVector unique_patches_str;
  StringVector unique_poppatches_str;
  StringVector unique_individuals_str;
  IntegerVector max_limit_int;
  
  CharacterVector df_class_lel = {"data.frame", "hfvdata"};
  
  if (by_pop.isNotNull()) {
    RObject by_pop_RO (by_pop);
    if (is<LogicalVector>(by_pop_RO)) {
      LogicalVector by_pop_log = as<LogicalVector>(by_pop_RO);
      by_pop_bool = by_pop_log(0);
    } else LefkoUtils::pop_error("by_pop", "a logical value", "", 1);
  }
  
  if (by_patch.isNotNull()) {
    RObject by_patch_RO (by_patch);
    if (is<LogicalVector>(by_patch_RO)) {
      LogicalVector by_patch_log = as<LogicalVector>(by_patch_RO);
      by_patch_bool = by_patch_log(0);
    } else LefkoUtils::pop_error("by_patch", "a logical value", "", 1);
  }
  
  if (by_indiv.isNotNull()) {
    RObject by_indiv_RO (by_indiv);
    if (is<LogicalVector>(by_indiv_RO)) {
      LogicalVector by_indiv_log = as<LogicalVector>(by_indiv_RO);
      by_indiv_bool = by_indiv_log(0);
    } else LefkoUtils::pop_error("by_indiv", "a logical value", "", 1);
  }
  
  if (prop_size.isNotNull()) {
    RObject prop_size_RO (prop_size);
    if (is<LogicalVector>(prop_size_RO)) {
      LogicalVector prop_size_log = as<LogicalVector>(prop_size_RO);
      prop_size_bool = prop_size_log(0);
    } else LefkoUtils::pop_error("prop_size", "a logical value", "", 1);
  }
  
  if (max_limit.isNotNull()) {
    RObject max_limit_RO (max_limit);
    if (is<IntegerVector>(max_limit_RO) || is<NumericVector>(max_limit_RO)) {
      max_limit_int = as<IntegerVector>(max_limit_RO);
      max_limit_bool = true;
      //sample_limit = max_limit_int(0);
    } else LefkoUtils::pop_error("max_limit", "an integer", "", 1);
  }
  
  if (reps.isNotNull()) {
    RObject reps_RO (reps);
    if (is<IntegerVector>(reps_RO) || is<NumericVector>(reps_RO)) {
      IntegerVector reps_int = as<IntegerVector>(reps_RO);
      reps_true = reps_int(0);
    } else LefkoUtils::pop_error("reps", "an integer", "", 1);
  }
  
  List hfv_list (reps_true);
  
  if (is<DataFrame>(data)) {
    true_data = as<DataFrame>(data);
    total_rows = static_cast<int>(true_data.nrows());
    
    CharacterVector data_names = true_data.attr("names");
    int df_var_num = static_cast<int>(data_names.length());
    
    if (true_data.hasAttribute("class")) {
      CharacterVector true_data_classes = wrap(true_data.attr("class"));
      bool found_hfv {false};
      for (int i = 0; i < true_data_classes.length(); i++) {
        if (stringcompare_hard(String(true_data_classes(i)), "hfvdata")) found_hfv = true;
      }
      
      if (!found_hfv) {
        LefkoUtils::pop_error("data", "an hfvdata object", "", 1);
      }
      
      if (popcol.isNotNull()) {
        RObject popcol_RO (popcol);
        if (is<StringVector>(popcol_RO)) {
          popcol_str = as<StringVector>(popcol_RO);
          
          for (int i = 0; i < data_names.length(); i++) {
            if (stringcompare_hard(as<std::string>(data_names(i)), String(popcol_str(0)))) {
              popcol_int = i;
              found_pop_col = true;
            }
          }
        } else if (is<IntegerVector>(popcol_RO) || is<NumericVector>(popcol_RO)) {
          IntegerVector popcol_intvec = as<IntegerVector>(popcol_RO);
          popcol_int = popcol_intvec(0);
          if (popcol_int >= 0 && popcol_int < df_var_num) found_pop_col = true;
        } else LefkoUtils::pop_error("popcol", "a string value", "", 1);
        
        if (is<IntegerVector>(true_data(popcol_int))) {
          IntegerVector pops_alldata = as<IntegerVector>(true_data(popcol_int));
          pops_allrows_str = StringVector(pops_alldata);
          
        } else if (is<NumericVector>(true_data(popcol_int))) {
          NumericVector pops_alldata = as<NumericVector>(true_data(popcol_int));
          pops_allrows_str = StringVector(pops_alldata);
          
        } else if (is<StringVector>(true_data(popcol_int))) {
          StringVector pops_alldata = as<StringVector>(true_data(popcol_int));
          pops_allrows_str = StringVector(pops_alldata);
          
        } else throw Rcpp::exception("Population identity is not in a recognized format.", false);
        
        unique_pops_str = sort_unique(pops_allrows_str);
        total_pops = static_cast<int>(unique_pops_str.length());
        
      } else if (by_pop_bool) {
        StringVector new_popcol_str = {"popid"};
        popcol_str = new_popcol_str;
        
        for (int i = 0; i < data_names.length(); i++) {
          if (stringcompare_hard(as<std::string>(data_names(i)), String(popcol_str(0)))) {
            popcol_int = i;
            found_pop_col = true;
          }
        }
        
        if (is<IntegerVector>(true_data(popcol_int))) {
          IntegerVector pops_alldata = as<IntegerVector>(true_data(popcol_int));
          pops_allrows_str = StringVector(pops_alldata);
          
        } else if (is<NumericVector>(true_data(popcol_int))) {
          NumericVector pops_alldata = as<NumericVector>(true_data(popcol_int));
          pops_allrows_str = StringVector(pops_alldata);
          
        } else if (is<StringVector>(true_data(popcol_int))) {
          StringVector pops_alldata = as<StringVector>(true_data(popcol_int));
          pops_allrows_str = StringVector(pops_alldata);
          
        } else throw Rcpp::exception("Population identity is not in a recognized format.", false);
        
        unique_pops_str = sort_unique(pops_allrows_str);
        total_pops = static_cast<int>(unique_pops_str.length());
      } else {
        StringVector new_popcol_str = {"popid"};
        popcol_str = new_popcol_str;
        
        for (int i = 0; i < data_names.length(); i++) {
          if (stringcompare_hard(as<std::string>(data_names(i)), String(popcol_str(0)))) {
            popcol_int = i;
            found_pop_col = true;
          }
        }
        
        StringVector pops_allrows_str_temp (total_rows);
        for (int i = 0; i < total_rows; i++) pops_allrows_str_temp(i) = "0";
        pops_allrows_str = pops_allrows_str_temp;
        unique_pops_str = sort_unique(pops_allrows_str);
      }
      
      if (patchcol.isNotNull()) {
        RObject patchcol_RO (patchcol);
        if (is<StringVector>(patchcol_RO)) {
          patchcol_str = as<StringVector>(patchcol_RO);
          
          for (int i = 0; i < data_names.length(); i++) {
            if (stringcompare_hard(as<std::string>(data_names(i)), String(patchcol_str(0)))) {
              patchcol_int = i;
              found_patch_col = true;
            }
          }
        } else if (is<IntegerVector>(patchcol_RO) || is<NumericVector>(patchcol_RO)) {
          IntegerVector patchcol_intvec = as<IntegerVector>(patchcol_RO);
          patchcol_int = patchcol_intvec(0);
          if (patchcol_int >= 0 && patchcol_int < df_var_num) found_patch_col = true;
        } else LefkoUtils::pop_error("patchcol", "a string value", "", 1);
        
        if (is<IntegerVector>(true_data(patchcol_int))) {
          IntegerVector patches_alldata = as<IntegerVector>(true_data(patchcol_int));
          patches_allrows_str = StringVector(patches_alldata);
          
        } else if (is<NumericVector>(true_data(patchcol_int))) {
          NumericVector patches_alldata = as<NumericVector>(true_data(patchcol_int));
          patches_allrows_str = StringVector(patches_alldata);
          
        } else if (is<StringVector>(true_data(patchcol_int))) {
          StringVector patches_alldata = as<StringVector>(true_data(patchcol_int));
          patches_allrows_str = StringVector(patches_alldata);
          
        } else throw Rcpp::exception("Patch identity is not in a recognized format.", false);
        
        unique_patches_str = sort_unique(patches_allrows_str); // Should be redone as pop-patches
        total_patches = static_cast<int>(unique_patches_str.length()); // Should be redone as pop-patches
        
      } else if (by_patch_bool) {
        StringVector new_patchcol_str = {"patchid"};
        patchcol_str = new_patchcol_str;
        
        for (int i = 0; i < data_names.length(); i++) {
          if (stringcompare_hard(as<std::string>(data_names(i)), String(patchcol_str(0)))) {
            patchcol_int = i;
            found_patch_col = true;
          }
        }
        
        if (is<IntegerVector>(true_data(patchcol_int))) {
          IntegerVector patches_alldata = as<IntegerVector>(true_data(patchcol_int));
          patches_allrows_str = StringVector(patches_alldata);
          
        } else if (is<NumericVector>(true_data(patchcol_int))) {
          NumericVector patches_alldata = as<NumericVector>(true_data(patchcol_int));
          patches_allrows_str = StringVector(patches_alldata);
          
        } else if (is<StringVector>(true_data(patchcol_int))) {
          StringVector patches_alldata = as<StringVector>(true_data(patchcol_int));
          patches_allrows_str = StringVector(patches_alldata);
          
        } else throw Rcpp::exception("Patch identity is not in a recognized format.", false);
        
        unique_patches_str = sort_unique(patches_allrows_str); // Should be redone as pop-patches
        total_patches = static_cast<int>(unique_patches_str.length()); // Should be redone as pop-patches
      } else {
        StringVector new_patchcol_str = {"patchid"};
        patchcol_str = new_patchcol_str;
        
        for (int i = 0; i < data_names.length(); i++) {
          if (stringcompare_hard(as<std::string>(data_names(i)), String(patchcol_str(0)))) {
            patchcol_int = i;
            found_patch_col = true;
          }
        }
        
        StringVector patches_allrows_str_temp (total_rows);
        for (int i = 0; i < total_rows; i++) patches_allrows_str_temp(i) = "0";
        patches_allrows_str = patches_allrows_str_temp; // Should be redone as pop-patches
        unique_patches_str = sort_unique(patches_allrows_str); // Should be redone as pop-patches
      }
      
      if (by_patch_bool && !found_patch_col) {
        LefkoUtils::pop_error("patch identity", "", "", 16);
      }
      
      StringVector poppatches_allrows_str (total_rows);    // By rowstr_temp (total_rows);
      for (int i = 0; i < total_rows; i++) {
        poppatches_allrows_str(i) = pops_allrows_str(i);
        poppatches_allrows_str(i) += " ";
        poppatches_allrows_str(i) += patches_allrows_str(i);
      }
      StringVector unique_poppatches_str = sort_unique(poppatches_allrows_str);
      total_poppatches = static_cast<int>(unique_poppatches_str.length());
      
      arma::uvec pop_by_row = zeros<uvec>(total_rows);
      arma::uvec poppatch_by_row = zeros<uvec>(total_rows);
      for (int i = 0; i < total_rows; i++) {
        for (int j = 0; j < total_pops; j++) {
          if (pops_allrows_str(i) == unique_pops_str(j)) pop_by_row(i) = j;
        }
        for (int j = 0; j < total_poppatches; j++) {
          if (poppatches_allrows_str(i) == unique_poppatches_str(j)) poppatch_by_row(i) = j;
        }
      }
      
      int length_of_max_limit_int = static_cast<int>(max_limit_int.length());
      if (!prop_size_bool && length_of_max_limit_int > 1) {
        if (max_limit_bool) {
          if (by_pop_bool && by_patch_bool) {
            if (length_of_max_limit_int != total_poppatches) {
              String eat_my_shorts = "Argument max_limit should include entries for all ";
              eat_my_shorts += total_poppatches;
              eat_my_shorts += " pop-patches.";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          } else if (by_pop_bool) {
            if (length_of_max_limit_int != total_pops) {
              String eat_my_shorts = "Argument max_limit should include entries for all ";
              eat_my_shorts += total_pops;
              eat_my_shorts += " populations.";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          } else if (by_patch_bool) {
            if (length_of_max_limit_int != total_poppatches) {
              String eat_my_shorts = "Argument max_limit should include entries for all ";
              eat_my_shorts += total_poppatches;
              eat_my_shorts += " patches.";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          }
        }
      }
      
      if (indivcol.isNotNull()) {
        RObject indivcol_RO (indivcol);
        if (is<StringVector>(indivcol_RO)) {
          indivcol_str = as<StringVector>(indivcol_RO);
        } else if (is<IntegerVector>(indivcol_RO) || is<NumericVector>(indivcol_RO)) {
          IntegerVector indivcol_intvec = as<IntegerVector>(indivcol_RO);
          indivcol_int = indivcol_intvec(0);
        } else LefkoUtils::pop_error("indivcol", "a string value", "", 1);
      } else {
        StringVector new_indivcol_str = {"individ"};
        indivcol_str = new_indivcol_str;
      }
      
      if (popcol_int >= 0 && popcol_int <= df_var_num) found_pop_col = true;
      if (patchcol_int >= 0 && patchcol_int <= df_var_num) found_patch_col = true;
      if (indivcol_int >= 0 && indivcol_int <= df_var_num) found_indiv_col = true;
      
      for (int i = 0; i < data_names.length(); i++) {
        if (stringcompare_hard(as<std::string>(data_names(i)), String(indivcol_str(0)))) {
          indivcol_int = i;
          found_indiv_col = true;
        }
      }
      
      if (by_pop_bool && !found_pop_col) {
        LefkoUtils::pop_error("population identity", "", "", 16);
      }
      if (by_patch_bool && !found_patch_col) {
        LefkoUtils::pop_error("patch identity", "", "", 16);
      }
      if (by_indiv_bool && !found_indiv_col) {
        LefkoUtils::pop_error("individual identity", "", "", 16);
      }
      
      if (by_indiv_bool) {
        if (is<IntegerVector>(true_data(indivcol_int))) {
          IntegerVector individuals = as<IntegerVector>(true_data(indivcol_int));
          individuals_allrows_str = StringVector(individuals);
          
        } else if (is<NumericVector>(true_data(indivcol_int))) {
          NumericVector individuals = as<NumericVector>(true_data(indivcol_int));
          individuals_allrows_str = StringVector(individuals);
          
        } else if (is<CharacterVector>(true_data(indivcol_int))) {
          CharacterVector individuals = as<CharacterVector>(true_data(indivcol_int));
          individuals_allrows_str = StringVector(individuals);
          
        }
        
        unique_individuals_str = sort_unique(individuals_allrows_str);
        total_indivs = unique_individuals_str.length();
        
        if (prop_size_bool) {
          // Equal sample proportions, by individual
          
          if (length_of_max_limit_int > 1) {
            String eat_my_shorts = "If prop_size is TRUE, then argument ";
            eat_my_shorts += "max_limit should be empty or set to a single ";
            eat_my_shorts += "value.";
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          
          if (by_patch_bool) {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_poppatches; j++) {
                arma::uvec current_poppatch_rows = find(poppatch_by_row == j);
                int found_indivs_selection = static_cast<int>(current_poppatch_rows.n_elem);
                StringVector current_indiv_selection (found_indivs_selection);
                for (int k = 0; k < found_indivs_selection; k++) {
                  current_indiv_selection(k) = individuals_allrows_str(current_poppatch_rows(k));
                }
                
                StringVector unique_indivs_current_poppatch = sort_unique(current_indiv_selection);
                
                int current_poppatch_sample_limit {0};
                if (max_limit_bool) {
                  double base_num = static_cast<double>(unique_indivs_current_poppatch.length());
                  double base_denom = static_cast<double>(total_indivs);
                  
                  double current_prop = base_num / base_denom;
                  double current_sample_double = current_prop * base_denom;
                  current_poppatch_sample_limit = static_cast<int>(round(current_sample_double));
                } else {
                  current_poppatch_sample_limit = static_cast<int>(unique_indivs_current_poppatch.length());
                }
                
                // Row determination
                CharacterVector sampled_individuals = Rcpp::sample(unique_indivs_current_poppatch,
                  current_poppatch_sample_limit, true);
                
                for (int k = 0; k < current_poppatch_sample_limit; k++) {
                  arma::uvec rows_identified_arma (total_rows, fill::zeros);
                  for (int l = 0; l < total_rows; l++) {
                    if (individuals_allrows_str(l) == sampled_individuals(k) && poppatch_by_row(l) == j) {
                      rows_identified_arma(l) = 1;
                    }
                  }
                  arma::uvec found_rows_for_indiv_arma = find(rows_identified_arma);
                  IntegerVector found_rows_for_indiv = wrap(found_rows_for_indiv_arma);
                  
                  if (j == 0 && k == 0) {
                    new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  } else {
                    DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  }
                }
                
              }
            }
          } else if (by_pop_bool) {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_pops; j++) {
                arma::uvec current_pop_rows = find(pop_by_row == j);
                int found_indivs_selection = static_cast<int>(current_pop_rows.n_elem);
                StringVector current_indiv_selection (found_indivs_selection);
                for (int k = 0; k < found_indivs_selection; k++) {
                  current_indiv_selection(k) = individuals_allrows_str(current_pop_rows(k));
                }
                
                StringVector unique_indivs_current_pop = sort_unique(current_indiv_selection);
                
                int current_pop_sample_limit {0};
                if (max_limit_bool) {
                  double base_num = static_cast<double>(unique_indivs_current_pop.length());
                  double base_denom = static_cast<double>(total_indivs);
                  
                  double current_prop = base_num / base_denom;
                  double current_sample_double = current_prop * base_denom;
                  current_pop_sample_limit = static_cast<int>(round(current_sample_double));
                } else {
                  current_pop_sample_limit = static_cast<int>(unique_indivs_current_pop.length());
                }
                
                // Row determination
                CharacterVector sampled_individuals = Rcpp::sample(unique_indivs_current_pop,
                  current_pop_sample_limit, true);
                
                for (int k = 0; k < current_pop_sample_limit; k++) {
                  arma::uvec rows_identified_arma (total_rows, fill::zeros);
                  for (int l = 0; l < total_rows; l++) {
                    if (individuals_allrows_str(l) == sampled_individuals(k) && pop_by_row(l) == j) {
                      rows_identified_arma(l) = 1;
                    }
                  }
                  arma::uvec found_rows_for_indiv_arma = find(rows_identified_arma);
                  IntegerVector found_rows_for_indiv = wrap(found_rows_for_indiv_arma);
                  
                  if (j == 0 && k == 0) {
                    new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  } else {
                    DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  }
                }
                
              }
            }
          } else {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              StringVector current_indiv_selection = individuals_allrows_str;
              
              StringVector unique_indivs_current = sort_unique(current_indiv_selection);
              int current_sample_limit = static_cast<int>(unique_indivs_current.length());
              
              // Row determination
              CharacterVector sampled_individuals = Rcpp::sample(unique_indivs_current,
                current_sample_limit, true);
              
              for (int k = 0; k < current_sample_limit; k++) {
                arma::uvec rows_identified_arma (total_rows, fill::zeros);
                for (int l = 0; l < total_rows; l++) {
                  if (individuals_allrows_str(l) == sampled_individuals(k)) {
                    rows_identified_arma(l) = 1;
                  }
                }
                arma::uvec found_rows_for_indiv_arma = find(rows_identified_arma);
                IntegerVector found_rows_for_indiv = wrap(found_rows_for_indiv_arma);
                
                if (k == 0) {
                  new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                  hfv_list(i) = clone(new_sampled_data_frame);
                } else {
                  DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                  new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                  hfv_list(i) = clone(new_sampled_data_frame);
                }
              }
              
            }
          }
        } else {
          // User-set sample number, by individual
          if (by_patch_bool) {
            if (length_of_max_limit_int > 1 && length_of_max_limit_int != total_poppatches) {
              String eat_my_shorts = "If prop_size is FALSE, then argument ";
              eat_my_shorts += "max_limit should be empty, set to a single ";
              eat_my_shorts += "value, or set to the number of pop-patches.";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
            
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_poppatches; j++) {
                arma::uvec current_poppatch_rows = find(poppatch_by_row == j);
                int found_indivs_selection = static_cast<int>(current_poppatch_rows.n_elem);
                StringVector current_indiv_selection (found_indivs_selection);
                for (int k = 0; k < found_indivs_selection; k++) {
                  current_indiv_selection(k) = individuals_allrows_str(current_poppatch_rows(k));
                }
                
                StringVector unique_indivs_current_poppatch = sort_unique(current_indiv_selection);
                
                int current_poppatch_sample_limit {0};
                if (max_limit_bool) {
                  if (length_of_max_limit_int > 1) {
                    current_poppatch_sample_limit = max_limit_int(j);
                  } else if (length_of_max_limit_int == 1) {
                    current_poppatch_sample_limit = max_limit_int(0);
                  } else {
                    current_poppatch_sample_limit = default_sample_set;
                  }
                } else {
                  current_poppatch_sample_limit = default_sample_set;
                }
                
                // Row determination
                CharacterVector sampled_individuals = Rcpp::sample(unique_indivs_current_poppatch,
                  current_poppatch_sample_limit, true);
                
                for (int k = 0; k < current_poppatch_sample_limit; k++) {
                  arma::uvec rows_identified_arma (total_rows, fill::zeros);
                  for (int l = 0; l < total_rows; l++) {
                    if (individuals_allrows_str(l) == sampled_individuals(k) && poppatch_by_row(l) == j) {
                      rows_identified_arma(l) = 1;
                    }
                  }
                  arma::uvec found_rows_for_indiv_arma = find(rows_identified_arma);
                  IntegerVector found_rows_for_indiv = wrap(found_rows_for_indiv_arma);
                  
                  if (j == 0 && k == 0) {
                    new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  } else {
                    DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  }
                }
                
              }
            }
          } else if (by_pop_bool) {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_pops; j++) {
                arma::uvec current_pop_rows = find(pop_by_row == j);
                int found_indivs_selection = static_cast<int>(current_pop_rows.n_elem);
                StringVector current_indiv_selection (found_indivs_selection);
                for (int k = 0; k < found_indivs_selection; k++) {
                  current_indiv_selection(k) = individuals_allrows_str(current_pop_rows(k));
                }
                
                StringVector unique_indivs_current_pop = sort_unique(current_indiv_selection);
                
                int current_pop_sample_limit {0};
                if (max_limit_bool) {
                  if (length_of_max_limit_int > 1) {
                    current_pop_sample_limit = max_limit_int(j);
                  } else if (length_of_max_limit_int == 1) {
                    current_pop_sample_limit = max_limit_int(0);
                  } else {
                    current_pop_sample_limit = default_sample_set;
                  }
                } else {
                  current_pop_sample_limit = default_sample_set;
                }
                
                // Row determination
                CharacterVector sampled_individuals = Rcpp::sample(unique_indivs_current_pop,
                  current_pop_sample_limit, true);
                
                for (int k = 0; k < current_pop_sample_limit; k++) {
                  arma::uvec rows_identified_arma (total_rows, fill::zeros);
                  for (int l = 0; l < total_rows; l++) {
                    if (individuals_allrows_str(l) == sampled_individuals(k) && pop_by_row(l) == j) {
                      rows_identified_arma(l) = 1;
                    }
                  }
                  arma::uvec found_rows_for_indiv_arma = find(rows_identified_arma);
                  IntegerVector found_rows_for_indiv = wrap(found_rows_for_indiv_arma);
                  
                  if (j == 0 && k == 0) {
                    new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  } else {
                    DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                    new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                    hfv_list(i) = clone(new_sampled_data_frame);
                  }
                }
                
              }
            }
          } else {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              StringVector current_indiv_selection = individuals_allrows_str;
              
              StringVector unique_indivs_current = sort_unique(current_indiv_selection);
              int current_sample_limit = static_cast<int>(unique_indivs_current.length());
              
              // Row determination
              CharacterVector sampled_individuals = Rcpp::sample(unique_indivs_current,
                current_sample_limit, true);
              
              for (int k = 0; k < current_sample_limit; k++) {
                arma::uvec rows_identified_arma (total_rows, fill::zeros);
                for (int l = 0; l < total_rows; l++) {
                  if (individuals_allrows_str(l) == sampled_individuals(k)) {
                    rows_identified_arma(l) = 1;
                  }
                }
                arma::uvec found_rows_for_indiv_arma = find(rows_identified_arma);
                IntegerVector found_rows_for_indiv = wrap(found_rows_for_indiv_arma);
                
                if (k == 0) {
                  new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                  hfv_list(i) = clone(new_sampled_data_frame);
                } else {
                  DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, found_rows_for_indiv);
                  new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                  hfv_list(i) = clone(new_sampled_data_frame);
                }
              }
              
            }
          }
        }
      } else {
        // By row
        if (prop_size_bool) {
          // Equal sample proportions, by row
          if (length_of_max_limit_int > 1) {
            String eat_my_shorts = "If prop_size is TRUE, then argument ";
            eat_my_shorts += "max_limit should be empty or set to a single ";
            eat_my_shorts += "value.";
            
            throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
          }
          
          if (by_patch_bool) {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_poppatches; j++) {
                arma::uvec current_poppatch_rows = find(poppatch_by_row == j);
                int found_rows_selection = static_cast<int>(current_poppatch_rows.n_elem);
                
                int current_poppatch_sample_limit {0};
                if (max_limit_bool) {
                  double base_num = static_cast<double>(found_rows_selection);
                  double base_denom = static_cast<double>(total_rows);
                  
                  double current_prop = base_num / base_denom;
                  double current_sample_double = current_prop * base_denom;
                  current_poppatch_sample_limit = static_cast<int>(round(current_sample_double));
                } else {
                  current_poppatch_sample_limit = found_rows_selection;
                }
                
                // Row determination
                IntegerVector current_poppatch_rows_int = wrap(current_poppatch_rows);
                IntegerVector sampled_rows = Rcpp::sample(current_poppatch_rows_int,
                  current_poppatch_sample_limit, true);
                
                if (j == 0) {
                  new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  hfv_list(i) = clone(new_sampled_data_frame);
                } else {
                  DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                  hfv_list(i) = clone(new_sampled_data_frame);
                }
              }
            }
          } else if (by_pop_bool) {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_pops; j++) {
                arma::uvec current_pop_rows = find(pop_by_row == j);
                int found_rows_selection = static_cast<int>(current_pop_rows.n_elem);
                
                int current_pop_sample_limit {0};
                if (max_limit_bool) {
                  double base_num = static_cast<double>(found_rows_selection);
                  double base_denom = static_cast<double>(total_rows);
                  
                  double current_prop = base_num / base_denom;
                  double current_sample_double = current_prop * base_denom;
                  current_pop_sample_limit = static_cast<int>(round(current_sample_double));
                } else {
                  current_pop_sample_limit = found_rows_selection;
                }
                
                // Row determination
                IntegerVector current_pop_rows_int = wrap(current_pop_rows);
                IntegerVector sampled_rows = Rcpp::sample(current_pop_rows_int,
                  current_pop_sample_limit, true);
                
                if (j == 0) {
                  new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  hfv_list(i) = clone(new_sampled_data_frame);
                } else {
                  DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                  hfv_list(i) = clone(new_sampled_data_frame);
                }
              }
            }
          } else {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              IntegerVector unique_rows_current = seq(0, (total_rows - 1));
              int current_sample_limit = total_rows;
              
              // Row determination
              IntegerVector sampled_rows = Rcpp::sample(unique_rows_current,
                current_sample_limit, true);
              
              new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
              hfv_list(i) = clone(new_sampled_data_frame);
            }
          }
        } else {
          // User-set sample number, by individual
          if (by_patch_bool) {
            if (length_of_max_limit_int > 1 && length_of_max_limit_int != total_poppatches) {
              String eat_my_shorts = "If prop_size is FALSE, then argument ";
              eat_my_shorts += "max_limit should be empty, set to a single ";
              eat_my_shorts += "value, or set to the number of pop-patches.";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
            
            
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_poppatches; j++) {
                arma::uvec current_poppatch_rows = find(poppatch_by_row == j);
                int found_rows_selection = static_cast<int>(current_poppatch_rows.n_elem);
                
                int current_poppatch_sample_limit {0};
                if (max_limit_bool) {
                  if (length_of_max_limit_int > 1) {
                    current_poppatch_sample_limit = max_limit_int(j);
                  } else if (length_of_max_limit_int == 1) {
                    current_poppatch_sample_limit = max_limit_int(0);
                  } else {
                    current_poppatch_sample_limit = default_sample_set;
                  }
                } else {
                  current_poppatch_sample_limit = default_sample_set;
                }
                
                // Row determination
                IntegerVector current_poppatch_rows_int = wrap(current_poppatch_rows);
                IntegerVector sampled_rows = Rcpp::sample(current_poppatch_rows_int,
                  current_poppatch_sample_limit, true);
                
                if (j == 0) {
                  new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  hfv_list(i) = clone(new_sampled_data_frame);
                } else {
                  DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                  hfv_list(i) = clone(new_sampled_data_frame);
                }
                
              }
            }
          } else if (by_pop_bool) {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              for (int j = 0; j < total_pops; j++) {
                arma::uvec current_pop_rows = find(pop_by_row == j);
                int found_rows_selection = static_cast<int>(current_pop_rows.n_elem);
                
                int current_pop_sample_limit {0};
                if (max_limit_bool) {
                  if (length_of_max_limit_int > 1) {
                    current_pop_sample_limit = max_limit_int(j);
                  } else if (length_of_max_limit_int == 1) {
                    current_pop_sample_limit = max_limit_int(0);
                  } else {
                    current_pop_sample_limit = default_sample_set;
                  }
                } else {
                  current_pop_sample_limit = default_sample_set;
                }
                
                // Row determination
                IntegerVector current_pop_rows_int = wrap(current_pop_rows);
                IntegerVector sampled_rows = Rcpp::sample(current_pop_rows_int,
                  current_pop_sample_limit, true);
                
                if (j == 0) {
                  new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  hfv_list(i) = clone(new_sampled_data_frame);
                } else {
                  DataFrame new_sampled_data_frame_next = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
                  new_sampled_data_frame = df_rbind(new_sampled_data_frame, new_sampled_data_frame_next);
                  hfv_list(i) = clone(new_sampled_data_frame);
                }
              }
            }
          } else {
            for (int i = 0; i < reps_true; i++) {
              DataFrame new_sampled_data_frame;
              
              int current_sample_limit {0};
              if (max_limit_bool) {
                if (length_of_max_limit_int > 0) {
                  current_sample_limit = max_limit_int(0);
                } else {
                  current_sample_limit = default_sample_set;
                }
              } else {
                current_sample_limit = default_sample_set;
              }
                
              // Row determination
              IntegerVector unique_rows_current = seq(0, (total_rows - 1));
              IntegerVector sampled_rows = Rcpp::sample(unique_rows_current,
                current_sample_limit, true);
              
              new_sampled_data_frame = LefkoUtils::df_subset_byrow(true_data, sampled_rows);
              hfv_list(i) = clone(new_sampled_data_frame);
            }
          }
        }
      }
    } else LefkoUtils::pop_error("data", "an hfvdata object", "", 1);
  } else LefkoUtils::pop_error("data", "an hfvdata object", "", 1);
  
  for (int i = 0; i < reps_true; i++) {
    DataFrame new_sampled_data_frame = as<DataFrame>(hfv_list(i));
    new_sampled_data_frame.attr("class") = df_class_lel;
  }
  CharacterVector hfv_list_class = {"hfvlist"};
  hfv_list.attr("class") = hfv_list_class;
  
  return hfv_list;
}

