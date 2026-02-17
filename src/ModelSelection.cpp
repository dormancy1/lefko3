#include <RcppArmadillo.h>
#include <LefkoUtils.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;


// Index of Functions
// 
// 1. List praxis  Workhorse Formula Creator for Function modelsearch()
// 2. List .stovokor  Main Formula Creation for Function modelsearch()
// 3. DataFrame create_pm  Creates a Skeleton Paramnames Object for Use in Function-based Modeling
// 4. NumericVector vrmf_inator  Convert modelextract Coefficient Vector to vrm_frame Vector
// 5. List miniMod  Minimize lefkoMod Object by Conversion to vrm_input Object


//' Workhorse Formula Creator for Function modelsearch()
//'
//' Function \code{praxis()} is the workhorse function used by function
//' \code{stovokor} to create individual vital rate model formulae, which are
//' then used as input in the global model calls used in function
//' \code{\link{modelsearch}()}.
//' 
//' @name praxis
//' 
//' @param surv A vector of strings indicating the names of the variables coding
//' survival.
//' @param obs A vector of strings indicating the names of the variables coding
//' observation status.
//' @param size A vector of strings indicating the names of the variables coding
//' primary size.
//' @param sizeb A vector of strings indicating the names of the variables
//' coding secondary size.
//' @param sizec A vector of strings indicating the names of the variables
//' coding tertiary size.
//' @param repst A vector of strings indicating the names of the variables
//' coding reproductive status.
//' @param fec A vector of strings indicating the names of the variables coding
//' fecundity.
//' @param matstat A vector of strings indicating the names of the variables
//' coding for maturity status.
//' @param historical A logical value indicating whether to create global models
//' with historical effects.
//' @param response An integer coding for the vital rate to be modeled. Codes
//' are as follows: 1) survival, 2) observation status, 3) primary size,
//' 4) secondary size, 5) tertiary size, 6) reproductive status, 7) fecundity,
//' and 8) maturity status.
//' @param suite A string indicating the scope of independent factors included
//' in the global models. Options include \code{"full"}, \code{"main"},
//' \code{"size"}, \code{"rep"}, and \code{"const"}.
//' @param approach A string indicating whether to use mixed model encoding 
//' (\code{"mixed"}) or GLM encoding (\code{"glm"}).
//' @param nojuvs A logical value indicating that juvenile rates should be
//' estimated (\code{FALSE}) or not (\code{TRUE}).
//' @param juvsize A logical value indicating whether to include size terms in
//' juvenile models.
//' @param indiv A string indicating the name of the variable coding individual
//' identity.
//' @param patch A string indicating the name of the variable coding patch
//' identity.
//' @param year A string indicating the name of the variable coding time
//' \emph{t}.
//' @param age A string indicating the name of the variable coding age.
//' @param densitycol The name of the density variable, or \code{"none"}.
//' @param indcova_raw A vector of strings indicating the names in times
//' \emph{t}+1, \emph{t}, and \emph{t}-1 of a specific individual covariate used
//' in the dataset.
//' @param indcovb_raw A vector of strings indicating the names in times
//' \emph{t}+1, \emph{t}, and \emph{t}-1 of a specific individual covariate used
//' in the dataset.
//' @param indcovc_raw A vector of strings indicating the names in times
//' \emph{t}+1, \emph{t}, and \emph{t}-1 of a specific individual covariate used
//' in the dataset.
//' @param sizebused A logical value denoting if secondary size variables are to
//' be used.
//' @param sizecused A logical value denoting if tertiary size variables are to
//' be used.
//' @param grouptest A logical value indicating whether to test for group
//' effect.
//' @param ageused A logical value indicating whether to test for age effect.
//' @param densityused A logical value indicating whether the density variable
//' is to be used.
//' @param indcovaused Logical value indicating whether individual covariate a
//' is used.
//' @param indcovbused Logical value indicating whether individual covariate b
//' is used.
//' @param indcovcused Logical value indicating whether individual covariate c
//' is used.
//' @param annucovaused Logical value indicating whether annual covariate a
//' is used.
//' @param annucovbused Logical value indicating whether annual covariate b
//' is used.
//' @param annucovcused Logical value indicating whether annual covariate c
//' is used.
//' @param pasrand A logical value indicating whether to treat patch as a random
//' variable within mixed models.
//' @param yasrand A logical value indicating whether to treat year as a random
//' variable within mixed models.
//' @param iaasrand A logical value indicating whether to treat indcova as
//' random.
//' @param ibasrand A logical value indicating whether to treat indcovb as
//' random.
//' @param icasrand A logical value indicating whether to treat indcovc as
//' random.
//' @param iaasfac A logical value indicating whether to treat indcova as a
//' factor variable.
//' @param ibasfac A logical value indicating whether to treat indcovb as a
//' factor variable.
//' @param icasfac A logical value indicating whether to treat indcovc as a
//' factor variable.
//' @param fectime An integer indicating whether to use reproductive output in
//' time \emph{t} (2) or time \emph{t}+1 (3) as the response for fecundity.
//' @param repstcheck A boolean variable denoting whether reproductive status is
//' being tested as a response within the suite of vital rates being estimated.
//' @param interactions A boolean value indicating whether all two-way
//' interactions between fixed factors should be created.
//' 
//' @return A list with four elements. The first and second are both one-element
//' string vectors, with the first coding for the adult vital rate global model,
//' and the second coding the juvenile vital rate global model. The third and
//' fourth code the number of terms tested in each of these models,
//' respectively.
//'
//' @keywords internal
//' @noRd
Rcpp::List praxis(const StringVector& surv, const StringVector& obs,
  const StringVector& size, const StringVector& sizeb,
  const StringVector& sizec, const StringVector& repst, const StringVector& fec,
  const StringVector& matstat, bool historical, int response,
  const String& suite, const String& approach, const bool nojuvs,
  bool juvsize, const String& indiv, const String& patch,
  const String& year, const String& age, const String& densitycol,
  const StringVector indcova_raw, const StringVector indcovb_raw,
  const StringVector indcovc_raw, const bool sizebused, const bool sizecused,
  const bool grouptest, const bool ageused, const bool densityused,
  const bool indcovaused, const bool indcovbused, const bool indcovcused,
  const bool annucovaused, const bool annucovbused, const bool annucovcused,
  bool pasrand, bool yasrand, bool iaasrand, bool ibasrand, bool icasrand,
  const bool iaasfac, const bool ibasfac, const bool icasfac, int fectime,
  bool repstcheck, const bool interactions) {
  
  if (nojuvs) juvsize = false;
  
  String randomtackonp {""};
  String randomtackony {""};
  String randomtackoni {""};
  String randomtackonia {""};
  String randomtackonib {""};
  String randomtackonic {""};
  String randomtackon {""};
  
  String fixedtackong {""};
  String fixedtackonp {""};
  String fixedtackony {""};
  String fixedtackonia {""};
  String fixedtackonib {""};
  String fixedtackonic {""};
  String fixedtackon {""};
  
  String fixedtackonaa {""};
  String fixedtackonab {""};
  String fixedtackonac {""};
  
  String sizesuffix;
  String jsizesuffix;
  String fecsuffix;
  
  String fullmainmodel;
  String juvmainmodel;
  String adult_out_model;
  String juv_out_model;
  
  int covcount {0};
  if (indcova_raw(1) != "none" && indcovaused) covcount += 1;
  if (indcovb_raw(1) != "none" && indcovbused) covcount += 1;
  if (indcovc_raw(1) != "none" && indcovcused) covcount += 1;
  
  int total_terms {0};
  int fixedcovcounter {0};
  int randomcovcounter {0};
  
  StringVector indcova (3);
  StringVector indcovb (3);
  StringVector indcovc (3);
  
  if (indcovaused) {
    if (!iaasrand && iaasfac) {
      String major_pain1 = "as.factor(";
      major_pain1 += indcova_raw(1);
      major_pain1 += ")";
      
      String major_pain2 = "as.factor(";
      major_pain2 += indcova_raw(2);
      major_pain2 += ")";
      
      indcova(1) = major_pain1;
      indcova(2) = major_pain2;
    } else indcova = indcova_raw;
  }
  if (indcovbused) {
    if (!ibasrand && ibasfac) {
      String major_pain1 = "as.factor(";
      major_pain1 += indcovb_raw(1);
      major_pain1 += ")";
      
      String major_pain2 = "as.factor(";
      major_pain2 += indcovb_raw(2);
      major_pain2 += ")";
      
      indcovb(1) = major_pain1;
      indcovb(2) = major_pain2;
    } else indcovb = indcovb_raw;
  }
  if (indcovcused) {
    if (!icasrand && icasfac) {
      String major_pain1 = "as.factor(";
      major_pain1 += indcovc_raw(1);
      major_pain1 += ")";
      
      String major_pain2 = "as.factor(";
      major_pain2 += indcovc_raw(2);
      major_pain2 += ")";
      
      indcovc(1) = major_pain1;
      indcovc(2) = major_pain2;
    } else indcovc = indcovc_raw;
  }
  
  // Year, patch, indiv, and group - all possible factor variables
  // Year, patch, and indiv may be random
  if (approach != "mixed") {
    yasrand = false;
    pasrand = false;
  }
  
  if (grouptest) {
    fixedtackong += "as.factor(group2)";
    fixedcovcounter += 1;
    total_terms++;
    
    if (historical) {
      fixedcovcounter += 1;
      fixedtackong += " + as.factor(group1)";
      total_terms++;
    }
    fixedtackon += fixedtackong;
  }
  
  if (year!= "none") {
    if (yasrand) {
      randomtackony += "(1 | ";
      randomtackony += year;
      randomtackony += ")";
      randomcovcounter += 1;
      total_terms++;
      
      randomtackon += randomtackony;
      
    } else {
      fixedtackony += "as.factor(";
      fixedtackony += year;
      fixedtackony += ")";
      total_terms++;
      fixedcovcounter += 1;
      
      if (fixedcovcounter > 1) fixedtackon += " + ";
      fixedtackon += fixedtackony;
    }
  }
  
  if (patch!= "none") {
    if (pasrand) {
      randomtackonp += "(1 | ";
      randomtackonp += patch;
      randomtackonp += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (randomcovcounter > 1) randomtackon += " + ";
      randomtackon += randomtackonp;
      
    } else {
      fixedtackonp += "as.factor(";
      fixedtackonp += patch;
      fixedtackonp += ")";
      total_terms++;
      fixedcovcounter += 1;
      
      if (fixedcovcounter > 1) fixedtackon += " + ";
      fixedtackon += fixedtackonp;
    }
  }
  
  if (indiv != "none" && approach == "mixed") {
    randomtackoni += "(1 | ";
    randomtackoni += indiv;
    randomtackoni += ")";
    total_terms++;
    randomcovcounter += 1;
    
    if (randomcovcounter > 1) randomtackon += " + ";
    randomtackon += randomtackoni;
  }
  
  // Add individual covariates to tacked-on sections
  if (indcova(1) != "none" && indcovaused) {
    if (!iaasrand) {
      fixedtackonia += indcova(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonia += " + ";
        fixedtackonia += indcova(2);
        fixedcovcounter += 1;
        total_terms++;
      }
      
      if ((fixedcovcounter > 2 && historical) || (fixedcovcounter > 1 && !historical))
        fixedtackon += " + ";
      fixedtackon += fixedtackonia;
      
    } else {
      randomtackonia += "(1 | ";
      randomtackonia += indcova(1);
      randomtackonia += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (historical) {
        randomtackonia += " + (1 | ";
        randomtackonia += indcova(2);
        randomtackonia += ")";
        randomcovcounter += 1;
        total_terms++;
      }
      
      if ((randomcovcounter > 2 && historical) || (randomcovcounter > 1 && !historical))
        randomtackon += " + ";
      randomtackon += randomtackonia;
    }
  }
  
  if (indcovb(1) != "none" && indcovbused) {
    if (!ibasrand) {
      fixedtackonib += indcovb(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonib += " + ";
        fixedtackonib += indcovb(2);
        fixedcovcounter += 1;
        total_terms++;
      }
      
      if ((fixedcovcounter > 2 && historical) || (fixedcovcounter > 1 && !historical))
        fixedtackon += " + ";
      fixedtackon += fixedtackonib;
      
    } else {
      randomtackonib += "(1 | ";
      randomtackonib += indcovb(1);
      randomtackonib += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (historical) {
        randomtackonib += " + (1 | ";
        randomtackonib += indcovb(2);
        randomtackonib += ")";
        randomcovcounter += 1;
        total_terms++;
      }
      
      if ((randomcovcounter > 2 && historical) || (randomcovcounter > 1 && !historical))
        randomtackon += " + ";
      randomtackon += randomtackonib;
    }
  }
  
  if (indcovc(1) != "none" && indcovcused) {
    if (!icasrand) {
      fixedtackonic += indcovc(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcovc(2);
        total_terms++;
      }
      
      if ((fixedcovcounter > 2 && historical) || (fixedcovcounter > 1 && !historical))
        fixedtackon += " + ";
      fixedtackon += fixedtackonic;
      
    } else {
      randomtackonic += "(1 | ";
      randomtackonic += indcovc(1);
      randomtackonic += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (historical) {
        randomtackonic += " + (1 | ";
        randomtackonic += indcovc(2);
        randomtackonic += ")";
        randomcovcounter += 1;
        total_terms++;
      }
      
      if ((randomcovcounter > 2 && historical) || (randomcovcounter > 1 && !historical))
        randomtackon += " + ";
      randomtackon += randomtackonic;
    }
  }
  
  if (interactions && !iaasrand && !ibasrand) {
    if ((indcova(1) != "none" && indcovb(1) != "none") && (indcovaused && indcovbused)) {
      fixedtackonib += " + ";
      fixedtackonib += indcova(1);
      fixedtackonib += ":";
      fixedtackonib += indcovb(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonib += " + ";
        fixedtackonib += indcova(2);
        fixedtackonib += ":";
        fixedtackonib += indcovb(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonib += " + ";
        fixedtackonib += indcova(1);
        fixedtackonib += ":";
        fixedtackonib += indcovb(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonib += " + ";
        fixedtackonib += indcova(2);
        fixedtackonib += ":";
        fixedtackonib += indcovb(1);
        fixedcovcounter += 1;
        total_terms++;
      }
    }
  }
  if (interactions && !iaasrand && !icasrand) {
    if ((indcova(1) != "none" && indcovc(1) != "none") && (indcovaused && indcovcused)) {
      fixedtackonic += " + ";
      fixedtackonic += indcova(1);
      fixedtackonic += ":";
      fixedtackonic += indcovc(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcova(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcova(1);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcova(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(1);
        fixedcovcounter += 1;
        total_terms++;
      }
    }
  }
  if (interactions && !ibasrand && !icasrand) {
    if ((indcovb(1) != "none" && indcovc(1) != "none") && (indcovbused && indcovcused)) {
      fixedtackonic += " + ";
      fixedtackonic += indcovb(1);
      fixedtackonic += ":";
      fixedtackonic += indcovc(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcovb(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcovb(1);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcovb(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(1);
        fixedcovcounter += 1;
        total_terms++;
      }
    }
  }
  
  if (annucovaused || annucovbused || annucovcused) {
    if (annucovaused) {
      fixedtackonaa += " + annucova2";
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonaa += " + annucova1";
        fixedcovcounter += 1;
        total_terms++;
      }
    }
      
    if (annucovbused) {
      fixedtackonab += " + annucovb2";
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonab += " + annucovb1";
        fixedcovcounter += 1;
        total_terms++;
      }
    }
    
    if (annucovcused) {
      fixedtackonac += " + annucovc2";
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonac += " + annucovc1";
        fixedcovcounter += 1;
        total_terms++;
      }
    }
    
    if (interactions) {
      if (annucovaused && annucovbused) {
        fixedtackonaa += " + annucova2:annucovb2";
        fixedcovcounter += 1;
        total_terms++;
        
        if (historical) {
          fixedtackonaa += " + annucova2:annucovb1";
          fixedtackonaa += " + annucova1:annucovb2";
          fixedtackonaa += " + annucova1:annucovb1";
          
          fixedcovcounter += 3;
          total_terms += 3;
        }
      }
      
      if (annucovaused && annucovcused) {
        fixedtackonaa += " + annucova2:annucovc2";
        fixedcovcounter += 1;
        total_terms++;
        
        if (historical) {
          fixedtackonaa += " + annucova2:annucovc1";
          fixedtackonaa += " + annucova1:annucovc2";
          fixedtackonaa += " + annucova1:annucovc1";
          
          fixedcovcounter += 3;
          total_terms += 3;
        }
      }
      
      if (annucovbused && annucovcused) {
        fixedtackonab += " + annucovb2:annucovc2";
        fixedcovcounter += 1;
        total_terms++;
        
        if (historical) {
          fixedtackonab += " + annucovb2:annucovc1";
          fixedtackonab += " + annucovb1:annucovc2";
          fixedtackonab += " + annucovb1:annucovc1";
          
          fixedcovcounter += 3;
          total_terms += 3;
        }
      }
    }
  }
  
  fixedtackon += fixedtackonaa;
  fixedtackon += fixedtackonab;
  fixedtackon += fixedtackonac;
  
  // Main model patterns
  int modelcounter {0};
  int juvmodelcounter {0};
  
  // First the juvenile model pattern
  if (!nojuvs) {
    juvmainmodel = " ~ ";
    
    if (suite != "const") {
      if (juvsize && suite != "rep") {
        juvmainmodel += size(1);
        juvmodelcounter = 1;
        
        if (sizebused) {
          if (juvmodelcounter > 0) juvmainmodel += " + ";
          juvmainmodel += sizeb(1);
          juvmodelcounter += 1;
          
          if (suite == "full" || interactions) {
            juvmainmodel += " + ";
            juvmainmodel += size(1);
            juvmainmodel += ":";
            juvmainmodel += sizeb(1);
            juvmodelcounter += 1;
          }
        }
        if (sizecused) {
          if (juvmodelcounter > 0) juvmainmodel += " + ";
          juvmainmodel += sizec(1);
          juvmodelcounter += 1;
          
          if (suite == "full" || interactions) {
            juvmainmodel += " + ";
            juvmainmodel += size(1);
            juvmainmodel += ":";
            juvmainmodel += sizec(1);
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + ";
              juvmainmodel += sizeb(1);
              juvmainmodel += ":";
              juvmainmodel += sizec(1);
              juvmodelcounter += 1;
            }
          }
        }
        if (densityused) {
          if (juvmodelcounter > 0) juvmainmodel += " + ";
          juvmainmodel += densitycol;
          juvmodelcounter += 1;
          
          if (interactions) {
            juvmainmodel += " + ";
            juvmainmodel += size(1);
            juvmainmodel += ":";
            juvmainmodel += densitycol;
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + ";
              juvmainmodel += sizeb(1);
              juvmainmodel += ":";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
            if (sizecused) {
              juvmainmodel += " + ";
              juvmainmodel += sizec(1);
              juvmainmodel += ":";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
          }
        }
        
        if (interactions) {
          if (annucovaused) {
            juvmainmodel += " + annucova2:";
            juvmainmodel += size(1);
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + annucova2:";
              juvmainmodel += sizeb(1);
              juvmodelcounter += 1;
            }
            
            if (sizecused) {
              juvmainmodel += " + annucova2:";
              juvmainmodel += sizec(1);
              juvmodelcounter += 1;
            }
            
            if (densityused) {
              juvmainmodel += " + annucova2:";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
            
            if (historical) {
              juvmainmodel += " + annucova2:";
              juvmainmodel += size(2);
              
              juvmainmodel += " + annucova1:";
              juvmainmodel += size(1);
              
              juvmainmodel += " + annucova1:";
              juvmainmodel += size(2);
              juvmodelcounter += 3;
              
              if (densityused) {
                juvmainmodel += " + annucova1:";
                juvmainmodel += densitycol;
                juvmodelcounter += 1;
              }
              
              if (sizebused) {
                juvmainmodel += " + annucova2:";
                juvmainmodel += sizeb(2);
                
                juvmainmodel += " + annucova1:";
                juvmainmodel += sizeb(1);
                
                juvmainmodel += " + annucova1:";
                juvmainmodel += sizeb(2);
                juvmodelcounter += 3;
              }
              
              if (sizecused) {
                juvmainmodel += " + annucova2:";
                juvmainmodel += sizec(1);
                
                juvmainmodel += " + annucova1:";
                juvmainmodel += sizec(1);
                
                juvmainmodel += " + annucova1:";
                juvmainmodel += sizec(2);
                juvmodelcounter += 3;
              }
            }
          }
          
          if (annucovbused) {
            juvmainmodel += " + annucovb2:";
            juvmainmodel += size(1);
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + annucovb2:";
              juvmainmodel += sizeb(1);
              juvmodelcounter += 1;
            }
            
            if (sizecused) {
              juvmainmodel += " + annucovb2:";
              juvmainmodel += sizec(1);
              juvmodelcounter += 1;
            }
            
            if (densityused) {
              juvmainmodel += " + annucovb2:";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
            
            if (historical) {
              juvmainmodel += " + annucovb2:";
              juvmainmodel += size(2);
              
              juvmainmodel += " + annucovb1:";
              juvmainmodel += size(1);
              
              juvmainmodel += " + annucovb1:";
              juvmainmodel += size(2);
              juvmodelcounter += 3;
              
              if (densityused) {
                juvmainmodel += " + annucovb1:";
                juvmainmodel += densitycol;
                juvmodelcounter += 1;
              }
              
              if (sizebused) {
                juvmainmodel += " + annucovb2:";
                juvmainmodel += sizeb(2);
                
                juvmainmodel += " + annucovb1:";
                juvmainmodel += sizeb(1);
                
                juvmainmodel += " + annucovb1:";
                juvmainmodel += sizeb(2);
                juvmodelcounter += 3;
              }
              
              if (sizecused) {
                juvmainmodel += " + annucovb2:";
                juvmainmodel += sizec(1);
                
                juvmainmodel += " + annucovb1:";
                juvmainmodel += sizec(1);
                
                juvmainmodel += " + annucovb1:";
                juvmainmodel += sizec(2);
                juvmodelcounter += 3;
              }
            }
          }
          
          if (annucovcused) {
            juvmainmodel += " + annucovc2:";
            juvmainmodel += size(1);
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + annucovc2:";
              juvmainmodel += sizeb(1);
              juvmodelcounter += 1;
            }
            
            if (sizecused) {
              juvmainmodel += " + annucovc2:";
              juvmainmodel += sizec(1);
              juvmodelcounter += 1;
            }
            
            if (densityused) {
              juvmainmodel += " + annucovc2:";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
            
            if (historical) {
              juvmainmodel += " + annucovc2:";
              juvmainmodel += size(2);
              
              juvmainmodel += " + annucovc1:";
              juvmainmodel += size(1);
              
              juvmainmodel += " + annucovc1:";
              juvmainmodel += size(2);
              juvmodelcounter += 3;
              
              if (densityused) {
                juvmainmodel += " + annucovc1:";
                juvmainmodel += densitycol;
                juvmodelcounter += 1;
              }
              
              if (sizebused) {
                juvmainmodel += " + annucovc2:";
                juvmainmodel += sizeb(2);
                
                juvmainmodel += " + annucovc1:";
                juvmainmodel += sizeb(1);
                
                juvmainmodel += " + annucovc1:";
                juvmainmodel += sizeb(2);
                juvmodelcounter += 3;
              }
              
              if (sizecused) {
                juvmainmodel += " + annucovc2:";
                juvmainmodel += sizec(1);
                
                juvmainmodel += " + annucovc1:";
                juvmainmodel += sizec(1);
                
                juvmainmodel += " + annucovc1:";
                juvmainmodel += sizec(2);
                juvmodelcounter += 3;
              }
            }
          }
        }
      } else if (densityused) {
        if (juvmodelcounter > 0) juvmainmodel += " + ";
        juvmainmodel += densitycol;
        juvmodelcounter += 1;
      } else if (fixedcovcounter == 0) {
        juvmainmodel += "1";
        juvmodelcounter += 1;
      }
    } else  if (fixedcovcounter == 0) {
      juvmainmodel += "1";
      juvmodelcounter += 1;
    }
    
    if (juvmodelcounter == 0) {
      juvmainmodel += "1";
      juvmodelcounter += 1;
    }
    if (fixedcovcounter > 0) {
      if (juvmodelcounter > 0) juvmainmodel += " + ";
      juvmainmodel += fixedtackon;
    }
    if (randomcovcounter > 0) {
      if (juvmodelcounter > 0) juvmainmodel += " + ";
      juvmainmodel += randomtackon;
    }
    
    juvmodelcounter += total_terms;
  }
  
  // Adult model pattern
  fullmainmodel = " ~ ";
  
  if (age != "none") {
    fullmainmodel += age;
    modelcounter = 1;
  }
  
  if (densityused) {
    if (modelcounter > 0) fullmainmodel += " + ";
    fullmainmodel += densitycol;
    modelcounter += 1;
  }
  
  if (suite == "full" || suite == "main" || suite == "size" || suite == "rep") {
    if (suite != "rep") {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += size(1);
      modelcounter += 1;
      
      if (interactions) {
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":annucova2";
          modelcounter += 1;
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
        }
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += size(2);
        modelcounter += 1;

        if (interactions) {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
          
          if (annucovaused) {
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":annucova2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":annucova2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":annucova1";
            modelcounter += 1;
          }
          if (annucovbused) {
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":annucovb2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":annucovb2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
          }
          if (annucovcused) {
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":annucovc2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":annucovc2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(2);
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
          }
        }
      }
      
      if (sizebused) {
        if (modelcounter > 0) fullmainmodel += " + ";
        fullmainmodel += sizeb(1);
        modelcounter += 1;
        
        if (interactions) {
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
          
          if (annucovaused) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":annucova2";
            modelcounter += 1;
          }
          if (annucovbused) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":annucovb2";
            modelcounter += 1;
          }
          if (annucovcused) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":annucovc2";
            modelcounter += 1;
          }
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(2);
          modelcounter += 1;

          if (interactions) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += size(2);
            modelcounter += 1;
            
            if (densityused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":";
              fullmainmodel += densitycol;
              modelcounter += 1;
            }
            
            if (annucovaused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(1);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":annucova1";
              modelcounter += 1;
            }
            if (annucovbused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(1);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":annucovb1";
              modelcounter += 1;
            }
            if (annucovcused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(1);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizeb(2);
              fullmainmodel += ":annucovc1";
              modelcounter += 1;
            }
          }
        }
      }
      
      if (sizecused) {
        if (modelcounter > 0) fullmainmodel += " + ";
        fullmainmodel += sizec(1);
        modelcounter += 1;
        
        if (interactions) {
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
          
          if (annucovaused) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":annucova2";
            modelcounter += 1;
          }
          if (annucovbused) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":annucovb2";
            modelcounter += 1;
          }
          if (annucovcused) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":annucovc2";
            modelcounter += 1;
          }
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(2);
          modelcounter += 1;
          
          if (interactions) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += size(2);
            modelcounter += 1;
            
            if (sizebused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(1);
              fullmainmodel += ":";
              fullmainmodel += sizec(2);
              modelcounter += 1;
            
              fullmainmodel += " + ";
              fullmainmodel += sizec(1);
              fullmainmodel += ":";
              fullmainmodel += sizeb(2);
              modelcounter += 1;
            }
            
            if (densityused) {
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":";
              fullmainmodel += densitycol;
              modelcounter += 1;
            }
            
            if (annucovaused) {
              fullmainmodel += " + ";
              fullmainmodel += sizec(1);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":annucova1";
              modelcounter += 1;
            }
            if (annucovbused) {
              fullmainmodel += " + ";
              fullmainmodel += sizec(1);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":annucovb1";
              modelcounter += 1;
            }
            if (annucovcused) {
              fullmainmodel += " + ";
              fullmainmodel += sizec(1);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(2);
              fullmainmodel += ":annucovc1";
              modelcounter += 1;
            }
          }
        }
      }
    }
    
    if (suite != "size") {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += repst(1);
      modelcounter += 1;
      
      if (interactions) {
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":annucova2";
          modelcounter += 1;
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
        }
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += repst(2);
        modelcounter += 1;
        
        if (interactions) {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
          
          if (annucovaused) {
            fullmainmodel += " + ";
            fullmainmodel += repst(1);
            fullmainmodel += ":annucova2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":annucova2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":annucova1";
            modelcounter += 1;
          }
          if (annucovbused) {
            fullmainmodel += " + ";
            fullmainmodel += repst(1);
            fullmainmodel += ":annucovb2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":annucovb2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
          }
          if (annucovcused) {
            fullmainmodel += " + ";
            fullmainmodel += repst(1);
            fullmainmodel += ":annucovc2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":annucovc2";
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += repst(2);
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
          }
        }
      }
    }
    
    if (suite == "full" || (suite == "main" && interactions)) {
      fullmainmodel += " + ";
      fullmainmodel += size(1);
      fullmainmodel += ":";
      fullmainmodel += repst(1);
      modelcounter += 1;
      
      if (sizebused) {
        fullmainmodel += " + ";
        fullmainmodel += sizeb(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
      }
      
      if (sizecused) {
        fullmainmodel += " + ";
        fullmainmodel += sizec(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += repst(2);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += sizeb(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(2);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += sizec(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(2);
          modelcounter += 1;
        }
      }
    }
    
    if (interactions) {
      if (age != "none" && ageused) {
        if (suite == "full" || suite == "main" || suite == "size") {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += size(1);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        }
        
        if (suite == "full" || suite == "main" || suite == "repst") {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += repst(1);
          modelcounter += 1;
        }
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":annucova2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":annucova1";
            modelcounter += 1;
          }
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
          }
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
          }
        }
        
        if (suite == "full" || suite == "main" || suite == "size") {
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += size(2);
            modelcounter += 1;
            
            if (sizebused) {
              fullmainmodel += " + ";
              fullmainmodel += age;
              fullmainmodel += ":";
              fullmainmodel += sizeb(2);
              modelcounter += 1;
            }
            if (sizecused) {
              fullmainmodel += " + ";
              fullmainmodel += age;
              fullmainmodel += ":";
              fullmainmodel += sizec(2);
              modelcounter += 1;
            }
          }
          if (suite == "full" || suite == "main" || suite == "rep") {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += repst(2);
            modelcounter += 1;
          }
        }
      }
      
      if (densityused) {
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += densitycol;
          fullmainmodel += ":annucova2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += densitycol;
            fullmainmodel += ":annucova1";
            modelcounter += 1;
          }
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += densitycol;
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += densitycol;
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
          }
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += densitycol;
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += densitycol;
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
          }
        }
      }
      
      if (indcova(1) != "none" && !iaasrand && indcovaused) {
        fullmainmodel += " + ";
        fullmainmodel += indcova(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcova(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical && indcova(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcova(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcova(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
        
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":annucova2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":annucova1";
            modelcounter += 1;
            
            if (indcova(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcova(2);
              fullmainmodel += ":annucova1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcova(2);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
            }
          }
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
            
            if (indcova(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcova(2);
              fullmainmodel += ":annucovb1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcova(2);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
            }
          }
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
            
            if (indcova(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcova(2);
              fullmainmodel += ":annucovc1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcova(2);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
            }
          }
        }
      }
      if (indcovb(1) != "none" && !ibasrand && indcovbused) {
        fullmainmodel += " + ";
        fullmainmodel += indcovb(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcovb(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical && indcovb(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcovb(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
        
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":annucova2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":annucova1";
            modelcounter += 1;
            
            if (indcovb(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcovb(2);
              fullmainmodel += ":annucova1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcovb(2);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
            }
          }
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
            
            if (indcovb(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcovb(2);
              fullmainmodel += ":annucovb1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcovb(2);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
            }
          }
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
            
            if (indcovb(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcovb(2);
              fullmainmodel += ":annucovc1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcovb(2);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
            }
          }
        }
      }
      
      if (indcovc(1) != "none" && !icasrand && indcovcused) {
        fullmainmodel += " + ";
        fullmainmodel += indcovc(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcovc(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical && indcovc(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcovc(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
        
        if (annucovaused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":annucova2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":annucova1";
            modelcounter += 1;
            
            if (indcovc(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcovc(2);
              fullmainmodel += ":annucova1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcovc(2);
              fullmainmodel += ":annucova2";
              modelcounter += 1;
            }
          }
        }
        if (annucovbused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":annucovb2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":annucovb1";
            modelcounter += 1;
            
            if (indcovc(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcovc(2);
              fullmainmodel += ":annucovb1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcovc(2);
              fullmainmodel += ":annucovb2";
              modelcounter += 1;
            }
          }
        }
        if (annucovcused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":annucovc2";
          modelcounter += 1;
          
          if (historical) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":annucovc1";
            modelcounter += 1;
            
            if (indcovc(2) != "none") {
              fullmainmodel += " + ";
              fullmainmodel += indcovc(2);
              fullmainmodel += ":annucovc1";
              modelcounter += 1;
              
              fullmainmodel += " + ";
              fullmainmodel += indcovc(2);
              fullmainmodel += ":annucovc2";
              modelcounter += 1;
            }
          }
        }
      }
    }
    
    if (fixedcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += fixedtackon;
    }
    
    if (randomcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += randomtackon;
    }
    
  } else if (suite == "rep") {
    if (age != "none" || covcount > 0) fullmainmodel += " + ";
    fullmainmodel += repst(1);
    modelcounter += 1;
    
    if (historical && repstcheck) {
      fullmainmodel += " + ";
      fullmainmodel += repst(2);
      modelcounter += 1;
      
      if (interactions) {
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        modelcounter += 1;
      }
    }
    
    if (fixedcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += fixedtackon;
    }
    
    if (randomcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += randomtackon;
    }
    
  } else if (suite == "cons") {
    if (fixedcovcounter == 0 && modelcounter == 0) fullmainmodel += "1";
    modelcounter += 1;
    
    if (fixedcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += fixedtackon;
    }
    
    if (randomcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += randomtackon;
    }
    
  } else {
    fullmainmodel = "none";
  }
  
  // Build the global models
  if (response == 1) {
    String surv0_first(surv[0]);
    adult_out_model = surv0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = surv0_first;
      juv_out_model += juvmainmodel;
      
    } else {
      juv_out_model = "none";
    }
    
  } else if (response == 2) {
    String obs0_first(obs[0]);
    adult_out_model = obs0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = obs0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 3) {
    String size0_first(size[0]);
    adult_out_model = size0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = size0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 4) {
    String sizeb0_first(sizeb[0]);
    adult_out_model = sizeb0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = sizeb0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 5) {
    String sizec0_first(sizec[0]);
    adult_out_model = sizec0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = sizec0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 6) {
    String repst0_first(repst[0]);
    adult_out_model = repst0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = repst0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 7) {
    if (fectime == 3) {
      String fec0_first(fec[0]);
      adult_out_model = fec0_first;
    } else {
      String fec1_first(fec[1]);
      adult_out_model = fec1_first;
    }
    adult_out_model += fullmainmodel;
  } else  if (response == 8) {
    if (!nojuvs) {
      String matstat0_first(matstat[0]);
      juv_out_model = matstat0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  }
  modelcounter += total_terms;
  
  List output(4);
  output(0) = adult_out_model;
  output(1) = juv_out_model;
  output(2) = modelcounter;
  output(3) = juvmodelcounter;
  
  CharacterVector out_names = {"adult_model", "juv_model", "total_terms",
    "juv_total_terms"};
  
  if (fullmainmodel == "none") {
    output["adult_model"] = 1;
    
    if (!nojuvs) {
      output["juv_model"] = 1;
    }
  }
  if (nojuvs) {
    output["juv_model"] = 0;
  }
  
  return output;
}

//' Main Formula Creation for Function modelsearch()
//'
//' Function \code{stovokor()} creates the list of formulae to be used as input
//' in the global model calls used in function \code{\link{modelsearch}()}.
//' 
//' @name stovokor
//' 
//' @param surv A vector of strings indicating the names of the variables coding
//' survival.
//' @param obs A vector of strings indicating the names of the variables coding
//' observation status.
//' @param size A vector of strings indicating the names of the variables coding
//' primary size.
//' @param sizeb A vector of strings indicating the names of the variables
//' coding secondary size.
//' @param sizec A vector of strings indicating the names of the variables
//' coding tertiary size.
//' @param repst A vector of strings indicating the names of the variables
//' coding reproductive status.
//' @param fec A vector of strings indicating the names of the variables coding
//' fecundity.
//' @param matstat A vector of strings indicating the names of the variables
//' coding for maturity status.
//' @param vitalrates A vector of strings indicating which vital rates will be
//' estimated.
//' @param historical A logical value indicating whether to create global models
//' with historical effects.
//' @param suite A string indicating the scope of independent factors included
//' in the global models. Options include \code{"full"}, \code{"main"},
//' \code{"size"}, \code{"rep"}, and \code{"const"}.
//' @param approach A string indicating whether to use mixed model encoding 
//' (\code{"mixed"}) or GLM encoding (\code{"glm"}).
//' @param nojuvs A logical value indicating that juvenile rates should be
//' estimated (\code{FALSE}) or not (\code{TRUE}).
//' @param juvsize A logical value indicating whether to include size terms in
//' juvenile models.
//' @param indiv A string indicating the name of the variable coding individual
//' identity.
//' @param patch A string indicating the name of the variable coding patch
//' identity.
//' @param year A string indicating the name of the variable coding time
//' \emph{t}.
//' @param age A string indicating the name of the variable coding age.
//' @param densitycol The name of the density variable, or \code{"none"}.
//' @param indcova A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indcovb A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indcovc A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param sizebused A logical vector indicating if secondary size variables are
//' to be tested in each of the 14 vital rate models.
//' @param sizecused A logical vector indicating if tertiary size variables are
//' to be tested in each of the 14 vital rate models.
//' @param grouptest A logical vector indicating if group is to be tested in
//' each of the 14 vital rate models.
//' @param ageused A logical vector indicating if age is to be tested in each of
//' the 14 vital rate models.
//' @param densityused A logical vector indicating if density is to be tested in
//' each of the 14 vital rate models.
//' @param indcovaused A logical vector indicating if individual covariate a is
//' to be tested in each of the 14 vital rate models.
//' @param indcovbused A logical vector indicating if individual covariate b is
//' to be tested in each of the 14 vital rate models.
//' @param indcovcused A logical vector indicating if individual covariate c is
//' to be tested in each of the 14 vital rate models.
//' @param annucovaused A logical vector indicating if annual covariate a is to
//' be tested in each of the 14 vital rate models.
//' @param annucovbused A logical vector indicating if annual covariate b is to
//' be tested in each of the 14 vital rate models.
//' @param annucovcused A logical vector indicating if annual covariate c is to
//' be tested in each of the 14 vital rate models.
//' @param pasrand A logical value indicating whether to treat patch as a random
//' variable within mixed models.
//' @param yasrand A logical value indicating whether to treat year as a random
//' variable within mixed models.
//' @param iaasrand A logical value indicating whether to treat indcova as
//' random.
//' @param ibasrand A logical value indicating whether to treat indcovb as
//' random.
//' @param icasrand A logical value indicating whether to treat indcovc as
//' random.
//' @param iaasfac A logical value indicating whether to treat indcova as a
//' factor variable.
//' @param ibasfac A logical value indicating whether to treat indcovb as a
//' factor variable.
//' @param icasfac A logical value indicating whether to treat indcovc as a
//' factor variable.
//' @param fectime An integer indicating whether to use reproductive output in
//' time \emph{t} (2) or time \emph{t}+1 (3) as the response for fecundity.
//' @param size_zero A boolean variable indicating whether the primary size
//' model is zero-inflated.
//' @param sizeb_zero A boolean variable indicating whether the secondary size
//' model is zero-inflated.
//' @param sizec_zero A boolean variable indicating whether the tertiary size
//' model is zero-inflated.
//' @param jsize_zero A boolean variable indicating whether the juvenile primary
//' size model is zero-inflated.
//' @param jsizeb_zero A boolean variable indicating whether the juvenile
//' secondary size model is zero-inflated.
//' @param jsizec_zero A boolean variable indicating whether the juvenile
//' tertiary size model is zero-inflated.
//' @param interactions A boolean variable indicating whether to include two-way
//' interactions among all fixed variables.
//' 
//' @return A list of four lists. The first list includes the 14 main global
//' models covering all 14 vital rate responses, followed by an associated
//' \code{paramnames} object, and a vector of 14 integers showing the number of
//' terms tested in each respective model. The next three lists repeat this
//' structure without a new paramnames object, but for a reduced set of models
//' (\code{alternate}), a glm-only version of the models (\code{glm.alternate}),
//' and a version without any individual covariates (\code{nocovs.alternate}).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.stovokor)]]
Rcpp::List stovokor(const StringVector& surv, const StringVector& obs,
  const StringVector& size, const StringVector& sizeb,
  const StringVector& sizec, const StringVector& repst, const StringVector& fec,
  const StringVector& matstat, const StringVector& vitalrates, bool historical,
  StringVector& suite, const String& approach, const bool nojuvs,
  bool juvsize, const String& indiv, const String& patch,
  const String& year, const String& age, const String& densitycol,
  const StringVector& indcova, const StringVector& indcovb,
  const StringVector& indcovc, const bool sizebused, const bool sizecused,
  const LogicalVector& grouptest, const LogicalVector& ageused,
  const LogicalVector& densityused, const LogicalVector& indcovaused,
  const LogicalVector& indcovbused, const LogicalVector& indcovcused,
  const LogicalVector& annucovaused, const LogicalVector& annucovbused,
  const LogicalVector& annucovcused, const bool pasrand, const bool yasrand,
  const bool iaasrand, const bool ibasrand, const bool icasrand,
  const bool iaasfac, const bool ibasfac, const bool icasfac, const int fectime,
  const bool size_zero, const bool sizeb_zero, const bool sizec_zero,
  const bool jsize_zero, const bool jsizeb_zero, const bool jsizec_zero,
  const bool interactions) {
  
  if (nojuvs) juvsize = false;
  
  int nvitalrates = vitalrates.length();
  bool survcheck {false};
  bool obscheck {false};
  bool sizecheck {false};
  bool repstcheck {false};
  bool feccheck  {false};
  
  const int sizel {static_cast<int>(size.length())};
  const int repstl {static_cast<int>(repst.length())};
  
  String fullsurvmodel {"none"};
  String fullobsmodel {"none"};
  String fullsizemodel {"none"};
  String fullsizebmodel {"none"};
  String fullsizecmodel {"none"};
  String fullrepstmodel {"none"};
  String fullfecmodel {"none"};
  
  String juvsurvmodel {"none"};
  String juvobsmodel {"none"};
  String juvsizemodel {"none"};
  String juvsizebmodel {"none"};
  String juvsizecmodel {"none"};
  String juvrepstmodel {"none"};
  String juvmatstmodel {"none"};
  
  String alt_fullsurvmodel {"none"};
  String alt_fullobsmodel {"none"};
  String alt_fullsizemodel {"none"};
  String alt_fullsizebmodel {"none"};
  String alt_fullsizecmodel {"none"};
  String alt_fullrepstmodel {"none"};
  String alt_fullfecmodel {"none"};
  
  String alt_juvsurvmodel {"none"};
  String alt_juvobsmodel {"none"};
  String alt_juvsizemodel {"none"};
  String alt_juvsizebmodel {"none"};
  String alt_juvsizecmodel {"none"};
  String alt_juvrepstmodel {"none"};
  String alt_juvmatstmodel {"none"};
  
  String nocovs_fullsurvmodel {"none"};
  String nocovs_fullobsmodel {"none"};
  String nocovs_fullsizemodel {"none"};
  String nocovs_fullsizebmodel {"none"};
  String nocovs_fullsizecmodel {"none"};
  String nocovs_fullrepstmodel {"none"};
  String nocovs_fullfecmodel {"none"};
  
  String nocovs_juvsurvmodel {"none"};
  String nocovs_juvobsmodel {"none"};
  String nocovs_juvsizemodel {"none"};
  String nocovs_juvsizebmodel {"none"};
  String nocovs_juvsizecmodel {"none"};
  String nocovs_juvrepstmodel {"none"};
  String nocovs_juvmatstmodel {"none"};
  
  String glm_fullsurvmodel {"none"};
  String glm_fullobsmodel {"none"};
  String glm_fullsizemodel {"none"};
  String glm_fullsizebmodel {"none"};
  String glm_fullsizecmodel {"none"};
  String glm_fullrepstmodel {"none"};
  String glm_fullfecmodel {"none"};
  
  String glm_juvsurvmodel {"none"};
  String glm_juvobsmodel {"none"};
  String glm_juvsizemodel {"none"};
  String glm_juvsizebmodel {"none"};
  String glm_juvsizecmodel {"none"};
  String glm_juvrepstmodel {"none"};
  String glm_juvmatstmodel {"none"};
  
  IntegerVector total_terms (14);
  IntegerVector alt_total_terms (14);
  IntegerVector glm_total_terms (14);
  IntegerVector nocovs_total_terms (14);
  
  // Determine which vital rates need global model formulae
  for (int i = 0; i < nvitalrates; i++) {
    if (vitalrates(i) == "surv") survcheck = 1;
    if (vitalrates(i) == "obs") obscheck = 1;
    if (vitalrates(i) == "size") sizecheck = 1;
    if (vitalrates(i) == "repst") repstcheck = 1;
    if (vitalrates(i) == "fec") feccheck = 1;
  }
  
  String suite_element;
  
  for (int i = 0; i < suite.length(); i++) {
    suite_element = String(suite(i));
    // Tests if inputs are appropriate for suite
    if (suite_element == "full" || suite_element == "main") {
      if (historical) {
        if (sizel != 3 && repstl != 3) {
          if (sizel == 2 && repstl == 2) {
            historical = false;
          } else if (repstl == 3) {
            suite(i) = "rep";
          } else if (sizel == 3) {
            suite(i) = "size";
          }
        }
      } else {
        if (sizel != 2 && repstl != 2) {
          if (sizel != 3 && repstl != 3) {
            suite(i) = "const";
          } else if (repstl > 1 && sizel < 2) {
            suite(i) = "rep";
          } else if (sizel > 1 && repstl < 2) {
            suite(i) = "size";
          }
        }
      }
    } else if (suite_element == "rep") {
      if (historical) {
        if (repstl != 3) {
          if (repstl == 2) {
            historical = false;
          } else {
            suite(i) = "const";
          }
        }
      } else {
        if (repstl != 2) {
          if (repstl != 3) {
            suite(i) = "const";
          }
        }
      }
    } else if (suite_element == "size") {
      if (historical) {
        if (sizel != 3) {
          if (sizel == 2) {
            historical = false;
          } else {
            suite(i) = "const";
          }
        }
      } else {
        if (sizel != 2) {
          if (sizel != 3) {
            suite(i) = "const";
          }
        }
      }
    }
  }
  
  StringVector alt_suite (14);
  
  for (int i = 0; i < 14; i++) {
    if (suite(i) == "full") alt_suite(i) = "main";
    if (suite(i) == "main") alt_suite(i) = "size";
    if (suite(i) == "size") alt_suite(i) = "const";
    if (suite(i) == "rep") alt_suite(i) = "const";
    if (suite(i) == "const") alt_suite(i) = "const";
  }
  
  // Build global models
  if (survcheck) {
    List surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 1, String(suite(0)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(0), ageused(0), densityused(0), indcovaused(0), indcovbused(0),
      indcovcused(0), annucovaused(0), annucovbused(0), annucovcused(0),
      pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
      fectime, repstcheck, interactions);
    
    String sp0_proxy(as<StringVector>(surv_prax[0]));
    fullsurvmodel = sp0_proxy;
    total_terms(0) = static_cast<int>(surv_prax(2));
    
    List alt_surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 1, String(alt_suite(0)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(0), ageused(0), densityused(0), indcovaused(0), indcovbused(0),
      indcovcused(0), annucovaused(0), annucovbused(0), annucovcused(0),
      pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
      fectime, repstcheck, interactions);
    
    String alt_sp0_proxy(as<StringVector>(alt_surv_prax[0]));
    alt_fullsurvmodel = alt_sp0_proxy;
    alt_total_terms(0) = static_cast<int>(alt_surv_prax(2));
    
    if (!nojuvs) {
      String sp1_proxy(as<StringVector>(surv_prax[1]));
      juvsurvmodel = sp1_proxy;
      total_terms(7) = static_cast<int>(surv_prax(3));
      
      String alt_sp1_proxy(as<StringVector>(alt_surv_prax[1]));
      alt_juvsurvmodel = sp1_proxy;
      alt_total_terms(7) = static_cast<int>(alt_surv_prax(3));
      
      List matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 8, String(suite(13)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
        indcovaused(13), indcovbused(13), indcovcused(13), annucovaused(13),
        annucovbused(13), annucovcused(13), pasrand, yasrand, iaasrand,
        ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck,
        interactions);
        
      String mt1_proxy(as<StringVector>(matstat_prax[1]));
      juvmatstmodel = mt1_proxy;
      total_terms(13) = static_cast<int>(matstat_prax(3));
      
      List alt_matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 8, String(alt_suite(13)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
        indcovaused(13), indcovbused(13), indcovcused(13), annucovaused(13),
        annucovbused(13), annucovcused(13), pasrand, yasrand, iaasrand,
        ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck,
        interactions);
        
      String alt_mt1_proxy(as<StringVector>(alt_matstat_prax[1]));
      alt_juvmatstmodel = alt_mt1_proxy;
      alt_total_terms(13) = static_cast<int>(alt_matstat_prax(3));
    }
    
    if (approach == "mixed") {
      List glm_surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 1, String(suite(0)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(0), ageused(0), densityused(0), indcovaused(0),
        indcovbused(0), indcovcused(0), annucovaused(0), annucovbused(0),
        annucovcused(0), false, false, false, false, false, iaasfac, ibasfac,
        icasfac, fectime, repstcheck, interactions);
      
      String glm_sp0_proxy(as<StringVector>(glm_surv_prax[0]));
      glm_fullsurvmodel = glm_sp0_proxy;
      glm_total_terms(0) = static_cast<int>(glm_surv_prax(2));
      
      if (!nojuvs) {
        String glm_sp1_proxy(as<StringVector>(glm_surv_prax[1]));
        glm_juvsurvmodel = glm_sp1_proxy;
        glm_total_terms(7) = static_cast<int>(glm_surv_prax(3));
        
        List glm_matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 8, String(suite(13)), "glm", nojuvs, juvsize,
          indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
          sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
          indcovaused(13), indcovbused(13), indcovcused(13), annucovaused(13),
          annucovbused(13), annucovcused(13), false, false, false, false, false,
          iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
          
        String glm_mt1_proxy(as<StringVector>(glm_matstat_prax[1]));
        glm_juvmatstmodel = glm_mt1_proxy;
        glm_total_terms(13) = static_cast<int>(glm_matstat_prax(3));
      }
    }
    
    int surv_check_indivs = static_cast<int>(indcovaused(0)) + 
      static_cast<int>(indcovbused(0)) + static_cast<int>(indcovcused(0));
    
    if (surv_check_indivs > 0) { 
      List nocovs_surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 1, String(suite(0)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(0), ageused(0), densityused(0), false,
        false, false, false, false, false, pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String nocovs_sp0_proxy(as<StringVector>(nocovs_surv_prax[0]));
      nocovs_fullsurvmodel = nocovs_sp0_proxy;
      nocovs_total_terms(0) = static_cast<int>(nocovs_surv_prax(2));
      
       if (!nojuvs) {
        String nocovs_sp1_proxy(as<StringVector>(nocovs_surv_prax[1]));
        nocovs_juvsurvmodel = nocovs_sp1_proxy;
        nocovs_total_terms(7) = static_cast<int>(nocovs_surv_prax(3));
      }
    }
    
    int matstat_check_indivs = static_cast<int>(indcovaused(13)) + 
      static_cast<int>(indcovbused(13)) + static_cast<int>(indcovcused(13));
    
    if (matstat_check_indivs > 0 && !nojuvs) { 
      List nocovs_matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 8, String(suite(13)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
        false, false, false, false, false, false, pasrand, yasrand, iaasrand,
        ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck,
        interactions);
        
      String nocovs_mt1_proxy(as<StringVector>(nocovs_matstat_prax[1]));
      nocovs_juvmatstmodel = nocovs_mt1_proxy;
      nocovs_total_terms(13) = static_cast<int>(nocovs_matstat_prax(3));
    }
  }
  
  if (obscheck) {
    List obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 2, String(suite(1)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(1), ageused(1), densityused(1), indcovaused(1), indcovbused(1),
      indcovcused(1), annucovaused(1), annucovbused(1), annucovcused(1),
      pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
      fectime, repstcheck, interactions);
    
    String ob0_proxy(as<StringVector>(obs_prax[0]));
    fullobsmodel = ob0_proxy;
    total_terms(1) = static_cast<int>(obs_prax(2));
    
    List alt_obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 2, String(alt_suite(1)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(1), ageused(1), densityused(1), indcovaused(1),
      indcovbused(1), indcovcused(1), annucovaused(1), annucovbused(1),
      annucovcused(1), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck, interactions);
    
    String alt_ob0_proxy(as<StringVector>(alt_obs_prax[0]));
    alt_fullobsmodel = alt_ob0_proxy;
    alt_total_terms(1) = static_cast<int>(alt_obs_prax(2));
    
    if (!nojuvs) {
      String ob1_proxy(as<StringVector>(obs_prax[1]));
      juvobsmodel = ob1_proxy;
      total_terms(8) = static_cast<int>(obs_prax(3));
      
      String alt_ob1_proxy(as<StringVector>(alt_obs_prax[1]));
      alt_juvobsmodel = alt_ob1_proxy;
      alt_total_terms(8) = static_cast<int>(alt_obs_prax(3));
    }
    
    if (approach == "mixed") { 
      List glm_obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 2, String(suite(1)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(1), ageused(1), densityused(1), indcovaused(1),
        indcovbused(1), indcovcused(1), annucovaused(1), annucovbused(1),
        annucovcused(1), false, false, false, false, false, iaasfac, ibasfac,
        icasfac, fectime, repstcheck, interactions);
      
      String glm_ob0_proxy(as<StringVector>(glm_obs_prax[0]));
      glm_fullobsmodel = glm_ob0_proxy;
      glm_total_terms(1) = static_cast<int>(glm_obs_prax(2));
      
      if (!nojuvs) {
        String glm_ob1_proxy(as<StringVector>(glm_obs_prax[1]));
        glm_juvobsmodel = glm_ob1_proxy;
        glm_total_terms(8) = static_cast<int>(glm_obs_prax(3));
      }
    }
    
    int obs_check_indivs = static_cast<int>(indcovaused(1)) + 
      static_cast<int>(indcovbused(1)) + static_cast<int>(indcovcused(1));
    
    if (obs_check_indivs > 0) { 
      List nocovs_obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 2, String(suite(1)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(1), ageused(1), densityused(1), false,
        false, false, false, false, false, pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String nocovs_ob0_proxy(as<StringVector>(nocovs_obs_prax[0]));
      nocovs_fullobsmodel = nocovs_ob0_proxy;
      nocovs_total_terms(1) = static_cast<int>(nocovs_obs_prax(2));
      
       if (!nojuvs) {
        String nocovs_ob1_proxy(as<StringVector>(nocovs_obs_prax[1]));
        nocovs_juvobsmodel = nocovs_ob1_proxy;
        nocovs_total_terms(8) = static_cast<int>(nocovs_obs_prax(3));
      }
    }
  }
  
  if (sizecheck) {
    String current_suite(suite(2));
    
    List size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 3, current_suite, approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(2), ageused(2), densityused(2), indcovaused(2), indcovbused(2),
      indcovcused(2), annucovaused(2), annucovbused(2), annucovcused(2),
      pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
      fectime, repstcheck, interactions);
    
    total_terms(2) = static_cast<int>(size_prax(2));
    
    List alt_size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 3, String(alt_suite(2)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(2), ageused(2), densityused(2), indcovaused(2),
      indcovbused(2), indcovcused(2), annucovaused(2), annucovbused(2),
      annucovcused(2), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck, interactions);
    
    alt_total_terms(2) = static_cast<int>(alt_size_prax(2));
    
    if (stringcompare_hard(current_suite, "full") && total_terms(2) > 18 && size_zero) {
      suite(2) = "main";
      alt_suite(2) = "size";
      
      List size_prax_main = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 3, String(suite(2)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(2), ageused(2), densityused(2),
        indcovaused(2), indcovbused(2), indcovcused(2), annucovaused(2),
        annucovbused(2), annucovcused(2), pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String sz0_proxy(as<StringVector>(size_prax_main[0]));
      fullsizemodel = sz0_proxy;
      total_terms(2) = static_cast<int>(size_prax_main(2));
      
      List alt_size_prax_main = praxis(surv, obs, size, sizeb, sizec, repst,
        fec, matstat, historical, 3, String(alt_suite(2)), approach, nojuvs,
        juvsize, indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(2), ageused(2), densityused(2),
        indcovaused(2), indcovbused(2), indcovcused(2), annucovaused(2),
        annucovbused(2), annucovcused(2), pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String alt_sz0_proxy(as<StringVector>(alt_size_prax_main[0]));
      alt_fullsizemodel = alt_sz0_proxy;
      alt_total_terms(2) = static_cast<int>(alt_size_prax_main(2));
      
      if (!nojuvs && jsize_zero) {
        String sz1_proxy(as<StringVector>(size_prax_main[1]));
        juvsizemodel = sz1_proxy;
        total_terms(9) = static_cast<int>(size_prax_main(3));
        
        String alt_sz1_proxy(as<StringVector>(alt_size_prax_main[1]));
        alt_juvsizemodel = alt_sz1_proxy;
        alt_total_terms(9) = static_cast<int>(alt_size_prax_main(3));
      }
    } else {
      String sz0_proxy(as<StringVector>(size_prax[0]));
      fullsizemodel = sz0_proxy;
      
      String alt_sz0_proxy(as<StringVector>(alt_size_prax[0]));
      alt_fullsizemodel = alt_sz0_proxy;
      
      if (!nojuvs) {
        String sz1_proxy(as<StringVector>(size_prax[1]));
        juvsizemodel = sz1_proxy;
        total_terms(9) = static_cast<int>(size_prax(3));
        
        String alt_sz1_proxy(as<StringVector>(alt_size_prax[1]));
        alt_juvsizemodel = alt_sz1_proxy;
        alt_total_terms(9) = static_cast<int>(alt_size_prax(3));
      }
    }
    
    if (approach == "mixed") {
      List glm_size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 3, current_suite, "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(2), ageused(2), densityused(2), indcovaused(2),
        indcovbused(2), indcovcused(2), annucovaused(2), annucovbused(2),
        annucovcused(2), false, false, false, false, false, iaasfac, ibasfac,
        icasfac, fectime, repstcheck, interactions);
      
      glm_total_terms(2) = static_cast<int>(glm_size_prax(2));
      
      if (!nojuvs) {
        String glm_sz1_proxy(as<StringVector>(glm_size_prax[1]));
        glm_juvsizemodel = glm_sz1_proxy;
        glm_total_terms(9) = static_cast<int>(glm_size_prax(3));
      }
    }
    
    int size_check_indivs = static_cast<int>(indcovaused(2)) + 
      static_cast<int>(indcovbused(2)) + static_cast<int>(indcovcused(2));
    
    if (size_check_indivs > 0) { 
      List nocovs_size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 3, String(suite(2)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(2), ageused(2), densityused(2), false,
        false, false, false, false, false, pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String nocovs_sz0_proxy(as<StringVector>(nocovs_size_prax[0]));
      nocovs_fullsizemodel = nocovs_sz0_proxy;
      nocovs_total_terms(2) = static_cast<int>(nocovs_size_prax(2));
      
       if (!nojuvs) {
        String nocovs_sz1_proxy(as<StringVector>(nocovs_size_prax[1]));
        nocovs_juvsizemodel = nocovs_sz1_proxy;
        nocovs_total_terms(9) = static_cast<int>(nocovs_size_prax(3));
      }
    }
    
    if (sizebused) { 
      String current_suite(suite(3));
      
      List sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 4, current_suite, approach, nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(3), ageused(3), densityused(3), indcovaused(3),
        indcovbused(3), indcovcused(3), annucovaused(3), annucovbused(3),
        annucovcused(3), pasrand, yasrand, iaasrand, ibasrand, icasrand,
        iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      total_terms(3) = static_cast<int>(sizeb_prax(2));
      
      List alt_sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 4, String(alt_suite(3)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(3), ageused(3), densityused(3),
        indcovaused(3), indcovbused(3), indcovcused(3), annucovaused(3),
        annucovbused(3), annucovcused(3), pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      alt_total_terms(3) = static_cast<int>(alt_sizeb_prax(2));
      
      if (stringcompare_hard(current_suite, "full") && total_terms(3) > 18 && sizeb_zero) { 
        suite(3) = "main";
        alt_suite(3) = "size";
        
        List sizeb_prax_main = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 4, String(suite(3)), approach, nojuvs, juvsize,
          indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
          sizebused, sizecused, grouptest(3), ageused(3), densityused(3),
          indcovaused(3), indcovbused(3), indcovcused(3), annucovaused(3),
          annucovbused(3), annucovcused(3), pasrand, yasrand, iaasrand,
          ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck,
          interactions);
        
        String szb0_proxy(as<StringVector>(sizeb_prax_main[0]));
        fullsizebmodel = szb0_proxy;
        total_terms(3) = static_cast<int>(sizeb_prax_main(2));
        
        List alt_sizeb_prax_main = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 4, String(alt_suite(3)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(3), ageused(3),
          densityused(3), indcovaused(3), indcovbused(3), indcovcused(3),
          annucovaused(3), annucovbused(3), annucovcused(3), pasrand, yasrand,
          iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
          repstcheck, interactions);
        
        String alt_szb0_proxy(as<StringVector>(alt_sizeb_prax_main[0]));
        alt_fullsizebmodel = alt_szb0_proxy;
        alt_total_terms(3) = static_cast<int>(alt_sizeb_prax_main(2));
        
        if (!nojuvs && jsizeb_zero) {
          String szb1_proxy(as<StringVector>(sizeb_prax_main[1]));
          juvsizebmodel = szb1_proxy;
          total_terms(10) = static_cast<int>(sizeb_prax_main(3));
          
          String alt_szb1_proxy(as<StringVector>(alt_sizeb_prax_main[1]));
          alt_juvsizebmodel = alt_szb1_proxy;
          alt_total_terms(10) = static_cast<int>(alt_sizeb_prax_main(3));
        }
      } else {
        String szb0_proxy(as<StringVector>(sizeb_prax[0]));
        fullsizebmodel = szb0_proxy;
        
        String alt_szb0_proxy(as<StringVector>(alt_sizeb_prax[0]));
        alt_fullsizebmodel = alt_szb0_proxy;
        
        if (!nojuvs) {
          String szb1_proxy(as<StringVector>(sizeb_prax[1]));
          juvsizebmodel = szb1_proxy;
          total_terms(10) = static_cast<int>(sizeb_prax(3));
          
          String alt_szb1_proxy(as<StringVector>(alt_sizeb_prax[1]));
          alt_juvsizebmodel = alt_szb1_proxy;
          alt_total_terms(10) = static_cast<int>(alt_sizeb_prax(3));
        }
      }
      
      if (approach == "mixed") {
        List glm_sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 4, current_suite, "glm", nojuvs, juvsize, indiv,
          patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
          sizecused, grouptest(3), ageused(3), densityused(3), indcovaused(3),
          indcovbused(3), indcovcused(3), annucovaused(3), annucovbused(3),
          annucovcused(3), false, false, false, false, false, iaasfac, ibasfac,
          icasfac, fectime, repstcheck, interactions);
        
        glm_total_terms(3) = static_cast<int>(glm_sizeb_prax(3));
        
        if (!nojuvs) {
          String glm_szb1_proxy(as<StringVector>(glm_sizeb_prax[1]));
          glm_juvsizebmodel = glm_szb1_proxy;
          glm_total_terms(10) = static_cast<int>(glm_sizeb_prax(3));
        }
      }
      
      int sizeb_check_indivs = static_cast<int>(indcovaused(3)) + 
        static_cast<int>(indcovbused(3)) + static_cast<int>(indcovcused(3));
      
      if (sizeb_check_indivs > 0) { 
        List nocovs_sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 4, String(suite(3)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(3), ageused(3),
          densityused(3), false, false, false, false, false, false, pasrand,
          yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
          fectime, repstcheck, interactions);
        
        String nocovs_szb0_proxy(as<StringVector>(nocovs_sizeb_prax[0]));
        nocovs_fullsizebmodel = nocovs_szb0_proxy;
        nocovs_total_terms(3) = static_cast<int>(nocovs_sizeb_prax(2));
        
         if (!nojuvs) {
          String nocovs_szb1_proxy(as<StringVector>(nocovs_sizeb_prax[1]));
          nocovs_juvsizebmodel = nocovs_szb1_proxy;
          nocovs_total_terms(10) = static_cast<int>(nocovs_sizeb_prax(3));
        }
      }
    }
    
    if (sizecused) { 
      String current_suite(suite(4));
        
      List sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 5, current_suite, approach, nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(4), ageused(4), densityused(4), indcovaused(4),
        indcovbused(4), indcovcused(4), annucovaused(4), annucovbused(4),
        annucovcused(4), pasrand, yasrand, iaasrand, ibasrand, icasrand,
        iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      total_terms(4) = static_cast<int>(sizec_prax(2));
      
      List alt_sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 5, String(alt_suite(4)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(4), ageused(4), densityused(4),
        indcovaused(4), indcovbused(4), indcovcused(4), annucovaused(4),
        annucovbused(4), annucovcused(4), pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      alt_total_terms(4) = static_cast<int>(alt_sizec_prax(2));
      
      if (stringcompare_hard(current_suite, "full") && total_terms(4) > 18 && sizec_zero) {
        suite(4) = "main";
        alt_suite(4) = "size";
      
        List sizec_prax_main = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 5, String(suite(4)), approach, nojuvs, juvsize,
          indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
          sizebused, sizecused, grouptest(4), ageused(4), densityused(4),
          indcovaused(4), indcovbused(4), indcovcused(4), annucovaused(4),
          annucovbused(4), annucovcused(4), pasrand, yasrand, iaasrand,
          ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck,
          interactions);
        
        String szc0_proxy(as<StringVector>(sizec_prax_main[0]));
        fullsizecmodel = szc0_proxy;
        total_terms(4) = static_cast<int>(sizec_prax_main(2));
        
        List alt_sizec_prax_main = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 5, String(alt_suite(4)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(4), ageused(4),
          densityused(4), indcovaused(4), indcovbused(4), indcovcused(4),
          annucovaused(4), annucovbused(4), annucovcused(4), pasrand, yasrand,
          iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
          repstcheck, interactions);
        
        String alt_szc0_proxy(as<StringVector>(alt_sizec_prax_main[0]));
        alt_fullsizecmodel = alt_szc0_proxy;
        alt_total_terms(4) = static_cast<int>(alt_sizec_prax_main(2));
        
        if (!nojuvs && jsizec_zero) {
          String szc1_proxy(as<StringVector>(sizec_prax_main[1]));
          juvsizecmodel = szc1_proxy;
          total_terms(11) = static_cast<int>(sizec_prax_main(3));
          
          String alt_szc1_proxy(as<StringVector>(alt_sizec_prax_main[1]));
          alt_juvsizecmodel = alt_szc1_proxy;
          alt_total_terms(11) = static_cast<int>(alt_sizec_prax_main(3));
        }
      } else {
        String szc0_proxy(as<StringVector>(sizec_prax[0]));
        fullsizecmodel = szc0_proxy;
        
        String alt_szc0_proxy(as<StringVector>(alt_sizec_prax[0]));
        alt_fullsizecmodel = alt_szc0_proxy;
        
        if (!nojuvs) {
          String szc1_proxy(as<StringVector>(sizec_prax[1]));
          juvsizecmodel = szc1_proxy;
          total_terms(11) = static_cast<int>(sizec_prax(3));
          
          String alt_szc1_proxy(as<StringVector>(alt_sizec_prax[1]));
          alt_juvsizecmodel = alt_szc1_proxy;
          alt_total_terms(11) = static_cast<int>(alt_sizec_prax(3));
        }
      }
      
      if (approach == "mixed") {
        List glm_sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 5, current_suite, "glm", nojuvs, juvsize, indiv,
          patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
          sizecused, grouptest(4), ageused(4), densityused(4), indcovaused(4),
          indcovbused(4), indcovcused(4), annucovaused(4), annucovbused(4),
          annucovcused(4), false, false, false, false, false, iaasfac, ibasfac,
          icasfac, fectime, repstcheck, interactions);
        
        glm_total_terms(4) = static_cast<int>(glm_sizec_prax(4));
        
        if (!nojuvs) {
          String glm_szc1_proxy(as<StringVector>(glm_sizec_prax[1]));
          glm_juvsizecmodel = glm_szc1_proxy;
          glm_total_terms(11) = static_cast<int>(glm_sizec_prax(3));
        }
      }
      
      int sizec_check_indivs = static_cast<int>(indcovaused(4)) + 
        static_cast<int>(indcovbused(4)) + static_cast<int>(indcovcused(4));
      
      if (sizec_check_indivs > 0) { 
        List nocovs_sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 5, String(suite(4)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(4), ageused(4),
          densityused(4), false, false, false, false, false, false, pasrand,
          yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
          fectime, repstcheck, interactions);
        
        String nocovs_szc0_proxy(as<StringVector>(nocovs_sizec_prax[0]));
        nocovs_fullsizecmodel = nocovs_szc0_proxy;
        nocovs_total_terms(4) = static_cast<int>(nocovs_sizec_prax(2));
        
         if (!nojuvs) {
          String nocovs_szc1_proxy(as<StringVector>(nocovs_sizec_prax[1]));
          nocovs_juvsizecmodel = nocovs_szc1_proxy;
          nocovs_total_terms(11) = static_cast<int>(nocovs_sizec_prax(3));
        }
      }
    }
  }
  
  if (repstcheck) {
    List repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 6, String(suite(5)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(5), ageused(5), densityused(5), indcovaused(5), indcovbused(5),
      indcovcused(5), annucovaused(5), annucovbused(5), annucovcused(5),
      pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
      fectime, repstcheck, interactions);
    
    String rp0_proxy(as<StringVector>(repst_prax[0]));
    fullrepstmodel = rp0_proxy;
    total_terms(5) = static_cast<int>(repst_prax(2));
    
    List alt_repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 6, String(alt_suite(5)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(5), ageused(5), densityused(5), indcovaused(5),
      indcovbused(5), indcovcused(5), annucovaused(5), annucovbused(5),
      annucovcused(5), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck, interactions);
    
    String alt_rp0_proxy(as<StringVector>(alt_repst_prax[0]));
    alt_fullrepstmodel = alt_rp0_proxy;
    alt_total_terms(5) = static_cast<int>(alt_repst_prax(2));
    
    if (!nojuvs) {
      String rp1_proxy(as<StringVector>(repst_prax[1]));
      juvrepstmodel = rp1_proxy;
      total_terms(12) = static_cast<int>(repst_prax(3));
      
      String alt_rp1_proxy(as<StringVector>(alt_repst_prax[1]));
      alt_juvrepstmodel = alt_rp1_proxy;
      alt_total_terms(12) = static_cast<int>(alt_repst_prax(3));
    }
    
    if (approach == "mixed") { 
      List glm_repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 6, String(suite(5)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(5), ageused(5), densityused(5), indcovaused(5),
        indcovbused(5), indcovcused(5), annucovaused(5), annucovbused(5),
        annucovcused(5), false, false, false, false, false, iaasfac, ibasfac,
        icasfac, fectime, repstcheck, interactions);
      
      String glm_rp0_proxy(as<StringVector>(glm_repst_prax[0]));
      glm_fullrepstmodel = glm_rp0_proxy;
      glm_total_terms(5) = static_cast<int>(glm_repst_prax(2));
      
      if (!nojuvs) {
        String glm_rp1_proxy(as<StringVector>(glm_repst_prax[1]));
        glm_juvrepstmodel = glm_rp1_proxy;
        glm_total_terms(12) = static_cast<int>(glm_repst_prax(3));
      }
    }
    
    int repst_check_indivs = static_cast<int>(indcovaused(5)) + 
      static_cast<int>(indcovbused(5)) + static_cast<int>(indcovcused(5));
    
    if (repst_check_indivs > 0) { 
      List nocovs_repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 6, String(suite(5)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(5), ageused(5), densityused(5), false,
        false, false, false, false, false, pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String nocovs_rp0_proxy(as<StringVector>(nocovs_repst_prax[0]));
      nocovs_fullrepstmodel = nocovs_rp0_proxy;
      nocovs_total_terms(5) = static_cast<int>(nocovs_repst_prax(2));
      
       if (!nojuvs) {
        String nocovs_rp1_proxy(as<StringVector>(nocovs_repst_prax[1]));
        nocovs_juvrepstmodel = nocovs_rp1_proxy;
        nocovs_total_terms(12) = static_cast<int>(nocovs_repst_prax(3));
      }
    }
  }
  
  if (feccheck) {
    List fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 7, String(suite(6)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(6), ageused(6), densityused(6), indcovaused(6), indcovbused(6),
      indcovcused(6), annucovaused(6), annucovbused(6), annucovcused(6),
      pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac,
      fectime, repstcheck, interactions);
    
    String fc0_proxy(as<StringVector>(fec_prax[0]));
    fullfecmodel = fc0_proxy;
    total_terms(6) = static_cast<int>(fec_prax(2));
    
    List alt_fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 7, String(alt_suite(6)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(6), ageused(6), densityused(6), indcovaused(6),
      indcovbused(6), indcovcused(6), annucovaused(6), annucovbused(6),
      annucovcused(6), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck, interactions);
    
    String alt_fc0_proxy(as<StringVector>(alt_fec_prax[0]));
    alt_fullfecmodel = alt_fc0_proxy;
    alt_total_terms(6) = static_cast<int>(alt_fec_prax(2));
    
    if (approach == "mixed") { 
      List glm_fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 7, String(suite(6)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(6), ageused(6), densityused(6), indcovaused(6),
        indcovbused(6), indcovcused(6), annucovaused(6), annucovbused(6),
        annucovcused(6), false, false, false, false, false, iaasfac, ibasfac,
        icasfac, fectime, repstcheck, interactions);
      
      String glm_fc0_proxy(as<StringVector>(glm_fec_prax[0]));
      glm_fullfecmodel = glm_fc0_proxy;
      glm_total_terms(6) = static_cast<int>(glm_fec_prax(2));
    }
    
    int fec_check_indivs = static_cast<int>(indcovaused(6)) + 
      static_cast<int>(indcovbused(6)) + static_cast<int>(indcovcused(6));
    
    if (fec_check_indivs > 0) { 
      List nocovs_fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 7, String(suite(6)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(6), ageused(6), densityused(6), false,
        false, false, false, false, false, pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck, interactions);
      
      String nocovs_fc0_proxy(as<StringVector>(nocovs_fec_prax[0]));
      nocovs_fullfecmodel = nocovs_fc0_proxy;
      nocovs_total_terms(6) = static_cast<int>(nocovs_fec_prax(2));
    }
  }
  
  // New paramnames object
  Rcpp::DataFrame paramnames = paramnames_skeleton(false);
  int pm_varno = static_cast<int>(paramnames.nrows());
  StringVector mainparams = as<StringVector>(paramnames["mainparams"]);
  
  StringVector modelparams (37);
  
  bool annca_used = is_true( any(annucovaused));
  bool anncb_used = is_true( any(annucovbused));
  bool anncc_used = is_true( any(annucovcused));
  
  for (int i = 0; i < pm_varno; i++) {
    modelparams(i) = "none"; // Default value
    
    if (stringcompare_hard(as<std::string>(mainparams(i)), "year2")) modelparams(i) = year;
    if (stringcompare_hard(as<std::string>(mainparams(i)), "individ")) modelparams(i) = indiv;
    if (stringcompare_hard(as<std::string>(mainparams(i)), "patch")) modelparams(i) = patch;
    if (stringcompare_hard(as<std::string>(mainparams(i)), "surv3")) modelparams(i) = surv(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "obs3")) modelparams(i) = obs(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "size3")) modelparams(i) = size(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "repst3")) modelparams(i) = repst(0);
    if (fectime == 3 && stringcompare_hard(as<std::string>(mainparams(i)), "fec3")) modelparams(i) = fec(0);
    if (fectime == 2 && stringcompare_hard(as<std::string>(mainparams(i)), "fec2")) modelparams(i) = fec(1);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "size2")) modelparams(i) = size(1);
    
    if (sizebused) {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb3")) modelparams(i) = sizeb(0);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb2")) modelparams(i) = sizeb(1);
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizeb1")) modelparams(i) = sizeb(2);
      
    } else {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb3")) modelparams(i) = "none";
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb2")) modelparams(i) = "none";
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizeb1")) modelparams(i) = "none";
    }
    
    if (sizecused) {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec3")) modelparams(i) = sizec(0);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec2")) modelparams(i) = sizec(1);
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizec1")) modelparams(i) = sizec(2);
      
    } else {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec3")) modelparams(i) = "none";
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec2")) modelparams(i) = "none";
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizec1")) modelparams(i) = "none";
    }
    
    if (stringcompare_hard(as<std::string>(mainparams(i)), "repst2")) modelparams(i) = repst(1);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "matst3")) modelparams(i) = matstat(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "matstat2")) modelparams(i) = matstat(1);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "age")) modelparams(i) = age;
    
    if (densityused && stringcompare_hard(as<std::string>(mainparams(i)), "density")) modelparams(i) = densitycol;
    if (grouptest && stringcompare_hard(as<std::string>(mainparams(i)), "group2")) modelparams(i) = "group2";
    if (indcovaused && stringcompare_hard(as<std::string>(mainparams(i)), "indcova2")) modelparams(i) = indcova(1);
    if (indcovbused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovb2")) modelparams(i) = indcovb(1);
    if (indcovcused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovc2")) modelparams(i) = indcovc(1);
    
    if (annca_used && stringcompare_hard(as<std::string>(mainparams(i)), "annucova2")) modelparams(i) = "annucova2";
    if (anncb_used && stringcompare_hard(as<std::string>(mainparams(i)), "annucovb2")) modelparams(i) = "annucovb2";
    if (anncc_used && stringcompare_hard(as<std::string>(mainparams(i)), "annucovc2")) modelparams(i) = "annucovc2";
    
    // Further historical terms, if used
    if (historical) {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "size1")) modelparams(i) = size(2);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "repst1")) modelparams(i) = repst(2);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "group1")) modelparams(i) = "group1";
      
      if (indcovaused && stringcompare_hard(as<std::string>(mainparams(i)), "indcova1")) modelparams(i) = indcova(2);
      if (indcovbused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovb1")) modelparams(i) = indcovb(2);
      if (indcovcused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovc1")) modelparams(i) = indcovc(2);
      
      if (annca_used && stringcompare_hard(as<std::string>(mainparams(i)), "annucova1")) modelparams(i) = "annucova1";
      if (anncb_used && stringcompare_hard(as<std::string>(mainparams(i)), "annucovb1")) modelparams(i) = "annucovb1";
      if (anncc_used && stringcompare_hard(as<std::string>(mainparams(i)), "annucovc1")) modelparams(i) = "annucovc1";
    }
  }
  paramnames["modelparams"] = modelparams;
  
  Rcpp::List main = List::create(Named("full.surv.model") = fullsurvmodel,
    _["full.obs.model"] = fullobsmodel, _["full.size.model"] = fullsizemodel,
    _["full.sizeb.model"] = fullsizebmodel, _["full.sizec.model"] = fullsizecmodel,
    _["full.repst.model"] = fullrepstmodel, _["full.fec.model"] = fullfecmodel,
    _["juv.surv.model"] = juvsurvmodel, _["juv.obs.model"] = juvobsmodel,
    _["juv.size.model"] = juvsizemodel, _["juv.sizeb.model"] = juvsizebmodel,
    _["juv.sizec.model"] = juvsizecmodel, _["juv.repst.model"] = juvrepstmodel,
    _["juv.matst.model"] = juvmatstmodel, _["paramnames"] = paramnames,
    _["total_terms"] = total_terms);
  
  Rcpp::List alternate = List::create(Named("full.surv.model") = alt_fullsurvmodel,
    _["full.obs.model"] = alt_fullobsmodel, _["full.size.model"] = alt_fullsizemodel,
    _["full.sizeb.model"] = alt_fullsizebmodel, _["full.sizec.model"] = alt_fullsizecmodel,
    _["full.repst.model"] = alt_fullrepstmodel, _["full.fec.model"] = alt_fullfecmodel,
    _["juv.surv.model"] = alt_juvsurvmodel, _["juv.obs.model"] = alt_juvobsmodel,
    _["juv.size.model"] = alt_juvsizemodel, _["juv.sizeb.model"] = alt_juvsizebmodel,
    _["juv.sizec.model"] = alt_juvsizecmodel, _["juv.repst.model"] = alt_juvrepstmodel,
    _["juv.matst.model"] = alt_juvmatstmodel, _["total_terms"] = alt_total_terms);
  
  Rcpp::List glm_alternate = List::create(Named("full.surv.model") = glm_fullsurvmodel,
    _["full.obs.model"] = glm_fullobsmodel, _["full.size.model"] = glm_fullsizemodel,
    _["full.sizeb.model"] = glm_fullsizebmodel, _["full.sizec.model"] = glm_fullsizecmodel,
    _["full.repst.model"] = glm_fullrepstmodel, _["full.fec.model"] = glm_fullfecmodel,
    _["juv.surv.model"] = glm_juvsurvmodel, _["juv.obs.model"] = glm_juvobsmodel,
    _["juv.size.model"] = glm_juvsizemodel, _["juv.sizeb.model"] = glm_juvsizebmodel,
    _["juv.sizec.model"] = glm_juvsizecmodel, _["juv.repst.model"] = glm_juvrepstmodel,
    _["juv.matst.model"] = glm_juvmatstmodel, _["total_terms"] = glm_total_terms);
  
  Rcpp::List nocovs_alternate = List::create(Named("full.surv.model") = nocovs_fullsurvmodel,
    _["full.obs.model"] = nocovs_fullobsmodel, _["full.size.model"] = nocovs_fullsizemodel,
    _["full.sizeb.model"] = nocovs_fullsizebmodel, _["full.sizec.model"] = nocovs_fullsizecmodel,
    _["full.repst.model"] = nocovs_fullrepstmodel, _["full.fec.model"] = nocovs_fullfecmodel,
    _["juv.surv.model"] = nocovs_juvsurvmodel, _["juv.obs.model"] = nocovs_juvobsmodel,
    _["juv.size.model"] = nocovs_juvsizemodel, _["juv.sizeb.model"] = nocovs_juvsizebmodel,
    _["juv.sizec.model"] = nocovs_juvsizecmodel, _["juv.repst.model"] = nocovs_juvrepstmodel,
    _["juv.matst.model"] = nocovs_juvmatstmodel, _["total_terms"] = nocovs_total_terms);
  
  List output = List::create(_["main"] = main, _["alternate"] = alternate,
    _["glm.alternate"] = glm_alternate, _["nocovs.alternate"] = nocovs_alternate);
  
  return output;
}

//' Creates a Skeleton Paramnames Object for Use in Function-based Modeling
//' 
//' Creates a simple skeleton \code{paramnames} object that can be entered as
//' input in functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
//' \code{\link{aflefko2}()}.
//' 
//' @name create_pm
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
//' @examples 
//' our_pm <- create_pm()
//' our_pm
//' 
//' @export
// [[Rcpp::export(create_pm)]]
Rcpp::DataFrame create_pm(bool name_terms = false) {
  
  DataFrame output = paramnames_skeleton(name_terms);
    
  return output;
}

//' Convert modelextract Coefficient Vector to vrm_frame Vector
//' 
//' @name vrmf_inator
//' 
//' This function reorders the coefficient vector developed by function
//' \code{modelextract} in the \code{lefko3} C++ header file \code{main_utils.h}
//' to match the order required to create the \code{vrm_frame} element within
//' a new \code{vrm_input} object.
//' 
//' @param coef_vec The main \code{coefficients} vector from the model proxy
//' object used.
//' @param zi A logical value indicating whether the coefficients to be handled
//' are from the zero process model in a zero-inflated size or fecundity model.
//' 
//' @return A NumericVector with all coefficients in the proper order for a
//' \code{vrm_frame} data frame.
//' 
//' @keywords internal
//' @noRd
Rcpp::NumericVector vrmf_inator (NumericVector coef_vec, bool zi = false) {
  
  NumericVector main_out (128);
  
  if (!zi) {
    main_out(0) = coef_vec(0); // y-intercept
    main_out(1) = coef_vec(4); // size2
    main_out(2) = coef_vec(3); // size1
    main_out(3) = coef_vec(100); // sizeb2
    main_out(4) = coef_vec(101); // sizeb1
    main_out(5) = coef_vec(102); // sizec2
    main_out(6) = coef_vec(103); // sizec1
    main_out(7) = coef_vec(2); // repst2
    main_out(8) = coef_vec(1); // repst1
    main_out(9) = coef_vec(11); // age
    main_out(10) = coef_vec(104); // density
    main_out(11) = coef_vec(16); // indcova2
    main_out(12) = coef_vec(19); // indcova1
    main_out(13) = coef_vec(17); // indcovb2
    main_out(14) = coef_vec(20); // indcovb1
    main_out(15) = coef_vec(18); // indcovc2
    main_out(16) = coef_vec(21); // indcovc1
    main_out(17) = coef_vec(5); // repst1 repst2
    main_out(18) = coef_vec(6); // size1 size2
    main_out(19) = coef_vec(7); // size1 repst1
    main_out(20) = coef_vec(8); // size2 repst2
    main_out(21) = coef_vec(9); // size2 repst1
    main_out(22) = coef_vec(10); // size1 repst2
    main_out(23) = coef_vec(12); // age size1
    main_out(24) = coef_vec(13); // age size2
    main_out(25) = coef_vec(14); // age repst1
    main_out(26) = coef_vec(15); // age repst2
    main_out(27) = coef_vec(22); // indcova2 size2
    main_out(28) = coef_vec(23); // indcovb2 size2
    main_out(29) = coef_vec(24); // indcovc2 size2
    main_out(30) = coef_vec(25); // indcova2 repst2
    main_out(31) = coef_vec(26); // indcovb2 repst2
    main_out(32) = coef_vec(27); // indcovc2 repst2
    main_out(33) = coef_vec(28); // indcova1 size1
    main_out(34) = coef_vec(29); // indcovb1 size1
    main_out(35) = coef_vec(30); // indcovc1 size1
    main_out(36) = coef_vec(31); // indcova1 repst1
    main_out(37) = coef_vec(32); // indcovb1 repst1
    main_out(38) = coef_vec(33); // indcovc1 repst1
    main_out(39) = coef_vec(34); // indcova2 indcovb2
    main_out(40) = coef_vec(35); // indcova2 indcovc2
    main_out(41) = coef_vec(36); // indcovb2 indcovc2
    main_out(42) = coef_vec(37); // indcova1 indcovb1
    main_out(43) = coef_vec(38); // indcova1 indcovc1
    main_out(44) = coef_vec(39); // indcovb1 indcovc1
    main_out(45) = coef_vec(40); // indcova2 indcovb1
    main_out(46) = coef_vec(41); // indcova1 indcovb2
    main_out(47) = coef_vec(42); // indcova2 indcovc1
    main_out(48) = coef_vec(43); // indcova1 indcovc2
    main_out(49) = coef_vec(44); // indcovb2 indcovc1
    main_out(50) = coef_vec(45); // indcovb1 indcovc2
    main_out(51) = coef_vec(105); // sizeb2 sizeb1
    main_out(52) = coef_vec(106); // sizec2 sizec1
    main_out(53) = coef_vec(107); // size1 sizeb1
    main_out(54) = coef_vec(108); // size1 sizec1
    main_out(55) = coef_vec(109); // sizeb1 sizec1
    main_out(56) = coef_vec(110); // size2 sizeb2
    main_out(57) = coef_vec(111); // size2 sizec2
    main_out(58) = coef_vec(112); // sizeb2 sizec2
    main_out(59) = coef_vec(113); // size1 sizeb2
    main_out(60) = coef_vec(114); // size1 sizec2
    main_out(61) = coef_vec(115); // sizeb1 sizec2
    main_out(62) = coef_vec(116); // size2 sizeb1
    main_out(63) = coef_vec(117); // size2 sizec1
    main_out(64) = coef_vec(118); // sizeb2 sizec1
    main_out(65) = coef_vec(119); // density size2
    main_out(66) = coef_vec(120); // density sizeb2
    main_out(67) = coef_vec(121); // density sizec2
    main_out(68) = coef_vec(122); // density size1
    main_out(69) = coef_vec(123); // density sizeb1
    main_out(70) = coef_vec(124); // density sizec1
    main_out(71) = coef_vec(125); // density repst2
    main_out(72) = coef_vec(126); // density repst1
    main_out(73) = coef_vec(127); // sizeb2 repst2
    main_out(74) = coef_vec(128); // sizec2 repst2
    main_out(75) = coef_vec(130); // sizeb1 repst1
    main_out(76) = coef_vec(131); // sizeb2 repst1
    main_out(77) = coef_vec(132); // sizeb1 repst2
    main_out(78) = coef_vec(133); // sizec1 repst1
    main_out(79) = coef_vec(134); // sizec2 repst1
    main_out(80) = coef_vec(135); // sizec1 repst2
    main_out(81) = coef_vec(136); // sizeb2 age
    main_out(82) = coef_vec(137); // sizec2 age
    main_out(83) = coef_vec(138); // density age
    main_out(84) = coef_vec(139); // sizeb1 age
    main_out(85) = coef_vec(140); // sizec1 age
    main_out(86) = coef_vec(141); // indcova2 sizeb2
    main_out(87) = coef_vec(142); // indcova2 sizec2
    main_out(88) = coef_vec(143); // indcova2 density
    main_out(89) = coef_vec(144); // indcova1 sizeb1
    main_out(90) = coef_vec(145); // indcova1 sizec1
    main_out(91) = coef_vec(146); // indcova1 sizeb2
    main_out(92) = coef_vec(147); // indcova1 sizec2
    main_out(93) = coef_vec(148); // indcova2 sizeb1
    main_out(94) = coef_vec(149); // indcova2 sizec1
    main_out(95) = coef_vec(150); // indcova1 density
    main_out(96) = coef_vec(151); // indcovb2 sizeb2
    main_out(97) = coef_vec(152); // indcovb2 sizec2
    main_out(98) = coef_vec(153); // indcovb2 density
    main_out(99) = coef_vec(154); // indcovb1 sizeb1
    main_out(100) = coef_vec(155); // indcovb1 sizec1
    main_out(101) = coef_vec(156); // indcovb1 sizeb2
    main_out(102) = coef_vec(157); // indcovb1 sizec2
    main_out(103) = coef_vec(158); // indcovb2 sizeb1
    main_out(104) = coef_vec(159); // indcovb2 sizec1
    main_out(105) = coef_vec(160); // indcovb1 density
    main_out(106) = coef_vec(161); // indcovc2 sizeb2
    main_out(107) = coef_vec(162); // indcovc2 sizec2
    main_out(108) = coef_vec(163); // indcovc2 density
    main_out(109) = coef_vec(164); // indcovc1 sizeb1
    main_out(110) = coef_vec(165); // indcovc1 sizec1
    main_out(111) = coef_vec(166); // indcovc1 sizeb2
    main_out(112) = coef_vec(167); // indcovc1 sizec2
    main_out(113) = coef_vec(168); // indcovc2 sizeb1
    main_out(114) = coef_vec(169); // indcovc2 sizec1
    main_out(115) = coef_vec(170); // indcovc1 density
    main_out(116) = coef_vec(171); // indcova2 size1
    main_out(117) = coef_vec(172); // indcovb2 size1
    main_out(118) = coef_vec(173); // indcovc2 size1
    main_out(119) = coef_vec(174); // indcova1 size2
    main_out(120) = coef_vec(175); // indcovb1 size2
    main_out(121) = coef_vec(176); // indcovc1 size2
    main_out(122) = coef_vec(177); // indcova2 repst1
    main_out(123) = coef_vec(178); // indcovb2 repst1
    main_out(124) = coef_vec(179); // indcovc2 repst1
    main_out(125) = coef_vec(180); // indcova1 repst2
    main_out(126) = coef_vec(181); // indcovb1 repst2
    main_out(127) = coef_vec(182); // indcovc1 repst2
    
  } else {
    
    main_out(0) = coef_vec(46); // y-intercept
    main_out(1) = coef_vec(50); // size2
    main_out(2) = coef_vec(49); // size1
    main_out(3) = coef_vec(200); // sizeb2
    main_out(4) = coef_vec(201); // sizeb1
    main_out(5) = coef_vec(202); // sizec2
    main_out(6) = coef_vec(203); // sizec1
    main_out(7) = coef_vec(48); // repst2
    main_out(8) = coef_vec(47); // repst1
    main_out(9) = coef_vec(57); // age
    main_out(10) = coef_vec(204); // density
    main_out(11) = coef_vec(62); // indcova2
    main_out(12) = coef_vec(65); // indcova1
    main_out(13) = coef_vec(63); // indcovb2
    main_out(14) = coef_vec(66); // indcovb1
    main_out(15) = coef_vec(64); // indcovc2
    main_out(16) = coef_vec(67); // indcovc1
    main_out(17) = coef_vec(51); // repst1 repst2
    main_out(18) = coef_vec(52); // size1 size2
    main_out(19) = coef_vec(53); // size1 repst1
    main_out(20) = coef_vec(54); // size2 repst2
    main_out(21) = coef_vec(55); // size2 repst1
    main_out(22) = coef_vec(56); // size1 repst2
    main_out(23) = coef_vec(58); // age size1
    main_out(24) = coef_vec(59); // age size2
    main_out(25) = coef_vec(60); // age repst1
    main_out(26) = coef_vec(61); // age repst2
    main_out(27) = coef_vec(68); // indcova2 size2
    main_out(28) = coef_vec(69); // indcovb2 size2
    main_out(29) = coef_vec(70); // indcovc2 size2
    main_out(30) = coef_vec(71); // indcova2 repst2
    main_out(31) = coef_vec(72); // indcovb2 repst2
    main_out(32) = coef_vec(73); // indcovc2 repst2
    main_out(33) = coef_vec(74); // indcova1 size1
    main_out(34) = coef_vec(75); // indcovb1 size1
    main_out(35) = coef_vec(76); // indcovc1 size1
    main_out(36) = coef_vec(77); // indcova1 repst1
    main_out(37) = coef_vec(78); // indcovb1 repst1
    main_out(38) = coef_vec(79); // indcovc1 repst1
    main_out(39) = coef_vec(80); // indcova2 indcovb2
    main_out(40) = coef_vec(81); // indcova2 indcovc2
    main_out(41) = coef_vec(82); // indcovb2 indcovc2
    main_out(42) = coef_vec(83); // indcova1 indcovb1
    main_out(43) = coef_vec(84); // indcova1 indcovc1
    main_out(44) = coef_vec(85); // indcovb1 indcovc1
    main_out(45) = coef_vec(86); // indcova2 indcovb1
    main_out(46) = coef_vec(87); // indcova1 indcovb2
    main_out(47) = coef_vec(88); // indcova2 indcovc1
    main_out(48) = coef_vec(89); // indcova1 indcovc2
    main_out(49) = coef_vec(90); // indcovb2 indcovc1
    main_out(50) = coef_vec(91); // indcovb1 indcovc2
    main_out(51) = coef_vec(205); // sizeb2 sizeb1
    main_out(52) = coef_vec(206); // sizec2 sizec1
    main_out(53) = coef_vec(207); // size1 sizeb1
    main_out(54) = coef_vec(208); // size1 sizec1
    main_out(55) = coef_vec(209); // sizeb1 sizec1
    main_out(56) = coef_vec(210); // size2 sizeb2
    main_out(57) = coef_vec(211); // size2 sizec2
    main_out(58) = coef_vec(212); // sizeb2 sizec2
    main_out(59) = coef_vec(213); // size1 sizeb2
    main_out(60) = coef_vec(214); // size1 sizec2
    main_out(61) = coef_vec(215); // sizeb1 sizec2
    main_out(62) = coef_vec(216); // size2 sizeb1
    main_out(63) = coef_vec(217); // size2 sizec1
    main_out(64) = coef_vec(218); // sizeb2 sizec1
    main_out(65) = coef_vec(219); // density size2
    main_out(66) = coef_vec(220); // density sizeb2
    main_out(67) = coef_vec(221); // density sizec2
    main_out(68) = coef_vec(222); // density size1
    main_out(69) = coef_vec(223); // density sizeb1
    main_out(70) = coef_vec(224); // density sizec1
    main_out(71) = coef_vec(225); // density repst2
    main_out(72) = coef_vec(226); // density repst1
    main_out(73) = coef_vec(227); // sizeb2 repst2
    main_out(74) = coef_vec(228); // sizec2 repst2
    main_out(75) = coef_vec(230); // sizeb1 repst1
    main_out(76) = coef_vec(231); // sizeb2 repst1
    main_out(77) = coef_vec(232); // sizeb1 repst2
    main_out(78) = coef_vec(233); // sizec1 repst1
    main_out(79) = coef_vec(234); // sizec2 repst1
    main_out(80) = coef_vec(235); // sizec1 repst2
    main_out(81) = coef_vec(236); // sizeb2 age
    main_out(82) = coef_vec(237); // sizec2 age
    main_out(83) = coef_vec(238); // density age
    main_out(84) = coef_vec(239); // sizeb1 age
    main_out(85) = coef_vec(240); // sizec1 age
    main_out(86) = coef_vec(241); // indcova2 sizeb2
    main_out(87) = coef_vec(242); // indcova2 sizec2
    main_out(88) = coef_vec(243); // indcova2 density
    main_out(89) = coef_vec(244); // indcova1 sizeb1
    main_out(90) = coef_vec(245); // indcova1 sizec1
    main_out(91) = coef_vec(246); // indcova1 sizeb2
    main_out(92) = coef_vec(247); // indcova1 sizec2
    main_out(93) = coef_vec(248); // indcova2 sizeb1
    main_out(94) = coef_vec(249); // indcova2 sizec1
    main_out(95) = coef_vec(250); // indcova1 density
    main_out(96) = coef_vec(251); // indcovb2 sizeb2
    main_out(97) = coef_vec(252); // indcovb2 sizec2
    main_out(98) = coef_vec(253); // indcovb2 density
    main_out(99) = coef_vec(254); // indcovb1 sizeb1
    main_out(100) = coef_vec(255); // indcovb1 sizec1
    main_out(101) = coef_vec(256); // indcovb1 sizeb2
    main_out(102) = coef_vec(257); // indcovb1 sizec2
    main_out(103) = coef_vec(258); // indcovb2 sizeb1
    main_out(104) = coef_vec(259); // indcovb2 sizec1
    main_out(105) = coef_vec(260); // indcovb1 density
    main_out(106) = coef_vec(261); // indcovc2 sizeb2
    main_out(107) = coef_vec(262); // indcovc2 sizec2
    main_out(108) = coef_vec(263); // indcovc2 density
    main_out(109) = coef_vec(264); // indcovc1 sizeb1
    main_out(110) = coef_vec(265); // indcovc1 sizec1
    main_out(111) = coef_vec(266); // indcovc1 sizeb2
    main_out(112) = coef_vec(267); // indcovc1 sizec2
    main_out(113) = coef_vec(268); // indcovc2 sizeb1
    main_out(114) = coef_vec(269); // indcovc2 sizec1
    main_out(115) = coef_vec(270); // indcovc1 density
    main_out(116) = coef_vec(271); // indcova2 size1
    main_out(117) = coef_vec(272); // indcovb2 size1
    main_out(118) = coef_vec(273); // indcovc2 size1
    main_out(119) = coef_vec(274); // indcova1 size2
    main_out(120) = coef_vec(275); // indcovb1 size2
    main_out(121) = coef_vec(276); // indcovc1 size2
    main_out(122) = coef_vec(277); // indcova2 repst1
    main_out(123) = coef_vec(278); // indcovb2 repst1
    main_out(124) = coef_vec(279); // indcovc2 repst1
    main_out(125) = coef_vec(280); // indcova1 repst2
    main_out(126) = coef_vec(281); // indcovb1 repst2
    main_out(127) = coef_vec(282); // indcovc1 repst2
  }
 
  return main_out;
}

//' Minimize lefkoMod Object by Conversion to vrm_input Object
//' 
//' This function takes a \code{lefkoMod} object, which consists of vital rate
//' models, their associated \code{dredge} model tables, and related metadata,
//' and converts them to minimal data frame lists useable in MPM creation and
//' projection. The main advantage to using this approach is in memory savings.
//' 
//' @name miniMod
//' 
//' @param lMod A \code{lefkoMod} object.
//' @param hfv_data The \code{hfv_data} formatted data frame used to develop
//' object \code{lMod}.
//' @param stageframe The stageframe used to develop object \code{lMod}.
//' @param all_years A vector giving the times / years used to develop object
//' \code{lMod}, exactly as used in the latter. Only needed if object
//' \code{hfv_data} not provided.
//' @param all_patches A vector giving the patch names used to develop object
//' \code{lMod}, exactly as used in the latter. Only needed if object
//' \code{hfv_data} not provided.
//' @param all_groups A vector giving the stage groups used to develop object
//' \code{lMod}, exactly as used in the latter. Only needed if object
//' \code{stageframe} not provided.
//' @param all_indcova The name of individual covariate a if quantitative and
//' non-categorical, or of the categories used if the covariate is a factor
//' variable. Only needed if object \code{hfv_data} not provided but individual
//' covariates used in vital rate models.
//' @param all_indcovb The name of individual covariate a if quantitative and
//' non-categorical, or of the categories used if the cvoariate is a factor
//' variable. Only needed if object \code{hfv_data} not provided but individual
//' covariates used in vital rate models.
//' @param all_indcovc The name of individual covariate a if quantitative and
//' non-categorical, or of the categories used if the covariate is a factor
//' variable. Only needed if object \code{hfv_data} not provided but individual
//' covariates used in vital rate models.
//' 
//' @return An object of class \code{vrm_input}. See function
//' \code{\link{vrm_import}()} for details.
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
//' lathmodelsln3 <- modelsearch(lathvertln, historical = TRUE, 
//'   approach = "mixed", suite = "main", 
//'   vitalrates = c("surv", "obs", "size", "repst", "fec"), juvestimate = "Sdl",
//'   bestfit = "AICc&k", sizedist = "gaussian", fecdist = "poisson", 
//'   indiv = "individ", patch = "patchid", year = "year2",
//'   year.as.random = TRUE, patch.as.random = TRUE, show.model.tables = TRUE,
//'   quiet = "partial")
//' 
//' lathmodels_mini <- miniMod(lathmodelsln3, hfv_data = lathvertln,
//'   stageframe = lathframeln)
//' lathmodels_mini
//' }
//' 
//' @export miniMod
// [[Rcpp::export(miniMod)]]
Rcpp::List miniMod (RObject lMod, Nullable<RObject> hfv_data = R_NilValue,
  Nullable<RObject> stageframe = R_NilValue, Nullable<RObject> all_years  = R_NilValue,
  Nullable<RObject> all_patches  = R_NilValue, Nullable<RObject> all_groups  = R_NilValue,
  Nullable<RObject> all_indcova  = R_NilValue, Nullable<RObject> all_indcovb  = R_NilValue,
  Nullable<RObject> all_indcovc  = R_NilValue) {
  
  String ms_class;
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
  
  DataFrame hfv_main;
  DataFrame sf_main;
  CharacterVector mainyears;
  CharacterVector mainpatches;
  CharacterVector maingroups;
  CharacterVector mainindcova;
  CharacterVector mainindcovb;
  CharacterVector mainindcovc;
  
  if (hfv_data.isNotNull()) { 
    if (is<DataFrame>(hfv_data)) { 
      hfv_main = as<DataFrame>(hfv_data);
      
      IntegerVector year2_vec = as<IntegerVector>(hfv_main["year2"]);
      IntegerVector year2_unique = sort_unique(year2_vec);
      
      mainyears = as<CharacterVector>(year2_unique);
      
      CharacterVector patch_vec = as<CharacterVector>(hfv_main["patchid"]);
      CharacterVector patch_unique = sort_unique(patch_vec);
      
      mainpatches = patch_unique;
      
      CharacterVector all_hfv_vars = as<CharacterVector>(hfv_main.attr("names"));
      
      IntegerVector ica2_indices = index_l3(all_hfv_vars, "indcova2");
      if (static_cast<int>(ica2_indices.length() > 0)) {
        if (is<NumericVector>(hfv_main[static_cast<int>(ica2_indices(0))]) ||
          is<IntegerVector>(hfv_main[static_cast<int>(ica2_indices(0))])) {
            NumericVector toasted_almonds_a2 = as<NumericVector>(hfv_main[static_cast<int>(ica2_indices(0))]);
            if (toasted_almonds_a2.hasAttribute("levels")) {
              mainindcova = as<CharacterVector>(toasted_almonds_a2.attr("levels"));
            } else mainindcova = {NA_STRING};
          
        } else if (is<CharacterVector>(hfv_main[static_cast<int>(ica2_indices(0))])) {
          CharacterVector found_indcova2 = as<CharacterVector>(hfv_main[static_cast<int>(ica2_indices(0))]);
          mainindcova = sort_unique(found_indcova2);
        }
      }
      IntegerVector icb2_indices = index_l3(all_hfv_vars, "indcovb2");
      if (static_cast<int>(icb2_indices.length() > 0)) {
        if (is<NumericVector>(hfv_main[static_cast<int>(icb2_indices(0))]) ||
          is<IntegerVector>(hfv_main[static_cast<int>(icb2_indices(0))])) {
            NumericVector toasted_almonds_b2 = as<NumericVector>(hfv_main[static_cast<int>(icb2_indices(0))]);
            if (toasted_almonds_b2.hasAttribute("levels")) {
              mainindcovb = as<CharacterVector>(toasted_almonds_b2.attr("levels"));
            } else mainindcovb = {NA_STRING};
          
        } else if (is<CharacterVector>(hfv_main[static_cast<int>(icb2_indices(0))])) {
          CharacterVector found_indcovb2 = as<CharacterVector>(hfv_main[static_cast<int>(icb2_indices(0))]);
          mainindcovb = sort_unique(found_indcovb2);
        }
      }
      IntegerVector icc2_indices = index_l3(all_hfv_vars, "indcovc2");
      if (static_cast<int>(icc2_indices.length() > 0)) {
        if (is<NumericVector>(hfv_main[static_cast<int>(icc2_indices(0))]) ||
          is<IntegerVector>(hfv_main[static_cast<int>(icc2_indices(0))])) {
            NumericVector toasted_almonds_c2 = as<NumericVector>(hfv_main[static_cast<int>(icc2_indices(0))]);
            if (toasted_almonds_c2.hasAttribute("levels")) {
              mainindcovc = as<CharacterVector>(toasted_almonds_c2.attr("levels"));
            } else mainindcovc = {NA_STRING};
          
        } else if (is<CharacterVector>(hfv_main[static_cast<int>(icc2_indices(0))])) {
          CharacterVector found_indcovc2 = as<CharacterVector>(hfv_main[static_cast<int>(icc2_indices(0))]);
          mainindcovc = sort_unique(found_indcovc2);
        }
      }
      
    }
  } else { 
    if (all_years.isNotNull()) { 
      mainyears = as<CharacterVector>(all_years);
    } else {
      throw Rcpp::exception("Argument all_years is required.", false);
    }
    
    if (all_patches.isNotNull()) { 
      mainpatches = as<CharacterVector>(all_patches);
    } else { 
      mainpatches = {NA_STRING};
    }
  }
  
  if (stageframe.isNotNull()) {
    if (is<DataFrame>(stageframe)) {
      sf_main = as<DataFrame>(stageframe);
      
      IntegerVector groups_vec = as<IntegerVector>(sf_main["group"]);
      IntegerVector groups_unique = sort_unique(groups_vec);
      
      maingroups = as<CharacterVector>(groups_unique);
    }
  } else {
    if (all_groups.isNotNull()) {
      maingroups = as<CharacterVector>(all_groups);
    } else {
      maingroups = {"0"};
    }
  }
  
  if (all_indcova.isNotNull() && !hfv_data.isNotNull()) {
    mainindcova = as<CharacterVector>(all_indcova);
  } else {
    mainindcova = {"0"};
  }
  
  if (all_indcovb.isNotNull() && !hfv_data.isNotNull()) {
    mainindcovb = as<CharacterVector>(all_indcovb);
  } else {
    mainindcovb = {"0"};
  }
  
  if (all_indcovc.isNotNull() && !hfv_data.isNotNull()) {
    mainindcovc = as<CharacterVector>(all_indcovc);
  } else {
    mainindcovc = {"0"};
  }
  
  if (is<List>(lMod)) {
    Rcpp::List ms_intermediate(lMod);
    
    if (ms_intermediate.hasAttribute("class")) {
      CharacterVector ms_class_vector = ms_intermediate.attr("class");
      ms_class = ms_class_vector(0);
    } else {
      throw Rcpp::exception("Object lMod is List of unknown class.", false);
    }
    
    if (stringcompare_simple(ms_class, "lefkoMo", false)) {
      surmodl = ms_intermediate["survival_model"];
      obsmodl = ms_intermediate["observation_model"];
      sizmodl = ms_intermediate["size_model"];
      sibmodl = ms_intermediate["sizeb_model"];
      sicmodl = ms_intermediate["sizec_model"];
      repmodl = ms_intermediate["repstatus_model"];
      fecmodl = ms_intermediate["fecundity_model"];
      jsurmodl = ms_intermediate["juv_survival_model"];
      jobsmodl = ms_intermediate["juv_observation_model"];
      jsizmodl = ms_intermediate["juv_size_model"];
      jsibmodl = ms_intermediate["juv_sizeb_model"];
      jsicmodl = ms_intermediate["juv_sizec_model"];
      jrepmodl = ms_intermediate["juv_reproduction_model"];
      jmatmodl = ms_intermediate["juv_maturity_model"];
      
      pmnames = as<DataFrame>(ms_intermediate["paramnames"]);
    } else {
      throw Rcpp::exception("Object lMod is List of unknown class.", false);
    }
    
  } else {
    throw Rcpp::exception("Object lMod is object of unknown class.", false);
  }
  
  List surv_proxy = modelextract(surmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List obs_proxy = modelextract(obsmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List size_proxy = modelextract(sizmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List sizeb_proxy = modelextract(sibmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List sizec_proxy = modelextract(sicmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List repst_proxy = modelextract(repmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List fec_proxy = modelextract(fecmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jsurv_proxy = modelextract(jsurmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jobs_proxy = modelextract(jobsmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jsize_proxy = modelextract(jsizmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jsizeb_proxy = modelextract(jsibmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jsizec_proxy = modelextract(jsicmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jrepst_proxy = modelextract(jrepmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  List jmatst_proxy = modelextract(jmatmodl, pmnames, mainyears, mainpatches,
    maingroups, mainindcova, mainindcovb, mainindcovc, false);
  
  // vrm_frame creation
  CharacterVector main_effects = {"intercept", "size2", "size1", "sizeb2",
    "sizeb1", "sizec2", "sizec1", "repst2", "repst1", "age", "density",
    "indcova2", "indcova1", "indcovb2", "indcovb1", "indcovc2", "indcovc1"}; // 17 total terms
  
  CharacterVector main_effect_1 = {main_effects[0], main_effects[1],
    main_effects[2], main_effects[3], main_effects[4], main_effects[5],
    main_effects[6], main_effects[7], main_effects[8], main_effects[9],
    main_effects[10], main_effects[11], main_effects[12], main_effects[13],
    main_effects[14], main_effects[15], main_effects[16], main_effects[8],
    main_effects[2], main_effects[2], main_effects[1], main_effects[1],
    main_effects[2], main_effects[9], main_effects[9], main_effects[9],
    main_effects[9], main_effects[11], main_effects[13], main_effects[15],
    main_effects[11], main_effects[13], main_effects[15], main_effects[12],
    main_effects[14], main_effects[16], main_effects[12], main_effects[14],
    main_effects[16], main_effects[11], main_effects[11], main_effects[13],
    main_effects[12], main_effects[12], main_effects[14], main_effects[11],
    main_effects[12], main_effects[11], main_effects[12], main_effects[13],
    main_effects[14], main_effects[3], main_effects[5], main_effects[2],
    main_effects[2], main_effects[4], main_effects[1], main_effects[1],
    main_effects[3], main_effects[2], main_effects[2], main_effects[4],
    main_effects[1], main_effects[1], main_effects[3], main_effects[10],
    main_effects[10], main_effects[10], main_effects[10], main_effects[10],
    main_effects[10], main_effects[10], main_effects[10], main_effects[3],
    main_effects[5], main_effects[4], main_effects[3], main_effects[4],
    main_effects[6], main_effects[5], main_effects[6], main_effects[3],
    main_effects[5], main_effects[10], main_effects[4], main_effects[6],
    main_effects[11], main_effects[11], main_effects[11], main_effects[12],
    main_effects[12], main_effects[12], main_effects[12], main_effects[11],
    main_effects[11], main_effects[12], main_effects[13], main_effects[13],
    main_effects[13], main_effects[14], main_effects[14], main_effects[14],
    main_effects[14], main_effects[13], main_effects[13], main_effects[14],
    main_effects[15], main_effects[15], main_effects[15], main_effects[16],
    main_effects[16], main_effects[16], main_effects[16], main_effects[15],
    main_effects[15], main_effects[16], main_effects[11], main_effects[13],
    main_effects[15], main_effects[12], main_effects[14], main_effects[16],
    main_effects[11], main_effects[13], main_effects[15], main_effects[12],
    main_effects[14], main_effects[16]};
  
  CharacterVector main_effect_2 = {"", "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", main_effects[7], main_effects[1], main_effects[8],
    main_effects[7], main_effects[8], main_effects[7], main_effects[2],
    main_effects[1], main_effects[8], main_effects[7], main_effects[1],
    main_effects[1], main_effects[1], main_effects[7], main_effects[7],
    main_effects[7], main_effects[2], main_effects[2], main_effects[2],
    main_effects[8], main_effects[8], main_effects[8], main_effects[13],
    main_effects[15], main_effects[15], main_effects[14], main_effects[16],
    main_effects[16], main_effects[14], main_effects[13], main_effects[16],
    main_effects[15], main_effects[16], main_effects[15], main_effects[4],
    main_effects[6], main_effects[4], main_effects[6], main_effects[6],
    main_effects[3], main_effects[5], main_effects[5], main_effects[3],
    main_effects[5], main_effects[5], main_effects[4], main_effects[6],
    main_effects[6], main_effects[1], main_effects[3], main_effects[5],
    main_effects[2], main_effects[4], main_effects[6], main_effects[7],
    main_effects[8], main_effects[7], main_effects[7], main_effects[8],
    main_effects[8], main_effects[7], main_effects[8], main_effects[8],
    main_effects[7], main_effects[9], main_effects[9], main_effects[9],
    main_effects[9], main_effects[9], main_effects[3], main_effects[5],
    main_effects[10], main_effects[4], main_effects[6], main_effects[3],
    main_effects[5], main_effects[4], main_effects[6], main_effects[10],
    main_effects[3], main_effects[5], main_effects[10], main_effects[4],
    main_effects[6], main_effects[3], main_effects[5], main_effects[4],
    main_effects[6], main_effects[10], main_effects[3], main_effects[5],
    main_effects[10], main_effects[4], main_effects[6], main_effects[3],
    main_effects[5], main_effects[4], main_effects[6], main_effects[10],
    main_effects[2], main_effects[2], main_effects[2], main_effects[1],
    main_effects[1], main_effects[1], main_effects[8], main_effects[8],
    main_effects[8], main_effects[7], main_effects[7], main_effects[7]};
  
  CharacterVector main_defined = {"y-intercept", "sizea in time t",
    "sizea in time t-1", "sizeb in time t", "sizeb in time t-1",
    "sizec in time t", "sizec in time t-1", "reproductive status in time t",
    "reproductive status in time t-1", "age in time t", "density in time t",
    "individual covariate a in time t", "individual covariate a in time t-1",
    "individual covariate b in time t", "individual covariate b in time t-1",
    "individual covariate c in time t", "individual covariate c in time t-1"};
  
  CharacterVector main_defined_1 = {main_defined[0], main_defined[1],
    main_defined[2], main_defined[3], main_defined[4], main_defined[5],
    main_defined[6], main_defined[7], main_defined[8], main_defined[9],
    main_defined[10], main_defined[11], main_defined[12], main_defined[13],
    main_defined[14], main_defined[15], main_defined[16], main_defined[8],
    main_defined[2], main_defined[2], main_defined[1], main_defined[1],
    main_defined[2], main_defined[9], main_defined[9], main_defined[9],
    main_defined[9], main_defined[11], main_defined[13], main_defined[15],
    main_defined[11], main_defined[13], main_defined[15], main_defined[12],
    main_defined[14], main_defined[16], main_defined[12], main_defined[14],
    main_defined[16], main_defined[11], main_defined[11], main_defined[13],
    main_defined[12], main_defined[12], main_defined[14], main_defined[11],
    main_defined[12], main_defined[11], main_defined[12], main_defined[13],
    main_defined[14], main_defined[3], main_defined[5], main_defined[2],
    main_defined[2], main_defined[4], main_defined[1], main_defined[1],
    main_defined[3], main_defined[2], main_defined[2], main_defined[4],
    main_defined[1], main_defined[1], main_defined[3], main_defined[10],
    main_defined[10], main_defined[10], main_defined[10], main_defined[10],
    main_defined[10], main_defined[10], main_defined[10], main_defined[3],
    main_defined[5], main_defined[4], main_defined[3], main_defined[4],
    main_defined[6], main_defined[5], main_defined[6], main_defined[3],
    main_defined[5], main_defined[10], main_defined[4], main_defined[6],
    main_defined[11], main_defined[11], main_defined[11], main_defined[12],
    main_defined[12], main_defined[12], main_defined[12], main_defined[11],
    main_defined[11], main_defined[12], main_defined[13], main_defined[13],
    main_defined[13], main_defined[14], main_defined[14], main_defined[14],
    main_defined[14], main_defined[13], main_defined[13], main_defined[14],
    main_defined[15], main_defined[15], main_defined[15], main_defined[16],
    main_defined[16], main_defined[16], main_defined[16], main_defined[15],
    main_defined[15], main_defined[16], main_defined[11], main_defined[13],
    main_defined[15], main_defined[12], main_defined[14], main_defined[16],
    main_defined[11], main_defined[13], main_defined[15], main_defined[12],
    main_defined[14], main_defined[16]};
  
  CharacterVector main_defined_2 = {"", "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", main_defined[7], main_defined[1], main_defined[8],
    main_defined[7], main_defined[8], main_defined[7], main_defined[2],
    main_defined[1], main_defined[8], main_defined[7], main_defined[1],
    main_defined[1], main_defined[1], main_defined[7], main_defined[7],
    main_defined[7], main_defined[2], main_defined[2], main_defined[2],
    main_defined[8], main_defined[8], main_defined[8], main_defined[13],
    main_defined[15], main_defined[15], main_defined[14], main_defined[16],
    main_defined[16], main_defined[14], main_defined[13], main_defined[16],
    main_defined[15], main_defined[16], main_defined[15], main_defined[4],
    main_defined[6], main_defined[4], main_defined[6], main_defined[6],
    main_defined[3], main_defined[5], main_defined[5], main_defined[3],
    main_defined[5], main_defined[5], main_defined[4], main_defined[6],
    main_defined[6], main_defined[1], main_defined[3], main_defined[5],
    main_defined[2], main_defined[4], main_defined[6], main_defined[7],
    main_defined[8], main_defined[7], main_defined[7], main_defined[8],
    main_defined[8], main_defined[7], main_defined[8], main_defined[8],
    main_defined[7], main_defined[9], main_defined[9], main_defined[9],
    main_defined[9], main_defined[9], main_defined[3], main_defined[5],
    main_defined[10], main_defined[4], main_defined[6], main_defined[3],
    main_defined[5], main_defined[4], main_defined[6], main_defined[10],
    main_defined[3], main_defined[5], main_defined[10], main_defined[4],
    main_defined[6], main_defined[3], main_defined[5], main_defined[4],
    main_defined[6], main_defined[10], main_defined[3], main_defined[5],
    main_defined[10], main_defined[4], main_defined[6], main_defined[3],
    main_defined[5], main_defined[4], main_defined[6], main_defined[10],
    main_defined[2], main_defined[2], main_defined[2], main_defined[1],
    main_defined[1], main_defined[1], main_defined[8], main_defined[8],
    main_defined[8], main_defined[7], main_defined[7], main_defined[7]};
  
  NumericVector surv_slopes = vrmf_inator(as<NumericVector>(surv_proxy["coefficients"]), false);
  NumericVector obs_slopes = vrmf_inator(as<NumericVector>(obs_proxy["coefficients"]), false);
  NumericVector sizea_slopes = vrmf_inator(as<NumericVector>(size_proxy["coefficients"]), false);
  NumericVector sizeb_slopes = vrmf_inator(as<NumericVector>(sizeb_proxy["coefficients"]), false);
  NumericVector sizec_slopes = vrmf_inator(as<NumericVector>(sizec_proxy["coefficients"]), false);
  NumericVector repst_slopes = vrmf_inator(as<NumericVector>(repst_proxy["coefficients"]), false);
  NumericVector fec_slopes = vrmf_inator(as<NumericVector>(fec_proxy["coefficients"]), false);
  NumericVector jsurv_slopes = vrmf_inator(as<NumericVector>(jsurv_proxy["coefficients"]), false);
  NumericVector jobs_slopes = vrmf_inator(as<NumericVector>(jobs_proxy["coefficients"]), false);
  NumericVector jsizea_slopes = vrmf_inator(as<NumericVector>(jsize_proxy["coefficients"]), false);
  NumericVector jsizeb_slopes = vrmf_inator(as<NumericVector>(jsizeb_proxy["coefficients"]), false);
  NumericVector jsizec_slopes = vrmf_inator(as<NumericVector>(jsizec_proxy["coefficients"]), false);
  NumericVector jrepst_slopes = vrmf_inator(as<NumericVector>(jrepst_proxy["coefficients"]), false);
  NumericVector jmatst_slopes = vrmf_inator(as<NumericVector>(jmatst_proxy["coefficients"]), false);
  NumericVector sizea_zi_slopes = vrmf_inator(as<NumericVector>(size_proxy["coefficients"]), true);
  NumericVector sizeb_zi_slopes = vrmf_inator(as<NumericVector>(sizeb_proxy["coefficients"]), true);
  NumericVector sizec_zi_slopes = vrmf_inator(as<NumericVector>(sizec_proxy["coefficients"]), true);
  NumericVector fec_zi_slopes = vrmf_inator(as<NumericVector>(fec_proxy["coefficients"]), true);
  NumericVector jsizea_zi_slopes = vrmf_inator(as<NumericVector>(jsize_proxy["coefficients"]), true);
  NumericVector jsizeb_zi_slopes = vrmf_inator(as<NumericVector>(jsizeb_proxy["coefficients"]), true);
  NumericVector jsizec_zi_slopes = vrmf_inator(as<NumericVector>(jsizec_proxy["coefficients"]), true);
  
  List vrm_frame (25);
  vrm_frame(0) = main_effect_1;
  vrm_frame(1) = main_defined_1;
  vrm_frame(2) = main_effect_2;
  vrm_frame(3) = main_defined_2;
  vrm_frame(4) = surv_slopes;
  vrm_frame(5) = obs_slopes;
  vrm_frame(6) = sizea_slopes;
  vrm_frame(7) = sizeb_slopes;
  vrm_frame(8) = sizec_slopes;
  vrm_frame(9) = repst_slopes;
  vrm_frame(10) = fec_slopes;
  vrm_frame(11) = jsurv_slopes;
  vrm_frame(12) = jobs_slopes;
  vrm_frame(13) = jsizea_slopes;
  vrm_frame(14) = jsizeb_slopes;
  vrm_frame(15) = jsizec_slopes;
  vrm_frame(16) = jrepst_slopes;
  vrm_frame(17) = jmatst_slopes;
  vrm_frame(18) = sizea_zi_slopes;
  vrm_frame(19) = sizeb_zi_slopes;
  vrm_frame(20) = sizec_zi_slopes;
  vrm_frame(21) = fec_zi_slopes;
  vrm_frame(22) = jsizea_zi_slopes;
  vrm_frame(23) = jsizeb_zi_slopes;
  vrm_frame(24) = jsizec_zi_slopes;
  
  CharacterVector vrm_frame_names = {"main_effect_1", "main_1_defined",
    "main_effect_2", "main_2_defined", "surv", "obs", "sizea", "sizeb", "sizec",
    "repst", "fec", "jsurv", "jobs", "jsizea", "jsizeb", "jsizec", "jrepst",
    "jmatst", "sizea_zi", "sizeb_zi", "sizec_zi", "fec_zi", "jsizea_zi",
    "jsizeb_zi", "jsizec_zi"};
  
  vrm_frame.attr("names") = vrm_frame_names;
  vrm_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, 128);
  StringVector needed_classes {"data.frame", "vrm_input"};
  vrm_frame.attr("class") = needed_classes;
  
  // year_frame creation
  List year_frame (22);
  year_frame(0) = mainyears;
  year_frame(1) = as<NumericVector>(surv_proxy["years"]);
  year_frame(2) = as<NumericVector>(obs_proxy["years"]);
  year_frame(3) = as<NumericVector>(size_proxy["years"]);
  year_frame(4) = as<NumericVector>(sizeb_proxy["years"]);
  year_frame(5) = as<NumericVector>(sizec_proxy["years"]);
  year_frame(6) = as<NumericVector>(repst_proxy["years"]);
  year_frame(7) = as<NumericVector>(fec_proxy["years"]);
  year_frame(8) = as<NumericVector>(jsurv_proxy["years"]);
  year_frame(9) = as<NumericVector>(jobs_proxy["years"]);
  year_frame(10) = as<NumericVector>(jsize_proxy["years"]);
  year_frame(11) = as<NumericVector>(jsizeb_proxy["years"]);
  year_frame(12) = as<NumericVector>(jsizec_proxy["years"]);
  year_frame(13) = as<NumericVector>(jrepst_proxy["years"]);
  year_frame(14) = as<NumericVector>(jmatst_proxy["years"]);
  year_frame(15) = as<NumericVector>(size_proxy["zeroyear"]);
  year_frame(16) = as<NumericVector>(sizeb_proxy["zeroyear"]);
  year_frame(17) = as<NumericVector>(sizec_proxy["zeroyear"]);
  year_frame(18) = as<NumericVector>(fec_proxy["zeroyear"]);
  year_frame(19) = as<NumericVector>(jsize_proxy["zeroyear"]);
  year_frame(20) = as<NumericVector>(jsizeb_proxy["zeroyear"]);
  year_frame(21) = as<NumericVector>(jsizec_proxy["zeroyear"]);
  
  CharacterVector yf_names = {"years", "surv", "obs", "sizea", "sizeb", "sizec",
    "repst", "fec", "jsurv", "jobs", "jsizea", "jsizeb", "jsizec", "jrepst",
    "jmatst", "sizea_zi", "sizeb_zi", "sizec_zi", "fec_zi", "jsizea_zi",
    "jsizeb_zi", "jsizec_zi"};
  year_frame.attr("names") = yf_names;
  year_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainyears.length()));
  StringVector yf_class {"data.frame"};
  year_frame.attr("class") = yf_class;
  
  // patch_frame creation
  List patch_frame (22);
  patch_frame(0) = mainpatches;
  patch_frame(1) = as<NumericVector>(surv_proxy["patches"]);
  patch_frame(2) = as<NumericVector>(obs_proxy["patches"]);
  patch_frame(3) = as<NumericVector>(size_proxy["patches"]);
  patch_frame(4) = as<NumericVector>(sizeb_proxy["patches"]);
  patch_frame(5) = as<NumericVector>(sizec_proxy["patches"]);
  patch_frame(6) = as<NumericVector>(repst_proxy["patches"]);
  patch_frame(7) = as<NumericVector>(fec_proxy["patches"]);
  patch_frame(8) = as<NumericVector>(jsurv_proxy["patches"]);
  patch_frame(9) = as<NumericVector>(jobs_proxy["patches"]);
  patch_frame(10) = as<NumericVector>(jsize_proxy["patches"]);
  patch_frame(11) = as<NumericVector>(jsizeb_proxy["patches"]);
  patch_frame(12) = as<NumericVector>(jsizec_proxy["patches"]);
  patch_frame(13) = as<NumericVector>(jrepst_proxy["patches"]);
  patch_frame(14) = as<NumericVector>(jmatst_proxy["patches"]);
  patch_frame(15) = as<NumericVector>(size_proxy["zeropatch"]);
  patch_frame(16) = as<NumericVector>(sizeb_proxy["zeropatch"]);
  patch_frame(17) = as<NumericVector>(sizec_proxy["zeropatch"]);
  patch_frame(18) = as<NumericVector>(fec_proxy["zeropatch"]);
  patch_frame(19) = as<NumericVector>(jsize_proxy["zeropatch"]);
  patch_frame(20) = as<NumericVector>(jsizeb_proxy["zeropatch"]);
  patch_frame(21) = as<NumericVector>(jsizec_proxy["zeropatch"]);
  
  CharacterVector pf_names = clone(yf_names);
  pf_names(0) = "patches";
  patch_frame.attr("names") = pf_names;
  patch_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainpatches.length()));
  patch_frame.attr("class") = yf_class;
  
  // group2_frame creation
  List group2_frame (22);
  group2_frame(0) = maingroups;
  group2_frame(1) = as<NumericVector>(surv_proxy["groups2"]);
  group2_frame(2) = as<NumericVector>(obs_proxy["groups2"]);
  group2_frame(3) = as<NumericVector>(size_proxy["groups2"]);
  group2_frame(4) = as<NumericVector>(sizeb_proxy["groups2"]);
  group2_frame(5) = as<NumericVector>(sizec_proxy["groups2"]);
  group2_frame(6) = as<NumericVector>(repst_proxy["groups2"]);
  group2_frame(7) = as<NumericVector>(fec_proxy["groups2"]);
  group2_frame(8) = as<NumericVector>(jsurv_proxy["groups2"]);
  group2_frame(9) = as<NumericVector>(jobs_proxy["groups2"]);
  group2_frame(10) = as<NumericVector>(jsize_proxy["groups2"]);
  group2_frame(11) = as<NumericVector>(jsizeb_proxy["groups2"]);
  group2_frame(12) = as<NumericVector>(jsizec_proxy["groups2"]);
  group2_frame(13) = as<NumericVector>(jrepst_proxy["groups2"]);
  group2_frame(14) = as<NumericVector>(jmatst_proxy["groups2"]);
  group2_frame(15) = as<NumericVector>(size_proxy["zerogroups2"]);
  group2_frame(16) = as<NumericVector>(sizeb_proxy["zerogroups2"]);
  group2_frame(17) = as<NumericVector>(sizec_proxy["zerogroups2"]);
  group2_frame(18) = as<NumericVector>(fec_proxy["zerogroups2"]);
  group2_frame(19) = as<NumericVector>(jsize_proxy["zerogroups2"]);
  group2_frame(20) = as<NumericVector>(jsizeb_proxy["zerogroups2"]);
  group2_frame(21) = as<NumericVector>(jsizec_proxy["zerogroups2"]);
  
  CharacterVector g2f_names = clone(yf_names);
  g2f_names(0) = "groups";
  group2_frame.attr("names") = g2f_names;
  group2_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(maingroups.length()));
  group2_frame.attr("class") = yf_class;
  
  // group1_frame creation
  List group1_frame (22);
  group1_frame(0) = maingroups;
  group1_frame(1) = as<NumericVector>(surv_proxy["groups1"]);
  group1_frame(2) = as<NumericVector>(obs_proxy["groups1"]);
  group1_frame(3) = as<NumericVector>(size_proxy["groups1"]);
  group1_frame(4) = as<NumericVector>(sizeb_proxy["groups1"]);
  group1_frame(5) = as<NumericVector>(sizec_proxy["groups1"]);
  group1_frame(6) = as<NumericVector>(repst_proxy["groups1"]);
  group1_frame(7) = as<NumericVector>(fec_proxy["groups1"]);
  group1_frame(8) = as<NumericVector>(jsurv_proxy["groups1"]);
  group1_frame(9) = as<NumericVector>(jobs_proxy["groups1"]);
  group1_frame(10) = as<NumericVector>(jsize_proxy["groups1"]);
  group1_frame(11) = as<NumericVector>(jsizeb_proxy["groups1"]);
  group1_frame(12) = as<NumericVector>(jsizec_proxy["groups1"]);
  group1_frame(13) = as<NumericVector>(jrepst_proxy["groups1"]);
  group1_frame(14) = as<NumericVector>(jmatst_proxy["groups1"]);
  group1_frame(15) = as<NumericVector>(size_proxy["zerogroups1"]);
  group1_frame(16) = as<NumericVector>(sizeb_proxy["zerogroups1"]);
  group1_frame(17) = as<NumericVector>(sizec_proxy["zerogroups1"]);
  group1_frame(18) = as<NumericVector>(fec_proxy["zerogroups1"]);
  group1_frame(19) = as<NumericVector>(jsize_proxy["zerogroups1"]);
  group1_frame(20) = as<NumericVector>(jsizeb_proxy["zerogroups1"]);
  group1_frame(21) = as<NumericVector>(jsizec_proxy["zerogroups1"]);
  
  group1_frame.attr("names") = g2f_names;
  group1_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(maingroups.length()));
  group1_frame.attr("class") = yf_class;
  
  // dist_frame creation
  CharacterVector distf_response = {"surv", "obs", "sizea", "sizeb", "sizec",
    "repst", "fec", "jsurv", "jobs", "jsizea", "jsizeb", "jsizec", "jrepst",
    "jmatst"};
  CharacterVector distf_dist (14);
  
  distf_dist(0) = String(surv_proxy["family"]);
  distf_dist(1) = String(obs_proxy["family"]);
  distf_dist(2) = String(size_proxy["family"]);
  distf_dist(3) = String(sizeb_proxy["family"]);
  distf_dist(4) = String(sizec_proxy["family"]);
  distf_dist(5) = String(repst_proxy["family"]);
  distf_dist(6) = String(fec_proxy["family"]);
  distf_dist(7) = String(jsurv_proxy["family"]);
  distf_dist(8) = String(jobs_proxy["family"]);
  distf_dist(9) = String(jsize_proxy["family"]);
  distf_dist(10) = String(jsizeb_proxy["family"]);
  distf_dist(11) = String(jsizec_proxy["family"]);
  distf_dist(12) = String(jrepst_proxy["family"]);
  distf_dist(13) = String(jmatst_proxy["family"]);
  
  for (int i = 0; i < 14; i++) {
    if (distf_dist(i) == "binomial") distf_dist(i) = "binom";
  }
  
  DataFrame dist_frame = DataFrame::create(_["response"] = distf_response,
    _["dist"] = distf_dist);
  
  // indcova2_frame creation
  List indcova2_frame (22);
  indcova2_frame(0) = maingroups;
  indcova2_frame(1) = as<NumericVector>(surv_proxy["indcova2s"]);
  indcova2_frame(2) = as<NumericVector>(obs_proxy["indcova2s"]);
  indcova2_frame(3) = as<NumericVector>(size_proxy["indcova2s"]);
  indcova2_frame(4) = as<NumericVector>(sizeb_proxy["indcova2s"]);
  indcova2_frame(5) = as<NumericVector>(sizec_proxy["indcova2s"]);
  indcova2_frame(6) = as<NumericVector>(repst_proxy["indcova2s"]);
  indcova2_frame(7) = as<NumericVector>(fec_proxy["indcova2s"]);
  indcova2_frame(8) = as<NumericVector>(jsurv_proxy["indcova2s"]);
  indcova2_frame(9) = as<NumericVector>(jobs_proxy["indcova2s"]);
  indcova2_frame(10) = as<NumericVector>(jsize_proxy["indcova2s"]);
  indcova2_frame(11) = as<NumericVector>(jsizeb_proxy["indcova2s"]);
  indcova2_frame(12) = as<NumericVector>(jsizec_proxy["indcova2s"]);
  indcova2_frame(13) = as<NumericVector>(jrepst_proxy["indcova2s"]);
  indcova2_frame(14) = as<NumericVector>(jmatst_proxy["indcova2s"]);
  indcova2_frame(15) = as<NumericVector>(size_proxy["zeroindcova2s"]);
  indcova2_frame(16) = as<NumericVector>(sizeb_proxy["zeroindcova2s"]);
  indcova2_frame(17) = as<NumericVector>(sizec_proxy["zeroindcova2s"]);
  indcova2_frame(18) = as<NumericVector>(fec_proxy["zeroindcova2s"]);
  indcova2_frame(19) = as<NumericVector>(jsize_proxy["zeroindcova2s"]);
  indcova2_frame(20) = as<NumericVector>(jsizeb_proxy["zeroindcova2s"]);
  indcova2_frame(21) = as<NumericVector>(jsizec_proxy["zeroindcova2s"]);
  
  CharacterVector ica2f_names = clone(yf_names);
  ica2f_names(0) = "indcova";
  indcova2_frame.attr("names") = ica2f_names;
  indcova2_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainindcova.length()));
  indcova2_frame.attr("class") = yf_class;
  
  // indcova1_frame creation
  List indcova1_frame (22);
  indcova1_frame(0) = maingroups;
  indcova1_frame(1) = as<NumericVector>(surv_proxy["indcova1s"]);
  indcova1_frame(2) = as<NumericVector>(obs_proxy["indcova1s"]);
  indcova1_frame(3) = as<NumericVector>(size_proxy["indcova1s"]);
  indcova1_frame(4) = as<NumericVector>(sizeb_proxy["indcova1s"]);
  indcova1_frame(5) = as<NumericVector>(sizec_proxy["indcova1s"]);
  indcova1_frame(6) = as<NumericVector>(repst_proxy["indcova1s"]);
  indcova1_frame(7) = as<NumericVector>(fec_proxy["indcova1s"]);
  indcova1_frame(8) = as<NumericVector>(jsurv_proxy["indcova1s"]);
  indcova1_frame(9) = as<NumericVector>(jobs_proxy["indcova1s"]);
  indcova1_frame(10) = as<NumericVector>(jsize_proxy["indcova1s"]);
  indcova1_frame(11) = as<NumericVector>(jsizeb_proxy["indcova1s"]);
  indcova1_frame(12) = as<NumericVector>(jsizec_proxy["indcova1s"]);
  indcova1_frame(13) = as<NumericVector>(jrepst_proxy["indcova1s"]);
  indcova1_frame(14) = as<NumericVector>(jmatst_proxy["indcova1s"]);
  indcova1_frame(15) = as<NumericVector>(size_proxy["zeroindcova1s"]);
  indcova1_frame(16) = as<NumericVector>(sizeb_proxy["zeroindcova1s"]);
  indcova1_frame(17) = as<NumericVector>(sizec_proxy["zeroindcova1s"]);
  indcova1_frame(18) = as<NumericVector>(fec_proxy["zeroindcova1s"]);
  indcova1_frame(19) = as<NumericVector>(jsize_proxy["zeroindcova1s"]);
  indcova1_frame(20) = as<NumericVector>(jsizeb_proxy["zeroindcova1s"]);
  indcova1_frame(21) = as<NumericVector>(jsizec_proxy["zeroindcova1s"]);
  
  indcova1_frame.attr("names") = ica2f_names;
  indcova1_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainindcova.length()));
  indcova1_frame.attr("class") = yf_class;
  
  // indcovb2_frame creation
  List indcovb2_frame (22);
  indcovb2_frame(0) = maingroups;
  indcovb2_frame(1) = as<NumericVector>(surv_proxy["indcovb2s"]);
  indcovb2_frame(2) = as<NumericVector>(obs_proxy["indcovb2s"]);
  indcovb2_frame(3) = as<NumericVector>(size_proxy["indcovb2s"]);
  indcovb2_frame(4) = as<NumericVector>(sizeb_proxy["indcovb2s"]);
  indcovb2_frame(5) = as<NumericVector>(sizec_proxy["indcovb2s"]);
  indcovb2_frame(6) = as<NumericVector>(repst_proxy["indcovb2s"]);
  indcovb2_frame(7) = as<NumericVector>(fec_proxy["indcovb2s"]);
  indcovb2_frame(8) = as<NumericVector>(jsurv_proxy["indcovb2s"]);
  indcovb2_frame(9) = as<NumericVector>(jobs_proxy["indcovb2s"]);
  indcovb2_frame(10) = as<NumericVector>(jsize_proxy["indcovb2s"]);
  indcovb2_frame(11) = as<NumericVector>(jsizeb_proxy["indcovb2s"]);
  indcovb2_frame(12) = as<NumericVector>(jsizec_proxy["indcovb2s"]);
  indcovb2_frame(13) = as<NumericVector>(jrepst_proxy["indcovb2s"]);
  indcovb2_frame(14) = as<NumericVector>(jmatst_proxy["indcovb2s"]);
  indcovb2_frame(15) = as<NumericVector>(size_proxy["zeroindcovb2s"]);
  indcovb2_frame(16) = as<NumericVector>(sizeb_proxy["zeroindcovb2s"]);
  indcovb2_frame(17) = as<NumericVector>(sizec_proxy["zeroindcovb2s"]);
  indcovb2_frame(18) = as<NumericVector>(fec_proxy["zeroindcovb2s"]);
  indcovb2_frame(19) = as<NumericVector>(jsize_proxy["zeroindcovb2s"]);
  indcovb2_frame(20) = as<NumericVector>(jsizeb_proxy["zeroindcovb2s"]);
  indcovb2_frame(21) = as<NumericVector>(jsizec_proxy["zeroindcovb2s"]);
  
  CharacterVector icb2f_names = clone(yf_names);
  icb2f_names(0) = "indcovb";
  indcovb2_frame.attr("names") = icb2f_names;
  indcovb2_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainindcovb.length()));
  indcovb2_frame.attr("class") = yf_class;
  
  // indcovb1_frame creation
  List indcovb1_frame (22);
  indcovb1_frame(0) = maingroups;
  indcovb1_frame(1) = as<NumericVector>(surv_proxy["indcovb1s"]);
  indcovb1_frame(2) = as<NumericVector>(obs_proxy["indcovb1s"]);
  indcovb1_frame(3) = as<NumericVector>(size_proxy["indcovb1s"]);
  indcovb1_frame(4) = as<NumericVector>(sizeb_proxy["indcovb1s"]);
  indcovb1_frame(5) = as<NumericVector>(sizec_proxy["indcovb1s"]);
  indcovb1_frame(6) = as<NumericVector>(repst_proxy["indcovb1s"]);
  indcovb1_frame(7) = as<NumericVector>(fec_proxy["indcovb1s"]);
  indcovb1_frame(8) = as<NumericVector>(jsurv_proxy["indcovb1s"]);
  indcovb1_frame(9) = as<NumericVector>(jobs_proxy["indcovb1s"]);
  indcovb1_frame(10) = as<NumericVector>(jsize_proxy["indcovb1s"]);
  indcovb1_frame(11) = as<NumericVector>(jsizeb_proxy["indcovb1s"]);
  indcovb1_frame(12) = as<NumericVector>(jsizec_proxy["indcovb1s"]);
  indcovb1_frame(13) = as<NumericVector>(jrepst_proxy["indcovb1s"]);
  indcovb1_frame(14) = as<NumericVector>(jmatst_proxy["indcovb1s"]);
  indcovb1_frame(15) = as<NumericVector>(size_proxy["zeroindcovb1s"]);
  indcovb1_frame(16) = as<NumericVector>(sizeb_proxy["zeroindcovb1s"]);
  indcovb1_frame(17) = as<NumericVector>(sizec_proxy["zeroindcovb1s"]);
  indcovb1_frame(18) = as<NumericVector>(fec_proxy["zeroindcovb1s"]);
  indcovb1_frame(19) = as<NumericVector>(jsize_proxy["zeroindcovb1s"]);
  indcovb1_frame(20) = as<NumericVector>(jsizeb_proxy["zeroindcovb1s"]);
  indcovb1_frame(21) = as<NumericVector>(jsizec_proxy["zeroindcovb1s"]);
  
  indcovb1_frame.attr("names") = icb2f_names;
  indcovb1_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainindcovb.length()));
  indcovb1_frame.attr("class") = yf_class;
  
  // indcovc2_frame creation
  List indcovc2_frame (22);
  indcovc2_frame(0) = maingroups;
  indcovc2_frame(1) = as<NumericVector>(surv_proxy["indcovc2s"]);
  indcovc2_frame(2) = as<NumericVector>(obs_proxy["indcovc2s"]);
  indcovc2_frame(3) = as<NumericVector>(size_proxy["indcovc2s"]);
  indcovc2_frame(4) = as<NumericVector>(sizeb_proxy["indcovc2s"]);
  indcovc2_frame(5) = as<NumericVector>(sizec_proxy["indcovc2s"]);
  indcovc2_frame(6) = as<NumericVector>(repst_proxy["indcovc2s"]);
  indcovc2_frame(7) = as<NumericVector>(fec_proxy["indcovc2s"]);
  indcovc2_frame(8) = as<NumericVector>(jsurv_proxy["indcovc2s"]);
  indcovc2_frame(9) = as<NumericVector>(jobs_proxy["indcovc2s"]);
  indcovc2_frame(10) = as<NumericVector>(jsize_proxy["indcovc2s"]);
  indcovc2_frame(11) = as<NumericVector>(jsizeb_proxy["indcovc2s"]);
  indcovc2_frame(12) = as<NumericVector>(jsizec_proxy["indcovc2s"]);
  indcovc2_frame(13) = as<NumericVector>(jrepst_proxy["indcovc2s"]);
  indcovc2_frame(14) = as<NumericVector>(jmatst_proxy["indcovc2s"]);
  indcovc2_frame(15) = as<NumericVector>(size_proxy["zeroindcovc2s"]);
  indcovc2_frame(16) = as<NumericVector>(sizeb_proxy["zeroindcovc2s"]);
  indcovc2_frame(17) = as<NumericVector>(sizec_proxy["zeroindcovc2s"]);
  indcovc2_frame(18) = as<NumericVector>(fec_proxy["zeroindcovc2s"]);
  indcovc2_frame(19) = as<NumericVector>(jsize_proxy["zeroindcovc2s"]);
  indcovc2_frame(20) = as<NumericVector>(jsizeb_proxy["zeroindcovc2s"]);
  indcovc2_frame(21) = as<NumericVector>(jsizec_proxy["zeroindcovc2s"]);
  
  CharacterVector icc2f_names = clone(yf_names);
  icc2f_names(0) = "indcovc";
  indcovc2_frame.attr("names") = icc2f_names;
  indcovc2_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainindcovc.length()));
  indcovc2_frame.attr("class") = yf_class;
  
  // indcovc1_frame creation
  List indcovc1_frame (22);
  indcovc1_frame(0) = maingroups;
  indcovc1_frame(1) = as<NumericVector>(surv_proxy["indcovc1s"]);
  indcovc1_frame(2) = as<NumericVector>(obs_proxy["indcovc1s"]);
  indcovc1_frame(3) = as<NumericVector>(size_proxy["indcovc1s"]);
  indcovc1_frame(4) = as<NumericVector>(sizeb_proxy["indcovc1s"]);
  indcovc1_frame(5) = as<NumericVector>(sizec_proxy["indcovc1s"]);
  indcovc1_frame(6) = as<NumericVector>(repst_proxy["indcovc1s"]);
  indcovc1_frame(7) = as<NumericVector>(fec_proxy["indcovc1s"]);
  indcovc1_frame(8) = as<NumericVector>(jsurv_proxy["indcovc1s"]);
  indcovc1_frame(9) = as<NumericVector>(jobs_proxy["indcovc1s"]);
  indcovc1_frame(10) = as<NumericVector>(jsize_proxy["indcovc1s"]);
  indcovc1_frame(11) = as<NumericVector>(jsizeb_proxy["indcovc1s"]);
  indcovc1_frame(12) = as<NumericVector>(jsizec_proxy["indcovc1s"]);
  indcovc1_frame(13) = as<NumericVector>(jrepst_proxy["indcovc1s"]);
  indcovc1_frame(14) = as<NumericVector>(jmatst_proxy["indcovc1s"]);
  indcovc1_frame(15) = as<NumericVector>(size_proxy["zeroindcovc1s"]);
  indcovc1_frame(16) = as<NumericVector>(sizeb_proxy["zeroindcovc1s"]);
  indcovc1_frame(17) = as<NumericVector>(sizec_proxy["zeroindcovc1s"]);
  indcovc1_frame(18) = as<NumericVector>(fec_proxy["zeroindcovc1s"]);
  indcovc1_frame(19) = as<NumericVector>(jsize_proxy["zeroindcovc1s"]);
  indcovc1_frame(20) = as<NumericVector>(jsizeb_proxy["zeroindcovc1s"]);
  indcovc1_frame(21) = as<NumericVector>(jsizec_proxy["zeroindcovc1s"]);
  
  indcovc1_frame.attr("names") = icc2f_names;
  indcovc1_frame.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, static_cast<int>(mainindcovc.length()));
  indcovc1_frame.attr("class") = yf_class;
  
  // st_frame creation
  NumericVector st_frame (14);
  st_frame(0) = static_cast<double>(as<NumericVector>(surv_proxy["sigma"])(0));
  st_frame(1) = static_cast<double>(as<NumericVector>(obs_proxy["sigma"])(0));
  
  if (distf_dist(2) == "negbin") {
    st_frame(2) = static_cast<double>(as<NumericVector>(size_proxy["theta"])(0));
  } else {
    st_frame(2) = static_cast<double>(as<NumericVector>(size_proxy["sigma"])(0));
  }
  if (distf_dist(3) == "negbin") {
    st_frame(3) = static_cast<double>(as<NumericVector>(sizeb_proxy["theta"])(0));
  } else {
    st_frame(3) = static_cast<double>(as<NumericVector>(sizeb_proxy["sigma"])(0));
  }
  if (distf_dist(4) == "negbin") {
    st_frame(4) = static_cast<double>(as<NumericVector>(sizec_proxy["theta"])(0));
  } else {
    st_frame(4) = static_cast<double>(as<NumericVector>(sizec_proxy["sigma"])(0));
  }
  
  st_frame(5) = static_cast<double>(as<NumericVector>(repst_proxy["sigma"])(0));
  
  if (distf_dist(6) == "negbin") {
    st_frame(6) = static_cast<double>(as<NumericVector>(fec_proxy["theta"])(0));
  } else {
    st_frame(6) = static_cast<double>(as<NumericVector>(fec_proxy["sigma"])(0));
  }
  
  st_frame(7) = static_cast<double>(as<NumericVector>(jsurv_proxy["sigma"])(0));
  st_frame(8) = static_cast<double>(as<NumericVector>(jobs_proxy["sigma"])(0));
  
  if (distf_dist(9) == "negbin") {
    st_frame(9) = static_cast<double>(as<NumericVector>(jsize_proxy["theta"])(0));
  } else {
    st_frame(9) = static_cast<double>(as<NumericVector>(jsize_proxy["sigma"])(0));
  }
  if (distf_dist(10) == "negbin") {
    st_frame(10) = static_cast<double>(as<NumericVector>(jsizeb_proxy["theta"])(0));
  } else {
    st_frame(10) = static_cast<double>(as<NumericVector>(jsizeb_proxy["sigma"])(0));
  }
  if (distf_dist(11) == "negbin") {
    st_frame(11) = static_cast<double>(as<NumericVector>(jsizec_proxy["theta"])(0));
  } else {
    st_frame(11) = static_cast<double>(as<NumericVector>(jsizec_proxy["sigma"])(0));
  }
  
  st_frame(12) = static_cast<double>(as<NumericVector>(jrepst_proxy["sigma"])(0));
  st_frame(13) = static_cast<double>(as<NumericVector>(jmatst_proxy["sigma"])(0));
  
  st_frame.attr("names") = distf_response;
  
  // Main output
  List output = List::create(_["vrm_frame"] = vrm_frame,
    _["year_frame"] = year_frame, _["patch_frame"] = patch_frame,
    _["group2_frame"] = group2_frame, _["group1_frame"] = group1_frame,
    _["dist_frame"] = dist_frame, _["indcova2_frame"] = indcova2_frame,
    _["indcova1_frame"] = indcova1_frame, _["indcovb2_frame"] = indcovb2_frame,
    _["indcovb1_frame"] = indcovb1_frame, _["indcovc2_frame"] = indcovc2_frame,
    _["indcovc1_frame"] = indcovc1_frame, _["st_frame"] = st_frame);
  output.attr("class") = "vrm_input";
  
  return output;
}

