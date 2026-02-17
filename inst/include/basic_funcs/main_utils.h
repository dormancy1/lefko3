#ifndef LEFKOUTILS_main_utils_H
#define LEFKOUTILS_main_utils_H

#include <RcppArmadillo.h>
#define BOOST_DISABLE_ASSERTS

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using namespace arma;

// Index of functions
// 1. NumericVector concat_dbl  Append NumericVector to the End of Another NumericVector
// 2. IntegerVector concat_int  Append IntegerVector to the End of Another IntegerVector
// 3. std::string stringremove  Remove One String From Another
// 4. StringVector concat_str  Append StringVector to the End of Another StringVector
// 
// 5. bool stringcompare_hard  Compares Two Strings Literally
// 6. List stringcompare_soft  Compares Two Strings, Assessing Inclusion
// 7. bool stringcompare_simple  Compares Two Strings, Assessing Inclusion
// 8. bool stringcompare_x  Compares Three Strings for Interaction Notation
// 9. CharacterVector stringsort  Sort String Elements
// 
// 10. IntegerVector int_sort  Sort Integer Elements
// 11. IntegerMatrix refsort_num  Index a Numeric Vector According to a Reference Vector
// 12. IntegerVector refsort_str  Index a String Vector According to a Reference Vector
// 13. IntegerMatrix refsort_str_m  Index a String Matrix According to a Reference Vector
// 
// 14. arma::sp_mat spmat_log  Calculate Logarithms of Non-Zero Elements of Sparse Matrix
// 15. Rcpp::IntegerVector shrink  Resize an IntegerVector
// 16. IntegerVector index_l3  Find Indices of a Matching String in a StringVector
// 
// 17. List df_subset  Subset Data Frame By A Specific Variable
// 18. bool df_compare  Check If Two Data Frames Are Equal
// 19. List df_subset_byrow  Subset Data Frame By Integer Vector of Row Numbers
// 20. List df_remove  Remove Rows With Specific Index Values From Data Frame
// 21. List df_shedrows  Shrink Data Frame According to Index Vector
// 22. bool df_duplicates  Search for Duplicate Data Frame Values
// 23. DataFrame df_rbind  Bind Two Data Frames By Row
// 
// 24. List exp_grd  Repeat First Vector for Each Element of Second Vector
// 
// 25. List numeric_extractor  Extract Key Components From Simple Numerical Model
// 26. List glm_extractor  Extract Key Components of lm/glm/negbin Objects
// 27. List vglm_extractor  Extract Key Components of vglm Objects
// 28. List zeroinfl_extractor  Extract Key Components of zeroinfl Objects
// 29. List lme4_extractor  Extract Key Components of merMod Objects
// 30. List glmmTMB_extractor  Extract Key Components of glmmTMB Objects
// 31. List S3_extractor  Extract Core Components From S3 Vital Rate Models
// 32. List S4_extractor  Extract Core Components From S4 Vital Rate Models
// 33. List vrm_extractor  Extract Core Components of vrm_input Models
// 
// 34. NumericMatrix revelations  Create Matrices of Year and Patch Terms in Models
// 35. double rimeotam  Create a Summation of Most Terms Needed in Vital Rate Calculation
// 36. arma::ivec foi_counter  Count Elements in Each Random Individual Covariate Portion of Model
// 37. NumericVector flightoficarus  Create Vector of Random Individual Covariate Terms
// 38. StringVector bootson  Create Concatenated Vector of Random Individual Covariate Term Names
// 39. NumericVector zero_flightoficarus  Create Vector of Random Individual Covariate Terms for Zero-Inflated Models
// 40. StringVector zero_bootson  Create Concatenated Vector of Random Individual Covariate Term Names from a Zero-Inflated Model
// 41. arma::imat foi_index  Create Index of Element Numbers for Random Individual Covariate Terms
// 42. NumericMatrix revelations_leslie  Create Matrices of Year and Patch Terms in Models in Leslie Models
// 43. arma::imat foi_index_leslie  Create Index of Element Numbers for Random Individual Covariate Terms in Leslie Models
// 44. List modelextract  Extract Coefficients from Linear Vital Rate Models
// 
// 45. double preouterator  Estimate Value for Vital Rate Based on Inputs
// 46. List jerzeibalowski  Estimate All Elements of Function-based Population Projection Matrix
// 47. List motherbalowski  Estimate All Elements of Function-based Leslie Population Projection Matrix
// 
// 48. DataFrame loy_inator  Converts Labels Element to LOY Data Frame
// 49. void matrix_reducer  Reduces Matrices In A Function-based lefkoMat Object
// 50. int whichbrew  Assess if MPM is ahistorical, historical, age-by-stage, or Leslie
// 51. void pop_error  Standardized Error Messages
// 
// 52. void density_prep  Format All Density-related Variables Based on Density Inputs
// 53. void equivalence_prep  Format All Equivalence Weight-related Variables Based on Input





namespace LefkoUtils {
  //' Append NumericVector to the End of Another NumericVector
  //' 
  //' This function appends one NumericVector fully to another.
  //' 
  //' @name concat_dbl
  //' 
  //' @param A Any NumericVector.
  //' @param B Any other NumericVector.
  //' 
  //' @return Returns a new NumericVector with elements of vector A followed by
  //' elements of vector B.
  //'
  //' @keywords internal
  //' @noRd
  inline NumericVector concat_dbl(const NumericVector& x,
    const NumericVector& y) {
    
    std::vector<double> xconv = as<std::vector<double> >(x);
    std::vector<double> yconv = as<std::vector<double> >(y);
    std::vector<double> z(xconv.size() + yconv.size());
    
    std::copy(xconv.begin(), xconv.end(), z.begin());
    std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
    
    Rcpp::NumericVector zconv(z.begin(), z.end());
    
    return zconv;
  }
  
  //' Append IntegerVector to the End of Another IntegerVector
  //' 
  //' Returns a new IntegerVector with elements of vector A followed by
  //' elements of vector B.
  //' 
  //' @name concat_int
  //' 
  //' @param A Any IntegerVector.
  //' @param B Any other IntegerVector.
  //' 
  //' @return Returns a new IntegerVector with elements of vector A followed by
  //' elements of vector B.
  //' 
  //' @keywords internal
  //' @noRd
  inline IntegerVector concat_int(const IntegerVector& x,
    const IntegerVector& y) {
    
    std::vector<long long int> xconv = as<std::vector<long long int> >(x);
    std::vector<long long int> yconv = as<std::vector<long long int> >(y);
    std::vector<long long int> z(x.size() + y.size());
    
    std::copy(xconv.begin(), xconv.end(), z.begin());
    std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
    
    Rcpp::IntegerVector zconv(z.begin(), z.end());
    
    return zconv ;
  }
  
  //' Remove One String From Another
  //' 
  //' This function takes the second string and removes it from the first, if
  //' the second string can be found within the first string. This is from
  //' https://www.codespeedy.com/remove-specific-substring-from-string-in-cpp/
  //' 
  //' @name stringremove
  //' 
  //' @param str1 The first string
  //' @param str2 The second string
  //' 
  //' @return A new shortened string.
  //' 
  //' @keywords internal
  //' @noRd
  inline std::string stringremove(std::string original, std::string substring) {
    std::size_t ind = original.find(substring);
    
    if(ind != std::string::npos) {
      original.erase(ind,substring.length());
      
    }
      
    return original;
  }
  
  //' Append StringVector to the End of Another StringVector
  //' 
  //' Returns a new StringVector with elements of vector A followed by
  //' elements of vector B.
  //' 
  //' @name concat_str
  //' 
  //' @param A Any StringVector.
  //' @param B Any other StringVector.
  //' 
  //' @return Returns a new StringVector with elements of vector A followed by
  //' elements of vector B.
  //' 
  //' @keywords internal
  //' @noRd
  inline StringVector concat_str(const StringVector& x, const StringVector& y) {
    
    std::vector<std::string> xconv = as<std::vector<std::string> >(x);
    std::vector<std::string> yconv = as<std::vector<std::string> >(y);
    std::vector<std::string> z(x.size() + y.size());
    
    std::copy(x.begin(), x.end(), z.begin());
    std::copy(y.begin(), y.end(), z.begin() + x.size());
    
    Rcpp::StringVector zconv(z.begin(), z.end());
    
    return zconv;
  }
  
  //' Compares Two Strings Literally
  //' 
  //' This function compares two strings element by element. Returns \code{FALSE}
  //' in case of any differences whatsoever.
  //' 
  //' @name stringcompare_hard
  //' 
  //' @param str1 The first string
  //' @param str2 The second string
  //' 
  //' @return A logical value. In case of any difference at all, it will return
  //' \code{FALSE}.
  //' 
  //' @keywords internal
  //' @noRd
  inline bool stringcompare_hard(std::string str1, std::string str2) {
    int str1_length = str1.size();
    int str2_length = str2.size();
    bool same = true;
    
    if (str1_length == str2_length && str1_length > 0) {
      for (int i = 0; i < str1_length; i++) {
        if (str1[i] != str2[i]) {
          same = false;
        }
      }
    } else if (str1_length != str2_length) {
      same = false;
    }
    
    return same;
  }
  
  //' Compares Two Strings, Assessing Inclusion
  //' 
  //' This function compares two strings, and will assess whether \code{str2} is
  //' contained within \code{str1}.
  //' 
  //' @name stringcompare_soft
  //' 
  //' @param str1 The first string
  //' @param str2 The second string
  //' 
  //' @return A list of two values. The first is a logical value indicating
  //' whether \code{str2} occurs within \code{str1}. The second element is an
  //' integer indicating at what element of \code{str1} \code{str2} begins.
  //' \code{FALSE}.
  //' 
  //' @keywords internal
  //' @noRd
  inline List stringcompare_soft(std::string str1, std::string str2) {
    int str1_length = str1.size();
    int str2_length = str2.size();
    int rem_check {0};
    bool same = false;
    unsigned int start_index {0};
    
    if (str1_length >= str2_length && str2_length > 0) {
      for (int i = 0; i < str1_length; i++) {
        if (str1[i] != str2[rem_check]) {
          rem_check = 0;
          // same = false;
        } else {
          if (rem_check == 0) start_index = i;
          rem_check += 1;
          if (rem_check >= str2_length) break;
        }
      }
      
      if (rem_check == str2_length) {
        same = true;
      }
    }
    
    Rcpp::List output = Rcpp::List::create(_["contains"] = same,
      _["start_index"] = start_index);
    
    return output;
  }
  
  //' Compares Two Strings, Assessing Inclusion
  //' 
  //' This function compares two strings, and will assess whether \code{str2} is
  //' contained within \code{str1}. It is a simpler version of 
  //' \code{stringcompare_soft()} that yields only the logical result.
  //' 
  //' @name stringcompare_simple
  //' 
  //' @param str1 The first string
  //' @param str2 The second string
  //' @param lower A logical value indicating whether to change all inputs to
  //' lower case before checking.
  //' 
  //' @return A logical value indicating whether \code{str2} occurs within
  //' \code{str1}.
  //' 
  //' @keywords internal
  //' @noRd
  inline bool stringcompare_simple(std::string str1, std::string str2,
    bool lower = false) {
    
    int str1_length = str1.size();
    int str2_length = str2.size();
    int rem_check {0};
    bool same = false;
    
    if (str1_length >= str2_length && str2_length > 0) {
      for (int i = 0; i < str1_length; i++) {
        if (!lower) {
          if (str1[i] != str2[rem_check]) {
            rem_check = 0;
          } else {
            rem_check += 1;
            if (rem_check >= str2_length) break;
          }
        } else {
          if (tolower(str1[i]) != tolower(str2[rem_check])) {
            rem_check = 0;
          } else {
            rem_check += 1;
            if (rem_check >= str2_length) break;
          }
        }
      }
      
      if (rem_check == str2_length) {
        same = true;
      }
    }
    
    return same;
  }
  
  //' Compares Three Strings for Interaction Notation
  //' 
  //' This function compares checks to see if one string is composed of the other
  //' two strings in R's interaction notation.
  //' 
  //' @name stringcompare_x
  //' 
  //' @param str1 The first string. Used for comparison.
  //' @param str2 The second string. Will be incorporated into interaction format.
  //' @param str3 The third string. Will be incorporated into interaction format.
  //' 
  //' @return A logical value. In case of any difference at all, it will return
  //' \code{FALSE}.
  //' 
  //' @keywords internal
  //' @noRd
  inline bool stringcompare_x(std::string str1, std::string str2,
    std::string str3) {
    
    int str1_length = str1.size();
    int str2_length = str2.size();
    int str3_length = str3.size();
    int combined_length = str2_length + str3_length + 1;
    bool same = false;
    bool same1 = true;
    bool same2 = true;
    
    if (str1_length == combined_length && str1_length > 0) {
      std::string x1 = str2;
      x1 += ":";
      x1 += str3;
      
      std::string x2 = str3;
      x2 += ":";
      x2 += str2;
      
      for (int i = 0; i < str1_length; i++) {
        if (str1[i] != x1[i]) {
          same1 = false;
        }
        if (str1[i] != x2[i]) {
          same2 = false;
        }
      }
    } else {
      same1 = false;
      same2 = false;
    }
    
    if (same1 || same2) same = true;
    
    return same;
  }
  
  //' Sort String Elements
  //' 
  //' This function is based on code obtained from R Bloggers
  //' (see https://www.r-bloggers.com/2013/01/handling-strings-with-rcpp/). It
  //' sorts the elements of a string vector in alphabetical order.
  //' 
  //' @name stringsort
  //' 
  //' @param string_input A string vector.
  //' 
  //' @return The sorted string vector.
  //' 
  //' @keywords internal
  //' @noRd
  inline CharacterVector stringsort(const CharacterVector& string_input ) {
    int len = string_input.size();
    
    std::vector<std::string> converted(len);
    for (int i=0; i < len; i++) converted[i] = as<std::string>(string_input(i));
    std::sort( converted.begin(), converted.end() );
    
    Rcpp::CharacterVector new_converted(len);
    new_converted = converted;
    
    return new_converted;
  }
  
  //' Sort Integer Elements
  //' 
  //' This function is based on code obtained from the Rcpp Gallery by Ross
  //' Bennett (see https://gallery.rcpp.org/articles/sorting/). It sorts the
  //' elements of an integer vector.
  //' 
  //' @name int_sort
  //' 
  //' @param int_input An integer vector.
  //' 
  //' @return The sorted integer vector.
  //' 
  //' @keywords internal
  //' @noRd
  inline IntegerVector int_sort(const IntegerVector& x) {
     Rcpp::IntegerVector y = clone(x);
     std::sort(y.begin(), y.end());
     
     return y;
  }
  
  //' Index a Numeric Matrix According to a Reference Vector
  //' 
  //' Function \code{refsort_num()} takes a numeric matrix and replaces it with an
  //' integer matrix showing the position of each element in the input matrix
  //' within the reference vector.
  //' 
  //' @name refsort_num
  //' 
  //' @param vec The matrix to index
  //' @param ref The vector to use as a reference
  //' 
  //' @return An integer matrix with integers referring to elements in vector
  //' \code{ref}.
  //' 
  //' @keywords internal
  //' @noRd
  inline IntegerMatrix refsort_num(const NumericMatrix& vec,
    const NumericVector& ref) {
    int vec_length = vec.length();
    int ref_length = ref.length();
    
    Rcpp::IntegerMatrix output(vec.nrow(), vec.ncol());
    
    for (int i = 0; i < vec_length; i++) {
      for (int j = 0; j < ref_length; j++) {
        if (vec[i] == ref[j]) output[i] = j + 1;
      }
    }
    
    return output;
  }
  
  //' Index a String Vector According to a Reference Vector
  //' 
  //' Function \code{refsort_str()} takes a string vector and replaces it with an
  //' integer vector showing the position of each element in the input vector
  //' within the reference vector.
  //' 
  //' @name refsort_str
  //' 
  //' @param vec The vector to index
  //' @param ref The vector to use as a reference
  //' 
  //' @return An integer vector with integers referring to elements in vector
  //' \code{ref}.
  //' 
  //' @keywords internal
  //' @noRd
  inline IntegerVector refsort_str(const CharacterVector& vec,
    const CharacterVector& ref) {
    int vec_length = vec.length();
    int ref_length = ref.length();
    
    Rcpp::IntegerVector output(vec_length);
    
    for (int i = 0; i < vec_length; i++) {
      for (int j = 0; j < ref_length; j++) {
        if (stringcompare_hard(as<std::string>(vec[i]), 
          as<std::string>(ref[j]))) output[i] = j + 1;
      }
    }
    
    return output;
  }
  
  //' Index a String Matrix According to a Reference Vector
  //' 
  //' Function \code{refsort_str_m()} takes a string matrix and replaces it with
  //' an integer matrix showing the position of each element in the input matrix
  //' within the reference vector.
  //' 
  //' @name refsort_str_m
  //' 
  //' @param vec The vector to index
  //' @param ref The vector to use as a reference
  //' 
  //' @return An integer matrix with integers referring to elements in vector
  //' \code{ref}.
  //' 
  //' @keywords internal
  //' @noRd
  inline IntegerMatrix refsort_str_m(const CharacterMatrix& vec,
    const CharacterVector& ref) {
    int vec_length = vec.length();
    int ref_length = ref.length();
    
    Rcpp::IntegerMatrix output(vec.nrow(), vec.ncol());
    
    for (int i = 0; i < vec_length; i++) {
      for (int j = 0; j < ref_length; j++) {
        if (stringcompare_hard(as<std::string>(vec[i]), 
          as<std::string>(ref[j]))) output[i] = j + 1;
      }
    }
    
    return output;
  }

  //' Calculate Logarithms of Non-Zero Elements of Sparse Matrix
  //' 
  //' Function \code{spmat_log} finds the non-zero elements in a sparse matrix,
  //' calculates their logs, and inserts them back into the matrix and returns it.
  //' Based on code developed by Coatless Professor and posted by him on
  //' StackOverflow.
  //' 
  //' @name spmat_log
  //' 
  //' @param B A sparse matrix. Note that this is assumed to be a population
  //' projection matrix, meaning that all values are either 0 or positive.
  //' 
  //' @return A sparse matrix with non-zero values as logs of the elements in the
  //' input matrix.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::sp_mat spmat_log(arma::sp_mat coremat)
  {
    arma::sp_mat::const_iterator start = coremat.begin();
    arma::sp_mat::const_iterator end   = coremat.end();
    arma::sp_mat::const_iterator it = start; 
    
    int n = std::distance(start, end);
    
    if (n > 0) {
      arma::umat locs(2, n);
      arma::uvec temp(2);
      arma::vec vals(n);
      arma::vec logvals(n);
      locs.zeros();
      temp.zeros();
      vals.zeros();
      logvals.zeros();
      
      for(int i = 0; i < n; ++i) {
        temp(0) = it.row();
        temp(1) = it.col();
        locs.col(i) = temp;
        vals(i) = coremat(temp(0), temp(1));
        logvals(i) = log(vals(i));
        ++it; // increment
      }
      
      coremat = arma::sp_mat(locs, logvals, coremat.n_rows, coremat.n_cols);
    }
    
    return coremat;
  }
  
  //' Resize an IntegerVector
  //' 
  //' This function resizes an IntegerVector. It is based on code provided by
  //' Dirk Eddelbuettel on StackExchange
  //' (https://stackoverflow.com/questions/13782943/how-to-resize-a-numericvector).
  //' 
  //' @name shrink
  //' 
  //' @param x_int An IntegerVector.
  //' 
  //' @return A new IntegerVector exactly one element smaller than \code{x_int}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::IntegerVector shrink(const IntegerVector& x_int) {
    arma::ivec x = as<arma::ivec>(x_int);
    
    arma::ivec y = x;
    y.resize(x.size() - 1);
    
    Rcpp::IntegerVector y_out = as<IntegerVector>(wrap(y));
    
    return y_out;
  }
  
  //' Find Indices of a Matching String in a StringVector
  //' 
  //' This function returns the C++ style vector indices of all elements in a
  //' StringVector matching a specific text.
  //' 
  //' @name index_l3
  //' 
  //' @param mainvec The StringVector to search through.
  //' @param target The string to find within \code{mainvec}.
  //' 
  //' @return An IntegerVector holding the indices of matching strings.
  //' 
  //' @keywords internal
  //' @noRd
  inline IntegerVector index_l3 (const StringVector& mainvec, String target) {
    int mainvec_length = static_cast<int>(mainvec.length());
    
    int matches {0};
    for (int i = 0; i < mainvec_length; i++) {
      if (stringcompare_hard(String(mainvec(i)), target)) matches++;
    }
    
    IntegerVector all_indices (matches);
    int  ai_elem_counter {0};
    for (int i = 0; i < mainvec_length; i++) {
      if (stringcompare_hard(String(mainvec(i)), target)) {
        all_indices(ai_elem_counter) = i;
        ai_elem_counter++;
        
        if (ai_elem_counter == matches) break;
      }
    }
    
    return all_indices;
  }
  
  //' Subset Data Frame By A Specific Variable
  //' 
  //' Returns a data frame subset of a single condition in a single variable.
  //' 
  //' @name df_subset
  //' 
  //' @param x The data frame to subset.
  //' @param level The value to search for of the variable to use for subsetting.
  //' @param cisNA A logical override for \code{level}, that tells R to use
  //' \code{NA} as the default subset value.
  //' @param equal Logical value indicating whether to test for equality. Defaults
  //' to \code{TRUE}.
  //' @param greater Logical value indicating whether to test for a greater value.
  //' Only used for integer and numeric variables. Defaults to \code{FALSE}.
  //' @param less Logical value indicating whether to test for a lower value. Only
  //' used for integer and numeric variables. Defaults to \code{FALSE}.
  //' @param c_ints A logical value indicating whether to use C++ encoding for
  //' indices. If \code{TRUE}, then indices start at 0. If \code{FALSE}, then
  //' they start at 1. Only applies to variable indices.
  //' @param var The name or number of the variable to use for subsetting.
  //' 
  //' @return A new data frame subset from the old.
  //' 
  //' @section Notes:
  //' Either of \code{greater} or \code{less} must always be set to \code{FALSE}.
  //' 
  //' String and logical variables are always subset by equality.
  //' 
  //' @keywords internal
  //' @noRd
  inline List df_subset(const DataFrame& x, RObject level, bool cisNA = false,
    bool equal = true, bool greater = false, bool less = false,
    bool c_ints = true, Nullable<RObject> var = R_NilValue) {
    
    if (greater && less) {
      throw Rcpp::exception("Criteria cannot include both less than and greater than.",
        false);
    }
    
    StringVector var_names = x.attr("names");
    StringVector df_class = x.attr("class");
    int no_vars = x.length();
    int no_rows = x.nrows();
    
    int chosen_var_i {-1};
    
    if (var.isNotNull()) {
      RObject var_ = as<RObject>(var);
      
      if (is<StringVector>(var_)) {
        StringVector full_var = as<StringVector>(var_);
        String var_s_ = String(full_var(0));
        
        for (int i = 0; i < no_vars; i++) {
          if (stringcompare_hard(var_s_, String(var_names(i)))) {
            chosen_var_i = i;
          }
        }
      } else if (is<IntegerVector>(var_) || is<NumericVector>(var_)) {
        IntegerVector var_i_ = as<IntegerVector>(var_);
        
        if (var_i_.length() != 1) {
          throw Rcpp::exception("Please enter only a single variable to subset on.",
            false);
        }
        
        if (c_ints) {
          chosen_var_i = var_i_(0);
        } else {
          chosen_var_i = var_i_(0) - 1;
        }
        
        if (chosen_var_i < 0 || chosen_var_i > (no_vars - 1)) {
          throw Rcpp::exception("Chosen variable number outside the range of the data frame.",
            false);
        }
      }
    }
    
    if (chosen_var_i == -1) {
      throw Rcpp::exception("No valid variable chosen for subsetting.", false);
    }
    
    RObject chosen_data = x[chosen_var_i];
    
    arma::uvec useable_indices;
    double chosen_level_d {0.0};
    String chosen_level_s;
    bool string_used {false};
    
    if (is<StringVector>(level)) {
      StringVector level_vec = as<StringVector>(level);
      
      if (!StringVector::is_na(level_vec(0))) {
        chosen_level_s = String(level_vec(0));
        string_used = true;
      } else {
        cisNA = true;
      }
      
    } else if (is<NumericVector>(level)) {
      NumericVector level_vec = as<NumericVector>(level);
      
      if (!NumericVector::is_na(level_vec(0))) {
        chosen_level_d = static_cast<double>(level_vec(0));
      } else {
        cisNA = true;
      }
      
    } else if (is<LogicalVector>(level)) {
      LogicalVector level_vec = as<LogicalVector>(level);
      
      if (!LogicalVector::is_na(level_vec(0))) {
        chosen_level_d = static_cast<double>(level_vec(0));
      } else {
        cisNA = true;
      }
      
    } else if (is<IntegerVector>(level)) {
      IntegerVector level_vec = as<IntegerVector>(level);
      
      if (!IntegerVector::is_na(level_vec(0))) {
        chosen_level_d = static_cast<double>(level_vec(0));
      } else {
        cisNA = true;
      }
    }
    
    if (is<NumericVector>(chosen_data) || is<IntegerVector>(chosen_data)) {
      arma::vec chosen_data_ = as<arma::vec>(chosen_data);
      
      if (string_used) {
        throw Rcpp::exception("String levels cannot be used to subset numeric/integer variables.",
          false);
      }
      
      if (cisNA) { 
        NumericVector chosen_data_n = as<NumericVector>(chosen_data);
        arma::uvec na_check (no_rows, fill::zeros);
        
        for (int i = 0; i < no_rows; i++) {
          if (NumericVector::is_na(chosen_data_n(i))) na_check(i) = 1;
        }
        useable_indices = find(na_check);
      } else if (equal && !greater && !less) {
        useable_indices = find(chosen_data_ == chosen_level_d);
      } else if (equal && greater && !less) {
        useable_indices = find(chosen_data_ >= chosen_level_d);
      } else if (equal && !greater && less) {
        useable_indices = find(chosen_data_ <= chosen_level_d);
      } else if (!equal && greater && !less) {
        useable_indices = find(chosen_data_ > chosen_level_d);
      } else if (!equal && !greater && less) {
        useable_indices = find(chosen_data_ < chosen_level_d);
      }
      
    } else if (is<LogicalVector>(chosen_data)) {
      arma::uvec chosen_data_ = as<arma::uvec>(chosen_data);
      
      if (string_used) {
        throw Rcpp::exception("String levels cannot be used to subset logical variables.",
          false);
      }
      
      arma::uvec level_vec = as<arma::uvec>(level);
      
      if (cisNA) { 
        LogicalVector chosen_data_l = as<LogicalVector>(chosen_data);
        arma::uvec na_check (no_rows, fill::zeros);
        
        for (int i = 0; i < no_rows; i++) {
          if (LogicalVector::is_na(chosen_data_l(i))) na_check(i) = 1;
        }
        useable_indices = find(na_check);
      } else {
        useable_indices = find(chosen_data_ == chosen_level_d);
      }
      
    } else if (is<StringVector>(chosen_data)) {
      StringVector chosen_data_ = as<StringVector>(chosen_data);
      
      arma::uvec cv_log (no_rows, fill::zeros);
      
      if (cisNA) { 
        for (int i = 0; i < no_rows; i++) {
          if (StringVector::is_na(chosen_data_(i))) cv_log(i) = 1;
        }
      } else {
        for (int i = 0; i < no_rows; i++) {
          if (stringcompare_hard(String(chosen_data_(i)), chosen_level_s)) {
            cv_log(i) = 1;
          }
        }
      }
      useable_indices = find(cv_log); 
      
    } else {
      throw Rcpp::exception("Chosen variable is not a recognized subsettable type.",
        false);
    }
    
    int new_rows = static_cast<int>(useable_indices.n_elem);
    List new_df (no_vars);
    
    for (int i = 0; i < no_vars; i++) {
      if (is<NumericVector>(x[i])) {
        NumericVector old_var_i = as<NumericVector>(x[i]);
        NumericVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<IntegerVector>(x[i])) {
        IntegerVector old_var_i = as<IntegerVector>(x[i]);
        IntegerVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        
        if (old_var_i.hasAttribute("levels")) {
          StringVector int_class = old_var_i.attr("class");
          if (stringcompare_simple(String(int_class(0)), "fact", false)) {
            new_var_i.attr("levels") = old_var_i.attr("levels");
            new_var_i.attr("class") = "factor";
          }
        }
        
        new_df(i) = new_var_i;
        
      } else if (is<LogicalVector>(x[i])) {
        LogicalVector old_var_i = as<LogicalVector>(x[i]);
        LogicalVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<StringVector>(x[i])) {
        StringVector old_var_i = as<StringVector>(x[i]);
        StringVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        new_df(i) = new_var_i;
        
      } else {
        throw Rcpp::exception("Variable found of unrecognized type.", false);
      }
    }
    
    new_df.attr("names") = var_names;
    new_df.attr("class") = df_class;
    
    StringVector row_names(new_rows);
    for (int i = 0; i < new_rows; i++) {
      row_names(i) = std::to_string(i+1);
    }
    new_df.attr("row.names") = row_names;
    
    return new_df;
  }
  
  //' Check If Two Data Frames Are Equal
  //' 
  //' Function \code{df_compare()} assesses whether two data frames are the same.
  //' 
  //' @name df_compare
  //' 
  //' @param x First data frame.
  //' @param y Second data frame.
  //' 
  //' @return Returns \code{TRUE} if the data frames are the same; otherwise,
  //' returns \code{FALSE}.
  //' @keywords internal
  //' @noRd
  inline bool df_compare(const DataFrame& x, const DataFrame& y) {
    // bool df_equal {false};
    bool all_equal {true};
    
    int x_length = x.length();
    int y_length = y.length();
    
    if (x_length == y_length) {
      CharacterVector x_varnames = x.attr("names");
      CharacterVector y_varnames = y.attr("names");
      
      for (int i = 0; i < x_length; i++) {
        RObject x_i = RObject(x(i));
        RObject y_i = RObject(y(i));
        
        if (x_varnames(i) != y_varnames(i)) {
          all_equal = false;
          break;
        }
        
        if (is<LogicalVector>(x_i) && !is<LogicalVector>(y_i)) {
          all_equal = false;
          break;
        }
        if (is<IntegerVector>(x_i) && !is<IntegerVector>(y_i)) {
          all_equal = false;
          break;
        }
        if (is<NumericVector>(x_i) && !is<NumericVector>(y_i)) {
          all_equal = false;
          break;
        }
        if (is<CharacterVector>(x_i) && !is<CharacterVector>(y_i)) {
          all_equal = false;
          break;
        }
        
        if (is<LogicalVector>(x_i)) {
          // Rcout << "LogicalVector" << endl;
          
          LogicalVector x_i_log = LogicalVector(x_i);
          LogicalVector y_i_log = LogicalVector(y_i);
          
          for (int j = 0; j < x_i_log.length(); j++) {
            if (LogicalVector::is_na(x_i_log(j)) && LogicalVector::is_na(y_i_log(j))) {
              // No worries here
            } else if (x_i_log(j) != y_i_log(j)) {
              // Rcout << "Different" << endl;
              all_equal = false;
              break;
            }
          }
          
        } else if (is<IntegerVector>(x_i)) {
          // Rcout << "IntegerVector" << endl;
          
          IntegerVector x_i_int = IntegerVector(x_i);
          IntegerVector y_i_int = IntegerVector(y_i);
          
          for (int j = 0; j < static_cast<int>(x_i_int.length()); j++) {
            if (IntegerVector::is_na(x_i_int(j)) && IntegerVector::is_na(y_i_int(j))) {
              // No worries here
            } else if (x_i_int(j) != y_i_int(j)) {
              // Rcout << "Different" << endl;
              all_equal = false;
              break;
            }
          }
          
          if (x_i_int.hasAttribute("levels") && !y_i_int.hasAttribute("levels")) {
            all_equal = false;
            break;
          }
          if (!x_i_int.hasAttribute("levels") && y_i_int.hasAttribute("levels")) {
            all_equal = false;
            break;
          }
          
          if (x_i_int.hasAttribute("levels")) {
            CharacterVector x_i_levels = CharacterVector(x_i_int.attr("levels"));
            CharacterVector y_i_levels = CharacterVector(y_i_int.attr("levels"));
            
            if (x_i_levels.length() != y_i_levels.length()) {
              all_equal = false;
              break;
            }
            
            for (int j = 0; j < static_cast<int>(x_i_levels.length()); j++) {
              if (x_i_levels(j) != y_i_levels(j)) {
                all_equal = false;
                break;
              }
            }
          }
          
        } else if (is<NumericVector>(x_i)) {
          // Rcout << "NumericVector" << endl;
          
          NumericVector x_i_num = NumericVector(x_i);
          NumericVector y_i_num = NumericVector(y_i);
          
          for (int j = 0; j < x_i_num.length(); j++) {
            if (NumericVector::is_na(x_i_num(j)) && NumericVector::is_na(y_i_num(j))) {
              // No worries here
            } else if (x_i_num(j) != y_i_num(j)) {
              // Rcout << "Different" << endl;
              all_equal = false;
              break;
            }
          }
          
        } else if (is<CharacterVector>(x_i)) {
          // Rcout << "CharacterVector" << endl;
          
          CharacterVector x_i_str = CharacterVector(x_i);
          CharacterVector y_i_str = CharacterVector(y_i);
          
          for (int j = 0; j < x_i_str.length(); j++) {
            if (CharacterVector::is_na(x_i_str(j)) && CharacterVector::is_na(y_i_str(j))) {
              // No worries here
            } else if (x_i_str(j) != y_i_str(j)) {
              // Rcout << "Different" << endl;
              all_equal = false;
              break;
            }
          }
          
        } else {
          throw Rcpp::exception("Unknown variable type in input data frame.", false);
        }
      }
    } else all_equal = false;
    
    return(all_equal);
  }
  
  //' Subset Data Frame By Integer Vector of Row Numbers
  //' 
  //' Returns a data frame subset of a single condition in a single variable.
  //' 
  //' @name df_subset_byrow
  //' 
  //' @param x The data frame to subset.
  //' @param rows_to_keep An IntegerVector giving the rows to include, in order.
  //' Entries do not need to be in order, nor do they need to be unique. Uses
  //' C++ indexing..
  //' 
  //' @return A new data frame subset from the old.
  //' 
  //' @keywords internal
  //' @noRd
  inline List df_subset_byrow(const DataFrame& x, IntegerVector rows_to_keep) {
    
    StringVector var_names = x.attr("names");
    StringVector df_class = x.attr("class");
    int no_vars = x.length();
    //int no_rows = x.nrows();
    
    int new_rows = static_cast<int>(rows_to_keep.length());
    List new_df (no_vars);
    
    for (int i = 0; i < no_vars; i++) {
      if (is<NumericVector>(x[i])) {
        NumericVector old_var_i = as<NumericVector>(x[i]);
        NumericVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(rows_to_keep(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<IntegerVector>(x[i])) {
        IntegerVector old_var_i = as<IntegerVector>(x[i]);
        IntegerVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(rows_to_keep(j));
        }
        
        if (old_var_i.hasAttribute("levels")) {
          StringVector int_class = old_var_i.attr("class");
          if (stringcompare_simple(String(int_class(0)), "fact", false)) {
            new_var_i.attr("levels") = old_var_i.attr("levels");
            new_var_i.attr("class") = "factor";
          }
        }
        
        new_df(i) = new_var_i;
        
      } else if (is<LogicalVector>(x[i])) {
        LogicalVector old_var_i = as<LogicalVector>(x[i]);
        LogicalVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(rows_to_keep(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<StringVector>(x[i])) {
        StringVector old_var_i = as<StringVector>(x[i]);
        StringVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(rows_to_keep(j));
        }
        new_df(i) = new_var_i;
        
      } else {
        throw Rcpp::exception("Variable found of unrecognized type.", false);
      }
    }
    
    new_df.attr("names") = var_names;
    new_df.attr("class") = df_class;
    
    StringVector row_names(new_rows);
    for (int i = 0; i < new_rows; i++) {
      row_names(i) = std::to_string(i+1);
    }
    new_df.attr("row.names") = row_names;
    
    return new_df;
  }
  
  //' Remove Rows With Specific Index Values From Data Frame
  //' 
  //' Returns a data frame subset of a single condition in a single variable.
  //' 
  //' @name df_remove
  //' 
  //' @param x The data frame to remove data rows from.
  //' @param level The value to search for of the variable to use for subsetting.
  //' @param cisNA A logical override for \code{level}, that tells R to use
  //' \code{NA} as the default subset value.
  //' @param equal Logical value indicating whether to test for equality. Defaults
  //' to \code{TRUE}.
  //' @param greater Logical value indicating whether to test for a greater value.
  //' Only used for integer and numeric variables. Defaults to \code{FALSE}.
  //' @param less Logical value indicating whether to test for a lower value. Only
  //' used for integer and numeric variables. Defaults to \code{FALSE}.
  //' @param c_ints A logical value indicating whether to use C++ encoding for
  //' indices. If \code{TRUE}, then indices start at 0. If \code{FALSE}, then
  //' they start at 1. Only applies to variable indices.
  //' @param var The name or number of the variable to use for subsetting.
  //' 
  //' @return A new data frame based on the old, but without data rows found
  //' corresponding to the chosen \code{level} of the variable used for
  //' indexing.
  //' 
  //' @section Notes:
  //' Either of \code{greater} or \code{less} must always be set to \code{FALSE}.
  //' 
  //' String and logical variables are always subset by equality.
  //' 
  //' @keywords internal
  //' @noRd
  inline List df_remove(const DataFrame& x, RObject level, bool cisNA = false,
    bool equal = true, bool greater = false, bool less = false,
    bool c_ints = true, Nullable<RObject> var = R_NilValue) {
    
    if (greater && less) {
      throw Rcpp::exception("Criteria cannot include both less and greater than.",
        false);
    }
    
    StringVector var_names = x.attr("names");
    StringVector df_class = x.attr("class");
    int no_vars = x.length();
    int no_rows = x.nrows();
    
    int chosen_var_i {-1};
    
    if (var.isNotNull()) {
      RObject var_ = as<RObject>(var);
      
      if (is<StringVector>(var_)) {
        StringVector full_var = as<StringVector>(var_);
        String var_s_ = String(full_var(0));
        
        for (int i = 0; i < no_vars; i++) {
          if (stringcompare_hard(var_s_, String(var_names(i)))) {
            chosen_var_i = i;
          }
        }
      } else if (is<IntegerVector>(var_) || is<NumericVector>(var_)) {
        IntegerVector var_i_ = as<IntegerVector>(var_);
        
        if (var_i_.length() != 1) {
          throw Rcpp::exception("Enter only a single variable to subset on.",
            false);
        }
        
        if (c_ints) {
          chosen_var_i = var_i_(0);
        } else {
          chosen_var_i = var_i_(0) - 1;
        }
        
        if (chosen_var_i < 0 || chosen_var_i > (no_vars - 1)) {
          throw Rcpp::exception("Chosen variable number outside the range of the data frame.",
            false);
        }
      }
    }
    
    if (chosen_var_i == -1) {
      throw Rcpp::exception("No valid variable chosen for subsetting.", false);
    }
    
    RObject chosen_data = x[chosen_var_i];
    
    arma::uvec useable_indices;
    double chosen_level_d {0.0};
    String chosen_level_s;
    bool string_used {false};
    
    if (is<StringVector>(level)) {
      StringVector level_vec = as<StringVector>(level);
      
      if (!StringVector::is_na(level_vec(0))) {
        chosen_level_s = String(level_vec(0));
        string_used = true;
      } else {
        cisNA = true;
      }
      
    } else if (is<NumericVector>(level)) {
      NumericVector level_vec = as<NumericVector>(level);
      
      if (!NumericVector::is_na(level_vec(0))) {
        chosen_level_d = static_cast<double>(level_vec(0));
      } else {
        cisNA = true;
      }
      
    } else if (is<LogicalVector>(level)) {
      LogicalVector level_vec = as<LogicalVector>(level);
      
      if (!LogicalVector::is_na(level_vec(0))) {
        chosen_level_d = static_cast<double>(level_vec(0));
      } else {
        cisNA = true;
      }
      
    } else if (is<IntegerVector>(level)) {
      IntegerVector level_vec = as<IntegerVector>(level);
      
      if (!IntegerVector::is_na(level_vec(0))) {
        chosen_level_d = static_cast<double>(level_vec(0));
      } else {
        cisNA = true;
      }
    }
    
    if (is<NumericVector>(chosen_data) || is<IntegerVector>(chosen_data)) {
      arma::vec chosen_data_ = as<arma::vec>(chosen_data);
      
      if (string_used) {
        throw Rcpp::exception("String levels cannot be used to subset numeric/integer variables.",
          false);
      }
      
      if (cisNA) { 
        NumericVector chosen_data_n = as<NumericVector>(chosen_data);
        arma::uvec na_check (no_rows, fill::ones);
        
        for (int i = 0; i < no_rows; i++) {
          if (NumericVector::is_na(chosen_data_n(i))) na_check(i) = 0;
        }
        useable_indices = find(na_check);
      } else if (equal && !greater && !less) {
        useable_indices = find(chosen_data_ != chosen_level_d);
      } else if (equal && greater && !less) {
        useable_indices = find(chosen_data_ < chosen_level_d);
      } else if (equal && !greater && less) {
        useable_indices = find(chosen_data_ > chosen_level_d);
      } else if (!equal && greater && !less) {
        useable_indices = find(chosen_data_ <= chosen_level_d);
      } else if (!equal && !greater && less) {
        useable_indices = find(chosen_data_ >= chosen_level_d);
      }
      
    } else if (is<LogicalVector>(chosen_data)) {
      arma::uvec chosen_data_ = as<arma::uvec>(chosen_data);
      
      if (string_used) {
        throw Rcpp::exception("String levels cannot be used to subset logical variables.",
          false);
      }
      
      arma::uvec level_vec = as<arma::uvec>(level);
      
      if (cisNA) { 
        LogicalVector chosen_data_l = as<LogicalVector>(chosen_data);
        arma::uvec na_check (no_rows, fill::ones);
        
        for (int i = 0; i < no_rows; i++) {
          if (LogicalVector::is_na(chosen_data_l(i))) na_check(i) = 0;
        }
        useable_indices = find(na_check);
      } else {
        useable_indices = find(chosen_data_ != chosen_level_d);
      }
      
    } else if (is<StringVector>(chosen_data)) {
      StringVector chosen_data_ = as<StringVector>(chosen_data);
      
      arma::uvec cv_log (no_rows, fill::ones);
      
      if (cisNA) { 
        for (int i = 0; i < no_rows; i++) {
          if (StringVector::is_na(chosen_data_(i))) cv_log(i) = 0;
        }
      } else {
        for (int i = 0; i < no_rows; i++) {
          if (stringcompare_hard(String(chosen_data_(i)), chosen_level_s)) {
            cv_log(i) = 0;
          }
        }
      }
      useable_indices = find(cv_log); 
      
    } else {
      throw Rcpp::exception("Chosen variable is not a recognized subsettable type.", false);
    }
    
    int new_rows = static_cast<int>(useable_indices.n_elem);
    List new_df (no_vars);
    
    for (int i = 0; i < no_vars; i++) {
      if (is<NumericVector>(x[i])) {
        NumericVector old_var_i = as<NumericVector>(x[i]);
        NumericVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<IntegerVector>(x[i])) {
        IntegerVector old_var_i = as<IntegerVector>(x[i]);
        IntegerVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        
        if (old_var_i.hasAttribute("levels")) {
          StringVector int_class = old_var_i.attr("class");
          if (stringcompare_simple(String(int_class(0)), "fact", false)) {
            new_var_i.attr("levels") = old_var_i.attr("levels");
            new_var_i.attr("class") = "factor";
          }
        }
        
        new_df(i) = new_var_i;
        
      } else if (is<LogicalVector>(x[i])) {
        LogicalVector old_var_i = as<LogicalVector>(x[i]);
        LogicalVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<StringVector>(x[i])) {
        StringVector old_var_i = as<StringVector>(x[i]);
        StringVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(useable_indices(j));
        }
        new_df(i) = new_var_i;
        
      } else {
        throw Rcpp::exception("Variable found of unrecognized type.", false);
      }
    }
    
    new_df.attr("names") = var_names;
    new_df.attr("class") = df_class;
    
    StringVector row_names(new_rows);
    for (int i = 0; i < new_rows; i++) {
      row_names(i) = std::to_string(i+1);
    }
    new_df.attr("row.names") = row_names;
    
    return new_df;
  }
  
  //' Shrink Data Frame According to Index Vector
  //' 
  //' This function will remove rows corresponding to the input vector.
  //' 
  //' @name df_shedrows
  //' 
  //' @param data The data frame to modify.
  //' @param index The vector of C++ format row indices to remove.
  //' 
  //' @return A data frame with rows removed from the original.
  //' 
  //' @keywords internal
  //' @noRd
  inline List df_shedrows (const DataFrame& data, const IntegerVector& index) {
    int no_vars = data.length();
    int old_rows = static_cast<int>(data.nrows());
    
    StringVector var_names = data.attr("names");
    StringVector df_class = data.attr("class");
    
    IntegerVector keep_index = seq(0, (old_rows - 1));
    keep_index = setdiff(keep_index, index);
    keep_index = sort_unique(keep_index);
    
    int new_rows = static_cast<int>(keep_index.length());
    List new_df (no_vars);
    
    for (int i = 0; i < no_vars; i++) {
      if (is<NumericVector>(data[i])) {
        NumericVector old_var_i = as<NumericVector>(data[i]);
        NumericVector new_var_i (new_rows);
        
        for (int j = 0; j < new_rows; j++) {
          new_var_i(j) = old_var_i(keep_index(j));
        }
        new_df(i) = new_var_i;
        
      } else if (is<IntegerVector>(data[i])) {
          IntegerVector old_var_i = as<IntegerVector>(data[i]);
          IntegerVector new_var_i (new_rows);
          
          for (int j = 0; j < new_rows; j++) {
            new_var_i(j) = old_var_i(keep_index(j));
          }
          
          if (old_var_i.hasAttribute("levels")) {
            StringVector int_class = old_var_i.attr("class");
            if (stringcompare_simple(String(int_class(0)), "fact", false)) {
              new_var_i.attr("levels") = old_var_i.attr("levels");
              new_var_i.attr("class") = "factor";
            }
          }
          
          new_df(i) = new_var_i;
          
        } else if (is<LogicalVector>(data[i])) {
          LogicalVector old_var_i = as<LogicalVector>(data[i]);
          LogicalVector new_var_i (new_rows);
          
          for (int j = 0; j < new_rows; j++) {
            new_var_i(j) = old_var_i(keep_index(j));
          }
          new_df(i) = new_var_i;
          
        } else if (is<StringVector>(data[i])) {
          StringVector old_var_i = as<StringVector>(data[i]);
          StringVector new_var_i (new_rows);
          
          for (int j = 0; j < new_rows; j++) {
            new_var_i(j) = old_var_i(keep_index(j));
          }
          new_df(i) = new_var_i;
          
        } else {
          throw Rcpp::exception("Variable found of unrecognized type.", false);
        }
      }
      
      new_df.attr("names") = var_names;
      new_df.attr("class") = df_class;
      
      StringVector row_names(new_rows);
      for (int i = 0; i < new_rows; i++) {
        row_names(i) = std::to_string(i+1);
      }
      new_df.attr("row.names") = row_names;
      
      return new_df;
  }
  
  //' Search for Duplicate Data Frame Values
  //' 
  //' This function returns provides a logical test of duplicated data points in
  //' a data frame.
  //' 
  //' @name df_duplicates
  //' 
  //' @param x A valid data frame.
  //' 
  //' @return \code{TRUE} if a single duplicate value is found; otherwise,
  //' \code{FALSE}.
  //' 
  //' @keywords internal
  //' @noRd
  inline bool df_duplicates(const DataFrame& x) {
    StringVector df_names = x.attr("names");
    int df_var_no = static_cast<int>(df_names.length());
    int df_rows_no = static_cast<int>(x.nrows());
    
    IntegerVector var_type (df_var_no);
    
    bool final_check {false};
    bool break_i {false};
    bool break_j {false};
    
    for (int i = 0; i < (df_rows_no - 1); i++) {
      arma::uvec old_check (df_rows_no, fill::zeros);
      arma::uvec old_ones (df_rows_no, fill::zeros);
      
      final_check = false;
      break_i = false;
      break_j = false;
      
      for (int j = 0; j < df_var_no; j++) {
        arma::uvec new_check (df_rows_no, fill::zeros);
        StringVector var_hold_str = as<StringVector>(x[j]);
        
        for (int k = i + 1; k < df_rows_no; k++) {
          if (j == 0) {
            if (StringVector::is_na(var_hold_str(i))) {
              if (StringVector::is_na(var_hold_str(k))) {
                old_check(k) = true;
              } else {
                old_check(k) = false;
              }
            } else {
              if (StringVector::is_na(var_hold_str(k))) {
                old_check(k) = false;
              } else {
                if (stringcompare_hard(String(var_hold_str(i)),
                  String(var_hold_str(k)))) {
                  old_check(k) = true;
                } else {
                  old_check(k) = false;
                }
              }
            }
          } else {
            if (StringVector::is_na(var_hold_str(i))) {
              if (StringVector::is_na(var_hold_str(k))) {
                new_check(k) = true;
              } else {
                new_check(k) = false;
              }
            } else {
              if (StringVector::is_na(var_hold_str(k))) {
                new_check(k) = false;
              } else {
                if (stringcompare_hard(String(var_hold_str(i)),
                  String(var_hold_str(k)))) {
                  new_check(k) = true;
                } else {
                  new_check(k) = false;
                }
              }
            }
          }
        }
        
        if (j == 0) {
          old_ones = find(old_check);
        } else {
          arma::uvec new_ones = find(new_check);
          arma::uvec old_new_intersect = intersect(old_ones, new_ones);
          
          int oni_length = static_cast<int>(old_new_intersect.n_elem);
          old_check.zeros();
          if (oni_length > 0) {
            old_check.elem(old_new_intersect) = ones<uvec>(oni_length);
            old_ones = find(old_check);
          } else {
            final_check = false;
            break_j = true;
            break;
          }
          
          if (j == (df_var_no - 1)) {
            if (oni_length > 0) {
              final_check = true;
              break_i = true;
            }
          }
          if (break_i || break_j) {
            break;
          }
        }
        if (break_i || break_j) break;
      }
      if (break_i) break;
    }
    
    return final_check;
  }
  
  //' Bind Two Data Frames By Row
  //' 
  //' This function takes two data frames, which must be composed of the same
  //' variables in the same order. Although the variables may have different
  //' names, the variables must be of the same type. The names of the variables in
  //' the new merged data frame will match those of the first data frame.
  //' Developed with the help of Microsoft Gemini AI.
  //' 
  //' @name df_rbind
  //' 
  //' @param df1 The first data frame, which will form the top of the new data
  //' frame and will be used as a reference for variable names.
  //' @param df2 The second data frame, which will be attached below data frame
  //' \code{df1}.
  //' 
  //' @return A new data frame composed of the merged data frames.
  //' 
  //' @keywords internal
  //' @noRd
  inline DataFrame df_rbind(DataFrame df1, DataFrame df2) {
    
    int nrows1 = df1.nrows();
    int nrows2 = df2.nrows();
    int ncols = df1.size();
  
    if (ncols != df2.size()) {
      stop("Data frames must have the same number of columns.");
    }
    
    CharacterVector df_names = df1.names();
  
    List out_df (ncols);
  
    for (int i = 0; i < ncols; ++i) {
      SEXP col1 = df1[i];
      SEXP col2 = df2[i];
  
      if (TYPEOF(col1) != TYPEOF(col2)) {
        stop("Columns must have the same type");
      }
  
      switch(TYPEOF(col1)){
        case INTSXP: {
          IntegerVector combinedCol(nrows1 + nrows2);
          IntegerVector col1Vector = as<IntegerVector>(col1);
          IntegerVector col2Vector = as<IntegerVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          bool current_int_class1 = col1Vector.hasAttribute("levels");
          if (current_int_class1) {
            CharacterVector col1lvls = col1Vector.attr("levels");
            CharacterVector col2lvls = col2Vector.attr("levels");
            
            CharacterVector col_mergedLevels = LefkoUtils::concat_str(col1lvls, col2lvls);
            CharacterVector unique_levels = unique(col_mergedLevels);
            
            combinedCol.attr("levels") = unique_levels;
            combinedCol.attr("class") = "factor";
          }
          out_df(i) = combinedCol;
          break;
        }
  
        case LGLSXP: {
          LogicalVector combinedCol(nrows1 + nrows2);
          LogicalVector col1Vector = as<LogicalVector>(col1);
          LogicalVector col2Vector = as<LogicalVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          out_df(i) = combinedCol;
          break;
        }
  
        case REALSXP: {
          NumericVector combinedCol(nrows1 + nrows2);
          NumericVector col1Vector = as<NumericVector>(col1);
          NumericVector col2Vector = as<NumericVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          out_df(i) = combinedCol;
          break;
        }
  
        case STRSXP: {
          CharacterVector combinedCol(nrows1 + nrows2);
          CharacterVector col1Vector = as<CharacterVector>(col1);
          CharacterVector col2Vector = as<CharacterVector>(col2);
          std::copy(col1Vector.begin(), col1Vector.end(), combinedCol.begin());
          std::copy(col2Vector.begin(), col2Vector.end(), combinedCol.begin() + nrows1);
          
          out_df(i) = combinedCol;
          break;
        }
        default:
          stop("Unsupported column type.");
      }
    }
    
    out_df.attr("names") = df_names;
    out_df.attr("class") = "data.frame";
    StringVector row_names(nrows1 + nrows2);
    for (int i = 0; i < (nrows1 + nrows2); i++) {
      row_names(i) = std::to_string(i+1);
    }
    out_df.attr("row.names") = row_names;
    
    return out_df;
  }
  
  //' Repeat First Vector for Each Element of Second Vector
  //' 
  //' This function is an Rcpp version of \code{expand.grid()}.
  //' 
  //' @param x Any R vector.
  //' @param y Any R vector.
  //' 
  //' @return A list with two elements. Each element is the expanded vector, and
  //' both vectors are of equal dimension (length of first vector multiplied by
  //' length of second vector).
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List exp_grd (RObject x, RObject y) {
    NumericVector x1d;
    IntegerVector x1i;
    LogicalVector x1l;
    StringVector x1s;
    
    NumericVector y1d;
    IntegerVector y1i;
    LogicalVector y1l;
    StringVector y1s;
    
    bool x_d {false};
    bool x_i {false};
    bool x_l {false};
    bool x_s {false};
    
    bool y_d {false};
    bool y_i {false};
    bool y_l {false};
    bool y_s {false};
    
    int x_length {0};
    if (is<NumericVector>(x)) {
      x1d = as<NumericVector>(x);
      x_d = true;
      x_length = static_cast<int>(x1d.length());
    } else if (is<IntegerVector>(x)) {
      x1i = as<IntegerVector>(x);
      x_i = true;
      x_length = static_cast<int>(x1i.length());
    } else if (is<LogicalVector>(x)) {
      x1l = as<LogicalVector>(x);
      x_l = true;
      x_length = static_cast<int>(x1l.length());
    } else if (is<StringVector>(x)) {
      x1s = as<StringVector>(x);
      x_s = true;
      x_length = static_cast<int>(x1s.length());
    } else {
      throw Rcpp::exception("Unrecognized variable type (exp_grd)", false);
    }
    
    int y_length {0};
    if (is<NumericVector>(y)) {
      y1d = as<NumericVector>(y);
      y_d = true;
      y_length = static_cast<int>(y1d.length());
    } else if (is<IntegerVector>(y)) {
      y1i = as<IntegerVector>(y);
      y_i = true;
      y_length = static_cast<int>(y1i.length());
    } else if (is<LogicalVector>(y)) {
      y1l = as<LogicalVector>(y);
      y_l = true;
      y_length = static_cast<int>(y1l.length());
    } else if (is<StringVector>(y)) {
      y1s = as<StringVector>(y);
      y_s = true;
      y_length = static_cast<int>(y1s.length());
    } else {
      throw Rcpp::exception("Unrecognized variable type (exp_grd)", false);
    }
    
    int new_length = x_length * y_length;
    
    List output (2);
    if (x_d && y_d) {
      NumericVector new_x (new_length);
      NumericVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1d(i);
          new_y( (i + (j * x_length)) ) = y1d(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_d && y_i) {
      NumericVector new_x (new_length);
      IntegerVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1d(i);
          new_y( (i + (j * x_length)) ) = y1i(j);
        }
      }
      
      if (y1i.hasAttribute("levels")) {
        StringVector y1i_levels = as<StringVector>(y1i);
        new_y.attr("levels") = y1i_levels;
        new_y.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_d && y_l) {
      NumericVector new_x (new_length);
      LogicalVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1d(i);
          new_y( (i + (j * x_length)) ) = y1l(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_d && y_s) {
      NumericVector new_x (new_length);
      StringVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1d(i);
          new_y( (i + (j * x_length)) ) = y1s(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_i && y_d) {
      IntegerVector new_x (new_length);
      NumericVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1i(i);
          new_y( (i + (j * x_length)) ) = y1d(j);
        }
      }
      
      if (x1i.hasAttribute("levels")) {
        StringVector x1i_levels = as<StringVector>(x1i);
        new_x.attr("levels") = x1i_levels;
        new_x.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_i && y_i) {
      IntegerVector new_x (new_length);
      IntegerVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1i(i);
          new_y( (i + (j * x_length)) ) = y1i(j);
        }
      }
      
      if (x1i.hasAttribute("levels")) {
        StringVector x1i_levels = as<StringVector>(x1i);
        new_x.attr("levels") = x1i_levels;
        new_x.attr("class") = "factor";
      }
      
      if (y1i.hasAttribute("levels")) {
        StringVector y1i_levels = as<StringVector>(y1i);
        new_y.attr("levels") = y1i_levels;
        new_y.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_i && y_l) {
      IntegerVector new_x (new_length);
      LogicalVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1i(i);
          new_y( (i + (j * x_length)) ) = y1l(j);
        }
      }
      
      if (x1i.hasAttribute("levels")) {
        StringVector x1i_levels = as<StringVector>(x1i);
        new_x.attr("levels") = x1i_levels;
        new_x.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_i && y_s) {
      IntegerVector new_x (new_length);
      StringVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1i(i);
          new_y( (i + (j * x_length)) ) = y1s(j);
        }
      }
      
      if (x1i.hasAttribute("levels")) {
        StringVector x1i_levels = as<StringVector>(x1i);
        new_x.attr("levels") = x1i_levels;
        new_x.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_l && y_d) {
      LogicalVector new_x (new_length);
      NumericVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1l(i);
          new_y( (i + (j * x_length)) ) = y1d(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_l && y_i) {
      LogicalVector new_x (new_length);
      IntegerVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1l(i);
          new_y( (i + (j * x_length)) ) = y1i(j);
        }
      }
      
      if (y1i.hasAttribute("levels")) {
        StringVector y1i_levels = as<StringVector>(y1i);
        new_y.attr("levels") = y1i_levels;
        new_y.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_l && y_l) {
      LogicalVector new_x (new_length);
      LogicalVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1l(i);
          new_y( (i + (j * x_length)) ) = y1l(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_l && y_s) {
      LogicalVector new_x (new_length);
      StringVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1l(i);
          new_y( (i + (j * x_length)) ) = y1s(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_s && y_d) {
      StringVector new_x (new_length);
      NumericVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1s(i);
          new_y( (i + (j * x_length)) ) = y1d(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_s && y_i) {
      StringVector new_x (new_length);
      IntegerVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1s(i);
          new_y( (i + (j * x_length)) ) = y1i(j);
        }
      }
      
      if (y1i.hasAttribute("levels")) {
        StringVector y1i_levels = as<StringVector>(y1i);
        new_y.attr("levels") = y1i_levels;
        new_y.attr("class") = "factor";
      }
      
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_s && y_l) {
      StringVector new_x (new_length);
      LogicalVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1s(i);
          new_y( (i + (j * x_length)) ) = y1l(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    } else if (x_s && y_s) {
      StringVector new_x (new_length);
      StringVector new_y (new_length);
      
      for (int i = 0; i < x_length; i++) {
        for (int j = 0; j < y_length; j++) {
          new_x( (i + (j * x_length)) ) = x1s(i);
          new_y( (i + (j * x_length)) ) = y1s(j);
        }
      }
      output(0) = new_x;
      output(1) = new_y;
      
    }
    
    return output;
  }
  
  //' Extract Key Components from Simple Numerical Model
  //' 
  //' This function creates a skeleton list needed for functions
  //' \code{jerzeibalowski()} and \code{motherbalowski()} when a vital rate model
  //' is simply a scalar.
  //' 
  //' @name numeric_extractor
  //' 
  //' @param object A numerical value, typical \code{1} or \code{0}.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}. Will generally be
  //' \code{numeric}.}
  //' \item{family}{The response distribution. Here, given as \code{constant}.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}. Here, given as \code{NULL}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.
  //' Here, given as \code{NULL}.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}. Here, given as
  //' \code{NULL}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables. Not used in \code{lm}/\code{glm}/\code{negbin} objects.
  //' Here, given as \code{NULL}.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}. Here, given as \code{NULL}.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model. Here, given as \code{NULL}.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}. Here, given as \code{NULL}.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to
  //' \code{1.0}.}
  //' \item{theta}{The estimated theta, if the response is negative binomial.
  //' Otherwise, will equal \code{1.0}.}
  //' 
  //' @keywords internal
  //' @noRd
  inline List numeric_extractor(const NumericVector& object) {
    String object_class = "numeric";
    String resp_family = "constant";
    int dist = 5;
    
    NumericVector coefs(1);
    coefs(0) = object(0);
    CharacterVector all_vars = {"Intercept"};
    
    double sigma {1.0};
    double theta {1.0};
    
    Rcpp::List output = List::create(_["class"] = object_class, _["family"] = resp_family,
      _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
      _["all_vars"] = all_vars, _["fixed_vars"] = all_vars, _["fixed_slopes"] = coefs,
      _["fixed_zi_vars"] = R_NilValue, _["fixed_zi_slopes"] = R_NilValue,
      _["random_vars"] = R_NilValue, _["random_slopes"] = R_NilValue,
      _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
      _["sigma"] = sigma, _["theta"] = theta);
  
    return output;
  }
  
  //' Extract Key Components of lm/glm/negbin Objects
  //' 
  //' This function extracts the components of an \code{lm}, \code{glm}, or
  //' \code{negbin} (function \code{glm.nb()}) object needed for functions
  //' \code{jerzeibalowski()} and \code{motherbalowski()}.
  //' 
  //' @name glm_extractor
  //' 
  //' @param object An \code{lm}, \code{glm}, or \code{negbin} object.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}. Will generally be either
  //' \code{lm} or \code{glm}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated. Not used in \code{lm}/\code{glm}/\code{negbin} objects.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated. Not used in \code{lm}/\code{glm}/\code{negbin} objects.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables. Not used in \code{lm}/\code{glm}/\code{negbin} objects.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
  //' objects.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model. Not used in \code{lm}/\code{glm}/
  //' \code{negbin} objects.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
  //' objects.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to
  //' \code{1.0}.}
  //' \item{theta}{The estimated theta, if the response is negative binomial.
  //' Otherwise, will equal \code{1.0}.}
  //' 
  //' @section Notes:
  //' Output from function \code{glm.nb()} is technically of class \code{negbin},
  //' but is treated as class \code{glm} here.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List glm_extractor(List object) {
    StringVector input_class = object.attr("class");
    std::string object_class = as<std::string>(input_class[0]);
    String resp_family = "gaussian";
    int dist = 2;
    double theta = 1.0;
    
    if (stringcompare_hard(object_class, "glm")) {
      List big_resp = object["family"];
      resp_family = as<String>(big_resp["family"]);
      
      List str_check_pois = stringcompare_soft(resp_family, "poisson");
      List str_check_negb = stringcompare_soft(resp_family, "negbin");
      List str_check_gamm = stringcompare_soft(resp_family, "gamma");
      List str_check_bin = stringcompare_soft(resp_family, "binomial");
      if (str_check_pois["contains"]) {
        dist = 0;
      } else if (str_check_negb["contains"]) {
        dist = 1;
        theta = object["theta"];
      } else if (str_check_gamm["contains"]) {
        dist = 3;
      } else if (str_check_bin["contains"]) {
        dist = 4;
      }
    } else if (stringcompare_hard(object_class, "negbin")) {
      resp_family = "negbin";
      dist = 1;
      theta = object["theta"];
    }
    
    NumericVector coefs = object["coefficients"];
    CharacterVector all_vars = coefs.attr("names");
    
    NumericVector residuals = object["residuals"];
    int no_residuals = residuals.length();
    double sum_squared_residuals {0.0};
    double sigma {1.0};
    
    for (int i = 0; i < no_residuals; i++) {
      sum_squared_residuals += (residuals(i) * residuals(i));
    }
    if (no_residuals > 2) {
      sigma = sum_squared_residuals / (static_cast<double>(no_residuals) - 2.0);
      sigma = sqrt(sigma);
    }
    
    Rcpp::List output = List::create(_["class"] = object_class, _["family"] = resp_family,
      _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
      _["all_vars"] = all_vars, _["fixed_vars"] = all_vars, _["fixed_slopes"] = coefs,
      _["fixed_zi_vars"] = R_NilValue, _["fixed_zi_slopes"] = R_NilValue,
      _["random_vars"] = R_NilValue, _["random_slopes"] = R_NilValue,
      _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
      _["sigma"] = sigma, _["theta"] = theta);
    
    return output;
  }
  
  //' Extract Key Components of vglm Objects
  //' 
  //' This function extracts the components of a \code{vglm} object needed for
  //' functions \code{jerzeibalowski()} and \code{motherbalowski()}.
  //' 
  //' @name vglm_extractor
  //' 
  //' @param object A \code{vglm} object.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated. Not used in \code{vglm} objects.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated. Always \code{TRUE} for \code{vglm} objects.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables. Not used in \code{vglm} objects.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}. Not used in \code{vglm} objects.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model. Not used in \code{vglm} objects.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}. Not used in \code{vglm} objects.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.}
  //' \item{theta}{The estimated theta, if the response is negative binomial.
  //' Otherwise, will equal \code{1.0}.}
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List vglm_extractor(S4 object) {
    int dist {0};
    double theta = object.slot("dispersion");
    
    S4 m_family = object.slot("family");
    String resp_family = m_family.slot("vfamily");
    
    List str_check_pois = stringcompare_soft(resp_family, "poisson");
    List str_check_negb = stringcompare_soft(resp_family, "negbin");
    if (str_check_pois["contains"]) {
      dist = 0;
    } else if (str_check_negb["contains"]) {
      dist = 1;
    } else {
      throw Rcpp::exception("Response distribution not recognized.", false);
    }
    
    NumericVector fixed_slopes = object.slot("coefficients");
    CharacterVector fixed_vars = fixed_slopes.attr("names");
    
    List all_terms = object.slot("terms");
    List all_terms_stuff = all_terms["terms"];
    NumericMatrix ats = all_terms_stuff.attr("factors");
    CharacterVector all_vars = rownames(ats);
    
    Rcpp::List output = List::create(_["class"] = "vglm", _["family"] = resp_family,
      _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = true,
      _["all_vars"] = all_vars, _["fixed_vars"] = fixed_vars,
      _["fixed_slopes"] = fixed_slopes, _["fixed_zi_vars"] = R_NilValue,
      _["fixed_zi_slopes"] = R_NilValue, _["random_vars"] = R_NilValue,
      _["random_slopes"] = R_NilValue, _["random_zi_vars"] = R_NilValue,
      _["random_zi_slopes"] = R_NilValue, _["sigma"] = 1.0, _["theta"] = theta);
    
    return output;
  }
  
  //' Extract Key Components of zeroinfl Objects
  //' 
  //' This function extracts the components of a \code{zeroinfl} object needed for
  //' functions \code{jerzeibalowski()} and \code{motherbalowski()} to work.
  //' 
  //' @name zeroinfl_extractor
  //' 
  //' @param object A \code{zeroinfl} object.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}. Will always be
  //' \code{glmmTMB}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated. Defaults to \code{FALSE} for \code{zeroinfl} objects.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}.}
  //' \item{random_vars}{A string vector holding the names of the random
  //' variables. Not used in \code{zeroinfl} objects.}
  //' \item{random_slopes}{A numeric vector holding the slope coefficients of the
  //' random variables, in the same order as \code{random_vars}. Not used in
  //' \code{zeroinfl} objects.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model. Not used in \code{zeroinfl} objects.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}. Not used in \code{zeroinfl} objects.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
  //' Equivalent output to lme4's \code{sigma()} function.}
  //' \item{theta}{The scale parameter theta used in the negative binomial
  //' distribution. Defaults to \code{1.0}.}
  //' 
  //' @section Notes:
  //' This function will only work in the case where random terms are given as
  //' \code{(1 | ranterm)}, where \code{ranterm} is the name of the random
  //' variable.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List zeroinfl_extractor(List object) {
    String model_family = object["dist"];
    int dist {0};
    double theta {1.0};
    bool zi = true;
    
    List str_check_pois = stringcompare_soft(model_family, "poisson");
    List str_check_negb = stringcompare_soft(model_family, "negbin");
    
    if (str_check_pois["contains"]) {
      dist = 0;
      
    } else if (str_check_negb["contains"]) {
      dist = 1;
      theta = object["theta"];
      
    } else {
      throw Rcpp::exception("Unrecognized response distribution.", false);
    }
    
    List model_part = object["model"];
    CharacterVector all_vars = model_part.attr("names");
    
    List all_coefs = object["coefficients"];
    NumericVector fixed_slopes = all_coefs["count"];
    NumericVector zi_slopes = all_coefs["zero"];
    CharacterVector fixed_terms = fixed_slopes.attr("names");
    CharacterVector zi_terms = zi_slopes.attr("names");
    
    List output = List::create(_["class"] = "zeroinfl", _["family"] = model_family,
      _["dist"] = dist, _["zero_inflated"] = zi, _["zero_truncated"] = false,
      _["all_vars"] = all_vars, _["fixed_vars"] = fixed_terms,
      _["fixed_slopes"] = fixed_slopes, _["fixed_zi_vars"] = zi_terms,
      _["fixed_zi_slopes"] = zi_slopes, _["random_vars"] = R_NilValue,
      _["random_slopes"] = R_NilValue, _["random_zi_vars"] = R_NilValue,
      _["random_zi_slopes"] = R_NilValue, _["sigma"] = 1.0, _["theta"] = theta);
    
    return output;
  }
  
  //' Extract Key Components of merMod Objects
  //' 
  //' This function extracts the components of a \code{merMod} object needed for
  //' functions \code{jerzeibalowski()} and \code{motherbalowski()} to work.
  //' 
  //' @name lme4_extractor
  //' 
  //' @param object A \code{merMod} object.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}. Will generally be either
  //' \code{lmerMod} or \code{glmerMod}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated. Not used in lme4.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated. Not used in lme4.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables. Not used in lme4 objects.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}. Not used in lme4 objects.}
  //' \item{random_vars}{A string vector holding the names of the random
  //' variables.}
  //' \item{random_slopes}{A numeric vector holding the slope coefficients of the
  //' random variables, in the same order as \code{random_vars}.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model. Not used in lme4.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}. Not used in lme4.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
  //' Equivalent output to lme4's \code{sigma()} function.}
  //' \item{theta}{Not used in lme4 output. Defaults to \code{1.0}.}
  //' 
  //' @section Notes:
  //' This function will only work in the case where random terms are given as
  //' \code{(1 | ranterm)}, where \code{ranterm} is the name of the random
  //' variable.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List lme4_extractor(S4 object) {
    std::string object_class = object.attr("class");
    NumericVector coefs = object.slot("beta");
    String resp_family = "gaussian";
    int dist = 2;
    List pp = object.slot("pp");
    List cnms = object.slot("cnms");
    StringVector ran_names = cnms.attr("names");
    List flist = object.slot("flist");
    int no_ran_terms = flist.length();
    double sigma {1.0};
    double theta {1.0};
    
    if (stringcompare_hard(object_class, "glmerMod")) {
      List big_resp = object.slot("resp");
      List resp_family_S3 = big_resp["family"];
      resp_family = as<String>(resp_family_S3["family"]);
      
      List str_check_pois = stringcompare_soft(resp_family, "poisson");
      List str_check_negb = stringcompare_soft(resp_family, "negbin");
      List str_check_gamm = stringcompare_soft(resp_family, "gamma");
      List str_check_bin = stringcompare_soft(resp_family, "binomial");
      if (str_check_pois["contains"]) {
        dist = 0;
      } else if (str_check_negb["contains"]) {
        dist = 1;
      } else if (str_check_gamm["contains"]) {
        dist = 3;
      } else if (str_check_bin["contains"]) {
        dist = 4;
      }
    }
    
    // Names of all variables
    DataFrame oframe = as<DataFrame>(object.slot("frame"));
    StringVector all_var_names = oframe.attr("names");
    
    // Fixed variables
    NumericMatrix ppX = pp["X"];
    StringVector coef_names = colnames(ppX);
    
    // Random variables
    arma::sp_mat lambdat_sp = as<sp_mat>(pp["Lambdat"]);
    arma::mat lambdat = arma::mat(lambdat_sp);
    NumericVector u = object.slot("u");
    arma::vec ucolvec = arma::vec(u);
    
    // b gives the random terms when all random terms are (1 | ranterm)
    arma::mat b = lambdat * ucolvec;
    NumericVector ran_slopes(b.begin(), b.end());
    
    CharacterVector ran_term_index(ran_slopes.length());
    int ran_term_index_counter {0};
    List ran_term_list(no_ran_terms);
    List ran_index_list(no_ran_terms);
    
    // Names of random variables
    for (int i = 0; i < no_ran_terms; i++) {
      IntegerVector new_term = flist[i];
      CharacterVector new_term_names = new_term.attr("levels");
      
      int new_term_names_length = new_term_names.length();
      NumericVector new_value(new_term_names_length);
      
      for (int j = 0; j < new_term_names_length; j++) {
        new_value(j) = ran_slopes(ran_term_index_counter);
        ran_term_index_counter++;
      }
      new_value.attr("names") = new_term_names;
      ran_term_list(i) = new_value;
      ran_index_list(i) = new_value.attr("names");
    }
    ran_index_list.attr("names") = ran_names;
    ran_term_list.attr("names") = ran_names;
    
    List devcomp = object.slot("devcomp");
    NumericVector cmp = devcomp["cmp"];
    IntegerVector dd = devcomp["dims"];
    CharacterVector cmp_names = cmp.attr("names");
    CharacterVector dd_names = dd.attr("names");
    int cmp_length = cmp.length();
    int dd_length = dd.length();
    int useSc_place = 0;
    int REML_place = 0;
    bool useSc = false;
    bool REML = false;
    
    for (int i = 0; i < dd_length; i++) {
      if (stringcompare_hard(as<std::string>(dd_names(i)), "useSc")) {
        useSc_place = i;
      } else if (stringcompare_hard(as<std::string>(dd_names(i)), "REML")) {
        REML_place = i;
      }
    }
    if (dd(useSc_place) > 0) {
      useSc = true;
      
      if (dd(REML_place) > 0) {
        REML = true;
      }
    }
    
    if (useSc) {
      for (int i = 0; i < cmp_length; i++) {
        if (REML) {
          if (stringcompare_hard(as<std::string>(cmp_names(i)), "sigmaREML")) {
            sigma = cmp(i);
          }
        } else {
          if (stringcompare_hard(as<std::string>(cmp_names(i)), "sigmaML")) {
            sigma = cmp(i);
          }
        }
      }
    }
    
    Rcpp::List output = List::create(_["class"] = object_class, _["family"] = resp_family,
      _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
      _["all_vars"] = all_var_names, _["fixed_vars"] = coef_names,
      _["fixed_slopes"] = coefs, _["fixed_zi_vars"] = R_NilValue,
      _["fixed_zi_slopes"] = R_NilValue, _["random_vars"] = ran_index_list,
      _["random_slopes"] = ran_term_list, _["random_zi_vars"] = R_NilValue,
      _["random_zi_slopes"] = R_NilValue, _["sigma"] = sigma, _["theta"] = theta);
    
    return output;
  }
  
  //' Extract Key Components of glmmTMB Objects
  //' 
  //' This function extracts the components of a \code{glmmTMB} object needed for
  //' functions \code{jerzeibalowski()} and \code{motherbalowski()} to work.
  //' 
  //' @name glmmTMB_extractor
  //' 
  //' @param object A \code{glmmTMB} object.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}. Will always be
  //' \code{glmmTMB}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}.}
  //' \item{random_vars}{A string vector holding the names of the random
  //' variables.}
  //' \item{random_slopes}{A numeric vector holding the slope coefficients of the
  //' random variables, in the same order as \code{random_vars}.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
  //' Equivalent output to lme4's \code{sigma()} function.}
  //' \item{theta}{The scale parameter theta used in the negative binomial
  //' distribution. Defaults to \code{1.0}.}
  //' 
  //' @section Notes:
  //' This function will only work in the case where random terms are given as
  //' \code{(1 | ranterm)}, where \code{ranterm} is the name of the random
  //' variable.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List glmmTMB_extractor(List object) {
    
    List the_green_flist;
    List the_green_ziflist;
    CharacterVector ran_zi_vars;
    
    List dissident_aggressor = object["fit"];
    NumericVector parameter_values = dissident_aggressor["par"];
    CharacterVector term_names = parameter_values.attr("names");
    
    NumericVector all_the_stuff = dissident_aggressor["parfull"];
    CharacterVector names_of_all_the_stuff = all_the_stuff.attr("names");
    int length_of_all_the_stuff = all_the_stuff.length();
    
    String class_thistime = {"glmmTMB"};
    int dist = 2;
    double sigma {1.0};
    double theta {1.0};
    
    List variance_crap = object["sdr"];
    NumericVector random_values = variance_crap["par.random"];
    CharacterVector random_values_tags = random_values.attr("names");
    int no_random_values = random_values.length(); // Has ALL random coefficients
    LogicalVector random_b(no_random_values);
    LogicalVector random_bzi(no_random_values);
    int no_ranb {0};
    int no_ranbzi {0};
    int no_ran_vars {0};
    int no_rz_vars {0};  
    
    for (int i = 0; i < no_random_values; i++) {
      if (stringcompare_hard(as<std::string>(random_values_tags(i)), "b")) {
        no_ranb++;
        random_b(i) = 1;
      } else if (stringcompare_hard(as<std::string>(random_values_tags(i)), "bzi")) {
        no_ranbzi++;
        random_bzi(i) = 1;
      }
    }
    
    DataFrame exploding_reloading = as<DataFrame>(object["frame"]);
    CharacterVector all_vars = exploding_reloading.attr("names");
    
    List the_green_manalishi = object["modelInfo"];
    CharacterVector ran_vars = the_green_manalishi["grpVar"];
    List the_green_reTrms = the_green_manalishi["reTrms"];
    List the_green_cond = the_green_reTrms["cond"];
    List the_green_zi = the_green_reTrms["zi"];
    
    if (no_ranb > 0) {
      the_green_flist = the_green_cond["flist"];
      no_ran_vars = the_green_flist.length();
    }
    if (no_ranbzi > 0) {
      the_green_ziflist = the_green_zi["flist"];
      no_rz_vars = the_green_ziflist.length();
      ran_zi_vars = the_green_ziflist.attr("names");
    }
    
    List obj = object["obj"];
    Environment env = obj["env"];
    List data = env["data"];
    NumericMatrix Xmat = data["X"];
    CharacterVector all_fixed_terms = colnames(Xmat);
    NumericMatrix Xmatzi = data["Xzi"];
    CharacterVector all_zi_terms = colnames(Xmatzi);
    
    int total_fixed_slopes {0};
    int total_zi_slopes {0};
    
    // Extract random slopes
    List ran_term_list(no_ran_vars);
    List ran_slope_list(no_ran_vars);
    List ran_zi_term_list;
    List ran_zi_slope_list;
    int ran_term_counter {0};
    
    // Extract cond random coefficients
    if (no_ranb > 0) {
      for (int i = 0; i < no_ran_vars; i++) {
        
        IntegerVector current_factor = the_green_flist[i];
        CharacterVector current_names = current_factor.attr("levels");
        
        ran_term_list(i) = current_names;
        int current_names_length = current_names.length();
        NumericVector current_slopes(current_names_length);
        
        for (int j = 0; j < current_names_length; j++) {
          if (random_b(ran_term_counter) > 0) {
            current_slopes(j) = random_values(ran_term_counter);
          }
          ran_term_counter++;
        }
        ran_slope_list(i) = current_slopes;
      }
      ran_term_list.attr("names") = ran_vars;
      ran_slope_list.attr("names") = ran_vars;
    }
    
    // Extract zi random coefficients
    if (no_ranbzi > 0) {
      List rztl(no_rz_vars);
      List rstl(no_rz_vars);
      
      for (int i = 0; i < no_rz_vars; i++) {
        IntegerVector current_factor = the_green_ziflist[i];
        CharacterVector current_names = current_factor.attr("levels");
        
        rztl(i) = current_names;
        int current_names_length = current_names.length();
        NumericVector current_slopes(current_names_length);
        
        for (int j = 0; j < current_names_length; j++) {
          if (random_bzi(ran_term_counter) > 0) {
            current_slopes(j) = random_values(ran_term_counter);
          }
          ran_term_counter++;
        }
        rstl(i) = current_slopes;
      }
      
      ran_zi_term_list = rztl;
      ran_zi_slope_list = rstl;
      ran_zi_term_list.attr("names") = ran_zi_vars;
      ran_zi_slope_list.attr("names") = ran_zi_vars;
    }
    
    // Extract dispersion parameter
    double dispersion {1.0};
    for (int i = 0; i < length_of_all_the_stuff; i++) {
      if (stringcompare_hard(as<std::string>(names_of_all_the_stuff(i)), "betad")) {
        dispersion = all_the_stuff(i);
      }
    }
    
    // Response distribution tests
    bool trunc = false;
    bool zi = false;
    
    List with_the_two_pronged_crown = the_green_manalishi["family"];
    String model_family = as<String>(with_the_two_pronged_crown["family"]);
    
    List str_check_pois = stringcompare_soft(model_family, "poisson");
    List str_check_negb = stringcompare_soft(model_family, "negbin");
    List str_check_nbin = stringcompare_soft(model_family, "nbinom");
    List str_check_gamm = stringcompare_soft(model_family, "Gamma");
    List str_check_bin = stringcompare_soft(model_family, "binomial");
    List str_check_trun = stringcompare_soft(model_family, "trunc");
    
    if (str_check_pois["contains"]) {
      dist = 0;
      
      if (str_check_trun["contains"]) {
        trunc = true;
      }
    } else if (str_check_negb["contains"] || str_check_nbin["contains"]) {
      dist = 1;
      theta = exp(dispersion);
      
      if (str_check_trun["contains"]) {
        trunc = true;
      }
    } else if (str_check_gamm["contains"]) {
      dist = 3;
      sigma = exp(-0.5 * dispersion);
      
      List str_check_zi = stringcompare_soft(model_family, "zi");
      if (str_check_zi["contains"]) {
        zi = true;
      }
    } else if (str_check_bin["contains"]) {
      dist = 4;
    } else if (dist == 2) {
      sigma = exp(0.5 * dispersion);
    } else {
      throw Rcpp::exception("Unrecognized response distribution.", false);
    }
    
    // Most zero-inflation tests (apart from ziGamma)
    int no_params = parameter_values.length();
    
    for (int i = 0; i < no_params; i++) {
      if (term_names(i) == "betazi" && zi == false) {
        //model_family += zi_addition;
        zi = true;
      }
      
      if (term_names(i) == "beta") {
        
        total_fixed_slopes++;
      } else if (term_names(i) == "betazi") {
        total_zi_slopes++;
      }
    }
    
    // Extract fixed and zero-inflated slopes
    NumericVector fixed_slopes(total_fixed_slopes);
    NumericVector zi_slopes(total_zi_slopes);
    CharacterVector fixed_terms(total_fixed_slopes);
    CharacterVector zi_terms(total_zi_slopes);
    int fixed_slope_counter {0};
    int zi_slope_counter {0};
    
    for (int i = 0; i < no_params; i++) {
      if (term_names(i) == "beta") {
        fixed_slopes(fixed_slope_counter) = parameter_values(i);
        fixed_terms(fixed_slope_counter) = all_fixed_terms(fixed_slope_counter);
        fixed_slope_counter++;
      } else if (term_names(i) == "betazi") {
        zi_slopes(zi_slope_counter) = parameter_values(i);
        zi_terms(zi_slope_counter) = all_zi_terms(zi_slope_counter);
        zi_slope_counter++;
      }
    }
    
    List output = List::create(_["class"] = class_thistime, _["family"] = model_family,
      _["dist"] = dist, _["zero_inflated"] = zi, _["zero_truncated"] = trunc,
      _["all_vars"] = all_vars, _["fixed_vars"] = fixed_terms,
      _["fixed_slopes"] = fixed_slopes, _["fixed_zi_vars"] = zi_terms,
      _["fixed_zi_slopes"] = zi_slopes, _["random_vars"] = ran_term_list,
      _["random_slopes"] = ran_slope_list, _["random_zi_vars"] = ran_zi_term_list,
      _["random_zi_slopes"] = ran_zi_slope_list, _["sigma"] = sigma,
      _["theta"] = theta);
    
    return output;
  }
  
  //' Extract Core Components From S3 Vital Rate Models
  //' 
  //' Function \code{S3_extractor()} extracts all needed terms from S3 objects
  //' used as vital rate models.
  //' 
  //' @name S3_extractor
  //' 
  //' @param object An S3 vital rate model other than objects of class
  //' \code{vrm_input}. Currently, this should be output from functions
  //' \code{lm()}, \code{glm()}, \code{glm.nb()}, \code{zeroinfl()}, and
  //' \code{glmmTMB()}.
  //' 
  //' @return A list describing the vital rate model in standard output required
  //' from function \code{modelextract()} to operate. Elements currently include:
  //' \item{class}{The exact class of \code{object}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}.}
  //' \item{random_vars}{A string vector holding the names of the random
  //' variables.}
  //' \item{random_slopes}{A numeric vector holding the slope coefficients of the
  //' random variables, in the same order as \code{random_vars}.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
  //' Equivalent output to lme4's \code{sigma()} function.}
  //' \item{theta}{The scale parameter theta used in the negative binomial
  //' distribution. Defaults to \code{1.0}.}
  //' 
  //' @section Notes:
  //' This function currently handles models developed with functions \code{lm()}
  //' and \code{glm()} from package \code{stats}, function \code{glm.nb()} from
  //' package \code{MASS}, function \code{zeroinfl()} from package \code{pscl},
  //' and function \code{glmmTMB()} from package \code{glmmTMB}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List S3_extractor(List object) {
    StringVector model_class = object.attr("class");
    int model_type {0}; // 0 = unknown, 1 = lm/glm/negbin, 2 = zeroinfl, 3 = glmmTMB
    
    List output;
    
    for (int i = 0; i < model_class.length(); i++) {
      if (stringcompare_hard(as<std::string>(model_class(i)), "lm")) {
        model_type = 1;
      } else if (stringcompare_hard(as<std::string>(model_class(i)), "zeroinfl")) {
        model_type = 2;
      } else if (stringcompare_hard(as<std::string>(model_class(i)), "glmmTMB")) {
        model_type = 3;
      }
    }
    
    if (model_type == 1) {
      output = glm_extractor(object);
    } else if (model_type == 2) {
      output = zeroinfl_extractor(object);
    } else if (model_type == 3) {
      output = glmmTMB_extractor(object);
    } else {
      throw Rcpp::exception("Model type unrecognized.", false);
    }
    
    return output;
  }
  
  //' Extract Core Components From S4 Vital Rate Models
  //' 
  //' Function \code{S4_extractor()} extracts all needed terms from S4 objects
  //' used as vital rate models.
  //' 
  //' @name S4_extractor
  //' 
  //' @param object An S4 vital rate model. Currently, this should be output from
  //' functions \code{vglm()}, \code{lmer()}, and \code{glmer()}.
  //' 
  //' @return A list describing the vital rate model in standard output required
  //' from function \code{modelextract()} to operate. Elements currently include:
  //' \item{class}{The exact class of \code{object}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}.}
  //' \item{random_vars}{A string vector holding the names of the random
  //' variables.}
  //' \item{random_slopes}{A numeric vector holding the slope coefficients of the
  //' random variables, in the same order as \code{random_vars}.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
  //' Equivalent output to lme4's \code{sigma()} function.}
  //' \item{theta}{The scale parameter theta used in the negative binomial
  //' distribution. Defaults to \code{1.0}.}
  //' 
  //' @section Notes:
  //' This function currently handles models developed with function \code{vglm()}
  //' from package \code{VGAM}, and functions \code{lmer()} and \code{glmer()}
  //' from package \code{lme4}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List S4_extractor(S4 object) {
    StringVector model_class = wrap(object.attr("class"));
    //String model_class = object.attr("class");
    
    List output;
    bool found {false};
    int mc_length = static_cast<int>(model_class.length());
    int counter {0};
    
    while (!found && counter < mc_length) {
      if (stringcompare_hard(String(model_class(counter)), "vglm")) {
        found = true;
        output = vglm_extractor(object);
        
      } else if (stringcompare_hard(String(model_class(counter)), "lmerMod") || 
          stringcompare_hard(String(model_class(counter)), "glmerMod")) {
        found = true;
        output = lme4_extractor(object);
        
      }
      counter++;
    }
    
    if (!found) {
      throw Rcpp::exception("Model type unrecognized.", false);
    }
    
    return output;
  }
  
  //' Extract Key Components of vrm_input Objects
  //' 
  //' This function extracts the components of a \code{vrm_input} object for
  //' functions \code{jerzeibalowski()} and \code{motherbalowski()}.
  //' 
  //' @name vrm_extractor
  //' 
  //' @param object A \code{vrm_input} object.
  //' 
  //' @return A list with the following elements:
  //' \item{class}{The exact class of \code{object}. Will be \code{vrm_input}.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated. Not used in \code{lm}/\code{glm}/\code{negbin} objects. Will
  //' equal \code{FALSE} if \code{dist = 5}.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated. Not used in \code{lm}/\code{glm}/\code{negbin} objects. Will
  //' equal \code{FALSE} if \code{dist = 5}.}
  //' \item{all_vars}{A vector holding the names of each variable used by
  //' \code{object}. Will equal only the intercept if \code{dist = 5}.}
  //' \item{fixed_vars}{A string vector holding the names of the fixed variables.
  //' Will equal only the intercept if \code{dist = 5}.}
  //' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
  //' fixed variables, in the same order as \code{fixed_vars}. Will equal only the
  //' intercept if \code{dist = 5}.}
  //' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
  //' fixed variables. Not used in \code{lm}/\code{glm}/\code{negbin} objects.
  //' Will equal \code{NULL} is \code{dist = 5}.}
  //' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the zero-inflated fixed variables, in the same order as
  //' \code{fixed_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
  //' objects. Will equal \code{NULL} is \code{dist = 5}.}
  //' \item{random_zi_vars}{A string vector holding the names of the random
  //' variables in the zero-inflation model. Not used in \code{lm}/\code{glm}/
  //' \code{negbin} objects. Will equal \code{NULL} is \code{dist = 5}.}
  //' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
  //' the random variables in the zero-inflation model, in the same order as
  //' \code{random_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
  //' objects. Will equal \code{NULL} is \code{dist = 5}.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to
  //' \code{1.0}.}
  //' \item{theta}{The estimated theta, if the response is negative binomial.
  //' Otherwise, will equal \code{1.0}.}
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List vrm_extractor(List object) {
    std::string object_class = "vrm_input";
    
    String distrib = as<String>(object["dist"]);
    double sigma_theta = object["sigma_theta"];
    
    String resp_family = "gaussian";
    int dist = 2;
    double theta {1.0};
    double sigma = sigma_theta;
    bool zero_inflation = false;
    bool zero_truncation = false;
    
    if (stringcompare_simple(distrib, "poisson", false)) {
      resp_family = "poisson";
      dist = 0;
    } else if (stringcompare_simple(distrib, "negbin", false)) {
      resp_family = "negbin";
      dist = 1;
      theta = sigma_theta;
    } else if (stringcompare_simple(distrib, "nbinom2", false)) {
      resp_family = "negbin";
      dist = 1;
      theta = sigma_theta;
    } else if (stringcompare_simple(distrib, "gamma", false)) {
      resp_family = "gamma";
      dist = 3;
    } else if (stringcompare_simple(distrib, "binom", false)) {
      resp_family = "binomial";
      dist = 4;
    } else if (stringcompare_simple(distrib, "cons", false)) {
      resp_family = "constant";
      dist = 5;
    }
    
    if (stringcompare_simple(distrib, "trunc", false)) {
      if (dist == 0 || dist == 1) zero_truncation = true; 
    }
    
    NumericVector coefs = as<NumericVector>(object["fixed_slopes"]);
    NumericVector zi_coefs = as<NumericVector>(object["fixed_zi"]);
    
    CharacterVector all_vars = as<CharacterVector>(object("effects_names"));
    
    NumericVector year_slopes = as<NumericVector>(object["year_slopes"]);
    NumericVector patch_slopes = as<NumericVector>(object["patch_slopes"]);
    NumericVector indcova2_slopes = as<NumericVector>(object["indcova2_slopes"]);
    NumericVector indcova1_slopes = as<NumericVector>(object["indcova1_slopes"]);
    NumericVector indcovb2_slopes = as<NumericVector>(object["indcovb2_slopes"]);
    NumericVector indcovb1_slopes = as<NumericVector>(object["indcovb1_slopes"]);
    NumericVector indcovc2_slopes = as<NumericVector>(object["indcovc2_slopes"]);
    NumericVector indcovc1_slopes = as<NumericVector>(object["indcovc1_slopes"]);
    
    NumericVector year_zi = as<NumericVector>(object["year_zi"]);
    NumericVector patch_zi = as<NumericVector>(object["patch_zi"]);
    NumericVector indcova2_zi = as<NumericVector>(object["indcova2_zi"]);
    NumericVector indcova1_zi = as<NumericVector>(object["indcova1_zi"]);
    NumericVector indcovb2_zi = as<NumericVector>(object["indcovb2_zi"]);
    NumericVector indcovb1_zi = as<NumericVector>(object["indcovb1_zi"]);
    NumericVector indcovc2_zi = as<NumericVector>(object["indcovc2_zi"]);
    NumericVector indcovc1_zi = as<NumericVector>(object["indcovc1_zi"]);
    
    double zi_year_sum = sum(year_zi);
    double zi_patch_sum = sum(patch_zi);
    double zi_indcova2_sum = sum(indcova2_zi);
    double zi_indcova1_sum = sum(indcova1_zi);
    double zi_indcovb2_sum = sum(indcovb2_zi);
    double zi_indcovb1_sum = sum(indcovb1_zi);
    double zi_indcovc2_sum = sum(indcovc2_zi);
    double zi_indcovc1_sum = sum(indcovc1_zi);
    
    double zi_sum = zi_year_sum + zi_patch_sum + zi_indcova2_sum + zi_indcova1_sum +
      zi_indcovb2_sum + zi_indcovb1_sum + zi_indcovc2_sum + zi_indcovc1_sum;
    
    if (zi_sum != 0.0) {
      zero_inflation = true;
    }
    
    if (zero_truncation && zero_inflation) {
      throw Rcpp::exception("Models cannot be both zero-inflated and zero-truncated.",
        false);
    }
    
    CharacterVector year_names = as<CharacterVector>(object["year_names"]);
    CharacterVector patch_names = as<CharacterVector>(object["patch_names"]);
    CharacterVector indcova_names = as<CharacterVector>(object["indcova_names"]);
    CharacterVector indcovb_names = as<CharacterVector>(object["indcovb_names"]);
    CharacterVector indcovc_names = as<CharacterVector>(object["indcovc_names"]);
    
    year_slopes.attr("names") = year_names;
    patch_slopes.attr("names") = patch_names;
    indcova2_slopes.attr("names") = indcova_names;
    indcova1_slopes.attr("names") = indcova_names;
    indcovb2_slopes.attr("names") = indcovb_names;
    indcovb1_slopes.attr("names") = indcovb_names;
    indcovc2_slopes.attr("names") = indcovc_names;
    indcovc1_slopes.attr("names") = indcovc_names;
    
    CharacterVector randomvar_names = {"year2", "patch", "indcova2", "indcova1",
      "indcovb2", "indcovb1", "indcovc2", "indcovc1"};
    
    List random_slopes (8);
    random_slopes(0) = year_slopes;
    random_slopes(1) = patch_slopes;
    random_slopes(2) = indcova2_slopes;
    random_slopes(3) = indcova1_slopes;
    random_slopes(4) = indcovb2_slopes;
    random_slopes(5) = indcovb1_slopes;
    random_slopes(6) = indcovc2_slopes;
    random_slopes(7) = indcovc1_slopes;
    
    List random_zi (8);
    random_zi(0) = year_zi;
    random_zi(1) = patch_zi;
    random_zi(2) = indcova2_zi;
    random_zi(3) = indcova1_zi;
    random_zi(4) = indcovb2_zi;
    random_zi(5) = indcovb1_zi;
    random_zi(6) = indcovc2_zi;
    random_zi(7) = indcovc1_zi;
    
    List random_names (8);
    random_names(0) = year_names;
    random_names(1) = patch_names;
    random_names(2) = indcova_names;
    random_names(3) = indcova_names;
    random_names(4) = indcovb_names;
    random_names(5) = indcovb_names;
    random_names(6) = indcovc_names;
    random_names(7) = indcovc_names;
    
    random_slopes.attr("names") = randomvar_names;
    random_zi.attr("names") = randomvar_names;
    random_names.attr("names") = randomvar_names;
    
    NumericVector group2_slopes = as<NumericVector>(object["group2_slopes"]);
    NumericVector group1_slopes = as<NumericVector>(object["group1_slopes"]);
    NumericVector group2_zi = as<NumericVector>(object["group2_zi"]);
    NumericVector group1_zi = as<NumericVector>(object["group1_zi"]);
    CharacterVector group2_names = clone(as<CharacterVector>(object["group_names"]));
    CharacterVector group1_names = clone(as<CharacterVector>(object["group_names"]));
    CharacterVector group2_zi_names = clone(as<CharacterVector>(object["group_names"]));
    CharacterVector group1_zi_names = clone(as<CharacterVector>(object["group_names"]));
    
    CharacterVector zi_vars = clone(all_vars);
    
    int no_group2_slopes = group2_slopes.length();
    int no_group1_slopes = group1_slopes.length();
    
    if (no_group2_slopes > 0) {
      for (int i = 0; i < no_group2_slopes; i++) {
        group2_names[i] = "group2" + group2_names[i];
      }
      all_vars = concat_str(all_vars, group2_names);
      coefs = concat_dbl(coefs, group2_slopes);
    }
    
    if (no_group1_slopes > 0) {
      for (int i = 0; i < no_group1_slopes; i++) {
        group1_names[i] = "group1" + group1_names[i];
      }
      all_vars = concat_str(all_vars, group1_names);
      coefs = concat_dbl(coefs, group1_slopes);
    }
    
    int no_group2_zi = group2_zi.length();
    int no_group1_zi = group1_zi.length();
    
    if (no_group2_zi > 0) {
      for (int i = 0; i < no_group2_zi; i++) {
        group1_zi_names[i] = "group2" + group2_zi_names[i];
      }
      zi_vars = concat_str(zi_vars, group2_zi_names);
      zi_coefs = concat_dbl(zi_coefs, group2_zi);
    }
    
    if (no_group1_zi > 0) {
      for (int i = 0; i < no_group1_zi; i++) {
        group1_zi_names[i] = "group1" + group1_zi_names[i];
      }
      zi_vars = concat_str(zi_vars, group1_zi_names);
      zi_coefs = concat_dbl(zi_coefs, group1_zi);
    }
    
    Rcpp::List output;
    if (dist == 5) {
      CharacterVector simple_vars = {"Intercept"};
      NumericVector simple_coefs (1);
      simple_coefs(0) = coefs(0);
      
      output = List::create(_["class"] = object_class, _["family"] = resp_family,
        _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
        _["all_vars"] = simple_vars, _["fixed_vars"] = simple_vars,
        _["fixed_slopes"] = simple_coefs, _["fixed_zi_vars"] = R_NilValue,
        _["fixed_zi_slopes"] = R_NilValue, _["random_vars"] = R_NilValue,
        _["random_slopes"] = R_NilValue, _["random_zi_vars"] = R_NilValue,
        _["random_zi_slopes"] = R_NilValue, _["sigma"] = sigma, _["theta"] = theta);
      
    } else  {
      output = List::create(_["class"] = object_class, _["family"] = resp_family,
        _["dist"] = dist, _["zero_inflated"] = zero_inflation,
        _["zero_truncated"] = zero_truncation, _["all_vars"] = all_vars,
        _["fixed_vars"] = all_vars, _["fixed_slopes"] = coefs,
        _["fixed_zi_vars"] = zi_vars, _["fixed_zi_slopes"] = zi_coefs,
        _["random_vars"] = random_names, _["random_slopes"] = random_slopes,
        _["random_zi_vars"] = random_names, _["random_zi_slopes"] = random_zi,
        _["sigma"] = sigma, _["theta"] = theta);
    }
    
    return output;
  }
  
  //' Create Matrices of Year and Patch Terms in Models
  //' 
  //' Function \code{revelations()} creates a matrix holding either the year or
  //' patch coefficients from all vital rate models. This reduces memory load in
  //' functions \code{\link{jerzeibalowski}()}, which may be important in some
  //' systems or compilers.
  //' 
  //' @name revelations
  //' 
  //' @param survproxy The proxy vital rate model covering survival from the main
  //' matrix estimator function.
  //' @param obsproxy The proxy vital rate model covering observation status from
  //' the main matrix estimator function.
  //' @param sizeproxy The proxy vital rate model covering primary size from the
  //' main matrix estimator function.
  //' @param sizebproxy The proxy vital rate model covering secondary size from
  //' the main matrix estimator function.
  //' @param sizecproxy The proxy vital rate model covering tertiary size from the
  //' main matrix estimator function.
  //' @param repstproxy The proxy vital rate model covering reproductive status
  //' from the main matrix estimator function.
  //' @param fecproxy The proxy vital rate model covering fecundity from the main
  //' matrix estimator function.
  //' @param jsurvproxy The proxy vital rate model covering juvenile survival from
  //' the main matrix estimator function.
  //' @param jobsproxy The proxy vital rate model covering juvenile observation
  //' status from the main matrix estimator function.
  //' @param jsizeproxy The proxy vital rate model covering juvenile primary size
  //' from the main matrix estimator function.
  //' @param jsizebproxy The proxy vital rate model covering juvenile secondary
  //' size from the main matrix estimator function.
  //' @param jsizecproxy The proxy vital rate model covering juvenile tertiary
  //' size from the main matrix estimator function.
  //' @param jrepstproxy The proxy vital rate model covering juvenile reproductive
  //' status from the main matrix estimator function.
  //' @param jmatstproxy The proxy vital rate model covering juvenile probability
  //' of becoming mature from the main matrix estimator function.
  //' @param mat_switch An integer coding for year (\code{1}) or patch (\code{2}).
  //' 
  //' @return A matrix with 14 columns corresponding to the number of vital rates
  //' and number of columns equal to the number of year or patches.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::NumericMatrix revelations(List survproxy, List obsproxy, List sizeproxy,
    List sizebproxy, List sizecproxy, List repstproxy, List fecproxy,
    List jsurvproxy, List jobsproxy, List jsizeproxy, List jsizebproxy,
    List jsizecproxy, List jrepstproxy,List jmatstproxy, int mat_switch) {
    
    NumericMatrix final_mat;
    
    if (mat_switch == 1) {
      NumericVector survyear = as<NumericVector>(survproxy["years"]);
      NumericVector obsyear = as<NumericVector>(obsproxy["years"]);
      NumericVector sizeyear = as<NumericVector>(sizeproxy["years"]);
      NumericVector sizebyear = as<NumericVector>(sizebproxy["years"]);
      NumericVector sizecyear = as<NumericVector>(sizecproxy["years"]);
      NumericVector repstyear = as<NumericVector>(repstproxy["years"]);
      NumericVector fecyear = as<NumericVector>(fecproxy["years"]);
      NumericVector jsurvyear = as<NumericVector>(jsurvproxy["years"]);
      NumericVector jobsyear = as<NumericVector>(jobsproxy["years"]);
      NumericVector jsizeyear = as<NumericVector>(jsizeproxy["years"]);
      NumericVector jsizebyear = as<NumericVector>(jsizebproxy["years"]);
      NumericVector jsizecyear = as<NumericVector>(jsizecproxy["years"]);
      NumericVector jrepstyear = as<NumericVector>(jrepstproxy["years"]);
      NumericVector jmatstyear = as<NumericVector>(jmatstproxy["years"]);
      
      int matrows = survyear.length();
      
      NumericMatrix year_mat(matrows, 14);
      year_mat(_, 0) = survyear;
      year_mat(_, 1) = obsyear;
      year_mat(_, 2) = sizeyear;
      year_mat(_, 3) = sizebyear;
      year_mat(_, 4) = sizecyear;
      year_mat(_, 5) = repstyear;
      year_mat(_, 6) = fecyear;
      year_mat(_, 7) = jsurvyear;
      year_mat(_, 8) = jobsyear;
      year_mat(_, 9) = jsizeyear;
      year_mat(_, 10) = jsizebyear;
      year_mat(_, 11) = jsizecyear;
      year_mat(_, 12) = jrepstyear;
      year_mat(_, 13) = jmatstyear;
      
      final_mat = year_mat;
      
    } else if (mat_switch == 2) {
      
      NumericVector survpatch = as<NumericVector>(survproxy["patches"]);
      NumericVector obspatch = as<NumericVector>(obsproxy["patches"]);
      NumericVector sizepatch = as<NumericVector>(sizeproxy["patches"]);
      NumericVector sizebpatch = as<NumericVector>(sizebproxy["patches"]);
      NumericVector sizecpatch = as<NumericVector>(sizecproxy["patches"]);
      NumericVector repstpatch = as<NumericVector>(repstproxy["patches"]);
      NumericVector fecpatch = as<NumericVector>(fecproxy["patches"]);
      NumericVector jsurvpatch = as<NumericVector>(jsurvproxy["patches"]);
      NumericVector jobspatch = as<NumericVector>(jobsproxy["patches"]);
      NumericVector jsizepatch = as<NumericVector>(jsizeproxy["patches"]);
      NumericVector jsizebpatch = as<NumericVector>(jsizebproxy["patches"]);
      NumericVector jsizecpatch = as<NumericVector>(jsizecproxy["patches"]);
      NumericVector jrepstpatch = as<NumericVector>(jrepstproxy["patches"]);
      NumericVector jmatstpatch = as<NumericVector>(jmatstproxy["patches"]);
      
      int matrows = survpatch.length();
      
      NumericMatrix patch_mat(matrows, 14);
      patch_mat(_, 0) = survpatch;
      patch_mat(_, 1) = obspatch;
      patch_mat(_, 2) = sizepatch;
      patch_mat(_, 3) = sizebpatch;
      patch_mat(_, 4) = sizecpatch;
      patch_mat(_, 5) = repstpatch;
      patch_mat(_, 6) = fecpatch;
      patch_mat(_, 7) = jsurvpatch;
      patch_mat(_, 8) = jobspatch;
      patch_mat(_, 9) = jsizepatch;
      patch_mat(_, 10) = jsizebpatch;
      patch_mat(_, 11) = jsizecpatch;
      patch_mat(_, 12) = jrepstpatch;
      patch_mat(_, 13) = jmatstpatch;
      
      final_mat = patch_mat;
    }
    
    return final_mat;
  }
  
  //' Create a Summation of Most Terms Needed in Vital Rate Calculation
  //' 
  //' Function \code{rimeotam()} provides the majority of the work in creating
  //' the linear model sum to be used in vital rate estimation in the MPM. Works
  //' specifically with functions \code{\link{jerzeibalowski}()} and
  //' \code{\link{motherbalowski}()}.
  //' 
  //' @name rimeotam
  //' 
  //' @param maincoefs The coefficients portion of the vital rate model proxy.
  //' @param fl1_i Reproductive status in time \emph{t}*-1.
  //' @param fl2n_i Reproductive status in time \emph{t}.
  //' @param sz1_i Primary size in time \emph{t}-1.
  //' @param sz2o_i Primary size in time \emph{t}.
  //' @param szb1_i Secondary size in time \emph{t}-1.
  //' @param szb2o_i Secondary size in time \emph{t}.
  //' @param szc1_i Tertiary size in time \emph{t}-1.
  //' @param szc2o_i Tertiary size in time \emph{t}.
  //' @param aage2_i Used age in time \emph{t}.
  //' @param inda_1 Value of numeric individual covariate a in time \emph{t}-1.
  //' @param inda_2 Value of numeric individual covariate a in time \emph{t}.
  //' @param indb_1 Value of numeric individual covariate b in time \emph{t}-1.
  //' @param indb_2 Value of numeric individual covariate b in time \emph{t}.
  //' @param indc_1 Value of numeric individual covariate c in time \emph{t}-1.
  //' @param indc_2 Value of numeric individual covariate c in time \emph{t}.
  //' @param anna_1 Value of numeric annual covariate a in time \emph{t}-1.
  //' @param anna_2 Value of numeric annual covariate a in time \emph{t}.
  //' @param annb_1 Value of numeric annual covariate b in time \emph{t}-1.
  //' @param annb_2 Value of numeric annual covariate b in time \emph{t}.
  //' @param annc_1 Value of numeric annual covariate c in time \emph{t}-1.
  //' @param annc_2 Value of numeric annual covariate c in time \emph{t}.
  //' @param used_dens Density value used.
  //' @param zi A logical value indicating whether model coefficients refer to the
  //' zero inflation portion of a model.
  //' 
  //' @return A single numeric value giving the sum of the products of the linear
  //' coefficients and the used status values.
  //' 
  //' @keywords internal
  //' @noRd
  inline double rimeotam(NumericVector maincoefs, const double fl1_i,
    const double fl2n_i, const double sz1_i, const double sz2o_i,
    const double szb1_i, const double szb2o_i, const double szc1_i,
    const double szc2o_i, const double aage2_i, const double inda_1,
    const double inda_2, const double indb_1, const double indb_2,
    const double indc_1, const double indc_2, const double anna_1,
    const double anna_2, const double annb_1, const double annb_2,
    const double annc_1, const double annc_2, const double used_dens, bool zi) {
    
    int add1 {0};
    int add2 {0};
    int add3 {0};
    
    if (zi) {
      add1 = 46;
      add2 = 100;
      add3 = 10;
    }
    
    double parti = maincoefs(0 + add1) + (maincoefs(1 + add1) * fl1_i) +
      (maincoefs(2 + add1) * fl2n_i) + (maincoefs(3 + add1) * sz1_i) +
      (maincoefs(4 + add1) * sz2o_i) + (maincoefs(5 + add1) * fl2n_i * fl1_i) + 
      (maincoefs(6 + add1) * sz2o_i * sz1_i) + (maincoefs(7 + add1) * sz1_i * fl1_i) +
      (maincoefs(8 + add1) * sz2o_i * fl2n_i) + (maincoefs(9 + add1) * sz2o_i * fl1_i) + 
      (maincoefs(10 + add1) * sz1_i * fl2n_i) + (maincoefs(11 + add1) * aage2_i) + 
      (maincoefs(12 + add1) * aage2_i * sz1_i) + (maincoefs(13 + add1) * aage2_i * sz2o_i) + 
      (maincoefs(14 + add1) * aage2_i * fl1_i) + (maincoefs(15 + add1) * aage2_i * fl2n_i) + 
      (maincoefs(16 + add1) * inda_2) + (maincoefs(17 + add1) * indb_2) +
      (maincoefs(18 + add1) * indc_2) + (maincoefs(19 + add1) * inda_1) +
      (maincoefs(20 + add1) * indb_1) + (maincoefs(21 + add1) * indc_1) + 
      (maincoefs(22 + add1) * inda_2 * sz2o_i) + (maincoefs(23 + add1) * indb_2 * sz2o_i) + 
      (maincoefs(24 + add1) * indc_2 * sz2o_i) + (maincoefs(25 + add1) * inda_2 * fl2n_i) + 
      (maincoefs(26 + add1) * indb_2 * fl2n_i) + (maincoefs(27 + add1) * indc_2 * fl2n_i) + 
      (maincoefs(28 + add1) * inda_1 * sz1_i) + (maincoefs(29 + add1) * indb_1 * sz1_i) +
      (maincoefs(30 + add1) * indc_1 * sz1_i) + (maincoefs(31 + add1) * inda_1 * fl1_i) + 
      (maincoefs(32 + add1) * indb_1 * fl1_i) + (maincoefs(33 + add1) * indc_1 * fl1_i) + 
      (maincoefs(34 + add1) * inda_2 * indb_2) + (maincoefs(35 + add1) * inda_2 * indc_2) + 
      (maincoefs(36 + add1) * indb_2 * indc_2) + (maincoefs(37 + add1) * inda_1 * indb_1) + 
      (maincoefs(38 + add1) * inda_1 * indc_1) + (maincoefs(39 + add1) * indb_1 * indc_1) + 
      (maincoefs(40 + add1) * inda_2 * indb_1) + (maincoefs(41 + add1) * inda_1 * indb_2) +
      (maincoefs(42 + add1) * inda_2 * indc_1) + (maincoefs(43 + add1) * inda_1 * indc_2) + 
      (maincoefs(44 + add1) * indb_2 * indc_1) + (maincoefs(45 + add1) * indb_1 * indc_2);
      
    double partii = (maincoefs(100 + add2) * szb2o_i) + (maincoefs(101 + add2) * szb1_i) + 
      (maincoefs(102 + add2) * szc2o_i) + (maincoefs(103 + add2) * szc1_i) + 
      (maincoefs(104 + add2) * used_dens) + (maincoefs(105 + add2) * szb1_i * szb2o_i) + 
      (maincoefs(106 + add2) * szc1_i * szc2o_i) + (maincoefs(107 + add2) * sz1_i * szb1_i) + 
      (maincoefs(108 + add2) * sz1_i * szc1_i) + (maincoefs(109 + add2) * szb1_i * szc1_i) + 
      (maincoefs(110 + add2) * sz2o_i * szb2o_i) + (maincoefs(111 + add2) * sz2o_i * szc2o_i) + 
      (maincoefs(112 + add2) * szb2o_i * szc2o_i) + (maincoefs(113 + add2) * sz1_i * szb2o_i) + 
      (maincoefs(114 + add2) * sz1_i * szc2o_i) + (maincoefs(115 + add2) * szb1_i * szc2o_i) + 
      (maincoefs(116 + add2) * sz2o_i * szb1_i) + (maincoefs(117 + add2) * sz2o_i * szc1_i) + 
      (maincoefs(118 + add2) * szb2o_i * szc1_i) + (maincoefs(119 + add2) * sz2o_i * used_dens) + 
      (maincoefs(120 + add2) * szb2o_i * used_dens) + (maincoefs(121 + add2) * szc2o_i * used_dens) + 
      (maincoefs(122 + add2) * sz1_i * used_dens) + (maincoefs(123 + add2) * szb1_i * used_dens) + 
      (maincoefs(124 + add2) * szc1_i * used_dens) + (maincoefs(125 + add2) * fl2n_i * used_dens) + 
      (maincoefs(126 + add2) * fl1_i * used_dens) + (maincoefs(127 + add2) * szb2o_i * fl2n_i) + 
      (maincoefs(128 + add2) * szc2o_i * fl2n_i) + (maincoefs(129 + add2) * inda_2 * aage2_i) + 
      (maincoefs(130 + add2) * szb1_i * fl1_i) + (maincoefs(131 + add2) * szb2o_i * fl1_i) + 
      (maincoefs(132 + add2) * szb1_i * fl2n_i) + (maincoefs(133 + add2) * szc1_i * fl1_i) + 
      (maincoefs(134 + add2) * szc2o_i * fl1_i) + (maincoefs(135 + add2) * szc1_i * fl2n_i) + 
      (maincoefs(136 + add2) * szb2o_i * aage2_i) + (maincoefs(137 + add2) * szc2o_i * aage2_i) + 
      (maincoefs(138 + add2) * used_dens * aage2_i) + (maincoefs(139 + add2) * szb1_i * aage2_i) + 
      (maincoefs(140 + add2) * szc1_i * aage2_i) + (maincoefs(141 + add2) * inda_2 * szb2o_i) + 
      (maincoefs(142 + add2) * inda_2 * szc2o_i) + (maincoefs(143 + add2) * inda_2 * used_dens) + 
      (maincoefs(144 + add2) * inda_1 * szb1_i) + (maincoefs(145 + add2) * inda_1 * szc1_i) + 
      (maincoefs(146 + add2) * inda_1 * szb2o_i) + (maincoefs(147 + add2) * inda_1 * szc2o_i) + 
      (maincoefs(148 + add2) * inda_2 * szb1_i) + (maincoefs(149 + add2) * inda_2 * szc1_i);
    
    double partiii = (maincoefs(150 + add2) * inda_1 * used_dens) +
      (maincoefs(151 + add2) * indb_2 * szb2o_i) + (maincoefs(152 + add2) * indb_2 * szc2o_i) + 
      (maincoefs(153 + add2) * indb_2 * used_dens) + (maincoefs(154 + add2) * indb_1 * szb1_i) + 
      (maincoefs(155 + add2) * indb_1 * szc1_i) + (maincoefs(156 + add2) * indb_1 * szb2o_i) + 
      (maincoefs(157 + add2) * indb_1 * szc2o_i) + (maincoefs(158 + add2) * indb_2 * szb1_i) + 
      (maincoefs(159 + add2) * indb_2 * szc1_i) + (maincoefs(160 + add2) * indb_1 * used_dens) + 
      (maincoefs(161 + add2) * indc_2 * szb2o_i) + (maincoefs(162 + add2) * indc_2 * szc2o_i) + 
      (maincoefs(163 + add2) * indc_2 * used_dens) + (maincoefs(164 + add2) * indc_1 * szb1_i) + 
      (maincoefs(165 + add2) * indc_1 * szc1_i) + (maincoefs(166 + add2) * indc_1 * szb2o_i) + 
      (maincoefs(167 + add2) * indc_1 * szc2o_i) + (maincoefs(168 + add2) * indc_2 * szb1_i) + 
      (maincoefs(169 + add2) * indc_2 * szc1_i) + (maincoefs(170 + add2) * indc_1 * used_dens) + 
      (maincoefs(171 + add2) * inda_2 * sz1_i) + (maincoefs(172 + add2) * indb_2 * sz1_i) + 
      (maincoefs(173 + add2) * indc_2 * sz1_i) + (maincoefs(174 + add2) * inda_1 * sz2o_i) + 
      (maincoefs(175 + add2) * indb_1 * sz2o_i) + (maincoefs(176 + add2) * indc_1 * sz2o_i) + 
      (maincoefs(177 + add2) * inda_2 * fl1_i) + (maincoefs(178 + add2) * indb_2 * fl1_i) + 
      (maincoefs(179 + add2) * indc_2 * fl1_i) + (maincoefs(180 + add2) * inda_1 * fl2n_i) + 
      (maincoefs(181 + add2) * indb_1 * fl2n_i) + (maincoefs(182 + add2) * indc_1 * fl2n_i) + 
      (maincoefs(183 + add2) * indc_1 * aage2_i) + (maincoefs(184 + add2) * inda_2 * inda_1) + 
      (maincoefs(185 + add2) * indb_2 * indb_1) + (maincoefs(186 + add2) * indc_2 * indc_1) + 
      (maincoefs(187 + add2) * sz2o_i * anna_2) + (maincoefs(188 + add2) * sz2o_i * anna_1) + 
      (maincoefs(189 + add2) * sz2o_i * annb_2) + (maincoefs(190 + add2) * sz2o_i * annb_1) + 
      (maincoefs(191 + add2) * sz2o_i * annc_2) + (maincoefs(192 + add2) * sz2o_i * annc_1) + 
      (maincoefs(193 + add2) * sz1_i * anna_2) + (maincoefs(194 + add2) * sz1_i * anna_1) + 
      (maincoefs(195 + add2) * sz1_i * annb_2) + (maincoefs(196 + add2) * sz1_i * annb_1) + 
      (maincoefs(197 + add2) * sz1_i * annc_2) + (maincoefs(198 + add2) * sz1_i * annc_1) + 
      (maincoefs(199 + add2) * szb2o_i * anna_2);
    
    double partiv = (maincoefs(300 + add2) * szb2o_i * anna_1) + 
      (maincoefs(301 + add2) * szb2o_i * annb_2) + (maincoefs(302 + add2) * szb2o_i * annb_1) + 
      (maincoefs(303 + add2) * szb2o_i * annc_2) + (maincoefs(304 + add2) * szb2o_i * annc_1) + 
      (maincoefs(305 + add2) * szb1_i * anna_2) + (maincoefs(306 + add2) * szb1_i * anna_1) + 
      (maincoefs(307 + add2) * szb1_i * annb_2) + (maincoefs(308 + add2) * szb1_i * annb_1) + 
      (maincoefs(309 + add2) * szb1_i * annc_2) + (maincoefs(310 + add2) * szb1_i * annc_1) + 
      (maincoefs(311 + add2) * szc2o_i * anna_2) + (maincoefs(312 + add2) * szc2o_i * anna_1) + 
      (maincoefs(313 + add2) * szc2o_i * annb_2) + (maincoefs(314 + add2) * szc2o_i * annb_1) + 
      (maincoefs(315 + add2) * szc2o_i * annc_2) + (maincoefs(316 + add2) * szc2o_i * annc_1) + 
      (maincoefs(317 + add2) * szc1_i * anna_2) + (maincoefs(318 + add2) * szc1_i * anna_1) + 
      (maincoefs(319 + add2) * szc1_i * annb_2) + (maincoefs(320 + add2) * szc1_i * annb_1) + 
      (maincoefs(321 + add2) * szc1_i * annc_2) + (maincoefs(322 + add2) * szc1_i * annc_1) + 
      (maincoefs(323 + add2) * fl2n_i * anna_2) + (maincoefs(324 + add2) * fl2n_i * anna_1) + 
      (maincoefs(325 + add2) * fl2n_i * annb_2) + (maincoefs(326 + add2) * fl2n_i * annb_1) + 
      (maincoefs(327 + add2) * fl2n_i * annc_2) + (maincoefs(328 + add2) * fl2n_i * annc_1) + 
      (maincoefs(329 + add2) * fl1_i * anna_2) + (maincoefs(330 + add2) * fl1_i * anna_1) + 
      (maincoefs(331 + add2) * fl1_i * annb_2) + (maincoefs(332 + add2) * fl1_i * annb_1) + 
      (maincoefs(333 + add2) * fl1_i * annc_2) + (maincoefs(334 + add2) * fl1_i * annc_1) + 
      (maincoefs(335 + add2) * aage2_i * anna_2) + (maincoefs(336 + add2) * aage2_i * anna_1) + 
      (maincoefs(337 + add2) * aage2_i * annb_2) + (maincoefs(338 + add2) * aage2_i * annb_1) + 
      (maincoefs(339 + add2) * aage2_i * annc_2) + (maincoefs(340 + add2) * aage2_i * annc_1) + 
      (maincoefs(341 + add2) * used_dens * anna_2) + (maincoefs(342 + add2) * used_dens * anna_1) + 
      (maincoefs(343 + add2) * used_dens * annb_2) + (maincoefs(344 + add2) * used_dens * annb_1) + 
      (maincoefs(345 + add2) * used_dens * annc_2) + (maincoefs(346 + add2) * used_dens * annc_1) + 
      (maincoefs(347 + add2) * inda_2 * anna_2) + (maincoefs(348 + add2) * inda_2 * anna_1) + 
      (maincoefs(349 + add2) * inda_2 * annb_2) + (maincoefs(350 + add2) * inda_2 * annb_1) + 
      (maincoefs(351 + add2) * inda_2 * annc_2) + (maincoefs(352 + add2) * inda_2 * annc_1);
    
    double partv = (maincoefs(353 + add2) * inda_1 * anna_2) + (maincoefs(354 + add2) * inda_1 * anna_1) + 
      (maincoefs(355 + add2) * inda_1 * annb_2) + (maincoefs(356 + add2) * inda_1 * annb_1) + 
      (maincoefs(357 + add2) * inda_1 * annc_2) + (maincoefs(358 + add2) * inda_1 * annc_1) + 
      (maincoefs(359 + add2) * indb_2 * anna_2) + (maincoefs(360 + add2) * indb_2 * anna_1) + 
      (maincoefs(361 + add2) * indb_2 * annb_2) + (maincoefs(362 + add2) * indb_2 * annb_1) + 
      (maincoefs(363 + add2) * indb_2 * annc_2) + (maincoefs(364 + add2) * indb_2 * annc_1) + 
      (maincoefs(365 + add2) * indb_1 * anna_2) + (maincoefs(366 + add2) * indb_1 * anna_1) + 
      (maincoefs(367 + add2) * indb_1 * annb_2) + (maincoefs(368 + add2) * indb_1 * annb_1) + 
      (maincoefs(369 + add2) * indb_1 * annc_2) + (maincoefs(370 + add2) * indb_1 * annc_1) + 
      (maincoefs(371 + add2) * indc_2 * anna_2) + (maincoefs(372 + add2) * indc_2 * anna_1) + 
      (maincoefs(373 + add2) * indc_2 * annb_2) + (maincoefs(374 + add2) * indc_2 * annb_1) + 
      (maincoefs(375 + add2) * indc_2 * annc_2) + (maincoefs(376 + add2) * indc_2 * annc_1) + 
      (maincoefs(377 + add2) * indc_1 * anna_2) + (maincoefs(378 + add2) * indc_1 * anna_1) + 
      (maincoefs(379 + add2) * indc_1 * annb_2) + (maincoefs(380 + add2) * indc_1 * annb_1) + 
      (maincoefs(381 + add2) * indc_1 * annc_2) + (maincoefs(382 + add2) * indc_1 * annc_1) + 
      (maincoefs(383 + add2) * anna_2) + (maincoefs(384 + add2) * anna_2* anna_1) + 
      (maincoefs(385 + add2) * anna_2 * annb_2) + (maincoefs(386 + add2) * anna_2 * annb_1) + 
      (maincoefs(387 + add2) * anna_2 * annc_2) + (maincoefs(388 + add2) * anna_2 * annc_1) + 
      (maincoefs(389 + add2) * anna_1) + (maincoefs(390 + add2) * anna_1 * annb_2) + 
      (maincoefs(391 + add2) * anna_1 * annb_1) + (maincoefs(392 + add2) * anna_1 * annc_2) + 
      (maincoefs(393 + add2) * anna_1 * annc_1) + (maincoefs(394 + add2) * annb_2) + 
      (maincoefs(395 + add2) * annb_2 * annb_1) + (maincoefs(396 + add2) * annb_2 * annc_2) + 
      (maincoefs(397 + add2) * annb_2 * annc_1) + (maincoefs(398 + add2) * annb_1) +
      (maincoefs(399 + add2) * annb_1 * annc_2);
    
    double partvi = (maincoefs(500 + add3) * annb_1 * annc_1) + (maincoefs(501 + add3) * annc_2) + 
      (maincoefs(502 + add3) * annc_2 * annc_1) + (maincoefs(503 + add3) * annc_1);
    
    double albatross = parti + partii + partiii + partiv + partv + partvi;
    
    return albatross;
  }
  
  //' Count Elements in Each Random Individual Covariate Portion of Model
  //' 
  //' Function \code{foi_counter()} counts the number of elements in each random
  //' individual covariate and returns that as a vector.
  //' 
  //' @name foi_counter
  //' 
  //' @param modelproxy A list holding the contents of a model processed with
  //' function \code{modelextract()}
  //' @param zi A logical value indicating whether to focus on the zero-inflation
  //' parameters.
  //' 
  //' @return A 6 element vector holding the numbers of elements in each random
  //' individual covariate in a model (either the cont portion or the zi portion).
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::ivec foi_counter(List modelproxy, bool zi) {
    
    arma::ivec return_vec(6, fill::zeros);
    
    if (!zi) {
      arma::vec modelinda2r = as<arma::vec>(modelproxy["indcova2s"]);
      arma::vec modelinda1r = as<arma::vec>(modelproxy["indcova1s"]);
      arma::vec modelindb2r = as<arma::vec>(modelproxy["indcovb2s"]);
      arma::vec modelindb1r = as<arma::vec>(modelproxy["indcovb1s"]);
      arma::vec modelindc2r = as<arma::vec>(modelproxy["indcovc2s"]);
      arma::vec modelindc1r = as<arma::vec>(modelproxy["indcovc1s"]);
      
      int v1_l = static_cast<int>(modelinda2r.n_elem);
      int v2_l = static_cast<int>(modelinda1r.n_elem);
      int v3_l = static_cast<int>(modelindb2r.n_elem);
      int v4_l = static_cast<int>(modelindb1r.n_elem);
      int v5_l = static_cast<int>(modelindc2r.n_elem);
      int v6_l = static_cast<int>(modelindc1r.n_elem);
      
      return_vec = {v1_l, v2_l, v3_l, v4_l, v5_l, v6_l};
    } else {
      arma::vec modelinda2r = as<arma::vec>(modelproxy["zeroindcova2s"]);
      arma::vec modelinda1r = as<arma::vec>(modelproxy["zeroindcova1s"]);
      arma::vec modelindb2r = as<arma::vec>(modelproxy["zeroindcovb2s"]);
      arma::vec modelindb1r = as<arma::vec>(modelproxy["zeroindcovb1s"]);
      arma::vec modelindc2r = as<arma::vec>(modelproxy["zeroindcovc2s"]);
      arma::vec modelindc1r = as<arma::vec>(modelproxy["zeroindcovc1s"]);
      
      int v1_l = static_cast<int>(modelinda2r.n_elem);
      int v2_l = static_cast<int>(modelinda1r.n_elem);
      int v3_l = static_cast<int>(modelindb2r.n_elem);
      int v4_l = static_cast<int>(modelindb1r.n_elem);
      int v5_l = static_cast<int>(modelindc2r.n_elem);
      int v6_l = static_cast<int>(modelindc1r.n_elem);
      
      return_vec = {v1_l, v2_l, v3_l, v4_l, v5_l, v6_l};
    }
    
    return return_vec;
  }
  
  //' Create Vector of Random Individual Covariate Terms
  //' 
  //' Function \code{flightoficarus()} creates vectors of random covariate
  //' terms.
  //' 
  //' @name flightoficarus
  //' 
  //' @param modelproxy A model proxy list extracted with function
  //' \code{modelextract()}.
  //' 
  //' @return A vector of numeric values for random categorical terms. The order
  //' is: 1) cov a time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1,
  //' 5) cov c time 2, and 6) cov c time 1. Rows may vary, but must be the same
  //' length for each model.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::NumericVector flightoficarus(List modelproxy) {
    NumericVector modelinda2r = as<NumericVector>(modelproxy["indcova2s"]);
    NumericVector modelinda1r = as<NumericVector>(modelproxy["indcova1s"]);
    NumericVector modelindb2r = as<NumericVector>(modelproxy["indcovb2s"]);
    NumericVector modelindb1r = as<NumericVector>(modelproxy["indcovb1s"]);
    NumericVector modelindc2r = as<NumericVector>(modelproxy["indcovc2s"]);
    NumericVector modelindc1r = as<NumericVector>(modelproxy["indcovc1s"]);
    
    int v1_l = modelinda2r.length();
    int v2_l = modelinda1r.length();
    int v3_l = modelindb2r.length();
    int v4_l = modelindb1r.length();
    int v5_l = modelindc2r.length();
    int v6_l = modelindc1r.length();
    int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
    
    NumericVector final_vec(vec_length);
    int all_counter {0};
    
    for (int i = 0; i < v1_l; i++) {
      final_vec(all_counter) = modelinda2r(i);
      all_counter++;
    }
    for (int i = 0; i < v2_l; i++) {
      final_vec(all_counter) = modelinda1r(i);
      all_counter++;
    }
    for (int i = 0; i < v3_l; i++) {
      final_vec(all_counter) = modelindb2r(i);
      all_counter++;
    }
    for (int i = 0; i < v4_l; i++) {
      final_vec(all_counter) = modelindb1r(i);
      all_counter++;
    }
    for (int i = 0; i < v5_l; i++) {
      final_vec(all_counter) = modelindc2r(i);
      all_counter++;
    }
    for (int i = 0; i < v6_l; i++) {
      final_vec(all_counter) = modelindc1r(i);
      all_counter++;
    }
    
    return final_vec;
  }
  
  //' Create Concatenated Vector of Random Individual Covariate Term Names
  //' 
  //' Function \code{bootson()} creates a concatenated string vector holding all
  //' covariate term names.
  //' 
  //' @name bootson
  //' 
  //' @param modelproxy A model proxy list extracted with function
  //' \code{modelextract()}.
  //' 
  //' @return A vector holding all covariate name terms. The order is: 1) cov a
  //' time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1, 5) cov c time 2,
  //' and 6) cov c time 1. Note that the element order is the same as in function
  //' \code{\link{.flightoficarus}()}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::StringVector bootson(List modelproxy) {
    NumericVector modelinda2r_df = as<NumericVector>(modelproxy["indcova2s"]);
    NumericVector modelinda1r_df = as<NumericVector>(modelproxy["indcova1s"]);
    NumericVector modelindb2r_df = as<NumericVector>(modelproxy["indcovb2s"]);
    NumericVector modelindb1r_df = as<NumericVector>(modelproxy["indcovb1s"]);
    NumericVector modelindc2r_df = as<NumericVector>(modelproxy["indcovc2s"]);
    NumericVector modelindc1r_df = as<NumericVector>(modelproxy["indcovc1s"]);
    
    StringVector modelinda2r_rownames = modelinda2r_df.attr("names");
    StringVector modelinda1r_rownames = modelinda1r_df.attr("names");
    StringVector modelindb2r_rownames = modelindb2r_df.attr("names");
    StringVector modelindb1r_rownames = modelindb1r_df.attr("names");
    StringVector modelindc2r_rownames = modelindc2r_df.attr("names");
    StringVector modelindc1r_rownames = modelindc1r_df.attr("names");
    
    int v1_l = modelinda2r_rownames.length();
    int v2_l = modelinda1r_rownames.length();
    int v3_l = modelindb2r_rownames.length();
    int v4_l = modelindb1r_rownames.length();
    int v5_l = modelindc2r_rownames.length();
    int v6_l = modelindc1r_rownames.length();
    int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
  
    StringVector final_vec(vec_length);
    int all_counter {0};
    
    for (int i = 0; i < v1_l; i++) {
      final_vec(all_counter) = modelinda2r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v2_l; i++) {
      final_vec(all_counter) = modelinda1r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v3_l; i++) {
      final_vec(all_counter) = modelindb2r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v4_l; i++) {
      final_vec(all_counter) = modelindb1r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v5_l; i++) {
      final_vec(all_counter) = modelindc2r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v6_l; i++) {
      final_vec(all_counter) = modelindc1r_rownames(i);
      all_counter++;
    }
    
    return final_vec;
  }
  
  //' Create Vector of Random Individual Covariate Terms for Zero-Inflated Models
  //' 
  //' Function \code{zero_flightoficarus()} creates vectors of random covariate
  //' terms from the binomial portion of a zero-inflated model.
  //' 
  //' @name zero_flightoficarus
  //' 
  //' @param modelproxy A model proxy list extracted with function
  //' \code{modelextract()}.
  //' 
  //' @return A vector of numeric values for random categorical terms. The order
  //' is: 1) cov a time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1,
  //' 5) cov c time 2, and 6) cov c time 1. Rows may vary, but must be the same
  //' length for each model.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::NumericVector zero_flightoficarus(List modelproxy) {
    NumericVector modelinda2r = as<NumericVector>(modelproxy["zeroindcova2s"]);
    NumericVector modelinda1r = as<NumericVector>(modelproxy["zeroindcova1s"]);
    NumericVector modelindb2r = as<NumericVector>(modelproxy["zeroindcovb2s"]);
    NumericVector modelindb1r = as<NumericVector>(modelproxy["zeroindcovb1s"]);
    NumericVector modelindc2r = as<NumericVector>(modelproxy["zeroindcovc2s"]);
    NumericVector modelindc1r = as<NumericVector>(modelproxy["zeroindcovc1s"]);
    
    int v1_l = modelinda2r.length();
    int v2_l = modelinda1r.length();
    int v3_l = modelindb2r.length();
    int v4_l = modelindb1r.length();
    int v5_l = modelindc2r.length();
    int v6_l = modelindc1r.length();
    int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
    
    NumericVector final_vec(vec_length);
    int all_counter {0};
    
    for (int i = 0; i < v1_l; i++) {
      final_vec(all_counter) = modelinda2r(i);
      all_counter++;
    }
    for (int i = 0; i < v2_l; i++) {
      final_vec(all_counter) = modelinda1r(i);
      all_counter++;
    }
    for (int i = 0; i < v3_l; i++) {
      final_vec(all_counter) = modelindb2r(i);
      all_counter++;
    }
    for (int i = 0; i < v4_l; i++) {
      final_vec(all_counter) = modelindb1r(i);
      all_counter++;
    }
    for (int i = 0; i < v5_l; i++) {
      final_vec(all_counter) = modelindc2r(i);
      all_counter++;
    }
    for (int i = 0; i < v6_l; i++) {
      final_vec(all_counter) = modelindc1r(i);
      all_counter++;
    }
    
    return final_vec;
  }
  
  //' Create Concatenated Vector of Random Individual Covariate Term Names from
  //' a Zero-Inflated Model
  //' 
  //' Function \code{zero_bootson()} creates a concatenated string vector holding
  //' all covariate term names from the binomial portion of a zero-inflated model.
  //' 
  //' @name zero_bootson
  //' 
  //' @param modelproxy A model proxy list extracted with function
  //' \code{modelextract()}.
  //' 
  //' @return A vector holding all covariate name terms. The order is: 1) cov a
  //' time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1, 5) cov c time 2,
  //' and 6) cov c time 1. Note that the element order is the same as in function
  //' \code{\link{.zero_flightoficarus}()}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::StringVector zero_bootson(List modelproxy) {
    NumericVector modelinda2r_df = as<NumericVector>(modelproxy["zeroindcova2s"]);
    NumericVector modelinda1r_df = as<NumericVector>(modelproxy["zeroindcova1s"]);
    NumericVector modelindb2r_df = as<NumericVector>(modelproxy["zeroindcovb2s"]);
    NumericVector modelindb1r_df = as<NumericVector>(modelproxy["zeroindcovb1s"]);
    NumericVector modelindc2r_df = as<NumericVector>(modelproxy["zeroindcovc2s"]);
    NumericVector modelindc1r_df = as<NumericVector>(modelproxy["zeroindcovc1s"]);
    
    StringVector modelinda2r_rownames = modelinda2r_df.attr("names");
    StringVector modelinda1r_rownames = modelinda1r_df.attr("names");
    StringVector modelindb2r_rownames = modelindb2r_df.attr("names");
    StringVector modelindb1r_rownames = modelindb1r_df.attr("names");
    StringVector modelindc2r_rownames = modelindc2r_df.attr("names");
    StringVector modelindc1r_rownames = modelindc1r_df.attr("names");
    
    int v1_l = modelinda2r_rownames.length();
    int v2_l = modelinda1r_rownames.length();
    int v3_l = modelindb2r_rownames.length();
    int v4_l = modelindb1r_rownames.length();
    int v5_l = modelindc2r_rownames.length();
    int v6_l = modelindc1r_rownames.length();
    int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
  
    StringVector final_vec(vec_length);
    int all_counter {0};
    
    for (int i = 0; i < v1_l; i++) {
      final_vec(all_counter) = modelinda2r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v2_l; i++) {
      final_vec(all_counter) = modelinda1r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v3_l; i++) {
      final_vec(all_counter) = modelindb2r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v4_l; i++) {
      final_vec(all_counter) = modelindb1r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v5_l; i++) {
      final_vec(all_counter) = modelindc2r_rownames(i);
      all_counter++;
    }
    for (int i = 0; i < v6_l; i++) {
      final_vec(all_counter) = modelindc1r_rownames(i);
      all_counter++;
    }
    
    return final_vec;
  }
  
  //' Create Index of Element Numbers for Random Individual Covariate Terms
  //' 
  //' Function \code{foi_index()} creates a matrix indexing the end points of
  //' each random individual covariate in the utilized vectors.
  //' 
  //' @name foi_index
  //' 
  //' @param surv_proxy Adult survival model proxy.
  //' @param obs_proxy Adult observation status model proxy.
  //' @param size_proxy Adult primary size model proxy.
  //' @param sizeb_proxy Adult secondary size model proxy.
  //' @param sizec_proxy Adult tertiary size model proxy.
  //' @param repst_proxy Adult reproductive status model proxy.
  //' @param fec_proxy Adult fecundity model proxy.
  //' @param jsurv_proxy Juvenile survival model proxy.
  //' @param jobs_proxy Juvenile observation status model proxy.
  //' @param jsize_proxy Juvenile primary size model proxy.
  //' @param jsizeb_proxy Juvenile secondary size model proxy.
  //' @param jsizec_proxy Juvenile tertiary size model proxy.
  //' @param jrepst_proxy Juvenile reproductive status model proxy.
  //' @param jmatst_proxy Juvenile maturity status model proxy.
  //' 
  //' @return An integer matrix with 6 rows and 20 columns. The columns contain
  //' the number of elements in each random individual covariate term, with the
  //' row order being: 1) cov a t2, 2) cov a t1, 3) cov b t2, 4) cov b t1,
  //' 5) cov c t2, and 6) cov c t1.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::imat foi_index(List surv_proxy, List obs_proxy, List size_proxy, 
    List sizeb_proxy, List sizec_proxy, List repst_proxy, List fec_proxy,
    List jsurv_proxy, List jobs_proxy, List jsize_proxy, List jsizeb_proxy,
    List jsizec_proxy, List jrepst_proxy, List jmatst_proxy) {
    
    arma::ivec surv_fc = foi_counter(surv_proxy, false);
    arma::ivec obs_fc = foi_counter(obs_proxy, false);
    arma::ivec size_fc = foi_counter(size_proxy, false);
    arma::ivec sizeb_fc = foi_counter(sizeb_proxy, false);
    arma::ivec sizec_fc = foi_counter(sizec_proxy, false);
    arma::ivec repst_fc = foi_counter(repst_proxy, false);
    arma::ivec fec_fc = foi_counter(fec_proxy, false);
    arma::ivec jsurv_fc = foi_counter(jsurv_proxy, false);
    arma::ivec jobs_fc = foi_counter(jobs_proxy, false);
    arma::ivec jsize_fc = foi_counter(jsize_proxy, false);
    arma::ivec jsizeb_fc = foi_counter(jsizeb_proxy, false);
    arma::ivec jsizec_fc = foi_counter(jsizec_proxy, false);
    arma::ivec jrepst_fc = foi_counter(jrepst_proxy, false);
    arma::ivec jmatst_fc = foi_counter(jmatst_proxy, false);
    arma::ivec size_fc_zi = foi_counter(size_proxy, true);
    arma::ivec sizeb_fc_zi = foi_counter(sizeb_proxy, true);
    arma::ivec sizec_fc_zi = foi_counter(sizec_proxy, true);
    arma::ivec fec_fc_zi = foi_counter(fec_proxy, true);
    arma::ivec jsize_fc_zi = foi_counter(jsize_proxy, true);
    arma::ivec jsizeb_fc_zi = foi_counter(jsizeb_proxy, true);
    arma::ivec jsizec_fc_zi = foi_counter(jsizec_proxy, true);
    
    arma::imat final_mat(6, 21, fill::zeros);
    
    for (int i = 0; i < 6; i++) {
      final_mat(i, 0) = surv_fc(i);
      final_mat(i, 1) = obs_fc(i);
      final_mat(i, 2) = size_fc(i);
      final_mat(i, 3) = sizeb_fc(i);
      final_mat(i, 4) = sizec_fc(i);
      final_mat(i, 5) = repst_fc(i);
      final_mat(i, 6) = fec_fc(i);
      final_mat(i, 7) = jsurv_fc(i);
      final_mat(i, 8) = jobs_fc(i);
      final_mat(i, 9) = jsize_fc(i);
      final_mat(i, 10) = jsizeb_fc(i);
      final_mat(i, 11) = jsizec_fc(i);
      final_mat(i, 12) = jrepst_fc(i);
      final_mat(i, 13) = jmatst_fc(i);
      final_mat(i, 14) = size_fc_zi(i);
      final_mat(i, 15) = sizeb_fc_zi(i);
      final_mat(i, 16) = sizec_fc_zi(i);
      final_mat(i, 17) = fec_fc_zi(i);
      final_mat(i, 18) = jsize_fc_zi(i);
      final_mat(i, 19) = jsizeb_fc_zi(i);
      final_mat(i, 20) = jsizec_fc_zi(i);
    }
    
    return final_mat;
  }
  
  //' Creates Matrices of Year and Patch Terms in Leslie Models
  //' 
  //' Function \code{revelations_leslie()} creates a matrix holding either the
  //' year or patch coefficients from Leslie vital rate models. This reduces
  //' memory load in function \code{\link{motherbalowski}()}.
  //' 
  //' @name revelations_leslie
  //' 
  //' @param survproxy The proxy vital rate model covering survival from the main
  //' matrix estimator function.
  //' @param fecproxy The proxy vital rate model covering fecundity from the main
  //' matrix estimator function.
  //' 
  //' @return A matrix with 2 columns corresponding to the number of vital rates
  //' and number of columns equal to the number of year or patches.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::NumericMatrix revelations_leslie(List survproxy, List fecproxy, int mat_switch) {
    
    NumericMatrix final_mat;
    
    if (mat_switch == 1) {
      NumericVector survyear = as<NumericVector>(survproxy["years"]);
      NumericVector fecyear = as<NumericVector>(fecproxy["years"]);
      
      int matrows = survyear.length();
      
      NumericMatrix year_mat(matrows, 2);
      year_mat(_, 0) = survyear;
      year_mat(_, 1) = fecyear;
      
      final_mat = year_mat;
      
    } else if (mat_switch == 2) {
      
      NumericVector survpatch = as<NumericVector>(survproxy["patches"]);
      NumericVector fecpatch = as<NumericVector>(fecproxy["patches"]);
      
      int matrows = survpatch.length();
      
      NumericMatrix patch_mat(matrows, 2);
      patch_mat(_, 0) = survpatch;
      patch_mat(_, 1) = fecpatch;
  
      final_mat = patch_mat;
    }
    
    return final_mat;
  }
  
  //' Create Index of Element Numbers for Random Individual Covariate Terms in
  //' Leslie Models
  //' 
  //' Function \code{foi_index_leslie()} creates a matrix indexing the end points
  //' of each random individual covariate in the utilized vectors. Used in
  //' function \code{\link{motherbalowski}()}.
  //' 
  //' @name foi_index_leslie
  //' 
  //' @param surv_proxy Adult survival model proxy.
  //' @param fec_proxy Adult fecundity model proxy.
  //' 
  //' @return An integer matrix with 6 rows and 3 columns. The columns contain the
  //' number of elements in each random individual covariate term, with the row
  //' order being: 1) cov a t2, 2) cov a t1, 3) cov b t2, 4) cov b t1,
  //' 5) cov c t2, and 6) cov c t1.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::imat foi_index_leslie(List surv_proxy, List fec_proxy) {
    
    arma::ivec surv_fc = foi_counter(surv_proxy, false);
    arma::ivec fec_fc = foi_counter(fec_proxy, false);
    arma::ivec fec_fc_zi = foi_counter(fec_proxy, true);
    
    arma::imat final_mat(6, 3, fill::zeros);
    
    for (int i = 0; i < 6; i++) {
      final_mat(i, 0) = surv_fc(i);
      final_mat(i, 1) = fec_fc(i);
      final_mat(i, 2) = fec_fc_zi(i);
    }
    
    return final_mat;
  }
  
  //' Extract Coefficients From Linear Vital Rate Models
  //' 
  //' Function \code{modelextract()} extracts coefficient values from linear
  //' models estimated through various linear modeling functions in R, to
  //' estimate vital rates in \code{lefko3}. Used to supply coefficients to
  //' \code{mpm_create()}, and through that function also to \code{flefko3()},
  //' \code{flefko2()}, \code{fleslie()}, and \code{aflefko2()}.
  //' 
  //' @name modelextract
  //' 
  //' @param object A linear model estimated through one of the methods used in
  //' function \code{modelsearch()}, or a \code{vrm_input} object.
  //' @param paramnames Data frame giving the names of standard coefficients
  //' required by matrix creation functions.
  //' @param mainyears A vector of the names of the monitoring occasions.
  //' @param mainpatches A vector of the names of the patches. Should be
  //' \code{NA} if no patches specified.
  //' @param maingroups A vector of the names of all stage groups.
  //' @param mainindcova A vector denoting values of individual covariate
  //' \code{a}, when that individual covariate is categorical.
  //' @param mainindcovb A vector denoting values of individual covariate
  //' \code{b}, when that individual covariate is categorical.
  //' @param mainindcovc A vector denoting values of individual covariate
  //' \code{c}, when that individual covariate is categorical.
  //' @param nodata A logical value used to determine whether to use
  //' \code{vrm_input} methods. Defaults to \code{FALSE}, in which case models
  //' were developed with function \code{modelsearch()}.
  //' 
  //' @return This function returns a list with the following elements:
  //' \item{coefficients}{Vector of fixed effect coefficients.}
  //' \item{years}{Vector of occasion coefficients, typically random.}
  //' \item{zeroyear}{Vector of zero-inflated occasion coefficients, typically
  //' random.}
  //' \item{patches}{Vector of patch coefficients, typically random.}
  //' \item{zeropatch}{Vector of zero-inflated patch coefficients, typically
  //' random.}
  //' \item{groups2}{Vector of group coefficients for time t.}
  //' \item{groups1}{Vector of group coefficients for time t-1.}
  //' \item{zerogroups2}{Vector of zero-inflated group coefficients for time
  //' t.}
  //' \item{zerogroups1}{Vector of zero-inflated group coefficients for time
  //' t-1.}
  //' \item{indcova2s}{Vector of individual covariate \code{a} values for time
  //' t.}
  //' \item{indcova1s}{Vector of individual covariate \code{a} values for time
  //' t-1.}
  //' \item{indcovb2s}{Vector of individual covariate \code{b} values for time
  //' t.}
  //' \item{indcovb1s}{Vector of individual covariate \code{b} values for time
  //' t-1.}
  //' \item{indcovc2s}{Vector of individual covariate \code{c} values for time
  //' t.}
  //' \item{indcovc1s}{Vector of individual covariate \code{c} values for time
  //' t-1.}
  //' \item{zeroindcova2s}{Vector of zero-inflated individual covariate \code{a}
  //' values for time t.}
  //' \item{zeroindcova1s}{Vector of zero-inflated individual covariate \code{a}
  //' values for time t-1.}
  //' \item{zeroindcovb2s}{Vector of zero-inflated individual covariate \code{b}
  //' values for time t.}
  //' \item{zeroindcovb1s}{Vector of zero-inflated individual covariate \code{b}
  //' values for time t-1.}
  //' \item{zeroindcovc2s}{Vector of zero-inflated individual covariate \code{c}
  //' values for time t.}
  //' \item{zeroindcovc1s}{Vector of zero-inflated individual covariate \code{c}
  //' values for time t-1.}
  //' \item{class}{The R class of the vital rate model.}
  //' \item{family}{The response distribution.}
  //' \item{dist}{An integer representing the response distribution. \code{0} = 
  //' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
  //' = binomial, and \code{5} = constant.}
  //' \item{zero_inflated}{A logical value indicating whether the distribution is
  //' zero-inflated.}
  //' \item{zero_truncated}{A logical value indicating whether the distribution is
  //' zero-truncated.}
  //' \item{sigma}{The residual standard deviation of the model. Defaults to
  //' \code{1.0}. Equivalent output to package lme4's \code{sigma()} function.}
  //' \item{theta}{The scale parameter theta used in the negative binomial
  //' distribution. Defaults to \code{1.0}.}
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List modelextract(RObject object, const DataFrame& paramnames,
    const CharacterVector& mainyears, const CharacterVector& mainpatches,
    RObject maingroups, RObject mainindcova, RObject mainindcovb,
    RObject mainindcovc, bool nodata = false) {
    
    CharacterVector fixed_zi_vars;
    NumericVector fixed_zi_slopes;
    List random_vars;
    List random_slopes;
    List random_zi_vars;
    List random_zi_slopes;
    
    NumericVector indcova2s;
    NumericVector indcova1s;
    NumericVector indcovb2s;
    NumericVector indcovb1s;
    NumericVector indcovc2s;
    NumericVector indcovc1s;
    NumericVector zeroindcova2s;
    NumericVector zeroindcova1s;
    NumericVector zeroindcovb2s;
    NumericVector zeroindcovb1s;
    NumericVector zeroindcovc2s;
    NumericVector zeroindcovc1s;
    
    List core_components;
    if (object.isS4()) {
      core_components = S4_extractor(as<S4>(object));
      
    } else if (is<List>(object) && !nodata) {
      core_components = S3_extractor(as<List>(object));
      
    } else if (nodata) {
      core_components = vrm_extractor(as<List>(object));
      
    } else {
      core_components = numeric_extractor(as<NumericVector>(object));
    }
    
    CharacterVector fixed_vars = core_components["fixed_vars"];
    NumericVector fixed_slopes = core_components["fixed_slopes"];
    int no_fixed_vars = static_cast<int>(fixed_vars.length());
    
    Nullable<CharacterVector> fixed_zi_vars_ = core_components["fixed_zi_vars"];
    Nullable<NumericVector> fixed_zi_slopes_ = core_components["fixed_zi_slopes"];
    
    if (fixed_zi_vars_.isNotNull()) {
      fixed_zi_vars = fixed_zi_vars_;
    } else {
      fixed_zi_vars = {"Intercept"};
    }
    if (fixed_zi_slopes_.isNotNull()) {
      fixed_zi_slopes = fixed_zi_slopes_;
    } else {
      fixed_zi_slopes = {0.0};
    }
    
    CharacterVector modelparam_names = paramnames["modelparams"];
    int pmnames_length = static_cast<int>(modelparam_names.length());
    
    std::string year2var = as<std::string>(modelparam_names(0));
    std::string individvar = as<std::string>(modelparam_names(1));
    std::string patchvar = as<std::string>(modelparam_names(2));
    std::string surv3var = as<std::string>(modelparam_names(3));
    std::string obs3var = as<std::string>(modelparam_names(4));
    std::string size3var = as<std::string>(modelparam_names(5));
    std::string sizeb3var = as<std::string>(modelparam_names(6));
    std::string sizec3var = as<std::string>(modelparam_names(7));
    std::string repst3var = as<std::string>(modelparam_names(8));
    std::string fec3var = as<std::string>(modelparam_names(9));
    std::string fec2var = as<std::string>(modelparam_names(10));
    std::string size2var = as<std::string>(modelparam_names(11));
    std::string size1var = as<std::string>(modelparam_names(12));
    std::string sizeb2var = as<std::string>(modelparam_names(13));
    std::string sizeb1var = as<std::string>(modelparam_names(14));
    std::string sizec2var = as<std::string>(modelparam_names(15));
    std::string sizec1var = as<std::string>(modelparam_names(16));
    std::string repst2var = as<std::string>(modelparam_names(17));
    std::string repst1var = as<std::string>(modelparam_names(18));
    std::string matst3var = as<std::string>(modelparam_names(19));
    std::string matst2var = as<std::string>(modelparam_names(20));
    std::string agevar = as<std::string>(modelparam_names(21));
    std::string densityvar = as<std::string>(modelparam_names(22));
    std::string indcova2var = as<std::string>(modelparam_names(23));
    std::string indcova1var = as<std::string>(modelparam_names(24));
    std::string indcovb2var = as<std::string>(modelparam_names(25));
    std::string indcovb1var = as<std::string>(modelparam_names(26));
    std::string indcovc2var = as<std::string>(modelparam_names(27));
    std::string indcovc1var = as<std::string>(modelparam_names(28));
    std::string group2var = as<std::string>(modelparam_names(29));
    std::string group1var = as<std::string>(modelparam_names(30));
    
    std::string annucova2var;
    std::string annucova1var;
    std::string annucovb2var;
    std::string annucovb1var;
    std::string annucovc2var;
    std::string annucovc1var;
    
    if (pmnames_length > 31) {
      annucova2var = as<std::string>(modelparam_names(31));
      annucova1var = as<std::string>(modelparam_names(32));
      annucovb2var = as<std::string>(modelparam_names(33));
      annucovb1var = as<std::string>(modelparam_names(34));
      annucovc2var = as<std::string>(modelparam_names(35));
      annucovc1var = as<std::string>(modelparam_names(36));
    } else {
      annucova2var = "none";
      annucova1var = "none";
      annucovb2var = "none";
      annucovb1var = "none";
      annucovc2var = "none";
      annucovc1var = "none";
    }
    
    int no_fixed_zi_slopes = fixed_zi_slopes.length();
    
    int no_years = mainyears.length();
    CharacterVector mainyears_text(mainyears);
    NumericVector year_coefs(no_years);
    NumericVector zeroyear_coefs(no_years);
    year_coefs.attr("names") = mainyears_text;
    zeroyear_coefs.attr("names") = mainyears_text;
    
    int no_patches = mainpatches.length();
    NumericVector patch_coefs(no_patches);
    NumericVector zeropatch_coefs(no_patches);
    
    if (no_patches < 2 && Rcpp::traits::is_na<STRSXP>(mainpatches(0))) {
      CharacterVector new_patch_names = {"pop"};
      patch_coefs.attr("names") = new_patch_names;
      zeropatch_coefs.attr("names") = new_patch_names;
    } else {
      patch_coefs.attr("names") = mainpatches;
      zeropatch_coefs.attr("names") = mainpatches;
    }
    
    CharacterVector maingroups_text(maingroups);
    int no_groups = maingroups_text.length();
    NumericVector group2_coefs(no_groups);
    NumericVector group1_coefs(no_groups);
    NumericVector zerogroup2_coefs(no_groups);
    NumericVector zerogroup1_coefs(no_groups);
    
    group2_coefs.attr("names") = maingroups_text;
    group1_coefs.attr("names") = maingroups_text;
    zerogroup2_coefs.attr("names") = maingroups_text;
    zerogroup1_coefs.attr("names") = maingroups_text;
    
    // Individual covariates
    CharacterVector indcova_names;
    CharacterVector indcovb_names;
    CharacterVector indcovc_names;
    int no_indcova_names {0};
    int no_indcovb_names {0};
    int no_indcovc_names {0};
    
    if (is<CharacterVector>(mainindcova) || is<NumericVector>(mainindcova)) {
      indcova_names = as<CharacterVector>(mainindcova);
      no_indcova_names = indcova_names.length();
      
      bool realnames = false;
      if (no_indcova_names > 0) {
        for (int i = 0; i < no_indcova_names; i++) {
          if (!stringcompare_hard(as<std::string>(indcova_names(i)), "none")) realnames = true;
        }
      }
      
      if (!realnames) {
        no_indcova_names = 0;
        
      } else {
        NumericVector indcova2s_inc(no_indcova_names);
        NumericVector indcova1s_inc(no_indcova_names);
        NumericVector zeroindcova2s_inc(no_indcova_names);
        NumericVector zeroindcova1s_inc(no_indcova_names);
        
        indcova2s = indcova2s_inc;
        indcova1s = indcova1s_inc;
        zeroindcova2s = zeroindcova2s_inc;
        zeroindcova1s = zeroindcova1s_inc;
        
        indcova2s.attr("names") = indcova_names;
        indcova1s.attr("names") = indcova_names;
        zeroindcova2s.attr("names") = indcova_names;
        zeroindcova1s.attr("names") = indcova_names;
      }
    }
    
    if (is<CharacterVector>(mainindcovb) || is<NumericVector>(mainindcovb)) {
      indcovb_names = as<CharacterVector>(mainindcovb);
      no_indcovb_names = indcovb_names.length();
      
      bool realnames = false;
      if (no_indcovb_names > 0) {
        for (int i = 0; i < no_indcovb_names; i++) {
          if (!stringcompare_hard(as<std::string>(indcovb_names(i)), "none")) realnames = true;
        }
      }
      
      if (!realnames) {
        no_indcovb_names = 0;
      } else {
        NumericVector indcovb2s_inc(no_indcovb_names);
        NumericVector indcovb1s_inc(no_indcovb_names);
        NumericVector zeroindcovb2s_inc(no_indcovb_names);
        NumericVector zeroindcovb1s_inc(no_indcovb_names);
        
        indcovb2s = indcovb2s_inc;
        indcovb1s = indcovb1s_inc;
        zeroindcovb2s = zeroindcovb2s_inc;
        zeroindcovb1s = zeroindcovb1s_inc;
        
        indcovb2s.attr("names") = indcovb_names;
        indcovb1s.attr("names") = indcovb_names;
        zeroindcovb2s.attr("names") = indcovb_names;
        zeroindcovb1s.attr("names") = indcovb_names;
      }
    }
    
    if (is<CharacterVector>(mainindcovc) || is<NumericVector>(mainindcovc)) {
      indcovc_names = as<CharacterVector>(mainindcovc);
      no_indcovc_names = indcovc_names.length();
      
      bool realnames = false;
      if (no_indcovc_names > 0) {
        for (int i = 0; i < no_indcovc_names; i++) {
          if (!stringcompare_hard(as<std::string>(indcovc_names(i)), "none")) realnames = true;
        }
      }
      
      if (!realnames) {
        no_indcovc_names = 0;
      } else {
        NumericVector indcovc2s_inc(no_indcovc_names);
        NumericVector indcovc1s_inc(no_indcovc_names);
        NumericVector zeroindcovc2s_inc(no_indcovc_names);
        NumericVector zeroindcovc1s_inc(no_indcovc_names);
        
        indcovc2s = indcovc2s_inc;
        indcovc1s = indcovc1s_inc;
        zeroindcovc2s = zeroindcovc2s_inc;
        zeroindcovc1s = zeroindcovc1s_inc;
        
        indcovc2s.attr("names") = indcovc_names;
        indcovc1s.attr("names") = indcovc_names;
        zeroindcovc2s.attr("names") = indcovc_names;
        zeroindcovc1s.attr("names") = indcovc_names;
      }
    }
    
    NumericVector coef_vec(514);
    
    for (int i = 0; i < no_fixed_vars; i++) {
      for (int j = 0; j < no_years; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)),
          as<std::string>(mainyears_text(j)), false)) {
          year_coefs(j) = fixed_slopes(i);
        }
      }
      
      for (int j = 0; j < no_patches; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), patchvar, false) ||
            stringcompare_simple(as<std::string>(fixed_vars(i)), "patch", false)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)),
            as<std::string>(mainpatches(j)), false)) {
            patch_coefs(j) = fixed_slopes(i);
          }
        }
      }
      
      for (int j = 0; j < no_groups; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), group2var, false)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)),
            as<std::string>(maingroups_text(j)), false)) {
            group2_coefs(j) = fixed_slopes(i);
          }
        }
      }
      for (int j = 0; j < no_groups; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), group1var, false)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)),
            as<std::string>(maingroups_text(j)), false)) {
            group1_coefs(j) = fixed_slopes(i);
          }
        }
      }
      
      if (no_indcova_names > 0) {
        for (int j = 0; j < no_indcova_names; j++) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcova2var, false)) {
            std::string new_check = stringremove(as<std::string>(fixed_vars(i)), indcova2var);
            
            if (!stringcompare_hard(as<std::string>(fixed_vars(i)), new_check) &&
                stringcompare_simple(new_check, as<std::string>(indcova_names(j)), false)) {
              indcova2s(j) = fixed_slopes(i);
            }
          }
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcova1var, false)) {
            std::string new_check = stringremove(as<std::string>(fixed_vars(i)), indcova1var);
            
            if (!stringcompare_hard(as<std::string>(fixed_vars(i)), new_check) &&
                stringcompare_simple(new_check, as<std::string>(indcova_names(j)), false)) {
              indcova1s(j) = fixed_slopes(i);
            }
          }
        }
      }
      
      if (no_indcovb_names > 0) {
        for (int j = 0; j < no_indcovb_names; j++) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovb2var, false)) {
            std::string new_check = stringremove(as<std::string>(fixed_vars(i)), indcovb2var);
            
            if (!stringcompare_hard(as<std::string>(fixed_vars(i)), new_check) &&
                stringcompare_simple(new_check, as<std::string>(indcovb_names(j)), false)) {
              indcovb2s(j) = fixed_slopes(i);
            }
          }
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovb1var, false)) {
            std::string new_check = stringremove(as<std::string>(fixed_vars(i)), indcovb1var);
            
            if (!stringcompare_hard(as<std::string>(fixed_vars(i)), new_check) &&
                stringcompare_simple(new_check, as<std::string>(indcovb_names(j)), false)) {
              indcovb1s(j) = fixed_slopes(i);
            }
          }
        }
      }
      
      if (no_indcovc_names > 0) {
        for (int j = 0; j < no_indcovc_names; j++) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovc2var, false)) {
            std::string new_check = stringremove(as<std::string>(fixed_vars(i)), indcovc2var);
            
            if (!stringcompare_hard(as<std::string>(fixed_vars(i)), new_check) &&
                stringcompare_simple(new_check, as<std::string>(indcovc_names(j)), false)) {
              indcovc2s(j) = fixed_slopes(i);
            }
          }
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovc1var, false)) {
            std::string new_check = stringremove(as<std::string>(fixed_vars(i)), indcovc1var);
            
            if (!stringcompare_hard(as<std::string>(fixed_vars(i)), new_check) &&
                stringcompare_simple(new_check, as<std::string>(indcovc_names(j)), false)) {
              indcovc1s(j) = fixed_slopes(i);
            }
          }
        }
      }
      
      if (stringcompare_simple(as<std::string>(fixed_vars(i)), "ntercep", false)) {
        coef_vec(0) = fixed_slopes(i);
      }
      
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), repst1var)) {
        coef_vec(1) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), repst2var)) {
        coef_vec(2) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), size1var)) {
        coef_vec(3) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), size2var)) {
        coef_vec(4) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), repst1var, repst2var)) {
        coef_vec(5) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, size2var)) {
        coef_vec(6) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, repst1var)) {
        coef_vec(7) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, repst2var)) {
        coef_vec(8) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, repst1var)) {
        coef_vec(9) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, repst2var)) {
        coef_vec(10) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), agevar)) {
        coef_vec(11) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, agevar)) {
        coef_vec(12) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, agevar)) {
        coef_vec(13) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), repst1var, agevar)) {
        coef_vec(14) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), repst2var, agevar)) {
        coef_vec(15) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcova2var)) {
        coef_vec(16) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovb2var)) {
        coef_vec(17) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovc2var)) {
        coef_vec(18) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcova1var)) {
        coef_vec(19) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovb1var)) {
        coef_vec(20) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovc1var)) {
        coef_vec(21) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, size2var)) {
        coef_vec(22) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, size2var)) {
        coef_vec(23) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, size2var)) {
        coef_vec(24) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, repst2var)) {
        coef_vec(25) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, repst2var)) {
        coef_vec(26) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, repst2var)) {
        coef_vec(27) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, size1var)) {
        coef_vec(28) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, size1var)) {
        coef_vec(29) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, size1var)) {
        coef_vec(30) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, repst1var)) {
        coef_vec(31) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, repst1var)) {
        coef_vec(32) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, repst1var)) {
        coef_vec(33) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovb2var)) {
        coef_vec(34) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovc2var)) {
        coef_vec(35) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, indcovc2var)) {
        coef_vec(36) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovb1var)) {
        coef_vec(37) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovc1var)) {
        coef_vec(38) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, indcovc1var)) {
        coef_vec(39) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovb1var)) {
        coef_vec(40) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovb2var)) {
        coef_vec(41) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovc1var)) {
        coef_vec(42) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovc2var)) {
        coef_vec(43) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, indcovc1var)) {
        coef_vec(44) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, indcovc2var)) {
        coef_vec(45) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, agevar)) {
        coef_vec(92) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, agevar)) {
        coef_vec(93) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, agevar)) {
        coef_vec(93) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, agevar)) {
        coef_vec(95) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, agevar)) {
        coef_vec(183) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcova1var)) {
        coef_vec(184) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, indcovb1var)) {
        coef_vec(185) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, indcovc1var)) {
        coef_vec(186) = fixed_slopes(i);
      }
      
      // New coefficients
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizeb2var)) {
        coef_vec(100) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizeb1var)) {
        coef_vec(101) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizec2var)) {
        coef_vec(102) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizec1var)) {
        coef_vec(103) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), densityvar)) {
        coef_vec(104) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, sizeb2var)) {
        coef_vec(105) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, sizec2var)) {
        coef_vec(106) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizeb1var)) {
        coef_vec(107) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizec1var)) {
        coef_vec(108) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, sizec1var)) {
        coef_vec(109) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizeb2var)) {
        coef_vec(110) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizec2var)) {
        coef_vec(111) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, sizec2var)) {
        coef_vec(112) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizeb2var)) {
        coef_vec(113) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizec2var)) {
        coef_vec(114) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, sizec2var)) {
        coef_vec(115) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizeb1var)) {
        coef_vec(116) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizec1var)) {
        coef_vec(117) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, sizec1var)) {
        coef_vec(118) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, densityvar)) {
        coef_vec(119) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, densityvar)) {
        coef_vec(120) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, densityvar)) {
        coef_vec(121) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, densityvar)) {
        coef_vec(122) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, densityvar)) {
        coef_vec(123) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, densityvar)) {
        coef_vec(124) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), repst2var, densityvar)) {
        coef_vec(125) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), repst1var, densityvar)) {
        coef_vec(126) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, repst2var)) {
        coef_vec(127) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, repst2var)) {
        coef_vec(128) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, agevar)) {
        coef_vec(129) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, repst1var)) {
        coef_vec(130) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, repst1var)) {
        coef_vec(131) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, repst2var)) {
        coef_vec(132) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, repst1var)) {
        coef_vec(133) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, repst1var)) {
        coef_vec(134) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, repst2var)) {
        coef_vec(135) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, agevar)) {
        coef_vec(136) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, agevar)) {
        coef_vec(137) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), densityvar, agevar)) {
        coef_vec(138) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, agevar)) {
        coef_vec(139) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, agevar)) {
        coef_vec(140) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizeb2var)) {
        coef_vec(141) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizec2var)) {
        coef_vec(142) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, densityvar)) {
        coef_vec(143) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizeb1var)) {
        coef_vec(144) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizec1var)) {
        coef_vec(145) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizeb2var)) {
        coef_vec(146) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizec2var)) {
        coef_vec(147) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizeb1var)) {
        coef_vec(148) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizec1var)) {
        coef_vec(149) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, densityvar)) {
        coef_vec(150) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizeb2var)) {
        coef_vec(151) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizec2var)) {
        coef_vec(152) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, densityvar)) {
        coef_vec(153) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizeb1var)) {
        coef_vec(154) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizec1var)) {
        coef_vec(155) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizeb2var)) {
        coef_vec(156) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizec2var)) {
        coef_vec(157) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizeb1var)) {
        coef_vec(158) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizec1var)) {
        coef_vec(159) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, densityvar)) {
        coef_vec(160) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizeb2var)) {
        coef_vec(161) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizec2var)) {
        coef_vec(162) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, densityvar)) {
        coef_vec(163) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizeb1var)) {
        coef_vec(164) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizec1var)) {
        coef_vec(165) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizeb2var)) {
        coef_vec(166) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizec2var)) {
        coef_vec(167) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizeb1var)) {
        coef_vec(168) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizec1var)) {
        coef_vec(169) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, densityvar)) {
        coef_vec(170) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, size1var)) {
        coef_vec(171) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, size1var)) {
        coef_vec(172) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, size1var)) {
        coef_vec(173) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, size2var)) {
        coef_vec(174) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, size2var)) {
        coef_vec(175) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, size2var)) {
        coef_vec(176) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, repst1var)) {
        coef_vec(177) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, repst1var)) {
        coef_vec(178) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, repst1var)) {
        coef_vec(179) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, repst2var)) {
        coef_vec(180) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, repst2var)) {
        coef_vec(181) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, repst2var)) {
        coef_vec(182) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, size2var)) {
        coef_vec(187) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, size2var)) {
        coef_vec(188) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, size2var)) {
        coef_vec(189) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, size2var)) {
        coef_vec(190) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, size2var)) {
        coef_vec(191) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, size2var)) {
        coef_vec(192) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, size1var)) {
        coef_vec(193) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, size1var)) {
        coef_vec(194) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, size1var)) {
        coef_vec(195) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, size1var)) {
        coef_vec(196) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, size1var)) {
        coef_vec(197) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, size1var)) {
        coef_vec(198) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, sizeb2var)) {
        coef_vec(199) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, sizeb2var)) {
        coef_vec(300) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, sizeb2var)) {
        coef_vec(301) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, sizeb2var)) {
        coef_vec(302) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, sizeb2var)) {
        coef_vec(303) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, sizeb2var)) {
        coef_vec(304) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, sizeb1var)) {
        coef_vec(305) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, sizeb1var)) {
        coef_vec(306) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, sizeb1var)) {
        coef_vec(307) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, sizeb1var)) {
        coef_vec(308) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, sizeb1var)) {
        coef_vec(309) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, sizeb1var)) {
        coef_vec(310) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, sizec2var)) {
        coef_vec(311) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, sizec2var)) {
        coef_vec(312) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, sizec2var)) {
        coef_vec(313) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, sizec2var)) {
        coef_vec(314) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, sizec2var)) {
        coef_vec(315) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, sizec2var)) {
        coef_vec(316) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, sizec1var)) {
        coef_vec(317) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, sizec1var)) {
        coef_vec(318) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, sizec1var)) {
        coef_vec(319) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, sizec1var)) {
        coef_vec(320) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, sizec1var)) {
        coef_vec(321) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, sizec1var)) {
        coef_vec(322) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, repst2var)) {
        coef_vec(323) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, repst2var)) {
        coef_vec(324) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, repst2var)) {
        coef_vec(325) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, repst2var)) {
        coef_vec(326) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, repst2var)) {
        coef_vec(327) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, repst2var)) {
        coef_vec(328) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, repst1var)) {
        coef_vec(329) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, repst1var)) {
        coef_vec(330) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, repst1var)) {
        coef_vec(331) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, repst1var)) {
        coef_vec(332) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, repst1var)) {
        coef_vec(333) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, repst1var)) {
        coef_vec(334) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, agevar)) {
        coef_vec(335) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, agevar)) {
        coef_vec(336) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, agevar)) {
        coef_vec(337) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, agevar)) {
        coef_vec(338) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, agevar)) {
        coef_vec(339) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, agevar)) {
        coef_vec(340) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, densityvar)) {
        coef_vec(341) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, densityvar)) {
        coef_vec(342) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, densityvar)) {
        coef_vec(343) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, densityvar)) {
        coef_vec(344) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, densityvar)) {
        coef_vec(345) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, densityvar)) {
        coef_vec(346) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, indcova2var)) {
        coef_vec(347) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, indcova2var)) {
        coef_vec(348) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, indcova2var)) {
        coef_vec(349) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, indcova2var)) {
        coef_vec(350) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, indcova2var)) {
        coef_vec(351) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, indcova2var)) {
        coef_vec(352) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, indcova1var)) {
        coef_vec(353) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, indcova1var)) {
        coef_vec(354) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, indcova1var)) {
        coef_vec(355) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, indcova1var)) {
        coef_vec(356) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, indcova1var)) {
        coef_vec(357) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, indcova1var)) {
        coef_vec(358) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, indcovb2var)) {
        coef_vec(359) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, indcovb2var)) {
        coef_vec(360) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, indcovb2var)) {
        coef_vec(361) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, indcovb2var)) {
        coef_vec(362) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, indcovb2var)) {
        coef_vec(363) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, indcovb2var)) {
        coef_vec(364) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, indcovb1var)) {
        coef_vec(365) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, indcovb1var)) {
        coef_vec(366) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, indcovb1var)) {
        coef_vec(367) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, indcovb1var)) {
        coef_vec(368) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, indcovb1var)) {
        coef_vec(369) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, indcovb1var)) {
        coef_vec(370) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, indcovc2var)) {
        coef_vec(371) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, indcovc2var)) {
        coef_vec(372) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, indcovc2var)) {
        coef_vec(373) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, indcovc2var)) {
        coef_vec(374) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, indcovc2var)) {
        coef_vec(375) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, indcovc2var)) {
        coef_vec(376) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, indcovc1var)) {
        coef_vec(377) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, indcovc1var)) {
        coef_vec(378) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, indcovc1var)) {
        coef_vec(379) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, indcovc1var)) {
        coef_vec(380) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, indcovc1var)) {
        coef_vec(381) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc1var, indcovc1var)) {
        coef_vec(382) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), annucova2var)) {
        coef_vec(383) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, annucova1var)) {
        coef_vec(384) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, annucovb2var)) {
        coef_vec(385) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, annucovb1var)) {
        coef_vec(386) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, annucovc2var)) {
        coef_vec(387) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova2var, annucovc1var)) {
        coef_vec(388) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), annucova1var)) {
        coef_vec(389) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, annucovb2var)) {
        coef_vec(390) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, annucovb1var)) {
        coef_vec(391) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, annucovc2var)) {
        coef_vec(392) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucova1var, annucovc1var)) {
        coef_vec(393) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), annucovb2var)) {
        coef_vec(394) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, annucovb1var)) {
        coef_vec(395) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, annucovc2var)) {
        coef_vec(396) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb2var, annucovc1var)) {
        coef_vec(397) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), annucovb1var)) {
        coef_vec(398) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, annucovc2var)) {
        coef_vec(399) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovb1var, annucovc1var)) {
        coef_vec(500) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), annucovc2var)) {
        coef_vec(501) = fixed_slopes(i);
      }
      if (stringcompare_x(as<std::string>(fixed_vars(i)), annucovc2var, annucovc1var)) {
        coef_vec(502) = fixed_slopes(i);
      }
      if (stringcompare_hard(as<std::string>(fixed_vars(i)), annucovc1var)) {
        coef_vec(503) = fixed_slopes(i);
      }
    }
    
    for (int i = 0; i < no_fixed_zi_slopes; i++) {
      for (int j = 0; j < no_years; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
          as<std::string>(mainyears_text(j)), false)) {
          zeroyear_coefs(j) = fixed_zi_slopes(i);
        }
      }
      for (int j = 0; j < no_patches; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), patchvar, false) ||
            stringcompare_simple(as<std::string>(fixed_zi_vars(i)), "patch", false)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
            as<std::string>(mainpatches(j)), false)) {
            zeropatch_coefs(j) = fixed_zi_slopes(i);
          }
        }
      }
      
      for (int j = 0; j < no_groups; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), group2var, false)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
            as<std::string>(maingroups_text(j)), false)) {
            zerogroup2_coefs(j) = fixed_zi_slopes(i);
          }
        }
      }
      for (int j = 0; j < no_groups; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), group1var, false)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
            as<std::string>(maingroups_text(j)), false)) {
            zerogroup1_coefs(j) = fixed_zi_slopes(i);
          }
        }
      }
      
      if (no_indcova_names > 0) {
        for (int j = 0; j < no_indcova_names; j++) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcova2var, false)) {
            if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
              as<std::string>(indcova_names(j)), false)) {
              zeroindcova2s(j) = fixed_zi_slopes(i);
            }
          }
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcova1var, false)) {
            if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
              as<std::string>(indcova_names(j)), false)) {
              zeroindcova1s(j) = fixed_zi_slopes(i);
            }
          }
        }
      }
      
      if (no_indcovb_names > 0) {
        for (int j = 0; j < no_indcovb_names; j++) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovb2var, false)) {
            if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
              as<std::string>(indcovb_names(j)), false)) {
              zeroindcovb2s(j) = fixed_zi_slopes(i);
            }
          }
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovb1var, false)) {
            if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
              as<std::string>(indcovb_names(j)), false)) {
              zeroindcovb1s(j) = fixed_zi_slopes(i);
            }
          }
        }
      }
      
      if (no_indcovc_names > 0) {
        for (int j = 0; j < no_indcovc_names; j++) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovc2var, false)) {
            if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
              as<std::string>(indcovc_names(j)), false)) {
              zeroindcovc2s(j) = fixed_zi_slopes(i);
            }
          }
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovc1var, false)) {
            if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)),
              as<std::string>(indcovc_names(j)), false)) {
              zeroindcovc1s(j) = fixed_zi_slopes(i);
            }
          }
        }
      }
      
      if (no_fixed_zi_slopes > 0) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), "ntercep", false)) {
          coef_vec(46) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), repst1var)) {
          coef_vec(47) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), repst2var)) {
          coef_vec(48) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), size1var)) {
          coef_vec(49) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), size2var)) {
          coef_vec(50) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst1var, repst2var)) {
          coef_vec(51) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, size2var)) {
          coef_vec(52) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, repst1var)) {
          coef_vec(53) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, repst2var)) {
          coef_vec(54) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, repst1var)) {
          coef_vec(55) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, repst2var)) {
          coef_vec(56) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), agevar)) {
          coef_vec(57) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, agevar)) {
          coef_vec(58) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, agevar)) {
          coef_vec(59) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst1var, agevar)) {
          coef_vec(60) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst2var, agevar)) {
          coef_vec(61) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcova2var)) {
          coef_vec(62) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovb2var)) {
          coef_vec(63) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovc2var)) {
          coef_vec(64) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcova1var)) {
          coef_vec(65) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovb1var)) {
          coef_vec(66) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovc1var)) {
          coef_vec(67) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, size2var)) {
          coef_vec(68) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, size2var)) {
          coef_vec(69) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, size2var)) {
          coef_vec(70) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, repst2var)) {
          coef_vec(71) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, repst2var)) {
          coef_vec(72) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, repst2var)) {
          coef_vec(73) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, size1var)) {
          coef_vec(74) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, size1var)) {
          coef_vec(75) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, size1var)) {
          coef_vec(76) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, repst1var)) {
          coef_vec(77) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, repst1var)) {
          coef_vec(78) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, repst1var)) {
          coef_vec(79) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovb2var)) {
          coef_vec(80) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovc2var)) {
          coef_vec(81) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, indcovc2var)) {
          coef_vec(82) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovb1var)) {
          coef_vec(83) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovc1var)) {
          coef_vec(84) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, indcovc1var)) {
          coef_vec(85) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovb1var)) {
          coef_vec(86) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovb2var)) {
          coef_vec(87) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovc1var)) {
          coef_vec(88) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovc2var)) {
          coef_vec(89) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, indcovc1var)) {
          coef_vec(90) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, indcovc2var)) {
          coef_vec(91) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, agevar)) {
          coef_vec(96) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, agevar)) {
          coef_vec(97) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, agevar)) {
          coef_vec(98) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, agevar)) {
          coef_vec(99) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, agevar)) {
          coef_vec(283) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcova1var)) {
          coef_vec(284) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, indcovb1var)) {
          coef_vec(285) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, indcovc1var)) {
          coef_vec(286) = fixed_zi_slopes(i);
        }
        
        // New coefficients
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizeb2var)) {
          coef_vec(200) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizeb1var)) {
          coef_vec(201) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizec2var)) {
          coef_vec(202) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizec1var)) {
          coef_vec(203) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), densityvar)) {
          coef_vec(204) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, sizeb2var)) {
          coef_vec(205) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, sizec2var)) {
          coef_vec(206) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizeb1var)) {
          coef_vec(207) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizec1var)) {
          coef_vec(208) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, sizec1var)) {
          coef_vec(209) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizeb2var)) {
          coef_vec(210) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizec2var)) {
          coef_vec(211) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, sizec2var)) {
          coef_vec(212) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizeb2var)) {
          coef_vec(213) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizec2var)) {
          coef_vec(214) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, sizec2var)) {
          coef_vec(215) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizeb1var)) {
          coef_vec(216) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizec1var)) {
          coef_vec(217) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, sizec1var)) {
          coef_vec(218) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, densityvar)) {
          coef_vec(219) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, densityvar)) {
          coef_vec(220) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, densityvar)) {
          coef_vec(221) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, densityvar)) {
          coef_vec(222) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, densityvar)) {
          coef_vec(223) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, densityvar)) {
          coef_vec(224) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst2var, densityvar)) {
          coef_vec(225) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst1var, densityvar)) {
          coef_vec(226) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, repst2var)) {
          coef_vec(227) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, repst2var)) {
          coef_vec(228) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, agevar)) {
          coef_vec(229) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, repst1var)) {
          coef_vec(230) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, repst1var)) {
          coef_vec(231) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, repst2var)) {
          coef_vec(232) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, repst1var)) {
          coef_vec(233) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, repst1var)) {
          coef_vec(234) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, repst2var)) {
          coef_vec(235) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, agevar)) {
          coef_vec(236) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, agevar)) {
          coef_vec(237) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), densityvar, agevar)) {
          coef_vec(238) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, agevar)) {
          coef_vec(239) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, agevar)) {
          coef_vec(240) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizeb2var)) {
          coef_vec(241) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizec2var)) {
          coef_vec(242) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, densityvar)) {
          coef_vec(243) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizeb1var)) {
          coef_vec(244) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizec1var)) {
          coef_vec(245) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizeb2var)) {
          coef_vec(246) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizec2var)) {
          coef_vec(247) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizeb1var)) {
          coef_vec(248) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizec1var)) {
          coef_vec(249) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, densityvar)) {
          coef_vec(250) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizeb2var)) {
          coef_vec(251) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizec2var)) {
          coef_vec(252) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, densityvar)) {
          coef_vec(253) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizeb1var)) {
          coef_vec(254) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizec1var)) {
          coef_vec(255) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizeb2var)) {
          coef_vec(256) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizec2var)) {
          coef_vec(257) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizeb1var)) {
          coef_vec(258) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizec1var)) {
          coef_vec(259) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, densityvar)) {
          coef_vec(260) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizeb2var)) {
          coef_vec(261) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizec2var)) {
          coef_vec(262) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, densityvar)) {
          coef_vec(263) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizeb1var)) {
          coef_vec(264) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizec1var)) {
          coef_vec(265) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizeb2var)) {
          coef_vec(266) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizec2var)) {
          coef_vec(267) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizeb1var)) {
          coef_vec(268) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizec1var)) {
          coef_vec(269) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, densityvar)) {
          coef_vec(270) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, size1var)) {
          coef_vec(271) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, size1var)) {
          coef_vec(272) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, size1var)) {
          coef_vec(273) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, size2var)) {
          coef_vec(274) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, size2var)) {
          coef_vec(275) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, size2var)) {
          coef_vec(276) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, repst1var)) {
          coef_vec(277) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, repst1var)) {
          coef_vec(278) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, repst1var)) {
          coef_vec(279) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, repst2var)) {
          coef_vec(280) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, repst2var)) {
          coef_vec(281) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, repst2var)) {
          coef_vec(282) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, size2var)) {
          coef_vec(287) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, size2var)) {
          coef_vec(288) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, size2var)) {
          coef_vec(289) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, size2var)) {
          coef_vec(290) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, size2var)) {
          coef_vec(291) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, size2var)) {
          coef_vec(292) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, size1var)) {
          coef_vec(293) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, size1var)) {
          coef_vec(294) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, size1var)) {
          coef_vec(295) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, size1var)) {
          coef_vec(296) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, size1var)) {
          coef_vec(297) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, size1var)) {
          coef_vec(298) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, sizeb2var)) {
          coef_vec(299) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, sizeb2var)) {
          coef_vec(400) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, sizeb2var)) {
          coef_vec(401) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, sizeb2var)) {
          coef_vec(402) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, sizeb2var)) {
          coef_vec(403) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, sizeb2var)) {
          coef_vec(404) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, sizeb1var)) {
          coef_vec(405) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, sizeb1var)) {
          coef_vec(406) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, sizeb1var)) {
          coef_vec(407) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, sizeb1var)) {
          coef_vec(408) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, sizeb1var)) {
          coef_vec(409) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, sizeb1var)) {
          coef_vec(410) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, sizec2var)) {
          coef_vec(411) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, sizec2var)) {
          coef_vec(412) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, sizec2var)) {
          coef_vec(413) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, sizec2var)) {
          coef_vec(414) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, sizec2var)) {
          coef_vec(415) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, sizec2var)) {
          coef_vec(416) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, sizec1var)) {
          coef_vec(417) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, sizec1var)) {
          coef_vec(418) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, sizec1var)) {
          coef_vec(419) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, sizec1var)) {
          coef_vec(420) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, sizec1var)) {
          coef_vec(421) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, sizec1var)) {
          coef_vec(422) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, repst2var)) {
          coef_vec(423) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, repst2var)) {
          coef_vec(424) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, repst2var)) {
          coef_vec(425) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, repst2var)) {
          coef_vec(426) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, repst2var)) {
          coef_vec(427) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, repst2var)) {
          coef_vec(428) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, repst1var)) {
          coef_vec(429) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, repst1var)) {
          coef_vec(430) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, repst1var)) {
          coef_vec(431) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, repst1var)) {
          coef_vec(432) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, repst1var)) {
          coef_vec(433) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, repst1var)) {
          coef_vec(434) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, agevar)) {
          coef_vec(435) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, agevar)) {
          coef_vec(436) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, agevar)) {
          coef_vec(437) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, agevar)) {
          coef_vec(438) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, agevar)) {
          coef_vec(439) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, agevar)) {
          coef_vec(440) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, densityvar)) {
          coef_vec(441) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, densityvar)) {
          coef_vec(442) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, densityvar)) {
          coef_vec(443) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, densityvar)) {
          coef_vec(444) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, densityvar)) {
          coef_vec(445) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, densityvar)) {
          coef_vec(446) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, indcova2var)) {
          coef_vec(447) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, indcova2var)) {
          coef_vec(448) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, indcova2var)) {
          coef_vec(449) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, indcova2var)) {
          coef_vec(450) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, indcova2var)) {
          coef_vec(451) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, indcova2var)) {
          coef_vec(452) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, indcova1var)) {
          coef_vec(453) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, indcova1var)) {
          coef_vec(454) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, indcova1var)) {
          coef_vec(455) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, indcova1var)) {
          coef_vec(456) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, indcova1var)) {
          coef_vec(457) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, indcova1var)) {
          coef_vec(458) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, indcovb2var)) {
          coef_vec(459) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, indcovb2var)) {
          coef_vec(460) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, indcovb2var)) {
          coef_vec(461) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, indcovb2var)) {
          coef_vec(462) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, indcovb2var)) {
          coef_vec(463) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, indcovb2var)) {
          coef_vec(464) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, indcovb1var)) {
          coef_vec(465) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, indcovb1var)) {
          coef_vec(466) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, indcovb1var)) {
          coef_vec(467) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, indcovb1var)) {
          coef_vec(468) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, indcovb1var)) {
          coef_vec(469) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, indcovb1var)) {
          coef_vec(470) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, indcovc2var)) {
          coef_vec(471) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, indcovc2var)) {
          coef_vec(472) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, indcovc2var)) {
          coef_vec(473) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, indcovc2var)) {
          coef_vec(474) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, indcovc2var)) {
          coef_vec(475) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, indcovc2var)) {
          coef_vec(476) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, indcovc1var)) {
          coef_vec(477) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, indcovc1var)) {
          coef_vec(478) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, indcovc1var)) {
          coef_vec(479) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, indcovc1var)) {
          coef_vec(480) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, indcovc1var)) {
          coef_vec(481) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc1var, indcovc1var)) {
          coef_vec(482) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), annucova2var)) {
          coef_vec(483) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, annucova1var)) {
          coef_vec(484) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, annucovb2var)) {
          coef_vec(485) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, annucovb1var)) {
          coef_vec(486) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, annucovc2var)) {
          coef_vec(487) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova2var, annucovc1var)) {
          coef_vec(488) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), annucova1var)) {
          coef_vec(489) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, annucovb2var)) {
          coef_vec(490) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, annucovb1var)) {
          coef_vec(491) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, annucovc2var)) {
          coef_vec(492) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucova1var, annucovc1var)) {
          coef_vec(493) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), annucovb2var)) {
          coef_vec(494) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, annucovb1var)) {
          coef_vec(495) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, annucovc2var)) {
          coef_vec(496) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb2var, annucovc1var)) {
          coef_vec(497) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), annucovb1var)) {
          coef_vec(498) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, annucovc2var)) {
          coef_vec(499) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovb1var, annucovc1var)) {
          coef_vec(510) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), annucovc2var)) {
          coef_vec(511) = fixed_zi_slopes(i);
        }
        if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), annucovc2var, annucovc1var)) {
          coef_vec(512) = fixed_zi_slopes(i);
        }
        if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), annucovc1var)) {
          coef_vec(513) = fixed_zi_slopes(i);
        }
      }
    }
    
    LogicalVector coef_vec_nacheck_long = is_na(coef_vec);
    LogicalVector coef_vec_nacheck = any(coef_vec_nacheck_long);
    if (coef_vec_nacheck(0) != 0) {
      Rf_warningcall(R_NilValue, "Some model coefficients are NA values. Matrices may contain NAs.");
    }
    
    // Random slopes
    Nullable<List> random_vars_ = core_components["random_vars"];
    Nullable<List> random_zi_vars_ = core_components["random_zi_vars"];
    Nullable<List> random_slopes_ = core_components["random_slopes"];
    Nullable<List> random_zi_slopes_ = core_components["random_zi_slopes"];
    
    if (random_vars_.isNotNull()) {
     
      random_vars = random_vars_;
      random_slopes = random_slopes_;
      
      CharacterVector random_names = random_vars.attr("names");
      int no_random_vars = random_slopes.length();
      
      for (int i = 0; i < no_random_vars; i++) {
        if (stringcompare_hard(as<std::string>(random_names(i)), year2var)) {
          CharacterVector ran_year_names = random_vars[i];
          NumericVector ran_year_slopes = random_slopes[i];
          int no_ran_year_slopes = ran_year_names.length();
          
          for (int j = 0; j < no_ran_year_slopes; j++) {
            for (int k = 0; k < no_years; k++) {
              if (stringcompare_hard(as<std::string>(ran_year_names(j)),
                as<std::string>(mainyears_text(k)))) {
                year_coefs(k) = ran_year_slopes(j);
              }
            }
          }
        }
        
        if (stringcompare_hard(as<std::string>(random_names(i)), patchvar) ||
            stringcompare_hard(as<std::string>(random_names(i)), "patch")) {
          CharacterVector ran_patch_names = random_vars[i];
          NumericVector ran_patch_slopes = random_slopes[i];
          int no_ran_patch_slopes = ran_patch_names.length();
          
          for (int j = 0; j < no_ran_patch_slopes; j++) {
            for (int k = 0; k < no_patches; k++) {
              if (stringcompare_hard(as<std::string>(ran_patch_names(j)),
                as<std::string>(mainpatches(k)))) {
                patch_coefs(k) = ran_patch_slopes(j);
              }
            }
          }
        }
        
        if (stringcompare_hard(as<std::string>(random_names(i)), indcova2var)) {
          CharacterVector ran_inda2_names = random_vars[i];
          NumericVector ran_inda2_slopes = random_slopes[i];
          int no_ran_inda2_slopes = ran_inda2_names.length();
        
          for (int j = 0; j < no_ran_inda2_slopes; j++) {
            for (int k = 0; k < no_indcova_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_inda2_names(j)),
                as<std::string>(indcova_names(k)))) {
                indcova2s(k) = ran_inda2_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_names(i)), indcova1var)) {
          CharacterVector ran_inda1_names = random_vars[i];
          NumericVector ran_inda1_slopes = random_slopes[i];
          int no_ran_inda1_slopes = ran_inda1_names.length();
        
          for (int j = 0; j < no_ran_inda1_slopes; j++) {
            for (int k = 0; k < no_indcova_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_inda1_names(j)),
                as<std::string>(indcova_names(k)))) {
                indcova1s(k) = ran_inda1_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_names(i)), indcovb2var)) {
          CharacterVector ran_indb2_names = random_vars[i];
          NumericVector ran_indb2_slopes = random_slopes[i];
          int no_ran_indb2_slopes = ran_indb2_names.length();
        
          for (int j = 0; j < no_ran_indb2_slopes; j++) {
            for (int k = 0; k < no_indcovb_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_indb2_names(j)),
                as<std::string>(indcovb_names(k)))) {
                indcovb2s(k) = ran_indb2_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_names(i)), indcovb1var)) {
          CharacterVector ran_indb1_names = random_vars[i];
          NumericVector ran_indb1_slopes = random_slopes[i];
          int no_ran_indb1_slopes = ran_indb1_names.length();
        
          for (int j = 0; j < no_ran_indb1_slopes; j++) {
            for (int k = 0; k < no_indcovb_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_indb1_names(j)),
                as<std::string>(indcovb_names(k)))) {
                indcovb1s(k) = ran_indb1_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_names(i)), indcovc2var)) {
          CharacterVector ran_indc2_names = random_vars[i];
          NumericVector ran_indc2_slopes = random_slopes[i];
          int no_ran_indc2_slopes = ran_indc2_names.length();
        
          for (int j = 0; j < no_ran_indc2_slopes; j++) {
            for (int k = 0; k < no_indcovc_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_indc2_names(j)),
                as<std::string>(indcovc_names(k)))) {
                indcovc2s(k) = ran_indc2_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_names(i)), indcovc1var)) {
          CharacterVector ran_indc1_names = random_vars[i];
          NumericVector ran_indc1_slopes = random_slopes[i];
          int no_ran_indc1_slopes = ran_indc1_names.length();
        
          for (int j = 0; j < no_ran_indc1_slopes; j++) {
            for (int k = 0; k < no_indcovc_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_indc1_names(j)),
                as<std::string>(indcovc_names(k)))) {
                indcovc1s(k) = ran_indc1_slopes(j);
              }
            }
          }
        }
      }
    }
    
    if (random_zi_vars_.isNotNull() && random_zi_slopes_.isNotNull()) {
      random_zi_vars = random_zi_vars_;
      random_zi_slopes = random_zi_slopes_;
      
      if (random_zi_slopes.length() > 0) {
        CharacterVector random_zi_names = random_zi_vars.attr("names");
        int no_random_zi_vars = random_zi_slopes.length();
        
        for (int i = 0; i < no_random_zi_vars; i++) {
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), year2var)) {
            CharacterVector ran_zi_year_names = random_zi_vars[i];
            NumericVector ran_zi_year_slopes = random_zi_slopes[i];
            int no_ran_zi_year_slopes = ran_zi_year_slopes.length();
            
            for (int j = 0; j < no_ran_zi_year_slopes; j++) {
              for (int k = 0; k < no_years; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_year_names(j)), 
                    as<std::string>(mainyears_text(k)))) {
                  zeroyear_coefs(k) = ran_zi_year_slopes(j);
                }
              }
            }
          }
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), patchvar) ||
              stringcompare_hard(as<std::string>(random_zi_names(i)), "patch")) {
            CharacterVector ran_zi_patch_names = random_zi_vars[i];
            NumericVector ran_zi_patch_slopes = random_zi_slopes[i];
            int no_ran_zi_patch_slopes = ran_zi_patch_slopes.length();
            
            for (int j = 0; j < no_ran_zi_patch_slopes; j++) {
              for (int k = 0; k < no_patches; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_patch_names(j)),
                    as<std::string>(mainpatches(k)))) {
                  zeropatch_coefs(k) = ran_zi_patch_slopes(j);
                }
              }
            }
          }
          
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcova2var)) {
            CharacterVector ran_zi_inda2_names = random_zi_vars[i];
            NumericVector ran_zi_inda2_slopes = random_zi_slopes[i];
            int no_ran_zi_inda2_slopes = ran_zi_inda2_slopes.length();
          
            for (int j = 0; j < no_ran_zi_inda2_slopes; j++) {
              for (int k = 0; k < no_indcova_names; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_inda2_names(j)),
                  as<std::string>(indcova_names(k)))) {
                  zeroindcova2s(k) = ran_zi_inda2_slopes(j);
                }
              }
            }
          }
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcova1var)) {
            CharacterVector ran_zi_inda1_names = random_zi_vars[i];
            NumericVector ran_zi_inda1_slopes = random_zi_slopes[i];
            int no_ran_zi_inda1_slopes = ran_zi_inda1_slopes.length();
          
            for (int j = 0; j < no_ran_zi_inda1_slopes; j++) {
              for (int k = 0; k < no_indcova_names; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_inda1_names(j)),
                  as<std::string>(indcova_names(k)))) {
                  zeroindcova1s(k) = ran_zi_inda1_slopes(j);
                }
              }
            }
          }
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovb2var)) {
            CharacterVector ran_zi_indb2_names = random_zi_vars[i];
            NumericVector ran_zi_indb2_slopes = random_zi_slopes[i];
            int no_ran_zi_indb2_slopes = ran_zi_indb2_slopes.length();
          
            for (int j = 0; j < no_ran_zi_indb2_slopes; j++) {
              for (int k = 0; k < no_indcovb_names; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_indb2_names(j)),
                  as<std::string>(indcovb_names(k)))) {
                  zeroindcovb2s(k) = ran_zi_indb2_slopes(j);
                }
              }
            }
          }
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovb1var)) {
            CharacterVector ran_zi_indb1_names = random_zi_vars[i];
            NumericVector ran_zi_indb1_slopes = random_zi_slopes[i];
            int no_ran_zi_indb1_slopes = ran_zi_indb1_slopes.length();
          
            for (int j = 0; j < no_ran_zi_indb1_slopes; j++) {
              for (int k = 0; k < no_indcovb_names; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_indb1_names(j)),
                  as<std::string>(indcovb_names(k)))) {
                  zeroindcovb1s(k) = ran_zi_indb1_slopes(j);
                }
              }
            }
          }
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovc2var)) {
            CharacterVector ran_zi_indc2_names = random_zi_vars[i];
            NumericVector ran_zi_indc2_slopes = random_zi_slopes[i];
            int no_ran_zi_indc2_slopes = ran_zi_indc2_slopes.length();
          
            for (int j = 0; j < no_ran_zi_indc2_slopes; j++) {
              for (int k = 0; k < no_indcovc_names; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_indc2_names(j)),
                  as<std::string>(indcovc_names(k)))) {
                  zeroindcovc2s(k) = ran_zi_indc2_slopes(j);
                }
              }
            }
          }
          if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovc1var)) {
            CharacterVector ran_zi_indc1_names = random_zi_vars[i];
            NumericVector ran_zi_indc1_slopes = random_zi_slopes[i];
            int no_ran_zi_indc1_slopes = ran_zi_indc1_slopes.length();
          
            for (int j = 0; j < no_ran_zi_indc1_slopes; j++) {
              for (int k = 0; k < no_indcovc_names; k++) {
                if (stringcompare_hard(as<std::string>(ran_zi_indc1_names(j)),
                  as<std::string>(indcovc_names(k)))) {
                  zeroindcovc1s(k) = ran_zi_indc1_slopes(j);
                }
              }
            }
          }
        }
      }
    }
    
    CharacterVector noneslot {"none"};
    if (no_indcova_names == 0) {
      indcova2s = {0.};
      indcova1s = {0.};
      indcova2s.attr("names") = noneslot;
      indcova1s.attr("names") = noneslot;
      
      zeroindcova2s = {0.};
      zeroindcova1s = {0.};
      zeroindcova2s.attr("names") = noneslot;
      zeroindcova1s.attr("names") = noneslot;
    }
    if (no_indcovb_names == 0) {
      indcovb2s = {0.};
      indcovb1s = {0.};
      indcovb2s.attr("names") = noneslot;
      indcovb1s.attr("names") = noneslot;
      
      zeroindcovb2s = {0.};
      zeroindcovb1s = {0.};
      zeroindcovb2s.attr("names") = noneslot;
      zeroindcovb1s.attr("names") = noneslot;
    }
    if (no_indcovc_names == 0) {
      indcovc2s = {0.};
      indcovc1s = {0.};
      indcovc2s.attr("names") = noneslot;
      indcovc1s.attr("names") = noneslot;
      
      zeroindcovc2s = {0.};
      zeroindcovc1s = {0.};
      zeroindcovc2s.attr("names") = noneslot;
      zeroindcovc1s.attr("names") = noneslot;
    }
    
    List output(28);
    
    output(0) = coef_vec;
    output(1) = year_coefs;
    output(2) = zeroyear_coefs;
    output(3) = patch_coefs;
    output(4) = zeropatch_coefs;
    output(5) = group2_coefs;
    output(6) = group1_coefs;
    output(7) = zerogroup2_coefs;
    output(8) = zerogroup1_coefs;
    output(9) = indcova2s;
    output(10) = indcova1s;
    output(11) = indcovb2s;
    output(12) = indcovb1s;
    output(13) = indcovc2s;
    output(14) = indcovc1s;
    output(15) = zeroindcova2s;
    output(16) = zeroindcova1s;
    output(17) = zeroindcovb2s;
    output(18) = zeroindcovb1s;
    output(19) = zeroindcovc2s;
    output(20) = zeroindcovc1s;
    output(21) = as<CharacterVector>(core_components["class"]);
    output(22) = as<CharacterVector>(core_components["family"]);
    output(23) = as<IntegerVector>(core_components["dist"]);
    output(24) = as<LogicalVector>(core_components["zero_inflated"]);
    output(25) = as<LogicalVector>(core_components["zero_truncated"]);
    output(26) = as<NumericVector>(core_components["sigma"]);
    output(27) = as<NumericVector>(core_components["theta"]);
  
    CharacterVector output_names = {"coefficients", "years", "zeroyear",
      "patches", "zeropatch", "groups2", "groups1", "zerogroups2",
      "zerogroups1", "indcova2s", "indcova1s", "indcovb2s", "indcovb1s",
      "indcovc2s", "indcovc1s", "zeroindcova2s", "zeroindcova1s",
      "zeroindcovb2s", "zeroindcovb1s", "zeroindcovc2s", "zeroindcovc1s",
      "class", "family", "dist", "zero_inflated", "zero_truncated", "sigma",
      "theta"};
    output.attr("names") = output_names;
    
    return output;
  }
  
  //' Estimate Value for Vital Rate Based on Inputs
  //' 
  //' Function \code{preouterator()} calculates the value of the vital rate called
  //' for by the function \code{jerzeibalowski()}.
  //' 
  //' @name preouterator
  //' 
  //' @param modelproxy A model_proxy object derived from function
  //' \code{modelextract()}.
  //' @param maincoefs The coefficients portion of the vital rate model proxy.
  //' @param randindex An integer matrix indexing all random covariates for all
  //' vital rates.
  //' @param dev_terms A numeric vector containing the deviations to the linear
  //' models input by the user. The order is: survival, observation status, size,
  //' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
  //' observation status, juvenile size, juvenile size_b, juvenile size_c,
  //' and juvenile reproductive status.
  //' @param vitalyear A matrix with year coefficients for all vital rates.
  //' @param vitalpatch A matrix with patch coefficients for all vital rates.
  //' @param chosen_r2inda A string identifying random covariate a in time t.
  //' @param chosen_r1inda A string identifying random covariate a in time t-1.
  //' @param chosen_r2indb A string identifying random covariate b in time t.
  //' @param chosen_r1indb A string identifying random covariate b in time t-1.
  //' @param chosen_r2indc A string identifying random covariate c in time t.
  //' @param chosen_r1indc A string identifying random covariate c in time t-1.
  //' @param chosen_f2inda_cat A string identifying fixed factor a in time t.
  //' @param chosen_f1inda_cat A string identifying fixed factor a in time t-1.
  //' @param chosen_f2indb_cat A string identifying fixed factor b in time t.
  //' @param chosen_f1indb_cat A string identifying fixed factor b in time t-1.
  //' @param chosen_f2indc_cat A string identifying fixed factor c in time t.
  //' @param chosen_f1indc_cat A string identifying fixed factor c in time t-1.
  //' @param status_terms A NumericVector containing, in order: fl1_i, fl2n_i,
  //' sz1_i, sz2o_i, szb1_i, szb2o_i, szc1_i, szc2o_i, aage2_i, inda_1, inda_2,
  //' indb_1, indb_2, indc_1, indc_2, used_dens, sz3_i, szb3_i, szc3_i,
  //' binwidth3_i, binbwidth3_i, bincwidth3_i, anna_2, anna_1, annb_2, annb_1,
  //' annc_2, and annc_1.
  //' @param modelgroups2 A vector of group slope coefficients for time t.
  //' @param modelgroups1 A vector of group slope coefficients for time t-1.
  //' @param modelgroups2zi A vector of zero-inflation model group slope
  //' coefficients for time t.
  //' @param modelgroups1zi A vector of zero-inflation model group slope
  //' coefficients for time t-1.
  //' @param modelyearzi A vector of zero-inflation model time slope coefficients.
  //' @param modelpatchzi A vector of zero-inflation model patch slope coefficients.
  //' @param modelind A vector of individual covariate slope coefficients.
  //' @param modelind_rownames A string vector with the names of the individual
  //' covariate coefficients.
  //' @param modelindzi A vector of individual covariate slope coefficients.
  //' @param modelind_rownames_zi A string vector with the names of the individual
  //' covariate coefficients.
  //' @param zi A logical value indicating whether model coefficients refer to the
  //' zero inflation portion of a model.
  //' @param sigma The sigma term in the \code{modelproxy} object.
  //' @param grp2o_i Stage group number in time \emph{t}.
  //' @param grp1_i Stage group number in time \emph{t}-1.
  //' @param patchnumber An integer index for pop-patch.
  //' @param yearnumber An integer index for monitoring occasion in time \emph{t}.
  //' @param vitaldist A parameter specifying the distribution of the vital rate.
  //' Current options are: Poisson (0), negative binomial (1), Gaussian (2),
  //' Gamma (3), and binomial (4).
  //' @param vitalrate An integer specifying the vital rate. 1 = surv, 2 = obs,
  //' 3 = size, 4 = sizeb, 5 = sizec, 6 = repst, 7 = fec, 8 = jsurv, 9 = jobs,
  //' 10 = jsize, 11 = jsizeb, 12 = jsizec, 13 = jrepst, 14 = jmatst.
  //' @param exp_tol A numeric value indicating the maximum limit for the
  //' \code{exp()} function to be used in vital rate calculations. Defaults to
  //' \code{700.0}.
  //' @param theta_tol A numeric value indicating a maximum value for theta in
  //' negative binomial probability density estimation. Defaults to
  //' \code{100000000.0}.
  //' @param ipm_cdf A logical value indicating whether to use the cumulative
  //' density function to estimate size transitions in continuous distributions
  //' (\code{true}), or the midpoint method (\code{false}).
  //' @param matrixformat An integer representing the style of matrix to develop.
  //' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
  //' and age-by-stage MPM (4).
  //' @param fecmod A scalar multiplier for fecundity.
  //' @param repentry_i Rep entry value for time t+1.
  //' @param negfec A logical value denoting whether to change negative estimated
  //' fecundity to 0.
  //' @param stage2n_i Numeric index of stage in time t.
  //' @param nostages The total number of stages in the stageframe.
  //' @param modeltrunc An integer coding for zero-truncation status.
  //' 
  //' @return A class double numeric value for the vital rate being estimated.
  //' 
  //' @keywords internal
  //' @noRd
  inline double preouterator(List modelproxy, NumericVector maincoefs, arma::imat randindex,
    NumericVector dev_terms, NumericMatrix vitalyear, NumericMatrix vitalpatch,
    String chosen_r2inda, String chosen_r1inda, String chosen_r2indb,
    String chosen_r1indb, String chosen_r2indc, String chosen_r1indc,
    String chosen_f2inda_cat, String chosen_f1inda_cat, String chosen_f2indb_cat,
    String chosen_f1indb_cat, String chosen_f2indc_cat, String chosen_f1indc_cat,
    NumericVector status_terms, NumericVector modelgroups2,
    NumericVector modelgroups1, NumericVector modelgroups2zi,
    NumericVector modelgroups1zi, NumericVector modelyearzi,
    NumericVector modelpatchzi, NumericVector modelind,
    StringVector modelind_rownames, NumericVector modelindzi,
    StringVector modelind_rownames_zi, bool zi, double sigma, double grp2o_i,
    double grp1_i, int patchnumber, int yearnumber, int vitaldist, int vitalrate,
    double exp_tol, double theta_tol, bool ipm_cdf, int matrixformat,
    double fecmod, double repentry_i, bool negfec, double stage2n_i, int nostages,
    int modeltrunc) {
    
    double preout {0.0};
    double all_out {0.0};
    double all_out_zi {0.0};
    
    int placeholder = vitalrate - 1;
    int placeholder_zi = placeholder + 12;
    int vitaltype {0}; // Binomial vital rates
    if (vitalrate == 3 || vitalrate == 4 || vitalrate == 5) {
      vitaltype = 1; // Size
    } else if (vitalrate == 10 || vitalrate == 11 || vitalrate == 12) {
      vitaltype = 1; // Juv size
      placeholder_zi = placeholder + 9;
    } else if (vitalrate == 7) {
      vitaltype = 2; // Fecundity
      placeholder_zi = placeholder + 11;
    }
    
    // For all / conditional models
    double mainsum = rimeotam(maincoefs, status_terms(0), status_terms(1),
      status_terms(2), status_terms(3), status_terms(4), status_terms(5),
      status_terms(6), status_terms(7), status_terms(8), status_terms(9),
      status_terms(10), status_terms(11), status_terms(12), status_terms(13),
      status_terms(14), status_terms(22), status_terms(23), status_terms(24),
      status_terms(25), status_terms(26), status_terms(27), status_terms(15),
      false);
    
    bool zi_processing = false;
    
    if (vitaltype == 1) {
      if (vitalrate == 3 || vitalrate == 10) {
        if (zi) zi_processing = true;
      } else if (vitalrate == 4 || vitalrate == 11) {
        if (zi) zi_processing = true;
      } else if (vitalrate == 5 || vitalrate == 12) {
        if (zi) zi_processing = true;
      } 
    } else if (vitaltype == 2) {
      if (zi && vitaldist < 2) zi_processing = true;  
    }
    
    // Creates covariate numerics for all models
    // Random covariate processing
    double chosen_randcova2 {0.0};
    if (chosen_r2inda != "none") {
      for (int indcount = 0; indcount < randindex(0, placeholder); indcount++) {
        if (chosen_r2inda == modelind_rownames(indcount)) {
          chosen_randcova2 = modelind(indcount);
        }
      }
    }
    double chosen_randcova1 {0.0};
    if (chosen_r1inda != "none") {
      int delectable_sum = randindex(0, placeholder);
      for (int indcount = 0; indcount < randindex(1, placeholder); indcount++) {
        if (chosen_r1inda == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcova1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb2 {0.0};
    if (chosen_r2indb != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder);
      for (int indcount = 0; indcount < randindex(2, placeholder); indcount++) {
        if (chosen_r2indb == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovb2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb1 {0.0};
    if (chosen_r1indb != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder);
      for (int indcount = 0; indcount < randindex(3, placeholder); indcount++) {
        if (chosen_r1indb == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovb1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc2 {0.0};
    if (chosen_r2indc != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder);
      for (int indcount = 0; indcount < randindex(4, placeholder); indcount++) {
        if (chosen_r2indc == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovc2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc1 {0.0};
    if (chosen_r1indc != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder) + randindex(4, placeholder);
      for (int indcount = 0; indcount < randindex(5, placeholder); indcount++) {
        if (chosen_r1indc == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovc1 = modelind(indcount + delectable_sum);
        }
      }
    }
    
    // Fixed factor covariate processing (all / conditional)
    double chosen_fixcova2 {0.0};
    if (chosen_f2inda_cat != "none") {
      for (int indcount = 0; indcount < randindex(0, placeholder); indcount++) {
        if (chosen_f2inda_cat == modelind_rownames(indcount)) {
          chosen_fixcova2 = modelind(indcount);
        }
      }
    }
    double chosen_fixcova1 {0.0};
    if (chosen_f1inda_cat != "none") {
      int delectable_sum = randindex(0, placeholder);
      for (int indcount = 0; indcount < randindex(1, placeholder); indcount++) {
        if (chosen_f1inda_cat == modelind_rownames(indcount + delectable_sum)) {
          chosen_fixcova1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_fixcovb2 {0.0};
    if (chosen_f2indb_cat != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder);
      for (int indcount = 0; indcount < randindex(2, placeholder); indcount++) {
        if (chosen_f2indb_cat == modelind_rownames(indcount + delectable_sum)) {
          chosen_fixcovb2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_fixcovb1 {0.0};
    if (chosen_f1indb_cat != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder);
      for (int indcount = 0; indcount < randindex(3, placeholder); indcount++) {
        if (chosen_f1indb_cat == modelind_rownames(indcount + delectable_sum)) {
          chosen_fixcovb1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_fixcovc2 {0.0};
    if (chosen_f2indc_cat != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder);
      for (int indcount = 0; indcount < randindex(4, placeholder); indcount++) {
        if (chosen_f2indc_cat == modelind_rownames(indcount + delectable_sum)) {
          chosen_fixcovc2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_fixcovc1 {0.0};
    if (chosen_f1indc_cat != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder) + randindex(4, placeholder);
      for (int indcount = 0; indcount < randindex(5, placeholder); indcount++) {
        if (chosen_f1indc_cat == modelind_rownames(indcount + delectable_sum)) {
          chosen_fixcovc1 = modelind(indcount + delectable_sum);
        }
      }
    }
    
    preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
      chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
      chosen_randcovc1 + chosen_fixcova2 + chosen_fixcova1 + chosen_fixcovb2 +
      chosen_fixcovb1 + chosen_fixcovc2 + chosen_fixcovc1 +
      modelgroups2(grp2o_i) + modelgroups1(grp1_i) + 
      vitalpatch(patchnumber, placeholder) +
      vitalyear(yearnumber, placeholder) + dev_terms(placeholder));
      
    if (preout > exp_tol && vitaldist < 2) preout = exp_tol;
    
    
    double preout_zi {0.0};
    
    if (zi_processing) {
      double mainsum_zi = rimeotam(maincoefs, status_terms(0), status_terms(1),
        status_terms(2), status_terms(3), status_terms(4), status_terms(5),
        status_terms(6), status_terms(7), status_terms(8), status_terms(9),
        status_terms(10), status_terms(11), status_terms(12), status_terms(13),
        status_terms(14), status_terms(22), status_terms(23), status_terms(24),
        status_terms(25), status_terms(26), status_terms(27), status_terms(15),
        true);
        
      // Only for size and fec
      double chosen_randcova2zi {0.0};
      if (chosen_r2inda != "none") {
        for (int indcount = 0; indcount < randindex(0, placeholder_zi); indcount++) {
          if (chosen_r2inda == modelind_rownames_zi(indcount)) {
            chosen_randcova2zi = modelindzi(indcount);
          }
        }
      }
      double chosen_randcova1zi {0.0};
      if (chosen_r1inda != "none") {
        int delectable_sum = randindex(0, placeholder_zi);
        for (int indcount = 0; indcount < randindex(1, placeholder_zi); indcount++) {
          if (chosen_r1inda == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_randcova1zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovb2zi {0.0};
      if (chosen_r2indb != "none") {
        int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi);
        for (int indcount = 0; indcount < randindex(2, placeholder_zi); indcount++) {
          if (chosen_r2indb == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_randcovb2zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovb1zi {0.0};
      if (chosen_r1indb != "none") {
        int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
          randindex(2, placeholder_zi);
        for (int indcount = 0; indcount < randindex(3, placeholder_zi); indcount++) {
          if (chosen_r1indb == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_randcovb1zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovc2zi {0.0};
      if (chosen_r2indc != "none") {
        int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
          randindex(2, placeholder_zi) + randindex(3, placeholder_zi);
        for (int indcount = 0; indcount < randindex(4, placeholder_zi); indcount++) {
          if (chosen_r2indc == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_randcovc2zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovc1zi {0.0};
      if (chosen_r1indc != "none") {
        int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
          randindex(2, placeholder_zi) + randindex(3, placeholder_zi) + randindex(4, placeholder_zi);
        for (int indcount = 0; indcount < randindex(5, placeholder_zi); indcount++) {
          if (chosen_r1indc == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_randcovc1zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      
      double chosen_fixcova2zi {0.0};
      if (chosen_f2inda_cat != "none") {
        for (int indcount = 0; indcount < randindex(0, 16); indcount++) {
          if (chosen_f2inda_cat == modelind_rownames_zi(indcount)) {
            chosen_fixcova2zi = modelindzi(indcount);
          }
        }
      }
      double chosen_fixcova1zi {0.0};
      if (chosen_f1inda_cat != "none") {
        int delectable_sum = randindex(0, 16);
        for (int indcount = 0; indcount < randindex(1, 16); indcount++) {
          if (chosen_f1inda_cat == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_fixcova1zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_fixcovb2zi {0.0};
      if (chosen_f2indb_cat != "none") {
        int delectable_sum = randindex(0, 16) + randindex(1, 16);
        for (int indcount = 0; indcount < randindex(2, 16); indcount++) {
          if (chosen_f2indb_cat == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_fixcovb2zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_fixcovb1zi {0.0};
      if (chosen_f1indb_cat != "none") {
        int delectable_sum = randindex(0, 16) + randindex(1, 16) + randindex(2, 16);
        for (int indcount = 0; indcount < randindex(3, 16); indcount++) {
          if (chosen_f1indb_cat == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_fixcovb1zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_fixcovc2zi {0.0};
      if (chosen_f2indc_cat != "none") {
        int delectable_sum = randindex(0, 16) + randindex(1, 16) + randindex(2, 16) +
          randindex(3, 16);
        for (int indcount = 0; indcount < randindex(4, 16); indcount++) {
          if (chosen_f2indc_cat == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_fixcovc2zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      double chosen_fixcovc1zi {0.0};
      if (chosen_f1indc_cat != "none") {
        int delectable_sum = randindex(0, 16) + randindex(1, 16) + randindex(2, 16) +
          randindex(3, 16) + randindex(4, 16);
        for (int indcount = 0; indcount < randindex(5, 16); indcount++) {
          if (chosen_f1indc_cat == modelind_rownames_zi(indcount + delectable_sum)) {
            chosen_fixcovc1zi = modelindzi(indcount + delectable_sum);
          }
        }
      }
      
      preout_zi = (mainsum_zi + chosen_randcova2zi + chosen_randcova1zi +
        chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
        chosen_randcovc1zi + chosen_fixcova2zi + chosen_fixcova1zi +
        chosen_fixcovb2zi + chosen_fixcovb1zi + chosen_fixcovc2zi +
        chosen_fixcovc1zi + modelgroups2zi(grp2o_i) + modelgroups1zi(grp1_i) + 
        modelpatchzi(patchnumber) + modelyearzi(yearnumber) +
        dev_terms(placeholder));
    }
    
    if (vitaltype == 0) {
      // Binomial vital rates only
      if (preout > exp_tol) preout = exp_tol;
        
      double pre_exp = exp(preout);
      all_out = pre_exp / (1.0 + pre_exp);
      
      // Rcout << "Binomial: preout: " << preout << " pre_exp: " << pre_exp <<
      //   " all_out: " << all_out << "\n";
      
    } else if (vitaltype == 1) {
      // Size vital rates
      
      double Used_size3 = status_terms(16);
      double Used_binwidth3 = status_terms(19);
      
      if (vitalrate == 4) {
        Used_size3 = status_terms(17);
        Used_binwidth3 = status_terms(20);
        
      } else if (vitalrate == 5) {
        Used_size3 = status_terms(18);
        Used_binwidth3 = status_terms(21);
        
      } else if (vitalrate == 11) {
        Used_size3 = status_terms(17);
        Used_binwidth3 = status_terms(20);
        
      } else if (vitalrate == 12) {
        Used_size3 = status_terms(18);
        Used_binwidth3 = status_terms(21);
      }
      
      if (zi_processing) {
        // Development of estimate of pi (probability of 0 in binomial zi model)
        if (preout_zi > exp_tol) preout_zi = exp_tol;
        
        double pre_exp_zi = exp(preout_zi);
        all_out_zi = pre_exp_zi / (1.0 + pre_exp_zi);
        
        // Rcout << "ZI Binomial: preout_zi: " << preout_zi << " pre_exp_zi: " << pre_exp_zi <<
        //   " all_out_zi: " << all_out_zi << "\n";
      }
      
      if (vitaldist == 0) {
        // Poisson distribution
        
        if (preout > exp_tol) preout = exp_tol;
        double lambda = exp(preout);
        
        double upper_boundary = (Used_size3 + (Used_binwidth3 / 2));
        double upper_boundary_int = floor(upper_boundary);
        
        double lower_boundary = (Used_size3 - (Used_binwidth3 / 2));
        double lower_boundary_int = floor(lower_boundary);
        
        if (ipm_cdf) {
          if (lower_boundary_int < 0.0) lower_boundary_int = 0.0;
          
          double sizefac {1.0};
          if (upper_boundary_int > 0.0) {
            sizefac = upper_boundary_int * tgamma(upper_boundary_int);
          }
          double main_out = boost::math::tgamma((upper_boundary_int + 1), lambda) /
            sizefac;
          
          if (upper_boundary_int > lower_boundary_int) {
            double sizefac_low {1.0};
            if (lower_boundary_int > 0.0) {
              sizefac_low = lower_boundary_int * tgamma(lower_boundary_int);
            }
            all_out = main_out - boost::math::tgamma((lower_boundary_int + 1), lambda) /
              sizefac_low;
          } else {
            all_out = main_out;
          }
          if (all_out < 0.0) all_out = 0.0; // Eliminates issues in some versions of Linux
          
          if (modeltrunc == 1) {
            double den_corr = (1.0 - (exp(-1 * lambda)));
            if (den_corr == 0.0 || NumericVector::is_na(den_corr)) {
              den_corr = 1 / (exp_tol * exp_tol);
            }
            all_out = all_out / den_corr;
          }
          // Rcout << "Poisson cdf: upper_boundary_int: " << upper_boundary_int << 
          //   " lower_boundary_int: " << lower_boundary_int << " lambda: " << lambda << 
          //   " all_out: " << all_out << "\n";
          
        } else {
          int y = static_cast<int>(upper_boundary_int);
          int y0 = static_cast<int>(lower_boundary_int);
          if (y0 < -1) y0 = -1;
          
          double current_prob {0.0};
          
          for (int summed_size = (y0 + 1); summed_size <= y; summed_size++) {
            double sizefac {1.0};
            if (Used_size3 > 0.0) {
              sizefac = Used_size3 * tgamma(Used_size3);
            }
            
            double den_corr {1.0};
            if (modeltrunc == 1) den_corr = (1.0 - (exp(-1 * lambda)));
            if (den_corr == 0.0 || NumericVector::is_na(den_corr)) {
              den_corr = 1.0 / (exp_tol * exp_tol);
            }
            
            current_prob += ((pow(lambda, Used_size3) * exp(-1.0 * lambda)) / sizefac) /
              den_corr;
          }
          all_out = current_prob;
          
          // Rcout << "Poisson mid: upper_boundary_int: " << upper_boundary_int <<
          //   " lower_boundary_int: " << lower_boundary_int << " lambda: " << lambda << 
          //   " current_prob: " << current_prob << "\n";
        }
        
        if (zi_processing) {
          double current_pi = all_out_zi;
          double current_chi = all_out;
          
          if (Used_size3 == 0) {
            all_out = current_pi + ((1.0 - current_pi) * current_chi);
          } else {
            all_out = (1.0 - current_pi) * current_chi;
          }
        }
        
      } else if (vitaldist == 1) {
        // Negative binomial
        
        double mu = exp(preout);
        
        double theta = modelproxy["sigma"];
        if (NumericVector::is_na(theta)) theta = 1.0;
        if (theta > theta_tol) theta = theta_tol;
        double alpha = 1.0 / theta;
        
        double upper_boundary = (Used_size3 + (Used_binwidth3 / 2));
        double upper_boundary_int = floor(upper_boundary);
        int y = static_cast<int>(upper_boundary_int);
        
        double lower_boundary = (Used_size3 - (Used_binwidth3 / 2));
        double lower_boundary_int = floor(lower_boundary);
        int y0 = static_cast<int>(lower_boundary_int);
        if (y0 < -1) y0 = -1;
        
        double log_amu = log(alpha) + log(mu);
        double log_mid = -1.0 * theta * log(1.0 + (alpha * mu));
        double den_corr {1.0};
        if (modeltrunc == 1) den_corr = 1.0 - exp(log_mid);
        if (den_corr == 0.0 || NumericVector::is_na(den_corr)) {
          den_corr = 1 / (exp_tol * exp_tol);
        }
        
        double current_prob {0.0};
        
        for (int summed_size = (y0 + 1); summed_size <= y; summed_size++) {
          double log_leftie = 0.0;
          for (int j = 0; j < summed_size; j++) {
            log_leftie = log(static_cast<double>(j) + theta) -
              log(static_cast<double>(j) + 1.0) + log_leftie;
          }
          double log_rightie = static_cast<double>(summed_size) *
            (log_amu - log(1.0 + (alpha * mu)));
          
          double raw_prob = log_leftie + log_mid + log_rightie;
          
          current_prob += exp(raw_prob) / den_corr;
        }
        all_out = current_prob;
        
        if (zi_processing) {
          double current_pi = all_out_zi;
          double current_chi = all_out;
          
          if (Used_size3 == 0) {
            all_out = current_pi + ((1.0 - current_pi) * current_chi);
          } else {
            all_out = (1.0 - current_pi) * current_chi;
          }
        }
        
        if (all_out < 0.0) all_out = 0.0; // Eliminates issues in some versions of Linux
        
        // Rcout << "Negbin: y: " << y << " y0: " << y0 << " alpha: " << alpha <<
        //   " mu: " << mu << " current_prob: " << current_prob << "\n";
        
      } else if (vitaldist == 2) {
        // Gaussian size distribution, assuming midpoint
        
        if (ipm_cdf) {
          double lower_size = Used_size3 - (0.5 * Used_binwidth3);
          double upper_size = Used_size3 + (0.5 * Used_binwidth3);
          
          double lower_prob = normcdf(lower_size, preout, sigma);
          double upper_prob = normcdf(upper_size, preout, sigma);
          
          all_out = upper_prob - lower_prob;
          
          // Rcout << "Gaussian cdf: upper_size: " << upper_size << " lower_size: " <<
          //   lower_size << " Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " preout: " << preout << " sigma: " <<
          //   sigma << " upper_prob: " << upper_prob << " lower_prob: " <<
          //   lower_prob << " all_out: " << all_out << "\n";
        } else {
          double sigma2 = sigma * sigma;
          
          all_out = (exp(-1 * (pow((Used_size3 - preout), 2) / (2.0 * sigma2))) / 
            ((pow((2 * M_PI), 0.5)) * sigma));
          all_out = all_out * Used_binwidth3; // Midpoint integration
          
          // Rcout << "Gaussian mid: Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " sigma: " << sigma << " preout: " <<
          //   preout << " all_out: " << all_out << "\n";
        }
      } else if (vitaldist == 3) {
        // Gamma size distribution, assuming midpoint
        
        double E_y = 1 / preout;
        double sigma2 = sigma * sigma;
        double alpha = 1.0 / sigma2;
        double beta = (alpha / E_y);
        
        if (ipm_cdf) {
          double lower_size = Used_size3 - (0.5 * Used_binwidth3);
          double upper_size = Used_size3 + (0.5 * Used_binwidth3);
          
          double lower_prob = boost::math::gamma_p(alpha, (beta * lower_size));
          double upper_prob = boost::math::gamma_p(alpha, (beta * upper_size));
          
          all_out = upper_prob - lower_prob;
          
          // Rcout << "Gamma cdf: upper_size: " << upper_size << " lower_size: " <<
          //   lower_size << " alpha: " << alpha << " beta: " << beta << " upper_prob: " <<
          //   upper_prob << " lower_prob: " << lower_prob << " all_out: " << all_out << "\n";
        } else {
          
          all_out = pow(beta, alpha) * (1.0 / tgamma(alpha)) * 
            pow(Used_size3, (alpha - 1.0)) * exp(-1.0 * beta * Used_size3);
          all_out = all_out * Used_binwidth3; // Midpoint integration
          
          // Rcout << "Gamma mid: Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " alpha: " << alpha << " beta: " << beta <<
          //   " all_out: " << all_out << "\n";
        }
      }
      
    } else if (vitaltype == 2) {
      // Fecundity
      if (matrixformat != 2 || stage2n_i != static_cast<double>(nostages+1)) {
        if (vitaldist == 0 || vitaldist == 1) {
          // Poisson and negative binomial fecundity
          if (preout > exp_tol) preout = exp_tol;
          
          if (zi_processing) {
            
            all_out = (1.0 - (exp(preout_zi) / (1.0 + exp(preout_zi)))) *
              (exp(preout) * fecmod * repentry_i);
            
          } else {
            
            all_out = exp(preout) * fecmod * repentry_i;
          }
        } else if (vitaldist == 2) {
          // Gaussian fecundity
          all_out = preout * fecmod * repentry_i;
          
          if (negfec && all_out < 0.0) all_out = 0.0;
          
        } else if (vitaldist == 3) {
          // Gamma fecundity
          all_out = (1.0 / preout) * fecmod * repentry_i;
        } else {
          all_out = maincoefs(0);
        }
      } else if (stage2n_i == static_cast<double>(nostages+1)) {
        // Fecundity in deVries-formatted hMPMs
        if (vitaldist == 0 || vitaldist == 1) {
          // Poisson and negative binomial fecundity
          
          if (preout > exp_tol) preout = exp_tol;
              
          if (zi_processing) {
            all_out = (1.0 - (exp(preout_zi) / (1.0 + exp(preout_zi)))) *
              (exp(preout) * fecmod * repentry_i);
            
          } else {
            all_out = exp(preout) * fecmod * repentry_i;
          }
        } else if (vitaldist == 2) {
          // Gaussian fecundity
          all_out = preout * fecmod * repentry_i;
          
          if (negfec && all_out < 0.0) {
            all_out = 0.0;
          }
        } else if (vitaldist == 3) {
          // Gamma fecundity
          all_out = (1.0 / preout) * fecmod * repentry_i;
        }
      } else {
        all_out = maincoefs(0);
      }
    }
    
    return all_out;
  }
  
  //' Estimate All Elements of Function-based Population Projection Matrix
  //' 
  //' Function \code{jerzeibalowski()} calculates matrix elements in
  //' function-based population projection matrices involving stages. Used in
  //' \code{mpm_create()}, and through that function in \code{flefko3()},
  //' \code{flefko2()}, and \code{aflefko2()}.
  //' 
  //' @name jerzeibalowski
  //' 
  //' @param AllStages Data frame with all required inputs for vital rate
  //' estimation other than the vital rate model coefficients themselves.
  //' Contains a row for each matrix element.
  //' @param stageframe Modified stageframe used in matrix calculations.
  //' @param matrixformat Integer representing the style of matrix to develop.
  //' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
  //' and age-by-stage MPM (4).
  //' @param survproxy List of coefficients estimated in model of survival.
  //' @param obsproxy List of coefficients estimated in model of observation.
  //' @param sizeproxy List of coefficients estimated in model of primary size.
  //' @param sizebproxy List of coefficients estimated in model of secondary
  //' size.
  //' @param sizecproxy List of coefficients estimated in model of tertiary
  //' size.
  //' @param repstproxy List of coefficients estimated in model of reproductive 
  //' status.
  //' @param fecproxy List of coefficients estimated in model of fecundity.
  //' @param jsurvproxy List of coefficients estimated in model of juvenile
  //' survival.
  //' @param jobsproxy List of coefficients estimated in model of juvenile
  //' observation.
  //' @param jsizeproxy List of coefficients estimated in model of juvenile
  //' primary size.
  //' @param jsizebproxy List of coefficients estimated in model of juvenile
  //' secondary size.
  //' @param jsizecproxy List of coefficients estimated in model of juvenile
  //' tertiary size.
  //' @param jrepstproxy List of coefficients estimated in model of juvenile
  //' reproductive status.
  //' @param jmatstproxy List of coefficients estimated in model of juvenile
  //' maturity probability.
  //' @param f2_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_anna A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{a} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_anna A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{a} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_annb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{b} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_annb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{b} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_annc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{c} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_annc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{c} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_inda_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{a} at each time \emph{t}
  //' to be used in analysis.
  //' @param f1_inda_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{a} at each time \emph{t-1}
  //' to be used in analysis.
  //' @param f2_indb_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{b} at each time \emph{t}
  //' to be used in analysis.
  //' @param f1_indb_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{b} at each time \emph{t-1}
  //' to be used in analysis.
  //' @param f2_indc_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{c} at each time \emph{t}
  //' to be used in analysis.
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
  //' models input by the user. The order is: survival, observation status, size,
  //' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
  //' observation status, juvenile size, juvenile size_b, juvenile size_c,
  //' juvenile reproductive status, and juvenile maturity status.
  //' @param dens_vr A logical value indicating whether any vital rates are
  //' density dependent.
  //' @param dvr_yn A logical vector indicating whether each vital rate is density
  //' dependent.
  //' @param dvr_style An integer vector indicating the style of density
  //' dependence for each vital rate.
  //' @param dvr_alpha A numeric vector indicating the value of alpha to use in
  //' density dependence for each vital rate.
  //' @param dvr_beta A numeric vector indicating the value of beta to use in
  //' density dependence for each vital rate.
  //' @param dens_n A numeric vector corresponding to the population size to use
  //' in vital rate density dependence calculations.
  //' @param dens A numeric value equal to the spatial density to be used in
  //' calculations.
  //' @param fecmod A scalar multiplier for fecundity.
  //' @param maxsize The maximum primary size to be used in element estimation.
  //' @param maxsizeb The maximum secondary size to be used in element estimation.
  //' @param maxsizec The maximum tertiary size to be used in element estimation.
  //' @param firstage The first age to be included in age-by-stage MPM estimation.
  //' @param finalage The final age to be included in age-by-stage MPM estimation.
  //' @param negfec A logical value denoting whether to change negative estimated
  //' fecundity to 0.
  //' @param yearnumber An integer specifying which time at time \emph{t} to
  //' develop matrices for. Must be in reference to the \code{listofyears} object
  //' developed in the \code{R} matrix estimator function.
  //' @param patchnumber An integer specifying which patch to develop matrices
  //' for. Must be in reference to the \code{listofyears} object developed in the
  //' \code{R} matrix estimator function.
  //' @param exp_tol A numeric value indicating the maximum limit for the
  //' \code{exp()} function to be used in vital rate calculations. Defaults to
  //' \code{700.0}.
  //' @param theta_tol A numeric value indicating a maximum value for theta in
  //' negative binomial probability density estimation. Defaults to
  //' \code{100000000.0}.
  //' @param ipm_cdf A logical value indicating whether to estimate size
  //' transitions using the cumulative density function in cases with continuous
  //' distributions. Defaults to \code{TRUE}, with the midpoint method used if
  //' \code{FALSE}.
  //' @param err_check A logical value indicating whether to export a matrix of
  //' conditional probabilities used to develop the \code{U} matrix. Defaults to
  //' \code{FALSE}.
  //' @param simplicity A logical value indicating whether to output all three
  //' matrices (\code{FALSE}), or just matrices \code{U} and \code{F}
  //' (\code{TRUE}). Defaults to \code{FALSE}.
  //' @param sparse A logical value indicating whether to output matrices in
  //' standard Armadillo format or in sparse matrix format. Defaults to
  //' \code{FALSE}.
  //' @param proj_only A logical value indicating whether to output all
  //' matrices outined in \code{simplicity}, or just the \code{A} matrix.
  //' 
  //' @return A list with 2, 3, or 4 elements. If \code{simplicity} is set to
  //' \code{FALSE}, then first 3 elements are matrices, including the main MPM
  //' (A), the survival-transition matrix (U), and a fecundity matrix (F). If
  //' simplicity is set to \code{TRUE}, then only the survival-transition matrix
  //' (U) and fecundity matrix (F) are output. If \code{err_check} is set to
  //' \code{TRUE}, then another element is added, which is a 7 column matrix
  //' showing survival probability, observation probability, reproduction
  //' probability, sizea transition probability, sizeb transition probability,
  //' sizec transition probability, and juvenile transition probability to
  //' maturity for each element of the final MPM. It is possible that, due to
  //' evolving development strategy, further columns are output, as well. Note
  //' that if \code{sparse = TRUE}, then output matrices are in sparse format.
  //' 
  //' @section Notes:
  //' The data frame AllStages introduces variables used in size and fecundity
  //' calculations. This DataFrame is broken up into long vectors composed of
  //' input sizes and related variables for these calculations. The "model" Lists
  //' bring in the vital rate models, and include random coefficients where
  //' needed. We also have a number of extra variables, that include such info as
  //' whether to use the Poisson, negative binomial, Gamma, or Gaussian
  //' distributions for size and fecundity calculations. If \code{sizedist},
  //' \code{sizebdist}, \code{sizecdist}, or \code{fecdist} equals 0, 1, 2, or 3,
  //' then the Poisson, negative binomial, Gaussian, or Gamma is used,
  //' respectively.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List jerzeibalowski(const DataFrame& AllStages,
    const DataFrame& stageframe, int matrixformat, const List& survproxy,
    const List& obsproxy, const List& sizeproxy, const List& sizebproxy,
    const List& sizecproxy, const List& repstproxy, const List& fecproxy,
    const List& jsurvproxy, const List& jobsproxy, const List& jsizeproxy,
    const List& jsizebproxy, const List& jsizecproxy, const List& jrepstproxy,
    const List& jmatstproxy, NumericVector f2_inda, NumericVector f1_inda,
    NumericVector f2_indb, NumericVector f1_indb, NumericVector f2_indc,
    NumericVector f1_indc, NumericVector f2_anna, NumericVector f1_anna,
    NumericVector f2_annb, NumericVector f1_annb, NumericVector f2_annc,
    NumericVector f1_annc, StringVector f2_inda_cat, StringVector f1_inda_cat,
    StringVector f2_indb_cat, StringVector f1_indb_cat, StringVector f2_indc_cat,
    StringVector f1_indc_cat, StringVector r2_inda, StringVector r1_inda,
    StringVector r2_indb, StringVector r1_indb, StringVector r2_indc,
    StringVector r1_indc, const NumericVector& dev_terms, bool dens_vr,
    LogicalVector dvr_yn, IntegerVector dvr_style, NumericVector dvr_alpha,
    NumericVector dvr_beta, NumericVector dens_n, double dens, double fecmod,
    double maxsize, double maxsizeb, double maxsizec,
    unsigned int firstage, unsigned int finalage, bool negfec, int yearnumber,
    int patchnumber, double exp_tol = 700.0, double theta_tol = 100000000.0,
    bool ipm_cdf = true, bool err_check = false, bool simplicity = false,
    bool sparse = false, bool proj_only = false) {
    
    NumericMatrix out;
    
    StringVector stagenames = stageframe["stage"];
    int nostages = stagenames.length();
    unsigned long matrixdim {0};
    
    int nostages_counter = nostages;
    for (int i = 0; i < nostages_counter; i++) {
      if (stringcompare_hard(as<std::string>(stagenames(i)), "AlmostBorn")) nostages -= 1;  
      if (stringcompare_hard(as<std::string>(stagenames(i)), "Dead")) nostages -= 1;
    }
    
    if (matrixformat == 1) { // Ehrlen-format hMPM
      matrixdim = nostages * nostages;
    } else if (matrixformat == 2) { // deVries-format hMPM
      matrixdim = nostages * (nostages + 1);
    } else if (matrixformat == 3) { // ahMPM
      matrixdim = nostages;
    } else if (matrixformat == 4) { // age-by-stage MPM
      matrixdim = nostages * (finalage - firstage + 1);
    }
    
    // Proxy model imports and settings
    bool sizezero = as<bool>(sizeproxy["zero_inflated"]);
    bool sizebzero = as<bool>(sizebproxy["zero_inflated"]);
    bool sizeczero = as<bool>(sizecproxy["zero_inflated"]);
    bool feczero = as<bool>(fecproxy["zero_inflated"]);
    bool jsizezero = as<bool>(jsizeproxy["zero_inflated"]);
    bool jsizebzero = as<bool>(jsizebproxy["zero_inflated"]);
    bool jsizeczero = as<bool>(jsizecproxy["zero_inflated"]);
    
    bool sizetrunc = as<bool>(sizeproxy["zero_truncated"]);
    bool sizebtrunc = as<bool>(sizebproxy["zero_truncated"]);
    bool sizectrunc = as<bool>(sizecproxy["zero_truncated"]);
    bool fectrunc = as<bool>(fecproxy["zero_truncated"]);
    bool jsizetrunc = as<bool>(jsizeproxy["zero_truncated"]);
    bool jsizebtrunc = as<bool>(jsizebproxy["zero_truncated"]);
    bool jsizectrunc = as<bool>(jsizecproxy["zero_truncated"]);
    
    NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
    NumericVector obscoefs = as<NumericVector>(obsproxy["coefficients"]);
    NumericVector sizecoefs = as<NumericVector>(sizeproxy["coefficients"]);
    NumericVector sizebcoefs = as<NumericVector>(sizebproxy["coefficients"]);
    NumericVector sizeccoefs = as<NumericVector>(sizecproxy["coefficients"]);
    NumericVector repstcoefs = as<NumericVector>(repstproxy["coefficients"]);
    NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
    NumericVector jsurvcoefs = as<NumericVector>(jsurvproxy["coefficients"]);
    NumericVector jobscoefs = as<NumericVector>(jobsproxy["coefficients"]);
    NumericVector jsizecoefs = as<NumericVector>(jsizeproxy["coefficients"]);
    NumericVector jsizebcoefs = as<NumericVector>(jsizebproxy["coefficients"]);
    NumericVector jsizeccoefs = as<NumericVector>(jsizecproxy["coefficients"]);
    NumericVector jrepstcoefs = as<NumericVector>(jrepstproxy["coefficients"]);
    NumericVector jmatstcoefs = as<NumericVector>(jmatstproxy["coefficients"]);
    
    double survsigma = as<double>(survproxy["sigma"]);
    double obssigma = as<double>(obsproxy["sigma"]);
    double sizesigma = as<double>(sizeproxy["sigma"]);
    double sizebsigma = as<double>(sizebproxy["sigma"]);
    double sizecsigma = as<double>(sizecproxy["sigma"]);
    double repstsigma = as<double>(repstproxy["sigma"]);
    double fecsigma = as<double>(fecproxy["sigma"]);
    double jsurvsigma = as<double>(jsurvproxy["sigma"]);
    double jobssigma = as<double>(jobsproxy["sigma"]);
    double jsizesigma = as<double>(jsizeproxy["sigma"]);
    double jsizebsigma = as<double>(jsizebproxy["sigma"]);
    double jsizecsigma = as<double>(jsizecproxy["sigma"]);
    double jrepstsigma = as<double>(jrepstproxy["sigma"]);
    double jmatstsigma = as<double>(jmatstproxy["sigma"]);
    
    int survdist = as<int>(survproxy["dist"]);
    int obsdist = as<int>(obsproxy["dist"]);
    int sizedist = as<int>(sizeproxy["dist"]);
    int sizebdist = as<int>(sizebproxy["dist"]);
    int sizecdist = as<int>(sizecproxy["dist"]);
    int repstdist = as<int>(repstproxy["dist"]);
    int fecdist = as<int>(fecproxy["dist"]);
    int jsurvdist = as<int>(jsurvproxy["dist"]);
    int jobsdist = as<int>(jobsproxy["dist"]);
    int jsizedist = as<int>(jsizeproxy["dist"]);
    int jsizebdist = as<int>(jsizebproxy["dist"]);
    int jsizecdist = as<int>(jsizecproxy["dist"]);
    int jrepstdist = as<int>(jrepstproxy["dist"]);
    int jmatstdist = as<int>(jmatstproxy["dist"]);
    
    if (NumericVector::is_na(sizesigma)) {
      if (sizedist == 1) {
        sizesigma = 1.0;
      } else {
        sizesigma = 0.0;
      }
    }
    if (NumericVector::is_na(sizebsigma)) {
      if (sizebdist == 1) {
        sizebsigma = 1.0;
      } else {
        sizebsigma = 0.0;
      }
    }
    if (NumericVector::is_na(sizecsigma)) {
      if (sizecdist == 1) {
        sizecsigma = 1.0;
      } else {
        sizecsigma = 0.0;
      }
    }
    if (NumericVector::is_na(fecsigma)) {
      if (fecdist == 1) {
        fecsigma = 1.0;
      } else {
        fecsigma = 0.0;
      }
    }
    if (NumericVector::is_na(jsizesigma)) {
      if (sizedist == 1) {
        jsizesigma = 1.0;
      } else {
        jsizesigma = 0.0;
      }
    }
    if (NumericVector::is_na(jsizebsigma)) {
      if (sizebdist == 1) {
        jsizebsigma = 1.0;
      } else {
        jsizebsigma = 0.0;
      }
    }
    if (NumericVector::is_na(jsizecsigma)) {
      if (sizecdist == 1) {
        jsizecsigma = 1.0;
      } else {
        jsizecsigma = 0.0;
      }
    }
    
    NumericMatrix vital_year = revelations(survproxy, obsproxy, sizeproxy,
      sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
      jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 1);
    
    NumericMatrix vital_patch = revelations(survproxy, obsproxy, sizeproxy,
      sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
      jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 2);
    
    arma::imat rand_index = foi_index(survproxy, obsproxy, sizeproxy,
      sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
      jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy);
    
    // NumericVector imports from model_proxy objects
    NumericVector sizeyearzi = as<NumericVector>(sizeproxy["zeroyear"]);
    NumericVector sizebyearzi = as<NumericVector>(sizebproxy["zeroyear"]);
    NumericVector sizecyearzi = as<NumericVector>(sizecproxy["zeroyear"]);
    NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
    NumericVector jsizeyearzi = as<NumericVector>(jsizeproxy["zeroyear"]);
    NumericVector jsizebyearzi = as<NumericVector>(jsizebproxy["zeroyear"]);
    NumericVector jsizecyearzi = as<NumericVector>(jsizecproxy["zeroyear"]);
    
    NumericVector dud_yearzi(sizeyearzi.length());
    
    NumericVector unisyzi = unique(sizeyearzi);
    NumericVector unisyzbi = unique(sizebyearzi);
    NumericVector unisyzci = unique(sizecyearzi);
    NumericVector unijsyzi = unique(jsizeyearzi);
    NumericVector unijsyzbi = unique(jsizebyearzi);
    NumericVector unijsyzci = unique(jsizecyearzi);
    NumericVector unifeci = unique(fecyearzi);
    
    NumericVector sizepatchzi = as<NumericVector>(sizeproxy["zeropatch"]);
    NumericVector sizebpatchzi = as<NumericVector>(sizebproxy["zeropatch"]);
    NumericVector sizecpatchzi = as<NumericVector>(sizecproxy["zeropatch"]);
    NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
    NumericVector jsizepatchzi = as<NumericVector>(jsizeproxy["zeropatch"]);
    NumericVector jsizebpatchzi = as<NumericVector>(jsizebproxy["zeropatch"]);
    NumericVector jsizecpatchzi = as<NumericVector>(jsizecproxy["zeropatch"]);
    
    NumericVector dud_patchzi(sizepatchzi.length());
    
    NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
    NumericVector obsgroups2 = as<NumericVector>(obsproxy["groups2"]);
    NumericVector sizegroups2 = as<NumericVector>(sizeproxy["groups2"]);
    NumericVector sizebgroups2 = as<NumericVector>(sizebproxy["groups2"]);
    NumericVector sizecgroups2 = as<NumericVector>(sizecproxy["groups2"]);
    NumericVector repstgroups2 = as<NumericVector>(repstproxy["groups2"]);
    NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
    NumericVector jsurvgroups2 = as<NumericVector>(jsurvproxy["groups2"]);
    NumericVector jobsgroups2 = as<NumericVector>(jobsproxy["groups2"]);
    NumericVector jsizegroups2 = as<NumericVector>(jsizeproxy["groups2"]);
    NumericVector jsizebgroups2 = as<NumericVector>(jsizebproxy["groups2"]);
    NumericVector jsizecgroups2 = as<NumericVector>(jsizecproxy["groups2"]);
    NumericVector jrepstgroups2 = as<NumericVector>(jrepstproxy["groups2"]);
    NumericVector jmatstgroups2 = as<NumericVector>(jmatstproxy["groups2"]);
    
    NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
    NumericVector obsgroups1 = as<NumericVector>(obsproxy["groups1"]);
    NumericVector sizegroups1 = as<NumericVector>(sizeproxy["groups1"]);
    NumericVector sizebgroups1 = as<NumericVector>(sizebproxy["groups1"]);
    NumericVector sizecgroups1 = as<NumericVector>(sizecproxy["groups1"]);
    NumericVector repstgroups1 = as<NumericVector>(repstproxy["groups1"]);
    NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
    NumericVector jsurvgroups1 = as<NumericVector>(jsurvproxy["groups1"]);
    NumericVector jobsgroups1 = as<NumericVector>(jobsproxy["groups1"]);
    NumericVector jsizegroups1 = as<NumericVector>(jsizeproxy["groups1"]);
    NumericVector jsizebgroups1 = as<NumericVector>(jsizebproxy["groups1"]);
    NumericVector jsizecgroups1 = as<NumericVector>(jsizecproxy["groups1"]);
    NumericVector jrepstgroups1 = as<NumericVector>(jrepstproxy["groups1"]);
    NumericVector jmatstgroups1 = as<NumericVector>(jmatstproxy["groups1"]);
    
    NumericVector sizegroups2zi = as<NumericVector>(sizeproxy["zerogroups2"]);
    NumericVector sizebgroups2zi = as<NumericVector>(sizebproxy["zerogroups2"]);
    NumericVector sizecgroups2zi = as<NumericVector>(sizecproxy["zerogroups2"]);
    NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
    NumericVector jsizegroups2zi = as<NumericVector>(jsizeproxy["zerogroups2"]);
    NumericVector jsizebgroups2zi = as<NumericVector>(jsizebproxy["zerogroups2"]);
    NumericVector jsizecgroups2zi = as<NumericVector>(jsizecproxy["zerogroups2"]);
    
    NumericVector dud_groups2zi(jsizecyearzi.length());
    
    NumericVector sizegroups1zi = as<NumericVector>(sizeproxy["zerogroups1"]);
    NumericVector sizebgroups1zi = as<NumericVector>(sizebproxy["zerogroups1"]);
    NumericVector sizecgroups1zi = as<NumericVector>(sizecproxy["zerogroups1"]);
    NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
    NumericVector jsizegroups1zi = as<NumericVector>(jsizeproxy["zerogroups1"]);
    NumericVector jsizebgroups1zi = as<NumericVector>(jsizebproxy["zerogroups1"]);
    NumericVector jsizecgroups1zi = as<NumericVector>(jsizecproxy["zerogroups1"]);
    
    NumericVector dud_groups1zi(jsizecyearzi.length());
    
    NumericVector survind = flightoficarus(survproxy);
    NumericVector obsind = flightoficarus(obsproxy);
    NumericVector sizeind = flightoficarus(sizeproxy);
    NumericVector sizebind = flightoficarus(sizebproxy);
    NumericVector sizecind = flightoficarus(sizecproxy);
    NumericVector repstind = flightoficarus(repstproxy);
    NumericVector fecind = flightoficarus(fecproxy);
    NumericVector jsurvind = flightoficarus(jsurvproxy);
    NumericVector jobsind = flightoficarus(jobsproxy);
    NumericVector jsizeind = flightoficarus(jsizeproxy);
    NumericVector jsizebind = flightoficarus(jsizebproxy);
    NumericVector jsizecind = flightoficarus(jsizecproxy);
    NumericVector jrepstind = flightoficarus(jrepstproxy);
    NumericVector jmatstind = flightoficarus(jmatstproxy);
    
    NumericVector sizeindzi = zero_flightoficarus(sizeproxy);
    NumericVector sizebindzi = zero_flightoficarus(sizebproxy);
    NumericVector sizecindzi = zero_flightoficarus(sizecproxy);
    NumericVector fecindzi = zero_flightoficarus(fecproxy);
    NumericVector jsizeindzi = zero_flightoficarus(jsizeproxy);
    NumericVector jsizebindzi = zero_flightoficarus(jsizebproxy);
    NumericVector jsizecindzi = zero_flightoficarus(jsizecproxy);
    
    StringVector survind_rownames = bootson(survproxy);
    StringVector obsind_rownames = bootson(obsproxy);
    StringVector sizeind_rownames = bootson(sizeproxy);
    StringVector sizebind_rownames = bootson(sizebproxy);
    StringVector sizecind_rownames = bootson(sizecproxy);
    StringVector repstind_rownames = bootson(repstproxy);
    StringVector fecind_rownames = bootson(fecproxy);
    StringVector jsurvind_rownames = bootson(jsurvproxy);
    StringVector jobsind_rownames = bootson(jobsproxy);
    StringVector jsizeind_rownames = bootson(jsizeproxy);
    StringVector jsizebind_rownames = bootson(jsizebproxy);
    StringVector jsizecind_rownames = bootson(jsizecproxy);
    StringVector jrepstind_rownames = bootson(jrepstproxy);
    StringVector jmatstind_rownames = bootson(jmatstproxy);
    
    StringVector sizeind_rownames_zi = zero_bootson(sizeproxy);
    StringVector sizebind_rownames_zi = zero_bootson(sizebproxy);
    StringVector sizecind_rownames_zi = zero_bootson(sizecproxy);
    StringVector fecind_rownames_zi = zero_bootson(fecproxy);
    StringVector jsizeind_rownames_zi = zero_bootson(jsizeproxy);
    StringVector jsizebind_rownames_zi = zero_bootson(jsizebproxy);
    StringVector jsizecind_rownames_zi = zero_bootson(jsizecproxy);
    
    // AllStages import and settings
    Rcpp::NumericVector stage3_num = as<NumericVector>(AllStages["stage3"]);
    Rcpp::NumericVector stage2n_num = as<NumericVector>(AllStages["stage2n"]);
    Rcpp::NumericVector stage2o_num = as<NumericVector>(AllStages["stage2o"]);
    arma::vec stage3 = as<arma::vec>(stage3_num);
    arma::vec stage2n = as<arma::vec>(stage2n_num);
    arma::vec stage2o = as<arma::vec>(stage2o_num);
    
    Rcpp::NumericVector sz3 = as<NumericVector>(AllStages["size3"]);
    Rcpp::NumericVector sz2n = as<NumericVector>(AllStages["size2n"]);
    Rcpp::NumericVector sz2o = as<NumericVector>(AllStages["size2o"]);
    Rcpp::NumericVector sz1 = as<NumericVector>(AllStages["size1"]);
    Rcpp::NumericVector szb3 = as<NumericVector>(AllStages["sizeb3"]);
    Rcpp::NumericVector szb2n = as<NumericVector>(AllStages["sizeb2n"]);
    Rcpp::NumericVector szb2o = as<NumericVector>(AllStages["sizeb2o"]);
    Rcpp::NumericVector szb1 = as<NumericVector>(AllStages["sizeb1"]);
    Rcpp::NumericVector szc3 = as<NumericVector>(AllStages["sizec3"]);
    Rcpp::NumericVector szc2n = as<NumericVector>(AllStages["sizec2n"]);
    Rcpp::NumericVector szc2o = as<NumericVector>(AllStages["sizec2o"]);
    Rcpp::NumericVector szc1 = as<NumericVector>(AllStages["sizec1"]);
    Rcpp::NumericVector ob3 = as<NumericVector>(AllStages["obs3"]);
    Rcpp::NumericVector fl3 = as<NumericVector>(AllStages["rep3"]);
    Rcpp::NumericVector fl2n = as<NumericVector>(AllStages["rep2n"]);
    Rcpp::NumericVector fl2o = as<NumericVector>(AllStages["rep2o"]);
    Rcpp::NumericVector fl1 = as<NumericVector>(AllStages["rep1"]);
    Rcpp::NumericVector mat3 = as<NumericVector>(AllStages["mat3"]);
    Rcpp::NumericVector mat2n = as<NumericVector>(AllStages["mat2n"]);
    Rcpp::NumericVector mat2o = as<NumericVector>(AllStages["mat2o"]);
    Rcpp::NumericVector mat1 = as<NumericVector>(AllStages["mat1"]);
    Rcpp::NumericVector immat2n = as<NumericVector>(AllStages["imm2n"]);
    Rcpp::NumericVector immat2o = as<NumericVector>(AllStages["imm2o"]);
    Rcpp::NumericVector immat1 = as<NumericVector>(AllStages["imm1"]);
    
    Rcpp::NumericVector repentry = as<NumericVector>(AllStages["repentry3"]);
    Rcpp::NumericVector indata2n = as<NumericVector>(AllStages["indata2n"]);
    Rcpp::NumericVector indata2o = as<NumericVector>(AllStages["indata2o"]);
    Rcpp::NumericVector binwidth3 = as<NumericVector>(AllStages["binwidth"]);
    Rcpp::NumericVector binbwidth3 = as<NumericVector>(AllStages["binbwidth"]);
    Rcpp::NumericVector bincwidth3 = as<NumericVector>(AllStages["bincwidth"]);
    Rcpp::NumericVector actualage2 = as<NumericVector>(AllStages["actualage"]);
    
    Rcpp::NumericVector grp3 = as<NumericVector>(AllStages["group3"]);
    Rcpp::NumericVector grp2n = as<NumericVector>(AllStages["group2n"]);
    Rcpp::NumericVector grp2o = as<NumericVector>(AllStages["group2o"]);
    Rcpp::NumericVector grp1 = as<NumericVector>(AllStages["group1"]);
    
    Rcpp::NumericVector indata = as<NumericVector>(AllStages["indata"]);
    
    Rcpp::NumericVector ovgivent = as<NumericVector>(AllStages["ovgiven_t"]);
    Rcpp::NumericVector ovgivenf = as<NumericVector>(AllStages["ovgiven_f"]);
    Rcpp::NumericVector ovoffsett = as<NumericVector>(AllStages["ovoffset_t"]);
    Rcpp::NumericVector ovoffsetf = as<NumericVector>(AllStages["ovoffset_f"]);
    arma::vec ovostt = as<arma::vec>(ovoffsett);
    arma::vec ovostf = as<arma::vec>(ovoffsetf);
    
    Rcpp::NumericVector ovestt_num = as<NumericVector>(AllStages["ovest_t"]);
    arma::vec ovestt = as<arma::vec>(ovestt_num);
    
    Rcpp::NumericVector ovestf_num = as<NumericVector>(AllStages["ovest_f"]);
    arma::vec ovestf = as<arma::vec>(ovestf_num);
    
    Rcpp::NumericVector ovsurvmult = as<NumericVector>(AllStages["ovsurvmult"]);
    Rcpp::NumericVector ovfecmult = as<NumericVector>(AllStages["ovfecmult"]);
    arma::vec ovsurvmult_arma = as<arma::vec>(ovsurvmult);
    arma::vec ovfecmult_arma = as<arma::vec>(ovfecmult);
    
    Rcpp::IntegerVector index321_int = as<IntegerVector>(AllStages["index321"]);
    arma::uvec index321 = as<arma::uvec>(index321_int);
    
    Rcpp::IntegerVector aliveandequal = as<IntegerVector>(AllStages["aliveandequal"]);
    
    int n = static_cast<int>(stage3.n_elem);
    
    arma::uvec offsettvec = find(ovostt != 0.);
    arma::uvec offsetfvec = find(ovostf != 0.);
    int offsetst = static_cast<int>(offsettvec.n_elem);
    int offsetsf = static_cast<int>(offsetfvec.n_elem);
    
    arma::uvec replacetvec = find(ovestt != -1.);
    arma::uvec replacefvec = find(ovestf != -1.);
    int replacementst = static_cast<int>(replacetvec.n_elem);
    int replacementsf = static_cast<int>(replacefvec.n_elem);
    
    arma::uvec tmults = find(ovsurvmult_arma != -1.0);
    arma::uvec fmults = find(ovfecmult_arma != -1.0);
    arma::uvec noreplacetvec = find(ovestt == -1.0);
    arma::uvec noreplacefvec = find(ovestf == -1.0);
    arma::uvec tmults_only = intersect(tmults, noreplacetvec);
    arma::uvec fmults_only = intersect(fmults, noreplacefvec);
    int tmults_only_st = static_cast<int>(tmults_only.n_elem);
    int fmults_only_st = static_cast<int>(fmults_only.n_elem);
    
    int repindex {0};
    int properindex {0};
    int proxyindex {0};
    
    // Determination of choices of fixed and random individual covariates,
    // and annual covariates
    double inda1 = f1_inda(yearnumber);
    double indb1 = f1_indb(yearnumber);
    double indc1 = f1_indc(yearnumber);
    double inda2 = f2_inda(yearnumber);
    double indb2 = f2_indb(yearnumber);
    double indc2 = f2_indc(yearnumber);
    
    String chosen_f2inda_cat = f2_inda_cat(yearnumber);
    String chosen_f1inda_cat = f1_inda_cat(yearnumber);
    String chosen_f2indb_cat = f2_indb_cat(yearnumber);
    String chosen_f1indb_cat = f1_indb_cat(yearnumber);
    String chosen_f2indc_cat = f2_indc_cat(yearnumber);
    String chosen_f1indc_cat = f1_indc_cat(yearnumber);
    
    String chosen_r2inda = r2_inda(yearnumber);
    String chosen_r1inda = r1_inda(yearnumber);
    String chosen_r2indb = r2_indb(yearnumber);
    String chosen_r1indb = r1_indb(yearnumber);
    String chosen_r2indc = r2_indc(yearnumber);
    String chosen_r1indc = r1_indc(yearnumber);
    
    double anna1 = f1_anna(yearnumber);
    double annb1 = f1_annb(yearnumber);
    double annc1 = f1_annc(yearnumber);
    double anna2 = f2_anna(yearnumber);
    double annb2 = f2_annb(yearnumber);
    double annc2 = f2_annc(yearnumber);
    
    // Density corrections under vital rate density dependence
    double vr1_dcorr {1.0};
    double vr2_dcorr {1.0};
    double vr3_dcorr {1.0};
    double vr4_dcorr {1.0};
    double vr5_dcorr {1.0};
    double vr6_dcorr {1.0};
    double vr7_dcorr {1.0};
    double vr8_dcorr {1.0};
    double vr9_dcorr {1.0};
    double vr10_dcorr {1.0};
    double vr11_dcorr {1.0};
    double vr12_dcorr {1.0};
    double vr13_dcorr {1.0};
    double vr14_dcorr {1.0};
    
    if (dens_vr) {
      // Adult survival
      if (dvr_yn(0)) {
        if (dvr_style(0) == 1) {
          vr1_dcorr = dvr_alpha(0) * exp(-1 * dvr_beta(0) * dens_n(0));
        } else if (dvr_style(0) == 2) {
          vr1_dcorr = dvr_alpha(0) / (1 + dvr_beta(0) * dens_n(0));
        } else if (dvr_style(0) == 3) {
          vr1_dcorr = 1 / (1 + exp(dvr_alpha(0) * dens_n(0) + dvr_beta(0)));
        } else if (dvr_style(0) == 4) {
          vr1_dcorr = 1 - (dens_n(0) / dvr_alpha(0));
          if (dvr_beta(0) != 0 && vr1_dcorr < -1.0) vr1_dcorr = 0.0; 
        }
      }
      
      // Adult observation
      if (dvr_yn(1)) {
        if (dvr_style(1) == 1) {
          vr2_dcorr = dvr_alpha(1) * exp(-1 * dvr_beta(1) * dens_n(1));
        } else if (dvr_style(1) == 2) {
          vr2_dcorr = dvr_alpha(1) / (1 + dvr_beta(1) * dens_n(1));
        } else if (dvr_style(1) == 3) {
          vr2_dcorr = 1 / (1 + exp(dvr_alpha(1) * dens_n(1) + dvr_beta(1)));
        } else if (dvr_style(1) == 4) {
          vr2_dcorr = 1 - (dens_n(1) / dvr_alpha(1));
          if (dvr_beta(1) != 0 && vr2_dcorr < -1.0) vr2_dcorr = 0.0; 
        }
      }
      
      // Adult sizea
      if (dvr_yn(2)) {
        if (dvr_style(2) == 1) {
          vr3_dcorr = dvr_alpha(2) * exp(-1 * dvr_beta(2) * dens_n(2));
        } else if (dvr_style(2) == 2) {
          vr3_dcorr = dvr_alpha(2) / (1 + dvr_beta(2) * dens_n(2));
        } else if (dvr_style(2) == 3) {
          vr3_dcorr = 1 / (1 + exp(dvr_alpha(2) * dens_n(2) + dvr_beta(2)));
        } else if (dvr_style(2) == 4) {
          vr3_dcorr = 1 - (dens_n(2) / dvr_alpha(2));
          if (dvr_beta(2) != 0 && vr3_dcorr < -1.0) vr3_dcorr = 0.0; 
        }
      }
      
      // Adult sizeb
      if (dvr_yn(3)) {
        if (dvr_style(3) == 1) {
          vr4_dcorr = dvr_alpha(3) * exp(-1 * dvr_beta(3) * dens_n(3));
        } else if (dvr_style(3) == 2) {
          vr4_dcorr = dvr_alpha(3) / (1 + dvr_beta(3) * dens_n(3));
        } else if (dvr_style(3) == 3) {
          vr4_dcorr = 1 / (1 + exp(dvr_alpha(3) * dens_n(3) + dvr_beta(3)));
        } else if (dvr_style(3) == 4) {
          vr4_dcorr = 1 - (dens_n(3) / dvr_alpha(3));
          if (dvr_beta(3) != 0 && vr4_dcorr < -1.0) vr4_dcorr = 0.0; 
        }
      }
      
      // Adult sizec
      if (dvr_yn(4)) {
        if (dvr_style(4) == 1) {
          vr5_dcorr = dvr_alpha(4) * exp(-1 * dvr_beta(4) * dens_n(4));
        } else if (dvr_style(4) == 2) {
          vr5_dcorr = dvr_alpha(4) / (1 + dvr_beta(4) * dens_n(4));
        } else if (dvr_style(4) == 3) {
          vr5_dcorr = 1 / (1 + exp(dvr_alpha(4) * dens_n(4) + dvr_beta(4)));
        } else if (dvr_style(4) == 4) {
          vr5_dcorr = 1 - (dens_n(4) / dvr_alpha(4));
          if (dvr_beta(4) != 0 && vr5_dcorr < -1.0) vr5_dcorr = 0.0; 
        }
      }
      
      // Adult reproduction
      if (dvr_yn(5)) {
        if (dvr_style(5) == 1) {
          vr6_dcorr = dvr_alpha(5) * exp(-1 * dvr_beta(5) * dens_n(5));
        } else if (dvr_style(5) == 2) {
          vr6_dcorr = dvr_alpha(5) / (1 + dvr_beta(5) * dens_n(5));
        } else if (dvr_style(5) == 3) {
          vr6_dcorr = 1 / (1 + exp(dvr_alpha(5) * dens_n(5) + dvr_beta(5)));
        } else if (dvr_style(5) == 4) {
          vr6_dcorr = 1 - (dens_n(5) / dvr_alpha(5));
          if (dvr_beta(5) != 0 && vr6_dcorr < -1.0) vr6_dcorr = 0.0; 
        }
      }
      
      // Adult fecundity
      if (dvr_yn(6)) {
        if (dvr_style(6) == 1) {
          vr7_dcorr = dvr_alpha(6) * exp(-1 * dvr_beta(6) * dens_n(6));
        } else if (dvr_style(6) == 2) {
          vr7_dcorr = dvr_alpha(6) / (1 + dvr_beta(6) * dens_n(6));
        } else if (dvr_style(6) == 3) {
          vr7_dcorr = 1 / (1 + exp(dvr_alpha(6) * dens_n(6) + dvr_beta(6)));
        } else if (dvr_style(6) == 4) {
          vr7_dcorr = 1 - (dens_n(6) / dvr_alpha(6));
          if (dvr_beta(6) != 0 && vr7_dcorr < -1.0) vr7_dcorr = 0.0; 
        }
      }
      
      // Juvenile survival
      if (dvr_yn(7)) {
        if (dvr_style(7) == 1) {
          vr8_dcorr = dvr_alpha(7) * exp(-1 * dvr_beta(7) * dens_n(7));
        } else if (dvr_style(7) == 2) {
          vr8_dcorr = dvr_alpha(7) / (1 + dvr_beta(7) * dens_n(7));
        } else if (dvr_style(7) == 3) {
          vr8_dcorr = 1 / (1 + exp(dvr_alpha(7) * dens_n(7) + dvr_beta(7)));
        } else if (dvr_style(7) == 4) {
          vr8_dcorr = 1 - (dens_n(7) / dvr_alpha(7));
          if (dvr_beta(7) != 0 && vr8_dcorr < -1.0) vr8_dcorr = 0.0; 
        }
      }
      
      // Juvenile observation
      if (dvr_yn(8)) {
        if (dvr_style(8) == 1) {
          vr9_dcorr = dvr_alpha(8) * exp(-1 * dvr_beta(8) * dens_n(8));
        } else if (dvr_style(8) == 2) {
          vr9_dcorr = dvr_alpha(8) / (1 + dvr_beta(8) * dens_n(8));
        } else if (dvr_style(8) == 3) {
          vr9_dcorr = 1 / (1 + exp(dvr_alpha(8) * dens_n(8) + dvr_beta(8)));
        } else if (dvr_style(8) == 4) {
          vr9_dcorr = 1 - (dens_n(8) / dvr_alpha(8));
          if (dvr_beta(8) != 0 && vr9_dcorr < -1.0) vr9_dcorr = 0.0; 
        }
      }
      
      // Juvenile sizea
      if (dvr_yn(9)) {
        if (dvr_style(9) == 1) {
          vr10_dcorr = dvr_alpha(9) * exp(-1 * dvr_beta(9) * dens_n(9));
        } else if (dvr_style(9) == 2) {
          vr10_dcorr = dvr_alpha(9) / (1 + dvr_beta(9) * dens_n(9));
        } else if (dvr_style(9) == 3) {
          vr10_dcorr = 1 / (1 + exp(dvr_alpha(9) * dens_n(9) + dvr_beta(9)));
        } else if (dvr_style(9) == 4) {
          vr10_dcorr = 1 - (dens_n(9) / dvr_alpha(9));
          if (dvr_beta(9) != 0 && vr10_dcorr < -1.0) vr10_dcorr = 0.0; 
        }
      }
      
      // Juvenile sizeb
      if (dvr_yn(10)) {
        if (dvr_style(10) == 1) {
          vr11_dcorr = dvr_alpha(10) * exp(-1 * dvr_beta(10) * dens_n(10));
        } else if (dvr_style(10) == 2) {
          vr11_dcorr = dvr_alpha(10) / (1 + dvr_beta(10) * dens_n(10));
        } else if (dvr_style(10) == 3) {
          vr11_dcorr = 1 / (1 + exp(dvr_alpha(10) * dens_n(10) + dvr_beta(10)));
        } else if (dvr_style(10) == 4) {
          vr11_dcorr = 1 - (dens_n(10) / dvr_alpha(10));
          if (dvr_beta(10) != 0 && vr11_dcorr < -1.0) vr11_dcorr = 0.0; 
        }
      }
      
      // Juvenile sizec
      if (dvr_yn(11)) {
        if (dvr_style(11) == 1) {
          vr12_dcorr = dvr_alpha(11) * exp(-1 * dvr_beta(11) * dens_n(11));
        } else if (dvr_style(11) == 2) {
          vr12_dcorr = dvr_alpha(11) / (1 + dvr_beta(11) * dens_n(11));
        } else if (dvr_style(11) == 3) {
          vr12_dcorr = 1 / (1 + exp(dvr_alpha(11) * dens_n(11) + dvr_beta(11)));
        } else if (dvr_style(11) == 4) {
          vr12_dcorr = 1 - (dens_n(11) / dvr_alpha(11));
          if (dvr_beta(11) != 0 && vr12_dcorr < -1.0) vr12_dcorr = 0.0; 
        }
      }
      
      // Juvenile reproduction
      if (dvr_yn(12)) {
        if (dvr_style(12) == 1) {
          vr13_dcorr = dvr_alpha(12) * exp(-1 * dvr_beta(12) * dens_n(12));
        } else if (dvr_style(12) == 2) {
          vr13_dcorr = dvr_alpha(12) / (1 + dvr_beta(12) * dens_n(12));
        } else if (dvr_style(12) == 3) {
          vr13_dcorr = 1 / (1 + exp(dvr_alpha(12) * dens_n(12) + dvr_beta(12)));
        } else if (dvr_style(12) == 4) {
          vr13_dcorr = 1 - (dens_n(12) / dvr_alpha(12));
          if (dvr_beta(12) != 0 && vr13_dcorr < -1.0) vr13_dcorr = 0.0; 
        }
      }
      
      // Juvenile maturity
      if (dvr_yn(13)) {
        if (dvr_style(13) == 1) {
          vr14_dcorr = dvr_alpha(13) * exp(-1 * dvr_beta(13) * dens_n(13));
        } else if (dvr_style(13) == 2) {
          vr14_dcorr = dvr_alpha(13) / (1 + dvr_beta(13) * dens_n(13));
        } else if (dvr_style(13) == 3) {
          vr14_dcorr = 1 / (1 + exp(dvr_alpha(13) * dens_n(13) + dvr_beta(13)));
        } else if (dvr_style(13) == 4) {
          vr14_dcorr = 1 - (dens_n(13) / dvr_alpha(13));
          if (dvr_beta(13) != 0 && vr14_dcorr < -1.0) vr14_dcorr = 0.0; 
        }
      }
    }
    
    // Matrix out collects conditional probabilities
    // Zero matrix with n rows and 7 columns: 0 surv, 1 obs, 2 repst,
    // 3 size, 4 size_b, 5 size_c, 6 matst, >6 are test variables
    if (err_check) {
      NumericMatrix zeroform(n, 7);
      std::fill(zeroform.begin(), zeroform.end(), 0.0); // Added for Linux issues
      out = zeroform;
      CharacterVector out_names = {"surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst"};
      colnames(out) = out_names;
    }
    NumericVector out_vec(7, 0.0);
    
    arma::mat survtransmat;
    arma::mat fectransmat;
    arma::sp_mat survtransmat_sp;
    arma::sp_mat fectransmat_sp;
    
    if (!sparse) {
      arma::mat survtransmat_pre(matrixdim, matrixdim, fill::zeros);
      arma::mat fectransmat_pre(matrixdim, matrixdim, fill::zeros);
      
      survtransmat = survtransmat_pre;
      fectransmat = fectransmat_pre;
    } else {
      arma::sp_mat survtransmat_sp_pre(matrixdim, matrixdim);
      arma::sp_mat fectransmat_sp_pre(matrixdim, matrixdim);
      
      survtransmat_sp = survtransmat_sp_pre;
      fectransmat_sp = fectransmat_sp_pre;
    }
    
    double fec_addedcoefs = sum(feccoefs);
    double jsurv_coefsadded = sum(jsurvcoefs);
    double mat_predicted {0.0};
    unsigned int k {0};
    
    // Loop runs through each line of AllStages, calculates each estimable matrix element
    for(int i = 0; i < n; i++) {
      out_vec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
      k = aliveandequal(i);
      mat_predicted = 0.0;
      
      if (err_check) out(i, 6) = 1.0; // Initialization of maturity status probability for typical case
      
      Rcpp::NumericVector statusterms = {fl1(i), fl2n(i), sz1(i), sz2o(i),
        szb1(i), szb2o(i), szc1(i), szc2o(i), actualage2(i), inda1, inda2, indb1,
        indb2, indc1, indc2, dens, sz3(i), szb3(i), szc3(i), binwidth3(i),
        binbwidth3(i), bincwidth3(i), anna1, anna2, annb1, annb2, annc1, annc2};
      
      if (ovgivent(i) == -1 && indata(i) == 1 && stage2n(i) == stage2o(i)) {
        if ((mat2n(i) == 1 && mat3(i) == 1) || (mat2o(i) == 1 && mat3(i) == 1)) {
          
          // Adult survival transitions
          if (survdist < 5) {
            
            out_vec(0) = preouterator(survproxy, survcoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, chosen_f2inda_cat,
              chosen_f1inda_cat, chosen_f2indb_cat, chosen_f1indb_cat,
              chosen_f2indc_cat, chosen_f1indc_cat, statusterms, survgroups2,
              survgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
              survind, survind_rownames, sizeindzi, sizeind_rownames_zi, false,
              survsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 1, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0);
            out_vec(0) = out_vec(0) * vr1_dcorr;
          } else {
            out_vec(0) = survcoefs(0);
            out_vec(0) = out_vec(0) * vr1_dcorr;
          }
          if (err_check) out(i, 0) = out_vec(0);
          
          if (obsdist < 5) {
            out_vec(1) = preouterator(obsproxy, obscoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, chosen_f2inda_cat,
              chosen_f1inda_cat, chosen_f2indb_cat, chosen_f1indb_cat,
              chosen_f2indc_cat, chosen_f1indc_cat, statusterms, obsgroups2,
              obsgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
              obsind, obsind_rownames, sizeindzi, sizeind_rownames_zi, false,
              obssigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 2, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0);
            out_vec(1) = out_vec(1) * vr2_dcorr;
            
          } else {
            out_vec(1) = obscoefs(0);
            out_vec(1) = out_vec(1) * vr2_dcorr;
          }
          if (err_check) out(i, 1) = out_vec(1);
          
          if (ob3(i) == 1 || obsdist == 5) {
            
            if (sizedist < 5) {
              bool used_sizezero = false;
              if (sizezero && sz3(i) == 0) used_sizezero = sizezero;
              
              out_vec(3) = preouterator(sizeproxy, sizecoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, sizegroups2, sizegroups1, sizegroups2zi, sizegroups1zi,
                sizeyearzi, sizepatchzi, sizeind, sizeind_rownames, sizeindzi,
                sizeind_rownames_zi, used_sizezero, sizesigma, grp2o(i), grp1(i),
                patchnumber, yearnumber, sizedist, 3, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
                sizetrunc);
              out_vec(3) = out_vec(3) * vr3_dcorr;
                
            } else {
              out_vec(3) = 1.0;
              out_vec(3) = out_vec(3) * vr3_dcorr;
            }
            if (err_check) out(i, 3) = out_vec(3);
            
            if (sizebdist < 5) {
              bool used_sizebzero = false;
              if (sizebzero && szb3(i) == 0) used_sizebzero = sizebzero;
              
              out_vec(4) = preouterator(sizebproxy, sizebcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, sizebgroups2, sizebgroups1, sizebgroups2zi,
                sizebgroups1zi, sizebyearzi, sizebpatchzi, sizebind,
                sizebind_rownames, sizebindzi, sizebind_rownames_zi, used_sizebzero,
                sizebsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist,
                4, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
                negfec, stage2n(i), nostages, sizebtrunc);
              out_vec(4) = out_vec(4) * vr4_dcorr;
            } else {
              out_vec(4) = 1.0;
              out_vec(4) = out_vec(4) * vr4_dcorr;
            }
            if (err_check) out(i, 4) = out_vec(4);
            
            if (sizecdist < 5) {
              bool used_sizeczero = false;
              if (sizeczero && szc3(i) == 0) used_sizeczero = sizeczero;
              
              out_vec(5) = preouterator(sizecproxy, sizeccoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, sizecgroups2, sizecgroups1, sizecgroups2zi,
                sizecgroups1zi, sizecyearzi, sizecpatchzi, sizecind,
                sizecind_rownames, sizecindzi, sizecind_rownames_zi, used_sizeczero,
                sizecsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist,
                5, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
                negfec, stage2n(i), nostages, sizectrunc);
              out_vec(5) = out_vec(5) * vr5_dcorr;
            } else {
              out_vec(5) = 1.0;
              out_vec(5) = out_vec(5) * vr5_dcorr;
            }
            if (err_check) out(i, 5) = out_vec(5);
            
            if (repstdist < 5) {
              out_vec(2) = preouterator(repstproxy, repstcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, repstgroups2, repstgroups1, dud_groups2zi, dud_groups1zi,
                dud_yearzi, dud_patchzi, repstind, repstind_rownames, sizeindzi,
                sizeind_rownames_zi, false, repstsigma, grp2o(i), grp1(i),
                patchnumber, yearnumber, 4, 6, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages, 0);
              out_vec(2) = out_vec(2) * vr6_dcorr;
                
              if (fl3(i) == 0) {
                out_vec(2) = 1.0 - out_vec(2);
              }
            } else {
              if (fl3(i) == 0) {
                out_vec(2) = repstcoefs(0);
                out_vec(2) = out_vec(2) * vr6_dcorr;
                out_vec(2) = 1.0 - out_vec(2);
              } else if (fl3(i) == 1) {
                out_vec(2) = repstcoefs(0);
                out_vec(2) = out_vec(2) * vr6_dcorr;
              } else {
                out_vec(2) = 0.0;
              }
            }
            if (err_check) out(i, 2) = out_vec(2);
            
          } else {
            out_vec(1) = 1.0 - out_vec(1);
            out_vec(2) = 1.0;
            out_vec(3) = 1.0;
            out_vec(4) = 1.0;
            out_vec(5) = 1.0;
            out_vec(6) = 1.0;
            
            if (err_check) {
              out(i, 1) = out_vec(1);
              out(i, 2) = out_vec(2);
              out(i, 3) = out_vec(3);
              out(i, 4) = out_vec(4);
              out(i, 5) = out_vec(5);
              out(i, 6) = out_vec(6);
            }
          }
          if (!sparse) {
            survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
            out_vec(4) * out_vec(5) * out_vec(6);
          } else {
            survtransmat_sp(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
            out_vec(4) * out_vec(5) * out_vec(6);
          }
          
        } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurv_coefsadded != 0.0) {
          // Juvenile to adult transitions
          if (jmatstdist < 5) {
            mat_predicted = preouterator(jmatstproxy, jmatstcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, jmatstgroups2, jmatstgroups1, dud_groups2zi,
              dud_groups1zi, dud_yearzi, dud_patchzi, jmatstind,
              jmatstind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              jmatstsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 21,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, 0);
            mat_predicted = mat_predicted * vr14_dcorr;
            
            if (mat3(i) > 0.5) {
              out_vec(6) = mat_predicted;
            } else {
              out_vec(6) = 1 - mat_predicted;
            }
          } else {
            if (mat3(i) > 0.5) {
              out_vec(6) = vr14_dcorr;
            } else {
              out_vec(6) = 1 - vr14_dcorr;
            }
          }
          if (err_check) out(i, 6) = out_vec(6);
          
          if (jsurvdist < 5) {
            out_vec(0) = preouterator(jsurvproxy, jsurvcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, jsurvgroups2, jsurvgroups1, dud_groups2zi,
              dud_groups1zi, dud_yearzi, dud_patchzi, jsurvind,
              jsurvind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              jsurvsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 8,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, 0);
            out_vec(0) = out_vec(0) * vr8_dcorr;
          } else {
            out_vec(0) = jsurvcoefs(0);
            out_vec(0) = out_vec(0) * vr8_dcorr;
          }
          if (err_check) out(i, 0) = out_vec(0);
          
          if (jobsdist < 5) {
            out_vec(1) = preouterator(jobsproxy, jobscoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, jobsgroups2, jobsgroups1, dud_groups2zi,
              dud_groups1zi, dud_yearzi, dud_patchzi, jobsind, jobsind_rownames,
              jsizeindzi, jsizeind_rownames_zi, false, jobssigma, grp2o(i),
              grp1(i), patchnumber, yearnumber, 4, 9, exp_tol, theta_tol,
              ipm_cdf, matrixformat, fecmod, repentry(i), negfec, stage2n(i),
              nostages, 0);
            out_vec(1) = out_vec(1) * vr9_dcorr;
          } else {
            out_vec(1) = jobscoefs(0);
            out_vec(1) = out_vec(1) * vr9_dcorr;
          }
          if (err_check) out(i, 1) = out_vec(1);
          
          if (ob3(i) == 1 || jobsdist == 5) {
            if (jsizedist < 5) {
              out_vec(3) = preouterator(jsizeproxy, jsizecoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, jsizegroups2, jsizegroups1, jsizegroups2zi,
                jsizegroups1zi, jsizeyearzi, jsizepatchzi, jsizeind,
                jsizeind_rownames, jsizeindzi, jsizeind_rownames_zi, jsizezero,
                jsizesigma, grp2o(i), grp1(i), patchnumber, yearnumber,
                sizedist, 10, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
                repentry(i), negfec, stage2n(i), nostages, jsizetrunc);
              out_vec(3) = out_vec(3) * vr10_dcorr;
            } else {
              out_vec(3) = 1.0;
              out_vec(3) = out_vec(3) * vr10_dcorr;
            }
            if (err_check) out(i, 3) = out_vec(3);
            
            if (jsizebdist < 5) {
              out_vec(4) = preouterator(jsizebproxy, jsizebcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, jsizebgroups2, jsizebgroups1, jsizebgroups2zi,
                jsizebgroups1zi, jsizebyearzi, jsizebpatchzi, jsizebind,
                jsizebind_rownames, jsizebindzi, jsizebind_rownames_zi,
                jsizebzero, jsizebsigma, grp2o(i), grp1(i), patchnumber,
                yearnumber, sizebdist, 11, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
                jsizebtrunc);
              out_vec(4) = out_vec(4) * vr11_dcorr;
            } else {
              out_vec(4) = 1.0;
              out_vec(4) = out_vec(4) * vr11_dcorr;
            }
            if (err_check) out(i, 4) = out_vec(4);
            
            if (jsizecdist < 5) {
              out_vec(5) = preouterator(jsizecproxy, jsizeccoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda,chosen_r1inda,
                chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
                chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
                chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
                statusterms, jsizecgroups2, jsizecgroups1, jsizecgroups2zi,
                jsizecgroups1zi, jsizecyearzi, jsizecpatchzi, jsizecind,
                jsizecind_rownames, jsizecindzi, jsizecind_rownames_zi,
                jsizeczero, jsizecsigma, grp2o(i), grp1(i), patchnumber,
                yearnumber, sizecdist, 12, exp_tol, theta_tol, ipm_cdf,
                matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
                jsizectrunc);
              out_vec(5) = out_vec(5) * vr12_dcorr;
            } else {
              out_vec(5) = 1.0;
              out_vec(5) = out_vec(5) * vr12_dcorr;
            }
            if (err_check) out(i, 5) = out_vec(5);
            
            if (jrepstdist < 5) {
              out_vec(2) = preouterator(jrepstproxy, jrepstcoefs, rand_index,
                dev_terms, vital_year, vital_patch, chosen_r2inda,
                chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
                chosen_r1indc, chosen_f2inda_cat, chosen_f1inda_cat,
                chosen_f2indb_cat, chosen_f1indb_cat, chosen_f2indc_cat,
                chosen_f1indc_cat, statusterms, jrepstgroups2, jrepstgroups1,
                dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
                jrepstind, jrepstind_rownames, jsizeindzi, jsizeind_rownames_zi,
                false, jrepstsigma, grp2o(i), grp1(i), patchnumber, yearnumber,
                4, 13, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
                repentry(i), negfec, stage2n(i), nostages, 0);
              out_vec(2) = out_vec(2) * vr13_dcorr;
                
              if (fl3(i) == 0) {
                out_vec(2) = 1.0 - out_vec(2);
              }
            } else {
              if (fl3(i) == 0) {
                out_vec(2) = jrepstcoefs(0);
                out_vec(2) = out_vec(2) * vr13_dcorr;
                out_vec(2) = 1.0 - out_vec(2);
              } else if (fl3(i) == 1) {
                out_vec(2) = jrepstcoefs(0);
                out_vec(2) = out_vec(2) * vr13_dcorr;
              } else {
                out_vec(2) = 0.0;
              }
            }
            if (err_check) out(i, 2) = out_vec(2);
            
          } else {
            out_vec(1) = 1.0 - out_vec(1);
            out_vec(2) = 1.0;
            out_vec(3) = 1.0;
            out_vec(4) = 1.0;
            out_vec(5) = 1.0;
            out_vec(6) = 1.0;
            
            if (err_check) {
              out(i, 1) = out_vec(1);
              out(i, 2) = out_vec(2);
              out(i, 3) = out_vec(3);
              out(i, 4) = out_vec(4);
              out(i, 5) = out_vec(5);
              out(i, 6) = out_vec(6);
            }
          }
          
          if (!sparse) {
            survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
              out_vec(4) * out_vec(5) * out_vec(6);
          } else {
            survtransmat_sp(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
              out_vec(4) * out_vec(5) * out_vec(6);
          }
        }
      } else if (ovgivent(i) != -1.0) {
        // All other transitions
        if (!sparse) {
          survtransmat(k) = ovgivent(i);
        } else {
          survtransmat_sp(k) = ovgivent(i);
        }
      }
      
      // Fecundity calculation
      if (indata2n(i) == 1 && fec_addedcoefs != 0.0 && repentry(i) > 0) {
        if (fl2o(i) > 0.0 && ovgivenf(i) == -1.0) {
          
          if (!sparse) {
            fectransmat(k) = preouterator(fecproxy, feccoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, fecgroups2, fecgroups1, fecgroups2zi, fecgroups1zi,
              fecyearzi, fecpatchzi, fecind, fecind_rownames, fecindzi,
              fecind_rownames_zi, feczero, fecsigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, fecdist, 7, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
              fectrunc);
            fectransmat(k) = fectransmat(k) * vr7_dcorr;
          } else {
            fectransmat_sp(k) = preouterator(fecproxy, feccoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              chosen_f2inda_cat, chosen_f1inda_cat, chosen_f2indb_cat,
              chosen_f1indb_cat, chosen_f2indc_cat, chosen_f1indc_cat,
              statusterms, fecgroups2, fecgroups1, fecgroups2zi, fecgroups1zi,
              fecyearzi, fecpatchzi, fecind, fecind_rownames, fecindzi,
              fecind_rownames_zi, feczero, fecsigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, fecdist, 7, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
              fectrunc);
            fectransmat_sp(k) = fectransmat_sp(k) * vr7_dcorr;
          }
          
        } else if (ovgivenf(i) != -1.0) {
          if (!sparse) {
            fectransmat(k) = ovgivenf(i);
            fectransmat(k) = fectransmat(k) * vr7_dcorr;
          } else {
            fectransmat_sp(k) = ovgivenf(i);
            fectransmat_sp(k) = fectransmat_sp(k) * vr7_dcorr;
          }
        }
      } else if (ovgivenf(i) != -1.0) {
        if (!sparse) {
          fectransmat(k) = ovgivenf(i);
          fectransmat(k) = fectransmat(k) * vr7_dcorr;
        } else {
          fectransmat_sp(k) = ovgivenf(i);
          fectransmat_sp(k) = fectransmat_sp(k) * vr7_dcorr;
        }
      }
    }
    
    double ov_mult {0.0};
    if (replacementst > 0) {
      for (int i = 0; i < replacementst; i++) {
        
        repindex = replacetvec(i); // AllStages index
        properindex = aliveandequal(repindex);
        arma::uvec rightindex = find(index321 == ovestt(repindex));
        
        if (rightindex.n_elem > 0) {
          proxyindex = aliveandequal(rightindex(0));
          
          ov_mult = ovsurvmult(repindex);
          if (ov_mult < 0.0) ov_mult = 1.0;
          
          if (!sparse) {
            survtransmat(properindex) = survtransmat(proxyindex) * ov_mult;
          } else {
            survtransmat_sp(properindex) = survtransmat_sp(proxyindex) * ov_mult;
          }
        }
      }
    }
    
    if (replacementsf > 0) {
      for (int i = 0; i < replacementsf; i++) {
        
        repindex = replacefvec(i); // AllStages index
        properindex = aliveandequal(repindex);
        arma::uvec rightindex = find(index321 == ovestf(repindex));
        
        if (rightindex.n_elem > 0) {
          proxyindex = aliveandequal(rightindex(0));
          
          ov_mult = ovfecmult(repindex);
          if (ov_mult < 0.0) ov_mult = 1.0;
          if (!sparse) {
            fectransmat(properindex) = fectransmat(proxyindex) * ov_mult;
          } else {
            fectransmat_sp(properindex) = fectransmat_sp(proxyindex) * ov_mult;
          }
        }
      }
    }
    
    if (offsetst > 0) {
      for(int i = 0; i < n; i++) {
        if (ovoffsett(i) != 0.) {
          k = aliveandequal(i);
          if (!sparse) {
            survtransmat(k) = survtransmat(k) + ovoffsett(i);
          } else {
            survtransmat_sp(k) = survtransmat_sp(k) + ovoffsett(i);
          }
        }
      }
    }
    
    if (offsetsf > 0) {
      for(int i = 0; i < n; i++) {
        if (ovoffsetf(i) != 0.) {
          k = aliveandequal(i);
          if (!sparse) {
            fectransmat(k) = fectransmat(k) + ovoffsetf(i);
          } else {
            fectransmat_sp(k) = fectransmat_sp(k) + ovoffsetf(i);
          }
        }
      }
    }
    
    if (tmults_only_st > 0) {
      for (int i = 0; i < tmults_only_st; i++) {
        repindex = tmults_only(i);
        properindex = aliveandequal(repindex);
        ov_mult = ovsurvmult(repindex);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        if (!sparse) {
          survtransmat(properindex) = survtransmat(properindex) * ov_mult;
        } else {
          survtransmat_sp(properindex) = survtransmat_sp(properindex) * ov_mult;
        }
      }
    }
    
    if (fmults_only_st > 0) {
      for (int i = 0; i < fmults_only_st; i++) {
        repindex = fmults_only(i);
        properindex = aliveandequal(repindex);
        ov_mult = ovfecmult(repindex);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        if (!sparse) {
          fectransmat(properindex) = fectransmat(properindex) * ov_mult;
        } else {
          fectransmat_sp(properindex) = fectransmat_sp(properindex) * ov_mult;
        }
      }
    }
    
    // Final output
    List output_final;
    
    if (!proj_only) {
      List output(4);
      
      if (!sparse) {
        if (!simplicity) {
          arma::mat amatrix = survtransmat + fectransmat;
          output(0) = amatrix;
        } else {
          output(0) = R_NilValue;
        }
        
        output(1) = survtransmat;
        output(2) = fectransmat;
        
        if (err_check) {
          output(3) = out;
        } else {
          output(3) = R_NilValue;
        }
      } else {
        if (!simplicity) {
          arma::sp_mat amatrix = survtransmat_sp + fectransmat_sp;
          output(0) = amatrix;
        } else {
          output(0) = R_NilValue;
        }
        
        output(1) = survtransmat_sp;
        output(2) = fectransmat_sp;
        
        if (err_check) {
          output(3) = out;
        } else {
          output(3) = R_NilValue;
        }
      }
      output_final = output;
      CharacterVector output_names = {"A", "U", "F", "out"};
      output_final.attr("names") = output_names;
    } else {
      List output (2);
      
      if (!sparse) {
        arma::mat amatrix = survtransmat + fectransmat;
        output(0) = amatrix;
        
        if (err_check) {
          output(1) = out;
        } else {
          output(1) = R_NilValue;
        }
      } else {
        arma::sp_mat amatrix = survtransmat_sp + fectransmat_sp;
        output(0) = amatrix;
        
        if (err_check) {
          output(1) = out;
        } else {
          output(1) = R_NilValue;
        }
      }
      output_final = output;
      CharacterVector output_names = {"A", "out"};
      output_final.attr("names") = output_names;
    }
    
    return output_final;
  }
  
  //' Estimate All Elements of Function-based Leslie Projection Matrix
  //' 
  //' Function \code{motherbalowski()} swiftly calculates matrix elements in
  //' function-based Leslie population projection matrices. Used in
  //' \code{mpm_create()}, and through this function in \code{fleslie()}.
  //' 
  //' @name motherbalowski
  //' 
  //' @param actualages An integer vector of all ages to be included in the
  //' matrices, in order.
  //' @param ageframe The modified stageframe used in matrix calculations.
  //' @param survproxy List of coefficients estimated in model of survival.
  //' @param fecproxy List of coefficients estimated in model of fecundity.
  //' @param f2_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_inda A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{a} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{b} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_indc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of individual covariate \code{c} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_anna A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{a} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_anna A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{a} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_annb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{b} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_annb A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{b} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_annc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{c} at
  //' each time \emph{t} to be used in analysis.
  //' @param f1_annc A numeric vector of length equal to the number of years,
  //' holding values equal to the mean value of annual covariate \code{c} at
  //' each time \emph{t}-1 to be used in analysis.
  //' @param f2_inda_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{a} at each time \emph{t}
  //' to be used in analysis.
  //' @param f1_inda_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{a} at each time \emph{t-1}
  //' to be used in analysis.
  //' @param f2_indb_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{b} at each time \emph{t}
  //' to be used in analysis.
  //' @param f1_indb_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{b} at each time \emph{t-1}
  //' to be used in analysis.
  //' @param f2_indc_cat A string vector of length equal to the number of years,
  //' holding categories of individual covariate \code{c} at each time \emph{t}
  //' to be used in analysis.
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
  //' @param surv_dev A numeric value indicating the deviation to the linear
  //' model of survival input by the user.
  //' @param fec_dev A numeric value indicating the deviation to the linear
  //' model of fecundity input by the user.
  //' @param dens A numeric value equal to the density to be used in calculations.
  //' @param fecmod A scalar multiplier for fecundity.
  //' @param finalage The final age to be included in Leslie MPM estimation.
  //' @param negfec A logical value denoting whether to change negative estimated
  //' fecundity to \code{0}.
  //' @param yearnumber An integer specifying which time at time \emph{t} to
  //' develop matrices for. Must be in reference to the \code{listofyears} object
  //' developed in the \code{R} matrix estimator function.
  //' @param patchnumber An integer specifying which patch to develop matrices
  //' for. Must be in reference to the \code{listofyears} object developed in the
  //' \code{R} matrix estimator function.
  //' @param dens_vr A logical value indicating whether any vital rates are
  //' density dependent.
  //' @param dvr_yn A logical vector indicating whether each vital rate is density
  //' dependent.
  //' @param dvr_style An integer vector indicating the style of density
  //' dependence for each vital rate.
  //' @param dvr_alpha A numeric vector indicating the value of alpha to use in
  //' density dependence for each vital rate.
  //' @param dvr_beta A numeric vector indicating the value of beta to use in
  //' density dependence for each vital rate.
  //' @param dens_n A numeric vector corresponding to the population size to use
  //' in vital rate density dependence calculations.
  //' @param exp_tol A numeric value indicating the maximum limit for the
  //' \code{exp()} function to be used in vital rate calculations. Defaults to
  //' \code{700.0}.
  //' @param theta_tol A numeric value indicating a maximum value for theta in
  //' negative binomial probability density estimation. Defaults to
  //' \code{100000000.0}.
  //' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
  //' \code{F}, rather than also outputting matrix \code{A}. Defaults to
  //' \code{FALSE}.
  //' @param sparse If \code{TRUE}, then only outputs matrices in sparse format.
  //' Defaults to \code{FALSE}.
  //' @param supplement An optional data frame edited and age-expanded showing
  //' supplemental transition information.
  //' 
  //' @return A list of 3 matrices, including the main MPM (A), the survival-
  //' transition matrix (U), and a fecundity matrix (F).
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List motherbalowski(IntegerVector actualages,
    const DataFrame& ageframe, const List& survproxy, const List& fecproxy,
    NumericVector f2_inda, NumericVector f1_inda, NumericVector f2_indb,
    NumericVector f1_indb, NumericVector f2_indc, NumericVector f1_indc,
    NumericVector f2_anna, NumericVector f1_anna, NumericVector f2_annb,
    NumericVector f1_annb, NumericVector f2_annc, NumericVector f1_annc,
    StringVector f2_inda_cat, StringVector f1_inda_cat, StringVector f2_indb_cat,
    StringVector f1_indb_cat, StringVector f2_indc_cat, StringVector f1_indc_cat,
    StringVector r2_inda, StringVector r1_inda, StringVector r2_indb,
    StringVector r1_indb, StringVector r2_indc, StringVector r1_indc,
    double surv_dev, double fec_dev, double dens, double fecmod,
    unsigned int finalage, bool negfec, int yearnumber, int patchnumber,
    bool dens_vr, LogicalVector dvr_yn, IntegerVector dvr_style,
    NumericVector dvr_alpha, NumericVector dvr_beta, NumericVector dens_n,
    double exp_tol = 700.0, double theta_tol = 100000000.0,
    bool simplicity = false, bool sparse = false,
    Nullable<DataFrame> supplement = R_NilValue) {
    
    // Determines the size of the matrix
    StringVector sf_agenames = as<StringVector>(ageframe["stage"]);
    IntegerVector sf_minage = as<IntegerVector>(ageframe["min_age"]);
    IntegerVector sf_maxage = as<IntegerVector>(ageframe["max_age"]);
    IntegerVector sf_repstatus = as<IntegerVector>(ageframe["repstatus"]);
    int noages = actualages.length();
    int start_age = min(actualages);
    //int last_age = max(actualages);
    
    bool cont = false;
    if (sf_maxage(noages - 1) == NA_INTEGER) {
      cont = true;
    }
    
    // Supplement processing
    IntegerVector ov_age2;
    IntegerVector ov_estage2;
    NumericVector ov_givenrate;
    NumericVector ov_offset;
    NumericVector ov_multiplier;
    IntegerVector ov_convtype;
    int supp_length {0};
    
    DataFrame supplement_;
    if (supplement.isNotNull()){
      supplement_ = as<DataFrame>(supplement);
      
      if (supplement_.containsElementNamed("age2")) {
        ov_age2 = clone(as<IntegerVector>(supplement_["age2"]));
        ov_estage2 = clone(as<IntegerVector>(supplement_["estage2"]));
        ov_givenrate = as<NumericVector>(supplement_["givenrate"]);
        ov_offset = as<NumericVector>(supplement_["offset"]);
        ov_multiplier = as<NumericVector>(supplement_["multiplier"]);
        ov_convtype = as<IntegerVector>(supplement_["convtype"]);
        
        supp_length = static_cast<int>(ov_givenrate.length());
        
        for (int i = 0; i < supp_length; i++) {
          ov_age2(i) = ov_age2(i) - start_age;
          
          if (!IntegerVector::is_na(ov_estage2(i))) {
            ov_estage2(i) = ov_estage2(i) - start_age;
          }
        }
      }
    }
    
    // Proxy model imports and settings
    NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
    NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
    
    bool feczero = as<bool>(fecproxy["zero_inflated"]);
    int survdist = as<int>(survproxy["dist"]);
    int fecdist = as<int>(fecproxy["dist"]);
    double fecsigma = as<double>(fecproxy["sigma"]);
    
    if (NumericVector::is_na(fecsigma)) {
      if (fecdist == 1) {
        fecsigma = 1.0;
      } else {
        fecsigma = 0.0;
      }
    }
  
    NumericMatrix vital_year = revelations_leslie(survproxy, fecproxy, 1);
    NumericMatrix vital_patch = revelations_leslie(survproxy, fecproxy, 2);
    
    NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
    NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
    NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
    NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
    NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
    NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
    NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
    NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
    
    NumericVector survind = flightoficarus(survproxy);
    NumericVector fecind = flightoficarus(fecproxy);
    NumericVector fecindzi = zero_flightoficarus(fecproxy);
    
    arma::imat rand_index = foi_index_leslie(survproxy, fecproxy);
    
    StringVector survind_rownames = bootson(survproxy);
    StringVector fecind_rownames = bootson(fecproxy);
    StringVector fecind_rownames_zi = zero_bootson(fecproxy);
    
    // Determination of choices of fixed and random ind covs, and annual covs
    double inda1 = f1_inda(yearnumber);
    double indb1 = f1_indb(yearnumber);
    double indc1 = f1_indc(yearnumber);
    double inda2 = f2_inda(yearnumber);
    double indb2 = f2_indb(yearnumber);
    double indc2 = f2_indc(yearnumber);
    
    String chosen_f2inda_cat = f2_inda_cat(yearnumber);
    String chosen_f1inda_cat = f1_inda_cat(yearnumber);
    String chosen_f2indb_cat = f2_indb_cat(yearnumber);
    String chosen_f1indb_cat = f1_indb_cat(yearnumber);
    String chosen_f2indc_cat = f2_indc_cat(yearnumber);
    String chosen_f1indc_cat = f1_indc_cat(yearnumber);
    
    String chosen_r2inda = r2_inda(yearnumber);
    String chosen_r1inda = r1_inda(yearnumber);
    String chosen_r2indb = r2_indb(yearnumber);
    String chosen_r1indb = r1_indb(yearnumber);
    String chosen_r2indc = r2_indc(yearnumber);
    String chosen_r1indc = r1_indc(yearnumber);
    
    double anna1 = f1_anna(yearnumber);
    double annb1 = f1_annb(yearnumber);
    double annc1 = f1_annc(yearnumber);
    double anna2 = f2_anna(yearnumber);
    double annb2 = f2_annb(yearnumber);
    double annc2 = f2_annc(yearnumber);
    
    // The output matrices
    arma::mat survtransmat;
    arma::mat fectransmat;
    arma::sp_mat survtransmat_sp;
    arma::sp_mat fectransmat_sp;
    
    if (!sparse) {
      arma::mat survtransmat_chuck (noages, noages, fill::zeros);
      arma::mat fectransmat_chuck (noages, noages, fill::zeros);
      
      survtransmat = survtransmat_chuck;
      fectransmat = fectransmat_chuck;
    } else { 
      arma::sp_mat survtransmat_chuck (noages, noages);
      arma::sp_mat fectransmat_chuck (noages, noages);
      
      survtransmat_sp = survtransmat_chuck;
      fectransmat_sp = fectransmat_chuck;
    }
    
    // Following loop runs through each age, and runs through
    // each estimable element in the matrix
    double fec_addedcoefs = sum(feccoefs);
    for(int i = 0; i < noages; i++) {
      // Adult survival transitions
      
      double preout {0.0};
      
      if (survdist < 5) {
        
        // Random covariate processing
        double chosen_randcova2 {0.0};
        if (chosen_r2inda != "none") {
          for (int indcount = 0; indcount < rand_index(0, 0); indcount++) {
            if (chosen_r2inda == survind_rownames(indcount)) {
              chosen_randcova2 = survind(indcount);
            }
          }
        }
        double chosen_randcova1 {0.0};
        if (chosen_r1inda != "none") {
          int delectable_sum = rand_index(0, 0);
          for (int indcount = 0; indcount < rand_index(1, 0); indcount++) {
            if (chosen_r1inda == survind_rownames(indcount + delectable_sum)) {
              chosen_randcova1 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb2 {0.0};
        if (chosen_r2indb != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0);
          for (int indcount = 0; indcount < rand_index(2, 0); indcount++) {
            if (chosen_r2indb == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovb2 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb1 {0.0};
        if (chosen_r1indb != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0);
          for (int indcount = 0; indcount < rand_index(3, 0); indcount++) {
            if (chosen_r1indb == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovb1 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc2 {0.0};
        if (chosen_r2indc != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
            rand_index(3, 0);
          for (int indcount = 0; indcount < rand_index(4, 0); indcount++) {
            if (chosen_r2indc == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovc2 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc1 {0.0};
        if (chosen_r1indc != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
            rand_index(3, 0) + rand_index(4, 0);
          for (int indcount = 0; indcount < rand_index(5, 0); indcount++) {
            if (chosen_r1indc == survind_rownames(indcount + delectable_sum)) {
              chosen_randcovc1 = survind(indcount + delectable_sum);
            }
          }
        }
        
        // Fixed factor covariate processing
        double chosen_fixcova2 {0.0};
        if (chosen_f2inda_cat != "none") {
          for (int indcount = 0; indcount < rand_index(0, 0); indcount++) {
            if (chosen_f2inda_cat == survind_rownames(indcount)) {
              chosen_fixcova2 = survind(indcount);
            }
          }
        }
        double chosen_fixcova1 {0.0};
        if (chosen_f1inda_cat != "none") {
          int delectable_sum = rand_index(0, 0);
          for (int indcount = 0; indcount < rand_index(1, 0); indcount++) {
            if (chosen_f1inda_cat == survind_rownames(indcount + delectable_sum)) {
              chosen_fixcova1 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_fixcovb2 {0.0};
        if (chosen_f2indb_cat != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0);
          for (int indcount = 0; indcount < rand_index(2, 0); indcount++) {
            if (chosen_f2indb_cat == survind_rownames(indcount + delectable_sum)) {
              chosen_fixcovb2 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_fixcovb1 {0.0};
        if (chosen_f1indb_cat != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0);
          for (int indcount = 0; indcount < rand_index(3, 0); indcount++) {
            if (chosen_f1indb_cat == survind_rownames(indcount + delectable_sum)) {
              chosen_fixcovb1 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_fixcovc2 {0.0};
        if (chosen_f2indc_cat != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
            rand_index(3, 0);
          for (int indcount = 0; indcount < rand_index(4, 0); indcount++) {
            if (chosen_f2indc_cat == survind_rownames(indcount + delectable_sum)) {
              chosen_fixcovc2 = survind(indcount + delectable_sum);
            }
          }
        }
        double chosen_fixcovc1 {0.0};
        if (chosen_f1indc_cat != "none") {
          int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
            rand_index(3, 0) + rand_index(4, 0);
          for (int indcount = 0; indcount < rand_index(5, 0); indcount++) {
            if (chosen_f1indc_cat == survind_rownames(indcount + delectable_sum)) {
              chosen_fixcovc1 = survind(indcount + delectable_sum);
            }
          }
        }
        
        double mainsum = rimeotam(survcoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, static_cast<double>(actualages(i)), inda1, inda2, indb1, indb2,
          indc1, indc2, anna1, anna2, annb1, annb2, annc1, annc2, dens, false);
        
        preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
          chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
          chosen_randcovc1 + chosen_fixcova2 + chosen_fixcova1 +
          chosen_fixcovb2 + chosen_fixcovb1 + chosen_fixcovc2 +
          chosen_fixcovc1 + survgroups2(0) + survgroups1(0) + 
          vital_patch(patchnumber, 0) + vital_year(yearnumber, 0) + surv_dev);
  
        if (preout > exp_tol) preout = exp_tol; // Catches numbers too high
        if (i < (noages - 1)) {
          if (!sparse) {
            survtransmat(i+1, i) = exp(preout) / (1.0 + exp(preout));
          } else {
            survtransmat_sp(i+1, i) = exp(preout) / (1.0 + exp(preout));
          }
        } else {
          if (cont) {
            if (!sparse) {
              survtransmat(i, i) = exp(preout) / (1.0 + exp(preout));
            } else {
              survtransmat_sp(i, i) = exp(preout) / (1.0 + exp(preout));
            }
          }
        }
      } else {
        if (i < (noages - 1)) {
          if (!sparse) {
            survtransmat(i+1, i) = survcoefs(0);
          } else {
            survtransmat_sp(i+1, i) = survcoefs(0);
          }
        } else {
          if (!sparse) {
            survtransmat(i, i) = survcoefs(0);
          } else {
            survtransmat_sp(i, i) = survcoefs(0);
          }
        }
      }
      
      // This next block calculates fecundity
      if (fec_addedcoefs != 0.0) {
        if (sf_repstatus(i) == 1) {
          
          // Random covariate processing
          double chosen_randcova2 {0.0};
          if (chosen_r2inda != "none") {
            for (int indcount = 0; indcount < rand_index(0, 6); indcount++) {
              if (chosen_r2inda == fecind_rownames(indcount)) {
                chosen_randcova2 = fecind(indcount);
              }
            }
          }
          double chosen_randcova1 {0.0};
          if (chosen_r1inda != "none") {
            int delectable_sum = rand_index(0, 6);
            for (int indcount = 0; indcount < rand_index(1, 6); indcount++) {
              if (chosen_r1inda == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcova1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb2 {0.0};
          if (chosen_r2indb != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6);
            for (int indcount = 0; indcount < rand_index(2, 6); indcount++) {
              if (chosen_r2indb == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovb2 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb1 {0.0};
          if (chosen_r1indb != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6);
            for (int indcount = 0; indcount < rand_index(3, 6); indcount++) {
              if (chosen_r1indb == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovb1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc2 {0.0};
          if (chosen_r2indc != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
              rand_index(3, 6);
            for (int indcount = 0; indcount < rand_index(4, 6); indcount++) {
              if (chosen_r2indc == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovc2 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc1 {0.0};
          if (chosen_r1indc != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
              rand_index(3, 6) + rand_index(4, 6);
            for (int indcount = 0; indcount < rand_index(5, 6); indcount++) {
              if (chosen_r1indc == fecind_rownames(indcount + delectable_sum)) {
                chosen_randcovc1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcova2zi {0.0};
          if (chosen_r2inda != "none") {
            for (int indcount = 0; indcount < rand_index(0, 16); indcount++) {
              if (chosen_r2inda == fecind_rownames_zi(indcount)) {
                chosen_randcova2zi = fecindzi(indcount);
              }
            }
          }
          double chosen_randcova1zi {0.0};
          if (chosen_r1inda != "none") {
            int delectable_sum = rand_index(0, 16);
            for (int indcount = 0; indcount < rand_index(1, 16); indcount++) {
              if (chosen_r1inda == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcova1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb2zi {0.0};
          if (chosen_r2indb != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16);
            for (int indcount = 0; indcount < rand_index(2, 16); indcount++) {
              if (chosen_r2indb == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovb2zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovb1zi {0.0};
          if (chosen_r1indb != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16);
            for (int indcount = 0; indcount < rand_index(3, 16); indcount++) {
              if (chosen_r1indb == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovb1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc2zi {0.0};
          if (chosen_r2indc != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
              rand_index(3, 16);
            for (int indcount = 0; indcount < rand_index(4, 16); indcount++) {
              if (chosen_r2indc == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovc2zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_randcovc1zi {0.0};
          if (chosen_r1indc != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
              rand_index(3, 16) + rand_index(4, 16);
            for (int indcount = 0; indcount < rand_index(5, 16); indcount++) {
              if (chosen_r1indc == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_randcovc1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
              
          // Fixed factor covariate processing
          double chosen_fixcova2 {0.0};
          if (chosen_f2inda_cat != "none") {
            for (int indcount = 0; indcount < rand_index(0, 6); indcount++) {
              if (chosen_f2inda_cat == fecind_rownames(indcount)) {
                chosen_fixcova2 = fecind(indcount);
              }
            }
          }
          double chosen_fixcova1 {0.0};
          if (chosen_f1inda_cat != "none") {
            int delectable_sum = rand_index(0, 6);
            for (int indcount = 0; indcount < rand_index(1, 6); indcount++) {
              if (chosen_f1inda_cat == fecind_rownames(indcount + delectable_sum)) {
                chosen_fixcova1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovb2 {0.0};
          if (chosen_f2indb_cat != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6);
            for (int indcount = 0; indcount < rand_index(2, 6); indcount++) {
              if (chosen_f2indb_cat == fecind_rownames(indcount + delectable_sum)) {
                chosen_fixcovb2 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovb1 {0.0};
          if (chosen_f1indb_cat != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6);
            for (int indcount = 0; indcount < rand_index(3, 6); indcount++) {
              if (chosen_f1indb_cat == fecind_rownames(indcount + delectable_sum)) {
                chosen_fixcovb1 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovc2 {0.0};
          if (chosen_f2indc_cat != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
              rand_index(3, 6);
            for (int indcount = 0; indcount < rand_index(4, 6); indcount++) {
              if (chosen_f2indc_cat == fecind_rownames(indcount + delectable_sum)) {
                chosen_fixcovc2 = fecind(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovc1 {0.0};
          if (chosen_f1indc_cat != "none") {
            int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
              rand_index(3, 6) + rand_index(4, 6);
            for (int indcount = 0; indcount < rand_index(5, 6); indcount++) {
              if (chosen_f1indc_cat == fecind_rownames(indcount + delectable_sum)) {
                chosen_fixcovc1 = fecind(indcount + delectable_sum);
              }
            }
          }
          
          double chosen_fixcova2zi {0.0};
          if (chosen_f2inda_cat != "none") {
            for (int indcount = 0; indcount < rand_index(0, 16); indcount++) {
              if (chosen_f2inda_cat == fecind_rownames_zi(indcount)) {
                chosen_fixcova2zi = fecindzi(indcount);
              }
            }
          }
          double chosen_fixcova1zi {0.0};
          if (chosen_f1inda_cat != "none") {
            int delectable_sum = rand_index(0, 16);
            for (int indcount = 0; indcount < rand_index(1, 16); indcount++) {
              if (chosen_f1inda_cat == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_fixcova1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovb2zi {0.0};
          if (chosen_f2indb_cat != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16);
            for (int indcount = 0; indcount < rand_index(2, 16); indcount++) {
              if (chosen_f2indb_cat == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_fixcovb2zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovb1zi {0.0};
          if (chosen_f1indb_cat != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16);
            for (int indcount = 0; indcount < rand_index(3, 16); indcount++) {
              if (chosen_f1indb_cat == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_fixcovb1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovc2zi {0.0};
          if (chosen_f2indc_cat != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
              rand_index(3, 16);
            for (int indcount = 0; indcount < rand_index(4, 16); indcount++) {
              if (chosen_f2indc_cat == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_fixcovc2zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          double chosen_fixcovc1zi {0.0};
          if (chosen_f1indc_cat != "none") {
            int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
              rand_index(3, 16) + rand_index(4, 16);
            for (int indcount = 0; indcount < rand_index(5, 16); indcount++) {
              if (chosen_f1indc_cat == fecind_rownames_zi(indcount + delectable_sum)) {
                chosen_fixcovc1zi = fecindzi(indcount + delectable_sum);
              }
            }
          }
          
          double preoutx {0.0};
          
          if (fecdist < 4) {
            if (feczero) {
              
              double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, static_cast<double>(actualages(i)), inda1, inda2,
                indb1, indb2, indc1, indc2, anna1, anna2, annb1, annb2, annc1,
                annc2, dens, true);
              
              preoutx = (mainsum + chosen_randcova2zi + chosen_randcova1zi +
                chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
                chosen_randcovc1zi + chosen_fixcova2zi + chosen_fixcova1zi +
                chosen_fixcovb2zi + chosen_fixcovb1zi + chosen_fixcovc2zi +
                chosen_fixcovc1zi + fecgroups2zi(0) + fecgroups1zi(0) + 
                fecpatchzi(patchnumber) + fecyearzi(yearnumber) + fec_dev);
              
            } else {
              
              double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, static_cast<double>(actualages(i)), inda1, inda2,
                indb1, indb2, indc1, indc2, anna1, anna2, annb1, annb2, annc1,
                annc2, dens, false);
              
              preoutx = (mainsum + chosen_randcova2 + chosen_randcova1 +
                chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
                chosen_randcovc1 + chosen_fixcova2 + chosen_fixcova1 +
                chosen_fixcovb2 + chosen_fixcovb1 + chosen_fixcovc2 +
                chosen_fixcovc1 + fecgroups2(0) + fecgroups1(0) + 
                vital_patch(patchnumber, 1) + vital_year(yearnumber, 1) +
                fec_dev);
            }
            
            if (fecdist == 0 || fecdist == 1) {
              // Poisson and negative binomial fecundity
              
              if (feczero) {
                if (preoutx > exp_tol) preoutx = exp_tol;
                if (!sparse) {
                  fectransmat(0, i) = (exp(preoutx) / (1.0 + exp(preoutx))) * fecmod;
                } else {
                  fectransmat_sp(0, i) = (exp(preoutx) / (1.0 + exp(preoutx))) * fecmod;
                }
              } else {
                if (preoutx > exp_tol) preoutx = exp_tol;
                
                if (!sparse) {
                  fectransmat(0, i) = exp(preoutx) * fecmod;
                } else {
                  fectransmat_sp(0, i) = exp(preoutx) * fecmod;
                }
              }
            } else if (fecdist == 2) {
              // Gaussian fecundity
              if (!sparse) {
                fectransmat(0, i) = preoutx * fecmod;
              
                if (negfec && fectransmat(0, i) < 0.0) {
                  fectransmat(0, i) = 0.0;
                }
              } else { 
                fectransmat_sp(0, i) = preoutx * fecmod;
              
                if (negfec && fectransmat_sp(0, i) < 0.0) {
                  fectransmat_sp(0, i) = 0.0;
                }
              }
            } else if (fecdist == 3) {
              // Gamma fecundity
              if (!sparse) {
                fectransmat(0, i) = (1.0 / preoutx) * fecmod;
              } else {
                fectransmat_sp(0, i) = (1.0 / preoutx) * fecmod;
              }
            }
            
          } else {
            if (!sparse) {
              fectransmat(0, i) = feccoefs(0);
            } else {
              fectransmat_sp(0, i) = feccoefs(0);
            }
          }
        }
      }
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
      } else if (target_col >= (noages - 1)) {
        throw Rcpp::exception("Some age2 values are too high.", false);
      }
      
      if (ov_convtype(l) == 1) {
        if (target_col >= (noages - 1) && cont) {
          target_row = target_col;
        } else {
          target_row = target_col + 1;
        }
        
        if (!NumericVector::is_na(ov_givenrate(l))) {
          if (!sparse) {
            survtransmat(target_row, target_col) = ov_givenrate(l);
          } else {
            survtransmat_sp(target_row, target_col) = ov_givenrate(l);
          }
        }
        
        if (!NumericVector::is_na(ov_offset(l))) {
          if (ov_offset(l) != 0.) {
            if (!sparse) {
              survtransmat(target_row, target_col) =
                survtransmat(target_row, target_col) + ov_offset(l);
            } else {
              survtransmat_sp(target_row, target_col) =
                survtransmat_sp(target_row, target_col) + ov_offset(l);
            }
          }
        }
        
        if (!IntegerVector::is_na(ov_estage2(l))) {
          proxy_col = ov_estage2(l);
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_col = noages - 1;
          } else if (proxy_col >= (noages - 1)) {
            throw Rcpp::exception("Some estage2 values are too high.", false);
          }
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_row = proxy_col; // Was target_col
          } else {
            proxy_row = proxy_col + 1;
          }
          
          if (!sparse) {
            survtransmat(target_row, target_col) = survtransmat(proxy_row, proxy_col);
          } else {
            survtransmat_sp(target_row, target_col) = survtransmat_sp(proxy_row, proxy_col);
          }
        }
        
        if (!NumericVector::is_na(ov_multiplier(l))) {
          if (!sparse) {
            survtransmat(target_row, target_col) *= ov_multiplier(l);
          } else {
            survtransmat_sp(target_row, target_col) *= ov_multiplier(l);
          }
        }
      } else if (ov_convtype(l) == 2) {
        target_row = 0;
        
        if (!NumericVector::is_na(ov_givenrate(l))) {
          if (!sparse) {
            fectransmat(target_row, target_col) = ov_givenrate(l);
          } else {
            fectransmat_sp(target_row, target_col) = ov_givenrate(l);
          }
        }
        
        if (!NumericVector::is_na(ov_offset(l))) {
          if (ov_offset(l) != 0.) {
            if (!sparse) {
              fectransmat(target_row, target_col) =
                fectransmat(target_row, target_col) + ov_offset(l);
            } else {
              fectransmat_sp(target_row, target_col) =
                fectransmat_sp(target_row, target_col) + ov_offset(l);
            }
          }
        }
        
        if (!IntegerVector::is_na(ov_estage2(l))) {
          proxy_col = ov_estage2(l);
          
          if (proxy_col >= (noages - 1) && cont) {
            proxy_col = noages - 1;
          } else if (proxy_col >= (noages - 1)) {
            throw Rcpp::exception("Some estage2 values are too high.", false);
          }
          
          proxy_row = 0;
          
          if (!sparse) {
            fectransmat(target_row, target_col) = fectransmat(proxy_row, proxy_col);
          } else {
            fectransmat_sp(target_row, target_col) = fectransmat_sp(proxy_row, proxy_col);
          }
        }
        
        if (!NumericVector::is_na(ov_multiplier(l))) {
          if (!sparse) {
            fectransmat(target_row, target_col) *= ov_multiplier(l);
          } else {
            fectransmat_sp(target_row, target_col) *= ov_multiplier(l);
          }
        }
      } else {
        target_row = 0;
        
        if (!NumericVector::is_na(ov_multiplier(l))) {
          if (!sparse) {
            fectransmat(target_row, target_col) *= ov_multiplier(l);
          } else {
            fectransmat_sp(target_row, target_col) *= ov_multiplier(l);
          }
        }
      }
    }
    
    List output;
    
    if (simplicity) {
      if (!sparse) {
        output = List::create(_["U"] = survtransmat, _["F"] = fectransmat);
      } else { 
        output = List::create(_["U"] = survtransmat_sp, _["F"] = fectransmat_sp);
      }
    } else {
      if (!sparse) {
        arma::mat amatrix = survtransmat + fectransmat;
        output = List::create(Named("A") = amatrix, _["U"] = survtransmat,
          _["F"] = fectransmat);
      } else {
        arma::sp_mat amatrix_sp = survtransmat_sp + fectransmat_sp;
        output = List::create(Named("A") = amatrix_sp, _["U"] = survtransmat_sp,
          _["F"] = fectransmat_sp);
      }
    }
    
    return output;
  }
  
  //' Converts Labels Element to LOY Data Frame
  //' 
  //' Function \code{loy_inator()} takes a \code{labels} element from a
  //' \code{lefkoMat} object as input, and produces an expanded version called
  //' a 'loy table', which is a data frame with multiple types of vector
  //' corresponding to the original input, and is used in indexing in various
  //' functions in \code{lefko3}.
  //' 
  //' @name loy_inator
  //' 
  //' @param labels A data frame with any two or all three of the \code{pop},
  //' \code{patch}, and \code{year2} string vectors used in \code{lefkoMat}
  //' objects.
  //' @param year2_trigger A logical value indicating whether to throw a warning
  //' if no \code{year2} vector exists within \code{labels}.
  //' 
  //' @return An expanded data frame including seven vectors.
  //' \item{pop}{The \code{pop} string vector in the input \code{labels} data
  //' frame. Denotes population identity.}
  //' \item{patch}{The \code{patch} string vector in the input \code{labels}
  //' data frame. Denotes patch identity.}
  //' \item{year2}{The \code{patch} string vector in the input \code{labels}
  //' data frame. Denotes the time in time \emph{t}.}
  //' \item{poppatch}{A string vector concatenating the identity of the
  //' associated \code{pop} and \code{patch} elements. Denotes the unique
  //' patch (technically pop-patch) identity.}
  //' \item{popc}{An integer vector denoting the population identity.}
  //' \item{poppatchc}{An integer vector denoting the unique patch (technically
  //' pop-patch) identity.}
  //' \item{year2c}{An integer vector denoting the unique identity of the time
  //' in time \emph{t}.}
  //' \item{patchesinpop}{}
  //' \item{yearsinpatch}{}
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::DataFrame loy_inator (DataFrame labels, bool year2_trigger) {
    
    StringVector poporder;
    StringVector patchorder;
    StringVector yearorder;
    
    if (!labels.hasAttribute("names")) {
      throw Rcpp::exception("This lefkoMat object lacks variable names in element labels.",
        false);
    }
    
    StringVector labels_vars = as<StringVector>(labels.attr("names"));
    int labels_novars = labels_vars.length();
    int loysize = labels.nrows();
    
    int found_vars {0};
    bool found_pop {false};
    bool found_patch {false};
    bool found_year2 {false};
    
    for (int i = 0; i < labels_novars; i++) {
      if (stringcompare_hard(as<std::string>(labels_vars(i)), "pop")) {
        poporder = as<StringVector>(labels["pop"]);
        found_pop = true;
        found_vars++;
      }
      if (stringcompare_hard(as<std::string>(labels_vars(i)), "patch")) {
        patchorder = as<StringVector>(labels["patch"]);
        found_patch = true;
        found_vars++;
      }
      if (stringcompare_hard(as<std::string>(labels_vars(i)), "year2")) {
        yearorder = as<StringVector>(labels["year2"]);
        found_year2 = true;
        found_vars++;
      }
    }
    
    if (year2_trigger && !found_year2) {
      Rf_warningcall(R_NilValue, "This lefkoMat object lacks annual matrices.\n");
    }
    if (found_vars < 2) {
      throw Rcpp::exception("Unusual labels element missing pop, patch, and/or year2 vectors.",
        false);
    }
    
    if (!found_pop) {
      poporder = rep("1", loysize);
    }
    if (!found_patch) {
      patchorder = rep("1", loysize);
    }
    if (!found_year2) {
      yearorder = rep("1", loysize);
    }
    
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
    
    return listofyears;
  }
  
  //' Reduces Matrices In A Function-based lefkoMat Object
  //' 
  //' This function reduces all of the matrices in a function-based
  //' \code{lefkoMat} object, looking for the common rows and columns that lack
  //' non-zero elements among all A matrices. It also calculates the matrix qc
  //' vector.
  //' 
  //' @name matrix_reducer
  //' 
  //' @param mat_list A list of matrices. In function-based encoding, there are
  //' three main elements, \code{A}, \code{U}, and \code{F}, each containing as
  //' many matrices as entries in the list-of-years.
  //' each of which is a list of square matrices of equal dimension.
  //' @param qc A three-element integer vector.
  //' @param ahstages The stageframe as processed by \code{sf_reassess} and
  //' the matrix creation function utilized, or the equivalent for Leslie MPMs.
  //' @param hstages The historical stage-pair index data frame.
  //' @param agestages The age-stage index data frame.
  //' @param age A logical value indicating whether the MPM is age-based.
  //' Defaults to \code{FALSE}.
  //' @param stage A logical value indicating whether the MPM is stage-based.
  //' Defaults to \code{TRUE}.
  //' @param historical A logical value indicating whether the MPM is
  //' historical. Defaults to \code{FALSE}.
  //' @param reduce A logical value indicating whether to reduce the matrices.
  //' Defaults to \code{FALSE}.
  //' @param simplicity A logical value indicating whether the MPM includes
  //' \code{A}, \code{U}, and \code{F} matrices, or only the latter. Defaults to
  //' \code{FALSE}, in which case all three are included.
  //' @param sparse_out A logical value indicating whether used matrices are in
  //' sparse format. Defaults to \code{FALSE}.
  //' 
  //' @return Reduces the matrices via reference, and creates a \code{matrixqc}
  //' object.
  //' 
  //' @keywords internal
  //' @noRd
  inline void matrix_reducer (List& mat_list, IntegerVector& qc,
    DataFrame& ahstages, DataFrame& hstages, DataFrame& agestages,
    bool age = false, bool stage = true, bool historical = false,
    bool reduce = false, bool simplicity = false, bool sparse_out = false) {
    
    List imported_A;
    if (!simplicity) imported_A = as<List>(mat_list["A"]);
    List imported_U = as<List>(mat_list["U"]);
    List imported_F = as<List>(mat_list["F"]);
    
    int no_matrices = imported_U.length();
    arma::uvec bits_to_remove;
    
    int u_transitions {0};
    int f_transitions {0};
    
    if (reduce) {
      if (!sparse_out) {
        for (int i = 0; i < no_matrices; i++) {
          arma::mat A_arma;
          
          if (simplicity) {
            if (!is<NumericMatrix>(imported_U(i))) continue;
            
            arma::mat U_arma = as<arma::mat>(imported_U(i));
            arma::mat F_arma = as<arma::mat>(imported_F(i));
            
            A_arma = U_arma + F_arma;
          } else {
            if (!is<NumericMatrix>(imported_U(i))) continue;
            A_arma = as<arma::mat>(imported_A(i));
          }
          
          arma::rowvec A_colsums = sum(A_arma, 0);
          arma::vec A_rowsums = sum(A_arma, 1);
          
          arma::uvec A_colzeros = find(A_colsums == 0.0);
          arma::uvec A_rowzeros = find(A_rowsums == 0.0);
          arma::uvec A_allzeros = intersect(A_colzeros, A_rowzeros);
          
          if (i == 0) {
            bits_to_remove = A_allzeros;
          } else {
            bits_to_remove = intersect(A_allzeros, bits_to_remove);
          }
        }
        
        if (bits_to_remove.n_elem > 0) {
          // Loop to reduce matrices and package them
          for (int i = 0; i < no_matrices; i++) {
            if (!is<NumericMatrix>(imported_U(i))) continue;
            arma::mat U_arma = as<arma::mat>(imported_U(i));
            arma::mat F_arma = as<arma::mat>(imported_F(i));
            
            U_arma.shed_rows(bits_to_remove);
            F_arma.shed_rows(bits_to_remove);
            
            U_arma.shed_cols(bits_to_remove);
            F_arma.shed_cols(bits_to_remove);
            
            arma::uvec U_nonzeros = find(U_arma);
            arma::uvec F_nonzeros = find(F_arma);
            u_transitions = u_transitions + static_cast<int>(U_nonzeros.n_elem);
            f_transitions = f_transitions + static_cast<int>(F_nonzeros.n_elem);
            
            imported_U(i) = U_arma;
            imported_F(i) = F_arma;
            
            if (!simplicity) {
              arma::mat A_arma = as<arma::mat>(imported_A(i));
              A_arma.shed_rows(bits_to_remove);
              A_arma.shed_cols(bits_to_remove);
              imported_A(i) = A_arma;
            }
          }
          
          if (!simplicity) mat_list(0) = imported_A;
          mat_list(1) = imported_U;
          mat_list(2) = imported_F;
          
          IntegerVector R_bits_to_remove = as<IntegerVector>(wrap(bits_to_remove));
          if (stage) {
            if (historical) {
              hstages = df_shedrows(hstages, R_bits_to_remove);
            } else if (age) {
              agestages = df_shedrows(agestages, R_bits_to_remove);
            } else {
              ahstages = df_shedrows(ahstages, R_bits_to_remove);
            }
          } else {
            ahstages = df_shedrows(ahstages, R_bits_to_remove);
          }
        } else {
          // No matrix reduction, with reduce == TRUE
          for (int i = 0; i < no_matrices; i++) {
            if (!is<NumericMatrix>(imported_U(i))) continue;
            arma::mat U_arma = as<arma::mat>(imported_U(i));
            arma::mat F_arma = as<arma::mat>(imported_F(i));
            
            arma::uvec U_nonzeros = find(U_arma);
            arma::uvec F_nonzeros = find(F_arma);
            
            u_transitions = u_transitions + static_cast<int>(U_nonzeros.n_elem);
            f_transitions = f_transitions + static_cast<int>(F_nonzeros.n_elem);
          }
        }
      } else {
        for (int i = 0; i < no_matrices; i++) {
          arma::sp_mat A_arma;
          
          if (simplicity) {
            if (!is<S4>(imported_U(i))) continue;
            
            arma::sp_mat U_arma = as<arma::sp_mat>(imported_U(i));
            arma::sp_mat F_arma = as<arma::sp_mat>(imported_F(i));
            
            A_arma = U_arma + F_arma;
          } else {
            if (!is<S4>(imported_U(i))) continue;
            A_arma = as<arma::sp_mat>(imported_A(i));
          }
          
          int A_cols = A_arma.n_cols;
          arma::uvec A_colrowzeros (A_cols, fill::zeros);
          
          for (int j = 0; j < A_cols; j++) {
            double col_sum = accu(A_arma.col(j));
            double row_sum = accu(A_arma.row(j));
            
            if ((col_sum + row_sum) == 0.0) A_colrowzeros(j) = 1;
          }
          
          arma::uvec A_allzeros = find(A_colrowzeros);
          
          if (i == 0) {
            bits_to_remove = A_allzeros;
          } else {
            bits_to_remove = intersect(A_allzeros, bits_to_remove);
          }
        }
        
        int remove_bits = static_cast<int>(bits_to_remove.n_elem);
        if (remove_bits > 0) {
          // Loop to reduce matrices and package them
          for (int i = 0; i < no_matrices; i++) {
            if (!is<S4>(imported_U(i))) continue;
            
            arma::sp_mat U_arma = as<arma::sp_mat>(imported_U(i));
            arma::sp_mat F_arma = as<arma::sp_mat>(imported_F(i));
            
            for (int j = 0; j < remove_bits; j++) {
              U_arma.shed_row(bits_to_remove(remove_bits - (j + 1)));
              F_arma.shed_row(bits_to_remove(remove_bits - (j + 1)));
              
              U_arma.shed_col(bits_to_remove(remove_bits - (j + 1)));
              F_arma.shed_col(bits_to_remove(remove_bits - (j + 1)));
            }
            arma::uvec U_nonzeros = find(U_arma);
            arma::uvec F_nonzeros = find(F_arma);
            u_transitions = u_transitions + static_cast<int>(U_nonzeros.n_elem);
            f_transitions = f_transitions + static_cast<int>(F_nonzeros.n_elem);
            
            imported_U(i) = U_arma;
            imported_F(i) = F_arma;
            
            if (!simplicity) {
              arma::sp_mat A_arma = as<arma::sp_mat>(imported_A(i));
              
              for (int j = 0; j < remove_bits; j++) {
                A_arma.shed_row(bits_to_remove(remove_bits - (j + 1)));
                A_arma.shed_col(bits_to_remove(remove_bits - (j + 1)));
              }
              imported_A(i) = A_arma;
            }
          }
          
          if (!simplicity) mat_list(0) = imported_A;
          mat_list(1) = imported_U;
          mat_list(2) = imported_F;
          
          IntegerVector R_bits_to_remove = as<IntegerVector>(wrap(bits_to_remove));
          if (stage) {
            if (historical) {
              hstages = df_shedrows(hstages, R_bits_to_remove);
            } else if (age) {
              agestages = df_shedrows(agestages, R_bits_to_remove);
            } else {
              ahstages = df_shedrows(ahstages, R_bits_to_remove);
            }
          } else {
            ahstages = df_shedrows(ahstages, R_bits_to_remove);
          }
        } else {
          // No matrix reduction, with reduce == TRUE
          for (int i = 0; i < no_matrices; i++) {
            if (!is<S4>(imported_U(i))) continue;
            arma::sp_mat U_arma = as<arma::sp_mat>(imported_U(i));
            arma::sp_mat F_arma = as<arma::sp_mat>(imported_F(i));
            
            arma::uvec U_nonzeros = find(U_arma);
            arma::uvec F_nonzeros = find(F_arma);
            
            u_transitions = u_transitions + static_cast<int>(U_nonzeros.n_elem);
            f_transitions = f_transitions + static_cast<int>(F_nonzeros.n_elem);
          }
        }
      }
      
    } else {
      // No matrix reduction
      if (!sparse_out) {
        for (int i = 0; i < no_matrices; i++) {
          if (!is<NumericMatrix>(imported_U(i))) continue;
          arma::mat U_arma = as<arma::mat>(imported_U(i));
          arma::mat F_arma = as<arma::mat>(imported_F(i));
          
          arma::uvec U_nonzeros = find(U_arma);
          arma::uvec F_nonzeros = find(F_arma);
          
          u_transitions = u_transitions + static_cast<int>(U_nonzeros.n_elem);
          f_transitions = f_transitions + static_cast<int>(F_nonzeros.n_elem);
        }
      } else {
        for (int i = 0; i < no_matrices; i++) {
          if (!is<S4>(imported_U(i))) continue;
          arma::sp_mat U_arma = as<arma::sp_mat>(imported_U(i));
          arma::sp_mat F_arma = as<arma::sp_mat>(imported_F(i));
          
          arma::uvec U_nonzeros = find(U_arma);
          arma::uvec F_nonzeros = find(F_arma);
          
          u_transitions = u_transitions + static_cast<int>(U_nonzeros.n_elem);
          f_transitions = f_transitions + static_cast<int>(F_nonzeros.n_elem);
        }
      }
    }
    
    qc(0) = u_transitions;
    qc(1) = f_transitions;
    qc(2) = no_matrices;
    
    return;
  }
  
  //' Assess if MPM is ahistorical, historical, age-by-stage, or Leslie
  //' 
  //' Function \code{whichbrew()} assesses whether a \code{lefkoMat} object is
  //' ahistorical, historical, age-based, or age-by-stage based.
  //' 
  //' @name whichbrew
  //' 
  //' @param ahstages The \code{ahstages} object from a \code{lefkoMat} object.
  //' @param hstages The \code{hstages} object from a \code{lefkoMat} object.
  //' @param agestages The \code{agestages} object from a \code{lefkoMat} object.
  //' 
  //' @return An integer corresponding to the type of \code{lefkoMat}object, with
  //' the following possible values: \code{0}: historical MPM, \code{1}:
  //' ahistorical MPM, \code{2}: age-by-stage MPM, and \code{3}: age-based MPM.
  //' 
  //' @keywords internal
  //' @noRd
  inline int whichbrew (Rcpp::DataFrame& ahstages, Rcpp::DataFrame& hstages,
    Rcpp::DataFrame& agestages) {
    
    int current_brew {1}; // 0 - hist; 1 - ahist; 2 - agestage; 3 - age
    
    int hst_cols = static_cast<int>(hstages.length());
    int ast_cols = static_cast<int>(agestages.length());
    
    if (hst_cols > 1) {
      current_brew = 0;
      
    } else if (ast_cols > 1) {
      current_brew = 2;
      
    } else {
      StringVector stage = as<StringVector>(ahstages["stage"]);
      if (stringcompare_simple(as<std::string>(stage(0)), "Age")) {
        current_brew = 3;
      } else current_brew = 1;
    }
    
    return current_brew;
  }
  
  //' Standardized Error Messages
  //' 
  //' Function \code{pop_error()} produces standardized error messages.
  //' 
  //' @name pop_error
  //' 
  //' @param input1 First string input.
  //' @param input2 Second string input.
  //' @param input3 Third string input.
  //' @param type Designates output message type.
  //' 
  //' @return Stops R and produces an error message.
  //' 
  //' @section Notes:
  //' Pop errors 3 and 6 were merged with 1.
  //' 
  //' @keywords internal
  //' @noRd
  inline void pop_error (String input1, String input2, String input3, int type = 1) {
    String eat_my_shorts;
    if (type == 1) { // Very useful
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " should be entered as ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 2) {
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be the same length as the number of ";
      eat_my_shorts += input2;
      eat_my_shorts += " entered in argument ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    } else if (type == 4) { // Very useful
      eat_my_shorts = "Matrix ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be square.";
      
    } else if (type == 5) {
      eat_my_shorts = "Matrix ";
      eat_my_shorts += input1;
      eat_my_shorts += " is not recognized.";
      
    } else if (type == 8) {
      eat_my_shorts = "Values input for ";
      eat_my_shorts += input1;
      eat_my_shorts += " in ";
      eat_my_shorts += input2;
      eat_my_shorts += " must be real numbers within set tolerance limits.";
      
    } else if (type == 9) {
      eat_my_shorts = "Variable names designating ";
      eat_my_shorts += input1;
      eat_my_shorts += " do not match variables names in input ";
      eat_my_shorts += input2;
      
    } else if (type == 10) { // Seems useful. Used many times
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be entered as a string vector showing ";
      eat_my_shorts += input2;
      eat_my_shorts += " in times t+1 and t, ";
      eat_my_shorts += "and time t-1 if a historical MPM is desired.";
      
    } else if (type == 11) {
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be entered if using ";
      eat_my_shorts += input2;
      eat_my_shorts += " object.";
      
    } else if (type == 13) { 
      eat_my_shorts = input1;
      eat_my_shorts += " is not recognized in arguments ";
      eat_my_shorts += input2;
      eat_my_shorts += " and ";
      eat_my_shorts += input3;
      
    } else if (type == 15) { // Very useful
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " is required if ";
      eat_my_shorts += input2;
      eat_my_shorts += " is not provided.";
      
    } else if (type == 16) { // Very useful
      eat_my_shorts = "Variable(s) coding for ";
      eat_my_shorts += input1;
      eat_my_shorts += " not found in dataset.";
      
    } else if (type == 17) { // Could this be removed? Used several times
      eat_my_shorts = "Some input ";
      eat_my_shorts += input1;
      eat_my_shorts += " values are not found in the ";
      eat_my_shorts += input2;
      eat_my_shorts += ".";
      
    } else if (type == 19) { // Could this be dropped for 29
      eat_my_shorts = "Vector ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be the same length as the number of ";
      eat_my_shorts += input2;
      eat_my_shorts += " in the ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    } else if (type == 20) {
      eat_my_shorts = "Matrix ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be the same number of columns as the number of ";
      eat_my_shorts += input2;
      eat_my_shorts += " in the ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    } else if (type == 21) { // Seems useful
      eat_my_shorts = "Some ";
      eat_my_shorts += input1;
      eat_my_shorts += " are not in an accepted style.";
      
    } else if (type == 22) {
      eat_my_shorts = "Argument ";
      eat_my_shorts += input1;
      eat_my_shorts += " must be a ";
      eat_my_shorts += input2;
      eat_my_shorts += " corresponding to the ";
      eat_my_shorts += input3;
      eat_my_shorts += ".";
      
    }
    
    throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
    
    return;
  }
  
  //' Format All Density-related Variables Based on Density Inputs
  //' 
  //' Function \code{density_prep()} takes the objects input in arguement
  //' \code{density} within the core projection functions and formats a number
  //' of control variables dictating how density will be interpreted during the
  //' projection.
  //' 
  //' @name density_prep
  //' 
  //' @param dens_index A reference to an empty List, which will be modified via
  //' this function. Will provide an index to the elements of the matrix to be
  //' altered during density adjustment.
  //' @param dyn_style An empty Armadillo unsigned integer vector, given as a
  //' reference to be modified. Will provide the density dependence styles to be
  //' used in density adjustment.
  //' @param dyn_alpha An empty Armadillo double floating point vector, given as
  //' a reference to be modified. Will provide the values of alpha to be used in
  //' density adjustment.
  //' @param dyn_beta An empty Armadillo double floating point vector, given as
  //' a reference to be modified. Will provide the values of beta to be used in
  //' density adjustment.
  //' @param dyn_gamma An empty Armadillo double floating point vector, given as
  //' a reference to be modified. Will provide the values of gamma to be used in
  //' density adjustment.
  //' @param dens_input The data frame input via the \code{density} argument, or
  //' a single data frame from the list input there.
  //' @param hstages The \code{hstages} data frame within the input lefkoMat
  //' object.
  //' @param stageframe The stageframe, typically within object \code{ahstages}
  //' within the lefkoMat object input.
  //' @param exp_tol The exponent maximmum limit being used.
  //' @param historical A Boolean value indicating whether the lefkoMat object
  //' is historical.
  //' 
  //' @return No objects are returned, though some inputs are modified.
  //' 
  //' @keywords internal
  //' @noRd
  inline void density_prep (List& dens_index, arma::uvec& dyn_style,
    arma::vec& dyn_alpha, arma::vec& dyn_beta, arma::vec& dyn_gamma,
    const DataFrame dens_input, const DataFrame hstages,
    const DataFrame stageframe, const double exp_tol, const bool historical) {
    
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
        
        if (static_cast<int>(find_hst3.n_elem) == 0 || static_cast<int>(find_hst2r.n_elem) == 0) {
          throw Rcpp::exception("Some stages entered in argument density do not exist in the life history model.",
            false);
        }
        if (static_cast<int>(find_hst1.n_elem) == 0 || static_cast<int>(find_hst2c.n_elem) == 0) {
          throw Rcpp::exception("Some stages entered in argument density do not exist in the life history model.",
            false);
        }
        
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
        
        if (static_cast<int>(find_ahst3.n_elem) == 0 || static_cast<int>(find_ahst2.n_elem) == 0) {
          throw Rcpp::exception("Some stages entered in argument density do not exist in the life history model.",
            false);
        }
        
        di_stage32_id(i) = find_ahst3(0);
        di_stage21_id(i) = find_ahst2(0);
        di_index(i) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
        
        ahst_3.zeros();
        ahst_2.zeros();
      }
      
      dens_index = Rcpp::List::create(_["index3"] = di_stage32_id,
        _["index2"] = di_stage21_id, _["index321"] = di_index);
    }
    
    dyn_style = as<arma::uvec>(dens_input["style"]);
    dyn_alpha = as<arma::vec>(dens_input["alpha"]);
    dyn_beta = as<arma::vec>(dens_input["beta"]);
    dyn_gamma = as<arma::vec>(dens_input["gamma"]);
    
    for (int i = 0; i < static_cast<int>(dyn_style.n_elem); i++) {
      if (dyn_style(i) < 1 || dyn_style(i) > 6) pop_error("density inputs", "", "", 21);
      
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
  
  //' Format All Equivalence Weight-related Variables Based on Input
  //' 
  //' @name equivalence_prep
  //' 
  //' @param equivalence_frame An empty data frame to be modified by this
  //' function. Will hold the core stage weights data.
  //' @param equivalence_vec An empty NumericVector to hold the stage-by-stage
  //' equivalence weights to be used.
  //' @param equivalence_vec_arma An empty Armadillo double floating point
  //' vector to hold the stage-by-stage equivalence weights to be used.
  //' @param eq_vec_length An empty integer to hold the length of the
  //' equivalence vectors.
  //' @param stage_weights The RObject version of the object input in argument
  //' \code{stage_weights}, or the first element of the list input there.
  //' @param format The MPM format code being used.
  //' @param used_matsize The matrix dimensions, as previously measured.
  //' 
  //' @return No objects are returned, though some inputs are modified.
  //' 
  //' @keywords internal
  //' @noRd
  inline void equivalence_prep (DataFrame& equivalence_frame,
    NumericVector& equivalence_vec, arma::vec& equivalence_vec_arma,
    int& eq_vec_length, const RObject stage_weights, const int format,
    const int used_matsize) {
    
    if (is<DataFrame>(stage_weights)) {
      equivalence_frame = as<DataFrame>(stage_weights);
      if (!equivalence_frame.hasAttribute("class")) {
        throw Rcpp::exception("Argument stage_weights should be a data frame of class lefkoEq.", 
          false);
      }
      
      CharacterVector eq_list_df_class = equivalence_frame.attr("class");
      bool found_lefkoEq {false};
      for (int j = 0; j < static_cast<int>(eq_list_df_class.length()); j++) {
        if (eq_list_df_class(j) == "lefkoEq" || eq_list_df_class(j) == "adaptEq") found_lefkoEq = true;
      }
      
      if (!found_lefkoEq) {
        throw Rcpp::exception("Argument stage_weights should be a data frame of class lefkoEq.",
          false);
      }
      
      IntegerVector eq_s2 = as<IntegerVector>(equivalence_frame["stage_id_2"]);
      IntegerVector eq_s1 = as<IntegerVector>(equivalence_frame["stage_id_1"]);
      IntegerVector eq_a2 = as<IntegerVector>(equivalence_frame["age2"]);
      IntegerVector eq_rn = clone(as<IntegerVector>(equivalence_frame["row_num"]));
      NumericVector eq_val = as<NumericVector>(equivalence_frame["value"]);
      
      eq_rn = eq_rn - 1;
      
      if (format < 3) {
        if (IntegerVector::is_na(eq_s1(0))) {
          throw Rcpp::exception("Enter stage pairs in lefkoEq objects used for historical MPMs.", 
            false);
        }
        if (IntegerVector::is_na(eq_s2(0))) {
          throw Rcpp::exception("Entries in column stage2 of lefkoEq objects cannot be empty except in Leslie MPMs.",
            false);
        }
      } else if (format > 3) {
        if (IntegerVector::is_na(eq_a2(0))) {
          throw Rcpp::exception("Enter ages in lefkoEq objects used for age-by-stage MPMs.",
            false);
        }
        if (format == 4) {
          if (IntegerVector::is_na(eq_s2(0))) {
            throw Rcpp::exception("Entries in column stage2 of lefkoEq objects cannot be empty except in Leslie MPMs.",
              false);
          }
        }
      } else {
        if (IntegerVector::is_na(eq_s2(0))) {
          throw Rcpp::exception("Entries in column stage2 of lefkoEq objects cannot be empty except in Leslie MPMs.",
            false);
        }
      }
      
      if (max(eq_rn) > used_matsize) {
        throw Rcpp::exception("Some row numbers in an entered lefkoEq object are too high.", 
          false);
      }
      
      if (min(eq_val) < 0.0) {
        throw Rcpp::exception("Values in argument stage weights may not be negative.", false);
      }
      
      NumericVector current_eq (used_matsize, 1.0);
      for (int j = 0; j < static_cast<int>(eq_rn.length()); j++) {
        current_eq(eq_rn(j)) = eq_val(j);
      }
      
      equivalence_vec = current_eq;
      equivalence_vec_arma = as<arma::vec>(equivalence_vec);
      eq_vec_length = static_cast<int>(equivalence_vec.length());
      
    } else if (is<NumericVector>(stage_weights) || is<IntegerVector>(stage_weights)) {
      equivalence_vec = as<NumericVector>(stage_weights);
      eq_vec_length = static_cast<int>(equivalence_vec.length());
      
      if (eq_vec_length != used_matsize) {
        throw Rcpp::exception("Numeric vector in argument stage_weights should be as long as there are rows in each matrix.", 
          false);
      }
      
      equivalence_vec_arma = as<arma::vec>(equivalence_vec);
      
    } else {
      throw Rcpp::exception("Argument stage_weights should be a data frame of class lefkoEq.", 
        false);
    }
  }
  
}
#endif
