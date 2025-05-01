#ifndef LEFKOUTILS_input_stuff_H
#define LEFKOUTILS_input_stuff_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// Function index:
// 1. bool stringcompare_input  Compares Two Strings, Assessing Inclusion
// 2. void RObj_TF_input_check  Take Generic RObject and Determine Value of 2 Boolean Variables With It
// 3. void RObj_DFr_input_check  Take Generic RObject and Extract a Data Frame, with Boolean Response
// 4. bool yesno_to_logic  Take Yes / No and Other Input to Yield a Boolean Value
// 5. void yesnoauto_to_logic  Take Yes / No and Other Input to Yield a Boolean Value
// 6. void integer_vectorizer  Create Standardized IntegerVectors Based on Non-Standard Input
// 7. void numeric_vectorizer  Create Standardized NumericVectors Based on Non-Standard Input
// 8. void integer_char_vectorizer  Create Standardized IntegerVectors Based on Non-Standard Integer and String Input



namespace LefkoInputs {
  
  //' Compares Two Strings, Assessing Inclusion
  //' 
  //' This function compares two strings, and will assess whether \code{str2} is
  //' contained within \code{str1}. It is a simpler version of 
  //' \code{stringcompare_soft()} that yields only the logical result. This
  //' function is the same as \code{stringcompare_simple()}.
  //' 
  //' @name stringcompare_input
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
  inline bool stringcompare_input(std::string str1, std::string str2,
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
  
  //' Take Generic RObject and Determine Value of 2 Boolean Variables With It
  //' 
  //' Function \code{RObj_TF_input_check()} determines the values of up to 2
  //' Boolean variables on the basis of a single Nullable RObject input, which
  //' can be entered in Logical, Integer, Numeric, or String formats.
  //' 
  //' @name RObj_TF_input_check
  //' 
  //' @param argument_name A String giving the name of the input argument.
  //' @param thirdvalue A String value to compare, for use in determining the
  //' second Boolean variable.
  //' @param output1 A Boolean variable based on yes/no or true/false inputs in
  //' argument \code{input}.
  //' @param output2 A Boolean variable determining if the value in \code{input}
  //' matches the value given in \code{thirdvalue}.
  //' @param thirdvalue_used A Boolean variable determining whether to search
  //' for \code{thirdvalue} in \code{input}.
  //' @param firstvalue_if_thirdtrue The Boolean value of \code{output1} if
  //' \code{output2} is true.
  //' @param input The RObject input by the user.
  //' 
  //' @return This function modifies two Boolean variables provided by
  //' reference. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void RObj_TF_input_check (String argument_name,
    String thirdvalue, bool& output1, bool& output2, bool thirdvalue_used = false,
    bool firstvalue_if_thirdtrue = true, Nullable<RObject> input = R_NilValue) {
    
    if (input.isNotNull()) {
      if (is<CharacterVector>(input)) {
        StringVector yesbits = {"y", "yes", "yea", "yeah", "t", "true", "ja", "tak"};
        StringVector nobits = {"n", "no", "non", "nah", "f", "false", "nein", "nie"};
        
        int yesbits_length = static_cast<int>(yesbits.length());
        int nobits_length = static_cast<int>(nobits.length());
        
        StringVector thirdbits;
        int thirdvalue_length {0};
        
        if (thirdvalue_used) {
          std::string thirdvalue_c = thirdvalue.get_cstring();
          thirdvalue_length = thirdvalue_c.length();
          bool tooclose {false};
          
          if (thirdvalue_c.substr(0) == "y" || thirdvalue_c.substr(0) == "Y" ||
              thirdvalue_c.substr(0) == "t" || thirdvalue_c.substr(0) == "T") {
            thirdvalue_length--;
            tooclose = true;
          }
          if (thirdvalue_c.substr(0) == "n" || thirdvalue_c.substr(0) == "N" ||
              thirdvalue_c.substr(0) == "f" || thirdvalue_c.substr(0) == "F") {
            thirdvalue_length--;
            tooclose = true;
          }
          
          if (thirdvalue_length < 1) throw Rcpp::exception("Argument thirdvalue not useable.", false);
          
          int start_index {0};
          if (tooclose) start_index++;
          
          StringVector thirdbits_pre (thirdvalue_length);
          for (int i = 0; i < thirdvalue_length; i++) {
            thirdbits_pre(i) = String(thirdvalue_c.substr(start_index, (start_index + i)));
          }
          thirdbits = thirdbits_pre;
        }
        
        StringVector input_check_vec = as<StringVector>(input);
        String input_check = String(input_check_vec(0));
        
        int third_check {0};
        int yes_check {0};
        int no_check {0};
        
        IntegerVector string_limits = {yesbits_length, nobits_length, thirdvalue_length};
        
        for (int i = 0; i < max(string_limits); i++) {
          if (i < yesbits_length) {
            if (stringcompare_input(input_check, String(yesbits(i)), true)) yes_check++;
          }
          if (i < nobits_length) {
            if (stringcompare_input(input_check, String(nobits(i)), true)) no_check++;
          }
          if (i < thirdvalue_length) {
            if (stringcompare_input(input_check, String(thirdbits(i)), true)) third_check++;
          }
        }
        
        if (third_check > 0) { 
          output1 = firstvalue_if_thirdtrue;
          output2 = true;
        } else if (yes_check > 0) {
          output1 = true;
          output2 = false;
        } else if (no_check > 0) {
          output1 = false;
          output2 = false;
        } else {
          String err_out = "Argument ";
          err_out += argument_name;
          err_out += " is invalid.";
          
          throw Rcpp::exception(err_out.get_cstring(), false);
        }
      } else if (is<LogicalVector>(input)) {
        LogicalVector input_check_vec = as<LogicalVector>(input);
        
        if (!LogicalVector::is_na(input_check_vec(0))) {
          output1 = static_cast<bool>(input_check_vec(0));
          output2 = false;
        } else {
          String err_out = "Argument ";
          err_out += argument_name;
          err_out += " may not equal NA.";
          
          throw Rcpp::exception(err_out.get_cstring(), false);
        }
      } else if (is<IntegerVector>(input)) {
        IntegerVector input_check_vec = as<IntegerVector>(input);
        
        if (!IntegerVector::is_na(input_check_vec(0))) {
          int input_first = static_cast<int>(input_check_vec(0));
          if (input_first > 0) {
            output1 = true;
            output2 = false;
          } else if (input_first == 0) {
            output1 = false;
            output2 = false;
          } else {
            String err_out = "Argument ";
            err_out += argument_name;
            err_out += " may not equal NA.";
            
            throw Rcpp::exception(err_out.get_cstring(), false);
          }
        }
      } else if (is<NumericVector>(input)) {
        NumericVector input_check_vec = as<NumericVector>(input);
        
        if (!NumericVector::is_na(input_check_vec(0))) {
          int input_first = static_cast<int>(input_check_vec(0));
          if (input_first > 0) {
            output1 = true;
            output2 = false;
          } else if (input_first == 0) {
            output1 = false;
            output2 = false;
          } else {
            String err_out = "Argument ";
            err_out += argument_name;
            err_out += " may not equal NA.";
            
            throw Rcpp::exception(err_out.get_cstring(), false);
          }
        }
      } else {
          String err_out = "Argument ";
          err_out += argument_name;
          err_out += " is invalid.";
          
          throw Rcpp::exception(err_out.get_cstring(), false);
      }
    }
  }
  
  //' Take Generic RObject and Extract a Data Frame, with Boolean Response
  //' 
  //' Function \code{RObj_DFr_input_check()} extracts a data frame from a
  //' single Nullable RObject input. It also sets a Boolean variable related to
  //' this data frame.
  //' 
  //' @name RObj_DFr_input_check
  //' 
  //' @param argument_name A String giving the name of the input argument.
  //' @param classname A String value to compare the name of the class of the
  //' input to.
  //' @param output1 The data frame to extract to, by reference.
  //' @param output2 A Boolean variable related to the value of the input.
  //' @param classname_check A Boolean variable stating whether to check the
  //' class of the object for the value of \code{classname}.
  //' @param output2value_if_notnull A Boolean variable determining what to
  //' assign to \code{output2} if \code{input} is null.
  //' @param input The RObject input by the user.
  //' 
  //' @return This function modifies a data frame and one Boolean variable by
  //' reference. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void RObj_DFr_input_check (String argument_name, String classname,
    DataFrame& output1, bool& output2, bool classname_check = false,
    bool output2value_if_notnull = true, Nullable<RObject> input = R_NilValue) {
    
    if (input.isNotNull()) {
      String added_bit = "Argument ";
      added_bit += argument_name;
      added_bit += " must be an object of class ";
      added_bit += classname;
      
      output2 = output2value_if_notnull;
      
      if (is<DataFrame>(input)) {
        output1 = as<DataFrame>(input);
        
        if (output1.hasAttribute("class") && classname_check) {
          CharacterVector chosen_input_class = output1.attr("class");
          bool found_class {false};
          
          for (int j = 0; j < static_cast<int>(chosen_input_class.length()); j++) {
            if (chosen_input_class(j) == classname) found_class = true;
          }
          if (!found_class) {
            throw Rcpp::exception(added_bit.get_cstring(), false);
          }
        }
      } else if (classname_check) {
        throw Rcpp::exception(added_bit.get_cstring(), false);
      } else {
        String new_bit = "Argument ";
        new_bit += argument_name;
        new_bit += " must be a data frame.";
        
        throw Rcpp::exception(new_bit.get_cstring(), false);
      }
    } else {
      if (output2value_if_notnull) {
        output2 = false;
      } else {
        output2 = true;
      }
    }
  }
  
  //' Take Yes / No and Other Input to Yield a Boolean Value
  //' 
  //' Function \code{yesno_to_logic()} takes a variety of inputs and interprets
  //' them, creating a Boolean response.
  //' 
  //' @name yesno_to_logic
  //' 
  //' @param input RObject to be interpreted.
  //' @param arg_name Name of argument that is being tested.
  //' 
  //' @return Returns a simple Boolean value, or produces an error for
  //' unintelligible input.
  //' 
  //' @keywords internal
  //' @noRd
  inline bool yesno_to_logic (RObject input, String arg_name) {
    bool final_result = false;
    
    if (is<StringVector>(input)) {
      StringVector yesbits = {"y", "yes", "yea", "yeah", "t", "true", "ja", "tak"};
      StringVector nobits = {"n", "no", "non", "nah", "f", "false", "nein", "nie"};
      
      StringVector input_check_vec = as<StringVector>(input);
      String input_check = String(input_check_vec(0));
      
      int yes_check {0};
      int no_check {0};
      
      for (int i = 0; i < 8; i++) {
        if (stringcompare_input(input_check, String(yesbits(i)), true)) yes_check++;
        if (stringcompare_input(input_check, String(nobits(i)), true)) no_check++;
      }
      
      if (yes_check > 0) {
        final_result = true;
      } else if (no_check > 0) {
        final_result = false;
      } else {
        String err_out = "Argument ";
        err_out += arg_name;
        err_out += " is invalid.";
        
        throw Rcpp::exception(err_out.get_cstring(), false);
      }
    } else if (is<LogicalVector>(input)) {
        LogicalVector input_check_vec = as<LogicalVector>(input);
        final_result = static_cast<bool>(input_check_vec(0));
    } else if (is<NumericVector>(input)) {
        IntegerVector input_check_vec = as<IntegerVector>(input);
        int input_first = static_cast<int>(input_check_vec(0));
        
        if (input_first == 1) final_result = true;
    } else {
      String err_out = "Argument ";
      err_out += arg_name;
      err_out += " is invalid.";
      
      throw Rcpp::exception(err_out.get_cstring(), false);
    }
    
    return final_result;
  }
  
  //' Take Yes / No and Other Input to Yield a Boolean Value
  //' 
  //' Function \code{yesnoauto_to_logic()} takes a variety of inputs and
  //' interprets them, altering two Boolean responses that are input as
  //' arguments.
  //' 
  //' @name yesnoauto_to_logic
  //' 
  //' @param input RObject to be interpreted.
  //' @param arg_name Name of the function argument being tested.
  //' @param yesno_only A Boolean variable holding the yes / no value.
  //' @param auto_only A Boolean variable holding whether "auto" was chosen.
  //' 
  //' @return Alters two Boolean variables on memory based on the input.
  //' 
  //' @keywords internal
  //' @noRd
  inline void yesnoauto_to_logic (RObject input, String arg_name,
    bool &yesno_only, bool &auto_only) {
    
    bool yesno_result = false;
    bool auto_result = false;
    
    if (is<StringVector>(input)) {
      StringVector yesbits = {"y", "yes", "yea", "yeah", "t", "true", "ja", "tak"};
      StringVector nobits = {"n", "no", "non", "nah", "f", "false", "nein", "nie"};
      StringVector autobits = {"au", "aut", "auto", "both", "jidou"};
      
      StringVector input_check_vec = as<StringVector>(input);
      String input_check = String(input_check_vec(0));
      
      int auto_check {0};
      int yes_check {0};
      int no_check {0};
      
      for (int i = 0; i < 8; i++) {
        if (i < 5) {
          if (stringcompare_input(input_check, String(autobits(i)), true)) auto_check++;
        }
        if (stringcompare_input(input_check, String(yesbits(i)), true)) yes_check++;
        if (stringcompare_input(input_check, String(nobits(i)), true)) no_check++;
      }
      
      if (auto_check > 0) { 
        auto_result = true;
      } else if (yes_check > 0) {
        yesno_result = true;
      } else if (no_check > 0) {
        yesno_result = false;
      } else {
        String err_out = "Argument ";
        err_out += arg_name;
        err_out += " is invalid.";
        
        throw Rcpp::exception(err_out.get_cstring(), false);
      }
    } else if (is<LogicalVector>(input)) {
        LogicalVector input_check_vec = as<LogicalVector>(input);
        yesno_result = static_cast<bool>(input_check_vec(0));
    } else if (is<NumericVector>(input)) {
        IntegerVector input_check_vec = as<IntegerVector>(input);
        int input_first = static_cast<int>(input_check_vec(0));
        
        if (input_first == 1) yesno_result = true;
    } else {
      String err_out = "Argument ";
      err_out += arg_name;
      err_out += " is invalid.";
      
      throw Rcpp::exception(err_out.get_cstring(), false);
    }
    
    yesno_only = yesno_result;
    auto_only = auto_result;
  }
  
  //' Create Standardized IntegerVectors Based on Non-Standard Input
  //' 
  //' @name integer_vectorizer
  //' 
  //' @param output The output reference, passed by reference.
  //' @param input The input vector, treated as an \code{RObject}.
  //' @param argument_name The name of the argument used as \code{input}, given as
  //' a String.
  //' @param stage_length An integer giving the length of the \code{stage2}
  //' vector.
  //' @param age_length An integer giving the length of the \code{age2} vector.
  //' @param min_limit The smallest integer to allow, if using limits.
  //' @param max_limit The largest integer to allow, if using limits.
  //' @param use_limits A Boolean variable indicating whether to limit allowable
  //' values.
  //' @param NAasOther A Boolean value indicating whether to treat \code{NA}
  //' values as the value specified in \code{change_value}.
  //' @param change_value The integer to set \code{NA} values to if
  //' \code{NAasOther = TRUE}.
  //' 
  //' @return This function modifies an input vector by reference, given as
  //' argument \code{output}. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void integer_vectorizer (IntegerVector& output, Nullable<RObject> input,
    String argument_name, int stage_length, int age_length, int min_limit,
    int max_limit, bool use_limits = false, bool NAasOther = false,
    int change_value = 0) {
    
    IntegerVector input_;
    int input_length {0};
    
    if (input.isNotNull()) {
      if (is<IntegerVector>(input)) {
        input_ = as<IntegerVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (IntegerVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (use_limits) {
          for (int i = 0; i < input_length; i++) { 
            if (!IntegerVector::is_na(input_(i))) {
              if (input_(i) < min_limit || input_(i) > max_limit) {
                String eat_my_shorts = "Argument ";
                eat_my_shorts += argument_name;
                eat_my_shorts += " must be an integer between ";
                eat_my_shorts += min_limit;
                eat_my_shorts += " and ";
                eat_my_shorts += max_limit;
                eat_my_shorts += ".";
                
                throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
              }
            }
          }
        }
        
        if (input_length != stage_length && input_length != age_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as vector stage2 or age2.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<NumericVector>(input)) {
        input_ = as<IntegerVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (NumericVector::is_na(input_(i)) || IntegerVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (use_limits) {
          for (int i = 0; i < input_length; i++) { 
            if (!IntegerVector::is_na(input_(i)) && !NumericVector::is_na(input_(i))) {
              if (input_(i) < min_limit || input_(i) > max_limit) {
                String eat_my_shorts = "Argument ";
                eat_my_shorts += argument_name;
                eat_my_shorts += " must be an integer between ";
                eat_my_shorts += min_limit;
                eat_my_shorts += " and ";
                eat_my_shorts += max_limit;
                eat_my_shorts += ".";
                
                throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
              }
            }
          }
        }
        
        if (input_length != stage_length && input_length != age_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as vector stage2 or age2.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<LogicalVector>(input)) {
          if (age_length != 0) {
            if (NAasOther) {
              IntegerVector input_temp (age_length, change_value);
              input_ = input_temp;
            } else {
              IntegerVector input_temp (age_length, NA_INTEGER);
              input_ = input_temp;
            }
          } else if (stage_length != 0) {
            if (NAasOther) {
              IntegerVector input_temp (stage_length, change_value);
              input_ = input_temp;
            } else {
              IntegerVector input_temp (stage_length, NA_INTEGER);
              input_ = input_temp;
            }
          }
      } else {
        String eat_my_shorts = "Please enter argument ";
        eat_my_shorts += argument_name;
        eat_my_shorts += " in integer format.";
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
    } else {
      if (NAasOther) {
        if (age_length != 0) {
          IntegerVector input_temp (age_length, change_value);
          input_ = input_temp;
        } else if (stage_length != 0) {
          IntegerVector input_temp (stage_length, change_value);
          input_ = input_temp;
        }
      } else {
        if (age_length != 0) {
          IntegerVector input_temp (age_length, NA_INTEGER);
          input_ = input_temp;
        } else if (stage_length != 0) {
          IntegerVector input_temp (stage_length, NA_INTEGER);
          input_ = input_temp;
        }
      }
    }
    output = input_;
  }
  
  //' Create Standardized NumericVectors Based on Non-Standard Input
  //' 
  //' @name numeric_vectorizer
  //' 
  //' @param output The output reference, passed by reference.
  //' @param input The input vector, treated as an \code{RObject}.
  //' @param argument_name The name of the argument used as \code{input}, given as
  //' a String.
  //' @param stage_length An integer giving the length of the \code{stage2}
  //' vector.
  //' @param age_length An integer giving the length of the \code{age2} vector.
  //' @param NAasOther A Boolean value indicating whether to treat \code{NA}
  //' values as value given in \code{change_value}.
  //' @param change_value The numeric value to change \code{NA}s to, if
  //' \code{NAasOther = TRUE}.
  //' 
  //' @return This function modifies an input vector by reference, given as
  //' argument \code{output}. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void numeric_vectorizer (NumericVector& output, Nullable<RObject> input,
    String argument_name, int stage_length, int age_length,
    bool NAasOther = false, double change_value = 0) {
    
    NumericVector input_;
    int input_length {0};
    
    if (input.isNotNull()) {
      if (is<NumericVector>(input)) {
        input_ = as<NumericVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (NumericVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        if (input_length != stage_length && input_length != age_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as vector stage2 or age2.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<LogicalVector>(input)) {
          if (age_length != 0) {
            if (NAasOther) {
              NumericVector input_temp (age_length, change_value);
              input_ = input_temp;
            } else {
              NumericVector input_temp (age_length, NA_REAL);
              input_ = input_temp;
            }
          } else if (stage_length != 0) {
            if (NAasOther) {
              NumericVector input_temp (stage_length, change_value);
              input_ = input_temp;
            } else {
              NumericVector input_temp (stage_length, NA_REAL);
              input_ = input_temp;
            }
          }
      } else {
        String eat_my_shorts = "Please enter argument ";
        eat_my_shorts += argument_name;
        eat_my_shorts += " in numeric format.";
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
    } else {
      if (NAasOther) {
        if (age_length != 0) {
          NumericVector input_temp (age_length, change_value);
          input_ = input_temp;
        } else if (stage_length != 0) {
          NumericVector input_temp (stage_length, change_value);
          input_ = input_temp;
        }
      } else {
        if (age_length != 0) {
          NumericVector input_temp (age_length, NA_REAL);
          input_ = input_temp;
        } else if (stage_length != 0) {
          NumericVector input_temp (stage_length, NA_REAL);
          input_ = input_temp;
        }
      }
    }
    output = input_;
  }
  
  //' Create Standardized IntegerVectors Based on Non-Standard Integer and String Input
  //' 
  //' @name integer_char_vectorizer
  //' 
  //' @param output The output reference, passed by reference.
  //' @param input The input vector, treated as an \code{RObject}.
  //' @param argument_name The name of the argument used as \code{input}, given as
  //' a String.
  //' @param stage_length An integer giving the length of the \code{stage2}
  //' vector.
  //' @param age_length An integer giving the length of the \code{age2} vector.
  //' @param int_limits An integer vector giving the possible integer values.
  //' @param char_limits A character vector giving the possible range of values.
  //' @param use_limits A Boolean variable indicating whether to limit allowable
  //' values.
  //' @param NAasOther A Boolean value indicating whether to treat \code{NA}
  //' values as the value specified in \code{change_value}.
  //' @param change_value The integer to set \code{NA} values to if
  //' \code{NAasOther = TRUE}.
  //' 
  //' @return This function modifies an input vector by reference, given as
  //' argument \code{output}. No real output is returned.
  //' 
  //' @keywords internal
  //' @noRd
  inline void integer_char_vectorizer (IntegerVector& output, Nullable<RObject> input,
    String argument_name, int stage_length, int age_length,
    IntegerVector int_limits, CharacterVector char_limits,
    bool NAasOther = false, int change_value = 0) {
    
    IntegerVector input_;
    int input_length {0};
    
    int allowable = static_cast<int>(int_limits.length());
    
    if (input.isNotNull()) {
      if (is<IntegerVector>(input) || is<NumericVector>(input)) {
        input_ = as<IntegerVector>(input);
        
        input_length = static_cast<int>(input_.length());
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (NumericVector::is_na(input_(i)) || IntegerVector::is_na(input_(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        for (int i = 0; i < input_length; i++) { 
          if (!IntegerVector::is_na(input_(i)) && !NumericVector::is_na(input_(i))) {
            bool found_value {false};
            
            for (int j = 0; j < allowable; j++) {
              if (input_(i) == int_limits(j)) found_value = true;
            }
            if (!found_value) {
              String eat_my_shorts = "Argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " must be an integer between ";
              eat_my_shorts += int_limits(0);
              eat_my_shorts += " and ";
              eat_my_shorts += int_limits(allowable - 1);
              eat_my_shorts += ".";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          }
        }
        
        if (input_length != stage_length && input_length != age_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as vector stage2 or age2.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else if (is<CharacterVector>(input)) {
        IntegerVector new_int_limits = clone(int_limits);
        CharacterVector char_int_limits = as<CharacterVector>(new_int_limits);
        
        CharacterVector CharInputVec = as<CharacterVector>(input);
        int input_length = static_cast<int>(CharInputVec.length());
        input_ = IntegerVector (input_length);
        
        if (NAasOther) {
          unsigned int na_count {0};
          for (int i = 0; i < input_length; i++) { 
            if (CharacterVector::is_na(CharInputVec(i))) {
              input_(i) = change_value;
              na_count++;
            }
            if (na_count == 1) {
              String eat_my_shorts = "NA values in argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " will be treated as ";
              eat_my_shorts += change_value;
              eat_my_shorts += " values.";
              
              Rf_warningcall(R_NilValue, "%s", eat_my_shorts.get_cstring());
            }
          }
        }
        
        for (int i = 0; i < input_length; i++) { 
          if (!CharacterVector::is_na(CharInputVec(i))) {
            bool found_value {false};
            
            for (int j = 0; j < allowable; j++) {
              if (CharInputVec(i) == char_limits(j)) {
                input_(i) = int_limits(j);
                found_value = true;
              } else if (CharInputVec(i) == char_int_limits(j)) {
                input_(i) = int_limits(j);
                found_value = true;
              }
            }
            
            if (!found_value) {
              String eat_my_shorts = "Argument ";
              eat_my_shorts += argument_name;
              eat_my_shorts += " must be an integer between ";
              eat_my_shorts += int_limits(0);
              eat_my_shorts += " and ";
              eat_my_shorts += int_limits(allowable - 1);
              eat_my_shorts += ".";
              
              throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
            }
          }
        }
        
        if (input_length != stage_length && input_length != age_length) {
          String eat_my_shorts = "Argument ";
          eat_my_shorts += argument_name;
          eat_my_shorts += " must be the same length as vector stage2 or age2.";
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
      } else {
        String eat_my_shorts = "Please enter argument ";
        eat_my_shorts += argument_name;
        eat_my_shorts += " in integer format.";
        
        throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
      }
    } else {
      if (NAasOther) {
        if (age_length != 0) {
          IntegerVector input_temp (age_length, change_value);
          input_ = input_temp;
        } else if (stage_length != 0) {
          IntegerVector input_temp (stage_length, change_value);
          input_ = input_temp;
        }
      } else {
        if (age_length != 0) {
          IntegerVector input_temp (age_length, NA_INTEGER);
          input_ = input_temp;
        } else if (stage_length != 0) {
          IntegerVector input_temp (stage_length, NA_INTEGER);
          input_ = input_temp;
        }
      }
    }
    output = input_;
  }
}
#endif