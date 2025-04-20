#ifndef LEFKOUTILS_input_stuff_H
#define LEFKOUTILS_input_stuff_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// Function index:
// 1. bool stringcompare_simple  Compares Two Strings, Assessing Inclusion
// 2. void RObj_TF_input_check  Take Generic RObject and Determine Value of 2 Boolean Variables With It
// 3. void RObj_DFr_input_check  Take Generic RObject and Extract a Data Frame, with Boolean Response



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
  
}
#endif