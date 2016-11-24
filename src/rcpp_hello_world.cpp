
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world()
{
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y = NumericVector::create(0.0, 1.0);
  List z = List::create(x, y);
  return z;
}

//' @title Sum two numbers
//' @description Sum two numbers
//' @param a A number
//' @param b Another number
//'
// [[Rcpp::export]]
NumericVector test_add(NumericVector a, NumericVector b)
{
  NumericVector out;
  out = a + b;
  return out;
}
