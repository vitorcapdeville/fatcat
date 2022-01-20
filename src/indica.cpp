#include <Rcpp.h>
using namespace Rcpp;

int indica(bool cond){
  if(cond){
    return 1;
  } else {
    return 0;
  }
}
