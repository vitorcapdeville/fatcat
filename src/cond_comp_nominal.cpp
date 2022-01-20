// #include <Rcpp.h>
// using namespace Rcpp;
//
// // This is a simple example of exporting a C++ function to R. You can
// // source this function into an R session using the Rcpp::sourceCpp
// // function (or via the Source button on the editor toolbar). Learn
// // more about Rcpp at:
// //
// //   http://www.rcpp.org/
// //   http://adv-r.had.co.nz/Rcpp.html
// //   http://gallery.rcpp.org/
// //
//
// // Só funciona para p = 1
// // [[Rcpp::export]]
// NumericVector probNominal(NumericMatrix beta, double f){
//   int k = beta.ncol() + 1;
//   NumericVector prob(k);
//   NumericVector aux(k-1);
//   aux[0] = exp(beta(0,0)*f);
//   aux[1] = exp(beta(0,1)*f);
//   aux[2] = exp(beta(0,2)*f);
//   double aux2 = sum(aux);
//   for(int s = 0; s <= k - 2; s++){
//     prob[s] = exp(beta(0,s)*f)/(1 + aux2);
//   }
//   prob[k-1] = 1/(1 + aux2);
//   return prob;
// }
//
//
// // [[Rcpp::export]]
// double logcondcompbetajkNominal(NumericVector f,NumericMatrix beta, NumericMatrix y,int j, int k){
//   double aux = 0;
//   int n = f.size();
//   for(int i = 0; i <= n-1; i++){
//     int ind = y(j-1,i) - 1;
//     aux += std::log(probNominal(beta,f[i])[ind]);
//   }
//   double aux2 = 0;
//   int t=0;
//   NumericVector aux3 = beta(_,k-1);
//   while(t<j-1 && t < aux3.size()-1){
//     aux2 -= (1/2000000)*std::pow(beta(t,k-1),2);
//     t += 1;
//   }
//   if(t==j-1){
//     aux2 -= (1/2000000)*std::pow(beta(t,k-1),2);
//   }
//   return aux+aux2;
// }
//
// // [[Rcpp::export]]
// double logcondcompfiNominal(double f, NumericMatrix beta,NumericMatrix y,int i){
//   double aux = 0;
//   int j = beta.nrow();
//   int k = beta.ncol();
//   int p = 1;
//   for(int l = 0; l <= j-1; l++){
//     int ind = y(l,i-1) - 1;
//     NumericMatrix aux3(1,k);
//     aux3(0,_) = beta(l,_);
//     aux += log(probNominal(aux3,f)[ind]);
//   }
//   double aux2 = 0;
//   for(int t=0; t<=p-1;t++){
//     aux2 -= std::pow(f,2)/2;
//   }
//   return(aux + aux2);
// }
