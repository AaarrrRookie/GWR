#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cross_vec_mat(Eigen::Map<Eigen::VectorXi> Match, Eigen::Map<Eigen::MatrixXd> Sv, Eigen::Map<Eigen::VectorXd> Y,int n, int nc) {

  NumericVector b(n);
  std::fill( b.begin(),  b.end(), 0);


  for(int iA=0; iA<nc; iA++){
    for(int i=0; i<n; i++){
      b(i)+=Sv(Match[i],iA)*Y[iA];
    }


    if(iA%10 == 0){
      Rcpp::Rcout << "Zeilen abgeschlossen: " << 100*iA/nc << "%" << std::endl;
    }

  }


  return b;

}
