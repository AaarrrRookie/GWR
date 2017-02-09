#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cross_vec(Eigen::Map<Eigen::VectorXi> Match, Eigen::Map<Eigen::MatrixXd> Sv,int n, int nc) {

  NumericMatrix b(nc,nc);
  std::fill( b.begin(),  b.end(), 0);


  for(int iA=0; iA<nc; iA++){
    for(int iB=0; iB<nc; iB++){
      for(int i=0; i<n; i++){
        b(iA,iB)+=Sv(Match[i],iB)*Sv(Match[i],iA);
      }
    }

    if(iA%10 == 0){
      Rcpp::Rcout << "Zeilen abgeschlossen: " << 100*iA/nc << "%" << std::endl;
    }

  }


  return b;

}
