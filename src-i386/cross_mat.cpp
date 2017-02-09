#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cross_mat(Eigen::Map<Eigen::VectorXi> Match, Eigen::Map<Eigen::VectorXi> ii, Eigen::Map<Eigen::VectorXi> p, Eigen::Map<Eigen::MatrixXd> Sv,int n, int nc) {

  NumericMatrix b(nc,nc);
  std::fill( b.begin(),  b.end(), 0);


  for(int iA=0; iA<nc; iA++){
    for(int iB=0; iB<nc; iB++){
      for(int i=p[iA]; i<p[iA+1]; i++){
        b(iA,iB)+=Sv(Match[ii[i]],iB);
      }
    }

    if(iA%10 == 0){
      Rcpp::Rcout << "Zeilen abgeschlossen: " << 100*iA/nc << "%" << std::endl;
    }

  }


  return b;

}
