#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector GWR_core2(Eigen::Map<Eigen::VectorXi> Match, Eigen::Map<Eigen::MatrixXd> D, Eigen::Map<Eigen::VectorXi> i, Eigen::Map<Eigen::VectorXi> p,int n, int nc) {

  int nr = D.rows();
  NumericMatrix b(nr,nc);
  NumericVector wwX(nc);

  for (int j=0; j<nr; ++j){

    std::fill( wwX.begin(),  wwX.end(), 0);
    double summe_w = 0;
    for (int k = 0; k < n; ++k) {
      Rcpp::checkUserInterrupt();
      if(D(Match[k],j) > 0){
        summe_w+= D(Match[k],j);
        for(int d = p[k]; d < p[k+1]; ++d){
          wwX[i[d]]+= D(Match[k],j);
        }
      }
    }

    for(int d = 0; d < nc; ++d){
      b(j,d)=wwX[d]/summe_w;
    }

    if(j%500 == 0){
      Rcpp::Rcout << "Zeilen abgeschlossen: " << 100*j/nr << "%" << std::endl;
    }

  }

  return b;

}
