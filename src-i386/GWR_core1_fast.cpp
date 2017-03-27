// Kern der GWR-Regression mit einer zustälzlichen additinven Variable
// Benötigt Bibliotheken Rcpp und RcppEigen

// Nur Auswertung der Koeffizienten ohne Varianz

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector GWR_core_time(Eigen::Map<Eigen::VectorXi> Match, Eigen::Map<Eigen::MatrixXd> D, Eigen::Map<Eigen::VectorXd> Y) {
  int n = Y.size();
  int nr = D.rows();
  NumericVector b(nr);

  for (int j=0; j<nr; ++j){

  Rcpp::checkUserInterrupt();
	double wwY = 0;
	double summe_w = 0;
    for (int k = 0; k < n; ++k) {
	  if(D(Match[k],j) > 0){
		  summe_w+= D(Match[k],j);
		  wwY+= D(Match[k],j)*Y[k];
	  }
    }
    wwY=wwY/summe_w;

	b[j]=wwY;

	if(j%500 == 0){
		Rcpp::Rcout << "Zeilen abgeschlossen: " << 100*j/nr << "%" << std::endl;
    }
  }

  return b;

}
