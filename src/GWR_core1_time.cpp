// Kern der GWR-Regression mit einer zustälzlichen additinven Variable
// Benötigt Bibliotheken Rcpp und RcppEigen

// Nur Auswertung der Koeffizienten ohne Varianz
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppEigen.h>
#include <progress.hpp>



using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector GWR_core_time(Eigen::Map<Eigen::VectorXi> Match,
                            Eigen::Map<Eigen::VectorXi> Match_T,
                            Eigen::Map<Eigen::MatrixXd> D,
                            Eigen::Map<Eigen::MatrixXd> D_T,
                            Eigen::Map<Eigen::VectorXd> Y) {


  int n = Y.size();
  int nr = D.rows();
  int nt = D_T.rows();
  NumericVector b(nr*nt);
  double temp;
  Progress p(nr*nt, true);

  for(int time=0; time<nt; ++time){
          if (Progress::check_abort() ){
            return -1.0;
          }
      for (int j=0; j<nr; ++j){
          p.increment();

	      double wwY = 0;
	      double summe_w = 0;

        for (int k = 0; k < n; ++k) {
	          if(D(Match[k],j) > 0){
	              temp = D(Match[k],j)*D_T(Match_T[k],time);
		            summe_w+= temp;
		            wwY+= temp*Y[k];
	           }
        }

      wwY=wwY/summe_w;

	    b[j+time*nr]=wwY;

  }
}

  return b;

}




NumericVector GWR_poisson_time(Eigen::Map<Eigen::VectorXi> Match,
                            Eigen::Map<Eigen::VectorXi> Match_T,
                            Eigen::Map<Eigen::MatrixXd> D,
                            Eigen::Map<Eigen::MatrixXd> D_T,
                            Eigen::Map<Eigen::VectorXd> Y) {


  int n = Y.size();
  int nr = D.rows();
  int nt = D_T.rows();
  NumericVector lambda(nr*nt);
  double temp;
  Progress p(nr*nt, true);

  for(int time=0; time<nt; ++time){
    if (Progress::check_abort() ){
      return -1.0;
    }
    for (int j=0; j<nr; ++j){
      p.increment();

      double wwY = 0;
      double summe_w = 0;

      for (int k = 0; k < n; ++k) {
        if(D(Match[k],j) > 0){
          temp = D(Match[k],j)*D_T(Match_T[k],time);
          summe_w+= temp;
          wwY+= temp*Y[k];
        }
      }

      wwY=wwY/summe_w;

      lambda[j+time*nr]=wwY;

    }
  }

  return lambda;

}
