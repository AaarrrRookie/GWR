#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace Rcpp;
using namespace RcppEigen;

// [[Rcpp::export]]
List MoranI(Eigen::Map<Eigen::VectorXi> Match,Eigen::Map<Eigen::VectorXi>  plz, Eigen::Map<Eigen::VectorXi>  Dplz, Eigen::Map<Eigen::MatrixXd> D, Eigen::Map<Eigen::VectorXd> Y,double TSS) {
  int n = plz.size();
  int nr = D.rows(), nplz=0;

  double zww4, s1=0, WSS=0;

  int a;
  NumericVector wY(n),SelbstGW(n), w(n), d(nr);

  for (int j=0; j<nr; ++j){

    Rcpp::checkUserInterrupt();
    a=Dplz[j];


    for (int k = 0; k < n; ++k) {
      w[k]=D(Match[k],j);
    }
    w=w/sum(w);

    nplz=0;
    double wwY=0;
    for(int i=0; i<n; ++i){
	  if(w[i]>0){
		  wwY+=w[i]*Y[i];
		  zww4+=4*w[i]*w[i];
	  }
      if(a==plz[i]){
        nplz+=1;
      }
    }

    s1+=zww4*nplz;
    zww4=0;


    for (int k = 0; k < n; ++k) {

      if(a==plz[k]){
        WSS+=(wwY-(w[k]*Y[k]))*Y[k];
      }

    }

  }

  double I=WSS/TSS;


  return Rcpp::List::create(Rcpp::Named("I") = I,
                            Rcpp::Named("s1") = s1,
							             Rcpp::Named("n") = n);
}
