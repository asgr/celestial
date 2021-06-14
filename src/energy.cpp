#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector potential_part(NumericMatrix part, Nullable<NumericMatrix> eval = R_NilValue,
                       NumericVector mass = 1, double soft = 1, double Munit = 1,
                       double Lunit = 1000, double Vunit = 1){
  
  double G = 6.67384e-11;
  double msol_to_kg = 1.98892e+30;
  double pc_to_m = 3.08568e+16;
  double g = G * msol_to_kg/(pc_to_m);
  g = g * Munit/(Lunit * Vunit * Vunit);
  // Final units are (Lunit*[pc]) * (Vunit*[m/s])^2/(Munit*[Msol])^2, so default is [kpc].[km/s]^2/[Msol]
  
  int npart = part.nrow();
  
  if(mass.length() == 1){
    mass = rep(mass, npart);
  }
  
  double soft2 = soft*soft;
  
  double newpot = 0;
  
  if(eval.isNull()){
    NumericVector pot_part(npart);
    for (int i = 0; i < npart; i++) {
      for (int j = 0; j < i; j++) {
        newpot = g / sqrt(
          pow(part(i,0) - part(j,0),2) +
          pow(part(i,1) - part(j,1),2) +
          pow(part(i,2) - part(j,2),2) +
          soft2
        );
        pot_part(i) -= newpot * mass(j);
        pot_part(j) -= newpot * mass(i);
      }
    }
    return pot_part;
  }else{
    NumericMatrix evalnonnull(eval); 
    int neval = evalnonnull.nrow();
    NumericVector pot_eval(neval);
    for (int i = 0; i < neval; i++) {
      for (int j = 0; j < npart; j++) {
        newpot =  g / sqrt(
          pow(evalnonnull(i,0) - part(j,0),2) +
          pow(evalnonnull(i,1) - part(j,1),2) +
          pow(evalnonnull(i,2) - part(j,2),2) +
          soft2
        );
        pot_eval(i) -= newpot * mass(j);
      }
    }
    return pot_eval;
  }
}

