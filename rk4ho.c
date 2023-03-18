#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NVARS 2

void rhs1(double *x, double *res) {
  res[0] = x[1]; //df/dx f(0) = x0
  res[1] = -x[0];
}

void rk4sys(double *res, double t0, double *x0, double h, int n){
  // f es la funcion evaluada en x0 en todo el sistema
  double f[NVARS], k[NVARS][4]; //k1_0,k2_0,k3_0,k4_0 - k1_1,k2_1,k3_1,k4_1
  double x[NVARS]; // Iniciales del rk sobre x
  //
  for(int i = 0; i<NVARS; i++){
  x[i] = x0[i];     //Chinga tu madre Mau :)
  res[i] = x[i];
  }
  
  for(int m = 0; m < n; m++) {  
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][0] = h * f[i]; //k1 para todos los f[i]
      x[i] = x[i] + 0.5 * k[i][0];
    }
    // printf("%lf\t%lf\t%lf\n", f[0], f[1], f[2]);
    rhs1(x, f); //t0
    for (int i = 0; i < NVARS; i++) {
      k[i][1] = h * f[i]; //k2 para todos los f[i]
      x[i] = x[i] + 0.5 * k[i][1];
    }
    // printf("%lf\t%lf\t%lf\n", f[0], f[1], f[2]);
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][2] = h * f[i]; //k3 para todos los f[i]
      x[i] = x[i] + k[i][2];
    }
    // printf("%lf\t%lf\t%lf\n", f[0], f[1], f[2]);
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][3] = h * f[i];
      res[i] = res[i] + (k[i][0] + 2*k[i][1] + + 2*k[i][2] + k[i][3])/6;
      x[i] = res[i];
    }
    printf("%lf\t%lf\t%lf\n", t0, x[0], x[1]);
    t0 = t0 + h; 
  }
}

int main (int argc, char *argv[]) {
  double t0, tf, h, n;
  double results[NVARS];
  double x0[NVARS] = {1,0};

  t0 = 0;
  tf = 6.28;
  n = 500;

  h = (tf - t0) / n;

  rk4sys(results, t0, x0, h, n);
  return 0;
}
