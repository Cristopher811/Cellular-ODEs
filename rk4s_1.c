#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NVARS 3

void rhs1(double *x, double *res) {
  res[0] = x[0]*x[1]*x[2];
  res[1] = x[1]+x[2]*x[0];
  res[2] = x[2]+x[0]*x[1];
}

void rk4sys(double *res, double t0, double x0, double h, int n){
  // f es la funcion evaluada en x0 en todo el sistema
  double f[NVARS], k[NVARS][4];
  double x[NVARS] = {x0, x0, x0}; // Iniciales del rk sobre x
  
  for(int m = 0; m < n; m++) {
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][0] = h * f[i];
      x[i] = x[i] + 0.5 * k[i][0];
    }
    // printf("%lf\t%lf\t%lf\n", f[0], f[1], f[2]);
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][1] = h * f[i];
      x[i] = x[i] + 0.5 * k[i][1];
    }
    // printf("%lf\t%lf\t%lf\n", f[0], f[1], f[2]);
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][2] = h * f[i];
      x[i] = x[i] + k[i][2];
    }
    // printf("%lf\t%lf\t%lf\n", f[0], f[1], f[2]);
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][3] = h * f[i];
      res[i] = res[i] + (k[i][0] + 2*k[i][1] + + 2*k[i][2] + k[i][3])/6;
      x[i] = res[i];
    }
    printf("%lf\t%lf\t%lf\t%lf\n", t0, x[0], x[1], x[2]); //? checar
    t0 = t0 + h;
  }
}

int main (int argc, char *argv[]) {
  double x0, t0, tf, h, n;
  double results[NVARS];

  x0 = 0.00128945;
  t0 = 0;
  tf = 10;
  n = 500;

  h = (tf - t0) / n;
  rk4sys(results, t0, x0, h, n);
  return 0;
}
