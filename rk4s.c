#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NVARS 3

double rhs1(double t, double *x, int i){

  double result;

  switch(i){
    case 0: result = 0-x[2]*x[0]-x[1]*x[0]-x[2]*x[0]+4.000000*x[2]*(0.200000-x[0]); break;
    case 1: result = 0-x[1]*x[1]-x[1]*x[1]-x[2]*x[1]+4.000000*x[2]*(0.200000-x[1]); break;
    case 2: result = 0-x[2]*x[1]-x[2]*x[0]-x[1]+4.000000*x[0]*(0.200000-x[2]); break;
  }
  printf("\t%lf\t%lf\t%lf\t\n",x[0],x[1],x[2]);
  return result;
}


double rk4(double (*rhs)(double, double), double t0, double x0, double h, int n, int verbose){
  double x1,f,k1,k2,k3,k4;
  int i;

  for (i = 1; i < n; i++) {
    f = rhs(t0,x0);
    k1 = h*f;
    f = rhs(t0+0.5*h,x0+0.5*k1);
    k2 = h*f;
    f = rhs(t0+0.5*h,x0+0.5*k2);
    k3 = h*f;
    f = rhs(t0+h,x0+k3);
    k4 = h*f;
    x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;

    if (verbose) {
      printf("\n\n k1 = %lf ", k1);
      printf("\n\n k2 = %lf ", k2);
      printf("\n\n k3 = %lf ", k3);
      printf("\n\n k4 = %lf ", k4);
      printf("\n\n x(%lf) = %lf\n ", t0+h,x1);
    }
    
    x0 = x1;
    t0 = t0+h;
  }
  return x1;
}


double rk4adpatative(double (*rhs)(double, double), double t0, double tf, double x0, double tol){
  double x1, x1halfstep, h;
  int i, n=1;

  h = tf-t0;
  x1 = rk4(rhs,t0,x0,h,n,0);
  n *=2;
  h *= 0.5;
  x1halfstep = rk4(rhs,t0,x0,h,n,0);
  
  while (fabs(x1halfstep-x1)>tol) {
    x1 = x1halfstep;
    n *= 2;
    h*= 0.5;
    x1halfstep = rk4(rhs,t0,x0,h,n,0);
  }
  return x1halfstep;
}


void rk4system(double (*f)(double ,double* ,int), double t, double *var, double step){
  double h = 0.5*step;
  double tvar1[NVARS], tvar2[NVARS], tvar3[NVARS];
  double k1[NVARS], k2[NVARS], k3[NVARS], k4[NVARS]; 
  int i;

  for (i=0;i<NVARS; i++) tvar1[i]=var[i]+0.5*(k1[i]=step*(*f)(t, var, i));
  for (i=0;i<NVARS; i++) tvar2[i]=var[i]+0.5*(k2[i]=step*(*f)(t, tvar1, i));
  for (i=0;i<NVARS; i++) tvar3[i]=var[i]+0.5*(k3[i]=step*(*f)(t, tvar2, i));
  for (i=0;i<NVARS; i++) k4[i]=step*(*f)(t+step,tvar3,i);
  for (i=0;i<NVARS; i++) var[i]+=(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
}

int main (int argc, char *argv[])
{
  double t0,x0,tn,h,xn;
  double (*f)(double,double);
  double *y,tmp;
  int i, n, v;

  sscanf(argv[1], "%lf", &t0);
  sscanf(argv[2], "%lf", &x0);
  sscanf(argv[3], "%lf", &tn);
  sscanf(argv[4], "%d", &n);
  sscanf(argv[5], "%d", &v);

  h = (tn - t0)/n;


  y = (double *)malloc(NVARS*sizeof(double));

  y[0] = x0;
  y[1] = x0;
  y[2] = x0;

  printf("\t%lf\t%lf\t%lf\n", y[0],y[1],y[2]);

  
  for(i=0;i<n;i++){
    tmp = t0+i*h;
    rk4system(rhs1,tmp,y,h);
    printf("\t%lf\t%lf\t%lf\t%lf\n", tmp, y[0],y[1],y[2]);
  }
  return 0;
}
