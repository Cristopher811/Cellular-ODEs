#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQR(X)  (X)*(X)
#define NVARS 2

double line(double t, double x)
{
  double m,b;
  double result;

  m = 1;
  b = 0;
  result = m*x + b; // result := function

  return result;
}

double rhs0(double t, double x)
{
  double result;
  result = SQR(t) - SQR(x);
  return result;
}

double rhs1(double t, double x) 
{
  double result;
  result = 1 + SQR(x);

  return result;
}

double rhs2(double t, double x)
{
  double result;
  result = SQR(x);
  
  return result;
}

double ho(double t, double *x, int i)
{
  double result; //result := function

  switch (i) {
    case 0: result = x[1]; break; // \dot{x} = y
    case 1: result = -x[0]; break; // \dot{y} = -x
  }
  return result;
}

//METHODS
double euler (double (*f)(double, double), double t0, double x0, double h, int n, int verbose)
{
  double slope, x1;
  int i;

  if (verbose) {
    printf("\nt0\tx0\tslope\txn\n");
    printf("--------------------\n");
  }
  for (i=0; i<n; i++) {
    slope = f(t0,x0);
    x1 = x0 + h*slope;
    if (verbose) {
      printf("%lf\t%lf\t%lf\t%lf\n",t0,x0,slope,x1);
    }

    x0 = x1;
    t0 = t0 + h;
  }
  return x1;
}

/*Runge Kutta  Method 4th order*/
double rk4 (double (*rhs)(double, double), double t0, double x0, double h, int n, int verbose)
{
  double  x1,f,k1,k2,k3,k4;
  int i;

  for (i=0; i<n; i++) {
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
      printf("\n\n x(%lf) = %lf\n", t0+h,x1);
    }
    x0 = x1;
    t0 = t0+h;
  }
  return x1;
}

/*Runge Kuta Method 4th order with adaptive stepsize*/
double rk4adaptive (double (*rhs)(double, double), double t0, double tf, double x0, double tol)
{
  double x1, x1halfStep, h;
  int i, n=1;

  h = tf-t0;
  x1 = rk4(rhs,t0,x0,h,n,0);
  while(fabs(x1halfStep-x1)>tol)
  {
    x1 = x1halfStep;
    n *= 2;
    h *= 0.5;
    x1halfStep = rk4(rhs,t0,x0,h,n,0);
  }
  return x1halfStep;
}


/*Runge Kuta Method 4th order for a ODE system*/
void rk4sys (double (*f)(double, double*,int), double t, double *var, double step)
{
  double h=0.5*step;
  double tvar1[NVARS], tvar2[NVARS], tvar3[NVARS]; // temporary storage arrays
  double k1[NVARS], k2[NVARS], k3[NVARS], k4[NVARS]; // for RK
  int i;

  for (i=0; i<NVARS; i++) {
    tvar1[i]=var[i]+0.5*(k1[i]=step*(*f)(t,var,i));
  }
  for (i=0; i<NVARS; i++) {
    tvar2[i]=var[i]+0.5*(k2[i]=step*(*f)(t+h,tvar1,i));
  }
  for (i=0; i<NVARS; i++) {
    tvar3[i]=var[i]+0.5*(k3[i]=step*(*f)(t+h,tvar2,i));
  }
  for (i=0; i<NVARS; i++) {
    k4[i] = step*(*f)(t+step,tvar3,i);
  }

  for (i=0; i<NVARS; i++) {
    var[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
  }
}

int main (int argc, char *argv[])
{
  double t0, x0, tn, h, xn;
  double (*f)(double, double);
  double *y, tmp;
  int i, n, ex, v;

  sscanf(argv[1], "%d", &ex); //example 
  sscanf(argv[2], "%lf", &t0); //initial time
  sscanf(argv[3], "%lf", &x0); //initial value
  sscanf(argv[4], "%lf", &tn); // final time
  sscanf(argv[5], "%d", &n);  // number of steps
  sscanf(argv[6], "%d", &v); // verbose 0,1

  //Calculating step size (h)
  h = (tn - t0)/n;

  switch (ex) {
    case 1:
    f = line;//2.71828182846
    break;
    case 2:
    f = rhs0;
    break;
    case 3:
    f = rhs1;
    break;
    case 4:
    f = rhs2;
    break;

    case 5:
    y = (double *)malloc(NVARS*sizeof(double)); //dy/dx

      y[0] = x0;
      if(x0)
        y[1] = 0.0;
      else
        y[1] = 1.0;

      for (i=0; i<n; i++) {
        tmp = t0 + i*h;
        rk4sys(ho, tmp, y, h); // (ho:=rhs, tmp:=final time, y:lhs, h:step)
        printf("%lf\t%lf\t%lf\n", tmp, y[0], y[1]);
      }
    break;
  }
  if(ex!=5)
  {
    xn = euler(f,t0,x0,h,n,v);
    printf("Euler: x(%lf) = %lf\n", tn, xn);
    xn = rk4(f,t0,x0,h,n,v);
    printf("RK4: x(%lf) = %lf\n", tn, xn);
    xn = rk4adaptive(f,t0,tn,x0,1.e-7);
    printf("RK$(adaptive stepsize): x(%lf) = %lf\n", tn, xn);
  }
  return 0;
}


