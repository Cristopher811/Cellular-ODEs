#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "./Sources/rhs.h"

//IDEA, CHANGE THE ARGUMENT OF RHS1 IN ORDER TO BE ABLE TO USE T0 AS LIMIT OF INTEGRATION.

int read_data_file(double *dataArray, FILE *file) {
	char buffer[100] = {'\0'};
	char *end = NULL;
	char c = 0;
	int i = 0;
	do {
		c = fgetc(file);
		if (c == EOF || c == '\n') {
			dataArray[i] = strtod(buffer, &end);
			buffer[0] = '\0';
			i++;
		} else
			strncat(buffer, &c, 1);
	} while(c != EOF);
	return 1;
}

void rhs1(double *x, double *res) {
  double sum;
  double rhs[NVARS];

  //r1(x,rhs);
/*   printf("%lu\t%lu\t\n", rhs,x); //esto ok */

/*   for (int i = 0; i<NVARS; i++) { 
    printf("%lf\t\n",x[i]); no ok -nan
  } */



  for(int i = 0; i<NVARS; i++){
    res[i] = x[i]; 
    rhs[i] = res[i];
    sum += rhs[i];
  /*   printf("SUM[%d]\t%lf\t\n\n", i , sum); */
  }
  
  for(int i = 0; i<NVARS; i++){
    //res[i] = x[i]; //checar *_*  CRIS RECUERDA QUE ESTO LO MODIFICASTE!!!!!
    res[i] = rhs[i] - x[i]*sum;
/*     printf("RES\t%lf\n", res[i]); */
  }
}

void rk4sys(double *res, double t0, double *x0, double h, int n){
  double f[NVARS], k[NVARS][4]; // f es la funcion evaluada en x0 en todo el sistema
  double x[NVARS]; // Iniciales del rk sobre x
  //
  for(int i = 0; i < NVARS; i++){
    x[i] = x0[i];
    res[i] = x[i]; // checar *_* CRIS RECUERDA QUE ESTO LO MODIFICASTE!!!
  }
/*   printf("%lu", x0); */
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
      printf("%lf\t%d\t%lf\t\n", t0, i, x[i]);
    }
    t0 = t0 + h;
  }
 }

int main (int argc, char *argv[]) {
  double x0[NVARS], t0, tf, h, n;
  double results[NVARS];

	FILE *cellData = NULL;
	// Abrimos y leemos el archivo
	cellData = fopen("./Sources/x_0.dat", "r");
	// Verificamos errores
	if (cellData == NULL) {
		fprintf(stderr, "Error al abrir el archivo x_0.dat");
		exit(EXIT_FAILURE);
	} 
	read_data_file(x0, cellData);
	fclose(cellData);
	for (int i = 0; i < NVARS; i++){
	//	printf("dx[%i]/dt = %lf\n", i, x0[i]);
  }

/*   printf("%lu\t\n", &x0); */
  t0 = 0;
  tf = 10;
  n = 50;

  h = (tf - t0) / n;

  rk4sys(results, t0, x0, h, n);
  return 0;
}
