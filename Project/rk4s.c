#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "./Sources/rhs.h"



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

double rhs1(double *x, double *ed){
  double sum;
  double Rx_i[NVARS];

  r0(x,Rx_i); //EVALUATION OF DATA in Rx_i  -> Change this depending of the cell you want to use


//Constructing the system of differential equations.
  for(int i = 0; i<NVARS; i++){
    sum += Rx_i[i]; // sum of each Rx_i (REMEMBER THAT WE ONLY NEED THE FINAL VALUE OF SUM TO CONSTRUCT THE SYSTEM OF DIFF. EQ.).
  }

  for(int i = 0; i<NVARS; i++){
    ed[i] = Rx_i[i] - x[i]*sum; //ith differential equation for the values x[i]; 
  }
  return sum;
}

void rk4sys(double *res, double t0, double *x0, double h){
  double f[NVARS], k[NVARS][4]; // f is the differential equation evaluated at each x0
  double x[NVARS]; // Initial conditions of the system (PROVIDED BY x_0.dat)
  double sum;
  double integral = 0;
  double dt = 0.0001;
  //
  for(int i = 0; i < NVARS; i++){
    x[i] = x0[i];
    res[i] = x[i];
  }

  while(integral<ln2){
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][0] = h * f[i];
      x[i] = x[i] + 0.5 * k[i][0];
    }
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][1] = h * f[i];
      x[i] = x[i] + 0.5 * k[i][1];
    }
    rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][2] = h * f[i];
      x[i] = x[i] + k[i][2];
    }
    sum = rhs1(x, f);
    for (int i = 0; i < NVARS; i++) {
      k[i][3] = h * f[i];
      res[i] = res[i] + (k[i][0] + 2*k[i][1] + 2*k[i][2] + k[i][3])/6;
      x[i] = res[i];
    }
    t0++;
    integral += sum*dt;
    printf("\n%lf\t",t0*0.00001);
    printf("%0.18lf\n", integral); // Change this depending of the data that we want to plot (i.e. either integral or sum)
    
    }
    printf("\t====FINAL====\n%lf\t",t0*0.00001);  
    printf("%0.18lf\n",integral);  
 }

int main (int argc, char *argv[]) {
  double x0[NVARS], t0, tf, h, n;
  double res[NVARS];

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


  t0 = 0;
  h = 0.001;

  rk4sys(res, t0, x0, h);
  return 0;
}
