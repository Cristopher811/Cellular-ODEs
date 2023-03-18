#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define SIZE 250

int read_data_file(double *dataArray, FILE *file);

int main() {
	FILE *cellData = NULL;
	double cellConcentration[SIZE];

	// Abrimos y leemos el archivo
	cellData = fopen("./x_0.dat", "r");
	// Verificamos errores
	if (cellData == NULL) {
		fprintf(stderr, "Error al abrir el archivo x_0.dat");
		exit(EXIT_FAILURE);
	}
	read_data_file(cellConcentration, cellData);
	fclose(cellData);
	for (int i = 0; i < SIZE; i++)
		printf("x[%i] = %lf\n", i, cellConcentration[i]);
	return 0;
}

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
