//FICHERO TEMPORAL HASTA REALIZAR LOS EXAMPLES
//Rellena la variable eopdata con la información del fichero
#include "EOPDATA.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void leerFichero(double *eopdata[13])
{
  int n = 20026;
	FILE *fich = fopen("PreliminaryOrbitDetermination/eop19620101.txt", "r");
	assert(fich != NULL);
	for (int i = 0; i < 13; i++)
	{
		eopdata[i] = malloc(sizeof(double)*n);
	}
	int j = 0;
	while(fscanf(fich, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
				&eopdata[0][j], &eopdata[1][j], &eopdata[2][j], &eopdata[3][j],
				&eopdata[4][j], &eopdata[5][j], &eopdata[6][j], &eopdata[7][j],
				&eopdata[8][j], &eopdata[9][j], &eopdata[10][j], &eopdata[11][j],
				&eopdata[12][j]) != EOF)
	{
		j++;
	}
	fclose(fich);
}

//USO:
//double *eopdata[13];
//leerFichero(eopdata);
