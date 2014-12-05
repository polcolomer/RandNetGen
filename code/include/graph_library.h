//Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

//--- INCLUDES ---------//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>
#include "../main.h"

//----- HEADERS --------//

GRAPH	read_network(char netNAME[]);
EDGE* 	create_edges(GRAPH G);
void	print_network (GRAPH G,char *nom,int init_node);
int		free_graph (GRAPH G);
int*    degree_distribution(GRAPH G);
double* clustering_spectrum (GRAPH G);
double  clustering_coeff (GRAPH G);
double  numOFtrianglesXnode (GRAPH G);
GRAPH   clustering(GRAPH G);
double* read_CK_fromFILE(GRAPH G,char fileNAME[]);
double* read_Knn_fromFILE(GRAPH G,char fileNAME[]);
double* Knn (GRAPH G);
