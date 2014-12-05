//Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

//--- INCLUDES ---------//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "./graph_library.h"
#include "../main.h"

//----- HEADERS --------//

int    swap_edges(GRAPH G, int s1,int s2,int r1,int r2);
int    rewiring_Pk(GRAPH G,int num_rewires,gsl_rng * randgsl);
int    choose_2_edges_random(GRAPH G,int *pos1,int *pos2,gsl_rng* randgsl);
int	   rewiring_Pkk(GRAPH G,int num_rewires,gsl_rng * randgsl);
int*   countEDGESwithK(GRAPH G);
int**  createPOSedgesK(GRAPH G,int *num_edges_k);
int    choose_2_edges_random_pkk(GRAPH G,int *pos_r,int *pos_s,int *num_edges_k,int **pos_edges_k,gsl_rng* randgsl);
