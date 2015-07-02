//Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

#ifndef ANEALING_PKKCBAR_H
#define ANEALING_PKKCBAR_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "graph_library.h"
#include "rewiring.h"

//----- HEADERS --------//

int    rewiring_PkkCbar_annealing(GRAPH* G,double B,double increment,double accmin,int rewires,gsl_rng* randgsl);
double calc_AH_PkkCbar(GRAPH* G,int s1,int s2,int r1,int r2,double C0,double *C1,double Caim);

#endif


