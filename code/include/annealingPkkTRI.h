//Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "./graph_library.h"
#include "./rewiring.h"
#include "./annealingTRI.h"
#include "../main.h"

//----- HEADERS --------//

int    rewiring_PkkTRI_annealing(GRAPH G,double B,double increment,double accmin,int rewires,gsl_rng* randgsl);
