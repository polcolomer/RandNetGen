//Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>


//----- INCLUDES --------//
#include "graph_library.h"
#include "annealingCk.h"
#include "annealingCbar.h"
#include "annealingTRI.h"
#include "annealingPkkCk.h"
#include "annealingPkkCbar.h"
#include "annealingPkkTRI.h"
#include "annealingKnn.h"

//----- HEADERS --------//
int read_arguments(char netNAME[],int*rewires,double* beta0,double* Abeta,double* accMIN,int* pkk,double* dk,char ck[],char cbar[],char tri[],char knn[],int* seed,int argc,char * const  argv[]);


#endif

