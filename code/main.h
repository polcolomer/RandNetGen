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

//**********************************************
//       STRUCTURES 
//**********************************************

//----- A NODE --------//
typedef struct NODE{
	
	int idnum;	/// ID from the original file
	int k;		/// degree
	int *out;	/// vector wth its neighbours
	int tri;    /// number of triangles or common neighbours

		
}NODE;

//----- AN EDGE --------//
typedef struct EDGE{
	
	int s;	  /// source
	int d;	  /// destination
	
}EDGE;

//----- A GRAPH --------//
typedef struct GRAPH{
	
	int N;				/// number of nodes
	NODE* node;			/// vector with the nodes
    int E;				/// number of edges
	EDGE* edge;			/// vector with all the edges
    int loops;			/// number of self edges
	int linksWrepeat;	/// number of links with repetitions
	int max_k;			/// maximum degree
	int* pk;            /// degree distribution
    double* ck;         /// clustering spectrum
    double  Ccoef;      /// clustering coefficient
    double triangles;   /// number of triangles divided by N
    double* Knn;        /// the average neighbour degree

}GRAPH;
//**********************************************
//**********************************************

//----- INCLUDES --------//

#include "./include/graph_library.h"
#include "./include/annealingCk.h"
#include "./include/annealingCbar.h"
#include "./include/annealingTRI.h"
#include "./include/annealingPkkCk.h"
#include "./include/annealingPkkCbar.h"
#include "./include/annealingPkkTRI.h"
#include "./include/annealingKnn.h"
#include "./include/rewiring.h"

//----- HEADERS --------//
int read_arguments(char netNAME[],int*rewires,double* beta0,double* Abeta,double* accMIN,int* pkk,char ck[],char cbar[],char tri[],char knn[],int* seed,int argc,char * const  argv[]);

#endif

