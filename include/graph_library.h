//Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

#ifndef GRAPH_LIBRARY_H
#define GRAPH_LIBRARY_H

//**********************************************
//       INCLUDES 
//**********************************************
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>

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
//       HEADERS 
//**********************************************

GRAPH	read_network(char netNAME[]);
void 	create_edges(GRAPH* G);
void	print_network (GRAPH* G,char *nom,int init_node);
int		free_graph (GRAPH* G);
void    degree_distribution(GRAPH* G);
double* clustering_spectrum (GRAPH* G);
double  clustering_coeff (GRAPH* G);
double  numOFtrianglesXnode (GRAPH* G);
double* read_CK_fromFILE(GRAPH* G,char fileNAME[]);
double* read_Knn_fromFILE(GRAPH* G,char fileNAME[]);
double* Knn (GRAPH* G);

#endif

