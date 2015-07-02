/************************************************************
 *
 *                   REWIRING LIBRARY
 *
 *		Functions useful when rewiring a network
 *
 * Copyright 2014 Pol Colomer de Simon. All rights reserved. 
 * Code under License GPLv3.
 *
 *************************************************************/



#include "rewiring.h"

//**********************************************************************
//**********************************************************************
///a function that do a rewiring preserving the degree distribution
int rewiring_Pk(GRAPH* G,int num_rewires,gsl_rng * randgsl){
	
    printf("Rewiring preserving P(k)...\n");fflush(stdout);
    
	int i,pos_r,pos_s,r1,r2,s1,s2;
	
	for(i=0;i<num_rewires*G->E;++i){

		while(!choose_2_edges_random(G,&pos_r,&pos_s,randgsl)){} ///we try to find two edges avoiding selfedges and multipledges

		r1 = G->edge[pos_r].s;		/// the nodes first link are
		r2 = G->edge[pos_r].d;
		     
		s1 = G->edge[pos_s].s;		/// the nodes second link are
		s2 = G->edge[pos_s].d;
		
		swap_edges(G,s1,s2,r1,r2);	/// we swap the two links

		G->edge[pos_r].d = s2;		/// we modify the edge vector
		G->edge[pos_s].d = r2;
		
	}
	
	return 1;
}
//**********************************************************************
//**********************************************************************
///a function that chooses two edgess at random that if they are swaped we avoid multipledges or selfedges.
/// we try to choose a second link E times. if we dont find it we return a 0
int choose_2_edges_random(GRAPH* G,int *pos_r,int *pos_s,gsl_rng* randgsl){
	
	int j,rand1,rand2=0,r1,r2,s1,s2,correct,tries;

	rand1 = gsl_rng_uniform_int(randgsl,G->E);	///we choose an edges at random
	
	r1 = G->edge[rand1].s;
	r2 = G->edge[rand1].d;
	
	if(r1==r2){printf("there is a selflink!!!\n");fflush(stdout);}
	
	correct=0;	///in order to be able to enter into the loop 
	
	tries = 0;	///a variable that controls how many times we have tried to choose another link

	while(!correct && tries < G->E){///we choose a second edges at random with some constrains
		
		rand2 = gsl_rng_uniform_int(randgsl,G->E);
		
		s1 = G->edge[rand2].s;
		s2 = G->edge[rand2].d;
		
		if(s1==s2){printf("there is a selflink!!!\n");fflush(stdout);}
		
		correct=1;											/// we assume that it is correct
		
		if(s1==r1 || s1==r2 || s2==r1 || s2==r2)correct=0;	/// we look if all the for nodes are different
		
		else {
			for(j=0;j<G->node[s1].k;++j){						 
				if(G->node[s1].out[j]==r2) correct=0;		/// we look if the node s1 is not already connected with r2 to avoid a multiple link
			}
			for(j=0;j<G->node[r1].k;++j){						
				if(G->node[r1].out[j]==s2) correct=0;		/// we look if the node r1 is not already connected with s2 to avoid a multiple link
			}
		}
		
		tries++;
		
	}

	*pos_r = rand1;
	*pos_s = rand2;
	
	if(tries == G->E) return 0;	///we could not fins a proper pair of links we return a 0. 
	else return 1;
	
}
//**********************************************************************
//**********************************************************************
///a function that swap to edges from the pairs s1,s2 and r1,r2 to s1,r2 and r1,s2
int swap_edges(GRAPH* G, int s1,int s2,int r1,int r2){
	
	int j,posr1=0,posr2=0,poss1=0,poss2=0;
	
	for(j=0;j<G->node[r2].k;++j){if(G->node[r2].out[j]==r1) {posr1=j;break;}}	///we look where the nodes were in the neighbours vector
	if(j==G->node[r2].k){printf("error swap_edges. neighbour not found\n");abort();}
	
	for(j=0;j<G->node[s2].k;++j){if(G->node[s2].out[j]==s1) {poss1=j;break;}}
	if(j==G->node[s2].k){printf("error swap_edges. neighbour not found\n");abort();}
	
	for(j=0;j<G->node[r1].k;++j){if(G->node[r1].out[j]==r2) {posr2=j;break;}}
	if(j==G->node[r1].k){printf("error swap_edges. neighbour not found\n");abort();}
	
	for(j=0;j<G->node[s1].k;++j){if(G->node[s1].out[j]==s2) {poss2=j;break;}}
	if(j==G->node[s1].k){printf("error swap_edges. neighbour not found\n");abort();}

	
	G->node[r1].out[posr2]=s2;
	G->node[r2].out[posr1]=s1;
	
	G->node[s1].out[poss2]=r2;
	G->node[s2].out[poss1]=r1;
	
	return 1;
}
//**********************************************************************
//**********************************************************************
///a function that do a rewiring preserving the joint degree distribution P(k,k')
int rewiring_Pkk(GRAPH* G,int num_rewires,gsl_rng * randgsl){
	
    int* numEDGESwithK = countEDGESwithK(G); /// we count how many edges with and node of degree k are

	int **pos_edges_k = createPOSedgesK(G,numEDGESwithK);
	
	int i,s1,s2,r1,r2,kr1,kr2,ks1,ks2,pos_r,pos_s,j=0,l=0;					/// rewiring variables that will store the proposed rewiring
	
    /******** we start the rewiring *************/
	
	printf("Rewiring preserving P(k,k')...\n");fflush(stdout);
	

	for(i=1; i<num_rewires*G->E ;++i){
		
		
		while(!choose_2_edges_random_pkk(G,&pos_r,&pos_s,numEDGESwithK,pos_edges_k,randgsl)){} /// we try to find two edges avoiding selfedges and multipledges
		
		r1 = G->edge[pos_r].s;		/// the nodes and its degree are
		r2 = G->edge[pos_r].d;
		
        kr1 = G->node[r1].k;
		kr2 = G->node[r2].k;
        
        s1 = G->edge[pos_s].s;
		s2 = G->edge[pos_s].d;
        
        ks1 = G->node[s1].k;
		ks2 = G->node[s2].k;
		
		swap_edges(G,s1,s2,r1,r2);		/// we make the proposed rewired
		
		
		G->edge[pos_r].d = s2;			/// we modify the edge vector
		G->edge[pos_s].d = r2;
		
				
		if(kr2!=ks2){					///we modify the vector pos_edges_k
			
			if(kr2!=kr1){
				for(j=0;j<numEDGESwithK[kr2];++j){if(pos_edges_k[kr2][j]==pos_r){break;}}
			}
			if(ks2!=ks1){
				for(l=0;l<numEDGESwithK[ks2];++l){if(pos_edges_k[ks2][l]==pos_s){break;}}
			}
			
			
			if     (kr2!=kr1 && ks2==ks1)  pos_edges_k[kr2][j] = pos_s;    
			else if(ks2!=ks1 && kr2==kr1)  pos_edges_k[ks2][l] = pos_r;
			else if(kr2!=kr1 && ks2!=ks1){									/// in case the two other nodes have different degree we just swap the positions of the edges
				pos_edges_k[kr2][j] = pos_s;  
				pos_edges_k[ks2][l] = pos_r;
			}
			
			j=0;
			l=0;
			
		}
		
		
	}
    
	free(numEDGESwithK);
	for(i=1;i<G->max_k+1;++i){free(pos_edges_k[i]);}
    free(pos_edges_k);
	
	return 0;
}
//**********************************************************************
//**********************************************************************
/// a function that counts how many edges are with one node with each degree
int* countEDGESwithK(GRAPH* G){
	
    int* numEDGESwithK = (int*)calloc(G->max_k+1,sizeof(int));
	int i,n1,n2,k1,k2;
	
	for(i=0;i<G->E;++i){			/// we look to all edges
		
		n1 = G->edge[i].s;		/// The nodes of the edge
		n2 = G->edge[i].d;
        
        k1 = G->node[n1].k;      /// we take both degree
        k2 = G->node[n2].k;		
			
        numEDGESwithK[k1] = numEDGESwithK[k1] + 1;	/// we add one edge with degree k
        
        if(k1!=k2) numEDGESwithK[k2] = numEDGESwithK[k2] + 1; ///we also add for the neighbour just in case they do not have the same degree
			
    }
	
	return numEDGESwithK;
}

//**********************************************************************
//**********************************************************************
///a function created the vector pos_edges_k. This vector tells us which edges have at least one node of degree equal to the first component.
///hence pos_edge_k[i][j] gives the name of the edge the has at least one node of degree i
int** createPOSedgesK(GRAPH* G,int *numEDGESwithK){
	
	int i;
	
    int **pos_edges_k = (int**)malloc(sizeof(int*)*(G->max_k+1)); /// we save space for each degree
    						
	for(i=1;i<G->max_k+1;++i){pos_edges_k[i] = (int*)malloc(sizeof(int)*numEDGESwithK[i]);}	/// we save space for each edge with at least one node of degree i
	
	int *pos_k = (int*)calloc(G->max_k+1,sizeof(int));										/// a vector that tell us where is the last edge that we putted in the pos_edges_k matrix
	
	int r1,r2,kr1,kr2;
	
	for(i=0;i<G->E;++i){
		
		r1  = G->edge[i].s;						/// we take the nodes and its degree of the edge we are loking for
		r2  = G->edge[i].d;
		kr1 = G->node[r1].k;
		kr2 = G->node[r2].k;
		
		pos_edges_k [ kr1 ] [ pos_k[kr1] ] = i;	/// we put the position of the edge to the matrix
		pos_k[kr1] = pos_k[kr1] + 1;			/// we add the number off edges that we have putted into the matrix
		
		if(kr1!=kr2){							/// in case they have different degree we also add the edge to the degree row of the degree of the second vector
			pos_edges_k[kr2][pos_k[kr2]] = i;
			pos_k[kr2] = pos_k[kr2] + 1;
		}
		
	}
	
	for(i=1;i<G->max_k+1;++i){
		
		if(pos_k[i] != numEDGESwithK[i]){ printf("problem at create_pos_edges_k\n");abort(); }
		
	}
	
	
	return pos_edges_k;
	
}

//**********************************************************************
//**********************************************************************
///a function that chooses two edgess at random that if they are swaped we avoid multipledges or selfedges.
/// two nodes ahve the same degree at least. we try to choose a second link E times. if we dont find it we return a 0
int choose_2_edges_random_pkk(GRAPH* G,int *pos_r,int *pos_s,int *numEDGESwithK,int **pos_edges_k,gsl_rng* randgsl){
	
	int j,rand1,rand2=0,r1,r2,s1,s2,kr1,ks1,ks2,correct,tries;
	
	rand1 = gsl_rng_uniform_int(randgsl,G->E);	///we choose an edges at random
	
	r1 = G->edge[rand1].s;		///the nodes are these ones
	r2 = G->edge[rand1].d;
	
	if(gsl_rng_uniform_int(randgsl,2)){	///with a 50% we took one node or the another. in case we choose the second one we flip them in order to be everything the same
		
		G->edge[rand1].s = r2;
		G->edge[rand1].d = r1;
		r1 = G->edge[rand1].s;
		r2 = G->edge[rand1].d;
		
	}
	
	
	kr1 = G->node[r1].k;		  /// we take the degrees of the nodes
		
	if(G->pk[kr1]==1) return 0; ///if we only have one node with this degree we go out and choose another edge
	
	
	correct=0;	///in order to be able to enter into the loop 
	
	tries = 0;	///a variable that controls how many times we have tried to choose another link. We will skip after E tries
	
	while(!correct && tries < G->E){///we choose a second edges at random with some constrains
		
		rand2 = rand1; ///to enter into the bucle. we can not choose the same edge
		while (rand2==rand1){				
			rand2 = gsl_rng_uniform_int(randgsl,numEDGESwithK[kr1]); ///we choose a random edge with at least one degree as kr1
		
			rand2 = pos_edges_k[kr1][rand2];/// we look which is the edge in the edge list
		}
		
		s1 = G->edge[rand2].s;				/// now we know the nodes of the edges
		s2 = G->edge[rand2].d;
		
		ks1 = G->node[s1].k;					/// we take the degrees of them
		ks2 = G->node[s2].k;					
			
		if((ks1!=kr1) || ((ks1==kr1 && ks2==kr1) && gsl_rng_uniform_int(randgsl,2))){	///in case the s1 is the node with the same degree we have to swap the edge. in case s1 s2 have the same with 50%
			
			G->edge[rand2].s = s2;			/// we swap the edge
			G->edge[rand2].d = s1;
			s1 = G->edge[rand2].s;			/// now the nodes are
			s2 = G->edge[rand2].d;
			ks1 = G->node[s1].k;				/// we take the degrees of them
			ks2 = G->node[s2].k;
			
		}
		
		
		correct=1;											/// we assume that it is correct
		
		if(s1==r1 || s1==r2 || s2==r1 || s2==r2) correct=0;	/// we look if all the for nodes are different
		
		else {
			
			for(j=0; j < ks1 ;++j){
				
				if(G->node[s1].out[j]==r2)correct=0;			/// we look if the node s1 is not already connected with r2 to avoid a multiple link
								
			}
			for(j=0; j < kr1 ;++j){						
				if(G->node[r1].out[j]==s2) correct=0;		/// we look if the node r1 is not already connected with s2 to avoid a multiple link
								
			}
		}

		tries++;
		
	}
	
	
	if(tries == G->E) return 0;	/// if we could not find a proper pair of links we return a 0. hope to be luckier next time
	
	*pos_r = rand1;				/// now we know the position of the two edges. 
	*pos_s = rand2;
	
	
	return 1;
	
}
//**********************************************************************
//**********************************************************************

