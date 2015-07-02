/************************************************************
 *
 *                    Graph Library
 *
 *		Functions useful when dealing with networks
 *
 * Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.
 *************************************************************/



#include "graph_library.h"

//**********************************************************************
//**********************************************************************
/// a function that reads and edge list and put it into the structures
GRAPH read_network(char netNAME[]){
	
	printf("Reading the network...\n");fflush(stdout);
	
    FILE *file = fopen(netNAME,"r");									///we open the file that was the argument
    
	int origen,desti,orig=0,dest=0;
	
	int N=0,nodes_saved=1000;
	GRAPH G;
	G.node   = (NODE*)malloc(sizeof(NODE)*nodes_saved);
	int *mem = (int*)malloc(sizeof(int)*nodes_saved);
	
	int org_new,dest_new;
	int j;
	int repeated;
	
	G.linksWrepeat=0;
	G.loops=0;
	G.max_k=0;
	
	
	while (!feof(file)){		///we start reading	
		
		repeated=0;
		if(fscanf(file,"%d\t%d\n",&origen,&desti)!=2)fprintf(stderr,"error reading a link form the edge list");///we read a link
		
		G.linksWrepeat = G.linksWrepeat+1;
		
		if(origen==desti){		///first we check if the link is a loop
			G.loops=G.loops+1;
		}
		
		else{
			
			if (origen < 0 ) printf("there is a negative G.node at line %i\n",G.linksWrepeat);
			if (desti  < 0 ) printf("there is a negative G.node at line %i\n",G.linksWrepeat);
			
			
			org_new=1; dest_new=1;	///primer assumim que el G.node es nou
			for(j=0;j<N;++j){
				
				if(G.node[j].idnum==origen){ org_new=0;  orig = j;}///Si no es nou ho anotem i apuntem el G.node que és
				if(G.node[j].idnum==desti) { dest_new=0; dest = j;}
				
				if (!dest_new && !org_new )break;	///si ja hem trobat els dos sortim
			}
			
			if(org_new){			///if the origin is a new G.node
				G.node[N].idnum = origen;
				mem   [N]       = 5;
				G.node[N].k     = 0;
				G.node[N].out   = (int*)malloc(sizeof(int)*mem[N]);
				origen = N;
				N++;
				if(N>nodes_saved-1){
					
					nodes_saved = nodes_saved*2;
					G.node = (NODE*)realloc(G.node,sizeof(NODE)*nodes_saved);
					mem    = (int*)realloc(mem,sizeof(int)*nodes_saved);
				}
				
			}
			else {origen=orig;}
			
			if(dest_new){		///if the destination is a new G.node
				G.node[N].idnum = desti;
				mem   [N]     = 5;
				G.node[N].k   = 0;
				G.node[N].out = (int*)malloc(sizeof(int)*mem[N]);
				desti = N;
				N++;
				if(N>nodes_saved-1){
					
					nodes_saved = nodes_saved*2;
					G.node = (NODE*)realloc(G.node,sizeof(NODE)*nodes_saved);
					mem    = (int*)realloc(mem,sizeof(int)*nodes_saved);
				}
				
			}
			else{desti=dest;}
			
			repeated=0;
			for (j=0;j<G.node[origen].k;++j){			///We check if the link was already there
				
				if (G.node[origen].out[j]==desti){	///if the link was already there we just get_out
					repeated=1;
					break;
				}
				
			}
			
			///if it is not a loop and it is not repated we save the link
			if(repeated==0){
				
				
				G.node[origen].k = G.node[origen].k +1;
				
				if(G.node[origen].k>mem[origen]-1){
					
					mem[origen]=mem[origen]*2;
					G.node[origen].out = (int*)realloc(G.node[origen].out,sizeof(int)*mem[origen]);
				}
				
				G.node[origen].out[G.node[origen].k-1] = desti;
				
				G.node[desti].k = G.node[desti].k +1;
				
				if(G.node[desti].k>mem[desti]-1){
					
					mem[desti]=mem[desti]*2;
					G.node[desti].out = (int*)realloc(G.node[desti].out,sizeof(int)*mem[desti]);
				}
				
				G.node[desti].out[G.node[desti].k-1] = origen;
				
			}
		}
		
	}
	
	G.E=0;
	for(j=0;j<N;++j){
		if(G.node[j].k>G.max_k) G.max_k = G.node[j].k;
		G.E = G.E + G.node[j].k;
	}
	G.E = G.E/2;
	G.N = N;
	
	free(mem);
    fclose(file);
	
	return G;
}
//**********************************************************************
//**********************************************************************
///a function that creates a matrix with all the edges without repetitions
void create_edges(GRAPH* G){
	
	G->edge = (EDGE*)malloc(sizeof(EDGE)*G->E);

	int i,j,l=0;
	
	for(i=0; i<G->N; ++i){
		for(j=0; j < G->node[i].k; ++j){
			
			if( i > G->node[i].out[j]){	///for not repeat
				
				G->edge[l].s = i;
				G->edge[l].d = G->node[i].out[j];
				l++;
			}
		}
	}

	if(l!=G->E){fprintf(stderr,"error at crear edges!! l=%i E=%i\n",l,G->E);exit(-1);}
	
	
	return;
	
}
//**********************************************************************
//**********************************************************************
///function to print the hole network in a link list of two columns
void print_network (GRAPH* G,char *nom,int init_node){
	
    FILE *fitxer = fopen(nom,"w");
	
    int i,j;
	for(i=0;i<G->N;++i){

		for (j=0;j<G->node[i].k;j++){
			
			fprintf(fitxer,"%i\t%i\n",i+init_node,G->node[i].out[j]+init_node);
		}
	}
	fclose(fitxer);
	
	return;
}

//**********************************************************************
//**********************************************************************
/// a function to free the memory saved for a network
int free_graph (GRAPH* G){
	
	int i;
	for( i=0; i<G->N; ++i){
		
		free(G->node[i].out);
	}
	
	free(G->node);
	
	return 1;
}
//**********************************************************************
//**********************************************************************
///function to compute the p(k). we will measure it as p(k)=N_k/N_tot
void degree_distribution(GRAPH* G){
	
	G->pk = (int*)calloc(G->max_k+1,sizeof(int));
	
    int i;
	for ( i=0; i<G->N; ++i){	///now we count how many nodes of a certain degree we have
		
		G->pk[G->node[i].k] = G->pk[G->node[i].k]+1;
		
	}
	
	return;
}
//**********************************************************************
//**********************************************************************
///function to measure the clustering spectrum C(k)
double* clustering_spectrum (GRAPH* G){
	
	double* ck = malloc(sizeof(double)*(G->max_k+1));
	int i=0,j=0,k=0,l=0,e=0,tri=0;
    double Ccoef=0;
	
	for(i=0;i<G->max_k+1;++i){	///we put to zero the vectors c_k and  number_k
		ck[i]=0;
	}
	
	///for each node i
	for(i=0; i<G->N; ++i){
		
		e=0;
		G->node[i].tri = 0;
		
		for(j=0; j<G->node[i].k; ++j){
			for(k=j+1; k<G->node[i].k; ++k){
				for (l=0;l<G->node[G->node[i].out[j]].k; ++l){
					if(G->node[G->node[i].out[j]].out[l]==G->node[i].out[k]) e=e+1; ///if there is a common neighbour we add one triagle to the node
						
					
				}			
			}
		}
		
		if(e!=0){
			tri = tri+e;
			G->node[i].tri = e;
            Ccoef = Ccoef + (2.*e)/(G->node[i].k*(G->node[i].k-1.));
			ck[G->node[i].k]= ck[G->node[i].k] + (double)e*2./G->node[i].k/(G->node[i].k-1.);
		}
		
	} 
	
		
	///and we take the mean without the graph with k=0,1
	for (i=0;i<G->max_k+1;++i){
		if (G->pk[i]>0) {
			ck[i] = ck[i]/G->pk[i];
		}
	}
	
	return ck;
}
//**********************************************************************
//**********************************************************************
///function to measure the LOCAL clustering coefficient c 
double clustering_coeff (GRAPH* G){
	
	int i=0,j=0,k=0,l=0,e=0,tri=0;
	double C=0;
	
	///for per cada node i
	for(i=0; i<G->N; ++i){
		
		e = 0;
		G->node[i].tri = 0;
		
		for(j=0; j<G->node[i].k; ++j){
			for(k=j+1; k<G->node[i].k; ++k){
				for (l=0;l<G->node[G->node[i].out[j]].k; ++l){
					if(G->node[G->node[i].out[j]].out[l]==G->node[i].out[k]) e=e+1; ///if there is a common neighbour we add one triagle to the node
						
					
				}			
			}
		}
		
		if(e!=0){
			tri = tri + e;
			G->node[i].tri = e;
			C = C + (2.*e)/(G->node[i].k*(G->node[i].k-1.));
		}
		
	} 

    return C/(G->N-G->pk[1]);
}
//**********************************************************************
//**********************************************************************
///function to measure the number of triangles
double numOFtrianglesXnode (GRAPH* G){
	
	int i=0,j=0,k=0,l=0,e=0,tri=0;
		
	///for per cada node i
	for(i=0; i<G->N; ++i){
		
		e = 0;
		G->node[i].tri = 0;
		
		for(j=0; j<G->node[i].k; ++j){
			for(k=j+1; k<G->node[i].k; ++k){
				for (l=0;l<G->node[G->node[i].out[j]].k; ++l){
					if(G->node[G->node[i].out[j]].out[l]==G->node[i].out[k]) e=e+1; ///if there is a common neighbour we add one triagle to the node
						
					
				}			
			}
		}
		
		if(e!=0){
			tri = tri + e;
			G->node[i].tri = e;
        }
		
	} 

    return (tri/3.)/G->N;
}

//**********************************************************************
//**********************************************************************
///function that reads the clustering spectrum froma file
double* read_CK_fromFILE(GRAPH* G,char fileNAME[]){
    
    FILE *file = fopen(fileNAME,"r");
    
    double* ckVEC  = calloc((G->max_k+1),sizeof(double));
    
    int k;
    double ck;
    
    while (!feof(file)){
        
        if(fscanf(file,"%d\t%lf\n",&k,&ck)!=2) printf("error reading the clustering spectrum from file");
        
        if(k>G->max_k) printf("the clustering spectrum file has larger degrees than our network. We omit them\n");
        else{  if(G->pk[k]>0) ckVEC[k] = ck;}
            
        
    }
    fclose(file);
        
    return ckVEC;
}
//**********************************************************************
//**********************************************************************
///function that read the average neighobur degree from a file
double* read_Knn_fromFILE(GRAPH* G,char fileNAME[]){
    
    FILE *file = fopen(fileNAME,"r");
    
    double* knnVEC  = calloc((G->max_k+1),sizeof(double));
    
    int k;
    double knn;
    
    while (!feof(file)){
        
        if(fscanf(file,"%d\t%lf\n",&k,&knn)!=2) printf("error reading the clustering spectrum from file");
        
        if(k>G->max_k) printf("the clustering spectrum file has larger degrees than our network. We omit them\n");
        else{  if(G->pk[k]>0) knnVEC[k] = knn;}
            
        
    }
    fclose(file);
        
    return knnVEC;
}
///function to measure the k_nn(k)
//**********************************************************************
// *********************************************************************
/// function that calculates the average neigbour degree fo a network
double* Knn (GRAPH* G){
	
		
	double *k_nn = (double*)calloc(G->max_k+1,sizeof(double));
	
	///we initialize the vectors
	int i,j;
	for(i=0;i<G->max_k+1;++i){k_nn[i]=0; }
	
	///we sum the degree all the neighbours of each node of degree k
	for (i=0;i<G->N;++i){
		
		for(j=0;j<G->node[i].k;++j){
			
			k_nn[G->node[i].k] = k_nn[G->node[i].k]+G->node[G->node[i].out[j]].k;
			
		}
	}
	
	///and we take the mean
	for (i=0;i<G->max_k+1;++i){
		
    	if(G->pk[i]>0){
		
			k_nn[i] = k_nn[i]/(G->pk[i]*i);
		}
	}
	
		
	return k_nn;
	
}
//**********************************************************************
// *********************************************************************
