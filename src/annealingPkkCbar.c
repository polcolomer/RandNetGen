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



#include "annealingPkkCbar.h"

//**********************************************************************
//**********************************************************************
int rewiring_PkkCbar_annealing(GRAPH* G,double B,double increment,double accmin,int rewires,gsl_rng* randgsl){
	
    
   	double Caim = G->Ccoef;  /// this is the clustering we want to achieve
	
	
 /***********************************************************************
	 we calculate the initial clustering of the random network
 ************************************************************************/
	
	double C  = clustering_coeff(G);
    
    double Cnew = C;
    
	printf("Caim %f Cinitial %f\n",Caim,C);
    
 /***********************************************************************
	       we do the rewiring preserving the P(k,k')
 ************************************************************************/	
    
    int* numEDGESwithK = countEDGESwithK(G); /// we count how many edges with and node of degree k are

	int **pos_edges_k = createPOSedgesK(G,numEDGESwithK);
	
	int s1,s2,r1,r2,kr1,kr2,ks1,ks2,pos_r,pos_s,j=0,l=0;					/// rewiring variables that will store the proposed rewiring
	
	double p,AH;
	
	int accepted=0,rewirestemp=0,pirem=0;					/// during the proces we count how many proposals rewirings are with AH>0, AH<o and AH=0
	double averAH = 0,averAHneg = 0,averAHpos = 0,oldacc=0;
	int numAH0 = 0,numAHneg = 0,numAHpos = 0;
	
    double H = fabs(C - Caim);                              /// initial energy
	
	
	/******** we start the rewiring *************/
	
	printf("Fixing the Clustering coefficient by an anneald rewiring preserving P(k,k')...\n");fflush(stdout);
	
	time_t start,end;										/// we will measure the time
	double dif;
	time (&start);
	
	FILE *ftemp = fopen("E_vs_T.dat","w");
	fprintf(ftemp,"#B\tEnergy\tacceptance\n");
	
    int i;
	for(i=1; oldacc>accmin || i<rewires*G->E+2 ;++i){
		
		
        /**** we propose a swap *********/			
		while(!choose_2_edges_random_pkk(G,&pos_r,&pos_s,numEDGESwithK,pos_edges_k,randgsl)){} /// we try to find two edges avoiding selfedges and multipledges
		
		r1 = G->edge[pos_r].s;		/// the nodes and its degree are
		r2 = G->edge[pos_r].d;
		
        kr1 = G->node[r1].k;
		kr2 = G->node[r2].k;
        
        s1 = G->edge[pos_s].s;
		s2 = G->edge[pos_s].d;
        
        ks1 = G->node[s1].k;
		ks2 = G->node[s2].k;
        
        /**** we calculate the increment of energy *********/
		AH = calc_AH_PkkCbar(G,s1,s2,r1,r2,C,&Cnew,Caim);	/// we calculate the increment of energy that would cause the rewiring
		
		averAH = averAH + fabs(AH);					    ///we also counbt the average AH of the proposals
		
		if(AH < 0.) {								    /// we count how many proposals have AH > 0
			
			numAHneg++;
			averAHneg = averAHneg + AH;
		}
		else if (AH > 0.) {							    /// we count how many proposals have AH < 0
			
			numAHpos++;
			averAHpos = averAHpos + AH;
			
		}
		else numAH0++;								/// we count how many proposals have AH = 0
		
		p =  gsl_rng_uniform(randgsl);				/// we throw a random number (0,1)
		
		/********** IF we acccept **************/
		
		if( p < exp(-B*AH) ){						
			
			
			swap_edges(G,s1,s2,r1,r2);			/// we make the proposed rewired
			
						
		    G->edge[pos_r].d = s2;						/// we modify the edge vector
            G->edge[pos_s].d = r2;
			
						
			if(kr2!=ks2){								///we modify the vector pos_edges_k
			
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
		    	
				
			C = Cnew;	/// we update the clustering vector
			
			if(fabs(AH)>0.)	accepted++;					/// we count how many changes we accept
			
			H = H + AH;
		}
		
		/********** IF we reject **************/
		
		else {								
			
			Cnew = C;	/// we recover the old clustering vector
		}
		
		rewirestemp++;
		
		/********** we reduce the temperature and we check the acceptance **************/
		
		if(rewirestemp > rewires*G->E ) {	///we try to find the appropiate temperature in order to have the desire acceptation
			
			printf("acceptance rate = %f "                ,(double)accepted/(numAHneg+numAHpos));				/// the acceptance
			//printf("AH = %f "                 ,averAH/(numAHneg+numAHpos));							/// the average energy of the proposed swaps
			//printf("numAHneg = %f AHneg = %e ",(double)numAHneg/rewirestemp,averAHneg/numAHneg);	/// the proportion of negative energy swaps and the average
			//printf("numAHpos = %f AHpos = %e ",(double)numAHpos/rewirestemp,averAHpos/numAHpos);	/// the proportion of positive energy swaps and the average
			//printf("numAH0=%f "               ,(double)numAH0/rewirestemp);							/// the proportion of proposals that do not change the energy
			printf("Beta=%e Energy=%e\n"              ,B,H);												/// the temperature and the energy
			fflush(stdout);
			
			fprintf(ftemp,"%f\t%f\t%f\n",B,H,(double)accepted/(numAHneg+numAHpos));
			
			if( ((double)accepted/(numAHneg+numAHpos)) > oldacc && i > rewires*G->E + 2 ) pirem++;	/// in case we havethe acceptance has increased 10 times the rewiring proces
			if(pirem>30) break;
			
			oldacc		= ((double)accepted/(numAHneg+numAHpos));								/// we save the old acceptance in order to compare with the next one
			accepted    = 0;												/// we put to zero all the acceptance counters
			rewirestemp = 0;
			averAH      = 0; numAHneg = 0;averAHneg = 0;
			numAH0      = 0; numAHpos = 0;averAHpos = 0;
			
			B = B*increment;												/// we reduce the temperature
			
		}
		
					
	}
	
	
	time (&end);									///we count the rewiring time and take conclusions
	dif = difftime (end,start);
	printf ("You rewired the entire network %.2f times and it took %.2lf seconds to run.\n",(double)i/G->E, dif );
	
 /***********************************************************************
		we print the network and we free the memory	
 ************************************************************************/
	
    printf("Cfinal %f\n",C);
		
	fclose(ftemp);
	
    free(numEDGESwithK);
	for(i=1;i<G->max_k+1;++i){free(pos_edges_k[i]);}
    free(pos_edges_k);
	
	
	
	return 0;
	
}
//**********************************************************************
//**********************************************************************
/// a function that calculates the difference in the clustering that the proposed change does
double calc_AH_PkkCbar(GRAPH* G,int s1,int s2,int r1,int r2,double C0,double *C1,double Caim){
	
    double Cnew = *C1;
    int norm = G->N-G->pk[1];
	int afectat,k; 			/// that variable stores the name of a common neighbour that also will change its cluster
	
	int ks1 = G->node[s1].k;	/// we look at the degrees of the afected nodes
	int ks2 = G->node[s2].k; 
	int kr1 = G->node[r1].k; 
	int kr2 = G->node[r2].k;
	
	int i,j;
	

    
	/******* we first discount the contribution of the links destroyed to the clustering ******/
	for(i=0; i<ks1; ++i){			
		for(j=0; j<ks2; ++j){
			
			if((G->node[s1].out[i] == G->node[s2].out[j])){   /// neighbours that had in common
				
				afectat  = G->node[s1].out[i];
                k        = G->node[afectat].k;
				
                Cnew     = Cnew - 2./(ks1*(ks1-1.))/norm;
                Cnew     = Cnew - 2./(ks2*(ks2-1.))/norm;
                Cnew     = Cnew - 2./(k*(k-1.))    /norm;
				
			}
			
		}
	}
	
		
	for(i=0; i<kr1; ++i){
		for(j=0; j<kr2; ++j){
			
			if(G->node[r1].out[i] == G->node[r2].out[j]){
				
                
                afectat  = G->node[r1].out[i];
                k        = G->node[afectat].k;
				
                Cnew     = Cnew - 2./(kr1*(kr1-1.))/norm;
                Cnew     = Cnew - 2./(kr2*(kr2-1.))/norm;
                Cnew     = Cnew - 2./(k*(k-1.))    /norm;
           
			}
			
		}
	}
	
	
	/*** and we add to the clustering the new common neighbours for the created links ******/
	
	for(i=0; i<ks1; ++i){			
		for(j=0; j<kr2; ++j){
			
			if(G->node[s1].out[i] == G->node[r2].out[j] && G->node[s1].out[i]!=s2 && G->node[s1].out[i]!=r1 ) {	/// we have to remember that the nodes s1 s2 and r1 r2 are no more neighbours!!

				afectat  = G->node[s1].out[i];
                k        = G->node[afectat].k;
				
                Cnew     = Cnew + 2./(ks1*(ks1-1.))/norm;
                Cnew     = Cnew + 2./(kr2*(kr2-1.))/norm;
                Cnew     = Cnew + 2./(k*(k-1.))    /norm;
			
			}
			
		}
	}
	
	for(i=0; i<kr1; ++i){
		for(j=0; j<ks2; ++j){
			
			if(G->node[r1].out[i] == G->node[s2].out[j] && G->node[r1].out[i]!=r2 && G->node[r1].out[i]!=s1){

				afectat  = G->node[r1].out[i];
                k        = G->node[afectat].k;
				
                Cnew     = Cnew + 2./(kr1*(kr1-1.))/norm;
                Cnew     = Cnew + 2./(ks2*(ks2-1.))/norm;
                Cnew     = Cnew + 2./(k*(k-1.))    /norm;
                				
			}
			
		}
	}
	
	/*** we calculate the energy increment ******/
    
    double AH;
    
    AH = fabs(Cnew-Caim) - fabs(C0-Caim);
    
    *C1 = Cnew;
		
	return AH;
}
