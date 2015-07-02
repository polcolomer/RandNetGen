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



#include "annealingCk.h"

//**********************************************************************
//**********************************************************************
int rewiring_Ck_annealing(GRAPH* G,double B,double increment,double accmin,int rewires,gsl_rng* randgsl){
	

 /***********************************************************************
	 we create a random network with the same degree sequence 
 ************************************************************************/

	double *ck  = clustering_spectrum(G);	/// we calculate the actual clustering spectrum

	double *ck2 = calloc(G->max_k+1,sizeof(double));
    
    double *ck_aim = G->ck;
    
	int i;/// the two are the same at the begining so we copy the first into the second one
	for (i=0; i<G->max_k+1; ++i){ ck2[i]    = ck[i]; }	
    
    		
 /***********************************************************************
	       we do the rewiring preserving the C(k)
 ************************************************************************/	
	
	int s1,s2,r1,r2,pos_r,pos_s,j=0;						/// rewiring variables that will store the proposed rewiring
	
	double p,AH;
	
	int accepted=0,rewirestemp=0,pirem=0;					/// during the proces we count how many proposals rewirings are with AH>0, AH<o and AH=0
	double averAH = 0,averAHneg = 0,averAHpos = 0,oldacc=0;
	int numAH0 = 0,numAHneg = 0,numAHpos = 0;
	double H =0;
	
	for(i=0;i<G->max_k +1;++i){								/// we calculate the initial energy
		
		H = H +fabs( ck[i] - ck_aim[i]);
	}
	
	
	/******** we start the rewiring *************/
	
	printf("Fixing C(k) by an annealed rewiring ...\n");fflush(stdout);
	
	time_t start,end;										/// we will measure the time
	double dif;
	time (&start);
	
		
	for(i=1; oldacc>accmin || i<rewires*G->E+2 ;++i){
		
			
		while(!choose_2_edges_random(G,&pos_r,&pos_s,randgsl)){} /// we try to find two edges avoiding selfedges and multipledges
		
		r1 = G->edge[pos_r].s;								/// the nodes are
		r2 = G->edge[pos_r].d; 
		
		s1 = G->edge[pos_s].s;
		s2 = G->edge[pos_s].d;
		
		AH = calc_AH(G,s1,s2,r1,r2,ck,ck2,ck_aim,H);	/// we calculate the increment of energy that would cause the rewiring
		
		averAH = averAH + fabs(AH);							///we also counbt the average AH of the proposals
		
		if(AH < 0.) {										/// we count how many proposals have AH > 0
			
			numAHneg++;
			averAHneg = averAHneg + AH;
		}
		else if (AH > 0.) {									/// we count how many proposals have AH < 0
			
			numAHpos++;
			averAHpos = averAHpos + AH;
			
		}
		else numAH0++;										/// we count how many proposals have AH = 0
		
		p =  gsl_rng_uniform(randgsl);						/// we throw a random number (0,1)
		
		/********** IF we acccept **************/
		
		if( p < exp(-B*AH) ){						
			
			
			swap_edges(G,s1,s2,r1,r2);					/// we make the proposed rewired
			
			G->edge[pos_r].d = s2;							/// we modify the edge vector
			G->edge[pos_s].d = r2;
			
			for( j=0; j<G->max_k+1; ++j){ck[j] = ck2[j];}		/// we update the clustering vector
			
			if(fabs(AH)>0.)	accepted++;						/// we coubt how many changes we accept
			
			H = H + AH;
		}
		
		/********** IF we reject **************/
		
		else {								
			
			for( j=0; j<G->max_k+1; ++j){ck2[j] = ck[j];}		/// we recover the old clustering vector
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
	
	
	time (&end);											///we count the rewiring time and take conclusions
	dif = difftime (end,start);
	printf ("You rewired the entire network %.2f times with %.2lf seconds.\n",(double)i/G->E, dif );
	
 /***********************************************************************
		 we free the memory
 ************************************************************************/
	
	free(ck);
	free(ck2);
	free(ck_aim);
	
	
	
	return 0;
	
}

//**********************************************************************
//**********************************************************************
/// a function that calculates the difference in the clustering that the proposed change does
double calc_AH(GRAPH* G,int s1,int s2,int r1,int r2,double *ck,double *ck2,double *ck_aim,double E0){
	
	int afectat,kaf; 			/// that variable stores the name of a common neighbour that also will change its cluster
	
	int ks1 = G->node[s1].k;	/// we look at the degrees of the afected nodes
	int ks2 = G->node[s2].k; 
	int kr1 = G->node[r1].k; 
	int kr2 = G->node[r2].k;
	
	int i,j;
	
	/******* we first discount the contribution of the links destroyed to the clustering ******/
	for(i=0; i<ks1; ++i){			
		for(j=0; j<ks2; ++j){
			
			if((G->node[s1].out[i] == G->node[s2].out[j])){
				
				afectat = G->node[s1].out[i];
				kaf     = G->node[afectat].k;
                
				ck2[ks1] = ck2[ks1] - 2./(ks1*(ks1-1.))/G->pk[ks1];
				ck2[ks2] = ck2[ks2] - 2./(ks2*(ks2-1.))/G->pk[ks2];
				ck2[G->node[afectat].k] = ck2[G->node[afectat].k] - 2./(kaf*(kaf-1.))/G->pk[kaf];
                
			}
		}
	}
	
		
	for(i=0; i<kr1; ++i){
		for(j=0; j<kr2; ++j){
			
			if(G->node[r1].out[i] == G->node[r2].out[j]){
				
				afectat	= G->node[r1].out[i];
                kaf     = G->node[afectat].k;
				
				ck2[kr1] = ck2[kr1] - 2./(kr1*(kr1-1.))/G->pk[kr1];
				ck2[kr2] = ck2[kr2] - 2./(kr2*(kr2-1.))/G->pk[kr2];
				ck2[G->node[afectat].k] = ck2[G->node[afectat].k] - 2./( kaf*( kaf-1.))/G->pk[kaf];
            
            }
        }
	}
	
	
	/*** and we add to the clustering the new common neighbours for the created links ******/
	
	for(i=0; i<ks1; ++i){			
		for(j=0; j<kr2; ++j){
			
			if(G->node[s1].out[i] == G->node[r2].out[j] && G->node[s1].out[i]!=s2 && G->node[s1].out[i]!=r1 ) {	/// we have to remember that the nodes s1 s2 and r1 r2 are no more neighbours!!
				
				afectat = G->node[s1].out[i];
				kaf     = G->node[afectat].k;
                
				ck2[ks1] = ck2[ks1] + 2./(ks1*(ks1-1.))/G->pk[ks1];
				ck2[kr2] = ck2[kr2] + 2./(kr2*(kr2-1.))/G->pk[kr2];
				ck2[G->node[afectat].k] = ck2[G->node[afectat].k] + 2./(kaf*(kaf-1.))/G->pk[kaf];
				
			}
		}
	}
	
	for(i=0; i<kr1; ++i){
		for(j=0; j<ks2; ++j){
			
			if(G->node[r1].out[i] == G->node[s2].out[j] && G->node[r1].out[i]!=r2 && G->node[r1].out[i]!=s1){
				
				afectat = G->node[r1].out[i];
                kaf     = G->node[afectat].k;
				
				ck2[kr1] = ck2[kr1] + 2./(kr1*(kr1-1.))/G->pk[kr1];
				ck2[ks2] = ck2[ks2] + 2./(ks2*(ks2-1.))/G->pk[ks2];
				ck2[G->node[afectat].k] = ck2[G->node[afectat].k] + 2./(kaf*(kaf-1.))/G->pk[kaf];
				                
			}
		}
	}
	
	/*** we calculate the energy increment ******/
	
	double AH,Enew = 0;
	
    for(i=0;i<G->max_k+1;++i){Enew = Enew + fabs(ck2[i] - ck_aim[i]);}
	
    AH = Enew - E0;
	
	return AH;
}

