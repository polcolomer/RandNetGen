/********************************************************************************
 * 
 *					         Random Network Generator
 *
 *  The program takes a network from and edge list file and generates a maximally
 *  random network with the same degree sequence and any possible combination of 
 *  list of other target properties.
 *
 *  The program can fix:
 *  1. The original joint degree distribution P(k,k')
 *  2. The average neighbour degree Knn(k)
 *  3. The clustering coefficient
 *  4. The number of triangles
 *  5. The clustering spectrum
 *
 * 
 * You need to have the GSL library installed
 *        
 * Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.
 * 
 * ******************************************************************************/

#include "main.h"


int main(int argc, char * const argv[]){

 /***********************************************************************
            Read arguments and initialize random generator	
 ************************************************************************/ 
	
	int seed,rewires,pkk;
    double beta0,Abeta,accMIN,dk;
	char netNAME [200],ck [200],cbar [200],tri [200],knn[200];
	
	read_arguments(netNAME,&rewires,&beta0,&Abeta,&accMIN,&pkk,&dk,ck,cbar,tri,knn,&seed,argc,argv);
	
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
	
 /**********************************************************************
                       we read the network 
 ***********************************************************************/ 	
	
	GRAPH G = read_network(netNAME); /// we read the network file. 
    
	create_edges(&G);                /// we create and edge list
    
    degree_distribution(&G);          /// we calculate how many nodes have a certain degree
    
 /**********************************************************************
                 we calculate the target properties
 ***********************************************************************/ 	

    if(!strcmp(ck  ,"original")) G.ck    = clustering_spectrum(&G);     /// IF want to fix the clustering spectrum 
    else if(strcmp(ck,"none"  )) G.ck    = read_CK_fromFILE(&G,ck);        
    
    if(!strcmp(cbar,"original")) G.Ccoef = clustering_coeff(&G);         /// IF want to fix the clustering coefficient
    else if(strcmp(cbar,"none")) G.Ccoef = atof(cbar);
    
    if(!strcmp(tri ,"original")) G.triangles = numOFtrianglesXnode(&G);  /// IF want to fix the number of triangles
    else if(strcmp(tri,"none" )) G.triangles = atof(tri);
    
    if(!strcmp(knn ,"original")) G.Knn = Knn(&G);                        /// IF want to fix the average neighbour degree
    else if(strcmp(knn,"none" )) G.Knn = read_Knn_fromFILE(&G,knn);
    
 
 /**********************************************************************
					we rewire the network	
 ***********************************************************************/
	
    char nom[200],output [200];                  /// output file name
    
    /****** we do NOT preserve the degree correlations ******/
    if(!pkk && !strcmp(knn,"none")){ 
        
        rewiring_Pk(&G,rewires,randgsl);
        
        sprintf(nom,"pk_%s" ,netNAME);               

        if(strcmp(ck  ,"none")){ rewiring_Ck_annealing  (&G,beta0,Abeta,accMIN,rewires,randgsl); sprintf(output,"ck%s"  ,nom);}
                                                                                                                     
        if(strcmp(cbar,"none")){ rewiring_Cbar_annealing(&G,beta0,Abeta,accMIN,rewires,randgsl); sprintf(output,"cbar%s",nom);}
                                                                                                                     
        if(strcmp(tri ,"none")){ rewiring_TRI_annealing (&G,beta0,Abeta,accMIN,rewires,randgsl); sprintf(output,"tri%s" ,nom) ;}
        
                
    }
        
    /****** IF we preserve the degree correlations ******/
    if( pkk && !strcmp(knn,"none") ){ 
        
        rewiring_Pkk(&G,rewires,randgsl);
        
        sprintf(nom,"pkk_%s" ,netNAME);               
       
        if(strcmp(ck  ,"none")){ rewiring_PkkCk_annealing(&G,beta0,Abeta,accMIN,rewires,randgsl)  ;sprintf(output,"ck%s"  ,nom);}
                                                                                                                       
        if(strcmp(cbar,"none")){ rewiring_PkkCbar_annealing(&G,beta0,Abeta,accMIN,rewires,randgsl);sprintf(output,"cbar%s",nom);}
                                                                                                                       
        if(strcmp(tri ,"none")){ rewiring_PkkTRI_annealing(&G,beta0,Abeta,accMIN,rewires,randgsl) ;sprintf(output,"tri%s" ,nom);}
        
        
    }
    
    /****** IF we want to fix the average neighbour degree ******/
    if(strcmp(knn,"none")){
        
        rewiring_Pk(&G,rewires,randgsl);
        
        rewiring_Knn_annealing (&G,beta0,Abeta,accMIN,rewires,randgsl); 
        
        sprintf(nom,"knn_%s" ,netNAME);
        
        if(strcmp(ck  ,"none")){ rewiring_PkkCk_annealing  (&G,beta0,Abeta,accMIN,rewires,randgsl); sprintf(output,"ck%s"  ,nom);}
                                                                                                                        
        if(strcmp(cbar,"none")){ rewiring_PkkCbar_annealing(&G,beta0,Abeta,accMIN,rewires,randgsl); sprintf(output,"cbar%s",nom);}
                                                                                                                        
        if(strcmp(tri ,"none")){ rewiring_PkkTRI_annealing (&G,beta0,Abeta,accMIN,rewires,randgsl); sprintf(output,"tri%s" ,nom) ;}
        
    }
 
    /*** we make sure the output name is the proper one ****/
    if( !strcmp(ck  ,"none") && !strcmp(cbar  ,"none") && !strcmp(tri  ,"none")) sprintf(output,"%s" ,nom);
    
    if(dk>0) sprintf(output,"dk%2.1f_%s",dk,netNAME);
   
 /**********************************************************************
					we print the results	
 ***********************************************************************/
    
    print_network (&G,output,0);
    
    free(G.edge);
    free(G.pk);
	free_graph(&G);	
    gsl_rng_free(randgsl);
	
    return 0;
	
}
//**********************************************************************
//**********************************************************************
//							FUNCTIONS
//**********************************************************************
//**********************************************************************
/// a function that read the arguments of the program
int read_arguments(char netNAME[],int*rewires,double* beta0,double* Abeta,double* accMIN,int* pkk,double* dk,char ck[],char cbar[],char tri[],char knn[],int* seed,int argc,char * const  argv[]){
	
  /******  Default values *******/
	
    *pkk     = 0;
    *rewires = 100;
    *beta0   = 100;
    *Abeta   = 1.4;
    *accMIN  = 0.00005;
    *seed    = time(NULL);
    sprintf(ck     ,"none");
    sprintf(cbar   ,"none");
    sprintf(tri    ,"none");
    sprintf(knn    ,"none");
    sprintf(netNAME,"none");
    
    *dk = 0;
    
    	
  /******  Reading arguments *******/

	int i;
	for (i=1; i<argc; ++i){
                    
		if      (!strcmp(argv[i],"-rewires")) *rewires = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-pkk"    )) *pkk     = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-seed"   )) *seed    = atoi(argv[++i]);
 		else if (!strcmp(argv[i],"-beta0"  )) *beta0   = atof(argv[++i]);
		else if (!strcmp(argv[i],"-Abeta"  )) *Abeta   = atof(argv[++i]);
		else if (!strcmp(argv[i],"-accMIN" )) *accMIN  = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dk"     )) *dk      = atof(argv[++i]);
        else if (!strcmp(argv[i],"-ck"     )) sprintf(ck     ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-cbar"   )) sprintf(cbar   ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-tri"    )) sprintf(tri    ,"%s",argv[++i]);
        else if (!strcmp(argv[i],"-knn"    )) sprintf(knn    ,"%s",argv[++i]);
       	else if (!strcmp(argv[i],"-net"    )) sprintf(netNAME,"%s",argv[++i]);
		
		else {
		  fprintf(stderr,"\nError:    '%s' not recognized as an argument type.\n",argv[i]);
		  fprintf(stderr,"          Before any argument you should first write the argument -key and then the value\n"); 
		  fprintf(stderr,"Example:  ./rng -net netFILE.net -cbar 0.25\n\n");
          fprintf(stderr,"        For more information read the README file \n\n");
		  exit(-1);
		}
       
	}

    /************** dk series *****************************/
    /// if the user uses the dk sries we have to rearrange the arguments
    if( *dk > 0.9 && *dk < 1.1 )  { /// dk = 1
    
        if(*pkk != 0) printf("-pkk 1 contradicts dk 1\n");
        if(strcmp(cbar,"none")) printf("-cbar %s contradicts dk 1\n",cbar);
        if(strcmp(ck  ,"none")) printf("-ck %s contradicts dk 1\n",ck);
        if(strcmp(tri ,"none")) printf("-tri %s contradicts dk 1\n",tri);
        if(strcmp(knn ,"none")) printf("-knn %s contradicts dk 1\n",knn);
        
        
    }
    else if( *dk > 1.9 && *dk < 2.05 ) { /// dk = 2
        
        if(*pkk != 0) printf("-pkk 1 is redundant with dk 2\n");
        *pkk = 1;
        
        if(strcmp(cbar,"none")) printf("-cbar %s contradicts dk 2\n",cbar);
        if(strcmp(ck  ,"none")) printf("-ck %s contradicts dk 2\n",ck);
        if(strcmp(tri ,"none")) printf("-tri %s contradicts dk 2\n",tri);
        if(strcmp(knn ,"none")) printf("-knn %s contradicts dk 2\n",knn);
        
    }
    else if( *dk > 2.05 && *dk < 2.11 ){ /// dk = 2.1
    
        if(*pkk != 0) printf("-pkk 1 is redundant with dk 2.1\n");
        *pkk = 1;
        
        if(strcmp(cbar,"none")) printf("-cbar option contradicts dk 2.1\n");
        sprintf(cbar,"original");
        
        if(strcmp(ck  ,"none")) printf("-ck option contradicts dk 2.1\n");
        if(strcmp(tri ,"none")) printf("-tri option contradicts dk 2.1\n");
        if(strcmp(knn ,"none")) printf("-knn option contradicts dk 2.1\n");
        
        
    }
    else if( *dk > 2.4 && *dk < 2.6 )  { /// dk = 2.5
        
        if(*pkk != 0) printf("-pkk 1 is redundant with dk 2.5\n");
        *pkk = 1;
        
        if(strcmp(ck  ,"none")) printf("-ck %s contradicts dk 2.5\n",ck);
        
        sprintf(ck,"original");
        
        
        if(strcmp(cbar,"none")) printf("-cbar option contradicts dk 2.5\n");
        if(strcmp(tri ,"none")) printf("-tri option contradicts dk 2.5\n");
        if(strcmp(knn ,"none")) printf("-knn option contradicts dk 2.5\n");   
        
    }
    else if(*dk==0){}
    else{                         /// dk = others
    
        printf("the dk value that you want is %2.1f. Are you sure this make sense?\n Possible values are 1,2,2.1,2.5\n Try again and good luck\n",*dk);
        abort();        
        
    }
    
    /******************************************************/
    
    
    
    
    printf("\nARGUMENTS:\n");
    printf("   Network : %s\n",netNAME);
    if  (*pkk) printf("   Preserving the P(k,k')\n");
    else      printf("   Preserving the P(k)\n");
    if(strcmp(cbar,"none")) printf("   Target cbar: %s\n",cbar);
    if(strcmp(ck  ,"none")) printf("   Target c(k): %s\n",ck);
    if(strcmp(tri ,"none")) printf("   Target triangles: %s\n",tri);
    if(strcmp(knn ,"none")) printf("   Target Knn(k): %s\n",knn);
    printf("   Num of rewires x step : %d*E\n   Initial beta : %f\n   Beta increment : %f\n   Min acceptation rate : %f\n",*rewires,*beta0,*Abeta,*accMIN);  
    printf("   Random seed %i\n\n",*seed);
    
    if(!strcmp(netNAME,"none")){
        
        fprintf(stderr,"You did not give any edgelist file name as an imput. Its is necessary\n\n");
        fprintf(stderr,"Example:  ./rng -net netFILE.net -cbar 0.25\n\n");
        fprintf(stderr,"        For more information read the README file \n\n");
        exit(-1);
    }
    
    if(strcmp(cbar,"none") && strcmp(ck,"none")){
        
        printf("You tried to fix the clustering coefficient and the clsuteirng spectrum. Redudant\n");
        printf("We keep the most restricted constrain. The clustering spectrum\n\n");
        sprintf(cbar   ,"%s","none");
    }
    
    if(strcmp(tri,"none") && strcmp(ck,"none")){
        
        printf("You tried to fix the number of triangles and the clusteirng spectrum. Redudant\n");
        printf("We keep the most restricted constrain. The clustering spectrum\n\n");
        sprintf(cbar   ,"%s","none");
    }
    
    if(strcmp(cbar,"none") && strcmp(tri,"none")){
        
        fprintf(stderr,"You tried to fix the clustering coefficient and the number of triangles.\n");
        fprintf(stderr,"Incompatible. Decide which one you prefer and run again the program.\n\n");
        exit(-1);
    }
    if(strcmp(knn,"none") && pkk){
    
        fprintf(stderr,"You tried to fix the Knn(k) preserving the P(k,k').\n");
        fprintf(stderr,"Incompatible. Decide which one you prefer and run again the program.\n\n");
        exit(-1);
    }
    
    
    
    
  
	
	return 1;
 }

//**********************************************************************
//**********************************************************************
//							END OF CODE		
//**********************************************************************
//**********************************************************************
