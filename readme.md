RandNetGen
========================================================================

 Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________

The program takes a network from an edge list file and generates a maximally random network with the same degree sequence and any possible combination of a list of other target properties.

The program can fix:
* The original joint degree distribution *P(k,k')*
* The average neighbour degree *Knn(k)*
* The clustering coefficient
* The number of triangles
* The clustering spectrum

Randomize a network fixing some properties is equivalent to generate dk-random graphs [4]. Increasing values of d capture progressively more properties of the network.
dk 1 is equivalent to randomize the network fixing only the degree sequence. dk 2 fixes additionally the degree correlations. dk 2.1 fixes also the clustering coefficient and 2.5 de full clustering spectrum.

Randomize a network keeping some properties of the original network can be useful for studying the effect of a certain feature on any topological or dynamical property of a network.
Compare real networks with randomized versions of them can also reveal formation patterns.
Code developed for the publication of articles in references [1,2,3].

## Requirements and Installation

  You only need to install the GSL library: http://www.gnu.org/software/gsl/
  
  How to tutorial: http://www.brianomeara.info/tutorials/brownie/gsl

## Compilation

  Simply run:

    $ make 

  The executable created is called **RandNetGen**


## Execution

The executable is called **RandNetGen**. The program has a set of options. To introduce on option you must introduce first its option key (preceded by a dash), one white space and then the argument value. The only compulsory option is the input network file (-net).
Arguments can appear in any order. If an argument does not appear the program gets the default value:

Examples
______________
 To randomize a network fixing its degree sequence and clustering spectrum
 
 	$ ./RandNetGen -net netFILEname -ck original
 

 To randomize a network in the dk 2.1 series
 
 	$ ./RandNetGen -net netFILEname -dk 2.1
   
 
 To get a network with the same joint degree distribution *P(k,k')* and a certain clustering coefficient (0.23).
 
 	$ ./RandNetGen -net netFILEname -pkk 1 -cbar 0.23    ## fixing 
 
 
 To get a network with the same degree sequence, with an average neighbour degree given from a file and the original number of triangles
 
 	$ ./RandNetGen -net netFILEname -knn knnFILEname -tri original
 
 
 
## Options and Arguments

#### Options:
```
	-net      Network input file
	-pkk      Rewiring method
	-knn      Average neighbour degree
	-ck       Clustering spectrum
	-cbar     Local clustering coefficient
	-tri      Number of triangles per node
	-dk       Dk serie randomization
	-rewires  Number of rewires each metropolis step
	-beta0    Initial Metropolis temperature 
	-Abeta    Increment ratio of the metropolis temperature 
	-accMIN   Minimum acceptance rate
	-seed     Random seed
```
#### Arguments
**Input network file**
```
-net <value>
```
Name of the input network file. It should be in the edge list format: 
A file with two columns with all the edges. With or without repetitions.
Undirected and unweighted networks only.

**Rewiring method** 
```
-pkk <value>
```
Two possible integer values: 0 or 1.
  * 0 if you only want to preserve the degree sequence
  * 1 if you want to preserve the joint degree distribution (so both, the degree sequence and the degree correlations.)
(Default: 0)

**Average neighbours degree (*Knn(k)*)**
```
-knn <value>:
```
The Average degree of the neighbours of nodes of degree *k*, *Knn(k)**. Three different values:
  * "original" : the program gets the *Knn(k)* of the original network as the target one.
  * "filename" : give the name of a file with the target *Knn(k)* you want. This file should have to columns. the first one is the degree and the second the *Knn(k)* of nodes of such degree
  * "none"     : The program does not fix the *Knn(k)*.
  * (Default: "none")

**Clustering spectrum**
```
-ck <value>
```
The clustering coefficient of nodes of degree *k*, *C(k)*. Three different values:
  * "original" : the program gets the clustering spectrum of the original network as the target one.
  * "filename" : give the name of a file with the target clustering spectrum you want. This file should have to columns. the first one is the degree and the second the clustering of nodes of such degree
  * "none"     : The program does not fix the clustering spectrum.
  * (Default: "none")

**Clustering coefficient**
```
-cbar <value>
```
LOCAL clustering coefficient. Three different values:
  * "original"   : the program gets the clustering coefficient of the original network as the target one.
  * float number : the numerical value of the target clustering coefficient.
  * "none"       : The program does not fix the clustering coefficient.
  * (Default: "none")

**Number of triangles of the network**
```
-tri <value>
```
Three different values:
  * "original"   : the program gets the number of triangles of the original network as the target one.
  * float number : the numerical value of the target number of triangles divided by the total number of nodes.
  * "none"       : The program does not fix the number of triangles.
  * (Default: "none")

**dk series**
```
-dk <value>
```
Three different values:
  * 0   : the program does not use the dk series terminology (default value)
  * 1   : randomize only preserving the degree sequence
  * 2   : randomize preserving the joint degree distribution (degree sequence and degree correlations). Equivalent to -pkk 1
  * 2.1 : randomize preserving the joint degree distribution and the clustering coefficient. Equivalent to -pkk 1 -cbar original
  * 2.5 : randomize preserving the joint degree distribution and clustering spectrum. Equivalent to -pkk 1 -ck original
  * (Default: 0)

**Number of rewires**
```
-rewires <value>
```
Number of rewires of each metropolis step for a given temperature. Its proportional to the total number of edges of the network E. So a value of 100 means that we do 100*E rewires each metropolis step. If the metropolis algorithm is not able to reach the aim value of the clustering try to increase this parameter.
(Default: 100)

**Initial metropolis temperature**
```
-beta0 <value>
```
Beta is the inverse of the temperature. Beta0 is the initial value of beta.
The optimal value depends on the size of the system.
I recommend to put a number that gives an initial acceptance rate above 90%.
(Default: 100)

**The temperature increment**
```
-Abeta <value>
```
The incremental factor of B each time we reduce the temperature during the metropolis algorithm.
*Bnew = Bold x Abeta*.
If the metropolis algorithm is not able to reach the aim value of the clustering try to reduce this parameter.
(Default: 1.4)

**The minimum acceptance rate**
```
-accMIN <value>
```
This parameter controls when the metropolis algorithm stops.
(Default: 0.00005)

**The random seed**
```
-seed <value>
```

(Default: Time of the CPU)

## Time performance
On an intel® Core™ i7-3770 CPU @ 3.40GHz × 8

Fixing the original clustering spectrum of the PGP network that has ~24000 edges takes approx 3.5 min

Fixing the original clustering spectrum of the PGP network preserving the joint degree distribution *P(k,k')* approx 11 min

Fixing a flat clustering spectrum of *C(k)=0.25* , preserving the joint degree distribution *P(k,k')* of a scale free network of 100.000 nodes and gamma exponent 3.1 (169518 edges) takes approx 86 hours.

|Dk-Series | PGP (24000 edges)  |  US-aiports (1087 edges)|
|:----------:|:------------------:|:-----------------------:| 
| 1        |        1 sec       |          1 sec          |
| 2        |       14 sec       |          1 sec          |
| 2.1      |       11 min       |         24 sec          |
| 2.5      |       11 min       |         24 sec          |




    

## Brief description of Random Network Generator

This program randomizes a undirected and unweighted network by rewiring all the links but preserving the degree sequence using a similar approach to the algorithm in ref. [1,2]. We use two different rewiring schemes. In the first one (-pkk 0), two different edges are chosen at random. Let these connect nodes A with B and C with D. Then, the two edges are swapped so that nodes A and D, on the one hand, and C and B, on the other, are now connected. We take care that no self-connections or multiple connections between the same pair of nodes are induced by this process. This rewiring scheme preserves the degree distribution of the original network but not degree-degree correlations. In the second rewiring scheme (-pkk 1), we first chose an edge at random and look at the degree of one of its attached nodes, k. Then, a second link attached to a node of the same degree k is chosen and the two links are swapped as before. Notice that this procedure preserves both the degree of each node and the actual nodes’ degrees at the end of the two original edges. Therefore, the procedure preserves the full degree-degree correlation structure encoded in the joint distribution *P(k k')*. Both procedures are ergodic and satisfy detailed balance.

In case the program is configured to fix any kind of clustering (spectrum,coefficient or number of triangles) or average neighbours degree (*Knn(k)*) it generates maximally random clustered networks by means of a biased rewiring procedure. Regardless of the rewiring scheme, explained above, at use, the process is biased so that generated graphs belong to an exponential ensemble of graphs. Here we consider ensembles where the Hamiltonian depends on the target property.

In case you are interested in randomizing weighted networks take a look at [multi edge randomizer](https://github.com/osagarra/Multi_edge_randomizer) by Oleguer Sagarra.

## Acknowledgements

I want to thank Sergio Oller, Chiara Orsini, Oriol Vilanova and Marian Boguñá for their very useful comments and suggestions.

## References 

[1] Pol Colomer-de-Simón, M Angeles Serrano, Mariano Beiró, 
    J. Ignacio Alvarez-Hamelin, and    Marián Boguñá,
    “Deciphering the global organization of clustering in real complex networks.” 
    [Scientific reports 3, 2517 (2013)](http://www.nature.com/srep/2013/130827/srep02517/full/srep02517.html).

[2] Pol Colomer-de-Simón, Marián Boguñá,
    "Double Percolation Phase Transition in Clustered Complex Networks"
    [Physical Review X, 4, 041020 (2014)](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.4.041020)

[3] Chiara Orsini, Marija Mitrović Dankulov, Almerima Jamakovic, Priya Mahadevan, Pol Colomer-de-Simón, Amin Vahdat, Kevin E. Bassler, Zoltán Toroczkai, Marián Boguñá, Guido Caldarelli, Santo Fortunato, Dmitri Krioukov,
    "How Random are Complex Networks"
    [Arxiv.org](http://arxiv.org/abs/1505.07503)

[4] P. Mahadevan, D. Krioukov, K. Fall, and A. Vahdat,
    Systematic Topology Analysis and Generation Using Degree Correlations
    [SIGCOMM 2006](http://dl.acm.org/citation.cfm?doid=1151659.1159930) 

## License

Copyright 2014 Pol Colomer de Simon.
All rights reserved. 
Code under License GPLv3.




