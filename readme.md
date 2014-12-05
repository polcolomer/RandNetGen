RandNetGen
========================================================================

 Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________

The program takes a network from an edge list file and generates a maximally random network with the same degree sequence and any possible combination of a list of other target properties.

The program can fix:
* The original joint degree distribution P(k,k')
* The average neighbour degree Knn(k)
* The clustering coefficient
* The number of triangles
* The clustering spectrum
 
Useful for studying the effect of a certain property on any topological or dynamical property of a network or find formation patterns in real networks. 

## Requirements and Installation

  You only need to install the GSL library: http://www.gnu.org/software/gsl/
  
  How to tutorial: http://www.brianomeara.info/tutorials/brownie/gsl

## Compilation

  Execute the bash script compile.sh

    ./compile.sh


## Execute Rewire

The executable is called rewire. To introduce an argument you must introduce first its key (preceded by a dash), one white space and then the argument value. 
Arguments can appear in any order. If an argument does not appear the program gets the default value:

Examples
______________
 To randomize a network fixing its degree sequence and clustering spectrum
 
 	./rng -net netFILEname -ck original
 
 
 To get a network with the same joint degree distribution P(k,k') and a certain clustering coefficient (0.23).
 
 	./rng -net netFILEname -pkk 1 -cbar 0.23
 
 
 To get a network with the same degree sequence, with an average neighbour degree given from a file and the original number of triangles
 
 	./rng -net netFILEname -knn knnFILEname -tri original
 
 
 
## Arguments

```
-net <value>
```
Name of the input network file. It should be in the edge list format: 
A file with two columns with all the edges. With or without repetitions.
Undirected and unweighted networks only.

```
-pkk <value>
```
Rewiring method. two possible integer values: 0 or 1.
  * 0 if you only want to preserve the degree sequence
  * 1 if you want to preserve the joint degree distribution (so both, the degree sequence and the degree correlations.)
(Default:0)

```
-knn <value>:
```
Average neighbours degree (Knn(k)). Three different values:
  * "original" : the program gets the Knn(k) of the original network as the target one.
  * "filename" : give the name of a file with the target Knn(k) you want. This file should have to columns. the first one is the degree and the second the Knn(k) of nodes of such degree
  * "none"     : The program does not fix the Knn(k).
  * (Default:"none")

```
-ck <value>
```
Clustering spectrum. Three different values:
  * "original" : the program gets the clustering spectrum of the original network as the target one.
  * "filename" : give the name of a file with the target clustering spectrum you want. This file should have to columns. the first one is the degree and the second the clustering of nodes of such degree
  * "none"     : The program does not fix the clustering spectrum.
  * (Default:"none")

```
-cbar <value>
```
LOCAL clustering coefficient. Three different values:
  * "original"   : the program gets the clustering coefficient of the original network as the target one.
  * float number : the numerical value of the target clustering coefficient.
  * "none"       : The program does not fix the clustering coefficient.
  * (Default:"none")

```
-tri <value>
```
Number of triangles of the network. Three different values:
  * "original"   : the program gets the number of triangles of the original network as the target one.
  * float number : the numerical value of the target number of triangles divided by the total number of nodes.
  * "none"       : The program does not fix the number of triangles.
  * (Default:"none")


```
-rewires <value>
```
Number of rewires of each metropolis step for a given temperature. Its proportional to the total number of edges of the network E. So a value of 100 means that we do 100*E rewires each metropolis step. If the metropolis algorithm is not able to reach the aim value of the clustering try to increase this parameter.
(Default:100)

```
-beta0 <value>
```
Initial beta value of the metropolis algorithm. beta is the inverse of the temperature.
The optimal value depends on the size of the system.
I recommend to put a number that gives an initial acceptance rate above 90%.
(Default: 100)

```
-Abeta <value>
```
The incremental factor of B each time we reduce the temperature during the metropolis algorithm.
Bnew = Bold*Abeta.
If the metropolis algorithm is not able to reach the aim value of the clustering try to reduce this parameter.
(Default: 1.4)

```
-accMIN <value>
```
The minimum acceptance rate.
This parameter controls when the metropolis algorithm stops.
(Default: 0.00005)

```
-seed <value>
```
The random seed
(Default: Time of the CPU)

## Time performance

Fixing the original clustering spectrum of the PGP network that has ~24000 edges takes approx 3.5 min

Fixing the original clustering spectrum of the PGP network preserving the joint degree distribution P(k,k') approx 13 min

    

## Brief description of Rewire

This program randomizes a undirected and unweighted network by rewiring all the links but preserving the degree sequence using a similar approach to the algorithm in ref. [1]. We use two different rewiring schemes. In the first one (-pkk 0), two different edges are chosen at random. Let these connect nodes A with B and C with D. Then, the two edges are swapped so that nodes A and D, on the one hand, and C and B, on the other, are now connected. We take care that no self-connections or multiple connections between the same pair of nodes are induced by this process. This rewiring scheme preserves the degree distribution of the original network but not degree-degree correlations. In the second rewiring scheme (-pkk 1), we first chose an edge at random and look at the degree of one of its attached nodes, k. Then, a second link attached to a node of the same degree k is chosen and the two links are swapped as before. Notice that this procedure preserves both the degree of each node and the actual nodes’ degrees at the end of the two original edges. Therefore, the procedure preserves the full degree-degree correlation structure encoded in the joint distribution P(k, k'). Both procedures are ergodic and satisfy detailed balance.

In case the program is configured to fix any kind of clustering (spectrum,coefficient or number of triangles) or average neighbours degree (Knn(k)) it generates maximally random clustered networks by means of a biased rewiring procedure. Regardless of the rewiring scheme, explained above, at use, the process is biased so that generated graphs belong to an exponential ensemble of graphs. Here we consider ensembles where the Hamiltonian depends on the target property.

In case you are interested in rewiring weighted networks take a look at [multi edge randomizer](https://github.com/osagarra/Multi_edge_randomizer) by Oleguer Sagarra.


## References 

[1] Pol Colomer-de-Simón, M Angeles Serrano, Mariano Beiró, 
    J. Ignacio Alvarez-Hamelin, and    Marián Boguñá,
    “Deciphering the global organization of clustering in real
    complex networks.” Scientific reports 3, 2517 (2013).

[2] Pol Colomer-de-Simón, Marián Boguñá,
    "Double Percolation Phase Transition in Clustered Complex Networks"
    Physical Review X, 4, 041020 (2014)

## License

Copyright 2014 Pol Colomer de Simon.
All rights reserved. 
Code under License GPLv3.




