GBP
===

Project containing implementations of inference in discrete probabilistic graphical models 
using several variants of Generalized Belief Propagation (GBP) message passing.

===
Usage is as follows:

  GeneralizedBeliefPropagation <uaifilename> <seed>

  Runs Genearlized Belief Propagation (GBP) on the network specified
  by <uaifilename>. The network must be in .uai format.
  Computes: 
      PR - The log10 value of the partition function (probability of evidence)
           and outputs to a file named uaifilename.PR
      BEL - Computes the belief approximations for all cliques in the model
           and outputs to a file named uaifilename.BEL
 Optional Arguments:
    -iters <value>  maximum number of iterations to run message passing for.
    -alpha <value>  step size used in message passing.
    -thresh <value>  convergence threshold.
    -dbl  runs a convergent, double-loop form of message passing.
    -iters_inner <value>  maximum number of iterations to run inner loop
          of double-loop algorithm for.
    -cvm  runs GBP on a region graph (RG) constructed using the 
		   cluster variation method (cvm).
    -ijgp <value>  Runs the mini-bucket schematic, where <value>
          specifies the iBound.
    -rg <rgfilename>
      Constructs the a region graph using the outer regions specified in
      the <rgfilename>. Inner regions are filled in using cvm.
    -verbose  Runs in verbose mode.
 

===
Sample input file is grid4x4.uai. Details of the file format can be 
found at: http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php

Sample region input file is grid4x4.rg, which specifies the 9 'faces'
of a 4-by-4 grid:

9
4 0 1 4 5
4 1 2 5 6
4 2 3 6 7
4 4 5 8 9
4 5 6 9 10
4 6 7 10 11
4 8 9 12 13
4 9 10 13 14
4 10 11 14 15

The first field is the number of regions (e.g. 9). Then for 
each region one first specifies the number of variables in that 
region (e.g. 4), followed by the variables in that region (e.g. 0 1 4 5).
