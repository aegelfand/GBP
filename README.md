GBP
===

Project containing implementations of inference in discrete probabilistic graphical models <br>
using several variants of Generalized Belief Propagation (GBP) message passing.

===
Usage:

  <b>GeneralizedBeliefPropagation</b> <i>&lt;uaifilename&gt; &lt;seed&gt;</i>

  Runs Genearlized Belief Propagation (GBP) on the network specified
  by in <i>uaifilename</i>. <br>The network must be in .uai format (see specification below). <br>
  Computes:  <br>
  &nbsp;&nbsp;&nbsp;&nbsp;    <b>PR</b> - The log10 value of the partition function (probability of evidence)
           and outputs to a file named <i>uaifilename.PR</i> <br>
  &nbsp;&nbsp;&nbsp;&nbsp;    <b>BEL</b> - Computes the belief approximations for all cliques in the model
           and outputs to a file named <i>uaifilename.BEL</i> <br>
 Optional Arguments: <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-iters</i> <i>&lt;value&gt;</i>  maximum number of iterations to run message passing for. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-alpha</i> <i>&lt;value&gt;</i>  step size used in message passing. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-thresh</i> <i>&lt;value&gt;</i>  convergence threshold. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-dbl</i>  runs a convergent, double-loop form of message passing. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-iters_inner</i> <i>&lt;value&gt;</i>  maximum number of iterations to run inner loop 
          of double-loop algorithm for. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-cvm</i>  runs GBP on a region graph (RG) constructed using the 
		   cluster variation method (cvm). <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-ijgp</i> <i>&lt;value&gt;</i>  Runs the mini-bucket schematic, where <value>
          specifies the iBound. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-rg</i> <i>&lt;rgfilename&gt;</i> 
      Constructs the a region graph using the outer regions specified in
      the file <i>rgfilename</i>. Inner regions are filled in using cvm. <br>
&nbsp;&nbsp;&nbsp;&nbsp;    <i>-verbose</i>  Runs in verbose mode.
 
===
Input File Specification:

Sample input file is <i>grid4x4.uai</i>. Details of the uai file format can be 
found at: http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php

Sample region input file is <i>grid4x4.rg</i>, which specifies the 9 'faces'
of a 4-by-4 grid: <br>

9 <br>
4 0 1 4 5 <br>
4 1 2 5 6 <br>
4 2 3 6 7 <br>
4 4 5 8 9 <br>
4 5 6 9 10 <br>
4 6 7 10 11 <br>
4 8 9 12 13 <br>
4 9 10 13 14 <br>
4 10 11 14 15 <br>
<br>
The first field is the number of regions (e.g. 9). Then for 
each region one first specifies the number of variables in that 
region (e.g. 4), followed by the variables in that region (e.g. 0 1 4 5).

===
Build Instructions:

The source code has no external dependencies, so you should simply be able to run 
<i>make all</i> in the <i>Release</i> folder to get a clean (and functional) build!
