# sasimplex

sasimplex implements "Simplex Simulated Annealing", as described in
Numerical Recipes, by Press et al. This implementation does not use
code from Numerical Recipes. It is based on simplex2.c in version 1.16
of the Gnu Scientific Library (GSL), rewritten to use the adaptive simplex
algorithm described by Fuchang Gao and Lixing Han. 2012. "Implementing
the Nelder-Mead simplex algorithm with adaptive parameters",
(Computational Optimization and Applications 51(1):259-277, 2012).

Each parameter can be constrained to lie between specified lower and upper bounds,
using the function sasimplex_set_bounds. 

There are no references to external variables, so multiple copies of
the minimizer can be started in different threads.

# DEPENDENCIES
The Gnu Scientific Library.

# FILES
sasimplex.c : source code for minimizer. 

sasimplex.h : header file for minimizer

annealsched.c : source code for annealing schedule

annealsched.h : header file for annealsched.h

Makefile : tells "make" how to compile xsasimplex.

xsasimplex.c : an example program, which deals with a bivariate
function with local minima at integers pairs such as (0,1) and
(2,3). The global minimum is at (0,0). However each parameter is
constrained within the interval [1,10]. The constrained minimum is at
[1,1].  

# COMPILING

First install the Gnu Scientific Library. Then the example program can
be compiled by typing "make xsasimplex" from the command. Without
using "make", the command line would look like this:

  gcc -o xsasimplex xsasimplex.c sasimplex.c annealsched.c -lgsl -lgslcblas

Execution generates the following output:

>  Using minimizer sasimplex.
>   0: size=0.0001 vscale=0.0007 aberr=2.0000 x=[+1.0000, +1.0000] converged
>   1: size=0.0001 vscale=0.0003 aberr=2.0014 x=[+1.0000, +1.0014] converged
>   2: size=0.0001 vscale=0.0004 aberr=2.0002 x=[+1.0001, +1.0001] converged
>   3: size=0.0001 vscale=0.0006 aberr=2.0002 x=[+1.0002, +1.0000] converged
>   4: size=0.0001 vscale=0.0007 aberr=2.0005 x=[+1.0004, +1.0000] converged
>   5: size=0.0001 vscale=0.0003 aberr=2.0002 x=[+1.0000, +1.0002] converged
>   6: size=0.0001 vscale=0.0008 aberr=2.0003 x=[+1.0002, +1.0001] converged
>   7: size=0.0001 vscale=0.0001 aberr=2.0001 x=[+1.0000, +1.0001] converged
>   8: size=0.0001 vscale=0.0017 aberr=2.0004 x=[+1.0004, +1.0000] converged
>   9: size=0.0001 vscale=0.0007 aberr=2.0004 x=[+1.0002, +1.0001] converged
>  10: size=0.0001 vscale=0.0004 aberr=2.0002 x=[+1.0002, +1.0000] converged
>  11: size=0.0001 vscale=0.0006 aberr=2.0001 x=[+1.0000, +1.0001] converged
>  12: size=0.0001 vscale=0.0000 aberr=2.0000 x=[+1.0000, +1.0000] converged
>  13: size=0.0001 vscale=0.0007 aberr=2.0001 x=[+1.0001, +1.0001] converged
>  14: size=0.0001 vscale=0.0001 aberr=2.0001 x=[+1.0001, +1.0000] converged
>  15: size=0.0001 vscale=0.0001 aberr=2.0002 x=[+1.0000, +1.0001] converged
>  16: size=0.0001 vscale=0.0006 aberr=2.0005 x=[+1.0000, +1.0005] converged
>  17: size=0.0000 vscale=0.0003 aberr=2.0001 x=[+1.0000, +1.0001] converged
>  18: size=0.0001 vscale=0.0009 aberr=2.0002 x=[+1.0001, +1.0001] converged
>  19: size=0.0001 vscale=0.0001 aberr=2.0001 x=[+1.0000, +1.0001] converged
>  Unconstrained minimum: (0.000000, 0.000000)

Each of the 20 lines is the output of a separate minimizer, started at
a random initial position.


