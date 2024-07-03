# Authors

The authors are: Angelamaria Cardone, Dajana Conte, Beatrice Paternoster (Department of Mathematics, University of Salerno, Italy).

Some information about the authors:

* Angelamaria Cardone webpage https://docenti.unisa.it/005020/home; email ancardone@unisa.it
* Dajana Conte webpage https://docenti.unisa.it/020280/home; email dajconte@unisa.it;
* Beatrice Paternoster webpage https://docenti.unisa.it/000793/home; email beapat@unisa.it.

# Description

The function fcoll.m solves an initial value problem for a fractional differential equation (FDE) by means of one step spline collocation methods.
One step spline collocation methods based are a generalization to FDEs of collocation methods. This code implements methods introduced in [3], where a suitable graded mesh is applied to obtain maximum order of convergence. Please, cite this code as [1,2] if you need.

[t,y,fval]=fcoll(f,b,gam,alpha,eta,r,N)

integrates the initial value problem for the scalar FDE of order alpha>0

     D^alpha y(t) = f(t,y(t)), t in [0,b]
     
     y^(k)(0) = gam(k+1),      k=0,...,K-1
     
where K is the smallest integer grater than alpha, and D^alpha is the fractional derivative according to the Caputo's definition. 

# Input and output arguments

INPUT
* f: function handle corresponding to the vector field of the FDE for a the scalar time variable t and the state variable y
* b: the right end of the interval [0,b]
* gam: set of initial conditions (vector of lenght K)
* alpha: order of the fractional derivative
* eta: vector of the collocation parameters. Constraint: 0 ≤ eta(1) ≤ · · · ≤eta(m).
* r>=1: grading exponent. For r=1 the mesh is uniform. For the optimal setting of r, compare [1,3].
* N: the number of mesh points.

OUTPUT
* t: vector of the graded mesh: t(i+1)=b(i/N)^r, i=0,...,N.
* y: vector of the approximate solution at the time steps t.
* fval: number of function evaluations. 

# References

[1] A. Cardone, D. Conte, B. Paternoster. (2021) A MATLAB Implementation of Spline Collocation Methods for Fractional Differential Equations. Lect. Notes Comput. Sci., vol 12949, 387-401. https://doi.org/10.1007/978-3-030-86653-2_29

[2] A. Cardone, D. Conte, Stability analysis of spline collocation methods for fractional differential equations, Math. Comput. Simulation 178 (2020), 501-514. https://doi.org/10.1016/j.matcom.2020.07.004

[3] A. Pedas, E. Tamme: Numerical solution of nonlinear fractional differential equations by spline collocation methods. J. Comput. Appl. Math. 255, 216–230 (2014). https://doi.org/10.1016/j.cam.2013.04.049
