function [t,y,fval]=fcoll(f,b,gam,alpha,eta,r,N)
%%  
%  FCOLL solves an initial value problem for a fractional differential equation 
%  (FDE) by means of one step spline collocation methods.
%
%  One step spline collocation methods based are a generalization to FDEs 
%  of collocation methods. This code implements methods introduced in [3], 
%  where a suitable graded mesh is applied to obtain maximum order of convergence.
%  Please, cite this code as [1,2] if you need.
%
%   [t,y,fval]=fcoll(f,b,gam,alpha,eta,r,N)
%   integrates the initial value problem for the scalar FDE of order alpha>0
%      D^alpha y(t) = f(t,y(t)), t in [0,b]
%      y^(k)(0) = gam(k+1),      k=0,...,K-1
%   where K is the smallest integer grater than alpha, and D^alpha is the
%   fractional derivative according to the Caputo's definition. f is a
%   function handle corresponding to the vector field of the FDE for a the
%   scalar time variable t and the state variable y.
%   The set of initial conditions gam is a vector of lenght K. 
%   eta(m) is equal to the vector of the collocation parameters. 
%   Constraint: 0 ≤ eta(1) ≤ · · · ≤eta(m).
%   r>=1 is  grading exponent. For r=1 the mesh is uniform. For the optimal
%   setting of r, compare [1,3].
%   N is the number of mesh points.
%
%  IN OUTPUT
%   t: vector of the graded mesh: t(i+1)=b(i/N)^r, i=0,...,N.
%   y: vector of the approximate solution at the time steps t.
%   fval is the number of the function evaluation. 
%
%  REFERENCES
%
%  [1] A. Cardone, D. Conte, B. Paternoster. (2021) A MATLAB Implementation
%      of Spline Collocation Methods for Fractional Differential Equations. 
% 	   Lect. Notes Comput. Sci., vol 12949, 387-401.
%      https://doi.org/10.1007/978-3-030-86653-2_29
%
%  [2] A. Cardone, D. Conte, Stability analysis of spline collocation methods 
%      for fractional differential equations, Math. Comput. Simulation 178 
%      (2020), 501-514. https://doi.org/10.1016/j.matcom.2020.07.004

%  [3] A. Pedas, E. Tamme: Numerical solution of nonlinear fractional 
%      differential equations by spline collocation methods. J. Comput. Appl. Math. 
%      255, 216–230 (2014). https://doi.org/10.1016/j.cam.2013.04.049
%%

t=b*((0:N)'/N).^r;
h=diff(t);

Alagr=matrix_Lagrange(eta);
A=matrix_A(alpha,eta,Alagr);
bvect=vector_b(alpha,Alagr); 
options=optimset('TolFun',1e-14,'TolX',1e-14,'Display','off');
m=length(eta);
Z=zeros(m,N); y=zeros(N+1,1); y(1)=gam(1);
iniz=feval(f,t(1)+eta*h(1),gam(1)*ones(m,1));
fval=0;
for j=1:N
    tj=t(j)+eta*h(j);
    Bj=lag(Z,Alagr,t,h,eta,alpha,j);
    Qj=Q(tj,alpha,gam);
    [Z(:,j),~,~,OUTPUT]=fsolve(@system_F,iniz,options,f,tj,A,Bj,Qj,h,j,alpha);
    fval=fval+OUTPUT.funcCount;
    
    Wj=lag_y(Z,Alagr,t,h,alpha,j); 
    Qj=Q(t(j+1),alpha,gam);
    y(j+1)=h(j)^(alpha)*bvect'*Z(:,j)+Wj+Qj;
    
    iniz=Z(:,j);
end




end