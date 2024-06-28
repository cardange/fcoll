function F=system_F(x,f,tj,A,Bj,Qj,h,j,alpha)
F=x-feval(f,tj,h(j)^(alpha)*A*x+Bj+Qj);
end
