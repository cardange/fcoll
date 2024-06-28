function z=Q(t,alpha,gam)
interv=0:ceil(alpha)-1;
z=sum(((t.^interv)./factorial(interv)).*(gam(interv+1).'),2);
end