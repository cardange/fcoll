function[b]=vector_b(alpha,Alagr)
m=size(Alagr,1);
interv=[0:m-1]';
w=gamma(interv+1)./gamma(interv+1+alpha);
b=Alagr*((1.^((alpha+interv))).*w);
end