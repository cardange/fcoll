function[A]=matrix_A(alpha,eta,Alagr)
m=length(eta);
interv=[0:m-1];
w=gamma(interv+1)./gamma(interv+1+alpha);
A=((eta.^((alpha+interv))).*w)*Alagr';
end