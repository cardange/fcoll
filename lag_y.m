function Wj=lag_y(Z,Alagr,t,h,alpha,j)
Wj=0;
m=size(Alagr,1);
interv=[0:m-1]';
w=gamma(interv+1)./gamma(interv+1+alpha);
for lambda=1:j-1
    d=(t(j+1)-t(lambda))/h(lambda);
    W=betainc(1/d,1+interv,alpha);
    g=Alagr*(((d.^((alpha+interv))).*w).*W);
    Wj=Wj+(h(lambda)^alpha)*g'*Z(:,lambda);
end
end