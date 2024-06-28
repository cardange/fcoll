function Bj=lag(Z,Alagr,t,h,eta,alpha,j)
m=length(eta);
Bj=zeros(m,1);
interv=[0:m-1];
w=gamma(interv+1)./gamma(interv+1+alpha);
for lambda=1:j-1
    d=(t(j)-t(lambda)+eta*h(j))/h(lambda);
    for k=1:m
        W(k,:)=betainc(1/d(k),1+interv,alpha);
    end
    E=(((d.^((alpha+interv))).*w).*W)*Alagr';
    Bj=Bj+(h(lambda)^alpha)*E*Z(:,lambda);
end
end