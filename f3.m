function z=f3(t,y)
global alpha3 lambda3 rho3 y03

for i=1:5
    sig(i)=i*alpha3;
    c=sig(i);
    dertc{i}=(c*t.^(c - alpha3)*gamma(c))./gamma(1 + c - alpha3);
end
sig(6)=2+alpha3;
c=sig(6);
dertc{6}=(c*t.^(c - alpha3)*gamma(c))./gamma(1 + c - alpha3);
    
u=y03*t.^0+t.^sig(1)+t.^sig(2)+t.^sig(3)+t.^sig(4)+t.^sig(5)+t.^sig(6);
g=dertc{1}+dertc{2}+dertc{3}+dertc{4}+dertc{5}+dertc{6}-lambda3*u-rho3*u.*(1-u.^2);

z=lambda3*y+rho3*y.*(1-y.^2)+g;
end
