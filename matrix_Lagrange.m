function[Alagr]=matrix_Lagrange(eta)
m=length(eta);
mu=1;
a=[1/(eta(mu)-eta(2)),-eta(2)/(eta(mu)-eta(2))];
for i=3:m
    b=[1/(eta(mu)-eta(i)),-eta(i)/(eta(mu)-eta(i))];
    a=conv(a,b);
end
Alagr(mu,:)=a;
for mu=2:m
    a=[1/(eta(mu)-eta(1)),-eta(1)/(eta(mu)-eta(1))];
    for i=[2:mu-1,mu+1:m]
        b=[1/(eta(mu)-eta(i)),-eta(i)/(eta(mu)-eta(i))];
        a=conv(a,b);
    end   
    Alagr(mu,:)=a;
end
Alagr=fliplr(Alagr);
end