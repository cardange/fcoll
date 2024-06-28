function z=sol3(t)
global alpha3 y03

for i=1:5
    sig(i)=i*alpha3;
end
sig(6)=2+alpha3;
z=y03+t.^sig(1)+t.^sig(2)+t.^sig(3)+t.^sig(4)+t.^sig(5)+t.^sig(6);
end