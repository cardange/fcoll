clear
f=@(t,y) y.^2+1/gamma(1.5)*t.^0.5-t.^2;
b=1;
gam=0;
alpha=1/2;
%eta=[(3-sqrt(3))/6; 1-((3-sqrt(3))/6)];
%r=2.6;
N=b/2^-5;
 eta=[1/2,1]'; 
 r=4;
[t,y]=fcoll(f,b,gam,alpha,eta,r,N);
sol=@(t) t;
err=abs(y(end)-sol(b))

subplot(2,1,1)
plot(t,sol(t),'b',t,y,'*r','LineWidth',1.5)
legend('Exact solution','Numerical solution','Location','southeast')
set(gca,'Fontsize',12)

subplot(2,1,2)
semilogy(t,abs(y-sol(t)),'LineWidth',1.5)
legend('Absolute error','Location','southeast')
set(gca,'Fontsize',12)

%saveas(gcf,'figura1.fig','fig')


