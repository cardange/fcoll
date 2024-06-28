clear
 close all
problem = input('problem ');
m=input('m (choose 2 or 3) ');


% gauss=[-sqrt(3/5);0;sqrt(3/5)];

% eta=1/2+gauss/2;
% eta=[(3-sqrt(3))/6; 1-((3-sqrt(3))/6)];
% eta=[0.5,1]';
switch m
    case 2
        eta=[1/3,2/3]';
    case 3
        eta=[1/2,3/4,1]';
    otherwise
        disp('error, set m=2 or 3')
end


m=length(eta);
if problem ==1
    % test n.1 ICCSA
    b=1;
    alpha=1/2;
    gam=[0];
    v=0.5;
    f=@f1;
    sol=@sol1;
elseif problem ==2
    % test problem n.2 ICCSA
    b=1;

    alpha=0.5;
    f=@(t,y) 40320/gamma(9-alpha)*t.^(8-alpha)-3*gamma(5+alpha/2)/gamma(5-alpha/2)*t.^(4-alpha/2)+...
        9/4*gamma(alpha+1)+(3/2*t.^(alpha/2)-t.^4).^3-y.^(3/2);

    gam=[0];
    v=0.5;
    sol=@(t) t.^8-3*t.^(4+alpha/2)+9/4*t.^(alpha);

elseif  problem==3
    % test problem n.3 iccsa
    b=8;
    global alpha3 lambda3 rho3 y03
    lambda3=-3; rho3=0.8;
    % alpha3=0.15; 
    alpha3=0.3; 
    alpha=alpha3;
    gam=[2];
    y03=gam;
    %r6=1.4; %c=[0 0.5]
    v=1-alpha;
    f=@f3;
    sol=@sol3;
elseif problem ==4
    % test 4.2 Babolian Vahidi Shoja 2014
    v=0.5;
    b=1;
    alpha=2.5;
    gam=[1;1;1];

    f=@(t,y) -y-y.^2+(1+erf(sqrt(t))).*exp(t)+exp(2*t);
    sol=@(t) exp(t);
elseif problem ==5
    % test 4.4 Babolian Vahidi Shoja 2014
    b=1;
    alpha=0.4;
    gam=[1];
    v=1-alpha;
    f=@(t,y) -(1+t.^2).*y.^2+...
        (t.^(1-alpha))/((1-alpha)*gamma(1-alpha))+(1+t.^2).*(1+t).^2;
    sol=@(t) 1+t;
end


% err=zeros(7,1);
% peff=zeros(7,1);
% fval=zeros(7,1);
% tempo=zeros(7,1);

n_run=6;
err=zeros(n_run,1);
peff=zeros(n_run,1);
fval=zeros(n_run,1);
tempo=zeros(n_run,1);

r=m/(1-v)

i=1;
% N=32;
% N=8;
 N=16;

[~,~,~]=fcoll(f,b,gam,alpha,eta,r,N);
now = tic();
[t,y,fval(i)]=fcoll(f,b,gam,alpha,eta,r,N);
tempo(i) = toc(now);
err(i)=abs(y(end)-sol(b));

for i=2:n_run
    N=N*2;

    now = tic();
    [t,y,fval(i)]=fcoll(f,b,gam,alpha,eta,r,N);
    tempo(i) = toc(now);
    err(i)=abs(y(end)-sol(b));
    peff(i)=(log(err(i-1))-log(err(i)))/log(2);
end
% format short g
vett_N=[ 16, 32, 64, 128, 256, 512]';
%vett_N=[ 32, 64, 128, 256, 512, 1024]';
% vett_N=[16, 32, 64, 128]';
% vett_N=[ 16 32, 64, 128, 256]';

const=err.*vett_N.^(m);

T = table(vett_N, err, peff, fval, tempo, const);
disp(T)

figure(1)
%subplot(2,1,1)
loglog(vett_N, err,'-*','LineWidth',1.5)
hold on
loglog(vett_N, const(end)*vett_N.^-(m),'--.','LineWidth',1.5)
grid
legend('OS m=2','slope p=m','Location','northeast')
xlabel('N');ylabel('err')
title(['Problem ', num2str(problem)])
set(gca,'Fontsize',12)

figure(2)
subplot(2,1,1)
semilogy(-log10(err),fval,'-*','LineWidth',1.5)
grid
xlabel('cd');ylabel('fval');
%title('work precision diagram')
hold on
%legend('Absolute error','Location','southeast')
set(gca,'Fontsize',12)
subplot(2,1,2)
semilogy(-log10(err),tempo,'-*','LineWidth',1.5)
% title('execution time')
grid
xlabel('cd');ylabel('time (sec)');
hold on
sgtitle(['Problem ', num2str(problem), '  m=', num2str(m)])
set(gca,'Fontsize',12)

figure(3)
subplot(2,1,1)
plot(vett_N, const,'-v','LineWidth',1.5)
title(['Problem ', num2str(problem)])
ylabel('error constant')
xlabel('N')
grid
set(gca,'Fontsize',12)

hold on





