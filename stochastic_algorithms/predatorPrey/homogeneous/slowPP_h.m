% CLASSICAL MONTE CARLO ALGORITHM
% HOMOGENEOUS PREDATOR-PREY 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.

clc
clear all
close all

Nrange = [100,1000,5000,7000,10000];

S = 1;

% parameters 
mu = 0.5;
b = 0.8; % preys birth rate 
d1 = 0.3; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.25; % competition rate 
p2 = 0.05; % competition rate     

% scaled parameters 
bt = mu*b;
d1t = (1-mu)*d1;
d2t = (1-mu)*d2;
p1t = mu*p1;
p2t = mu*p2;


T = 1000;
t = 0:T-1;

cc = 0;
for N=Nrange
    cc=cc+1;
    
    r = 2*bt-d2t;
    K = 1-d2t/(2*bt);


    % mean field equations 
    f = @(t,x)[x(1)*(2*p1t*x(2)-d1t);...
        x(2)*(r*(1-(x(2))/K)-2*(p1t+p2t+bt)*x(1))];


    n0 = N/4; % initial predators population 
    m0 = N/4; % initial preys population 

    y0 = [n0/N;m0/N]; % initial densities 

    % mean field solutions 
    [tout,yout] = ode45(f,t,y0);


    for s=1:S
        
        M=zeros(1,N);
        I=randperm(N);
        M(I(1:n0))=1;
        M(I(n0+1:n0+m0))=2;

        dens1n(1,s)=n0/N;  
        dens1m(1,s)=m0/N;
        pn=n0;
        pm=m0;    

        cp = cputime;
        ip = T*(mu*N+(1-mu)*N);
        dt = T/ip;
        i=1;
    
        for j = 1:dt:T
            i=i+1;

            pA=(M(3)==1).*(rand<d1*(1-mu));

            pBE=(M(1)==2 & M(2)==0).*(rand<b*mu);
            pEB=(M(1)==0 & M(2)==2).*(rand<b*mu);

            pB=(M(3)==2).*(rand<d2*(1-mu));

            pAB = (M(1)==1 & M(2)==2).*(rand<p1*mu);
            pAB1 =(M(1)==2 & M(2)==1).*(rand<p1*mu);


            pBA = (M(1)==1 & M(2)==2).*(rand<p2*mu);
            pBA1 =(M(1)==2 & M(2)==1).*(rand<p2*mu);



            pn=pn-pA+pAB+pAB1;
            pm=pm+pBE+pEB-pB-pBA-pBA1-pAB-pAB1;   

            pnM=max(pn);
            pmM=max(pm);
            N1=max(N,pnM+pmM);

            M=zeros(1,N1);
            I=randperm(N1);
            M(I(1:round(pn)))=1;
            M(I(1+round(pn):round(pn+pm)))=2;

            dens1n(i,s)=pn/N;
            dens1m(i,s)=pm/N;
        
        end
        time(s,cc)=cputime-cp
    end
        densN=mean(dens1n,2);
        youtN = yout(:,1);
        densM = mean(dens1m,2);
        youtM = yout(:,2);

end

