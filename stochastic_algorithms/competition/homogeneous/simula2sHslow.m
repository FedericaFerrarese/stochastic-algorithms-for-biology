% CLASSICAL MONTE CARLO ALGORITHM
% HOMOGENEOUS 2 SPECIES 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.

clc
clear all
close all

Nrange=[100,1000,3000,5000];

S = 1; 

% transition rates 
mu = 0.5;   
b1 = 0.5;   % birth rate population A
b2 = 0.5;   % birth rate population B
d1 = 0.5;   % death rate population A 
d2 = 0.5;   % death rate population B
c11 = 0.5;  % intraspecies competition rate A
c22 = 0.5;  % intraspecies competition rate B
c12 = 0.;  % interspecific competition rate AB
c21 = 0.;   % interspecific competition rate BA 
% scaled parameters 
b1t = mu*b1;
b2t = mu*b2;
c11t = mu*c11;
c12t = mu*c12;
c22t = mu*c22;
c21t = mu*c21;
d1t = (1-mu)*d1;
d2t = (1-mu)*d2;

T=500;
t=0:T-1;

cc=0;

for N=Nrange
    cc=cc+1;
    
    n0=N/2;
    m0=N/4;
    % initial densities
    y0 = [n0/N;m0/N];

    % mean field equations
    f = @(t,y)[2*b1t*y(1)*(1-y(1)-y(2))-(c11t*y(1)^2+d1t*y(1)+2*c12t*y(1)*y(2));...
    2*b2t*y(2)*(1-y(1)-y(2))-(c22t*y(2)^2+d2t*y(2)+2*c21t*y(1)*y(2))];

    % mean field solutions 
    [tout,yout] = ode45(f,t,y0);
    
    for s=1:S
        M = zeros(1,N);
        I = randperm(N);
        M(I(1:n0)) = 1;
        M(I(n0+1:n0+m0)) = 2;
        dens1n(1,s) = n0/N;  
        dens1m(1,s) = m0/N;
        pn = n0;
        pm = m0;
        ip = T*(mu*N+(1-mu)*N);
        dt = T/ip;
        i=1;
        cp=cputime;
        for j=1:dt:T
            i = i+1;

            pAA = (M(1)==1 & M(2)==1).*(rand<c11*mu);
            pAE = (M(1)==1 & M(2)==0).*(rand<b1*mu);
            pEA = (M(1)==0 & M(2)==1).*(rand<b1*mu);
            pA = (M(3)==1).*(rand<d1*(1-mu));

            pBB = (M(1)==2 & M(2)==2).*(rand<c22*mu);
            pBE = (M(1)==2 & M(2)==0).*(rand<b2*mu);
            pEB = (M(1)==0 & M(2)==2).*(rand<b2*mu);
            pB = (M(3)==2).*(rand<d2*(1-mu));

            pAB =(M(1)==1 & M(2)==2).*(rand<c12*mu);
            pAB1 = (M(1)==2 & M(2)==1).*(rand<c12*mu);
            pBA = (M(1)==1 & M(2)==2).*(rand<c21*mu);
            pBA1 = (M(1)==2 & M(2)==1).*(rand<c21*mu);

            pn = pn-pAA+pAE+pEA-pA-pAB-pAB1;
            pm = pm-pBB+pBE+pEB-pB-pBA-pBA1;   

            N1 = max(N,pn+pm);
            I = randperm(N1);
            M = zeros(1,N1);
            M(I(1:pn)) = 1;
            M(I(1+pn:pn+pm)) = 2;
            dens1n(i+1,s) = pn/N;
            dens1m(i+1,s) = pm/N;
    end
    time(s,cc) = cputime-cp;
    end
end

