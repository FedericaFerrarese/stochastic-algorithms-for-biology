% CLASSICAL MONTE CARLO ALGORITHM
% HOMOGENEOUS 1 SPECIES 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.

clc
clear all
close all

Nrange = [100,1000,3000,5000,7000,10000]; 

S = 1; % simulations

% paramteters
mu = 0.5;
b = 0.9; % birth 
d = 0.9; % death
c = 0.9; % competition

% scaled parameters
bt = b*mu;
ct = c*mu;
dt = d*(1-mu);

% time discretization
T = 1000;
t = 0:T-1;

cc = 0; 
for N=Nrange
    
    cc = cc+1; 
    n0=N/2;
    f=@(t,phi)2*bt*phi*(1-phi)-phi*(dt+ct*phi);
    
    [tout,yout]=ode45(f,t,n0/N);
    
    for s=1:S
        M=zeros(1,N);
        I=randperm(N);
        M(I(1:n0))=1;
        dens(1,s)=sum(M)/N;  
        p=n0;
        ip = T*(mu*N+(1-mu)*N);
        dt = T/ip;
        i=1;
        cp=cputime;
    for j=1:dt:T
        i=i+1;
        
        pAA=(M(1)==1 & M(2)==1).*(rand<c*mu);
        pAE=(M(1)==1 & M(2)==0).*(rand<b*mu);
        pEA=(M(1)==0 & M(2)==1).*(rand<b*mu);
        pA=(M(3)==1).*(rand<d*(1-mu));
        
        p=p-pAA+pAE+pEA-pA;
        
        I1=randperm(N);
        M=zeros(1,N);
        M(I1(1:p))=1;
        dens(i+1,s)=p/N;
        
    end
    timeS(s,cc)= cputime-cp
    end

end



