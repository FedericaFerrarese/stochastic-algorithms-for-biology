% NEW APPROACH TO THE MONTE CARLO ALGORITHM
% HOMOGENEOUS 1 SPECIES 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.

clc
clear all
close all

% parameters 
Nrange = [100,500,1000,5000,10000]; 

S = 1;
% time 
T = 1000;
t = 0:T;

% transition rates 
mu = 0.5;
b = 0.7; % birth rate
c = 0.5; % competition rate
d = 0.3; % death rate

% scaled parameters
bt = mu*b;
ct = mu*c;
dt = (1-mu)*d;


cc = 0; % counter 
for N=Nrange
    cc=cc+1;
    % initial population
    n0 = N/2;
    % intial density
    y0 = n0/N; 

    % mean field equation
    f = @(t,phi)2*bt*phi*(1-phi)-ct*phi^2-dt*phi;

    % mean field solution
    [tout,yout] = ode45(f,t,y0);
    youtF(:,cc) = yout;
    
    for s = 1:S
        M = zeros(1,N);
        I = randperm(N);
        p = n0;  % initial population
        M(I(1:p)) = 1;
        dens1f(1,s) = sum(M)/N; % initial density
        I = randperm(N);
        J = randperm(N);
        cp=cputime;
        for i=1:T

            Mmu1 = M(I(1:round(N*mu))); 
            Mmu2 = M(J(1:round(N*mu)));
            M1mmu = M(I(round(N*mu)+1:N));
            
            pAA = sum(Mmu1==1 & Mmu2==1);
            pAE = sum(Mmu1==1 & Mmu2==0);
            pEA = sum(Mmu1==0 & Mmu2==1);
            pA = sum(M1mmu==1);  


            rAA = rand(1,pAA);
            rAE = rand(1,pAE);
            rEA = rand(1,pEA);
            rA = rand(1,pA);

            p=p-sum(rA<d)+sum(rAE<b)+sum(rEA<b)-sum(rAA<c);
      
            % density at time i+1
            dens1f(i+1,s) = p/N; 
            I1=randperm(N);
            M = zeros(1,N);
            M(I1(1:p)) = 1;
        end
        timeF(s,cc)= cputime-cp;
    end
    dens(:,cc)=mean(dens1f,2);
end
