% NEW APPROACH TO THE MONTE CARLO ALGORITHM
% HOMOGENEOUS PREDATOR-PREY 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm
% and computes the error between the numerical solution and the
% simulations.
clc
clear all
close all

% parameters
Nrange = [100,1000,5000,7000,10000];

S = 1;

T = 1000;

% transition rates 
mu = 0.5;
b = 0.8; % preys birth rate 
d1 = 0.3; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.25; % competition rate (+1 predator -1 prey)
p2 = 0.05; % competition rate (-1 prey) 


% scaled parameters 
bt = mu*b;
d1t = (1-mu)*d1;
d2t = (1-mu)*d2;
p1t = mu*p1;
p2t = mu*p2;


% % time interval 
tf = linspace(0,T,T);

cc = 0;
for N=Nrange
    cc=cc+1;
    r = 2*bt-d2t;
    K = 1-d2t/(2*bt);


    % mean field equations 
    f = @(t,x)[x(1)*(2*p1t*x(2)-d1t);...
        x(2)*(r*(1-(x(2))/K)-2*(p1t+p2t+bt)*x(1))];


    n0 = round(N/4); % initial predators population 
    m0 = round(N/2); % initial preys population 

    y0 = [n0/N;m0/N]; % initial densities 

    % mean field solutions 
    [tout,yout] = ode45(f,tf,y0);
  
    for s=1:S
        M = zeros(1,N);
        I = randperm(N);
        pn = n0; % initial predators population 
        pm = m0; % initial preys population 
        M(I(1:pn)) = 1;
        M(I(pn+1:pn+pm)) = 2;

        dens1N(1,s) = pn/N; % initial predators density
        dens1M(1,s) = pm/N; % initial preys density
        I1 = randperm(N);
        J = randperm(N);

        tic
        for i=1:T-1  
            Mmu1 = M(I1(1:round(N*mu)));
            Mmu2 = M(J(1:round(N*mu)));
            M1mmu = M(I1(round(N*mu)+1:N));
            pBE = sum(Mmu1==2 & Mmu2==0);
            pEB = sum(Mmu1==0 & Mmu2==2);
            pA = sum(M1mmu==1);     
            pB = sum(M1mmu==2);
            pAB = sum(Mmu1==1 & Mmu2==2);
            pAB1 = sum(Mmu1==2 & Mmu2==1);
            rA=rand(1,pA);
            rAB=rand(1,pAB);
            rAB1=rand(1,pAB1);
            rBE=rand(1,pBE);
            rEB=rand(1,pEB);
            rB=rand(1,pB);

            pn=pn-sum(rA<d1)+sum(rAB<p1)+sum(rAB1<p1);
            pm=pm+sum(rBE<b)+sum(rEB<b)-(sum(rAB<p1)+sum(rAB1<p1)+...
                sum(rAB<p2)+sum(rAB1<p2)+sum(rB<d2));

            N1=max(N,pn+pm);

             I = randperm(N1);
            M = zeros(1,N1);
            M(I(1:pn)) = 1;
            M(I(pn+1:pn+pm)) = 2;

            dens1N(i+1,s) = pn/N;
            dens1M(i+1,s) = pm/N;
        end
           timeF(s,cc) =  toc
        end
        densN=mean(dens1N,2);
        youtN = yout(:,1);
        densM = mean(dens1M,2);
        youtM = yout(:,2);
        errN(cc) = norm(youtN-densN,inf);
        errM(cc) = norm(youtM-densM,inf);
end

 