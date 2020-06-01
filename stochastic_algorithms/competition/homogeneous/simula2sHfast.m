% NEW APPROACH TO THE MONTE CARLO ALGORITHM
% HOMOGENEOUS 2 SPECIES 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.
clc
clear all
close all

S = 1;
% parameters 
Nrange=[100,1000,3000,5000]; 

% time 
T = 1000;
t = 0:T;

% transition rates 
mu = 0.5;   
b1 = 0.5;   % birth rate population A
b2 = 0.5;   % birth rate population B
d1 = 0.1;   % death rate population A 
d2 = 0.2;   % death rate population B
c11 = 0.4;  % intraspecies competition rate A
c22 = 0.1;  % intraspecies competition rate B
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

cc = 0; % counter 
for N=Nrange

    cc = cc+1;
    
    % initial populations
    n0 = round(N/2);
    m0 = N/2;
    % initial densities
    y0 = [n0/N;m0/N];
    
    % mean field equations
    f = @(t,y)[2*b1t*y(1)*(1-y(1)-y(2))-(c11t*y(1)^2+d1t*y(1)+2*c12t*y(1)*y(2));...
    2*b2t*y(2)*(1-y(1)-y(2))-(c22t*y(2)^2+d2t*y(2)+2*c21t*y(1)*y(2))];
   
    % mean field solutions 
    [tout,yout] = ode45(f,t,y0);


     for s = 1:S
        pn = n0;
        pm = m0;
        % initial random distribution of the populations
        M = zeros(1,N);
        I = randperm(N);
        dens1N(1,s) = pn/N;
        dens1M(1,s) = pm/N;
        % initial populations 
        M(I(1:pn)) = 1;
        M(I(pn+1:pm+pn)) = 2;
        % initial densities 

        cp(s)=cputime;
        for i = 1:T
            I = randperm(N);
            J = randperm(N);
            Mmu1 = M(I(1:round(N*mu)));
            Mmu2 = M(J(1:round(N*mu)));

            M1mmu = M(I(round(N*mu)+1:N));

            pAA = sum(Mmu1==1 & Mmu2==1);
            pAE = sum(Mmu1==1 & Mmu2==0);
            pEA = sum(Mmu1==0 & Mmu2==1);
            pA = sum(M1mmu==1);  

            pBB = sum(Mmu1==2 & Mmu2==2);
            pBE = sum(Mmu1==2 & Mmu2==0);
            pEB = sum(Mmu1==0 & Mmu2==2);
            pB = sum(M1mmu==2); 

            pAB = sum(Mmu1==1 & Mmu2==2);
            pBA = sum(Mmu1==2 & Mmu2==1);

            rAA = rand(1,pAA);
            rAE = rand(1,pAE);
            rEA = rand(1,pEA);
            rA = rand(1,pA);
            rBB = rand(1,pBB);
            rBE = rand(1,pBE);
            rEB = rand(1,pEB);
            rB = rand(1,pB);
            rAB = rand(1,pAB);
            rBA = rand(1,pBA);

            pn=pn-sum(rA<d1)+sum(rAE<b1)+sum(rEA<b1)-sum(rAA<c11)-...
            sum(rAB<c12)-sum(rBA<c12);
            pm=pm-sum(rB<d2)+sum(rBE<b2)+sum(rEB<b2)-sum(rBB<c22)-...
            sum(rAB<c21)-sum(rBA<c21);
            
            pn(pn<0)=0;
            pm(pm<0)=0;
            
            N1=max(N,pn+pm);
            
            I1 = randperm(N1);
            M = zeros(1,N1);
            M(I1(1:pn)) = 1;
            M(I1(pn+1:pm+pn)) = 2;

            % densities at time i+1
            dens1N(i+1,s) = pn/N1;
            dens1M(i+1,s) = pm/N1;
        end
         time(s,cc) =  cputime-cp(s)
     end
     densN(:,cc)=mean(dens1N,2);
     densM(:,cc)=mean(dens1M,2);
 
end
