% CLASSICAL MONTE CARLO ALGORITHM
% HETEROGENEOUS 1 SPECIES 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.
clc
clear all
close all 

% parameters
global N mx m b c d n0 q1 q2

Nrange=[100,1000,3000,5000,7000,10000];
S = 1; % simulations

q1 = 1/3;
q2 = 1/3; 

% transition rates  
b = 0.5; % birth 
d = 0.5; % death
c = 0.5; % competition
m = 0.5; % migration

% spatial discretization
mx = 100;
ax = 0;
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';

% time discretization
ntot = 100;


cc = 0; % counter 

for N = Nrange
      bt2=N/3;
    mt1=N;
     ht1 = bt2/(mt1-1);
    cc = cc+1;
    
    n0 = N/2; % initial population in each cell 
    phi0 = n0/N; % initial density in each cell 
    % finite differences matrix 
    C = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));

    % reaction term
    g = @(y) 2*b*y.*(1-y)-c*y.^2-d*y;

    % Neumann boudary conditions
    C(1,1:2) = [-2,2]/hx^2;
    C(mx,mx-1:mx) = [2,-2]/hx^2;

    % piecewise continuous initial condition
    y0 = zeros(mx,1);
    y0(1:mx/2)=phi0;
    C1 = C*m;
    y(:,1) = y0;
    
    for s=1:S
        % initial population distribution 
        M = zeros(N,mx);
        pn = zeros(1,mx);     

        A = 1:N;

        I = repmat(A',[1,mx]);        
            for j=1:mx/2
                I(:,j) = I(randperm(N),j);
                M(I(1:n0,j),j) = 1;
                pn(j) = sum(M(:,j));
            end

            % initial density 
            dens1n(1,:,s) = pn/N;
            cp=cputime;
            
            ip = ntot*(q1*N+q2*N+(1-q1-q2)*N);
            dt = ntot/ip;
            i = 0;
            
            for k = 1: dt: ntot % time
                
                i = i+1;
              
                pn = sol1sNHslow(M,pn);

                M = zeros(N,mx);

                A = 1:N;

                I = repmat(A',[1,mx]);

                for j=1:mx
                    I(:,j) = I(randperm(N),j);      
                    M(I(1:pn(j),j),j) = 1;
                end
                % population in each cell at time i+1, simulation s
                dens1n(i+1,:,s)=pn/N;
            end
            for i=1:ntot
                % numerical solution: Euler method
                y(:,i+1) = y(:,i)+ht1*(C1*y(:,i)+g(y(:,i)));
            end
    end
        densNs(:,:,cc)=mean(dens1n,3);

    timeS(s,cc)=cputime-cp;
end


