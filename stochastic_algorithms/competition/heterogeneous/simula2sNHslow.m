% CLASSICAL MONTE CARLO ALGORITHM
% HETEROGENEOUS 2 SPECIES 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.

clc
clear all
close all 

% parameters
global N mx m1 m2 b1 b2 c11 c22 c12 c21 d1 d2 n0 q1 q2

S=1;

Nrange=[100,500,1000,5000,10000]; 

q1 = 1/3;
q2 = 1/3; 

% transition rates  
b1 = 0.5; % birth 
b2 = 0.5;
d1 = 0.5; % death
d2 = 0.5;
c11 = 0.5; % competition
c22 = 0.5;
c12 = 0;
c21 = 0;
m1 = 0.5; % migration
m2 = 0.5; 

% spatial discretization
mx = 100;
ax = 0;
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';

% time discretization
ntot=100;


cc = 0; % counter 
for N = Nrange
    cc = cc+1;
    bt2=round(N/3);
    mt2=N;
    ht1 = (bt2-at1)/(mt2-1);%dt
    n0 = N/2; % initial population in each cell ù
    m0 = N/2;
    phi0 = n0/N;
    psi0 = m0/N;
    
    % finite difference matrix without bc 
    C = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));

    % Neumann boundary conditions 
    C(1,1:2) = [-2,2]/hx^2;
    C(mx,mx-1:mx) = [2,-2]/hx^2;
    
    % piecewise continuos initial data
    y0phi = phi0* (x < mx/2)+  0* (x >=mx/2);
    y0psi = 0* (x < mx/2)+  psi0* (x >=mx/2);
    y0 = [y0phi;y0psi];

    % finite differences matrix with bc
    C1 = [m1*C,sparse(zeros(mx));sparse(zeros(mx)),m2*C];


    % cross diffusion term 
    C2 = @(y)[m1*y(1:mx);m2*y(mx+1:2*mx)]'.*C1*[y(mx+1:2*mx);y(1:mx)]-...
        [m1*y(mx+1:2*mx);m2*y(1:mx)]'.*C1*[y(1:mx);y(mx+1:2*mx)];

    % reaction term 
    g = @(y) [2*b1*y(1:mx).*(1-y(1:mx)-y(mx+1:2*mx))-c11*y(1:mx).^2-...
        2*c12*y(1:mx).*y(mx+1:2*mx)-d1*y(1:mx);...
        2*b2*y(mx+1:2*mx).*(1-y(1:mx)-y(mx+1:2*mx))-c22*y(mx+1:2*mx).^2-...
        2*c21*y(1:mx).*y(mx+1:2*mx)-d2*y(mx+1:2*mx)];


    for s=1:S
        M = zeros(N,mx);
        pn = zeros(1,mx);     
        pm = zeros(1,mx); 
        
        A = 1:N;
        
        I = repmat(A',[1,mx]);
        
        for j=1:mx/2
            I(:,j) = I(randperm(N),j);
            M(I(1:n0,j),j) = 1;
            pn(j) = sum(M(:,j));
        end
        
        for j=mx/2+1:mx
            I(:,j) = I(randperm(N),j);
            M(I(1:m0,j),j) = 2;
            pm(j) = sum(M(:,j)==2);
        end

        % initial density 
        dens1n(1,:,s) = pn/N;
        dens1m(1,:,s) = pm/N;
        
        cp=cputime;
        
        ip = ntot*(q1*N+q2*N+(1-q1-q2)*N);
        dt = ntot/ip;
        i = 0;
        for i = 1:ntot % time
        
        [pn,pm] = sol2sNHslow(M,pn,pm);
        
        M = zeros(N,mx);
        
        A = 1:N;
        
        I = repmat(A',[1,mx]);
        J = repmat(A',[1,mx]); 
        
        for j=1:mx
            I(:,j) = I(randperm(N),j);
            M(I(1:pn(j),j),j) = 1;
            M(I(pn(j)+1:pn(j)+pm(j),j),j) = 2;
        end
        
        % population in each cell at time i+1, simulation s
        dens1n(i+1,:,s)=pn/N;
        dens1m(i+1,:,s)=pm/N;
        end
        for i=1:ntot
         y(:,i+1) = y(:,i)+ht1*(C1*y(:,i)+g(y(:,i)));
        end
     end
    timef(s,cc)=cputime-cp;
    densNs(:,:,cc)=mean(dens1n,3);
    densMs(:,:,cc)=mean(dens1m,3); 
end


