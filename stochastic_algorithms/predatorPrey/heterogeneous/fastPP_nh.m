% NEW APPROACH TO THE MONTE CARLO ALGORITHM
% HETEROGENEOUS PREDATOR-PREY 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm 
% and computes the error between the numerical solution and the
% simulations.

clc
clear all
close all 

% parameters
global N mx m1 m2 p1 p2 b d1 d2 n0 m0 q1 q2

Nrange = [100,250,500,750,1000];

ntot = 1000;
ce = 4;

q1 = 1/3;
q2 = 1/3; 

% transition rates  

b = 0.1; % preys birth rate 
d1 = 0.1; % predators death rate 
d2 = 0.5; % preys death rate 
p1 = 0.25; % competition rate (+1 predator -1 prey)
p2 = 0.05; % competition rate (-1 prey)

m1 = 0.5; % migration
m2 = 0.5; 


r = 2*b-d2;
K = 1-d2/(2*b);

% spatial discretization
mx = 100;
ax = 0;
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';

cc = 0; % counter 
 for N=Nrange
    cc = cc+1;
    ht1 = q1;
    n0 = round(N/4); % predatori
    m0 = round(N/2); % prede
    phi0 = n0/N;
    psi0 = m0/N;
    
    % finite difference matrix without bc 
    C = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));

    % Neumann boundary conditions 
    C(1,1:2) = [-2,1]/hx^2;
    C(1,mx)=1/hx^2;
    C(mx,mx-1:mx) = [1,-2]/hx^2;
    C(mx,1)=1/hx^2;
    % piecewise continuos initial data
    y0phi=zeros(mx,1);
    y0psi=zeros(mx,1);
    y0phi(mx/2-ce:mx/2+ce+1)=n0/N;
    y0psi(mx/2-ce:mx/2+ce+1)=m0/N;
    y0=[y0phi;y0psi];

    % finite differences matrix with bc
    C1 = [m1*C,sparse(zeros(mx));sparse(zeros(mx)),m2*C];
    C11=[C,sparse(zeros(mx));sparse(zeros(mx)),C];

    % cross diffusion term 
    C2 = @(y)[m1*y(1:mx);m2*y(mx+1:2*mx)]'.*C11*[y(mx+1:2*mx);y(1:mx)]-...
        [m1*y(mx+1:2*mx);m2*y(1:mx)]'.*C11*[y(1:mx);y(mx+1:2*mx)];

    % reaction term 
    g = @(y) [2*p1*y(1:mx).*y(mx+1:2*mx)-d1*y(1:mx);...
        r*y(mx+1:2*mx).*(1-y(mx+1:2*mx)/K)-...
        2*(p1+p2+b)*y(1:mx).*y(mx+1:2*mx)];
    
    M1 = zeros(N,mx);
    M2 = zeros(N,mx);


    [i0,I]=sort(rand(N,ce*2+2));
    [i11,J]=sort(rand(N,ce*2+2));
    i1=N*[0:1:ce*2+1];
    I=I+i1;
    J=J+i1;
    Mt1=zeros(N,ce*2+2);
    Mt2=zeros(N,ce*2+2);
    Mt1(I(1:n0,:)) = 1;
    Mt2(J(1:n0,:)) = 1;

    Mt1(I(n0+1:n0+m0,:)) = 2;
    Mt2(J(n0+1:n0+m0,:)) = 2;

    M1(:,mx/2-ce:mx/2+ce+1)=Mt1;
    M2(:,mx/2-ce:mx/2+ce+1)=Mt2;

    pn = sum(M1==1);

    pm = sum(M1==2);

    % initial densities
    dens1n(1,:) = pn/N;
    dens1m(1,:) = pm/N;


    cp=cputime;
    y(:,1)=y0;

    for i = 1:ntot % time

        [pn,pm] = solPPfastNH(M1,M2,pn,pm);

        pnM=max(pn);
        pmM=max(pm);
        N1=max(N,pnM+pmM);
        M1 = zeros(N1,mx);
        M2 = zeros(N1,mx);


        [i0,I]=sort(rand(N1,mx));
        [i11,J]=sort(rand(N1,mx));

        for j=1:mx
            M1(I(1:pn(j),j),j)=1;
            M2(J(1:pn(j),j),j)=1;
            M1(I(pn(j)+1:pn(j)+pm(j),j),j)=2;
            M2(J(pn(j)+1:pn(j)+pm(j),j),j)=2;
        end
        % population in each cell at time i+1, simulation s
        dens1n(i+1,:)=pn/N1;
        dens1m(i+1,:)=pm/N1;

        y(:,i+1) = y(:,i)+ht1*(C1*y(:,i)+g(y(:,i))+C2(y(:,i)));

    end      
    timeF(cc)=cputime-cp
    errN(cc)=norm(dens1n(end,:)'-y(1:mx,end),inf);
    errM(cc)= norm(dens1m(end,:)'-y(mx+1:2*mx,end),inf);     
end
