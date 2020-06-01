% CLASSICAL MONTE CARLO ALGORITHM
% HETEROGENEOUS PREDATOR-PREY 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm 
clc
clear all
close all 

% parameters
global N mx m1 m2 p1 p2 b d1 d2 n0 m0 q1 q2

Nrange=[100,250,500,750,1000];

ntot =100;

ce=3;

S=1;
q1 = 1/3;
q2 = 1/3; 

% transition rates  
b = 0.1; % preys birth rate 
d1 = 0.1; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.25; % competition rate (+1 predator -1 prey)
p2 = 0.05; % competition rate (-1 prey)

m1 = 1; % migration
m2 =1; 


r = 2*b-d2;
K = 1-d2/(2*b);
% spatial discretization
mx =100;
ax = 0;
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';

% time discretization
% mt =50;
% at1 = 0;


% number of simulations
cc = 0; % counter 
for N=Nrange
    cc = cc+1;
    bt2=N/3;
    mt1=N;
    ht1 = bt2/(mt1-1);%dt
    n0 = round(N/4); 
    m0 = round(N/2); 
    phi0 = n0/N;
    psi0 = m0/N;
    
    % finite difference matrix without bc 
    C = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));

    % Neumann boundary conditions 
    C(1,1:2) = [-2,2]/hx^2;
    C(mx,mx-1:mx) = [2,-2]/hx^2;
    
    % initial data
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


    M = zeros(N,mx);

    [i0,I]=sort(rand(N,ce*2+2));

    i1=N*[0:1:ce*2+1];
    I=I+i1;

    Mt=zeros(N,ce*2+2);

    Mt(I(1:n0,:)) = 1;

    Mt(I(n0+1:n0+m0,:)) = 2; 

    M(:,mx/2-ce:mx/2+ce+1)=Mt;


    pn = sum(M==1);

    pm = sum(M==2);

    % initial density 
    dens1n(1,:) = pn/N;
    dens1m(1,:) = pm/N;

    cp=cputime;
    y(:,1)=y0;


    ip = ntot*(q1*N+q2*N+(1-q1-q2)*N);
    dt = ntot/ip;
    i=0;
    for k= 1:dt:ntot % time
        i=i+1;

        [pn,pm] = solPPslowNH(M,pn,pm);

        pnM=max(pn);
        pmM=max(pm);
        N1=max(N,pnM+pmM);

        M= zeros(N1,mx);

        [i0,I]=sort(rand(N1,mx));
       
        for j=1:mx
            M(I(1:pn(j),j),j)=1;
            
            M(I(pn(j)+1:pn(j)+pm(j),j),j)=2;
            
        end
            
        pnM=max(pn);
        pmM=max(pm);
        N1=max(N,pnM+pmM);
        
        M= zeros(N1,mx);
        
        [i0,I]=sort(rand(N1,mx));
       
        for j=1:mx
            M(I(1:pn(j),j),j)=1;
            
            M(I(pn(j)+1:pn(j)+pm(j),j),j)=2;
            
        end
             
%         % population in each cell at time i+1, simulation s
        dens1n(i+1,:)=pn/N;
        dens1m(i+1,:)=pm/N;
         end
     
    timeS(cc)=cputime-cp;
    densNf(:,:,cc)=mean(dens1n,3);
    densMf(:,:,cc)=mean(dens1m,3);

     for i=1:ntot
        y(:,i+1) = y(:,i)+ht1*(C1*y(:,i)+g(y(:,i))+C2(y(:,i)));
    end
    

end
