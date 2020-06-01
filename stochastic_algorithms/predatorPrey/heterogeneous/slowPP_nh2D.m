% CLASSICAL MONTE CARLO ALGORITHM
% HETEROGENEOUS PREDATOR-PREY 2 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm 
clc
clear all
close all 

% parameters
global N mx my m1 m2 p1 p2 b d1 d2 n0 m0 q1 q2

Nrange=[100,250,500,750,1000];

ntot = 100;
ce=2;

S=1;
q1 = 1/3;
q2 = 1/3; 

% transition rates  

b = 0.1; % preys birth rate 
d1 = 0.1; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.25; % competition rate (+1 predator -1 prey)
p2 = 0.05; % competition rate (-1 prey)

m1 = 0.5; % migration
m2 =0.5; 


r = 2*b-d2;
K = 1-d2/(2*b);
% spatial discretization
mx = 100;
ax = 0;
bx = 100;
hx = (bx-ax) /(mx-1); % dx
x = linspace(ax,bx,mx);
x = x';
my = 100;
ay = 0;
by = 100;
hy = (by-ay) /(my-1); % dx
y1 = linspace(ay,by,my);
y1 = y1';
[X,Y]=ndgrid(x,y1);

cc = 0; % counter 

for N=Nrange
    cc = cc+1;
    bt2=N/3;
    mt1=N;
    ht1 = bt2/(mt1-1);%dt

    n0 = round(N/4); % predatori
    m0 = round(N/2); % prede
    phi0 = n0/N;
    psi0 = m0/N;
    
    % finite differences
    Cx = toeplitz(sparse([1,1],[1,2],[-2,1]/hx^2,1,mx));
    Cy = toeplitz(sparse([1,1],[1,2],[-2,1]/hy^2,1,mx));
    idx = eye(mx);
    idy = eye(my);

    % Dirichelet boundary condition 
    Cx(1,1:2) = [-2,2]/hx^2;
    Cx(mx,mx-1:mx) = [2,-2]/hx^2;

    Cy(1,1:2) = [-2,2]/hy^2;
    Cy(my,my-1:my) = [2,-2]/hy^2;

    C = kron(idy,Cx)+kron(Cy,idy);

    C1 = [m1*C,sparse(zeros(mx*my));...
    sparse(zeros(mx*my)),m2*C];
   
    % cross diffusion term
    C2 = @(y)[y(1:mx*my);y(mx*my+1:2*mx*my)]'*C1*[y(mx*my+1:2*mx*my);y(1:mx*my)]-...
        [y(mx*my+1:2*mx*my);y(1:mx*my)]'*C1*[y(1:mx*my);y(mx*my+1:2*mx*my)];

    % reaction term 
    g = @(y) [2*p1*y(1:mx*my).*y(mx*my+1:2*mx*my)-d1*y(1:mx*my);...
        r*y(mx*my+1:2*mx*my).*(1-y(mx*my+1:2*mx*my)/K)-...
        2*(p1+p2+b)*y(1:mx*my).*y(mx*my+1:2*mx*my)];

    % initial data
    y0phi=zeros(mx,my);
    y0psi=zeros(mx,my);
    y0phi(mx/2-ce:mx/2+ce+1,my/2-ce:my/2+ce+1)=n0/N;
    y0psi(mx/2-ce:mx/2+ce+1,my/2-ce:my/2+ce+1)=m0/N;
    y0 = [y0phi(:);y0psi(:)];
    
    M = zeros(mx,my,N);
    
    [i0,I]=sort(rand(N,mx,my));
            
    for i=mx/2-ce:mx/2+ce+1
        for j=my/2-ce:my/2+ce+1
            M(i,j,I(1:n0,i,j))=1;
            M(i,j,I(n0+1:n0+m0,i,j))=2;
            
         end
    end
     
    pn=sum(M==1,3);
    pm=sum(M==2,3);

    % initial density 
    dens1n(:,:,1) = pn/N;
    dens1m(:,:,1) = pm/N;

    cp=cputime;

    ip = ntot*(q1*N+q2*N+(1-q1-q2)*N);
    dt = ntot/ip;
    i=0;
    for k1= 1:dt:ntot % time
        i=i+1;
        [pn,pm] = solPPslowNH2D(M,pn,pm);
        
        pnM=max(max(pn));
        pmM=max(max(pm));
        N1=max(N,pnM+pmM);
        M = zeros(mx,my,N1);
        
        [i0,I]=sort(rand(N1,mx,my));
            for k=1:mx
                for j=1:my
                    M(k,j,I(1:pn(k,j),k,j))=1;
                    M(k,j,I(pn(k,j)+1:pn(k,j)+pm(k,j),k,j))=2;
            
                end
            end
        
%         % population in each cell at time i+1, simulation s
        dens1n(:,:,i+1) = pn/N;
        dens1m(:,:,i+1) = pm/N;
     end        
        timeS(cc)=cputime-cp
        dn=reshape(dens1n(:,:,end),[mx*my,1]);
        dm=reshape(dens1m(:,:,end),[mx*my,1]);

    y(:,1)=y0;
    for k2=1:ntot
        y(:,k2+1) = y(:,k2)+ht1*(C1*y(:,k2)+C2(y(:,k2))+g(y(:,k2)));
        yphi(:,:,k2+1) = reshape(y(1:mx*my,k2+1),[mx,my]);
        ypsi(:,:,k2+1) = reshape(y(mx*my+1:2*mx*my,k2+1),[mx,my]);
    end
        
end

