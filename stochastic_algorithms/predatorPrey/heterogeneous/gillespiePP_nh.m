% GILLESPIE ALGORITHM
% HETEROGENEOUS PREDATOR-PREY 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm 
clc
clear all
close all


% parameters
ce = 49;
Nrange = [100,250,500,750,1000];
mx = 100;
t_end = 100;

S = 1;

% transition rates 
q1=1/3;
q2=1/3;
b = 0.9; % preys birth rate 
d1 = 0.9; % predators death rate 
d2 = 0.9; % preys death rate 
p1 = 0.9; % competition rate (+1 predator -1 prey)
p2 = 0.9; % competition rate (-1 prey) 
m1=0.9;
m2=0.9;
cc=0;

for N=Nrange

    cc=cc+1;
    t=0;

    r1 = 2*b-d2;
    K = 1-d2/(2*b);
    ax=0;
    bx=100;
    x=linspace(ax,bx,mx);
    hx=(bx-ax)/(mx-1);

    n0=N/4;
    m0=N/4;
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
    r1*y(mx+1:2*mx).*(1-y(mx+1:2*mx)/K)-...
    2*(p1+p2+b)*y(1:mx).*y(mx+1:2*mx)];

    A=n0;
    B=m0;
    E=N-n0-m0;


    j=0;
    t1(1)=t;
    
    A1(1,:)=zeros(1,mx);
    B1(1,:)=zeros(1,mx);
    A1(1,mx/2-ce:mx/2+ce+1)=A;
    B1(1,mx/2-ce:mx/2+ce+1)=B;
    E1=N-A1-B1;

    iter=0;
    tic
    y(:,1)=y0;

    while t<t_end 
        iter=iter+1;
        j=j+1;

        a(1,:)=2*b*q1*B1(j,:).*E1(j,:)/N;
        a(2,:)=2*p1*q1*A1(j,:).*B1(j,:)/N;
        a(3,:)=2*p2*q1*A1(j,:).*B1(j,:)/N;
        a(4,:)=d1*q1*A1(j,:);
        a(5,:)=d2*q1*B1(j,:);
        a(6,1:mx-1)=q1*(m1*A1(j,1:mx-1).*E1(j,2:mx))/N;
        a(7,2:mx)=q1*(m1*A1(j,2:mx).*E1(j,1:mx-1))/N;
        a(8,1:mx-1)=q1*(m2*B1(j,1:mx-1).*E1(j,2:mx))/N;
        a(9,2:mx)=q1*(m2*B1(j,2:mx).*E1(j,1:mx-1))/N;

        a0=sum(a);
        r=rand(1,mx+1);

        M=zeros(1,mx);

        for m=1:size(a,1)
            M(1,:)=M(1,:)+a(m,:);
            sel(m,:)=(M(1,:)>r(1:mx).*a0);
        end
        s=size(a,1)-sum(sel)+1;

        s6=(s==6);
        s61=[0,s6(1:mx-1)];
        s7=(s==7);
        s71=[s7(2:mx),0];

        s8=(s==8);
        s81=[0,s8(1:mx-1)];
        s9=(s==9);
        s91=[s9(2:mx),0];

        B1(j+1,:)=B1(j,:)+(s==1)-(s==2)-(s==3)-(s==5)-s8-s9+s81+s91;

        A1(j+1,:)=A1(j,:)+(s==2)-(s==4)-s6-s7+s61+s71;
        
        E1(j+1,:)=N-A1(j+1,:)-B1(j+1,:);

        tau=1/max(a0)*log(1/r(2));
        t=t+tau;
       
        t1(j)=t;
    end

    timeG(cc)=toc;
end
ht1=0.3;
y(:,1)=y0;

for i=1:t_end
    y(:,i+1) = y(:,i)+ht1*(C1*y(:,i)+C2(y(:,i))+g(y(:,i)));
end
