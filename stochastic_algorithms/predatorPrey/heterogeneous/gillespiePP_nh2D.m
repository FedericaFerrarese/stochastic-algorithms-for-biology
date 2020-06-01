% GILLESPIE ALGORITHM
% HETEROGENEOUS PREDATOR-PREY 2 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm 
clc 
clear all
close all


% parameters
ce=49;
Nrange=[100,250,500,750,1000];

mx=100;
my=100;
ax=0;
bx=100;
x=linspace(ax,bx,mx);
hx=(bx-ax)/(mx-1);
ay=0;
by=100;
y1=linspace(ay,by,my);
hy=(by-ay)/(my-1);

t_end=10;

S=1;

% transition rates 
q1=1/3;
q2=1/3;
b = 0.1; % preys birth rate 
d1 = 0.1; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.25; % competition rate (+1 predator -1 prey)
p2 = 0.05; % competition rate (-1 prey) 
m1=0.1;
m2=0.1;
cc=0;

for N=Nrange

    t=0;
    cc=cc+1;

    r1 = 2*b*q1-d2*q1;
    K = 1-d2*q1/(2*b*q1);

    n0=N/4;
    m0=N/2;
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
    g = @(y) [2*p1*q1*y(1:mx*my).*y(mx*my+1:2*mx*my)-d1*q1*y(1:mx*my);...
        r1*y(mx*my+1:2*mx*my).*(1-y(mx*my+1:2*mx*my)/K)-...
        2*q1*(p1+p2+b)*y(1:mx*my).*y(mx*my+1:2*mx*my)];

    % initial data
    y0phi=zeros(mx,my);
    y0psi=zeros(mx,my);
    y0phi(mx/2-ce:mx/2+ce+1,my/2-ce:my/2+ce+1)=n0/N;
    y0psi(mx/2-ce:mx/2+ce+1,my/2-ce:my/2+ce+1)=m0/N;
    y0 = [y0phi(:);y0psi(:)];

    y(:,1)=y0;


    A = n0;
    B = m0;
    E = N-n0-m0;


    j=0;
    t1(1)=t;
    A1(1,:,:)=zeros(1,mx,my);
    B1(1,:,:)=zeros(1,mx,my);
    A1(1,mx/2-ce:mx/2+ce+1,my/2-ce:my/2+ce+1)=A;
    B1(1,mx/2-ce:mx/2+ce+1,my/2-ce:my/2+ce+1)=B;
    E1=N-A1(1,:,:)-B1(1,:,:);
    cp=cputime;
    while t<t_end
        j=j+1;
        a(1,:,:)=2*b*q1*B1(j,:,:).*E1(j,:,:)/N;
        a(2,:,:)=2*p1*q1*A1(j,:,:).*B1(j,:,:)/N;
        a(3,:,:)=2*p2*q1*A1(j,:,:).*B1(j,:,:)/N;
        a(4,:,:)=d1*A1(j,:,:)*(1-q1-q2);
        a(5,:,:)=d2*B1(j,:,:)*(1-q1-q2);

        a(6,1:mx-1,:)=q1*(m1*A1(j,1:mx-1,:).*E1(j,2:mx,:))/N;
        a(7,2:mx,:)=q1*(m1*A1(j,2:mx,:).*E1(j,1:mx-1,:))/N;
        a(8,1:mx-1,:)=q1*(m2*B1(j,1:mx-1,:).*E1(j,2:mx,:))/N; % va in quella dopo
        a(9,2:mx,:)=q1*(m2*B1(j,2:mx,:).*E1(j,1:mx-1,:))/N; % va in quella prima

        a(10,:,1:my-1)=q1*(m1*A1(j,:,1:my-1).*E1(j,:,2:my))/N;
        a(11,:,2:my)=q1*(m1*A1(j,:,2:my).*E1(j,:,1:my-1))/N;
        a(12,:,1:my-1)=q1*(m2*B1(j,:,1:my-1).*E1(j,:,2:my))/N; % va in quella dopo
        a(13,:,2:my)=q1*(m2*B1(j,:,2:my).*E1(j,:,1:my-1))/N; % va in quella prima
        a0=sum(a);
        r=rand(1,mx+1,my+1);
        M=zeros(1,mx,my);

        for m=1:size(a,1)
            M(1,:,:)=M(1,:,:)+a(m,:,:);
            sel(m,:,:)=(M(1,:,:)>r(1,1:mx,1:my).*a0);
        end
        s=size(a,1)-sum(sel)+1;

        % righe
        s6=(s==6);
        s61=zeros(1,mx,my);
        s61(1,2:mx,:)=s6(1,1:mx-1,:);
        s7=(s==7);
        s71=zeros(1,mx,my);
        s71(1,1:mx-1,:)=s7(1,2:mx,:);


        s8=(s==8);
        s81=zeros(1,mx,my);
        s81(1,2:mx,:)=s8(1,1:mx-1,:);
        s9=(s==9);
        s91=zeros(1,mx,my);
        s91(1,1:mx-1,:)=s9(1,2:mx,:);

        % colonne
        s10=(s==10);
        s101=zeros(1,mx,my);
        s101(1,2:mx,:)=s10(1,1:mx-1,:);
        s11=(s==11);
        s111=zeros(1,mx,my);
        s111(1,1:mx-1,:)=s11(1,2:mx,:);


        s12=(s==12);
        s121=zeros(1,mx,my);
        s121(1,2:mx,:)=s12(1,1:mx-1,:);
        s13=(s==13);
        s131=zeros(1,mx,my);
        s131(1,1:mx-1,:)=s13(1,2:mx,:);

        B1(j+1,:,:)=B1(j,:,:)+(s==1)-(s==2)-(s==3)-(s==5)-s8-s9+s81+s91-s12-s13+s121+s131;

        A1(j+1,:,:)=A1(j,:,:)+(s==2)-(s==4)-s6-s7+s61+s71-s10-s11+s101+s111;

        E1(j+1,:,:)=N-A1(j+1,:,:)-B1(j+1,:,:);

        tau=1/max(max(a0))*log(1/r(1,mx+1,my+1));
        t=t+tau;

        t1(j)=t;
end
      timeG(cc)=cputime-cp
end

  ht1=0.3;
for i=1:t_end
    y(:,i+1) = y(:,i)+ht1*(C1*y(:,i)+C2(y(:,i))+g(y(:,i)));
    yphi(:,:,i+1) = reshape(y(1:mx*my,i+1),[mx,my]);
    ypsi(:,:,i+1) = reshape(y(mx*my+1:2*mx*my,i+1),[mx,my]);
end
