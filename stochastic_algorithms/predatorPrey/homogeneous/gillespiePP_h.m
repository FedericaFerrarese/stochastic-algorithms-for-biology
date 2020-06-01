% GILLESPIE ALGORITHM
% HOMOGENEOUS PREDATOR-PREY 1 DIMENSION
% This code implements the new approach to the Monte Carlo algorithm.

clc 
clear all
close all


% parameters
Nrange = [100,1000,5000,7000,10000];

S = 1;

t_end = 1000;

% transition rates 
mu = 0.5;
b = 0.8; % preys birth rate 
d1 = 0.3; % predators death rate 
d2 = 0.; % preys death rate 
p1 = 0.9; % competition rate (+1 predator -1 prey)
p2 = 0.9; % competition rate (-1 prey) 

 
cc = 0;
for N=Nrange 
    cc = cc+1;

    t = 0;

    n0 = N/4;
    m0 = N/4;

    A = n0;
    B = m0;
    E = N-n0-m0;

    for k = 1:S
        j = 1;
        t1(1) = t;
        A1(1) = A;
        B1(1) = B;
        tic
        while t<t_end
            j = j+1;

            a(1) = 2*b*mu*B*E/N;
            a(2) = 2*p1*mu*A*B/N;
            a(3) = 2*p2*mu*A*B/N;
            a(4) = d1*A*(1-mu);
            a(5) = d2*B*(1-mu);

            a0 = sum(a);

            r = rand(1,2);

            M = 0;
            i = 1;
            for m = 1:5
                M = M+a(m);
                if M >= r(1)*a0
                    sel(i) = m;
                    i = i+1;
                end
            end
            s = sel(1);

            if s==1
                B = B+1;
                E = N-A-B;
            end
            if s==2
                B = B-1;
                A = A+1;
                E = N-A-B;
            end
            if s==3
                B = B-1;
                E = N-A-B;
            end
            if s==4
                A = A-1;
                E = N-A-B;
            end
            if s==5
                B = B-1;
                E = N-A-B;
            end

            tau = 1/(a0)*log(1/r(2));
            t = t+tau;

            t1(j) = t;
            A1(j) = A;
            B1(j) = B;
        end
        timeG(k,cc) =  toc;
    end
 end
