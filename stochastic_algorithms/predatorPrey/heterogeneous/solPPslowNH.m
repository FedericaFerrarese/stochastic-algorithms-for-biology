% SIMULATIONS
% INPUT
% M = vector of initial random distribution of the population 
% pn = vector of predator population in each cell time i 
% pm = vector of prey population in each cell time i 
% OUTPUT
% pn = vector of predator population in each cell time i+1
% pm = vector of prey population in each cell time i+1
function [pn,pm] = solPPslowNH(M,pn,pm)
global mx m1 m2 p1 p2 b d1 d2 q1 q2 
    % q1: birth + competition

    pA=(M(3,:)==1).*(rand(1,mx)<d1*(1-q1-q2));


    pBE=(M(1,:)==2 & M(2,:)==0).*(rand(1,mx)<b*q1);
    pEB=(M(1,:)==0 & M(2,:)==2).*(rand(1,mx)<b*q1);

    pB=(M(3,:)==2).*(rand(1,mx)<d2*(1-q1-q2));

    pAB = (M(1,:)==2 & M(2,:)==1).*(rand(1,mx)<p1*q1);
    pAB1 = (M(1,:)==1 & M(2,:)==2).*(rand(1,mx)<p1*q1);

    pBA = (M(1,:)==2 & M(2,:)==1).*(rand(1,mx)<p2*q1);
    pBA1 =(M(1,:)==1 & M(2,:)==2).*(rand(1,mx)<p2*q1);

    pP1A=(M(4,1:mx-1)==1 & M(4,2:mx)==0).*(rand(1,mx-1)<m1*q2);
    pM1A=(M(4,1:mx-1)==0 & M(4,2:mx)==1).*(rand(1,mx-1)<m1*q2);
    pP1A(mx)=0;
    pM1A(mx)=0;

    pP1B=(M(4,1:mx-1)==2 & M(4,2:mx)==0).*(rand(1,mx-1)<m2*q2);
    pM1B=(M(4,1:mx-1)==0 & M(4,2:mx)==2).*(rand(1,mx-1)<m2*q2);
    pP1B(mx)=0;
    pM1B(mx)=0;


    rPcA = zeros(1,mx);
    rPcA(2:mx)=pP1A(1:mx-1);
    rMcA = zeros(1,mx);
    rMcA(2:mx)=pM1A(1:mx-1);


    rPcB = zeros(1,mx);
    rPcB(2:mx)=pP1B(1:mx-1);
    rMcB = zeros(1,mx);
    rMcB(2:mx)=pM1B(1:mx-1);

    pn=pn-pA+pAB+pAB1-pP1A+rPcA+pM1A-rMcA;
    pm=pm+pBE+pEB-pB-pBA-pBA1-pAB-pAB1-pP1B+rPcB+pM1B-rMcB;

    pn(pn<0)=0;
    pm(pm<0)=0;
end
