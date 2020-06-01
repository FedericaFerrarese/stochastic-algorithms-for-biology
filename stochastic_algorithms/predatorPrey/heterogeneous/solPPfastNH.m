% SIMULATIONS
% INPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% pn = vector of population A in each cell at time i 
% pm = vector of population B in each cell at time i
% OUTPUT
% pn = vector of population A in each cell at time i+1
% pm = vector of population B in each cell at time i+1

function [pn,pm] = solPPfastNH(M1,M2,pn,pm)
global N mx m1 m2 p1 p2 b d1 d2 q1 q2
% q1: birth + competition
        Mq1 = M1(1:round(N*q1),:);
        Mq1c = M2(1:round(N*q1),:);
        
        % q2: migration
        Mq2 = M1(round(N*q1)+1:round(N*q1+N*q2),:);

        % q3: death
        Mq3 = M1(round(N*q1+N*q2)+1:N,:);
       
        
        pEB = sum((Mq1==2 & Mq1c==0).*(rand(round(N*q1),mx)<b)); 
        pBE = sum((Mq1==0 & Mq1c==2).*(rand(round(N*q1),mx)<b));
        
        pAB = sum((Mq1==1 & Mq1c==2).*(rand(round(N*q1),mx)<p1));
        pAB1 = sum((Mq1==2 & Mq1c==1).*(rand(round(N*q1),mx)<p1));
        pBA = sum((Mq1==1 & Mq1c==2).*(rand(round(N*q1),mx)<p2));
        pBA1 = sum((Mq1==2 & Mq1c==1).*(rand(round(N*q1),mx)<p2));
        
        pA = sum((Mq3==1).*(rand(round(N*(1-q1-q2)),mx)<d1)); 
        pB = sum((Mq3==2).*(rand(round(N*(1-q1-q2)),mx)<d2)); 
        
        pP1A =  sum((Mq2(:,1:mx-1)==1 & Mq2(:,2:mx)==0).*(rand(size(Mq2,1),mx-1)<m1)); % m1
        pP1A(mx)=0;
     
        pM1A = sum((Mq2(:,2:mx)==1 & Mq2(:,1:mx-1)==0).*(rand(size(Mq2,1),mx-1)<m1));
        pM1A(mx)=0;
        
        pP1B =  sum((Mq2(:,1:mx-1)==2 & Mq2(:,2:mx)==0).*(rand(size(Mq2,1),mx-1)<m2)); % m2
        pP1B(mx)=0;
      
        pM1B = sum((Mq2(:,2:mx)==2 & Mq2(:,1:mx-1)==0).*(rand(size(Mq2,1),mx-1)<m2));
        pM1B(mx)=0;
        
        rPcA = [0,pP1A(1:mx-1)];
        rMcA = [0,pM1A(1:mx-1)];
        rPcB = [0,pP1B(1:mx-1)];
        rMcB = [0,pM1B(1:mx-1)];
        
        pn = pn-pA-pP1A+rPcA+pM1A-rMcA+pAB+pAB1;
        pm = pm+pBE+pEB-pB-pP1B+rPcB+pM1B-rMcB-pBA-pBA1-pAB-pAB1;
        
        pn(pn<0)=0;
        pm(pm<0)=0;
        pn(pn>N)=N;
        pm(pm>N)=N;
end
        