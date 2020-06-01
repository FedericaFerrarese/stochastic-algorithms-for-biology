% SIMULATIONS
% INPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% pn = vector of population A in each cell at time i 
% pm = vector of population B in each cell at time i
% OUTPUT
% pn = vector of population A in each cell at time i+1
% pm = vector of population B in each cell at time i+1

function [pn,pm] = sol2sNHfast(M1,M2,pn,pm)
global N mx m1 m2 b1 b2 c11 c22 c12 c21 d1 d2 q1 q2
        % q1: birth + competition
        Mq1 = M1(1:round(N*q1),:);
        Mq1c = M2(1:round(N*q1),:);
        
        % q2: migration
        Mq2 = M1(round(N*q1)+1:round(N*q1+N*q2),:);

        % q3: death
        Mq3 = M1(round(N*q1+N*q2)+1:N,:);
        

        % POPULATION A 
        % AjAj in AjEj
        pAA = sum((Mq1==1 & Mq1c==1).*(rand(round(N*q1),mx)<c11)); % competition
        
        % AjEj in AjAj
        pAE = sum((Mq1==1 & Mq1c==0).*(rand(round(N*q1),mx)<b1)); % birth
        
        % EjAj in AjAj
        pEA = sum((Mq1==0 & Mq1c==1).*(rand(round(N*q1),mx)<b1)); % birth
    
        % Aj in Ej
        pA = sum((Mq3==1).*(rand(round(N*(1-q1-q2)),mx)<d1)); % death
       
        % from cell j to j+1

        % EiAj in AiEj
        pP1A =  sum((Mq2(:,1:mx-1)==1 & Mq2(:,2:mx)==0).*(rand(size(Mq2,1),mx-1)<m1)); % migration
        pP1A(mx)=0;
      
        % from cell j to j-1

        % EiAj in AiEj
        pM1A = sum((Mq2(:,2:mx)==1 & Mq2(:,1:mx-1)==0).*(rand(size(Mq2,1),mx-1)<m1)); % migration
        pM1A(mx)=0;
        
           % POPULATION B 
        % BjBj in BjEj
        pBB = sum((Mq1==2 & Mq1c==2).*(rand(round(N*q1),mx)<c22)); % competition
        
        % AjEj in AjAj
        pBE = sum((Mq1==2 & Mq1c==0).*(rand(round(N*q1),mx)<b2)); % birth
        
        % EjAj in AjAj
        pEB = sum((Mq1==0 & Mq1c==2).*(rand(round(N*q1),mx)<b2)); % birth
    
        % Aj in Ej
        pB = sum((Mq3==2).*(rand(round(N*(1-q1-q2)),mx)<d2)); % death
       
        % from cell j to j+1

        % EiBj in BiEj
        pP1B =  sum((Mq2(:,1:mx-1)==2 & Mq2(:,2:mx)==0).*(rand(size(Mq2,1),mx-1)<m2)); % migration
        pP1B(mx)=0;
      
        % from cell j to j-1

        % EiAj in AiEj
        pM1B = sum((Mq2(:,2:mx)==2 & Mq2(:,1:mx-1)==0).*(rand(size(Mq2,1),mx-1)<m2)); % migration
        pM1B(mx)=0;
        
        % INTERACTION AB
        pAB = sum((Mq1==1 & Mq1c==2).*(rand(round(N*q1),mx)<c12)); % competition
        pAB1 = sum((Mq1==2 & Mq1c==1).*(rand(round(N*q1),mx)<c12)); % competition 
        
        pBA = sum((Mq1==1 & Mq1c==2).*(rand(round(N*q1),mx)<c21)); % competition
        pBA1 = sum((Mq1==2 & Mq1c==1).*(rand(round(N*q1),mx)<c21)); % competition 
       
        
        rPcA = [0,pP1A(1:mx-1)];
        rMcA = [0,pM1A(1:mx-1)];
        rPcB = [0,pP1B(1:mx-1)];
        rMcB = [0,pM1B(1:mx-1)];
        
        % update population sizes 
        pn = pn-pAA+pAE+pEA-pA-pP1A+rPcA+pM1A-rMcA-pAB-pAB1;
        pm = pm-pBB+pBE+pEB-pB-pP1B+rPcB+pM1B-rMcB-pBA-pBA1;
        
        pn(pn<0) = 0;
        pm(pm<0) = 0;
end
        