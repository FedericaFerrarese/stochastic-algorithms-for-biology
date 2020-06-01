% SIMULATIONS
% INPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% pn = vector of population in each cell at time i  
% OUTPUT
% pn = vector of population in each cell time i+1
function pn = sol1sNHfast(M1,M2,pn)
global N mx m b c d q1 q2
        % q1: birth + competition
        Mq1 = M1(1:round(N*q1),:);
        Mq1c = M2(1:round(N*q1),:);
        
        % q2: migration
        Mq2 = M1(round(N*q1)+1:round(N*q1+N*q2),:);

        % q3: death
        Mq3 = M1(round(N*q1+N*q2)+1:N,:);
        
        % AjAj in AjEj
        pAA = sum((Mq1==1 & Mq1c==1).*(rand(round(N*q1),mx)<c)); % competition
        
        % AjEj in AjAj
        pAE = sum((Mq1==1 & Mq1c==0).*(rand(round(N*q1),mx)<b)); % birth
        
        % EjAj in AjAj
        pEA = sum((Mq1==0 & Mq1c==1).*(rand(round(N*q1),mx)<b)); % birth 
    
        % Aj in Ej
        pA = sum((Mq3==1).*(rand(round(N*q1),mx)<d)); % death

        % from cell j to j+1
        % EiAj in AiEj
        pP1 = sum((Mq2(:,1:mx-1)==1 & Mq2(:,2:mx)==0).*(rand(size(Mq2,1),mx-1)<m)); % migration
        pP1(mx)=0;

        % from cell j to j-1
        % EiAj in AiEj
        pM1 = sum((Mq2(:,2:mx)==1 & Mq2(:,1:mx-1)==0).*(rand(size(Mq2,1),mx-1)<m)); % migration
        pM1(mx)=0;
        
        rPc = [0,pP1(1:mx-1)];
        rMc = [0,pM1(1:mx-1)];
        
        % update population size
        pn=pn-pAA+pAE+pEA-pA-pP1+rPc+pM1-rMc;
        pn(pn<0)=0;

end
        