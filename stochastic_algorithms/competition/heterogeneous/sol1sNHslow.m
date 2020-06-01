% SIMULATIONS
% INPUT
% M = vector of initial random distribution of the population 
% p = vector of population in each cell time i 
% OUTPUT
% p = vector of population in each cell time i+1
function p = sol1sNHslow(M,p)
global mx m b c d q1 q2

    % competition
    pAA = (M(1,:)==1 & M(2,:)==1).*(rand(1,mx)<c*q1); % competition
    
    % birth
    pAE = (M(1,:)==1 & M(2,:)==0).*(rand(1,mx)<b*q1); % birth
    pEA = (M(1,:)==0 & M(2,:)==1).*(rand(1,mx)<b*q1); % birth
    
    % death
    pA = (M(3,:)==1).*(rand(1,mx)<d*(1-q1-q2)); % death
    
    pP1A=(M(4,1:mx-1)==1 & M(4,2:mx)==0).*(rand(1,mx-1)<m*q2); % migration
    pM1A=(M(4,1:mx-1)==0 & M(4,2:mx)==1).*(rand(1,mx-1)<m*q2); % migration
    pP1A(mx)=0;
    pM1A(mx)=0;
   
    rPcA = zeros(1,mx);
    rPcA(2:mx)=pP1A(1:mx-1);
    rMcA = zeros(1,mx);
    rMcA(2:mx)=pM1A(1:mx-1);
    
    % update population size
    p = p-pAA+pAE+pEA-pA-pP1A+rPcA+pM1A-rMcA;    

    p(p<0)=0;
           
end
