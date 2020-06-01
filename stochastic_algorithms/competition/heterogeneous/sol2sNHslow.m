% SIMULATIONS
% INPUT
% M = vector of initial random distribution of the population  
% pn = vector of population A in each cell time i 
% pm = vector of population B in each cell time i
% OUTPUT
% pn = vector of population A in each cell time i+1
% pm = vector of population B in each cell time i+1
function [pn,pm] = sol2sNHslow(M,pn,pm)
global  mx m1 m2 b1 b2 c11 c22 c12 c21 d1 d2 q1 q2
      
        pAA=(M(1,:)==1 & M(2,:)==1).*(rand(1,mx)<c11*q1);
        pAE=(M(1,:)==0 & M(2,:)==1).*(rand(1,mx)<b1*q1);
        pEA=(M(1,:)==1 & M(2,:)==0).*(rand(1,mx)<b1*q1);
        pA=(M(3,:)==1).*(rand(1,mx)<d1*(1-q1-q2));
        
        pP1A=(M(4,1:mx-1)==1 & M(4,2:mx)==0).*(rand(1,mx-1)<m1*q2);
        pM1A=(M(4,1:mx-1)==0 & M(4,2:mx)==1).*(rand(1,mx-1)<m1*q2);
        pP1A(mx)=0;
        pM1A(mx)=0;
        
        
          
        pBB=(M(1,:)==2 & M(2,:)==2).*(rand(1,mx)<c22*q1);
        pBE=(M(1,:)==0 & M(2,:)==2).*(rand(1,mx)<b2*q1);
        pEB=(M(1,:)==2 & M(2,:)==0).*(rand(1,mx)<b2*q1);
        pB=(M(3,:)==2).*(rand(1,mx)<d2*(1-q1-q2));
        
        pP1B=(M(4,1:mx-1)==2 & M(4,2:mx)==0).*(rand(1,mx-1)<m2*q2);
        pM1B=(M(4,1:mx-1)==0 & M(4,2:mx)==2).*(rand(1,mx-1)<m2*q2);
        pP1B(mx)=0;
        pM1B(mx)=0;
        
        pAB = (M(1,:)==2 & M(2,:)==1).*(rand(1,mx)<c12*q1);
        pAB1 =(M(1,:)==1 & M(2,:)==2).*(rand(1,mx)<c12*q1);
        
        pBA = (M(1,:)==2 & M(2,:)==1).*(rand(1,mx)<c21*q1);
        pBA1 =(M(1,:)==1 & M(2,:)==2).*(rand(1,mx)<c21*q1);
        
        % update population sizes
        pn=pn-pAA+pAE+pEA-pA-pP1A+rPcA+pM1A-rMcA-pAB-pAB1;      
        pm=pm-pBB+pBE+pEB-pB-pP1B+rPcB+pM1B-rMcB-pBA-pBA1;
             
        pn(pn<0)=0;
        pm(pm<0)=0;
end
