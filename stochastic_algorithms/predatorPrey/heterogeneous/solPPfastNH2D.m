% SIMULATIONS
% INPUT
% M1 = vector of initial random distribution of the population 
% M2 = vector of initial random distribution of the population 
% pn = vector of population A in each cell at time i 
% pm = vector of population B in each cell at time i
% OUTPUT
% pn = vector of population A in each cell at time i+1
% pm = vector of population B in each cell at time i+1

function [pn,pm] = solPPfastNH2D(M1,M2,pn,pm)
global N mx my m1 m2 p1 p2 b d1 d2 q1 q2
        % q1: birth + competition
        Mq1 = M1(:,:,1:round(N*q1));
        Mq1c = M2(:,:,1:round(N*q1));
        
        % q2: migration
        Mq2 = M1(:,:,round(N*q1)+1:round(N*q1+N*q2));

        % q3: death
        Mq3 = M1(:,:,round(N*q1+N*q2)+1:N);
        
        
        pEB = sum((Mq1==2 & Mq1c==0).*(rand(mx,my,round(N*q1))<b),3); 
        pBE = sum((Mq1==0 & Mq1c==2).*(rand(mx,my,round(N*q1))<b),3);
        
        pAB = sum((Mq1==1 & Mq1c==2).*(rand(mx,my,round(N*q1))<p1),3); 
        pAB1 = sum((Mq1==2 & Mq1c==1).*(rand(mx,my,round(N*q1))<p1),3);     
        pBA = sum((Mq1==1 & Mq1c==2).*(rand(mx,my,round(N*q1))<p2),3);  
        pBA1 = sum((Mq1==2 & Mq1c==1).*(rand(mx,my,round(N*q1))<p2),3);

        pA = sum((Mq3==1).*(rand(mx,my,round(N*(1-q1-q2)))<d1),3);

        pB = sum((Mq3==2).*(rand(mx,my,round(N*(1-q1-q2)))<d2),3);
        
        pP1Ar =  sum((Mq2(1:mx-1,:,:)==1 & Mq2(2:mx,:,:)==0).*(rand(mx-1,my,size(Mq2,3))<m1),3); 
        pP1Ar(mx,:)=0;

        % EiAj in AiEj
        pM1Ar =  sum((Mq2(1:mx-1,:,:)==0 & Mq2(2:mx,:,:)==1).*(rand(mx-1,my,size(Mq2,3))<m1),3); 
        pM1Ar(mx,:)=0;

        pP1Ac =  sum((Mq2(:,1:my-1,:)==1 & Mq2(:,2:my,:)==0).*(rand(mx,my-1,size(Mq2,3))<m1),3); 
        pP1Ac(:,my)=0;

        pM1Ac =  sum((Mq2(:,1:my-1,:)==0 & Mq2(:,2:my,:)==1).*(rand(mx,my-1,size(Mq2,3))<m1),3); 
        pM1Ac(:,my)=0;
      
        pP1Br =  sum((Mq2(1:mx-1,:,:)==2 & Mq2(2:mx,:,:)==0).*(rand(mx-1,my,size(Mq2,3))<m2),3); 
        pP1Br(mx,:)=0;

        pM1Br =  sum((Mq2(1:mx-1,:,:)==0 & Mq2(2:mx,:,:)==2).*(rand(mx-1,my,size(Mq2,3))<m2),3); 
        pM1Br(mx,:)=0;

        pP1Bc =  sum((Mq2(:,1:my-1,:)==2 & Mq2(:,2:my,:)==0).*(rand(mx,my-1,size(Mq2,3))<m2),3); 
        pP1Bc(:,my)=0;
  
        pM1Bc =  sum((Mq2(:,1:my-1,:)==0 & Mq2(:,2:my,:)==2).*(rand(mx,my-1,size(Mq2,3))<m2),3); 
        pM1Bc(:,my)=0;
        
        rPcAr=zeros(mx,my);
        rMcAr=zeros(mx,my);
        rPcAc=zeros(mx,my);
        rMcAc=zeros(mx,my);
        
        rPcAr(2:mx,:) =pP1Ar(1:mx-1,:);
        rMcAr(2:mx,:) =pM1Ar(1:mx-1,:);

        rPcAc(:,2:my) =pP1Ac(:,1:my-1);
        rMcAc(:,2:my) =pM1Ac(:,1:my-1);
        
        rPcBr=zeros(mx,my);
        rMcBr=zeros(mx,my);
        rPcBc=zeros(mx,my);
        rMcBc=zeros(mx,my);
        
        rPcBr(2:mx,:) =pP1Br(1:mx-1,:);
        rMcBr(2:mx,:) =pM1Br(1:mx-1,:);
        rPcBc(:,2:my) =pP1Bc(:,1:my-1);
        rMcBc(:,2:my) =pM1Bc(:,1:my-1);
         
        pn=pn-pA+pAB+pAB1+round((-pP1Ac+rPcAc+pM1Ac-rMcAc-pP1Ar+rPcAr+pM1Ar-rMcAr));
        pm=pm+pBE+pEB-pB-pBA-pBA1-pAB-pAB1+round((-pP1Bc+rPcBc+pM1Bc-rMcBc-pP1Br+rPcBr+pM1Br-rMcBr));

        pn(pn<0)=0;
        pm(pm<0)=0;
        
end
        