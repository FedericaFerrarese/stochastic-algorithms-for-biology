% SIMULATIONS
% INPUT
% M = vector of initial random distribution of the population 
% pn = vector of predator population in each cell time i 
% pm = vector of prey population in each cell time i 
% OUTPUT
% pn = vector of predator population in each cell time i+1
% pm = vector of prey population in each cell time i+1
function [pn,pm] = solPPslowNH2D(M,pn,pm)
global mx my m1 m2 p1 p2 b d1 d2 q1 q2
    
pA=(M(:,:,3)==1).*(rand(mx,my,1)<d1*(1-q1-q2));
pB=(M(:,:,3)==2).*(rand(mx,my,1)<d2*(1-q1-q2));


pBE=(M(:,:,1)==2 & M(:,:,2)==0).*(rand(mx,my,1)<b*q1);
pEB=(M(:,:,1)==0 & M(:,:,2)==2).*(rand(mx,my,1)<b*q1);

pAB = (M(:,:,1)==2 & M(:,:,2)==1).*(rand(mx,my,1)<p1*q1);
pAB1 = (M(:,:,1)==1 & M(:,:,2)==2).*(rand(mx,my,1)<p1*q1);

pBA = (M(:,:,1)==2 & M(:,:,2)==1).*(rand(mx,my,1)<p2*q1);
pBA1 = (M(:,:,1)==1 & M(:,:,2)==2).*(rand(mx,my,1)<p2*q1);

% righe
pP1Ar=(M(1:mx-1,:,4)==1 & M(2:mx,:,4)==0).*(rand(mx-1,my,1)<m1*q2);
pM1Ar=(M(2:mx,:,4)==1 & M(1:mx-1,:,4)==0).*(rand(mx-1,my,1)<m1*q2);
pP1Ar(mx,:)=0;
pM1Ar(mx,:)=0;

% righe
pP1Br=(M(1:mx-1,:,4)==2 & M(2:mx,:,4)==0).*(rand(mx-1,my,1)<m2*q2);
pM1Br=(M(2:mx,:,4)==2 & M(1:mx-1,:,4)==0).*(rand(mx-1,my,1)<m2*q2);
pP1Br(mx,:)=0;
pM1Br(mx,:)=0;

% colonne
pP1Ac=(M(:,1:my-1,4)==1 & M(:,2:my,4)==0).*(rand(mx,my-1,1)<m1*q2);
pM1Ac=(M(:,2:my,4)==1 & M(:,1:my-1,4)==0).*(rand(mx,my-1,1)<m1*q2);
pP1Ac(:,my)=0;
pM1Ac(:,my)=0;

pP1Bc=(M(:,1:my-1,4)==2 & M(:,2:my,4)==0).*(rand(mx,my-1,1)<m2*q2);
pM1Bc=(M(:,2:my,4)==2 & M(:,1:my-1,4)==0).*(rand(mx,my-1,1)<m2*q2);
pP1Bc(:,my)=0;
pM1Bc(:,my)=0;



rPcAr = zeros(mx,my);
rPcAr(2:mx,:)=pP1Ar(1:mx-1,:);
rMcAr = zeros(mx,my);
rMcAr(2:mx,:)=pM1Ar(1:mx-1,:);

rPcAc = zeros(mx,my);
rPcAc(:,2:my)=pP1Ac(:,1:my-1);
rMcAc = zeros(mx,my);
rMcAc(:,2:my)=pM1Ac(:,1:my-1);

rPcBr = zeros(mx,my);
rPcBr(2:mx,:)=pP1Br(1:mx-1,:);
rMcBr = zeros(mx,my);
rMcBr(2:mx,:)=pM1Br(1:mx-1,:);


rPcBc = zeros(mx,my);
rPcBc(:,2:my)=pP1Bc(:,1:my-1);
rMcBc = zeros(mx,my);
rMcBc(:,2:my)=pM1Bc(:,1:my-1);

pn=pn-pA+pAB+pAB1+round((-pP1Ar+rPcAr+pM1Ar-rMcAr-pP1Ac+rPcAc+pM1Ac-rMcAc)/2);

pm=pm+pBE+pEB-pB-pBA-pBA1-pAB-pAB1+round((-pP1Br+rPcBr+pM1Br-rMcBr-pP1Bc+rPcBc+pM1Bc-rMcBc)/2);
   
pn(pn<0)=0;
pm(pm<0)=0;
end
