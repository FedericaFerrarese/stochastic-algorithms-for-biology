% HOMOGENEOUS PREDATOR-PREY MODEL: new approach, Gillespie and classical
% algoritms
clc
clear all
close all
% COMPARISON COMPUTATIONAL COSTS
% k=1--> comparison between the computational costs as the parameters vary
% for N fixed
% k=2 --> comparison between the computational costs as N varies for fixed
% parameters (=0.1)
% k=3 --> comparison between the computational costs as N varies for fixed
% parameters (=0.9)

% ORDER --> k=4

k=1;

if k==1
load('HS8.mat')
load('HF8.mat')
load('HG8.mat')

par=0.1:0.1:1;

figure
semilogy(par,timeF,'r-*','linewidth',3)
hold on
semilogy(par,timeG,'b-*','linewidth',3)
semilogy(par,time,'k-*','linewidth',3)
xlabel('Competition parameters')
ylabel('Time (s)')
legend('Our algorithm','Gillespie','Classical algorithm')
title('Computational cost')
set(gca,'FontSize',12,'FontWeight','bold')
axis([0.1 1 1e-3 1e3])

end

if k==2
load('HFN01.mat')
load('HGN01.mat')

figure
plot(Nrange,timeF,'r-*','linewidth',3)
hold on
plot(Nrange,timeG,'b-*','linewidth',3)
xlabel('N')
ylabel('Time')
legend('Our algorithm','Gillespie')
title('Parameters=0.1')
end

if k==3
load('HFN09.mat')
load('HGN09.mat')

figure
plot(Nrange,timeF,'r-*','linewidth',3)
hold on
plot(Nrange,timeG,'b-*','linewidth',3)
xlabel('N')
ylabel('Time')
legend('Our algorithm','Gillespie')
title('Parameters=0.9')
end

if k==4
   load('HFerror.mat')
    figure
loglog(Nrange,errN,'r*-','linewidth',3)
hold on
% loglog(Nrange,errN,'r*-','linewidth',3)
loglog(Nrange,errN(1).*(Nrange/Nrange(1)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
    legend({'Monte Carlo error',...
'N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Monte Carlo error predators')
set(gca,'FontSize',12,'FontWeight','bold')


 figure
loglog(Nrange,errM,'b*-','linewidth',3)
hold on
% loglog(Nrange,errN,'r*-','linewidth',3)
loglog(Nrange,errM(2).*(Nrange/Nrange(2)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
    legend({'Monte Carlo error',...
'N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Monte Carlo error preys')
set(gca,'FontSize',12,'FontWeight','bold')
end
