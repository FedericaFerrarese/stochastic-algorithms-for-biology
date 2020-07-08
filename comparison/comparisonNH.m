% HETEROGENEOUS PREDATOR-PREY MODEL (1D-case): new approach, Gillespie and
% classical algoritms
clc
clear all
close all
% COMPARISON COMPUTATIONAL COSTS
% k=1--> comparison between the computational costs as the parameters vary
% for N fixed
% k=2 --> comparison between the computational costs as N varies for fixed
% migration parameters (=0.1)
% k=3 --> comparison between the computational costs as N varies for fixed
% migration parameters (=0.9)

% ORDER --> k=4
k=1
if k==1
load('NHSpar.mat')
load('NHFpar.mat')
load('NHGpar.mat')

par=0.1:0.1:1;

figure
semilogy(par,timeF,'r-*','linewidth',3)
hold on
semilogy(par,time,'b-*','linewidth',3)
semilogy(par,timeS,'k-*','linewidth',3)
xlabel('Migration parameters')
ylabel('Time (s)')
legend('New algorithm','Gillespie algorithm','Classical algorithm')
title('Computational cost')
set(gca,'FontSize',12,'FontWeight','bold')
axis([0.1 1 1e-1 1e3])
end

if k==2
load('NHFN01.mat')
load('NHGN01.mat')

figure
semilogy(Nrange,timeF,'r-*','linewidth',3)
hold on
semilogy(Nrange,time,'b-*','linewidth',3)
xlabel('N')
ylabel('Time (s)')
legend('New algorithm','Gillespie algorithm')

title('Computational cost m_1 = m_2 = 0.1')
set(gca,'FontSize',12,'FontWeight','bold')
axis([100 1e3 1e-1 1e3])

end

if k==3
load('NHFN09.mat')
load('NHGN09.mat')

figure
semilogy(Nrange,timeF,'r-*','linewidth',3)
hold on
semilogy(Nrange,time,'b-*','linewidth',3)
xlabel('N')
ylabel('Time (s)')

legend('New algorithm','Gillespie algorithm')

title('Computational cost m_1 = m_2 = 0.9')
set(gca,'FontSize',12,'FontWeight','bold')
axis([100 1e3 1e-1 1e3])
end

if k==4
   load('NHerror.mat')
    figure
loglog(Nrange,errN,'r*-','linewidth',3)
hold on
loglog(Nrange,errN(2).*(Nrange/Nrange(2)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
    legend({'Monte Carlo error predators',...
'N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Monte Carlo error predators')
set(gca,'FontSize',12,'FontWeight','bold')
axis([10 10000 1e-3 10])

 figure
loglog(Nrange,errM,'b*-','linewidth',3)
hold on
loglog(Nrange,errM(3).*(Nrange/Nrange(3)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
    legend({'Monte Carlo error preys',...
'N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Monte Carlo error preys')
set(gca,'FontSize',12,'FontWeight','bold')
axis([10 10000 1e-3 10])
end
