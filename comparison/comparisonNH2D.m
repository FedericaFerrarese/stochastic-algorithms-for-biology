% HETEROGENEOUS PREDATOR-PREY MODEL (2D-case): new approach, Gillespie and
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
% migration parameters (=0.5)

% ORDER --> k=4
k=1

% k=1 --> N = 500, b=d1=d2=p1=p2=0.1:0.1:1, t = 1000 --> F+G+S

if k==1
load('NHSpar2D.mat')
load('NHFpar2D.mat')
load('NHGpar2D.mat')

par=0.1:0.1:1;

figure
semilogy(par,timeF,'r-*','linewidth',3)
hold on
semilogy(par,timeG,'b-*','linewidth',3)
semilogy(par,timeS,'k-*','linewidth',3)
xlabel('Migration parameters')
ylabel('Time (s)')
legend('New algorithm','Gillespie','Classical algorithm')
title('Computational cost')
set(gca,'FontSize',12,'FontWeight','bold')
axis([0.1 1 1 1e3])
end

if k==2
load('NH2DFN01.mat')
load('NH2DGN01.mat')

figure
semilogy(Nrange,timeF,'r-*','linewidth',3)
hold on
semilogy(Nrange,timeG,'b-*','linewidth',3)
xlabel('N')
% ylabel('Time')
legend('New algorithm','Gillespie algorithm')

title('Computational cost m_1 = m_2 = 0.1')
set(gca,'FontSize',12,'FontWeight','bold')
axis([100 1e3 1e-1 1e3])
end


if k==3
load('NH2DFN05.mat')
load('NH2DGN05.mat')
figure
semilogy(Nrange,timeF,'r-*','linewidth',3)
hold on
semilogy(Nrange,timeG,'b-*','linewidth',3)
xlabel('N')
% ylabel('Time')
legend('New algorithm','Gillespie algorithm')

title('Computational cost m_1 = m_2 = 0.5')
set(gca,'FontSize',12,'FontWeight','bold')
axis([100 1e3 1e-1 1e3])

end

if k==4
load('NH2Derror.mat', 'Nrange')
load('NH2Derror.mat', 'errM')
load('NH2Derror.mat', 'errN')
    figure
loglog(Nrange,errN,'r*-','linewidth',3)
hold on
% loglog(Nrange,errN,'r*-','linewidth',3)
loglog(Nrange,errN(end-7).*(Nrange/Nrange(end-7)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
    legend({'Monte Carlo error preys',...
'N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Monte Carlo error predators')
set(gca,'FontSize',12,'FontWeight','bold')
axis([10 10000 1e-3 10])


 figure
loglog(Nrange,errM,'b*-','linewidth',3)
hold on
% loglog(Nrange,errN,'r*-','linewidth',3)
loglog(Nrange,errM(end-7).*(Nrange/Nrange(end-7)).^(-1/2),'k','linewidth',3)
xlabel('N')
ylabel('Error')
    legend({'Monte Carlo error preys',...
'N^{(-1/2)}'},...
'FontSize',12,'FontWeight','bold')
title('Monte Carlo error preys')
set(gca,'FontSize',12,'FontWeight','bold')
axis([10 10000 1e-3 10])
end

