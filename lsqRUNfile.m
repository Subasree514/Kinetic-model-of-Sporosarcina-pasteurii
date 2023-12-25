close all
clear all
x_n=zeros(6,6);%simulation values
%% experiment data and std deviation data
load Data;%%experimental data
load sddata;%%Standard Deviation values
Expdata=[Data(:,1),Data(:,2),Data(:,3),Data(:,4),Data(:,5),Data(:,6)];%c_measured
sd=[sddata(:,1),sddata(:,2),sddata(:,3),sddata(:,4),sddata(:,5),sddata(:,6)]; %sigmanew
%% initial parameters and scaling factors
parameters=[0.82,0.72,0.75,0.82,0.61,0.5,0.8];%scaled parameters
%parameters=0.3*ones(7,1);
sf=[0.05 0.4 0.95 0.95 7.67 610 6000];%scaling factors
%% Bounds
lb=0.1*ones(7,1);%%lower bounds
ub=1*ones(7,1);%%upper bounds

%% sigma values(2)
sigmanew=[sd(:,1);sd(:,2);sd(:,3);sd(:,4);sd(:,5);sd(:,6)];%%standard deviation from experiments
%sigmanew=[15*ones(12,1);400*ones(6,1);3e3*ones(6,1);9e2*ones(6,1);100*ones(6,1)];%%Range of experimental data
%%
tspan=(0:1:5)';%%time points
x_n0=[Data(1,1),Data(1,2),Data(1,3),Data(1,4),Data(1,5),Data(1,6)];%Zeroth day data
options=optimoptions('lsqnonlin','StepTolerance',1e-8,'FunctionTolerance',1e-8,'MaxFunctionEvaluation',20000);
[values,resnorm,ex,flag]=lsqnonlin(@(para)lsqfun(para,Expdata,sd,tspan,sigmanew,sf),parameters,lb,ub,options);
[t_n,x_n] = ode15s(@(t,x)ODEfile(t,x,values,sf),tspan,x_n0);%Call ODE file
%% plots
figure() 
subplot(2,3,1)
plot(t_n,x_n(:,1),'k-',t_n,Expdata(:,1),'o','MarkerEdgeColor','black')
title('Biomass')
xlabel({'Time','(days)'})
ylabel({'Concentration','(gDW/l)'})
 subplot(2,3,2)
 plot(t_n,x_n(:,2),'k-',t_n,Expdata(:,2),'o','MarkerEdgeColor','black')
 title('Glucose')
 xlabel({'Time','(days)'})
ylabel({'Concentration','(g/l)'})
 subplot(2,3,3)
 plot(t_n,x_n(:,3),'k-',t_n,Expdata(:,3),'o','MarkerEdgeColor','black')
  title('Extracellular Urea')
  xlabel({'Time','(days)'})
ylabel({'Concentration','(mM)'})
 subplot(2,3,4)
 plot(t_n,x_n(:,4),'k-',t_n,Expdata(:,4),'o','MarkerEdgeColor','black')
 title('Intracellular Urea')
 xlabel({'Time','(days)'})
ylabel({'Concentration','(mM)'})
 subplot(2,3,5)
 plot(t_n,x_n(:,5),'k-',t_n,Expdata(:,5),'o','MarkerEdgeColor','black')
 title('Intracellular Ammonium')
 xlabel({'Time','(days)'})
ylabel({'Concentration','(mM)'})
 subplot(2,3,6)
 plot(t_n,x_n(:,6),'k-',t_n,Expdata(:,6),'o','MarkerEdgeColor','black')
title('Extracellular Ammonium')
xlabel({'Time','(days)'})
ylabel({'Concentration','(mM)'})


