function err=lsqfun(parameters,Data,sddata,tspan,sigmanew,s_f)%lsqfun file
load Data;%Experimental data-matfile
%% call ODE function
tspan=(0:1:5)';%time points
x_n=zeros(6,6);%simulation
x_n0=[Data(1,1),Data(1,2),Data(1,3),Data(1,4),Data(1,5),Data(1,6)];%Zeroth day data points
[t_n,x_n] = ode15s(@(t,x)ODEfile(t,x,parameters,s_f),tspan,x_n0);%Call ODE file
%% lsqnonlin function 
Expdata=[Data(:,1),Data(:,2),Data(:,3),Data(:,4),Data(:,5),Data(:,6)];%Experimental data points
C_measured=[Expdata(:,1);Expdata(:,2);Expdata(:,3);Expdata(:,4);Expdata(:,5);Expdata(:,6)];%Experimental
C_Predicted=[x_n(:,1);x_n(:,2);x_n(:,3);x_n(:,4);x_n(:,5);x_n(:,6)];%Simulated
error=(C_measured-C_Predicted)./sigmanew;%Objective function
err=error;
end
