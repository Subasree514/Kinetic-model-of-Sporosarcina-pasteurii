function f1 = ODEfile(t,x,parameters,sf)
f1=zeros(6,1);%ODE Simulation
x_n=zeros(6,6);
%% parameters=[k_m,v_max,ks,Umax,Y_XS,km,vmax];%%%
%% parameters=[0.041mmol/l 0.288mmol/gdw/d 0.716g/l 0.78/d 4.68gDW/g 305mmol/l 4800mmol/l/d];%% actual parameter values

%% parameters:
%% Kinetic Parameters from the literature to be estimated
k_m=parameters(1); %Michaelis Menten constant ammonium assimilation[mmol/l] 
v_max=parameters(2); %Michaelis Menten constant ammonium assimilation[mmol/gDW/d]
ks=parameters(3);%Saturation coefficient of glucose[mM]
%% Parameters from the experiment to be estimated
Umax=parameters(4);%specific growth rate constant[1/day]
Y_XS=parameters(5);%maximum yield of biomass on glucose[gDW/g] 
%% active transport constants
km=parameters(6);%Michaelis Menten constant ureolysis[mmol/l]
vmax=parameters(7);%Michaelis Menten constant ureolysis[mmol/l/d]
%% carrying capacity of the biomass equation
X_max=8.25;%gDW/l
%% States (Process variables)
X = x(1); % X: biomass  (gDW/l)
S = x(2);% S: substrate (g/l)
Urea_ex = x(3); % Urea extracellular(mM)
Urea_in = x(4); % Urea intracellular(mM)
Amm_in = x(5); % Ammonium intracellular(mM)
Amm_ex = x(6); % Ammonium extracellular(mM)
%% expressions of growth rate and some constants used
U_1=(Umax*sf(4)*S)/((ks*sf(3))+S);%(Growth rate)
n_cells =14.5e8*exp(U_1*t);%(per 100 ml)
%n_cells=[14500000 6800000 6966666.667 12233333.33 12233333.33 12233333.33 12233333.33]';%Cell numbers
surf_area=6e-12;%Surface area of cell membrane[m^2]
per_amm =4.1e-5;% permeability coefficient ammonia[m/d]
per_amm_cell=per_amm*surf_area; %Ammonia permeability coefficient*Area of the surface [m^3/d]
per_urea=1.6e-4;%permeability coefficient urea[m/d]
per_urea_cell=per_urea*surf_area;%Urea permeability coefficient*Area of the surface [m^3/d]
v_flask=1e-4;%Volume of flask in m^-3 to account for the uniiformity of the units (mmol/l=mol/m^3)
%% Rate expressions(4)
r0=(per_urea_cell*n_cells*(Urea_ex-Urea_in))./v_flask;%mol/m^3/d = mmol/l/d Urea diffusion
r1=(vmax*sf(7)*Urea_in)/((km*sf(6))+Urea_in);%% mmol/l/d Ureolysis
r2=(v_max*sf(2)*Amm_in)/((k_m*sf(1))+Amm_in);%%mmol/l/d Ammonium assimilation
r3=(per_amm_cell*n_cells*(Amm_in-Amm_ex))/v_flask;%%mol/m^3/d = mmol/l/d Passive transport of ammonium
%% Model Equations of 6 process variables
f1(1,1)=U_1*x(1)*(1-(x(1)/X_max));%Biomass gDW/l/d
f1(2,1)=-(U_1*x(1))/(Y_XS*sf(5));%Substrate g/l/d
f1(3,1)= -r0;%extracellular urea mmol/l/d
f1(4,1)=r0-r1;%Intracellular urea mmol/l/d
f1(5,1)= (2*r1)-(r2*X)-r3;%Intracellular ammonia mmol/l/d
f1(6,1)= r3;%extracellular ammonia mmol/l/d
end
