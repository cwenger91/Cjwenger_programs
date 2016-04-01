function [ ROC_args ] = RateOfClimb( Weight,Thrust,span,Sref,ct,target_alt,Aero)
%Calculates the rate of climb parameters for an aircraft in the climbing
%phase
%NOTE: Units are in english system, also if you need to change any
%constants or aerodynamic parameters do so in the sections below
%----------------------------Inputs---------------------------------------%
%       Weight --> Inital weight                                          %
%       span --> Wing span                                                %
%       Sref --> Wing reference area                                      %
%       Thrust --> Thrust for aircraft during climb                       %
%       ct --> Thrust specific fuel consumption                           %
%       target_alt --> Target altitude for climb                          %
%       Aero(1) --> CD0                                                   %
%       Aero(2) --> Cl                                                    %
%       Aero(3) --> Cd                                                    %
%-----------------------------Outputs-------------------------------------%
%       dw --> Change in weight                                           %
%       dx --> Change is distance                                         %
%       dt --> Change in time                                             %
%       ROC --> Rate of climb for aircraft                                %
%% CONSTANTS 
g = 32.2; %ft
TEMP_SL = 518.69; % R
DENSITY_SL = 2.3769E-3;
POWER_SL = 1000;
R = 1716; %ft*lb/(R*slug)
aT = -0.00356; %R/ft  
AR = span.^2./Sref; % Aspect Ratio
K = 1/(pi.*0.7.*AR);
%% AERODYNAMIC PARAMETERS
CD0 = Aero(1);
Cl = Aero(2);
Cd = Aero(3);
ClCd = Cl/Cd;
ClCd_half = 3/4*(1/2*K*CD0^3)^(1/4);
LDmax = sqrt(1/(4*CD0*K));
%% RATE OF CLIMB CALCULATIONS
T_W = Thrust/Weight;
W_S = Weight/Sref;
Z = 1 + sqrt(1 + 3/(LDmax^2*T_W^2));
height = [0:1000:target_alt];
power = 0;
for ii = 1:length(height)
    temp = TEMP_SL + aT.*height(ii);
    density = DENSITY_SL.*(temp./TEMP_SL).^(-g./(aT.*R)-1);
    ROCmax(ii) = sqrt(W_S*Z/(3*density*CD0))*T_W^(3/2)*(1-Z/6-3/(2*T_W^2*LDmax^2*Z));
    power = power + (density./DENSITY_SL).^0.7*POWER_SL;
end
dt = 0;
for ii = 1:length(height)
   dt = dt + 1000/ROCmax(ii);
end
dw = sum(ROCmax)/power;
dx = 2./ct*sqrt(2./(density.*span)).*ClCd_half*(sqrt(Weight)-sqrt(dw));
ROC_args = [dw;dx;dt];
