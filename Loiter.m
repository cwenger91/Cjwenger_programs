function [ Loiter_args ] = Loiter( Weight,span,Sref,ct,alt,time )
% Calculates the changes in weight and time for a given aircraft in loiter
%NOTE: Units are in english system, also if you need to change any
%constants or aerodynamic parameters do so in the sections below
%----------------------------Inputs---------------------------------------%
%       Weight --> Inital weight                                          %
%       span --> Wing span                                                %
%       Sref --> Wing reference area                                      %
%       ct --> Thrust specific fuel consumption                           %
%       alt --> Loiter altitude                                           %
%       time --> Loiter time                                              %
%-----------------------------Outputs-------------------------------------%
%       dw --> Change in weight                                           %
%       dx --> Change is distance                                         %
%       dt --> Change in time                                             %
%% CONSTANTS 
g = 32.2; %ft
TEMP_SL = 518.69; % R
DENSITY_SL = 2.3769E-3;
R = 1716; %ft*lb/(R*slug)
aT = -0.00356; %R/ft  
AR = span.^2./Sref; % Aspect Ratio
K = 1/(pi.*0.7.*AR);
%% AERODYNAMIC PARAMETERS
CD0 = 1;
Cl = 1;
Cd = 2;
ClCd_half = 3/4*(1/2*K*CD0^3)^(1/4);
ClCd = Cl/Cd;
%% ATMOSPHERE PROPERTIES
temp = TEMP_SL + aT.*alt;
density = DENSITY_SL.*(temp./TEMP_SL).^(-g./(aT.*R)-1);
%% LOITER CALCULATIONS
dw = Weight/exp(time/(1./ct.*ClCd))
dx = 2./ct*sqrt(2./(density.*span)).*ClCd_half*(sqrt(Weight)-sqrt(dw));
dt = time;
Loiter_args = [dw;dx;dt];
end

