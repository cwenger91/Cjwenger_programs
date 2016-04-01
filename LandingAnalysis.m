function [ Landing_args ] = LandingAnalysis( WL,span,Sref,Clmax,CDO )
% Calculates the landing distance in the horizonatal distance required to
% clear a 50 ft obstacle with free roll and then braking to a complete stop
% NOTE: All calculations are done in english units if you want to change to
% SI units the constants section of the code needs to be set up for the SI
% system
%-------------------------------Inputs------------------------------------%
%     WL --> Landing weight with 1/5 fuel weight                          %
%     span --> Wing span                                                  %
%     Sref --> refrence area of wing                                      %
%     CDO --> coefficent of drag at lift equals zero                      %
%     Clmax --> Maximum coefficent of lift                                %
%-------------------------------Outputs-----------------------------------%
%     total_distance --> Total landing distance                           %
%     sA --> Air distance                                                 %
%     sFR --> Free roll distance                                          %
%     sB --> Braking Distance                                             %
%% CONSTANTS DEFINED FOR PROBLEM STATEMENT
g = 32.2; %ft/s^2
density = 0.0023769; % slugs/ft^3
u = 0.015; % friction of the runway 
AR = span.^2./Sref; % Aspect Ratio
K = 1/(pi.*0.7.*AR);
CD_LandingGear = 0.014; % coeff. of drag for the landing gear
%% AERODYNAMICS CALCULATIONS
WTO = 551810; %Weight at takeoff
Vstall = sqrt((WTO./Sref).*2/(density.*Clmax)); % Stall velocity
V50 = 1.3.*Vstall; % Velocity over the 50-ft obstacle
VTD = 1.15.*Vstall; % Touchdown Velocity
CD = (CDO + CD_LandingGear + K.*Clmax.^2);
D = 0.5.*density.*VTD.^2*Sref.*CD;
L = WL;
%% LANDING DISTANCE
sA = L./D.*((V50.^2-VTD.^2)./(2.*g)+50);
% or sA = 50/tand(3) 
sFR = 3.*VTD;
sB = WL/(g.*u.*density.*Sref.*(CD./u-Clmax)).*log(1+density./2.*Sref./WL.*...
    (CD./u-Clmax).*VTD.^2);
total_distance = sA+sFR+sB;
Landing_args = [total_distance; sA; sFR; sB];
end

