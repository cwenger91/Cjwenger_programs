function [ TakeOff_args ] = TakeOffAnalysis( WTO,span,Sref,CDO,Clmax,T )
%TakeOffAnalsis gets the distance required for an aircraft to acclerate
%from velocity equal to zero to the belocity at takeoff speed to climb over
%a 35 or 50 ft obstacle
% NOTE: All calculations are done in english units if you want to change to
% SI units the constants section of the code needs to be set up for the SI
% system
%-------------------------------Inputs------------------------------------%
%     WTO --> Weight at takeoff                                           %
%     span --> Wing span                                                  %
%     Sref --> refrence area of wing                                      %
%     CDO --> coefficent of drag at lift equals zero                      %
%     Clmax --> Maximum coefficent of lift                                %
%     T --> Thrust at sea level                                           %
%-------------------------------Outputs-----------------------------------%
%     total_distance --> Total landing distance                           %
%     total_time --> Total time to takeoff                                %
%     sG --> Ground distance                                              %
%     sR --> Rotation distance                                            %
%     sTR --> Transition Distance                                         %
%     tG --> Ground time                                                  %
%     tTR --> Transtion time                                              %
%% CONSTANTS DEFINED FOR PROBLEM STATEMENT
g = 32.2; %ft/s^2
density = 0.0023769; % slugs/ft^3
u = 0.015; % friction 
hOB = 50; % ft
AR = span.^2./Sref;
K = 1/(pi.*0.7.*AR);
CD_LandingGear = 0.014; % coeff. of drag for the landing gear
%% AERODYNAMICS CALCULATIONS
Vstall = sqrt((WTO./Sref).*2/(density.*Clmax)); % Stall velocity
VTO = 1.2.*Vstall; % Takeoff velocity
V = 0.707.*VTO; % Runway velocity
D = 0.5.*density.*V.^2*Sref.*(CDO + CD_LandingGear+ K.*Clmax.^2);
L = 0.5.*density.*V.^2*Sref.*Clmax;
a = g./WTO.*(T-D-u.*(WTO-L));
%% RUNWAY DISTANCE AND TIME
Radius = VTO.^2/(0.15.*g);
htr = Radius-hOB;
theta_CL = 0.05; % Takeoff angle in rads
sG = 0.5*VTO.^2/a;
sR = 2.*VTO;
sTR = Radius.*sin(theta_CL);
if Radius > 50
    sCL = 0;
else
sCL = (50-htr)./tan(theta_CL);
end
tan(theta_CL)
total_distance = sG+sR+sTR+sCL;

tG = VTO./a;
tTR = (sTR+(sCL./cos(theta_CL)))./VTO;
total_time = tG + tTR;
TakeOff_args = [total_distance; total_time;sG;sR;sTR;sCL;tG;tTR];
end

