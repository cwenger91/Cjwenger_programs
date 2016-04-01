function [Ae,Ath,Patm1,Patm2,Patm3] = NozzleSolution( e_x,th_x,Pt,gama)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

g1 = (gama+1)/(2*(gama-1));
g2 = gama/(gama-1);

Ae = 3+2.8*e_x/6*(2+e_x/6);
Ath = 3+2.8*th_x/6*(2+th_x/6);
%Area ratio of Ae/A*
rA = Ae/Ath;

% 1st Critical Calculations
Me1 = (2/(gama+1))^g1/rA;
rPt1 = (1+(gama-1)/2*Me1^2)^g2;
Patm1 = Pt/rPt1;

% 3rd Critical Calculations
Me3 = ((gama+1)/(gama-1))^((gama+1)/4)*rA^((gama-1)/2);
rPt3 = (1+(gama-1)/2*Me3^2)^g2;
Patm3 = Pt/rPt3;

% 2nd Critical Calculations
M1 = ((gama+1)/(gama-1))^((gama+1)/4)*rA^((gama-1)/2);
rPt = (1+(gama-1)/2*M1^2)^g2;
M2 = sqrt((M1^2+2/(gama-1))/(2*gama/(gama-1)*M1^2-1));
rP2_1 = (1+gama*M1^2)/(1+gama*M2^2);
rPt2 = rPt/rP2_1;
Patm2 = Pt/rPt2;

end

