function [ M2,Pt2_Pt1 ] = RayleighFlow( T1,M1,mDot,mfDot,h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
gama = 1.4; % constant for air
R = 287;
Cp = gama*R/(gama-1);
f = mfDot/mDot;
QDot = f*h;
Tt1 = (1 + (gama -1)/2*M1^2)*T1;

Tt2_Tt1 = 1+(f*h)/(Cp*Tt1);
fRM1 = (1+(gama-1)/2*M1^2)/((1+gama*M1^2)^2)*M1^2;
f2 = Tt2_Tt1*fRM1
M(1) = sqrt(1-2*gama*f2 + sqrt((2*gama*f2-1)^2-4*f2*(f2*gama^2-(gama-1)/2))...
     /(2*f2*gama^2+1-gama))
M(2) = sqrt(1-2*gama*f2 - sqrt((2*gama*f2-1)^2-4*f2*(f2*gama^2-(gama-1)/2))...
     /(2*f2*gama^2+1-gama))
 
 if(M1 < 1)
     if(QDot > 0)
         M2 = max(M);
     elseif(QDot < 0)
         M2 = min(M)
     end
 elseif(M1>1)
     if(QDot > 0)
         M2 = min(M);
     elseif(QDot < 0)
         M2 = max(M)
     end
 end
Pt2_Pt1 = ((1+(gama-1)/2*M2^2)/(1+(gama-1)/2*M1^2))^(gama/(gama-1))...
    *((1+gama*M1^2)/(1+gama*M2^2))
end

