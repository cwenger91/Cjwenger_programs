function [ M2,PtRatio ] = NormalShock(M1,gama)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M2 = sqrt((M1^2+2/(gama-1)) / (2*gama/(gama-1)*M1^2-1));
PtRatio = (((gama+1)*M1^2)/((gama-1)*M1^2+2))^(gama/(gama-1))*((gama+1)/(2*gama*M1^2+1-gama))^(1/(gama-1));
end

