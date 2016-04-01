function [ M2, PtRatio ] = ObliqueShock(M1,delta,sigma,gama)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

M1n = M1*sind(sigma);
M2n = sqrt((M1n^2+2/(gama-1)) / (2*gama/(gama-1)*M1n^2-1));
M2 = M2n / sind(sigma-delta);

PtRatio = (((gama+1)*M1n^2)/((gama-1)*M1n^2+2))^(gama/(gama-1))*((gama+1)/(2*gama*M1n^2+1-gama))^(1/(gama-1));
end

