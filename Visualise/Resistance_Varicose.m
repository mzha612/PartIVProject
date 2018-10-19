function [R] = Resistance_Varicose(Vessel,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% V_radius = Vessel.Radius;
V_h = Vessel.h;

R = 1.066*exp(9.0449*(1-V_h)); %cs1

if param == 0.1
    R = 1.2832*exp(5.4793*(1-V_h));
end
end

