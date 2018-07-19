%% Change in Resistance with respect to radius for a constant egl thickness
%{
Simple script to show the effect of vessel radius and resistance for a
fixed egl thickness
%}
%{
Author = Michael Zhang
Date created = 10-07-18
%}
clc
clear all
close all
% digits(32)
%% Parameter
phi_f = 0.99;
mu_f = 10^-3; % viscosity of water, Pa.s
K = 10^10;

EGL = 1e-6; %thickness. m
num_R = 200; % Number of points
Radius = logspace(log10(0.05/2),log10(2.5e-6),num_R); % Range of radii

for i = 1:num_R
    h = 1 - (EGL)/Radius(i); % non-dim ratio of egl to radius
    chi = K*(Radius(i))^2/(phi_f * mu_f); % non-dim chi
    
    Resistance(i)=  mu_f /Radius(i).^3 *Resistance_Summets(h,chi); % Pa.s.m-3.m
    Resistance_Poiseuille(i) = PoiseuilleFlow(Radius(i));
end

figure
plot(Radius*10^6, Resistance./Resistance_Poiseuille)
xlabel('Radius (\mum)')
ylabel('increase in resistivity')

figure
semilogx(Radius*10^6, Resistance./Resistance_Poiseuille)
xlabel('Radius (\mum)')
ylabel('increase in resistivity')
