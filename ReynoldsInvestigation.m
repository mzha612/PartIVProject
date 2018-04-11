clear all
close all
clc

mu = 1.002; % kg m^-1 s^-1
rho = 1000; % kg m^-3

Re = @(u,L) rho * u * L / mu;

% Arterial Mean Velocity
x_art = linspace(20,100,81)*1e-6;
u = @(x) (0.054.*x/1e-6 + 5.1)/1000;

[x,y] = meshgrid(x_art,u(x_art));
ReMap = Re(u(x_art)',x_art);

figure
contourf(x,y,ReMap)
colorbar
xlabel('Arterial vessel Diameter (m)')
ylabel('Arterial mean velocity (m/s)')
% Venous Mean Velocity
x_ven = linspace(20,110,91)*1e-6;
u = @(x) (0.065.*x/1e-6 + 3.25)/1000;

[x,y] = meshgrid(x_ven,u(x_ven));
ReMap = Re(u(x_ven)',x_ven);

figure
contourf(x,y,ReMap)
colorbar
xlabel('Veinous vessel Diameter (m)')
ylabel('Veinous mean velocity (m/s)')
%