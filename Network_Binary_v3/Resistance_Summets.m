function [Resistivity] = Resistance_Summets(h)
%% Compute the constants
%{
Author = Michael Zhang
Date created = 03-06-18
%}

%{
Calculates the Resistance from summets2018.
%}

%{
Inputs:
    h = 1 - epsilon(normalised)
Outputs:
   	Resistivity = Resistance per length per area? TODO
%}

%% Parameters
%TODO manage parameters
% close all
% global chi
chi = 250;

phi_f = 0.99; % Fluid phase fraction
% phi_s = 1 - phi_f; % Solid fraction
% phi = phi_s/phi_f;
% h = .8;
dx = 0.01; % TODO.
x_l_range = 0:dx:h;
x_f_range =  h:dx:1;


%% Calculate constants
% [A,aBcs,D1,D2,E,M] = Calculate_Constants(chi,h,phi_f)
[A,aBcs,D1,D2,~,~] = Calculate_Constants(chi,h,phi_f);

%% Calculate flow rate

v_f = @(x) (D1*exp(-sqrt(chi).*x)+D2*exp(sqrt(chi).*x)-(A+aBcs)/chi).^2;
v_l = @(x) ((A*x.^2)/2 + 1).^2;

% u = @(x) -D1*exp(-sqrt(chi).*x)-D2*exp(sqrt(chi).*x) + (A*(phi+1)*x.^2)/2 + 1 + E*x+M;

% figure 
% plot(x_l_range,v_l(x_l_range))
% 
% figure 
% plot(x_f_range,v_f(x_f_range))
% % 
% % figure
% % plot(x_f_range,u(x_f_range))
% 
% x = [x_l_range, x_f_range];
% y = [v_l(x_l_range),v_f(x_f_range)];
% figure
% plot([fliplr(-x),x],[fliplr(y),y])

%% Calculate resistance

% Flowrate as the volume integral for the 1D graph.
Flowrate = (integral(v_f,h,1) + integral(v_l,0,h)) * 2 * pi;

% 
Resistivity = -A/Flowrate;

end

