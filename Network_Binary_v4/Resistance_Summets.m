function [Resistivity] = Resistance_Summets(h)
%% Resistance_Summets v4
%{
Calculates the Resistance from appendix A in summets2018.
Does this by finding the values of the constants, then using the constants
to build up a flow velocity profile. The profile is then integrated to
obtain the flow rate. Resisitvity is finally obtained by dividing the
pressure gradient A with the flow rate.
%}

%{
Inputs:
    h = 1 - epsilon(normalised)
Outputs:
   	Resistivity     Non-dimensional resisitivity (resistance per length)
%}

%{
Author = Michael Zhang
Date created = 09-06-18
%}
%% Parameters
% if chi is to change, then this would be better suited as an input or
chi = 250;

% similar to chi
phi_f = 0.99; % Fluid phase fraction

%% Constants
% E, M are to do with elastic displacement of the wall, which we are
% ignoring
% [A,aBcs,D1,D2,E,M] = Calculate_Constants(chi,h,phi_f)
[A,aBcs,D1,D2,~,~] = Calculate_Constants(chi,h,phi_f);

%% Flow Velocity

v_f = @(x) (D1*exp(-sqrt(chi).*x)+D2*exp(sqrt(chi).*x)-(A+aBcs)/chi).^2;
v_l = @(x) ((A*x.^2)/2 + 1).^2;

% u = @(x) -D1*exp(-sqrt(chi).*x)-D2*exp(sqrt(chi).*x) + (A*(phi+1)*x.^2)/2 + 1 + E*x+M;

% isPlotting = false
% if isPlotting
%     % non-dimensional radial discretisation for the velocity profile if we are
%     % plotting
%     dr = 0.01;
%     
%     r = [0:dr:h, h:dr:1];
%     v = [v_l(0:dr:h),v_f(h:dr:1)];
%     
%     figure
%     plot([fliplr(-r),r],[fliplr(v),v])
%     xlabel('r')
%     ylabel('v')
% end

%% Flow Rate
% Flowrate as the integral for the 1D graph.
Flowrate = integral(v_f,h,1) + integral(v_l,0,h);


%% Calculate resistance
% Division of the pressure gradient A, with the flow rate
Resistivity = -A/Flowrate;

end

