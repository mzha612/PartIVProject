function [Resistivity] = Resistance_Summets(h,chi)
%% Resistance_Summets v5
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
    chi = dimensionless parameter.
Outputs:
   	Resistivity, Non-dimensional resisitivity (resistance per length)
%}

%{
Author = Michael Zhang
Date created = 09-06-18
%}
%% Parameters
phi_f = 0.99; % Fluid phase fraction

%% Constants
% E, M are to do with elastic displacement of the wall, which we are
% ignoring
[A,aBcs,D1,D2,~,~] = Calculate_Constants(chi,h,phi_f);

%% Flow Velocity
% v_f = @(x) (D1*exp(-sqrt(chi).*x)+D2*exp(sqrt(chi).*x)-(A+aBcs)/chi);
% v_l = @(x) ((A*x.^2)/2 + 1);

%% Flow Rate
% Flowrate as the integral for the 1D graph.

int_v_l = A*h^3/6+h;
int_v_f = (-D1/sqrt(chi)*exp(-sqrt(chi))+D2/sqrt(chi)*exp(sqrt(chi))-(A+aBcs)/chi)...
    - (-D1/sqrt(chi)*exp(-sqrt(chi).*h)+D2/sqrt(chi)*exp(sqrt(chi).*h)-(A+aBcs)*h/chi);

Flowrate = 2*(int_v_f + int_v_l);
% Flowrate = 2*(integral(v_f,h,1) + integral(v_l,0,h));

%% Calculate resistance
% Division of the pressure gradient A, with the flow rate
Resistivity = -A/Flowrate;

end

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