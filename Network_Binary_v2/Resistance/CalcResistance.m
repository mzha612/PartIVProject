function [resistance, flow_rate, pressure_gradient] = ...
    CalcResistance(H, H_ref)
% This function calculates the resistivity for a given vessel radius by
% dividing the flow rate by the pressure drop gradient

% Input:
%   * H = vessel radius in units of meteres
%   * H_ref = reference vessel radius from the network in units of meters

% Output:
%   * resistivity = non-dimensional flow rate / pressure drop gradient
%     (assumed to be constant)
%   * flow_rate = non-dimensional, scaled flow rate
%   * A = non-dimensional, scaled pressure gradient

%% Parameters


epsilon_true = 1*10^(-6); % Thickness of EGL layer (metres)
epsilon = epsilon_true / H; % [Non-dimensional] thickness of EGL layer
K = 9.9e9; % Hydraulic resistivity (Pascal second per metres squared)
mu_f = 10^(-3); % Plasma dynamic viscosity (Pascal second)
phi_f = 0.99; % Fluid phase fraction
phi_s = 1 - phi_f; % Solud fraction
phi = phi_s/phi_f;
chi = (K*(H^2))/(phi_f*mu_f); % ! [Non-dimensional] Darcy permeability
% Non-dimesionalised, normalized with respect to the radius of the given vessel
V = 10^(-3); % Flow velocity (metres per second)
P = mu_f*V/H; % Characteristic Pressure

% Non-dimensional, normalized radius of lumen to the total radius
h = 1 - epsilon;




%% Compute the flow rate

a = H_ref/H; % Scaling factor

% In the lumen
x_start_l = 0;
x_end_l = h/a;

% In the EGL
x_start_EGL = x_end_l;
x_end_EGL = 1/a;

% In the lumen
% 3-D version
% volume_lumen = 2*pi*((0.25*A*((x_end_l*(H_ref/H))^4) + ...
%     ((x_end_l*(H_ref/H))^2)) - ...
%     (0.25*A*((x_start_l*(H_ref/H))^4) + ((x_start_l*(H_ref/H))^2)));

% Scale the parameters
A_s = (a^2)*A;

% 2-D version
volume_lumen = @(x) (A_s*(x^3))/6 + x;
volume_lumen_half = volume_lumen(x_end_l) - volume_lumen(x_start_l);

% In the EGL
% 3-D version
% volume_EGL_0 = @(x) -((sqrt(chi)).^(-1)).*D1*(x*(H_ref/H)).*exp...
%     (-sqrt(chi)*(x*(H_ref/H)));
% volume_EGL_1 = @(x) (chi.^(-1)).*D1.*exp(-sqrt(chi)*(x*(H_ref/H)));
% volume_EGL_2 = @(x) ((sqrt(chi)).^(-1)).*D2*(x*(H_ref/H)).*exp...
%     (sqrt(chi)*(x*(H_ref/H)));
% volume_EGL_3 = @(x) (chi.^(-1)).*D2.*exp(sqrt(chi)*(x*(H_ref/H)));
% volume_EGL_4 = @(x) (((A+aBcs)./((x*(H_ref/H))*chi)*(x*(H_ref/H))));

% volume_EGL = 2*pi*(volume_EGL_0(x_end_EGL) - volume_EGL_0(x_start_EGL) ...
%     -(volume_EGL_1(x_end_EGL) - volume_EGL_1(x_start_EGL)) ...
%     + volume_EGL_2(x_end_EGL) - volume_EGL_2(x_start_EGL) ...
%     -(volume_EGL_3(x_end_EGL) - volume_EGL_3(x_start_EGL)) ...
%     -(volume_EGL_4(x_end_EGL) - volume_EGL_4(x_start_EGL)));

% 2-D version
volume_EGL_0 = @(x) ((-1)/sqrt(chi))*(D1*(exp(-sqrt(chi)*x*(a-1))))*exp(-sqrt(chi)*x);
volume_EGL_1 = @(x) (1/sqrt(chi))*(D2*exp(sqrt(chi)*x*(a-1)))*exp(sqrt(chi)*x);
volume_EGL_2 = @(x) ((A+aBcs)/chi)*x;

volume_EGL_half = (volume_EGL_0(x_end_EGL) + volume_EGL_1(x_end_EGL) ...
    - volume_EGL_2(x_end_EGL)) - ...
    (volume_EGL_0(x_start_EGL) + volume_EGL_1(x_start_EGL) ...
    - volume_EGL_2(x_start_EGL));


flow_rate = 2*(volume_lumen_half + volume_EGL_half);
pressure_gradient = A_i*((H_ref^2)/(H^2));
resistance = abs(pressure_gradient*(H/H_ref)*length_factor)/flow_rate;


end