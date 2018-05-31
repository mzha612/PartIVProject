function [resistance, flow_rate, pressure_gradient] = ...
    CalcResistance_0(H, H_ref, length_factor)
% This function calculates the resistnace for a given vessel radius by
% dividing the flow rate by the pressure drop
% Input:
%   * H = vessel radius in units of meteres
%   * H_ref = reference vessel radius from the network in units of meters
%   * length_factor = length_factor * vessel radius = vessel length
% Output:
%   * resistance = non-dimensional flow rate / (pressure drop gradient
% (assumed to be constant)*length of vessel)
%   * flow_rate = non-dimensional flow rate
%   * A = non-dimensional pressure gradient

%% Parameters


epsilon_true = 1*10^(-6); % Thickness of EGL layer (metres)
epsilon = epsilon_true / H; % [Non-dimensional] thickness of EGL layer
K = 9.9e9; % Hydraulic resistivity (Pascal second per metres squared)
mu_f = 10^(-3); % Plasma dynamic viscosity (Pascal second)
phi_f = 0.99; % Fluid phase fraction
phi_s = 1 - phi_f; % Solud fraction
phi = phi_s/phi_f;
chi = (K*(H^2))/(phi_f*mu_f); % ! [Non-dimensional] Darcy permeability
% Non-mimesionalised with respect to the radius of the given vessel
V = 10^(-3); % Flow velocity (metres per second)
P = mu_f*V/H; % Characteristic Pressure

% [Non-dimensional] radius of lumen to the total radius
h = 1 - epsilon;

%% Compute the constants

Q1 = h;
Q2 = h*chi.*exp(-sqrt(chi));
Q3 = sqrt(chi).*exp(-sqrt(chi)*h);

Q = repelem(Q1, length(chi)) ./ (Q2 + Q3);


T1 = Q.*exp(-sqrt(chi)*h);
T2 = (h*sqrt(chi)/2) + (repelem(phi_f, length(chi))./(h*sqrt(chi)));
T3 = repelem(phi_f, length(chi))./chi;
T4 = phi_f*Q.*exp(-sqrt(chi)*h);

T = (T1 .* T2) - T3 + T4;


R1 = exp(sqrt(chi)*h);
R2 = (sqrt(chi)*h).*exp(sqrt(chi));
R3 = exp(-sqrt(chi)*h);
R4 = (sqrt(chi)*h).*exp(-sqrt(chi));

R = (R1 - R2)./(R3 + R4);


S1 = T2;
S2 = R1 - (R .* R3);
S3 = phi_f * (R1 + (R .* R3));

S = S1 .* S2 - S3;


aBcs1 = (S.*sqrt(chi)).^(-1);
aBcs2 = R .* (exp(-sqrt(chi)) - exp(-sqrt(chi)*h)) - exp(sqrt(chi))...
    + exp(sqrt(chi)*h);
aBcs3 = (1-h)./(S*sqrt(chi) .* h);
aBcs4 = -S2;
aBcs5 = (1-h)/chi;
aBcs6 = sqrt(chi).^(-1);
aBcs7 = (T./S).*(exp(sqrt(chi)) - exp(sqrt(chi)*h));
aBcs8 = ((R.*T)./S + Q).*(exp(-sqrt(chi)) - exp(-sqrt(chi)*h));
aBcs9 = (1-h)./(sqrt(chi)*h);
aBcs10 = (T./S).*exp(sqrt(chi)*h);
aBcs11 = ((R.*T)./S + Q).*exp(-sqrt(chi)*h);

aBcs = (aBcs1.*aBcs2 -aBcs3.*aBcs4)...
    .*((aBcs5 - aBcs6 .* ( aBcs7  - aBcs8)...
    + aBcs9 .* (aBcs10 - aBcs11)).^(-1));


A1 = repelem(sqrt(chi), length(h))./h;
A2 = (T.*exp(sqrt(chi)*h))./S;
A3 = exp(-sqrt(chi)*h).*((R .* T) ./ S + Q);
A4 = -S2./S;

A = A1.*(aBcs .* (A2 - A3) + A4);
A_i = A;

D2 = (aBcs.*T - 1)./S;
D1 = R.*D2 + aBcs .* Q;


E1 = sqrt(chi)*(phi + 1);
E2 = D2 .* exp(sqrt(chi)*h);
E3 = D1 .* exp(-sqrt(chi)*h);
E4 = A.*h*(phi+ 1);

E = E1.*(E2 - E3)-E4;


M1 = D1 .* exp(-sqrt(chi));
M2 = D2 .* exp(sqrt(chi));
M3 = A*(phi+1)/2;

M = M1 + M2 - M3 - 1;




%% Compute the flow rate

% In the lumen
x_start_l = 0;
x_end_l = h*(H/H_ref);

% In the EGL
x_start_EGL = x_end_l;
x_end_EGL = 1*(H/H_ref);

% In the lumen
% 3-D version
% volume_lumen = 2*pi*((0.25*A*((x_end_l*(H_ref/H))^4) + ...
%     ((x_end_l*(H_ref/H))^2)) - ...
%     (0.25*A*((x_start_l*(H_ref/H))^4) + ((x_start_l*(H_ref/H))^2)));

% 2-D version
volume_lumen = @(x) (A*(x^3))/6 + x;
volume_lumen_half = volume_lumen(x_end_l*(H_ref/H)) - volume_lumen(x_start_l*(H_ref/H));

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
volume_EGL_0 = @(x) ((-1)/sqrt(chi))*D1*exp(-sqrt(chi)*x);
volume_EGL_1 = @(x) (1/sqrt(chi))*D2*exp(sqrt(chi)*x);
volume_EGL_2 = @(x) ((A+aBcs)/chi)*x;

volume_EGL_half = (volume_EGL_0(x_end_EGL*(H_ref/H)) + volume_EGL_1(x_end_EGL*(H_ref/H)) ...
    - volume_EGL_2(x_end_EGL)*(H_ref/H)) - ...
    (volume_EGL_0(x_start_EGL*(H_ref/H)) + volume_EGL_1(x_start_EGL*(H_ref/H)) ...
    - volume_EGL_2(x_start_EGL)*(H_ref/H));


flow_rate = 2*(volume_lumen_half + volume_EGL_half);
pressure_gradient = A_i*((H_ref^2)/(H^2));
resistance = abs(pressure_gradient*(H/H_ref)*length_factor)/flow_rate;

%% Check

scaling_factor = H_ref/H;

% chi_i = chi;
% chi = chi_i*scaling_factor^2;

% In the lumen
x_start_l_i = 0;
x_end_l_i = h;

% In the EGL
x_start_EGL_i = x_end_l_i;
x_end_EGL_i = 1;

% In the lumen

% 2-D version
volume_lumen = @(x) (A*(x^3))/6 + x;
volume_lumen_half_i = volume_lumen(x_end_l_i) - volume_lumen(x_start_l_i);

% In the EGL
volume_EGL_0 = @(x) ((-1)/sqrt(chi))*D1*exp(-sqrt(chi)*x);
volume_EGL_1 = @(x) (1/sqrt(chi))*D2*exp(sqrt(chi)*x);
volume_EGL_2 = @(x) ((A+aBcs)/chi)*x;

% 2-D version

volume_EGL_half_i = (volume_EGL_0(x_end_EGL_i) + volume_EGL_1(x_end_EGL_i) ...
    - volume_EGL_2(x_end_EGL_i)) - ...
    (volume_EGL_0(x_start_EGL_i) + volume_EGL_1(x_start_EGL_i) ...
    - volume_EGL_2(x_start_EGL_i));


flow_rate_i = 2*(volume_lumen_half_i + volume_EGL_half_i);
flow_rate = flow_rate_i*(1/scaling_factor);
pressure_gradient_i = A_i;
pressure_gradient = pressure_gradient_i*(scaling_factor^2);

resistance = abs(pressure_gradient*(H/H_ref)*length_factor)/flow_rate;
resistance = abs(pressure_gradient*1)/flow_rate;

end