function [x_l_i, velocity_l_i, x_f_i, velocity_f_i, A_i...
    x_l, velocity_l, x_f, velocity_f, A] = ...
    CalcVelocityProfile(H, H_ref)
% This function calculates a velocity profile and a flow rate using
% the referenced values in the lieterature, but arbitarily sets the
% pressure gradient A positive if specified in the input to check that
% the equations are accounting for the pressure gradient

%% Parameters
% epsilon = 0.2;
% h = 1 - epsilon;
% chi = 250;
% phi_f = 0.99;
% phi_s = 1 - phi_f;
% phi = phi_s/phi_f;


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

%% Compute the velocity profile w.r.t. the current vessel radius

% In the lumen
x_start_l = 0;
x_end_l = h;
x_l_i = linspace(x_start_l, x_end_l, 300);
velocity_l_i = (A*(x_l_i.^2))/2 + 1;

% In the EGL
x_start_EGL = x_end_l;
x_end_EGL = 1;
x_f_i = linspace(x_start_EGL, x_end_EGL, 300);
velocity_f_i = D1*exp(-sqrt(chi)*x_f_i) + D2*exp(sqrt(chi)*x_f_i) - ...
    (A + aBcs)/chi;

%% Compute the velocity profile w.r.t. the reference vessel radius

% In the lumen
x_start_l = 0;
x_end_l = h*(H/H_ref);
x_l = linspace(x_start_l, x_end_l, 300);
velocity_l = (A*((x_l*(H_ref/H)).^2))/2 + 1;

% In the EGL
x_start_EGL = x_end_l;
x_end_EGL = 1*(H/H_ref);
x_f = linspace(x_start_EGL, x_end_EGL, 300);
velocity_f = D1*exp(-sqrt(chi)*x_f*(H_ref/H)) + D2*exp(sqrt(chi)*x_f*(H_ref/H)) - ...
    (A + aBcs)/chi;

A = A_i*((H_ref^2)/(H^2));

end