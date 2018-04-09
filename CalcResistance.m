function resistance = CalcResistance(H)
% This function calculates the resistnace for a given vessel radius by
% dividing the flow rate by the pressure drop
% Input:
%   * H = vessel radius in units of meteres (i.e. not scaled)
% Output:
%   * resistance = flow rate / pressure drop from inlet to outlet

%% Parameters

epsilon_true = 1*10^(-6); % Thickness of EGL layer (metres)
epsilon = epsilon_true / H; % [Non-dimensional] thickness of EGL layer
K = 9.9e9; % Hydraulic resistivity (Pascal second per metres squared)
mu_f = 10^(-3); % Plasma dynamic viscosity (Pascal second)
phi_f = 0.99; % Fluid phase fraction
phi_s = 1 - phi_f; % Solud fraction
phi = phi_s/phi_f;
chi = (K*(H^2))/(phi_f*mu_f); % [Non-dimensional] EGL permeability
V = 10^(-3); % Flow velocity (metres per second)

% h = (H - epsilon_true) / H;
h = 1 - epsilon; % [Non-dimensional] radius of lumen to the total radius
length_factor = 5; % length_factor * vessel radius = length of the vessel

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
aBcs3 = repelem((1-h), length(chi))./(S.*sqrt(chi) * h);
aBcs4 = -S2;
aBcs5 = repelem((1-h), length(chi))./chi;
aBcs6 = sqrt(chi).^(-1);
aBcs7 = (T./S).*(exp(sqrt(chi)) - exp(sqrt(chi)*h));
aBcs8 = ((R.*T)./S + Q).*(exp(-sqrt(chi)) - exp(-sqrt(chi)*h));
aBcs9 = repelem((1-h), length(chi))./(sqrt(chi)*h);
aBcs10 = (T./S).*exp(sqrt(chi)*h);
aBcs11 = ((R.*T)./S + Q).*exp(-sqrt(chi)*h);

aBcs = (aBcs1.*aBcs2 -aBcs3.*aBcs4)...
    .*((aBcs5 - aBcs6 .* ( aBcs7  - aBcs8)...
    + aBcs9 .* (aBcs10 - aBcs11)).^(-1));


A1 = sqrt(chi)/h;
A2 = (T.*exp(sqrt(chi)*h))./S;
A3 = exp(-sqrt(chi)*h).*((R .* T) ./ S + Q);
A4 = -S2./S;

% Pressure gradient (i.e. pressure change per length)
A = A1.*(aBcs .* (A2 - A3) + A4);

D2 = (aBcs.*T - 1)./S;
D1 = R.*D2 + aBcs .* Q;


E1 = sqrt(chi)*(phi + 1);
E2 = D2 .* exp(sqrt(chi)*h);
E3 = D1 .* exp(-sqrt(chi)*h);
E4 = A*h*(phi+ 1);

E = E1.*(E2 - E3)-E4;


M1 = D1 .* exp(-sqrt(chi));
M2 = D2 .* exp(sqrt(chi));
M3 = A*(phi+1)/2;

M = M1 + M2 - M3 - 1;
%% Compute the flow rate

% Convert back to the dimensional values

% In the lumen
x_start_l = 0;
x_end_l = H - epsilon_true;

volume_lumen = V*pi*((0.25*A*(x_end_l^4) + (x_end_l^2)) - ...
    (0.25*A*(x_start_l^4) + (x_start_l^2)));

% In the EGL
x_start_EGL = x_end_l;
x_end_EGL = H;

volume_EGL_0 = @(x) -((sqrt(chi)).^(-1)).*D1*x.*exp(-sqrt(chi)*x);
volume_EGL_1 = @(x) (chi.^(-1)).*D1.*exp(-sqrt(chi)*x);
volume_EGL_2 = @(x) ((sqrt(chi)).^(-1)).*D2*x.*exp(sqrt(chi)*x);
volume_EGL_3 = @(x) (chi.^(-1)).*D2.*exp(sqrt(chi)*x);
volume_EGL_4 = @(x) (((A+aBcs)./(x*chi)*x));

volume_EGL = V*2*pi*(volume_EGL_0(x_end_EGL) - volume_EGL_0(x_start_EGL) ...
    -(volume_EGL_1(x_end_EGL) - volume_EGL_1(x_start_EGL)) ...
    + volume_EGL_2(x_end_EGL) - volume_EGL_2(x_start_EGL) ...
    -(volume_EGL_3(x_end_EGL) - volume_EGL_3(x_start_EGL)) ...
    -(volume_EGL_4(x_end_EGL) - volume_EGL_4(x_start_EGL)));

flow_rate = volume_lumen + volume_EGL;
resistance = abs(A*length_factor)./flow_rate;

end