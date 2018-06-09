function [flow_rate, A, resistance] = CalcFlowRate_vary_h(h)
% Calculate a flow rate by integrating the velocity profile based on
% different chi values

%% Parameters
epsilon = 1 - h;
chi = 250;
phi_f = 0.99;
phi_s = 1 - phi_f;
phi = phi_s/phi_f;

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
x_end_l = h;

volume_lumen = pi*(((A.*(x_end_l.^4))/4 + (x_end_l.^2)) - ...
    ((A*(x_start_l^4))/4 + (x_start_l^2)));

% In the EGL
x_start_EGL = x_end_l;
x_end_EGL = 1;

volume_EGL_0 = @(x) -((sqrt(chi)).^(-1)).*D1.*x.*exp(-sqrt(chi)*x);
volume_EGL_1 = @(x) (chi.^(-1)).*D1.*exp(-sqrt(chi).*x);
volume_EGL_2 = @(x) ((sqrt(chi)).^(-1)).*D2.*x.*exp(sqrt(chi).*x);
volume_EGL_3 = @(x) (chi.^(-1)).*D2.*exp(sqrt(chi).*x);
volume_EGL_4 = @(x) (((A+aBcs)./(x.*chi).*x));

volume_EGL = 2*pi*(volume_EGL_0(x_end_EGL) - volume_EGL_0(x_start_EGL) ...
    -(volume_EGL_1(x_end_EGL) - volume_EGL_1(x_start_EGL)) ...
    + volume_EGL_2(x_end_EGL) - volume_EGL_2(x_start_EGL) ...
    -(volume_EGL_3(x_end_EGL) - volume_EGL_3(x_start_EGL)) ...
    -(volume_EGL_4(x_end_EGL) - volume_EGL_4(x_start_EGL)));

flow_rate = volume_lumen + volume_EGL;
resistance = abs(A)./flow_rate;
end