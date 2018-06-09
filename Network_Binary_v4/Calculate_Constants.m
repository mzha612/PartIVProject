function [A,aBcs,D1,D2,E,M] = Calculate_Constants(chi,h,phi_f)
%% Compute the constants v4

%{
Calculates the constant values from the Appendix A in \autocite{summets_2018}.
For Zero-order flow for the straight walled vessel
%}

%{
Inputs:
    chi         % EGL non-dimensional permeability
    h           % 1 - epsilon the non-dimensional EGL thickness
    phi_f       % Fluid phase fraction
Outputs:
    A,aBcs,D1,D2,E,M
%}

%{
Author = Michael Zhang
Date created = 09-06-18
%}
%% phi
phi_s = 1 - phi_f; % Solid phase fraction
phi = phi_s/phi_f; % Non-dimensional phase fraction?

%% Q
Q1 = h;
Q2 = h*chi*exp(-sqrt(chi));
Q3 = sqrt(chi)*exp(-sqrt(chi)*h);

Q = Q1 / (Q2 + Q3);

%% T
T1 = Q*exp(-sqrt(chi)*h);
T2 = (h*sqrt(chi)/2) + (phi_f/(h*sqrt(chi)));
T3 = phi_f/chi;
T4 = phi_f*Q*exp(-sqrt(chi)*h);

T = (T1 * T2) - T3 + T4;

%% R
R1 = exp(sqrt(chi)*h);
R2 = sqrt(chi)*h*exp(sqrt(chi));
R3 = exp(-sqrt(chi)*h);
R4 = sqrt(chi)* h* exp(-sqrt(chi));

R = (R1 - R2)/(R3 + R4);

%% S
S1 = T2;
S2 = R1 - (R * R3);
S3 = phi_f * (R1 + (R * R3));

S = S1 * S2 - S3;

%% \alphaBc_s
aBcs1 = 1/(S*sqrt(chi));
aBcs2 = R * (exp(-sqrt(chi)) - exp(-sqrt(chi)*h)) - exp(sqrt(chi)) + exp(sqrt(chi)*h);

aBcs3 = (1-h)/(S*sqrt(chi) * h);
aBcs4 = -S2;

aBcs5 = (1-h)/chi;
aBcs6 = 1/sqrt(chi);
aBcs7 = (T/S)*(exp(sqrt(chi)) - exp(sqrt(chi)*h));
aBcs8 = ((R*T)/S + Q)*(exp(-sqrt(chi)) - exp(-sqrt(chi)*h));

aBcs9 = (1-h)/(sqrt(chi)*h);
aBcs10 = (T/S)*exp(sqrt(chi)*h);
aBcs11 = ((R*T)/S + Q)*exp(-sqrt(chi)*h);

aBcs = (aBcs1*aBcs2 -aBcs3*aBcs4)...
    *((aBcs5 - aBcs6 * ( aBcs7  - aBcs8)...
    + aBcs9 * (aBcs10 - aBcs11))^-1);

%% A
% A, the pressure gradient here is calculated such that the maximum
% velocity of the fluid is 1.
A1 = sqrt(chi)/h;
A2 = (T*exp(sqrt(chi)*h))/S;
A3 = exp(-sqrt(chi)*h)*((R * T) / S + Q);
A4 = -S2/S;

A = A1*(aBcs * (A2 - A3) + A4);

%% D1, D2
D2 = (aBcs*T - 1)/S;
D1 = R*D2 + aBcs * Q;

%% E
E1 = sqrt(chi)*(phi + 1);
E2 = D2 * exp(sqrt(chi)*h);
E3 = D1 * exp(-sqrt(chi)*h);
E4 = A*h*(phi+ 1);

E = E1*(E2 - E3)-E4;

%% M
M1 = D1 * exp(-sqrt(chi));
M2 = D2 * exp(sqrt(chi));
M3 = A*(phi+1)/2;

M = M1 + M2 - M3 - 1;
end

