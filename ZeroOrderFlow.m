%%
% \autocite{Sumets_2018}
% Investigation of the zero oder flow for straight walled vessels.
%%
clear all
close all
clc

%% Parameters
chi = 250;
epsilon = 0.2;
h = 1 - epsilon;
phi_f = 0.99;
phi = 1;


%% Constants

Q1 = h;
Q2 = h*chi*exp(-sqrt(chi));
Q3 = sqrt(chi)*exp(-sqrt(chi)*h);

Q = Q1 / (Q2 + Q3);

T1 = Q*exp(-sqrt(chi)*h);
T2 = (h*sqrt(chi)/2) + (phi_f/(h*sqrt(chi)) );
T3 = phi_f/chi;
T4 = phi_f*Q*exp(-sqrt(chi)*h);

T = (T1 * T2) - T3 + T4;

R1 = exp(sqrt(chi)*h);
R2 = sqrt(chi)*h*exp(sqrt(chi));
R3 = exp(-sqrt(chi)*h);
R4 = sqrt(chi)* h* exp(-sqrt(chi));

R = (R1 - R2)/(R3 + R4);

S1 = T2;
S2 = R1 - (R * R3);
S3 = phi_f * (R1 + (R * R3));

S = S1 * S2 - S3;

aBcs1 = 1/(S*sqrt(chi));
aBcs2 = R * (exp(-sqrt(chi)) - exp(-sqrt(chi)*h)) - exp(sqrt(chi)) + exp(sqrt(chi)*h);
aBcs3 = (1-h)/(S*sqrt(chi) * h);
aBcs4 = -S2;
aBcs5 = (1-h)/chi;
aBcs6 = 1\sqrt(chi);
aBcs7 = T/S*(exp(sqrt(chi)) - exp(sqrt(chi)*h));
aBcs8 = (R*T/S + Q)*(exp(-sqrt(chi)) - exp(-sqrt(chi)*h));
aBcs9 = (1-h)/(sqrt(chi)*h);
aBcs10 = T/S*exp(sqrt(chi)*h);
aBcs11 = (R*T/S + Q)*exp(-sqrt(chi)*h);

aBcs = (aBcs1*aBcs2 -aBcs3*aBcs4)...
    *(aBcs5 - aBcs6 * ( aBcs7  - aBcs8)...
    + aBcs9 * (aBcs10 - aBcs11))^-1;

A1 = sqrt(chi)/h;
A2 = T*exp(sqrt(chi)*h)/S;
A3 = exp(-sqrt(chi)*h)*(R * T / S + Q);
A4 = -S2/S;

A = A1*(aBcs * (A2 - A3) + A4);

D2 = (aBcs*T - 1)/S;
D1 = R*D2 + aBcs * Q;

E1 = sqrt(chi)*(phi + 1);
E2 = D2 * exp(sqrt(chi)*h);
E3 = D1 * exp(-sqrt(chi)*h);
E4 = A*h*(phi+ 1);

E = E1*(E2 - E3)-E4;

M1 = D1 * exp(-sqrt(chi));
M2 = D2 * exp(sqrt(chi));
M3 = A*(phi+1)/2;

M = M1 + M2 - M3 - 1;
%% Governing Equations
% CL1 = 0;
% CF1 = 0;

V_l = @(x2, A) A .* x2.^2 / 2 + 1;
V_f = @(D1,D2,A,aBcs,x2) D1 * exp(-sqrt(chi).*x2) + D2 .* exp(sqrt(chi).*x2) - (A + aBcs)/chi;
U = 1; %% TODO
%% Plots

x2m = linspace(-1,-1+epsilon,201);
x2 = linspace(-1+epsilon,1-epsilon,2001);
x2p = linspace(1-epsilon,1,201);

% Reproduce figure 9 to make sure the above is correct.
% TODO
figure
% plot([x2m, x2, x2p],[V_f(D1,D2,A,aBcs,x2m), V_l(x2,A), V_f(D1,D2,A,aBcs,x2p)])
plot(x2p,V_f(D1,D2,A,aBcs,x2p))
figure
plot(x2m,V_f(D1,D2,A,aBcs,x2m))
figure
plot(x2,V_l(x2, A))

