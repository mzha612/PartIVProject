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
phi_s = 1 - phi_f;
phi = phi_s/phi_f;


%% Constants

Q1 = h;
Q2 = h*chi*exp(-sqrt(chi));
Q3 = sqrt(chi)*exp(-sqrt(chi)*h);

Q = Q1 / (Q2 + Q3);
%%
T1 = Q*exp(-sqrt(chi)*h);
T2 = (h*sqrt(chi)/2) + (phi_f/(h*sqrt(chi)));
T3 = phi_f/chi;
T4 = phi_f*Q*exp(-sqrt(chi)*h);

T = (T1 * T2) - T3 + T4;
%%
R1 = exp(sqrt(chi)*h);
R2 = sqrt(chi)*h*exp(sqrt(chi));
R3 = exp(-sqrt(chi)*h);
R4 = sqrt(chi)* h* exp(-sqrt(chi));

R = (R1 - R2)/(R3 + R4);

%%
S1 = T2;
S2 = R1 - (R * R3);
S3 = phi_f * (R1 + (R * R3));

S = S1 * S2 - S3;

%%
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

%%
A1 = sqrt(chi)/h;
A2 = (T*exp(sqrt(chi)*h))/S;
A3 = exp(-sqrt(chi)*h)*((R * T) / S + Q);
A4 = -S2/S;

A = A1*(aBcs * (A2 - A3) + A4);

%%
D2 = (aBcs*T - 1)/S;
D1 = R*D2 + aBcs * Q;

%%
E1 = sqrt(chi)*(phi + 1);
E2 = D2 * exp(sqrt(chi)*h);
E3 = D1 * exp(-sqrt(chi)*h);
E4 = A*h*(phi+ 1);

E = E1*(E2 - E3)-E4;

%%
M1 = D1 * exp(-sqrt(chi));
M2 = D2 * exp(sqrt(chi));
M3 = A*(phi+1)/2;

M = M1 + M2 - M3 - 1;
%% Governing Equations
% CL1 = 0;
% CF1 = 0;

% V_l = @(x2, A) A .* x2.^2 / 2 + 1;
% V_f = @(D1,D2,A,aBcs,x2) D1 * exp(-sqrt(chi).*x2) + D2 .* exp(sqrt(chi).*x2) - (A + aBcs)/chi;

V_l = @(x2, A) 0.5*A*(x2.^2) + 1;
V_f = @(D1,D2,A,aBcs,x2) D1 * exp(-sqrt(chi)*x2) + D2 * exp(sqrt(chi)*x2) - (A + aBcs)/chi;
U = @(D1, D2, A, E, M, x2) -D1*exp(-sqrt(chi)*x2) - D2*exp(sqrt(chi)*x2)...
    + ((A*(phi+1)*(x2.^2))/2) + E*x2 + M;

%% Plots

x2_EGL_neg = linspace(-1,-1+epsilon,201);
x2_lumen = linspace(-1+epsilon,1-epsilon,2001);
x2_EGL_pos = linspace(1-epsilon,1,201);

% 1) Reproduce the velocities of the fluid
figure

% cmap = colormap('parula');
cmap = colormap('colorcube');

% Fluid phase in negative coordinates
ax1 = subplot(2, 3, 1);
plot(ax1, x2_EGL_neg, V_f(D1,D2,A,aBcs,x2_EGL_neg), 'linewidth', 2.0, ...
    'color',cmap(10, :))
title(ax1, 'V_f')
xlabel('x_2')
ylabel('V_{0f}^{(0)}')
grid(ax1, 'on')

% Fluid phase in positive coordinates
ax2 = subplot(2, 3, 2);
plot(ax2, x2_EGL_pos, V_f(D1,D2,A,aBcs,x2_EGL_pos), 'linewidth', 2.0, ...
    'color',cmap(10, :))
title(ax2, 'V_f')
xlabel(ax2, 'x_2')
ylabel(ax2, 'V_{0f}^{(0)}')
grid(ax2, 'on')

% Lumen
ax3 = subplot(2, 3, 3);
plot(ax3, x2_lumen,V_l(x2_lumen, A), 'linewidth', 2.0, ...
    'color',cmap(30, :))
title(ax3, 'V_l')
xlabel(ax3, 'x_2')
ylabel(ax3, 'V_{0l}^{(0)}')
grid(ax3, 'on')

% Adjusted to reproduce the figure 9

ax4 = subplot(2, 3, [5,6]);
hold(ax4, 'on')
plot(ax4, x2_EGL_pos, V_f(D1,D2,A,aBcs,x2_EGL_pos), 'linewidth', 2.0, ...
    'color',cmap(10, :))
plot(ax4, x2_lumen, V_l(x2_lumen, A), 'linewidth', 2.0, ...
    'color',cmap(30, :))
plot(ax4, fliplr(x2_EGL_neg), V_f(D1,D2,A,aBcs,x2_EGL_pos), 'linewidth', 2.0, ...
    'color',cmap(10, :))
hold(ax4, 'off')
title(ax4, 'Adjusted velocity values to reproduce the Figure 9')
xlabel(ax4, 'x_2')
ylabel(ax4, 'V_{0\beta}^{(0)}')
grid(ax4, 'on')

linkaxes([ax3, ax2], 'xy')
ylim(ax3, [-0.1, 1.0])

set([ax1, ax2, ax3, ax4], 'FontSize', 15)

%% Plot U
figure
ax1 = subplot(1, 2, 1);
plot(ax1, x2_EGL_neg, U(D1, D2, A, E, M, x2_EGL_neg), 'linewidth', 2.0, ...
    'color',cmap(20, :))
title(ax1, 'Plotted for negative x_2 values')
xlabel('x_2')
ylabel('U_{0}^{(0)}')
grid(ax1, 'on')

ax2 = subplot(1, 2, 2);
plot(ax2, x2_EGL_pos, U(D1, D2, A, E, M, x2_EGL_pos), 'linewidth', 2.0, ...
    'color',cmap(40, :))
title(ax2, 'Plotted for positive x_2 values')
xlabel('x_2')
ylabel('U_{0}^{(0)}')
grid(ax2, 'on')
set([ax1, ax2], 'FontSize', 15)