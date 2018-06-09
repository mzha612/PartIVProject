close all
clear
clc

%% Normalized (subscript n)

H_i = 1*10^(-6)*5;
H_0 = H_i*0.4;
a = H_0/H_i; % Scaling factor

epsilon_true = 1*10^(-6); % Thickness of EGL layer (metres)
epsilon = epsilon_true / H_i; % [Non-dimensional] thickness of EGL layer

K = 9.9e9; % Hydraulic resistivity (Pascal second per metres squared)

mu_f = 10^(-3); % Plasma dynamic viscosity (Pascal second)
phi_f = 0.99; % Fluid phase fraction
phi_s = 1 - phi_f; % Solud fraction
phi = phi_s/phi_f;
chi = (K*(H_i^2))/(phi_f*mu_f); % [Non-dimensional] Darcy permeability
x_n = 1;

% Non-mimesionalised with respect to the radius of the given vessel


% [Non-dimensional] radius of lumen to the total radius
h_n = 1 - epsilon;

Q1 = h_n;
Q2 = h_n*chi.*exp(-sqrt(chi));
Q3 = sqrt(chi).*exp(-sqrt(chi)*h_n);

Q = repelem(Q1, length(chi)) ./ (Q2 + Q3);


T1 = Q.*exp(-sqrt(chi)*h_n);
T2 = (h_n*sqrt(chi)/2) + (repelem(phi_f, length(chi))./(h_n*sqrt(chi)));
T3 = repelem(phi_f, length(chi))./chi;
T4 = phi_f*Q.*exp(-sqrt(chi)*h_n);

T = (T1 .* T2) - T3 + T4;


R1 = exp(sqrt(chi)*h_n);
R2 = (sqrt(chi)*h_n).*exp(sqrt(chi));
R3 = exp(-sqrt(chi)*h_n);
R4 = (sqrt(chi)*h_n).*exp(-sqrt(chi));

R = (R1 - R2)./(R3 + R4);


S1 = T2;
S2 = R1 - (R .* R3);
S3 = phi_f * (R1 + (R .* R3));

S = S1 .* S2 - S3;


aBcs1 = (S.*sqrt(chi)).^(-1);
aBcs2 = R .* (exp(-sqrt(chi)) - exp(-sqrt(chi)*h_n)) - exp(sqrt(chi))...
    + exp(sqrt(chi)*h_n);
aBcs3 = (1-h_n)./(S*sqrt(chi) .* h_n);
aBcs4 = -S2;
aBcs5 = (1-h_n)/chi;
aBcs6 = sqrt(chi).^(-1);
aBcs7 = (T./S).*(exp(sqrt(chi)) - exp(sqrt(chi)*h_n));
aBcs8 = ((R.*T)./S + Q).*(exp(-sqrt(chi)) - exp(-sqrt(chi)*h_n));
aBcs9 = (1-h_n)./(sqrt(chi)*h_n);
aBcs10 = (T./S).*exp(sqrt(chi)*h_n);
aBcs11 = ((R.*T)./S + Q).*exp(-sqrt(chi)*h_n);

aBcs = (aBcs1.*aBcs2 -aBcs3.*aBcs4)...
    .*((aBcs5 - aBcs6 .* ( aBcs7  - aBcs8)...
    + aBcs9 .* (aBcs10 - aBcs11)).^(-1));


A1 = repelem(sqrt(chi), length(h_n))./h_n;
A2 = (T.*exp(sqrt(chi)*h_n))./S;
A3 = exp(-sqrt(chi)*h_n).*((R .* T) ./ S + Q);
A4 = -S2./S;

A = A1.*(aBcs .* (A2 - A3) + A4);

D2 = (aBcs.*T - 1)./S;
D1 = R.*D2 + aBcs .* Q;


E1 = sqrt(chi)*(phi + 1);
E2 = D2 .* exp(sqrt(chi)*h_n);
E3 = D1 .* exp(-sqrt(chi)*h_n);
E4 = A.*h_n*(phi+ 1);

E = E1.*(E2 - E3)-E4;


M1 = D1 .* exp(-sqrt(chi));
M2 = D2 .* exp(sqrt(chi));
M3 = A*(phi+1)/2;

M = M1 + M2 - M3 - 1;

velocity_lumen_n = @(x) 0.5*(A*x.^2)+1;
velocity_EGL_n = @(x) D1*exp(-sqrt(chi)*x) + D2*exp(sqrt(chi)*x)-...
    ((A+aBcs)/chi);

x_lumen_n = linspace(0, h_n, 100);
x_EGL_n = linspace(h_n, x_n, 100);

disp(sprintf('Lumen, Normalized result (x = %.2f): %.2f', x_n, velocity_lumen_n(x_n)))
disp(sprintf('EGL, Normalized result (x = %.2f): %.2f', (h_n+x_n)/2, velocity_EGL_n((h_n+x_n)/2)))
subplot(1,2,1)
hold on
plot(x_lumen_n, velocity_lumen_n(x_lumen_n))
plot(x_EGL_n, velocity_EGL_n(x_EGL_n))
yl = ylim;
% plot([h_n, h_n], [yl(1), yl(2)], 'r--', 'linewidth', 1.5)
title('Normalized')
grid on
hold off

%% Scaled (subscript s)

x_s = x_n/a;
h_s = h_n/a;
A_s = (a^2)*A;

velocity_lumen_s = @(x) 0.5*(A_s*x.^2)+1;

velocity_EGL_s = @(x) (D1*(exp(-sqrt(chi)*x*(a-1))))*...
    exp(-sqrt(chi)*x) + ...
    (D2*exp(sqrt(chi)*x*(a-1)))*exp(sqrt(chi)*x) - ...
    (A+aBcs)/chi;

x_lumen_s = linspace(0, h_s,100);
x_EGL_s = linspace(h_s, x_s,100);

for i = 1:length(x_EGL_s)
    velocity_EGL_s_values(i) = velocity_EGL_s(x_EGL_s(i));
end

disp(sprintf('Lumen, Scaled result (x'' = %.2f): %.2f', x_s, velocity_lumen_s(x_s)))
disp(sprintf('EGL, Scaled result (x = %.2f): %.2f', (h_s+x_s)/2, velocity_EGL_s((h_s+x_s)/2)))
subplot(1,2,2)
hold on
plot(x_lumen_s, velocity_lumen_s(x_lumen_s))
plot(x_EGL_s, velocity_EGL_s_values)
xl = xlim;
yl = ylim;

% plot([h_s, h_s], [yl(1), yl(2)], 'r--', 'linewidth', 1.5)

hold off
grid on
title('Scaled')


