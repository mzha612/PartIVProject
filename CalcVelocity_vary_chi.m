function [V_lumen, V_fluid, U] = CalcVelocity_vary_chi...
    (chi, x2_lumen, x2_EGL_pos)

%% Parameters
epsilon = 0.2;
h = 1 - epsilon;
phi_f = 0.99;
phi_s = 1 - phi_f;
phi = phi_s/phi_f;

%% Compute the constants

Q1 = h;
Q2 = h*chi.*exp(-sqrt(chi));
Q3 = sqrt(chi).*exp(-sqrt(chi)*h);

Q = repelem(Q1, length(chi)) / (Q2 + Q3);


T1 = Q.*exp(-sqrt(chi)*h);
T2 = (h*sqrt(chi)/2) + (repelem(phi_f, length(chi))/(h*sqrt(chi)));
T3 = repelem(phi_f, length(chi))/chi;
T4 = phi_f*Q*exp(-sqrt(chi)*h);

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
aBcs8 = ((R.*T)/S + Q).*(exp(-sqrt(chi)) - exp(-sqrt(chi)*h));
aBcs9 = repelem((1-h), length(chi))./(sqrt(chi)*h);
aBcs10 = (T./S).*exp(sqrt(chi)*h);
aBcs11 = ((R.*T)/S + Q).*exp(-sqrt(chi)*h);

aBcs = (aBcs1.*aBcs2 -aBcs3.*aBcs4)...
    .*((aBcs5 - aBcs6 .* ( aBcs7  - aBcs8)...
    + aBcs9 .* (aBcs10 - aBcs11)).^(-1));


A1 = sqrt(chi)/h;
A2 = (T.*exp(sqrt(chi)*h))/S;
A3 = exp(-sqrt(chi)*h).*((R .* T) / S + Q);
A4 = -S2./S;

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


%% Compute the velocities



V_lumen = (0.5*repmat(A', 1, length(x2_lumen))).* ...
    repmat(x2_lumen.^2, length(chi), 1) + 1;

V_fluid = zeros(length(chi), length(x2_EGL_pos));
U = zeros(length(chi), length(x2_EGL_pos));
for ind = 1:1:length(chi)
    V_fluid(ind, :) = D1(ind) .* exp(-sqrt(chi(ind))*x2_EGL_pos) ...
        + D2(ind) .*  exp(sqrt(chi(ind))*x2_EGL_pos) - (A(ind) + ...
        aBcs(ind))/chi(ind);
    
    U(ind, :) = -D1(ind)*exp(-sqrt(chi(ind))*x2_EGL_pos) - ...
        D2(ind)*exp(sqrt(chi(ind))*x2_EGL_pos)...
    + ((A(ind)*(phi+1)*(x2_EGL_pos.^2))/2) + E(ind)*x2_EGL_pos + M(ind);
end


end

