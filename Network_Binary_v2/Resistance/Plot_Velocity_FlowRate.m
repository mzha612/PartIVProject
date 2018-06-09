%% This is a script to check how the velocity profile is dependent on the
% pressure gradient
% Positive pressure gradient (arbitary) -> Flows from outlet to inlet
% Negative pressure gradient -> Flows from inlet to oulet

clear
close all
clc

%% Velocity profile w.r.t the current vessel radius

H = 10*10^(-6);
H_ref = 5*10^(-6);

[x_l_i, velocity_l_i, x_f_i, velocity_f_i, A_i, ...
    x_l, velocity_l, x_f, velocity_f, A] = CalcVelocityProfile(H, H_ref);


% cmap = colormap('parula');
cmap = colormap('colorcube');

% Plot velocity profile & dispaly flow rate value
ax1 = subplot(1, 2, 1);
hold(ax1, 'on')
plot(ax1, x_l_i, velocity_l_i, 'linewidth', 2.0, 'color', cmap(10, :))
plot(ax1, x_f_i, velocity_f_i, 'linewidth', 2.0, 'color', cmap(20, :))
% title(ax1, sprintf('Velocity profile at A = %.3f | Flow rate = %.3f', ...
%     A, flow_rate))
title(ax1, sprintf('Velocity profile at A = %.3f', A_i))
xlabel(ax1, 'x_2_i')
ylabel(ax1, '{V_{0\beta}}^{(0)}')
grid(ax1, 'on')
hold(ax1, 'off')
legend(ax1, '{V_l}^{(0)}', 'V_f^{(0)}')

%% Velocity profile w.r.t the reference vessel radius

% Plot velocity profile & dispaly flow rate value

ax2 = subplot(1, 2, 2);
hold(ax2, 'on')
plot(ax2, x_l, velocity_l, 'linewidth', 2.0, 'color', cmap(10, :))
plot(ax2, x_f, velocity_f, 'linewidth', 2.0, 'color', cmap(20, :))
% title(ax1, sprintf('Velocity profile at A = %.3f | Flow rate = %.3f', ...
%     A, flow_rate))
title(ax2, sprintf('Velocity profile at A = %.3f', A))
xlabel(ax2, 'x_2')
ylabel(ax2, '{V_{0\beta}}^{(0)}')
grid(ax2, 'on')
hold(ax2, 'off')
legend(ax2, '{V_l}^{(0)}', 'V_f^{(0)}')