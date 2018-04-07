clear
close all
clc

ref_epsilon = 0.2;
ref_h = 1-ref_epsilon;
epsilon = linspace(0.01, 0.7, 300);
h_values = 1 - epsilon;
% h_values = 250;

[Q, A, resistance] = CalcFlowRate_vary_h(h_values);

% cmap = colormap('parula');
cmap = colormap('colorcube');

% Plot flow rate vs h
figure
ax1 = subplot(3, 1, 1);
hold(ax1, 'on')
plot(ax1, h_values, Q, 'linewidth', 2.0, ...
    'color',cmap(10, :))
yl = ylim(ax1);
plot(ax1, [ref_h, ref_h], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax1, 'Flow rate vs h')
xlabel(ax1, 'h')
ylabel(ax1, 'Flow rate')
grid(ax1, 'on')
hold(ax1, 'off')

% Plot pressure gradient A vs h
ax2 = subplot(3, 1, 2);
hold(ax2, 'on')
plot(ax2, h_values, A, 'linewidth', 2.0, ...
    'color', cmap(15, :))
yl = ylim(ax2);
plot(ax2, [ref_h, ref_h], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax2, 'Pressure gradient A vs h')
xlabel(ax2, 'h')
ylabel(ax2, 'A')
grid(ax2, 'on')
hold(ax2, 'off')

% Plot resistance vs h
ax3 = subplot(3, 1, 3);
hold(ax3, 'on')
plot(ax3, h_values, resistance, 'linewidth', 2.0, ...
    'color',cmap(20, :))
yl = ylim(ax3);
plot(ax3, [ref_h, ref_h], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax3, '(Pressure gradient/Flow rate) = Resistance vs h')
xlabel(ax3, 'h')
ylabel(ax3, 'Resistance')
grid(ax3, 'on')
hold(ax3, 'off')

% linkaxes([ax1, ax3], 'x')
set([ax1, ax2, ax3], 'FontSize', 13)
