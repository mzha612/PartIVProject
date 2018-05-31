clear
% close all
clc

ref_chi = 250;
chi_values = linspace(ref_chi-50, ref_chi+50, 300);
% chi_values = 250;

[Q, A, resistance] = CalcFlowRate_vary_chi(chi_values);

% cmap = colormap('parula');
cmap = colormap('colorcube');

% Plot flow rate vs chi
figure
ax1 = subplot(3, 1, 1);
hold(ax1, 'on')
plot(ax1, chi_values, Q, 'linewidth', 2.0, ...
    'color',cmap(10, :))
yl = ylim(ax1);
plot(ax1, [ref_chi, ref_chi], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax1, 'Flow rate vs \chi')
xlabel(ax1, '\chi')
ylabel(ax1, 'Flow rate')
grid(ax1, 'on')
hold(ax1, 'off')

% Plot pressure gradient A vs chi
ax2 = subplot(3, 1, 2);
hold(ax2, 'on')
plot(ax2, chi_values, A, 'linewidth', 2.0, ...
    'color', cmap(15, :))
yl = ylim(ax2);
plot(ax2, [ref_chi, ref_chi], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax2, 'Pressure gradient A vs \chi')
xlabel(ax2, '\chi')
ylabel(ax2, 'A')
grid(ax2, 'on')
hold(ax2, 'off')

% Plot resistance vs chi
ax3 = subplot(3, 1, 3);
hold(ax3, 'on')
plot(ax3, chi_values, resistance, 'linewidth', 2.0, ...
    'color',cmap(20, :))
yl = ylim(ax3);
plot(ax3, [ref_chi, ref_chi], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax3, '(Pressure gradient/Flow rate) = Resistance vs \chi')
xlabel(ax3, '\chi')
ylabel(ax3, 'Resistance')
grid(ax3, 'on')
hold(ax3, 'off')

% linkaxes([ax1, ax3], 'x')
set([ax1, ax2, ax3], 'FontSize', 13)
