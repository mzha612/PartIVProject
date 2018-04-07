clear
close all
clc

%% Vary h
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
ax1 = subplot(3, 2, 1);
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
ax2 = subplot(3, 2, 3);
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
ax3 = subplot(3, 2, 5);
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

%% Vary chi


ref_chi = 250;
chi_values = linspace(ref_chi-50, ref_chi+50, 300);
% chi_values = 250;

[Q, A, resistance] = CalcFlowRate_vary_chi(chi_values);

% cmap = colormap('parula');
cmap = colormap('colorcube');

% Plot flow rate vs chi
ax4 = subplot(3, 2, 2);
hold(ax4, 'on')
plot(ax4, chi_values, Q, 'linewidth', 2.0, ...
    'color',cmap(10, :))
yl = ylim(ax4);
plot(ax4, [ref_chi, ref_chi], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax4, 'Flow rate vs \chi')
xlabel(ax4, '\chi')
ylabel(ax4, 'Flow rate')
grid(ax4, 'on')
hold(ax4, 'off')

% Plot pressure gradient A vs chi
ax5 = subplot(3, 2, 4);
hold(ax5, 'on')
plot(ax5, chi_values, A, 'linewidth', 2.0, ...
    'color', cmap(15, :))
yl = ylim(ax5);
plot(ax5, [ref_chi, ref_chi], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax5, 'Pressure gradient A vs \chi')
xlabel(ax5, '\chi')
ylabel(ax5, 'A')
grid(ax5, 'on')
hold(ax5, 'off')

% Plot resistance vs chi
ax6 = subplot(3, 2, 6);
hold(ax6, 'on')
plot(ax6, chi_values, resistance, 'linewidth', 2.0, ...
    'color',cmap(20, :))
yl = ylim(ax6);
plot(ax6, [ref_chi, ref_chi], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
    'color', 'r')
title(ax6, '(Pressure gradient/Flow rate) = Resistance vs \chi')
xlabel(ax6, '\chi')
ylabel(ax6, 'Resistance')
grid(ax6, 'on')
hold(ax6, 'off')

% linkaxes([ax1, ax4], 'y')
% linkaxes([ax2, ax5], 'y')
% linkaxes([ax3, ax6], 'y')
set([ax1, ax2, ax3, ax4, ax5, ax6], 'FontSize', 13)

