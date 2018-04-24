clear
close all
clc

%% Vary h
ref_epsilon = 0.2;
ref_h = 1-ref_epsilon;
epsilon = linspace(0, 0.7, 300);
h_values = linspace(10*10^(-6), 40*10^(-6), 1000);
h_values = 10*10^(-6):5*10^(-8):40*10^(-6);
ref_chi = 250;
H_ref = 5*10^(-6);
H_ref = h_values(floor(0.5*length(h_values)));
H_ref_ind = find(H_ref==h_values);
resistance = zeros(1, length(h_values));
Q = zeros(1, length(h_values));
A = zeros(1, length(h_values));
length_factor = 10;
for i = 1:1:length(resistance)
    [resistance(i), Q(i), A(i)] = CalcResistance(h_values(i), H_ref, length_factor);
end


% cmap = colormap('parula');
cmap = colormap('colorcube');

% Plot flow rate vs h
figure
ax1 = subplot(2, 6, 1:2);
hold(ax1, 'on')
plot(ax1, h_values, Q, 'linewidth', 2.0, ...
    'color',cmap(10, :))
plot(ax1, h_values(H_ref_ind), Q(H_ref_ind), 'linewidth', 2.0, ...
    'color', 'red', 'linestyle', 'none', 'marker', 'x', 'markersize', 8)
yl = ylim(ax1);
plot(ax1, [h_values(H_ref_ind), h_values(H_ref_ind)], [yl(1), Q(H_ref_ind)], 'linewidth', 1.0, ...
    'color', 'red', 'linestyle', ':')
title(ax1, 'Flow rate vs Vessel radius')
xlabel(ax1, 'Vessel radius (m)')
ylabel(ax1, 'Flow rate Q')
grid(ax1, 'on')
legend(ax1, '', sprintf('H_{ref} = %.2e m', H_ref), 'location', 'best')
hold(ax1, 'off')

% Plot pressure gradient A vs h
ax2 = subplot(2, 6, 3:4);
hold(ax2, 'on')
plot(ax2, h_values, A, 'linewidth', 2.0, ...
    'color', cmap(15, :))
plot(ax2, h_values(H_ref_ind), A(H_ref_ind), 'linewidth', 2.0, ...
    'color', 'red', 'linestyle', 'none', 'marker', 'x', 'markersize', 8)
yl = ylim(ax2);
plot(ax2, [h_values(H_ref_ind), h_values(H_ref_ind)], [yl(1), A(H_ref_ind)], 'linewidth', 1.0, ...
    'color', 'red', 'linestyle', ':')
title(ax2, 'Pressure gradient A vs Vessel radius')
xlabel(ax2, 'Vessel radius (m)')
ylabel(ax2, 'A')
grid(ax2, 'on')
legend(ax2, '', sprintf('H_{ref} = %.2e m', H_ref), 'location', 'best')
hold(ax2, 'off')

% Plot resistance vs h
ax3 = subplot(2, 6, 5:6);
hold(ax3, 'on')
a = 50;
scatter(ax3, h_values, resistance, a, 'x', 'linewidth', 1.0, ...
    'MarkerEdgeColor', cmap(20, :))
plot(ax3, h_values(H_ref_ind), resistance(H_ref_ind), 'linewidth', 2.0, ...
    'color', 'red', 'linestyle', 'none', 'marker', 'x', 'markersize', 8)
yl = ylim(ax3);
plot(ax3, [h_values(H_ref_ind), h_values(H_ref_ind)], [yl(1), resistance(H_ref_ind)], 'linewidth', 1.0, ...
    'color', 'red', 'linestyle', ':')
title(ax3, 'abs(Pressure gradient)/Flow rate = Resistance vs Vessel radius')
xlabel(ax3, 'Vessel radius (m)')
ylabel(ax3, 'Resistance')
grid(ax3, 'on')
legend(ax3, '', sprintf('H_{ref} = %.2e m', H_ref), 'location', 'best')
hold(ax3, 'off')

% Plot log(resistance) vs h
ax4 = subplot(2, 6, 7:9);
hold(ax4, 'on')
a = 50;
scatter(ax4, h_values, log(resistance), a, 'x', 'linewidth', 1.0, ...
    'MarkerEdgeColor', cmap(20, :))
y = log(resistance);
plot(ax4, h_values(H_ref_ind), y(H_ref_ind), 'linewidth', 2.0, ...
    'color', 'red', 'linestyle', 'none', 'marker', 'x', 'markersize', 8)
yl = ylim(ax4);
plot(ax4, [h_values(H_ref_ind), h_values(H_ref_ind)], [yl(1), y(H_ref_ind)], 'linewidth', 1.0, ...
    'color', 'red', 'linestyle', ':')
title(ax4, 'log(Resistance) vs Vessel radius')
xlabel(ax4, 'Vessel radius (m)')
ylabel(ax4, 'log(Resistance)')
grid(ax4, 'on')
legend(ax4, '', sprintf('H_{ref} = %.2e m', H_ref), 'location', 'best')
hold(ax4, 'off')


%% Fitting a curve to resistance vs h

% linear to log(resistance) vs h
log_R_NaN = logical(isnan(log(resistance)));
fit_y = log(resistance(~log_R_NaN));
fit_x = h_values(~log_R_NaN);
[linear_log, linear_log_S] = polyfit(fit_x, fit_y, 1);

% Polynomial of order 2
R_NaN = logical(isnan(resistance));
fit_y = resistance(~R_NaN);
fit_x = h_values(~R_NaN);
[poly_2, poly_2_S] = polyfit(fit_x, fit_y, 2);

% Polynomial of order 3
[poly_3, poly_3_S] = polyfit(fit_x, fit_y, 3);


ax5 = subplot(2, 6, 10:12);
cmap2 = colormap('lines');
hold(ax5, 'on')
a = 50;
scatter(ax5, h_values, resistance, a, 'x', 'linewidth', 1.0, ...
    'MarkerEdgeColor', cmap(20, :))
plot(ax5, h_values, exp(polyval(linear_log, h_values)), 'linewidth', 2.0, ...
    'color', cmap2(1, :))
plot(ax5, h_values, polyval(poly_2, h_values), 'linewidth', 2.0, ...
    'color', cmap2(2, :))
plot(ax5, h_values, polyval(poly_3, h_values), 'linewidth', 2.0, ...
    'color', cmap2(3, :))
% yl = ylim(ax5);
% plot(ax5, [ref_h, ref_h], [yl(1), yl(2)], 'linestyle', '--', 'linewidth', 1.5, ...
%     'color', 'r')
% title(ax5, sprintf...
%     ('Resistance vs h | chi = %.1f',...
%     ref_chi))
title(ax5, 'Resistance vs vessel radius')
xlabel(ax5, 'Vessel radius (m)')
ylabel(ax5, 'Resistance')
grid(ax5, 'on')
hold(ax5, 'off')
legend(ax5, 'From analytical expression', ...
    sprintf('Exponential fit | norm(Residuals) = %.4f', linear_log_S.normr),...
    sprintf('Quadratic fit | norm(Residuals) = %.4f', poly_2_S.normr),...
    sprintf('Cubic fit | norm(Residuals) = %.4f', poly_3_S.normr))

% linkaxes([ax1, ax4], 'y')
% linkaxes([ax2, ax5], 'y')
% linkaxes([ax3, ax6], 'y')
% set([ax1, ax2, ax3, ax4, ax5, ax6], 'FontSize', 13)

