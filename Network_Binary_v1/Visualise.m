function [] = Visualise(Nodes,Vessels,num_bif,R_0,u)
%% Co-ordinates of the nodes

num_bif = 4; % number of bifurations to the middle. ie expansion.
num_gen = num_bif * 2;  % Number of generations 
num_levels = num_bif * 2 + 1; % Number of levels

gen_vessels = Doubling(num_bif);  % Array, Number of vessels in each generation
lvl_nodes = [1, gen_vessels, 1]; % Array, Number of nodes in each level
lvl_nodes(num_bif + 1) = [];

num_vessels = sum(gen_vessels); % Number of vessels total
num_nodes = sum(lvl_nodes); % Number of nodes total.

gen_vessel_ID = getID(gen_vessels,1);
lvl_node_ID = getID(lvl_nodes,0);


theta_init = pi/4;
theta_end = pi/7;
theta = linspace(theta_init, theta_end, num_bif-1); % Angle of branching at the bifurcation point
length_factor = 5; % length_factor * vessel radius = vessel length
max_linewidth = 10; % Linewidth of the plot representing the maximum vessel
                   % radius in the network
                   
node_marker = (max_linewidth/0.5)*10;
node_txt = 12;    
origin_x = 0;
origin_y = 0;

% Store the (x, y) coordinates of the nodes

% Locating the first vessel and the last vessel at the origin
Nodes{1}.x = 0;
Nodes{1}.y = 0;

Nodes{2}.x = length_factor*Vessels{Nodes{2}.Parent_Vessel}.Radius;
Nodes{2}.y = 0;

Nodes{num_nodes-1}.y = 0;
Nodes{num_nodes}.y = 0;

% X co-ordinates
% Starting from the LHS
for i = 3:1:num_levels
    
    % Nodes in the current level
    nodes_at_lvl = lvl_node_ID{i}+1;
    
    for j = 1:1:length(nodes_at_lvl)
        % x coordinate : x of nodes at a lower level + length of the
        % parent vessel
        Nodes{nodes_at_lvl(j)}.x = Nodes{lvl_node_ID{i-1}(1)+1}.x + ...
            length_factor*...
            Vessels{Nodes{nodes_at_lvl(j)}.Parent_Vessel(1)}.Radius;
    end
end

% Y co-ordinates
% Use the symmetry
k = 1;
for i = 2:1:num_bif-1
    
    % Nodes in the current level
    nodes_at_lvl = lvl_node_ID{i}+1;
    
    % Nodes in the next level
    nodes_at_nxt_lvl = lvl_node_ID{i+1}+1;
    
    % Using the symmetry -> nodes on the LHS
    nodes_sym = lvl_node_ID{num_levels-(i)}+1;
    
    for j = 1:1:length(nodes_at_lvl)
        % Nodes in the next level; top
        Nodes{nodes_at_nxt_lvl((j-1)*2+1)}.y = Nodes{nodes_at_lvl(j)}.y...
            + tan(theta(k))*length_factor*...
           Vessels{Nodes{nodes_at_lvl(j)}.Daughter_Vessel(1)}.Radius;
        Nodes{nodes_sym((j-1)*2+1)}.y = ...
           Nodes{nodes_at_nxt_lvl((j-1)*2+1)}.y;
        
        % Nodes in the next level; bottom
        Nodes{nodes_at_nxt_lvl((j-1)*2+2)}.y = Nodes{nodes_at_lvl(j)}.y...
            - tan(theta(k))*length_factor*...
            Vessels{Nodes{nodes_at_lvl(j)}.Daughter_Vessel(1)}.Radius;
        Nodes{nodes_sym((j-1)*2+2)}.y = ...
           Nodes{nodes_at_nxt_lvl((j-1)*2+2)}.y;
    end
    k = k+1;
end

% Y co-ordinates for the nodes at the centre

nodes_at_lvl = lvl_node_ID{num_bif}+1;

% Nodes in the next level
nodes_at_nxt_lvl = lvl_node_ID{num_bif+1}+1;

for j = 1:1:length(nodes_at_lvl)
    
    % Nodes in the next level; top
    Nodes{nodes_at_nxt_lvl((j-1)*2+1)}.y = Nodes{nodes_at_lvl(j)}.y...
        + tan(theta(k))*length_factor*...
       Vessels{Nodes{nodes_at_lvl(j)}.Daughter_Vessel(1)}.Radius;
   
    % Nodes in the next level; bottom
    Nodes{nodes_at_nxt_lvl((j-1)*2+2)}.y = Nodes{nodes_at_lvl(j)}.y...
        - tan(theta(k))*length_factor*...
        Vessels{Nodes{nodes_at_lvl(j)}.Daughter_Vessel(1)}.Radius;
    
end

%% Draw the struacture of the network
ax1 = subplot(2,1,1);
title_font_sz = 16;
quant_font_sz = 1;
title(ax1, sprintf...
    ('Structure of the network | No. Bifurcations = %i', num_bif), ...
    'FontSize', title_font_sz)
hold(ax1, 'on')
for i = 1:1:num_vessels
    
    x = [Nodes{Vessels{i}.Start_Node+1}.x, Nodes{Vessels{i}.End_Node+1}.x];
    y = [Nodes{Vessels{i}.Start_Node+1}.y, Nodes{Vessels{i}.End_Node+1}.y];
    % Plot the vessels
    linewidth_scale = (1-(R_0 - Vessels{i}.Radius)/R_0)*max_linewidth;
    plot(ax1, x, y, 'linewidth', linewidth_scale, 'color', [0.8, 0.8, 0.8])

end
hold(ax1, 'off')


for i = 1:1:num_vessels
    % Vessel numbering
%     text(0.5*(Nodes{Vessels{i}.Start_Node+1}.x+Nodes{Vessels{i}.End_Node+1}.x),...
%     0.5*(Nodes{Vessels{i}.Start_Node+1}.y+Nodes{Vessels{i}.End_Node+1}.y)...
%     +(1-(R_0 - Vessels{i}.Radius)/R_0)*max_linewidth, sprintf('Vessel %i', i), ...
%         'HorizontalAlignment', 'center', 'color', [0, 0, 0], ...
%         'FontSize' ,node_txt*quant_font_sz, 'VerticalAlignment', 'bottom')
    % Vessel's resistance & Area
    text(0.5*(Nodes{Vessels{i}.Start_Node+1}.x+Nodes{Vessels{i}.End_Node+1}.x),...
    0.5*(Nodes{Vessels{i}.Start_Node+1}.y+Nodes{Vessels{i}.End_Node+1}.y),...
    sprintf('K_{%i} = %.3e\nR_{%i} = %.2f %sm', i, Vessels{i}.Resistance, i, Vessels{i}.Radius, '\mu'), ...
        'HorizontalAlignment', 'center', 'color', [0, 0, 0], ...
        'FontSize' ,node_txt*quant_font_sz, 'VerticalAlignment', 'middle')
end    
    
% Plot the nodes
hold(ax1, 'on')
for i = 1:1:num_nodes
    scatter(ax1, Nodes{i}.x, Nodes{i}.y, node_marker, ...
        'MarkerFacecolor', [1, 0, 0], 'MarkerEdgeColor', 'none')
    % Node numbering
    text(Nodes{i}.x, Nodes{i}.y, sprintf('%i', i-1), ...
        'HorizontalAlignment', 'center', 'color', [1, 1, 1], ...
        'FontSize' ,node_txt, 'FontWeight', 'bold')
end
hold(ax1, 'off')

%% Draw the network with the solved quantities
ax2 = subplot(2,1,2);
title(ax2, 'Solved network', 'FontSize', title_font_sz)

% Plot the flow rate
cmap_Q = flip(colormap('autumn'));
max_Q = max(u(num_nodes+1:end));
min_Q = min(u(num_nodes+1:end));
cmap_Q_ind = linspace(min_Q, max_Q, size(cmap_Q, 1));

hold(ax2, 'on')

for i = 1:1:num_vessels
    
    x = [Nodes{Vessels{i}.Start_Node+1}.x, Nodes{Vessels{i}.End_Node+1}.x];
    y = [Nodes{Vessels{i}.Start_Node+1}.y, Nodes{Vessels{i}.End_Node+1}.y];
    
    linewidth_scale = (1-(R_0 - Vessels{i}.Radius)/R_0)*max_linewidth;
    [~, cmap_i] = min(abs(cmap_Q_ind - u(num_nodes+i)));
    plot(ax2, x, y, 'linewidth', linewidth_scale, 'color', cmap_Q(cmap_i, :))
    text(0.5*(Nodes{Vessels{i}.Start_Node+1}.x+Nodes{Vessels{i}.End_Node+1}.x),...
        0.5*(Nodes{Vessels{i}.Start_Node+1}.y+Nodes{Vessels{i}.End_Node+1}.y),...
        sprintf('Q_{%i}=%.3e', i, u(num_nodes+i)), ...
        'HorizontalAlignment', 'center', 'color', [0, 0, 0], ...
        'BackgroundColor', 'none', 'FontSize' ,node_txt*quant_font_sz, ...
        'VerticalAlignment', 'middle', 'Margin', 1)
end
hold(ax2, 'off')


% Plot the pressure
cmap_P = flip(colormap('winter'));
cmap_P = flip(colormap('summer'));
max_P = max(u(1:num_nodes));
min_P = min(u(1:num_nodes));
cmap_P_ind = linspace(min_P, max_P, size(cmap_P, 1));

hold(ax2, 'on')
for i = 1:1:num_nodes
    [~, cmap_i] = min(abs(cmap_P_ind - u(i)));
    scatter(ax2, Nodes{i}.x, Nodes{i}.y, node_marker, ...
        'MarkerFacecolor', cmap_P(cmap_i, :), 'MarkerEdgeColor', 'none')
    
    text(Nodes{i}.x, Nodes{i}.y , sprintf('P_{%i}=%.2f', i-1, u(i)), ...
        'HorizontalAlignment', 'center', 'color', [0, 0, 0], ...
        'FontSize' ,node_txt*quant_font_sz, 'VerticalAlignment', 'middle')
end
hold(ax2, 'off')

yl = ylim(ax2);
xl = xlim(ax2);
% ylim(ax1, [yl(1)-(yl(2)-yl(1))*0.05, yl(2)+(yl(2)-yl(1))*0.05]);
xlim(ax1, [xl(1)-(xl(2)-xl(1))*0.05, xl(2)+(xl(2)-xl(1))*0.05]);
linkaxes([ax1, ax2], 'xy')

%% Draw the network with the solved quantities (separately)
figure
ax3 = subplot(1,1,1);
title(ax3, 'Solved network', 'FontSize', title_font_sz)

% Plot the flow rate
cmap_Q = flip(colormap('autumn'));
max_Q = max(u(num_nodes+1:end));
min_Q = min(u(num_nodes+1:end));
cmap_Q_ind = linspace(min_Q, max_Q, size(cmap_Q, 1));

hold(ax3, 'on')

for i = 1:1:num_vessels
    
    x = [Nodes{Vessels{i}.Start_Node+1}.x, Nodes{Vessels{i}.End_Node+1}.x];
    y = [Nodes{Vessels{i}.Start_Node+1}.y, Nodes{Vessels{i}.End_Node+1}.y];
    
    linewidth_scale = (1-(R_0 - Vessels{i}.Radius)/R_0)*max_linewidth;
    [~, cmap_i] = min(abs(cmap_Q_ind - u(num_nodes+i)));
    plot(ax3, x, y, 'linewidth', linewidth_scale, 'color', cmap_Q(cmap_i, :))
    text(0.5*(Nodes{Vessels{i}.Start_Node+1}.x+Nodes{Vessels{i}.End_Node+1}.x),...
        0.5*(Nodes{Vessels{i}.Start_Node+1}.y+Nodes{Vessels{i}.End_Node+1}.y),...
        sprintf('Q_{%i}=%.3e', i, u(num_nodes+i)), ...
        'HorizontalAlignment', 'center', 'color', [0, 0, 0], ...
        'BackgroundColor', 'none', 'FontSize' ,node_txt*quant_font_sz, ...
        'VerticalAlignment', 'middle', 'Margin', 1)
end
hold(ax3, 'off')


% Plot the pressure
cmap_P = flip(colormap('winter'));
cmap_P = flip(colormap('summer'));
max_P = max(u(1:num_nodes));
min_P = min(u(1:num_nodes));
cmap_P_ind = linspace(min_P, max_P, size(cmap_P, 1));

hold(ax3, 'on')
for i = 1:1:num_nodes
    [~, cmap_i] = min(abs(cmap_P_ind - u(i)));
    scatter(ax3, Nodes{i}.x, Nodes{i}.y, node_marker, ...
        'MarkerFacecolor', cmap_P(cmap_i, :), 'MarkerEdgeColor', 'none')
    
    text(Nodes{i}.x, Nodes{i}.y , sprintf('P_{%i}=%.2f', i-1, u(i)), ...
        'HorizontalAlignment', 'center', 'color', [0, 0, 0], ...
        'FontSize' ,node_txt*quant_font_sz, 'VerticalAlignment', 'middle')
end
hold(ax3, 'off')

yl = ylim(ax3);
xl = xlim(ax3);
% ylim(ax3, [yl(1)-(yl(2)-yl(1))*0.05, yl(2)+(yl(2)-yl(1))*0.05]);
xlim(ax3, [xl(1)-(xl(2)-xl(1))*0.05, xl(2)+(xl(2)-xl(1))*0.05]);

end