%% Rigid network
% Create a small network.

clear
close all
clc
%% Parameters

% Generation is a "continuous range between two levels".
% Level is a "discrete point"
P0 = 150;          % Inital Pressure
Pinf = 0.5*P0;          % Final Pressure  
R_0 = 50; % Input vessel radius in micronmetres

num_bif = 4; % number of bifurations to the middle. ie expansion.
num_gen = num_bif * 2;  % Number of generations 
num_levels = num_bif * 2 + 1; % Number of levels

gen_vessels = Doubling(num_bif);  % Array, Number of vessels in each generation
lvl_nodes = [1, gen_vessels, 1]; % Array, Number of nodes in each level
lvl_nodes(num_bif + 1) = [];

% Radii
gen_radii = zeros(1,num_gen);
gen_radii(1) = R_0; %% TODO more effecient.
for i = 1:num_bif
    [gen_radii(i+1)] = BifurcateRadius(gen_radii(i));
end
gen_radii(num_bif+1:end) = fliplr(gen_radii(1:num_bif));

%% Create Nodes and Vessels
num_vessels = sum(gen_vessels); % Number of vessels total
num_nodes = sum(lvl_nodes); % Number of nodes total.

gen_vessel_ID = getID(gen_vessels,1);
lvl_node_ID = getID(lvl_nodes,0);

k = 1;
for i = 1:num_gen
    for j = 1:gen_vessels(i)
        Vessels{k}.ID = k;
        Vessels{k}.Generation = i;
        Vessels{k}.Radius = gen_radii(i);
        Vessels{k}.Area = pi*gen_radii(i)^2;
        Vessels{k}.Resistance = K(Vessels{k}.Radius);
%         Vessels{k}.Resistance = CalcResistance(Vessels{k}.Radius*10^(-6)); % Micron -> metres
        
        % On the RHS with respect to the middle
        if (i > num_bif)
            Vessels{k}.Position = 1;
            if i < num_gen
                Vessels{k}.Start_Node = lvl_node_ID{i}(j);
                Vessels{k}.End_Node = lvl_node_ID{i+1}(ceil(j/2));
%                 Vessels{k}.End_Node = lvl_node_ID{i+1}(floor(j*0.4) + 1);
            else
                % Node numbering starts at 0
                Vessels{k}.Start_Node = num_nodes - 2;
                Vessels{k}.End_Node = num_nodes - 1;
            end
        % On the LHS with respect to the middle   
        elseif (i <= num_bif)
            Vessels{k}.Position = -1;    
            if i > 1
%                 Vessels{k}.Start_Node = lvl_node_ID{i}(floor(j*0.4) + 1);
                Vessels{k}.Start_Node = floor(k/2);
                Vessels{k}.End_Node = lvl_node_ID{i+1}(j);
            else
                Vessels{k}.Start_Node = 0;
                Vessels{k}.End_Node = 1;
            end
        end
        k = k + 1;
    end
end

% for i = 1:num_vessels
% disp(Vessels{i})
% end

Nodes{1}.ID = 0;
Nodes{1}.Level = 1;
Nodes{1}.Position = -1;
Nodes{1}.Parent_Vessel = [];
Nodes{1}.Daughter_Vessel = 1;

k = 2;
for i = 1:num_levels-2
    for j = 1:lvl_nodes(i+1)
        Nodes{k}.ID = k - 1;
        Nodes{k}.Level = i + 1;

        if i < num_bif
            Nodes{k}.Position = -1;
            Nodes{k}.Parent_Vessel = gen_vessel_ID{i}(j);
            Nodes{k}.Daughter_Vessel = gen_vessel_ID{i+1}(2 * floor(j*0.5)+1: 2 * floor(j*0.5)+ 2);
        elseif i > num_bif
            Nodes{k}.Position = 1;
            Nodes{k}.Daughter_Vessel = gen_vessel_ID{i+1}(j);
            Nodes{k}.Parent_Vessel = gen_vessel_ID{i}(2 * floor(j*0.5)+1: 2 * floor(j*0.5)+ 2);
        else
            Nodes{k}.Position = 0;
            Nodes{k}.Daughter_Vessel = gen_vessel_ID{i+1}(j);
            Nodes{k}.Parent_Vessel = gen_vessel_ID{i}(j);
        end
        k = k +1;
    end
end

Nodes{num_nodes}.ID = num_nodes;
Nodes{num_nodes}.Level = num_levels;
Nodes{num_nodes}.Position = 1;
Nodes{num_nodes}.Parent_Vessel = num_vessels;
Nodes{num_nodes}.Daughter_Vessel = [];


% for i = 1:num_nodes
% disp(Nodes{i})
% end
       

%% Formulate Matrix
 
AQ = zeros(num_nodes - 2, num_vessels);
AP = zeros(num_vessels, num_nodes + num_vessels);

% Mass Conservation
for i = 2:num_nodes - 1
    a = Nodes{i}.Parent_Vessel;
    b = Nodes{i}.Daughter_Vessel;
    
    AQ(i - 1, a(1)) = 1;
    AQ(i - 1, b(1)) = -1;
    
    if Nodes{i}.Position == -1
        AQ(i - 1, b(2)) = -1;
    elseif Nodes{i}.Position == 1
        AQ(i - 1, a(2)) = 1;
    end
end

clear a b

% Flow equation
for i = 1:num_vessels    
    AP(i,Vessels{i}.End_Node + 1) = 1;
    AP(i,Vessels{i}.Start_Node + 1) = -1;
    AP(i,num_nodes + Vessels{i}.ID) =  K(Vessels{i}.Radius);
%     AP(i,num_nodes + Vessels{i}.ID) =  CalcResistance(Vessels{k}.Radius*10^(-6));
end


% Inital and Final pressure Constraints
AR = zeros(2, num_nodes + num_vessels);
AR(1,1) = 1;
AR(2,num_nodes) = 1;

A = [AP; zeros(num_nodes - 2,num_nodes), AQ; AR];

clear AP AQ AR;
% RHS
b = [zeros(num_nodes + num_vessels - 2,1); P0; Pinf];

% Solve
u = A\b;

%% Assign Flow and Pressure values
for i = 1:num_nodes
    Nodes{i}.Pressure = u(i);
end

for i = 1:num_vessels
   Vessels{i}.Flow = u(num_nodes + i); 
end

clear i j k;

Visualise(Nodes,Vessels,num_bif,R_0,u)
