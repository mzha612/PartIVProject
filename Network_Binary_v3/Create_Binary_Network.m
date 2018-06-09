function [Vessels, Nodes] = Create_Binary_Network(num_bif,d_R_0)
%% Create_Binary_Network v3
%{
Function Creates the vessel segment infomation to be able to form a binary
tree network.
%}

%{
Inputs:
    num_bif % number of bifurations to the middle. ie expansion.
    d_R_0   % input vessel radius in micronmetres
    
Outputs:
    Vessels: Structure, vessel segment infomation
    Nodes: Structure, nodal infomation
%}

%{
Author = Michael Zhang
Date created = 02-06-18

v3 Derived from Networks_v1, and v2.
%}

%% Parameters

% Generation is a "continuous range between two levels".
% Level is a "discrete point"

global num_vessels num_nodes
art2ven = 1;    % Length ratio - venulouse hold 4times more vblood than artrioles./viens/artiers
num_gen = num_bif * 2;              % Number of generations

gen_vessels = Doubling(num_bif);    % Array, Number of vessels in each generation
lvl_nodes = [1, gen_vessels, 1];    % Array, Number of nodes in each level
lvl_nodes(num_bif + 1) = [];

num_vessels = sum(gen_vessels);         % Number of vessels total
gen_vessel_ID = getID(gen_vessels,1);   % ID for each vessel in each generation
lvl_node_ID = getIDnode(lvl_nodes,0);

% ID for each node in each level
lvl_node_ID{1} = 1;
lvl_node_ID{end} = num_nodes;
gen_vessel_ID = [0, gen_vessel_ID, 0];
num_nodes = sum(lvl_nodes);

% Creates the radii for each generation
gen_radii = zeros(1,num_gen);
gen_radii(1) = d_R_0; %% TODO more effecient.
for i = 1:num_bif-1
    [gen_radii(i+1)] = BifurcateRadius(gen_radii(i));
end
gen_radii(num_bif+1:end) = art2ven * fliplr(gen_radii(1:num_bif));

%% Create Nodes and Vessels

Vessels = cell(1,num_vessels); % Holds cell infomation
Nodes = cell(1,num_nodes);     % Holds nodal infomation

k = 1;
for i = 1:num_gen
    for j = 1:gen_vessels(i)
        Vessels{k}.ID = gen_vessel_ID{i+1}(j);
        Vessels{k}.Generation = i;
        
        Vessels{k}.Radius = gen_radii(i);
        Vessels{k}.Length = [];
        Vessels{k}.n_Radius = 1;
        
        %         Vessels{k}.n_Area = gen_radii(i)^2;

        %         Vessels{k}.Resistance = K(Vessels{k}.Radius);
        %         Vessels{k}.Resistance = CalcResistance(Vessels{k}.Radius*10^(-6)); % Micron -> metres
        %         Vessels{k}.BC = 0;
        %         if i <= (num_gen)/2 - 1         % Left side
        %             Vessels{k}.Parent_Node = lvl_node_ID{i}(ceil(j/2));
        %             Vessels{k}.Daughter_Node = lvl_node_ID{i+2}(2*j-1:2*j);
        %         elseif i > (num_gen )/2 + 1     % Right side
        %             Vessels{k}.Parent_Node= lvl_node_ID{i}(2*j-1:2*j);
        %             Vessels{k}.Daughter_Node = lvl_node_ID{i+2}(ceil(j/2));
        %         elseif i == (num_gen)/2         % Left side of middle
        %             Vessels{k}.Parent_Node = lvl_node_ID{i}(ceil(j/2));
        %             Vessels{k}.Daughter_Node = lvl_node_ID{i+2}(j);
        %         else                            % Right side of middle
        %             Vessels{k}.Parent_Node = lvl_node_ID{i}(j);
        %             Vessels{k}.Daughter_Node = lvl_node_ID{i+2}(ceil(j/2));
        %         end
        if i <= num_gen/2
            Vessels{k}.Parent_Node = lvl_node_ID{i}(ceil(j/2));
            Vessels{k}.Daughter_Node = lvl_node_ID{i+1}(j);
            Vessels{k}.Length = getLength(Vessels{k}.Radius);
        else
            Vessels{k}.Parent_Node = lvl_node_ID{i}(j);
            Vessels{k}.Daughter_Node = lvl_node_ID{i+1}(ceil(j/2));
            Vessels{k}.Length = 14.54*Vessels{k}.Radius;
        end
        Vessels{k}.n_Length = Vessels{k}.Length / Vessels{k}.Radius; 
        k = k + 1;
    end
end

for i = 1:num_nodes
    Nodes{i}.ID = i;
    Nodes{i}.Parent_Vessel=[];
    Nodes{i}.Daughter_Vessel=[];
    Nodes{i}.Parent_Node = [];
    Nodes{i}.Daughter_Node = [];
    Nodes{i}.num_Daughter_Nodes = 0;
    Nodes{i}.num_Parent_Nodes = 0;
end


for i = 1:num_vessels
    Nodes{Vessels{i}.Parent_Node}.Daughter_Vessel = [Nodes{Vessels{i}.Parent_Node}.Daughter_Vessel ,Vessels{i}.ID];
    Nodes{Vessels{i}.Daughter_Node}.Parent_Vessel = [Nodes{Vessels{i}.Daughter_Node}.Parent_Vessel, Vessels{i}.ID];
end

for i = 1:num_nodes
    for j = 1: num_vessels
        if Vessels{j}.Parent_Node == i
            Nodes{i}.Daughter_Node = [Nodes{i}.Daughter_Node ,Vessels{j}.Daughter_Node];
            Nodes{i}.num_Daughter_Nodes = Nodes{i}.num_Daughter_Nodes + 1; 
        end
        if Vessels{j}.Daughter_Node== i
            Nodes{i}.Parent_Node  = [Nodes{i}.Parent_Node ,Vessels{j}.Parent_Node ];
            Nodes{i}.num_Parent_Nodes = Nodes{i}.num_Parent_Nodes + 1; 
        end
    end
end


for i = 1:num_nodes
    Nodes{i}.ID = i;
    if isempty(Nodes{i}.Parent_Vessel)
        Nodes{i}.BC = 1;
    elseif isempty(Nodes{i}.Daughter_Vessel)
        Nodes{i}.BC = -1;
    else
        Nodes{i}.BC = 0;
    end
end


end

function [num] = Doubling(generations)
% Function that creates the bifurations
% Could be improved
num = 2*ones(1,2*generations);
for i = 0:generations
    num(i+1)= num(i+1)^i;
end
num(generations+1:end) = fliplr(num(1:generations));
end

function [new_radius] = BifurcateRadius(current_radius)
% Power law? for a the same radius, gamma = 2.7?
new_radius = current_radius*0.773;
end

function ID = getID(num_items, k)
% Creates the ID number
branches = 0:sum(num_items);
for i = 1:size(num_items,2)
    ID{i} = branches((1+k:k+num_items(i)));
    k = k + num_items(i);
end
end

function ID = getIDnode(num_items, k)
% Creates the ID number
branches = 0:sum(num_items);
for i = 1:size(num_items,2)
    ID{i} = branches((1+k:k+num_items(i))) + 1;
    k = k + num_items(i);
end
end

function length = getLength(radius)
if radius > 50
    length = (15.75*((radius/10000).^1.1))*10000/radius;
else
    length = (1.79*((radius/10000).^0.47))*10000/radius;
end

end