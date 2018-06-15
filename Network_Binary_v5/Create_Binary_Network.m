function [Vessel, Node] = Create_Binary_Network(num_bif,d_R_0)
%% Create_Binary_Network v4
%{
Function Creates the vessel segment infomation to be able to form a binary
tree network.
%}

%{
Inputs:
    num_bif % number of bifurations to the middle. ie expansion.
    d_R_0   % input vessel radius in micronmetres
Outputs:
    Vessel: Structure, vessel segment infomation
    Node: Structure, nodal infomation
%}

%{
Author = Michael Zhang
Date created = 09-06-18
%}

%% Parameters

% Generation is a "continuous range between two levels".
% Level is a "discrete point"

global num_vessels num_nodes
art2ven = 2;                        % Length ratio - venuoles hold 4 times more blood than artrioles./viens/artiers
num_gen = num_bif * 2;              % Number of generations

gen_vessels = Bifurcate(num_bif);    % Array, Number of vessels in each generation
lvl_nodes = [1, gen_vessels, 1];    % Array, Number of nodes in each level
lvl_nodes(num_bif + 1) = [];

num_vessels = sum(gen_vessels);         % Number of vessels total
gen_vessel_ID = getID(gen_vessels,1);   % ID for each vessel in each generation
lvl_node_ID = getIDnode(lvl_nodes,0);

% ID for each node in each level
lvl_node_ID{1} = 1;

gen_vessel_ID = [0, gen_vessel_ID, 0];
num_nodes = sum(lvl_nodes);
lvl_node_ID{end} = num_nodes;
% Creates the radii for each generation
gen_radii = zeros(1,num_gen);
gen_radii(1) = d_R_0; %% TODO more effecient.
for i = 1:num_bif-1
    [gen_radii(i+1)] = BifurcateRadius(gen_radii(i));
end
gen_radii(num_bif+1:end) = art2ven * fliplr(gen_radii(1:num_bif));

%% Create Nodes and Vessels
Vessel = cell(1,num_vessels); % Holds cell infomation
Node = cell(1,num_nodes);     % Holds nodal infomation

v = 1; % for each vessel
for i = 1:num_gen % for each generation
    for j = 1:gen_vessels(i)    % for each vessel in the generation
        Vessel{v}.ID = gen_vessel_ID{i+1}(j);
        Vessel{v}.Generation = i;
        
        Vessel{v}.Radius = gen_radii(i);
        Vessel{v}.Length = [];
        Vessel{v}.n_Radius = 1;
        Vessel{v}.n_Length = [];
        
        if i <= num_gen/2   % if on branching side of the binary tree
            Vessel{v}.Parent_Node = lvl_node_ID{i}(ceil(j/2));
            Vessel{v}.Daughter_Node = lvl_node_ID{i+1}(j);
            Vessel{v}.Length = getLength(Vessel{v}.Radius);
        else % if on the merging side of the binary tree
            Vessel{v}.Parent_Node = lvl_node_ID{i}(j);
            Vessel{v}.Daughter_Node = lvl_node_ID{i+1}(ceil(j/2));
            Vessel{v}.Length = 14.54*Vessel{v}.Radius; % Veins
        end
        
        Vessel{v}.n_Length = Vessel{v}.Length / Vessel{v}.Radius;
        v = v + 1;
    end
end

% Assign ndoal structure and temporary values
for n = 1:num_nodes 
    Node{n}.ID = n;
    Node{n}.Parent_Vessel=[];
    Node{n}.Daughter_Vessel=[];
    Node{n}.Parent_Node = [];
    Node{n}.Daughter_Node = [];
    Node{n}.num_Daughter_Nodes = 0;
    Node{n}.num_Parent_Nodes = 0;
end

% Assigns the appropriate parent and daugher vessels each node has
for v = 1:num_vessels 
    Node{Vessel{v}.Parent_Node}.Daughter_Vessel = [Node{Vessel{v}.Parent_Node}.Daughter_Vessel ,Vessel{v}.ID];
    Node{Vessel{v}.Daughter_Node}.Parent_Vessel = [Node{Vessel{v}.Daughter_Node}.Parent_Vessel, Vessel{v}.ID];
end

% Assigns the appropritte parent and daugher nodes that each
% node sees, used for connectivity
for n = 1:num_nodes
    for v = 1: num_vessels
        if Vessel{v}.Parent_Node == n
            Node{n}.Daughter_Node = [Node{n}.Daughter_Node ,Vessel{v}.Daughter_Node];
            Node{n}.num_Daughter_Nodes = Node{n}.num_Daughter_Nodes + 1;
        end
        if Vessel{v}.Daughter_Node== n
            Node{n}.Parent_Node  = [Node{n}.Parent_Node ,Vessel{v}.Parent_Node ];
            Node{n}.num_Parent_Nodes = Node{n}.num_Parent_Nodes + 1;
        end
    end
end

% Assigns whether a node is a entry and exit ie boundary
for n = 1:num_nodes
    if isempty(Node{n}.Parent_Vessel)
        Node{n}.BC = 1; % Input
    elseif isempty(Node{n}.Daughter_Vessel)
        Node{n}.BC = -1; % Output
    else
        Node{n}.BC = 0; % Otherwise
    end
end

end

function [bifurcations] = Bifurcate(num_bif)
% Creates the number of vessels for each generation level of the binary
% tree
% TODO: Could be improved
bifurcations = 2*ones(1,2*num_bif);
% 1,2,4,8,16 ...
for i = 0:num_bif
    bifurcations(i+1)= bifurcations(i+1)^i;
end
% Mirror
bifurcations(num_bif+1:end) = fliplr(bifurcations(1:num_bif));
end

function [r1] = BifurcateRadius(r0)
% Murrys power law for a the same radius,
% gamma = 2.7,  k = 0.774
% gamma = 3, k = 0.7937
% gamma = 2.333 , k == .7427
r1 = r0*0.774;
end

function ID = getID(num_items, k)
% Creates the ID number for each vessel in the network
branches = 0:sum(num_items);
for i = 1:size(num_items,2)
    ID{i} = branches((1+k:k+num_items(i)));
    k = k + num_items(i);
end
end

function ID = getIDnode(num_items, k)
% Creates the ID number for each node in the network
branches = 0:sum(num_items);
for i = 1:size(num_items,2)
    ID{i} = branches((1+k:k+num_items(i))) + 1;
    k = k + num_items(i);
end
end

function length = getLength(radius)
% Length to Radius ratio for pulmonary vessels arteries and veins
% Obtained from U. Qureshi 2013.
if radius > 50 %um arteries
    %formula is for radius in cm.
    length = (15.75*((radius/10000).^1.1))*10000/radius;
else
    length = (1.79*((radius/10000).^0.47))*10000/radius;
end
% Also for vein but currently in body above
end