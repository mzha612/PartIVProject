function [Vessel,Node] = Create_Real_Rat_Mesentery()
%% Create_Binary_Network v4
%{
Function Creates the vessel segment infomation to be able to solve the rat mesentery
data obtained from secomb. 1990 - 913 vessels.
%}

%{
Inputs:
    Implicit(Datafile) 'Rat_Mesentery_1990.Mat'
Outputs:
    Vessel: Structure, vessel segment infomation
    Node: Structure, nodal infomation
%}

%{
Author = Michael Zhang
Date created = 09-06-18

Most of the code is recylced from create binary network v4
%}
%% Parameters

global num_vessels num_nodes
num_vessels = 913;
num_boundary = 65;

%% Open and read in file
% Data first imported into matlab, then saved.
% Boundary, Diameter, ID, Length, NodeIn, NodeOut
load('Rat_Mesentery_1990.Mat');

%% Create and assign Vesssel and Nodal infomation

% Vessel names are not consistant with vessel ID numbers, ie certain names
% are missing. and not in order
vessel_name_pool = zeros(1000,1);
for v = 1:num_vessels
    vessel_name_pool(ID(v)) = ID(v);
end
vessel_name_pool(vessel_name_pool == 0) = [];

node_name_pool = zeros(1000,1);
for v = 1:num_vessels
    node_name_pool(NodeIn(v)) = NodeIn(v);
    node_name_pool(NodeOut(v)) = NodeOut(v);
end
node_name_pool(node_name_pool == 0) = [];
num_nodes = size(node_name_pool,1);

% Predefining size
Vessel = cell(num_vessels,1);
Node = cell(num_nodes,1);

% Initialise vessel infomation, from the data set
for v = 1:num_vessels
    Vessel{v}.Vessel_ID = v;
    Vessel{v}.Name = vessel_name_pool(v);
    Vessel{v}.Parent_Node = find(node_name_pool==NodeIn(v));
    Vessel{v}.Daughter_Node = find(node_name_pool==NodeOut(v));
    Vessel{v}.Radius = Diameter(v) / 2;
    Vessel{v}.Length = Length(v);
    Vessel{v}.n_Radius = 1;
    Vessel{v}.n_Length = Vessel{v}.Length / Vessel{v}.Radius;
end

% Initalise nodal information
for n = 1:num_nodes
    Node{n}.Node_ID = n;
    Node{n}.Name = node_name_pool(n);
    Node{n}.BC = 0;
    Node{n}.Parent_Vessel=[];
    Node{n}.Daughter_Vessel=[];
    Node{n}.Parent_Node = [];
    Node{n}.Daughter_Node = [];
    Node{n}.num_Daughter_Nodes = 0;
    Node{n}.num_Parent_Nodes = 0;
end

% Assigns/appends the appropriate parent and daugher vessels each node has
for v = 1:num_vessels
    Node{Vessel{v}.Parent_Node}.Daughter_Vessel = [Node{Vessel{v}.Parent_Node}.Daughter_Vessel ,Vessel{v}.Vessel_ID];
    Node{Vessel{v}.Daughter_Node}.Parent_Vessel = [Node{Vessel{v}.Daughter_Node}.Parent_Vessel, Vessel{v}.Vessel_ID];
end

% Assigns the appropriate paraent and daughter vessesls each vessel has
for v = 1: num_vessels
    Vessel{v}.Parent_Vessel = Node{Vessel{v}.Parent_Node}.Parent_Vessel;
    Vessel{v}.Daughter_Vessel = Node{Vessel{v}.Daughter_Node}.Daughter_Vessel;
end

% Assigns the appropriate parent and daugher nodes that each
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

% For each boundary node, find whether it is a input or an output node, and
% assign BC value
for b = 1:num_boundary
    for v = 1:num_vessels
        if ~isempty(find(Boundary(b) == NodeIn(v),1))
            Node{find(node_name_pool==NodeIn(v))}.BC = 1; % input
        end
        if ~isempty(find(Boundary(b) == NodeOut(v),1))          
            Node{find(node_name_pool==NodeOut(v))}.BC = -1; % Output
        end
    end
end

end