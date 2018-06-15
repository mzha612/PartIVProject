function [Vessel,Node] = CombineTree(tree1,tree2)
%% CombineTree v2
%{
Grows a RSB or RTB for left and right side
Inverts infomation of, one side
Adds cappillaries to the ends to "Connect" - from a distribution

Inputs:
    tree1,tree2L: Function handles that referes to the differnt topological
    models
Outputs:
    Vessel: Updated Structure,containing vessel segment infomation
    Node: Updated Structure, containing nodal infomation
%}

%{
Author = Michael Zhang
Date created = 15-06-18
%}
%% Grow left and right tree
Algorithm1 = tree1;
Algorithm2 = tree2;

global num_vessels num_nodes num_generations num_capillaries
%
% num_capillaries = 180;

Vessel = cell(1,1);
Node = cell(1,1);

% Create left tree
[Vessel1,Node1] = Algorithm1(Vessel,Node,"Arteriole");
Total_Nodes = numel(Node1);
Total_Vessels = numel(Vessel1);
num_gen1 = num_generations;

% Create right tree
[Vessel2,Node2] = Algorithm2(Vessel,Node,"Venuole");
Total_Nodes = Total_Nodes + numel(Node2);
Total_Vessels =Total_Vessels + numel(Vessel2);
num_gen2 = num_generations;
%% Flip right side

% TODO arterial generation and veunole generation number
for n = 1:numel(Node2)
    Node2{n}.Generation = num_gen2 - Node2{n}.Generation ;
    Node2{n}.xy = [-1,1].*Node2{n}.xy + [3 0];
    [Node2{n}.Daughter_Node,Node2{n}.Parent_Node] = deal(  Node2{n}.Parent_Node,Node2{n}.Daughter_Node);
    [Node2{n}.Daughter_Vessel,Node2{n}.Parent_Vessel] = deal(  Node2{n}.Parent_Vessel,Node2{n}.Daughter_Vessel);
    [Node2{n}.num_Daughter_Nodes,Node2{n}.num_Parent_Nodes] = deal(  Node2{n}.num_Parent_Nodes,Node2{n}.num_Daughter_Nodes);
    Node2{n}.BC = Node2{n}.BC*-1;
end

for v = 1:numel(Vessel2)
    Vessel2{v}.Generation = num_gen2 - Vessel2{v}.Generation ;
    Vessel2{v}.xy_Start = [-1,1].*Vessel2{v}.xy_Start + [3 0];
    Vessel2{v}.xy_End = [-1,1].*Vessel2{v}.xy_End + [3 0];
    [Vessel2{v}.Daughter_Node,Vessel2{v}.Parent_Node] = deal(  Vessel2{v}.Parent_Node,Vessel2{v}.Daughter_Node);
    [Vessel2{v}.Daughter_Vessel,Vessel2{v}.Parent_Vessel] = deal(  Vessel2{v}.Parent_Vessel,Vessel2{v}.Daughter_Vessel);
end

%% Combine Both sides with vessels.

t1 = 0;
t2 = 0;
for n = 1:numel(Node1)
    if Node1{n}.Terminal
        t1 = t1 + 1;
        terminal_node1(t1) = n;
    end
end
for n = 1:numel(Node2)
    if Node2{n}.Terminal
        t2 = t2 + 1;
        terminal_node2(t2) = n;
    end
end

% Randomly swap terminal end list, such that random connections can be made
terminal_node1 = terminal_node1(randperm(num_capillaries));
terminal_node2 = terminal_node2(randperm(num_capillaries));

%% Middle connections

Vessel3 = cell(1,1);
for c = 1:num_capillaries
    n1 = terminal_node1(c);
    n2 = terminal_node2(c);
    
    Vessel3{c}.ID = c + Total_Vessels;
    Vessel3{c}.Generation = num_gen1 + 1; %% Need to fix
    Vessel3{c}.xy_Start = Node1{n1}.xy;
    Vessel3{c}.xy_End = Node2{n2}.xy;
    Vessel3{c}.Daughter_Vessel = Node2{n2}.Daughter_Vessel + numel(Vessel1);
    Vessel3{c}.Daughter_Node = n2 + numel(Node1);
    Vessel3{c}.Parent_Vessel= Node1{n1}.Parent_Vessel;
    Vessel3{c}.Parent_Node = n1;
    Vessel3{c}.xy_Length = 1;
    
    Vessel1{Node1{n1}.Parent_Vessel}.Daughter_Vessel = c+ Total_Vessels;
    Vessel2{Node2{n2}.Daughter_Vessel}.Parent_Vessel = c+ Total_Vessels - numel(Vessel1);
    
    Node1{n1}.Daughter_Vessel = c + Total_Vessels;
    Node1{n1}.Daughter_Node = n2 + numel(Node1);
    Node2{n2}.Parent_Vessel = c + Total_Vessels - numel(Vessel1);
    Node2{n2}.Parent_Node = n1 - numel(Node1);
    
    Node1{n1}.num_Daughter_Nodes = 1;
    Node1{n1}.num_Parent_Nodes = 1;
    
    Node2{n2}.num_Daughter_Nodes = 1;
    Node2{n2}.num_Parent_Nodes = 1;
    
    Node1{n1}.BC = 0;
    Node1{n1}.Terminal = 0;
    
    Node2{n2}.BC = 0;
    Node2{n2}.Terminal = 0;
end

% Assign random vessel radius
for c = 1:num_capillaries
    Vessel3{c}.Radius = inf;
    % Make sure the random radius is smaller than both parent and daughter.
    while ~(Vessel3{c}.Radius < Vessel1{Vessel3{c}.Parent_Vessel}.Radius )...
            && ~(Vessel3{c}.Radius < Vessel2{Vessel3{c}.Daughter_Vessel - numel(Vessel1)}.Radius)
        Vessel3{c}.Radius = Sample_Diameter("Capillary")/2;
    end
    % aSsign random length from sample distribution
    Vessel3{c}.Length = Sample_Length("Capillary");
    Vessel3{c}.n_Length = Vessel3{c}.Length/Vessel3{c}.Radius;
end

%% ReAssign infomation.
% Right hand side
for n = 1:numel(Node2)
    Node2{n}.ID =  Node2{n}.ID + numel(Node1) ;
    Node2{n}.Daughter_Node = Node2{n}.Daughter_Node + numel(Node1);
    Node2{n}.Daughter_Vessel = Node2{n}.Daughter_Vessel + numel(Vessel1);
    Node2{n}.Parent_Vessel = Node2{n}.Parent_Vessel + numel(Vessel1);
    Node2{n}.Parent_Node = Node2{n}.Parent_Node + numel(Node1);
end

for v = 1:numel(Vessel2)
    Vessel2{v}.ID = Vessel2{v}.ID + numel(Vessel1);
    Vessel2{v}.Daughter_Node = Vessel2{v}.Daughter_Node + numel(Node1) ;
    Vessel2{v}.Daughter_Vessel = Vessel2{v}.Daughter_Vessel + numel(Vessel1) ;
    if Node2{Vessel2{v}.Daughter_Node - numel(Node1)}.BC == -1
        % if its the terminal vessel
        Vessel2{v}.Daughter_Vessel = [];
    end
    Vessel2{v}.Parent_Node = Vessel2{v}.Parent_Node + numel(Node1) ;
    Vessel2{v}.Parent_Vessel = Vessel2{v}.Parent_Vessel + numel(Vessel1);
end

%% Combine
Node = [Node1 Node2];
Vessel = [Vessel1, Vessel2, Vessel3];

num_nodes = numel(Node);
num_vessels = numel(Vessel);
end


