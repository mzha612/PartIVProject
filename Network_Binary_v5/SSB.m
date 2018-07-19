function [Vessel,Node] = SSB(Vessel,Node,Vessel_Type)
%% Symmetric Segment Branching
%{
Topological growth model
Function that performs symetric segemnt growth or binary tree,
Adapted from RTB.m
Growth occurs stocatically on terminal nodes into two segments.
%}
%{
Inputs:
    Vessel: Emtpy Cell that will hold vessel segment infomation
    Node: Empty cell that will hold nodal infomation
    Vessel_Type: String, indicating which distribution to sample from
Outputs:
    Vessel: Updated Structure,containing vessel segment infomation
    Node: Updated Structure, containingnodal infomation
%}
%{
Author = Michael Zhang
Date created = 11-06-18
%}
%% Parameters
global num_vessels num_nodes num_generations num_capillaries
num_vessels = 1;
num_nodes = 2;

%% Initialise Inital nodes and vessel

% Currently hardcoded inital segments
% TODO: possible some arbitary starting points,
Node{1}.ID = 1;
Node{1} = InitaliseNode(Node{1});
Node{1}.xy = [-0.2,0.5];
Node{1}.num_Daughter_Nodes =  1;
Node{1}.Terminal = false;
Node{1}.BC = 1;
Node{1}.Daughter_Vessel = 1;
Node{1}.Daughter_Node = 2;

Node{2}.ID = 2;
Node{2} = InitaliseNode(Node{2});
Node{2}.xy = [0,0.5];
Node{2}.num_Daughter_Nodes =  0;
Node{2}.Parent_Vessel = 1;
Node{2}.Parent_Node = 1;
Node{2}.Generation = 1;

Vessel{1}.ID = 1;
Vessel{1} = InitaliseVessel(Vessel{1});
Vessel{1}.Parent_Node = 1;
Vessel{1}.Daughter_Node = 2;

%% Growth
c = 0;
c_max = ceil(log2(num_capillaries));
num_capillaries = 2^c_max;

while (c < c_max)
    seed = [];
    num_seed = 0; % Number of possible sprouting points
    for n = 1:num_nodes % Only sprout from terminal nodes
        % Could be done with a constant list and pushing/popping
        if Node{n}.Terminal
            num_seed = num_seed + 1;
            seed(num_seed) = n; % Save the node number,
        end
    end
    %     TODO: have an updating list rather than make it every single time
    % Randomly choose a seed from the list
    for j = 1:2^c
        
        sprout = seed(j);
        
        %% Create and update nodal/vessel infomation
        Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
        Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
        
        Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
        Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
        
        Node{sprout}.num_Daughter_Nodes =  2;
        Node{sprout}.Terminal = false;
        Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
        Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
        
        %% Re-evaluate
        num_nodes = numel(Node);
        num_vessels = numel(Vessel);
    end
    c = c + 1;
end

% Update Nodal infomation.
for n = 1:num_nodes
    if Node{n}.Terminal %termination criteria > num_capillaries
        Node{n}.BC = -1;
    end
    Node{n}.num_Parent_Nodes = size(Node{n}.Parent_Node,2);
end

% Assigns vessel daughter infomation
for n = 2:num_nodes
    for d = 1:Node{n}.num_Daughter_Nodes
        Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel = [Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel, Node{n}.Daughter_Vessel(d)];
    end
    Vessel{Node{n}.Parent_Vessel}.Daughter_Node = n;
end

%%

Diameter_Data = load('Frequency_Diameter.mat');
for v = 1:num_vessels
    % Finds the number of generations of each vessel
    Gen(v) =  Vessel{v}.Generation;
    
    % Samples radii from distribution
    Vessel{v}.Radius = Sample_Diameter(Vessel_Type,Diameter_Data)/2;
end
num_generations = max(Gen);

count = 0;
while count < num_vessels-1
    count = 0;
    for v = 2:num_vessels
        if Vessel{v}.Radius >= Vessel{Vessel{v}.Parent_Vessel}.Radius
            temp = Vessel{v}.Radius;
            Vessel{v}.Radius = Vessel{Vessel{v}.Parent_Vessel}.Radius;
            Vessel{Vessel{v}.Parent_Vessel}.Radius = temp;
        else
            count = count + 1;
        end
    end
end
Length_Data = load('Frequency_Length.mat');
for v = 1:num_vessels
    Vessel{v}.Length = Sample_Length(Vessel_Type, Length_Data);
    Vessel{v}.n_Length = Vessel{v}.Length/Vessel{v}.Radius;
end
% Number of vessels in each generation
for g = 1:num_generations
    v_gen(g) = sum(Gen == g);
end

% x space
x_coord = linspace(0,1,num_generations);

% y space
y_coord = cell(num_generations-1,1);
for g = 2:num_generations
    y_coord{g-1} = linspace(0,1,v_gen(g)+2);
end

% Assigns xy coordinates to each node depending on generation number
for n = 3:num_nodes
    Node{n}.xy = [x_coord(Node{n}.Generation), y_coord{Node{n}.Generation-1}(1+v_gen(Node{n}.Generation))];
    v_gen(Node{n}.Generation) = v_gen(Node{n}.Generation) - 1;
end

% Assigns xy coordiniates of the start and end of each vessel from the
% nodal coordinates, and calculates length of the vessel.
for v = 1:num_vessels
    Vessel{v}.xy_Start = Node{Vessel{v}.Parent_Node}.xy;
    Vessel{v}.xy_End = Node{Vessel{v}.Daughter_Node}.xy;
    Vessel{v}.xy_Length = norm([Vessel{v}.xy_End,Vessel{v}.xy_Start]);
end

if Vessel_Type == "Arteriole"
    for v = 1:num_vessels
        Vessel{v}.Arterial_Generation = Vessel{v}.Generation;
        Vessel{v}.Venule_Generation = [];
        Vessel{v}.Type = "Arteriole";
    end
elseif Vessel_Type == "Venule"
    for v = 1:num_vessels
        Vessel{v}.Arterial_Generation = [];
        Vessel{v}.Venous_Generation = Vessel{v}.Generation;
        Vessel{v}.Type = "Venule";
    end
else
    disp('error, vesseltype')
end

% for n = 1:num_nodes
%     disp(Node{n})
% end
%
% for v = 1:num_vessels
%     disp(Vessel{v})
% end


end

function Node = InitaliseNode(Node)
% Initalises the nodal infomation
Node.Generation = 0;
Node.xy = [];

Node.Daughter_Vessel = [];
Node.Daughter_Node = [];
Node.num_Daughter_Nodes = 0;

Node.Parent_Vessel = [];
Node.Parent_Node = [];
Node.num_Parent_Nodes = 0;

Node.BC = 0;
Node.Terminal = true;
end

function Vessel = InitaliseVessel(Vessel)
% Initalises the Vessel infomation
Vessel.Generation = 1;
Vessel.xy_Start = [];
Vessel.xy_End = [];

Vessel.Daughter_Vessel = [];
Vessel.Daughter_Node = [];

Vessel.Parent_Vessel = [];
Vessel.Parent_Node = [];
end

function Node = UpdateNode(Node,n,v,sprout)
% Updates nodal infomation when grown
Node{n}.ID = n;
Node{n} = InitaliseNode(Node{n});
Node{n}.Generation  = Node{sprout}.Generation + 1;
Node{n}.Parent_Node = sprout;
Node{n}.Parent_Vessel = v;
Node{n}.num_Parent_Nodes =  1;
end

function Vessel = UpdateVessel(Vessel,Node,v,sprout)
% Updates nodal infomation when grown
Vessel{v}.ID = v;
Vessel{v} = InitaliseVessel(Vessel{v});
Vessel{v}.Parent_Node = sprout;
Vessel{v}.Parent_Vessel = Node{sprout}.Parent_Vessel;
Vessel{v}.Generation = Vessel{Vessel{v}.Parent_Vessel}.Generation + 1;
end