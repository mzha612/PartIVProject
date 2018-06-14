function [Vessel,Node] = RSB(Vessel,Node)
%% Random Segment Branching
%{
Topological growth model
Function that performs random segment growth.
segment Growth occurs stocatically at any node.
%}

%{
Author = Michael Zhang
Date created = 11-06-18
%}
%% Parameters
global num_vessels num_nodes num_generations num_capillaries

%% Initialise Inital nodes and vessel
num_vessels = 1;
num_nodes = 2;

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
Vessel{1}.Parent_Vessel = 0;

%% Growth
c = 0; % Count of capillaries/or terminal segments, currently
while (c < num_capillaries)
    c = 0;
    sprout = (randi([2,num_nodes],1)); % any segment (NODE) can bifurcate
    
    %% Different segment branching scenarios
    if Node{sprout}.Terminal % If bifurcating from a terminal node, do the same as the RTB
        Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
        Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
        
        Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
        Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
        
        Node{sprout}.num_Daughter_Nodes = 2;
        Node{sprout}.Terminal = false;
        
        Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
        Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
    else % If bifurcating from a non terminal node
        
        for d = 1:Node{sprout}.num_Daughter_Nodes
            % Store the previous daughter nodes and vessels
            old_d_Nodes(d) = Node{sprout}.Daughter_Node(d);
            old_d_Vessels(d) = Node{sprout}.Daughter_Vessel(d);
            
            % Update the daughter node/vessel with the new parent
            % infomation
            Node{Node{sprout}.Daughter_Node(d)}.Parent_Node = num_nodes + 2;
            % +2 refers to all the duaghters splitting to the second bifurcation
            % This doesn't matter unless it is affected by some thing that
            % depends on the way it splits.
            Vessel{Node{sprout}.Daughter_Vessel(d)}.Parent_Node = num_nodes + 2;
            Vessel{Node{sprout}.Daughter_Vessel(d)}.Parent_Vessel = num_vessels + 2;
            
            % Recursively update each of the sprouting node's current
            % daughters
            [Vessel,Node] = UpdateDaughters(Vessel,Node,Node{sprout}.Daughter_Node(d));
        end
        
        % Create New Daughter Node
        Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
        Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
        
        % Create New terminal Node
        Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
        Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
        
        % update the sprouting node with new daughter infomation
        Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
        Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
        
        % Update thew new nodes with infomation about the previous
        % daughters of the sprouting node ie (thier new daughters)
        Node{num_nodes + 1}.Terminal = true;
        Node{num_nodes + 2}.Terminal = false;
        % +2 here, hardcoded, branch adding onto from the second bifurcation.
        Node{num_nodes + 1}.Daughter_Node = [];
        Node{num_nodes + 2}.Daughter_Node = old_d_Nodes;
        
        Node{num_nodes + 1}.Daughter_Vessel = [];
        Node{num_nodes + 2}.Daughter_Vessel = old_d_Vessels;
        
        Node{num_nodes + 1}.num_Daughter_Nodes = 0;
        Node{num_nodes + 2}.num_Daughter_Nodes = 2;
    end
    
    %% Re-evaluate
    num_nodes = numel(Node);
    num_vessels = numel(Vessel);
    for n = 1:num_nodes
        if Node{n}.Terminal %terminatioion criteria > num_capillaries
            c = c + 1;
        end
    end
end

% Update Nodal infomation.
for n = 1:num_nodes
    if Node{n}.Terminal %termination criteria > num_capillaries
        Node{n}.BC = -1;
    end
    Node{n}.num_Parent_Nodes = size(Node{n}.Parent_Node,2);
end

% Assigns Vessel daughter infomation
for n = 2:num_nodes
    for d = 1:Node{n}.num_Daughter_Nodes
        Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel = ...
            [Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel, Node{n}.Daughter_Vessel(d)];
    end
    Vessel{Node{n}.Parent_Vessel}.Daughter_Node = n;
end

%% Coordinates for visualisation

% ----------------------------------------------------------------------- %
%% TODO make it sample from distribution


% Finds the number of generations of each vessel
for v = 1:num_vessels
    Gen(v) =  Vessel{v}.Generation;
end
num_generations = max(Gen);

radii = linspace(5,30,num_generations+1);
length = 300;
for v = 1:num_vessels
    Vessel{v}.Radius = radii(Vessel{v}.Generation);
    Vessel{v}.Length = 300;
    Vessel{v}.n_Length = 300/Vessel{v}.Radius;
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
    Node{n}.xy = [x_coord(Node{n}.Generation), y_coord{Node{n}.Generation-1}(end - v_gen(Node{n}.Generation))];
    v_gen(Node{n}.Generation) = v_gen(Node{n}.Generation) - 1;
end

% Assigns xy coordiniates of the start and end of each vessel from the
% nodal coordinates, and calculates length of the vessel.
for v = 1:num_vessels
    Vessel{v}.xy_Start = Node{Vessel{v}.Parent_Node}.xy;
    Vessel{v}.xy_End = Node{Vessel{v}.Daughter_Node}.xy;
    Vessel{v}.Length = norm([Vessel{v}.xy_End,Vessel{v}.xy_Start]);
end
% ----------------------------------------------------------------------- %
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

function [Vessel,Node] = UpdateDaughters(Vessel,Node,n)
% Recursively update each daughter with updated infomation when 'shifted
% down'
Vessel{Node{n}.Parent_Vessel}.Generation =  Vessel{Node{n}.Parent_Vessel}.Generation + 1;
Node{n}.Generation  = Node{n}.Generation + 1;
if ~Node{n}.Terminal
    for d = 1:Node{n}.num_Daughter_Nodes
        [Vessel,Node] = UpdateDaughters(Vessel,Node,Node{n}.Daughter_Node(d));
    end
end

end


