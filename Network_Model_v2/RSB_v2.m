function [Vessel,Node] = RSB(Vessel,Node,Vessel_Type)
%% Random Segment Branching v2
%{
Topological growth model
Function that performs random segment growth.
segment Growth occurs stocatically at any node.
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

%% Initialise Inital nodes and vessel
num_vessels = 1;
num_nodes = 2;

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
Node{1}.Generation = 0;

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
Vessel{1}.Generation = 1;

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

[Vessel,Node] = UpdateDaughters(Vessel,Node,Node{1}.Daughter_Node);


if (Vessel_Type == "Arteriole")
%     % No weights
%     P_d =[4.687488329463876;2.465387856758793]% gamma
%     P_l =   [1.002871132180170;2.888884270441188e+02]% gamma
%     
% 
%     P_d = [4.558909962071017;2.541429434656689];% gamma
%     P_l = [1.309264848264005;2.051390713688520e+02]; % gamma
%     
%     dist_l = 'Gamma';
%     
% %     P_d = [14.081118405943762;2.626613312178127]; % wei
% % %     P_l = [1.257530962313703;2.342497763481835e+02]; % wei
%     dist_d = 'Gamma';

%     Diameter_Data = readtable('art_dia_gamma.csv');
%     Diameter_Data = table2array(Diameter_Data);
%     Length_data = readtable('art_len_linear.csv');
%     Length_data = table2array(Length_data);
elseif (Vessel_Type == "Venule")
%     No weights
 P_d = [2.791278438217233;6.368312347248351];% gamma
   P_l = [1.198108031965868;2.425088293496245e+02];% gamma
   P_l=  [1.715404811576471;1.948040327217233e+02];
   
   
   P_d = [3.820106390430977;4.307168289097728];
% %    last half 5
%     P_d = [3.202898614399205;5.790492383117489]; % gamma
%     P_l = [1.555822060711918;1.922850913173208e+02]; % gamma
    dist_l = 'Gamma';
    dist_d = 'Gamma';
    
% %     
%     % First half 5
%     P_d = [2.939412230957511;5.951741024615105]; % gamma
%     P_l = [1.071743460163905;2.830265577626656e+02]; % gamma

% %     No weights
%         P_d = [19.685876560537350;1.766938522123263]; % weibull  
%     P_l = [3.011502415291832e+02;1.127054787436254]; % weibull
% % %     
% % % %     P_d = [21.051966968392160;1.895731923533652]; % weibull  
% % % % %     P_l = [3.064548114327770e+02;1.148573351456710]; % weibull
% %     
%     dist_d = 'Weibull';
%     dist_l = 'Weibull';
%     Diameter_Data = readtable('ven_dia_gamma.csv');
%     Diameter_Data = table2array(Diameter_Data);
%     Length_data = readtable('ven_len_gamma .csv');
%     Length_data = table2array(Length_data);
else
    disp('error sample diameter') 
end
%     Diameter_Data = load('Frequency_Diameter.mat');

% p = [39.468958081272575,-0.459462643821905,8.657837845701500];
% modelFun =  @(x) p(1)*x.^(p(2)) + p(3);
    
for v = 1:num_vessels
    % Finds the number of generations of each vessel
    Gen(v) =  Vessel{v}.Generation;
    
    % Samples radii from distribution
%     Vessel{v}.Radius = Sample_Diameter(Vessel_Type,Diameter_Data)/2;
%         Vessel{v}.Radius = Diameter_Data(randi(size(Diameter_Data,1),[1,1]))/2*10^-6;
%     Vessel{v}.Radius = min(icdf(dist_d,rand(1),P_d(1),P_d(2)) + 2.2,50);
%     Vessel{v}.Radius = Vessel{v}.Radius/2*10^-6;
%     Vessel{v}.Radius = 25/2*10^-6;
% 
    lll = icdf(dist_d,rand(1),P_d(1),P_d(2));
%     while lll > 50
%         lll = icdf(dist_d,rand(1),P_d(1),P_d(2));
%     end
    Vessel{v}.Radius = lll + 2.5;
    Vessel{v}.Radius =  Vessel{v}.Radius/2*10^-6;
%     Vessel{v}.Radius = modelFun(Vessel{v}.Generation)/2*10^-6;
end
num_generations = max(Gen);

count = 0;

while count < num_vessels-1
    count = 0;
    % TODO: inital segment radii
    for v = 2:num_vessels
        % Swaps radii depending if bigger or not.
        if Vessel{v}.Radius  >= Vessel{Vessel{v}.Parent_Vessel}.Radius+ 1e-6
            % Could use "deal"
            temp = Vessel{v}.Radius;
            Vessel{v}.Radius = Vessel{Vessel{v}.Parent_Vessel}.Radius;
            Vessel{Vessel{v}.Parent_Vessel}.Radius = temp;
        else
            count = count + 1;
        end
    end
end

% Sample vessel lengths
% Length_Data = load('Frequency_Length.mat');

for v = 1:num_vessels
%         Vessel{v}.Length = Length_data(randi(size(Length_data,1),[1,1]))*10^-6;
    

    Vessel{v}.Length = (icdf(dist_l,rand(1),P_l(1),P_l(2)) + 50)*10^-6;
%     Vessel{v}.Length = (500)*10^-6;
%     Vessel{v}.Length = Sample_Length(Vessel_Type,Length_Data);
    Vessel{v}.n_Length = Vessel{v}.Length/Vessel{v}.Radius;
end

%% XY Coordinates

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
    if n == 125
       disp('hi') 
    end
    Node{n}.xy = [x_coord(Node{n}.Generation), y_coord{Node{n}.Generation-1}(end - v_gen(Node{n}.Generation))];
    v_gen(Node{n}.Generation) = v_gen(Node{n}.Generation) - 1;
end

% Assigns xy coordiniates of the start and end of each vessel from the
% nodal coordinates, and calculates length of the vessel.
for v = 1:num_vessels
    Vessel{v}.xy_Start = Node{Vessel{v}.Parent_Node}.xy;
    Vessel{v}.xy_End = Node{Vessel{v}.Daughter_Node}.xy;
    %     Vessel{v}.xy_Length = norm([Vessel{v}.xy_End,Vessel{v}.xy_Start]);
end
%%
% for n = 1:num_nodes
%     disp(Node{n})
% end
%
% for v = 1:num_vessels
%     disp(Vessel{v})
% end
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
end

function Node = InitaliseNode(Node)
% Initalises the nodal infomation
Node.Generation = 0;
Node.xy = [];
Node.Growth_Angle = [];
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
% Vessel{Node{n}.Parent_Vessel}.Generation =  Vessel{Node{n}.Parent_Vessel}.Generation + 1;
Node{n}.Generation  = Node{n}.Generation + 1;
Node{n}.Generation  = Node{Node{n}.Parent_Node}.Generation + 1;
if ~Node{n}.Terminal
    for d = 1:Node{n}.num_Daughter_Nodes
        Vessel{Node{n}.Daughter_Vessel(d)}.Generation =  Vessel{Node{n}.Parent_Vessel}.Generation + 1;
        [Vessel,Node] = UpdateDaughters(Vessel,Node,Node{n}.Daughter_Node(d));
    end
end

end
