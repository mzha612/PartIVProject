function [Vessel,Node] = AR_Terminal(Vessel,Node,Vessel_Type)
%% AR_Terminal - Network_Model_v2
%{
    Area restricted growth algorthm for which terminal branching occurs
%}
%{
Input:

Output:

%}
%{
    Michael Zhang
    21-7-18
%}
% close
%% Parameters
% rng(2)
global num_vessels num_nodes num_capillaries num_generations
% Vessel_Type
% num_capillaries = 18;
num_vessels = 1;
num_nodes = 2;
%
% Vessel_Type = "Arteriole";

global x_min x_max y_min y_max


% Initalise the first nodes + vessel

Node{1}.ID = 1;
Node{1} = InitaliseNode(Node{1});
Node{1}.xy = [0,0];
Node{1}.num_Daughter_Nodes =  1;
Node{1}.Terminal = false;
Node{1}.BC = 1;
Node{1}.Daughter_Vessel = 1;
Node{1}.Daughter_Node = 2;
Node{1}.num_Parent_Nodes = 0;

Node{2}.ID = 2;
Node{2} = InitaliseNode(Node{2});
Node{2}.xy = [0,Sample_Length(Vessel_Type)];
Node{2}.num_Daughter_Nodes =  0;
Node{2}.Parent_Vessel = 1;
Node{2}.Parent_Node = 1;
Node{2}.Generation = 1;
Node{2}.num_Parent_Nodes = 1;

Vessel{1}.ID = 1;
Vessel{1} = InitaliseVessel(Vessel{1});
Vessel{1}.Parent_Node = 1;
Vessel{1}.Daughter_Node = 2;
Vessel{1}.xy_Start = Node{Vessel{1}.Parent_Node}.xy;
Vessel{1}.xy_End = Node{Vessel{1}.Daughter_Node}.xy;

new_nodes = [1,2];
new_vessels = 1;

figure
hold on
plot(x_min,y_min,x_max,y_max)
for nn = 1:2
    plot(Node{new_nodes(nn)}.xy,'bo')
end
plot([Vessel{new_vessels(1)}.xy_Start;Vessel{new_vessels(1)}.xy_End],'r')

c = 0;
while (c < num_capillaries)
    
    
    c = 0;
    %     isValid = 0;
    %     while(~isValid)
    %         seed = [];
    %         num_seed = 0; % Number of possible sprouting points
    %         for n = 1:num_nodes % Only sprout from terminal nodes
    %             % Could be done with a constant list and pushing/popping
    %             if Node{n}.Terminal
    %                 num_seed = num_seed + 1;
    %                 seed(num_seed) = n; % Save the node number,
    %             end
    %         end
    %         %     TODO: have an updating list rather than make it every single time
    %         % Randomly choose a seed from the list
    %         sprout = seed(randi([1,num_seed],1));
    %
    %
    %
    %         %% Create and update nodal/vessel infomation
    %         Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
    %         Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
    %
    %         Node{sprout}.num_Daughter_Nodes =  2;
    %         Node{sprout}.Terminal = false;
    %         Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
    %         Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
    %
    %         Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
    %         Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
    %
    %         Vessel{num_vessels + 1}.xy_Start =  Node{Vessel{num_vessels+1}.Parent_Node}.xy;
    %         Vessel{num_vessels + 1}.xy_End = Node{Vessel{num_vessels+1}.Daughter_Node(1)}.xy;
    %         Vessel{num_vessels + 2}.xy_Start =  Node{Vessel{num_vessels+2}.Parent_Node}.xy;
    %         Vessel{num_vessels + 2}.xy_End = Node{Vessel{num_vessels+2}.Daughter_Node(2)}.xy;
    %
    %         isValid = Check_Valid(Vessel,Node,sprout);
    %         if ~isValid
    %             Node{sprout}.Terminal = true;
    %         end
    %         %         for n = 1:num_nodes+2
    %         %             disp(Node{n})
    %         %         end
    
    
    isValid = 0;
    while(~isValid)
        N_Temp = Node;
        V_Temp = Vessel;
        seed = [];
        num_seed = 0; % Number of possible sprouting points
        for n = 1:num_nodes % Only sprout from terminal nodes
            % Could be done with a constant list and pushing/popping
            if N_Temp{n}.Terminal
                num_seed = num_seed + 1;
                seed(num_seed) = n; % Save the node number,
            end
        end
        %     TODO: have an updating list rather than make it every single time
        % Randomly choose a seed from the list
        sprout = seed(randi([1,num_seed],1));
        
        %% Create and update nodal/vessel infomation
        N_Temp = UpdateNode(N_Temp,num_nodes + 1,num_vessels + 1,sprout);
        N_Temp = UpdateNode(N_Temp,num_nodes + 2,num_vessels + 2,sprout);
        
        N_Temp{sprout}.num_Daughter_Nodes =  2;
        N_Temp{sprout}.Terminal = false;
        N_Temp{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
        N_Temp{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
        
        V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + 1,sprout);
        V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + 2,sprout);
        
        V_Temp{num_vessels + 1}.xy_Start =  N_Temp{V_Temp{num_vessels+1}.Parent_Node}.xy;
        V_Temp{num_vessels + 1}.xy_End = N_Temp{V_Temp{num_vessels+1}.Daughter_Node(1)}.xy;
        V_Temp{num_vessels + 2}.xy_Start =  N_Temp{V_Temp{num_vessels+2}.Parent_Node}.xy;
        V_Temp{num_vessels + 2}.xy_End = N_Temp{V_Temp{num_vessels+2}.Daughter_Node(2)}.xy;
        
        isValid = Check_Valid(V_Temp,N_Temp,sprout);
        if ~isValid
            N_Temp{sprout}.Terminal = true;
        end
        %         for n = 1:num_nodes+2
        %             disp(Node{n})
        %         end
    end
    Node = N_Temp ;
    Vessel= V_Temp;
    
    
    %% Re-evaluate
    
    
    
    new_nodes = [num_nodes + 1,num_nodes + 2];
    new_vessels = [num_vessels + 1,num_vessels + 2];
    
    for nn = 1:size(new_nodes,2)
        plot(Node{new_nodes(nn)}.xy(1),Node{new_nodes(nn)}.xy(2),'bo')
        %         pause(0.01)
    end
    for nv = 1:size(new_vessels,2)
        plot([Vessel{new_vessels(nv)}.xy_Start(1),Vessel{new_vessels(nv)}.xy_End(1)],[Vessel{new_vessels(nv)}.xy_Start(2),Vessel{new_vessels(nv)}.xy_End(2)],'r')
        %         pause(0.01)
    end
    
    num_nodes = numel(Node);
    num_vessels = numel(Vessel);
    for n = 1:num_nodes
        if Node{n}.Terminal %termination criteria > num_capillaries
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

% Assigns vessel daughter infomation
for n = 2:num_nodes
    for d = 1:Node{n}.num_Daughter_Nodes
        Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel = [Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel, Node{n}.Daughter_Vessel(d)];
    end
    Vessel{Node{n}.Parent_Vessel}.Daughter_Node = n;
end

% No weights
P_d =[4.687488329463876;2.465387856758793];% gamma
P_l =   [1.002871132180170;2.888884270441188e+02];% gamma
dist_l = 'Gamma';
dist_d = 'Gamma';
for v = 1:num_vessels
    % Finds the number of generations of each vessel
    Gen(v) =  Vessel{v}.Generation;
    
    % Samples radii from distribution
    Vessel{v}.Radius = min(icdf(dist_d,rand(1),P_d(1),P_d(2)) + 2.2,50);
    Vessel{v}.Radius = Vessel{v}.Radius/2*10^-6;
    Vessel{v}.Radius = Sample_Diameter(Vessel_Type)/2*10^-6;
end
num_generations = max(Gen);

count = 0;
while count < num_vessels-1
    count = 0;
    for v = 2:num_vessels
        if Vessel{v}.Radius >= Vessel{Vessel{v}.Parent_Vessel}.Radius + 1e-6
            temp = Vessel{v}.Radius;
            Vessel{v}.Radius = Vessel{Vessel{v}.Parent_Vessel}.Radius;
            Vessel{Vessel{v}.Parent_Vessel}.Radius = temp;
        else
            count = count + 1;
        end
    end
end
% Length_Data = load('Frequency_Length.mat');
for v = 1:num_vessels
    Vessel{v}.Length = norm(Vessel{v}.xy_End - Vessel{v}.xy_Start) * 10^-6;
    %     Vessel{v}.Length = Sample_Length(Vessel_Type, Length_Data);
    Vessel{v}.n_Length = Vessel{v}.Length/Vessel{v}.Radius;
end
% Number of vessels in each generation
for g = 1:num_generations
    v_gen(g) = sum(Gen == g);
end

% x space
% x_coord = linspace(0,1,num_generations);

% y space
% y_coord = cell(num_generations-1,1);
% for g = 2:num_generations
%     y_coord{g-1} = linspace(0,1,v_gen(g)+2);
% end

% % Assigns xy coordinates to each node depending on generation number
% for n = 3:num_nodes
%     Node{n}.xy = [x_coord(Node{n}.Generation), y_coord{Node{n}.Generation-1}(1+v_gen(Node{n}.Generation))];
%     v_gen(Node{n}.Generation) = v_gen(Node{n}.Generation) - 1;
% end

% Assigns xy coordiniates of the start and end of each vessel from the
% nodal coordinates, and calculates length of the vessel.
% for v = 1:num_vessels
%     Vessel{v}.xy_Start = Node{Vessel{v}.Parent_Node}.xy;
%     Vessel{v}.xy_End = Node{Vessel{v}.Daughter_Node}.xy;
%     Vessel{v}.xy_Length = norm([Vessel{v}.xy_End,Vessel{v}.xy_Start]);
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


%
% for n = 1:num_nodes
%     disp(Node{n})
% end
%
%
% for v = 1:num_vessels
%     disp(Vessel{v})
% end


end




%% function [Vessel,Node] = RTB(Vessel,Node,Vessel_Type)
%
%
% %% Growth
% c = 0;  % Count of capillaries/or terminal segments, currently
% while (c < num_capillaries)
%     c = 0;
%     seed = [];
%     num_seed = 0; % Number of possible sprouting points
%     for n = 1:num_nodes % Only sprout from terminal nodes
%         % Could be done with a constant list and pushing/popping
%         if Node{n}.Terminal
%             num_seed = num_seed + 1;
%             seed(num_seed) = n; % Save the node number,
%         end
%     end
%     %     TODO: have an updating list rather than make it every single time
%     % Randomly choose a seed from the list
%     sprout = seed(randi([1,num_seed],1));
%
%     %% Create and update nodal/vessel infomation
%     Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
%     Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
%
%     Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
%     Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
%
%     Node{sprout}.num_Daughter_Nodes =  2;
%     Node{sprout}.Terminal = false;
%     Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
%     Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
%
%     %% Re-evaluate
%     num_nodes = numel(Node);
%     num_vessels = numel(Vessel);
%     for n = 1:num_nodes
%         if Node{n}.Terminal %termination criteria > num_capillaries
%             c = c + 1;
%         end
%     end
% end
%
% % Update Nodal infomation.
% for n = 1:num_nodes
%     if Node{n}.Terminal %termination criteria > num_capillaries
%         Node{n}.BC = -1;
%     end
%     Node{n}.num_Parent_Nodes = size(Node{n}.Parent_Node,2);
% end
%
% % Assigns vessel daughter infomation
% for n = 2:num_nodes
%     for d = 1:Node{n}.num_Daughter_Nodes
%         Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel = [Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel, Node{n}.Daughter_Vessel(d)];
%     end
%     Vessel{Node{n}.Parent_Vessel}.Daughter_Node = n;
% end
%
% %%
% if (Vessel_Type == "Arteriole")
%     Diameter_Data = readtable('art_dia_gamma.csv');
%     Diameter_Data = table2array(Diameter_Data);
%     Length_data = readtable('art_len_linear.csv');
%     Length_data = table2array(Length_data);
% elseif (Vessel_Type == "Venule")
%     Diameter_Data = readtable('ven_dia_gamma.csv');
%     Diameter_Data = table2array(Diameter_Data);
%     Length_data = readtable('ven_len_gamma .csv');
%     Length_data = table2array(Length_data);
% else
%     disp('error sample diameter')
% end
% % Diameter_Data = load('Frequency_Diameter.mat');
% for v = 1:num_vessels
%     % Finds the number of generations of each vessel
%     Gen(v) =  Vessel{v}.Generation;
%
%     % Samples radii from distribution
%     Vessel{v}.Radius = Diameter_Data(randi(size(Diameter_Data,1),[1,1]))/2*10^-6;
% end
% num_generations = max(Gen);
%
% count = 0;
% while count < num_vessels-1
%     count = 0;
%     for v = 2:num_vessels
%         if Vessel{v}.Radius >= Vessel{Vessel{v}.Parent_Vessel}.Radius + 1e-6
%             temp = Vessel{v}.Radius;
%             Vessel{v}.Radius = Vessel{Vessel{v}.Parent_Vessel}.Radius;
%             Vessel{Vessel{v}.Parent_Vessel}.Radius = temp;
%         else
%             count = count + 1;
%         end
%     end
% end
% % Length_Data = load('Frequency_Length.mat');
% for v = 1:num_vessels
%     Vessel{v}.Length = Length_data(randi(size(Length_data,1),[1,1]))*10^-6;
%     %     Vessel{v}.Length = Sample_Length(Vessel_Type, Length_Data);
%     Vessel{v}.n_Length = Vessel{v}.Length/Vessel{v}.Radius;
% end
% % Number of vessels in each generation
% for g = 1:num_generations
%     v_gen(g) = sum(Gen == g);
% end
%
% % x space
% x_coord = linspace(0,1,num_generations);
%
% % y space
% y_coord = cell(num_generations-1,1);
% for g = 2:num_generations
%     y_coord{g-1} = linspace(0,1,v_gen(g)+2);
% end
%
% % Assigns xy coordinates to each node depending on generation number
% for n = 3:num_nodes
%     Node{n}.xy = [x_coord(Node{n}.Generation), y_coord{Node{n}.Generation-1}(1+v_gen(Node{n}.Generation))];
%     v_gen(Node{n}.Generation) = v_gen(Node{n}.Generation) - 1;
% end
%
% % Assigns xy coordiniates of the start and end of each vessel from the
% % nodal coordinates, and calculates length of the vessel.
% for v = 1:num_vessels
%     Vessel{v}.xy_Start = Node{Vessel{v}.Parent_Node}.xy;
%     Vessel{v}.xy_End = Node{Vessel{v}.Daughter_Node}.xy;
%     Vessel{v}.xy_Length = norm([Vessel{v}.xy_End,Vessel{v}.xy_Start]);
% end
%
% if Vessel_Type == "Arteriole"
%     for v = 1:num_vessels
%         Vessel{v}.Arterial_Generation = Vessel{v}.Generation;
%         Vessel{v}.Venule_Generation = [];
%         Vessel{v}.Type = "Arteriole";
%     end
% elseif Vessel_Type == "Venule"
%     for v = 1:num_vessels
%         Vessel{v}.Arterial_Generation = [];
%         Vessel{v}.Venous_Generation = Vessel{v}.Generation;
%         Vessel{v}.Type = "Venule";
%     end
% else
%     disp('error, vesseltype')
% end
%
%
%
% % for n = 1:num_nodes
% %     disp(Node{n})
% % end
% %
% % for v = 1:num_vessels
% %     disp(Vessel{v})
% % end
%
%
% end

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
Node{n}.Generation = Node{sprout}.Generation + 1;
Node{n}.Parent_Node = sprout;
Node{n}.Parent_Vessel = v;
Node{n}.num_Parent_Nodes =  1;

[Node] = Node_xy(Node,n);
end

function [Node] = Node_xy(Node,n)

% global Vessel_Type
Vessel_Type = "Arteriole";
Node{n}.Growth_Angle = Sample_Angle;
Length = Sample_Length(Vessel_Type);

points =  Node{Node{Node{n}.Parent_Node}.Parent_Node}.xy - Node{Node{n}.Parent_Node}.xy;
alpha = atan2d(points(2),points(1));
if (mod(n,2)== 1)
    theta = Node{n}.Growth_Angle + alpha;
    Node{n}.xy = Node{Node{n}.Parent_Node}.xy + [cosd(theta)*Length,sind(theta)*Length] ;
else
    theta = 360 - Node{n}.Growth_Angle + alpha;
    Node{n}.xy = Node{Node{n}.Parent_Node}.xy + [cosd(theta)*Length,sind(theta)*Length] ;
end



% sign = -1 + 2 * (mod(n,2)== 1);
% Node{n}.xy = Node{Node{n}.Parent_Node}.xy + sign*[cosd(theta)*Node{n}.Length,sind(theta)*Node{n}.Length ] ;
end

function Vessel = UpdateVessel(Vessel,Node,v,sprout)
% Updates nodal infomation when grown
Vessel{v}.ID = v;
Vessel{v} = InitaliseVessel(Vessel{v});
Vessel{v}.Parent_Node = sprout;
Vessel{v}.Parent_Vessel = Node{sprout}.Parent_Vessel;
Vessel{v}.Generation = Vessel{Vessel{v}.Parent_Vessel}.Generation + 1;

Vessel{v}.Daughter_Node = Node{sprout}.Daughter_Node;
end



%% Comparison Plots
