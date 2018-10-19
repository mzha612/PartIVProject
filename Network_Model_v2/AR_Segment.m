function [Vessel,Node] = AR_Segment(Vessel,Node,Vessel_Type)
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
% global Vessel_Type

Vessel_Type = "Arteriole";

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
Node{2}.Growth_Angle = 180;

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

c = 0; % Count of capillaries/or terminal segments, currently
while (c < num_capillaries)
    c = 0;
    sprout = (randi([2,num_nodes],1)); % any segment (NODE) can bifurcate
    Vessel_Temp = Vessel;
    Node_Temp = Node;
    %% Different segment branching scenarios
    if Node_Temp{sprout}.Terminal % If bifurcating from a terminal node, do the same as the RTB
        isValid = 0;
        %         while(~isValid)
        seed = [];
        num_seed = 0; % Number of possible sprouting points
        for n = 1:num_nodes % Only sprout from terminal nodes
            % Could be done with a constant list and pushing/popping
            if Node_Temp{n}.Terminal
                num_seed = num_seed + 1;
                seed(num_seed) = n; % Save the node number,
            end
        end
        %     TODO: have an updating list rather than make it every single time
        % Randomly choose a seed from the list
        sprout = seed(randi([1,num_seed],1));
        
        %% Create and update nodal/vessel infomation
        [Node_Temp,~,~] = UpdateNode(Node_Temp,num_nodes + 1,num_vessels + 1,sprout);
        [Node_Temp,~,~] = UpdateNode(Node_Temp,num_nodes + 2,num_vessels + 2,sprout);
        
        Node_Temp{sprout}.num_Daughter_Nodes =  2;
        Node_Temp{sprout}.Terminal = false;
        Node_Temp{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
        Node_Temp{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
        
        Vessel_Temp = UpdateVessel(Vessel_Temp,Node_Temp,num_vessels + 1,sprout);
        Vessel_Temp = UpdateVessel(Vessel_Temp,Node_Temp,num_vessels + 2,sprout);
        
        Vessel_Temp{num_vessels + 1}.xy_Start =  Node_Temp{Vessel_Temp{num_vessels+1}.Parent_Node}.xy;
        Vessel_Temp{num_vessels + 1}.xy_End = Node_Temp{Vessel_Temp{num_vessels+1}.Daughter_Node(1)}.xy;
        Vessel_Temp{num_vessels + 2}.xy_Start =  Node_Temp{Vessel_Temp{num_vessels+2}.Parent_Node}.xy;
        Vessel_Temp{num_vessels + 2}.xy_End = Node_Temp{Vessel_Temp{num_vessels+2}.Daughter_Node(2)}.xy;
        
        isValid = Check_Valid(Vessel_Temp,Node_Temp,sprout);
        if ~isValid
            %             Node{sprout}.Terminal = true;
            %             Node{num_nodes + 1} = [];
            %             Node{num_nodes + 2} = [];
            %             Vessel{num_vessels + 1} =[];
            %             Vessel{num_vessels + 2} =[];
        else
            Vessel = Vessel_Temp;
            Node = Node_Temp;
            num_nodes = numel(Node);
            num_vessels = numel(Vessel);
        end
    else % If bifurcating from a non terminal node
        isValid = 1;
        
        
        for d = 1:Node_Temp{sprout}.num_Daughter_Nodes
            % Store the previous daughter nodes and vessels
            old_d_Nodes(d) = Node_Temp{sprout}.Daughter_Node(d);
            old_d_Vessels(d) = Node_Temp{sprout}.Daughter_Vessel(d);
            
            % Update the daughter node/vessel with the new parent
            % infomation
            Node_Temp{Node_Temp{sprout}.Daughter_Node(d)}.Parent_Node = num_nodes + 2;
            % +2 refers to all the duaghters splitting to the second bifurcation
            % This doesn't matter unless it is affected by some thing that
            % depends on the way it splits.
            Vessel_Temp{Node_Temp{sprout}.Daughter_Vessel(d)}.Parent_Node = num_nodes + 2;
            Vessel_Temp{Node_Temp{sprout}.Daughter_Vessel(d)}.Parent_Vessel = num_vessels + 2;
            
            % Update the sprouting node with new daughter infomation
            Node_Temp{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
            Node_Temp{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
            
            % Update thew new nodes with infomation about the previous
            % daughters of the sprouting node ie (thier new daughters)
            Node_Temp{num_nodes + 1}.Terminal = true;
            Node_Temp{num_nodes + 2}.Terminal = false;
            % +2 here, hardcoded, branch adding onto from the second bifurcation.
            Node_Temp{num_nodes + 1}.Daughter_Node = [];
            Node_Temp{num_nodes + 2}.Daughter_Node = old_d_Nodes;
            
            Node_Temp{num_nodes + 1}.Daughter_Vessel = [];
            Node_Temp{num_nodes + 2}.Daughter_Vessel = old_d_Vessels;
            
            Node_Temp{num_nodes + 1}.num_Daughter_Nodes = 0;
            Node_Temp{num_nodes + 2}.num_Daughter_Nodes = 2;
            [Node_Temp,Angle,Length]  = UpdateNode(Node_Temp,num_nodes + d,num_vessels + d,sprout);
            Vessel_Temp = UpdateVessel(Vessel_Temp,Node_Temp,num_vessels + d,sprout);%
            
        end
        [Vessel_Temp,Node_Temp] = Transform(Vessel_Temp,Node_Temp,sprout,Angle,Length,sprout);
        %         [Vessel_Temp,Node_Temp,isValid] = UpdateDaughters(Vessel_Temp,Node_Temp,sprout);
        %             if ~isValid
        %                 break
        %             end
        % Recursively update each of the sprouting node's current
        % daughters
        %         end
        
        % Create New Daughter Node
        %             [Node_Temp,~,~]  = UpdateNode(Node_Temp,num_nodes + 1,num_vessels + 1,sprout);
        %             Vessel = UpdateVessel(Vessel,Node_Temp,num_vessels + 1,sprout);%
        % %             % Create New terminal Node
        %             [Node,~,~]  = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
        %             Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
        %
        
        
        
        %         for d = 1:Node_Temp{sprout}.num_Daughter_Nodes
        %             Vessel_Temp{num_vessels + 1}.xy_Start =  Node_Temp{Vessel_Temp{num_vessels+1}.Parent_Node}.xy;
        %             Vessel_Temp{num_vessels + 1}.xy_End = Node_Temp{Vessel_Temp{num_vessels+1}.Daughter_Node(1)}.xy;
        %             Vessel_Temp{num_vessels + 2}.xy_Start =  Node_Temp{Vessel_Temp{num_vessels+2}.Parent_Node}.xy;
        %             Vessel_Temp{num_vessels + 2}.xy_End = Node_Temp{Vessel_Temp{num_vessels+2}.Daughter_Node(2)}.xy;
        %         for n = 1:num_nodes+2
        %             disp(Node_Temp{n})
        %         end
        %
        %
        %         for v = 1:num_vessels+2
        %             disp(Vessel_Temp{v})
        %         end
        %
        %         %         end
        
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
        
        Vessel_Temp{num_vessels + 1}.xy_Start =  Node_Temp{Vessel_Temp{num_vessels+1}.Parent_Node}.xy;
        Vessel_Temp{num_vessels + 1}.xy_End = Node_Temp{Vessel_Temp{num_vessels+1}.Daughter_Node(1)}.xy;
        Vessel_Temp{num_vessels + 2}.xy_Start =  Node_Temp{Vessel_Temp{num_vessels+2}.Parent_Node}.xy;
        Vessel_Temp{num_vessels + 2}.xy_End = Node_Temp{Vessel_Temp{num_vessels+2}.Daughter_Node(2)}.xy;
        
        if ~isValid
            %
            %             Node{num_nodes + 1} = [];
            %             Node{num_nodes + 2} = [];
            %             Vessel{num_vessels + 1} =[];
            %             Vessel{num_vessels + 2} =[];
        else
            Vessel = Vessel_Temp;
            Node = Node_Temp;
            num_nodes = numel(Node_Temp);
            num_vessels = numel(Vessel_Temp);
        end
        if isValid
            new_nodes = [num_nodes,num_nodes-1];
            new_vessels = [num_vessels,num_vessels-1];
            figure
            hold on
            for nn = 1:num_nodes
                plot(Node{nn}.xy(1),Node{nn}.xy(2),'bo')
                %                     pause(1)
            end
            for nv = 1:num_vessels
                plot([Vessel{nv}.xy_Start(1),Vessel{nv}.xy_End(1)],[Vessel{nv}.xy_Start(2),Vessel{nv}.xy_End(2)],'r')
                %                     pause(1)
            end
            
            %         for nn = 1:size(new_nodes,2)
            %             plot(Node{new_nodes(nn)}.xy(1),Node{new_nodes(nn)}.xy(2),'bo')
            %             %                     pause(1)
            %         end
            %         for nv = 1:size(new_vessels,2)
            %             plot([Vessel{new_vessels(nv)}.xy_Start(1),Vessel{new_vessels(nv)}.xy_End(1)],[Vessel{new_vessels(nv)}.xy_Start(2),Vessel{new_vessels(nv)}.xy_End(2)],'r')
            %             %                     pause(1)
            %         end
            for n = 1:num_nodes
                disp(Node{n})
            end
            
            for v = 1:num_vessels
                disp(Vessel{v})
            end
            
        end
        
    end
    
    %% Re-evaluate
    %     num_nodes = numel(Node);
    %     num_vessels = numel(Vessel);
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

% Assigns vessel daughter infomation
for n = 2:num_nodes
    for d = 1:Node{n}.num_Daughter_Nodes
        Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel = [Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel, Node{n}.Daughter_Vessel(d)];
    end
    Vessel{Node{n}.Parent_Vessel}.Daughter_Node = n;
end

for v = 1:num_vessels
    % Finds the number of generations of each vessel
    Gen(v) =  Vessel{v}.Generation;
    % Samples radii from distribution
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
% [Node_Temp,Angle,Length]  = UpdateNode(Node_Temp,num_nodes + d,num_vessels + d,sprout);

function [Node,angle,length] = UpdateNode(Node,n,v,sprout)
% Updates nodal infomation when grown
Node{n}.ID = n;
Node{n} = InitaliseNode(Node{n});
Node{n}.Generation = Node{sprout}.Generation + 1;
Node{n}.Parent_Node = sprout;
Node{n}.Parent_Vessel = v;
Node{n}.num_Parent_Nodes =  1;

[Node,angle,length] = Node_xy(Node,n);
end

function [Node,theta,Length] = Node_xy(Node,n)

global Vessel_Type

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
%  Vessel_Temp = UpdateVessel(Vessel_Temp,Node_Temp,num_vessels + d,sprout);%

function Vessel = UpdateVessel(Vessel,Node,v,sprout)
% Updates nodal infomation when grown
Vessel{v}.ID = v;
Vessel{v} = InitaliseVessel(Vessel{v});
Vessel{v}.Parent_Node = sprout;
Vessel{v}.Parent_Vessel = Node{sprout}.Parent_Vessel;
Vessel{v}.Generation = Vessel{Vessel{v}.Parent_Vessel}.Generation + 1;

Vessel{v}.Daughter_Node = Node{sprout}.Daughter_Node;
end

% function [Vessel,Node] = RSB(Vessel,Node,Vessel_Type)
% %% Random Segment Branching v2
% %{
% Topological growth model
% Function that performs random segment growth.
% segment Growth occurs stocatically at any node.
% %}
% %{
% Inputs:
%     Vessel: Emtpy Cell that will hold vessel segment infomation
%     Node: Empty cell that will hold nodal infomation
%     Vessel_Type: String, indicating which distribution to sample from
% Outputs:
%     Vessel: Updated Structure,containing vessel segment infomation
%     Node: Updated Structure, containingnodal infomation
% %}
% %{
% Author = Michael Zhang
% Date created = 11-06-18
% %}
% %% Parameters
% global num_vessels num_nodes num_generations num_capillaries
%
% %% Initialise Inital nodes and vessel
% num_vessels = 1;
% num_nodes = 2;
%
% % Currently hardcoded inital segments
% % TODO: possible some arbitary starting points,
% Node{1}.ID = 1;
% Node{1} = InitaliseNode(Node{1});
% Node{1}.xy = [-0.2,0.5];
% Node{1}.num_Daughter_Nodes =  1;
% Node{1}.Terminal = false;
% Node{1}.BC = 1;
% Node{1}.Daughter_Vessel = 1;
% Node{1}.Daughter_Node = 2;
%
% Node{2}.ID = 2;
% Node{2} = InitaliseNode(Node{2});
% Node{2}.xy = [0,0.5];
% Node{2}.num_Daughter_Nodes =  0;
% Node{2}.Parent_Vessel = 1;
% Node{2}.Parent_Node = 1;
% Node{2}.Generation = 1;
%
% Vessel{1}.ID = 1;
% Vessel{1} = InitaliseVessel(Vessel{1});
% Vessel{1}.Parent_Node = 1;
% Vessel{1}.Daughter_Node = 2;
% Vessel{1}.Parent_Vessel = 0;
%
% %% Growth
%
% c = 0; % Count of capillaries/or terminal segments, currently
% while (c < num_capillaries)
%     c = 0;
%     sprout = (randi([2,num_nodes],1)); % any segment (NODE) can bifurcate
%
%     %% Different segment branching scenarios
%     if Node{sprout}.Terminal % If bifurcating from a terminal node, do the same as the RTB
%         Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
%         Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
%
%         Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
%         Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
%
%         Node{sprout}.num_Daughter_Nodes = 2;
%         Node{sprout}.Terminal = false;
%
%         Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
%         Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
%     else % If bifurcating from a non terminal node
%         for d = 1:Node{sprout}.num_Daughter_Nodes
%             % Store the previous daughter nodes and vessels
%             old_d_Nodes(d) = Node{sprout}.Daughter_Node(d);
%             old_d_Vessels(d) = Node{sprout}.Daughter_Vessel(d);
%
%             % Update the daughter node/vessel with the new parent
%             % infomation
%             Node{Node{sprout}.Daughter_Node(d)}.Parent_Node = num_nodes + 2;
%             % +2 refers to all the duaghters splitting to the second bifurcation
%             % This doesn't matter unless it is affected by some thing that
%             % depends on the way it splits.
%             Vessel{Node{sprout}.Daughter_Vessel(d)}.Parent_Node = num_nodes + 2;
%             Vessel{Node{sprout}.Daughter_Vessel(d)}.Parent_Vessel = num_vessels + 2;
%
%             % Recursively update each of the sprouting node's current
%             % daughters
%             [Vessel_Temp,Node_Temp,isValid] = UpdateDaughters(Vessel,Node,Node{sprout}.Daughter_Node(d));
%         end
%
%         % Create New Daughter Node
%         Node = UpdateNode(Node,num_nodes + 1,num_vessels + 1,sprout);
%         Vessel = UpdateVessel(Vessel,Node,num_vessels + 1,sprout);
%
%         % Create New terminal Node
%         Node = UpdateNode(Node,num_nodes + 2,num_vessels + 2,sprout);
%         Vessel = UpdateVessel(Vessel,Node,num_vessels + 2,sprout);
%
%         % update the sprouting node with new daughter infomation
%         Node{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
%         Node{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
%
%         % Update thew new nodes with infomation about the previous
%         % daughters of the sprouting node ie (thier new daughters)
%         Node{num_nodes + 1}.Terminal = true;
%         Node{num_nodes + 2}.Terminal = false;
%         % +2 here, hardcoded, branch adding onto from the second bifurcation.
%         Node{num_nodes + 1}.Daughter_Node = [];
%         Node{num_nodes + 2}.Daughter_Node = old_d_Nodes;
%
%         Node{num_nodes + 1}.Daughter_Vessel = [];
%         Node{num_nodes + 2}.Daughter_Vessel = old_d_Vessels;
%
%         Node{num_nodes + 1}.num_Daughter_Nodes = 0;
%         Node{num_nodes + 2}.num_Daughter_Nodes = 2;
%     end
%
%     %% Re-evaluate
%     num_nodes = numel(Node);
%     num_vessels = numel(Vessel);
%     for n = 1:num_nodes
%         if Node{n}.Terminal %terminatioion criteria > num_capillaries
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
% % Assigns Vessel daughter infomation
% for n = 2:num_nodes
%     for d = 1:Node{n}.num_Daughter_Nodes
%         Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel = ...
%             [Vessel{Node{n}.Parent_Vessel}.Daughter_Vessel, Node{n}.Daughter_Vessel(d)];
%     end
%     Vessel{Node{n}.Parent_Vessel}.Daughter_Node = n;
% end

function [Vessel,Node] = Transform(Vessel,Node,n,theta,length,og_n)
% function that rotates and translates all the nodes and segments that
% are daughters of the node n.
% if (mod(n,2)== 1)
%     theta = Node{og_n}.Growth_Angle + angle;
% else
%     theta = 360 - Node{og_n}.Growth_Angle + angle;
% end

% angle = 180 - Node{og_n}.Growth_Angle;
angle = 0;
sprout_node_xy = Node{og_n}.xy;

% Translate to origin
Node{n}.xy = Node{n}.xy - sprout_node_xy;

% Rotate about origin,
Node{n}.xy  = ([cosd(angle),-sind(angle);sind(angle),cosd(angle)] * Node{n}.xy')'  ;
% Shift back

Node{n}.xy = Node{n}.xy + sprout_node_xy;

% Apply Translation

Node{n}.xy = Node{n}.xy + [cosd(theta)*length,sind(theta)*length];

% Vessel{Node{n}.Daughter_Vessel}.xy_Start = Node{n}.xy;
% Vessel{Node{n}.Parent_Vessel}.xy_Start = [cosd(angle)*length,sind(angle)*length];

% Recursively update each daughter with updated infomation when 'shifted
% down'
if ~Node{n}.Terminal
    for d = 1:Node{n}.num_Daughter_Nodes
        [Vessel,Node] = Transform(Vessel,Node,Node{n}.Daughter_Node(d),theta,length,og_n);
    end
end

end

%%
function [Vessel,Node,isValid] = UpdateDaughters(Vessel,Node,n)
% Recursively update each daughter with updated infomation when 'shifted
% down'
Vessel{Node{n}.Parent_Vessel}.Generation =  Vessel{Node{n}.Parent_Vessel}.Generation + 1;
Node{n}.Generation  = Node{n}.Generation + 1;
% [Vessel,Node] = Transform(Vessel,Node,n,angle,length,og_n);
[isValid] = Check_Valid(Vessel,Node,Node{n}.Parent_Node);
if isValid
    if ~Node{n}.Terminal
        for d = 1:Node{n}.num_Daughter_Nodes
            [Vessel,Node,isValid] = UpdateDaughters(Vessel,Node,Node{n}.Daughter_Node(d));
        end
    end
end
end


