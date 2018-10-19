function [Vessel,Node] = AR_preTerminal(Vessel,Node,Vessel_Type)
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

% figure
% hold on
% plot(x_min,y_min,x_max,y_max)
% for nn = 1:2
%     plot(Node{new_nodes(nn)}.xy,'bo')
% end
% plot([Vessel{new_vessels(1)}.xy_Start;Vessel{new_vessels(1)}.xy_End],'r')
preterminal = 0;
terminal = 0;
c = 0;
while (c < num_capillaries)
    
    
    c = 0;
    sprout = (randi([2,num_nodes],1)); % any segment (NODE) can bifurcate
    
    N_Temp = Node;
    V_Temp = Vessel;
    
    
    if rand(1) < 0.9
        isValid = 0;
        while ~isValid
            N_Temp = Node;
            V_Temp = Vessel;
    
            sprout = (randi([2,num_nodes],1)); % any segment (NODE) can bifurcate
            if Node{sprout}.Terminal
                %% Create and update nodal/vessel infomation
                N_Temp = UpdateNode(N_Temp,num_nodes + 1,num_vessels + 1,sprout);
                N_Temp = UpdateNode(N_Temp,num_nodes + 2,num_vessels + 2,sprout);
                
                N_Temp{sprout}.num_Daughter_Nodes =  2;
                N_Temp{sprout}.Terminal = false;
                N_Temp{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
                N_Temp{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
                
                V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + 1,sprout);
                V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + 2,sprout);
                V_Temp{num_vessels + 1}.Daughter_Node = N_Temp{sprout}.Daughter_Node(1);
                
                V_Temp{num_vessels + 2}.Daughter_Node = N_Temp{sprout}.Daughter_Node(2);
                V_Temp{num_vessels + 1}.xy_Start =  N_Temp{V_Temp{num_vessels+1}.Parent_Node}.xy;
                V_Temp{num_vessels + 1}.xy_End = N_Temp{V_Temp{num_vessels+1}.Daughter_Node}.xy;
                V_Temp{num_vessels + 2}.xy_Start =  N_Temp{V_Temp{num_vessels+2}.Parent_Node}.xy;
                V_Temp{num_vessels + 2}.xy_End = N_Temp{V_Temp{num_vessels+2}.Daughter_Node}.xy;
                
                isValid = Check_Valid2(V_Temp,N_Temp,sprout);
                
                if ~isValid
                else
                    Node = N_Temp ;
                    Vessel= V_Temp;
                    terminal = terminal + 1;
                end
            end
        end
    else
        isValid = 0;
        while ~isValid
            N_Temp = Node;
             V_Temp = Vessel;
    
            if num_nodes == 2
                break
            end
            sprout = (randi([2,num_nodes],1));
            
            
            if Node{sprout}.Terminal
                
            elseif Node{Node{sprout}.Daughter_Node(1)}.Terminal || Node{Node{sprout}.Daughter_Node(2)}.Terminal
                
                for d = 1:N_Temp{sprout}.num_Daughter_Nodes
                    % Store the previous daughter nodes and vessels
                    old_d_Nodes(d) = N_Temp{sprout}.Daughter_Node(d);
                    old_d_Vessels(d) = N_Temp{sprout}.Daughter_Vessel(d);
                    
                    % Update the daughter node/vessel with the new parent
                    % infomation
                    N_Temp{N_Temp{sprout}.Daughter_Node(d)}.Parent_Node = num_nodes + 2;
                    % +2 refers to all the duaghters splitting to the second bifurcation
                    % This doesn't matter unless it is affected by some thing that
                    % depends on the way it splits.
                    V_Temp{N_Temp{sprout}.Daughter_Vessel(d)}.Parent_Node = num_nodes + 2;
                    V_Temp{N_Temp{sprout}.Daughter_Vessel(d)}.Parent_Vessel = num_vessels + 2;
                    
                    N_Temp = UpdateNode(N_Temp,num_nodes + d,num_vessels + d,sprout);
                    V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + d,sprout);
                    
                    
                end
                V_Temp{num_vessels + 1}.Daughter_Node = num_nodes + 1;
                V_Temp{num_vessels + 2}.Daughter_Node = num_nodes + 2;
                %             % Create New Daughter Node
                %             N_Temp = UpdateNode(N_Temp,num_nodes + 1,num_vessels + 1,sprout);
                %             V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + 1,sprout);
                %
                %             % Create New terminal Node
                %             N_Temp = UpdateNode(N_Temp,num_nodes + 2,num_vessels + 2,sprout);
                %             V_Temp = UpdateVessel(V_Temp,N_Temp,num_vessels + 2,sprout);
                
                % update the sprouting node with new daughter infomation
                N_Temp{sprout}.Daughter_Node =  [num_nodes + 1, num_nodes + 2];
                N_Temp{sprout}.Daughter_Vessel = [num_vessels + 1, num_vessels + 2];
                
                % Update thew new nodes with infomation about the previous
                % daughters of the sprouting node ie (thier new daughters)
                N_Temp{num_nodes + 1}.Terminal = true;
                N_Temp{num_nodes + 2}.Terminal = false;
                % +2 here, hardcoded, branch adding onto from the second bifurcation.
                N_Temp{num_nodes + 1}.Daughter_Node = [];
                N_Temp{num_nodes + 2}.Daughter_Node = old_d_Nodes;
                
                N_Temp{num_nodes + 1}.Daughter_Vessel = [];
                N_Temp{num_nodes + 2}.Daughter_Vessel = old_d_Vessels;
                
                N_Temp{num_nodes + 1}.num_Daughter_Nodes = 0;
                N_Temp{num_nodes + 2}.num_Daughter_Nodes = 2;
                
                V_Temp{num_vessels + 1}.xy_Start =  N_Temp{sprout}.xy;
                V_Temp{num_vessels + 1}.xy_End = N_Temp{N_Temp{sprout}.Daughter_Node(1)}.xy;
                V_Temp{num_vessels + 2}.xy_Start =  N_Temp{sprout}.xy;
                V_Temp{num_vessels + 2}.xy_End = N_Temp{N_Temp{sprout}.Daughter_Node(2)}.xy;
                
                V_Temp{old_d_Vessels(1)}.xy_Start =  N_Temp{N_Temp{sprout}.Daughter_Node(2)}.xy;
                %             V_Temp{old_d_Vessels(1)}.xy_End = N_Temp{N_Temp{sprout}.Daughter_Node(2)}.xy;
                
                V_Temp{old_d_Vessels(2)}.xy_Start =  N_Temp{N_Temp{sprout}.Daughter_Node(2)}.xy;
                %             V_Temp{old_d_Vessels(2)}.xy_End = N_Temp{N_Temp{sprout}.Daughter_Node(2)}.xy;
                
                
                %             for n = 1:num_nodes + 2
                %                 disp(N_Temp{n})
                %             end
                %
                %
                %             for v = 1:num_vessels + 2
                %                 disp(V_Temp{v})
                %             end
                for d = 1:2
                    xy = N_Temp{N_Temp{sprout}.Daughter_Node(d)}.xy - N_Temp{sprout}.xy;
                    %                 V_Temp{num_vessels + d}.xy_Start
                    %                 V_Temp{num_vessels + d}.xy_End
                    %                xy = -xy;
                    % Recursively update each of the sprouting node's current
                    % daughters
                    [V_Temp,N_Temp] = UpdateDaughters(V_Temp,N_Temp,N_Temp{sprout}.Daughter_Node(d),xy);
                end
                %             for n = 1:num_nodes + 2
                %                 disp(N_Temp{n})
                %             end
                
                
                for v = 1:num_vessels + 2
                    V_Temp{v}.xy_Start = N_Temp{V_Temp{v}.Parent_Node}.xy;
                    V_Temp{v}.xy_End = N_Temp{V_Temp{v}.Daughter_Node}.xy;
                    %                 disp(V_Temp{v})
                    
                end
                
                isValid = 0;
                count = 0;
                angle = 0;
                while(~isValid && count < 20)
                    
                    [V_Temp,N_Temp] = Rotate_Daughters(V_Temp,N_Temp,sprout,sprout,angle);
                    
                    
                    isValid = Check_Valid_Segment(V_Temp,N_Temp,sprout,1);
                    angle = Sample_Angle;
                    count = count + 1;
                end
                
                
                
                
                if isValid
                    Node = N_Temp ;
                    Vessel= V_Temp;
                    preterminal = preterminal + 1;
                    
                    
                    
                    
                end
                %     elseif
                
            end
        end
    end
    
    
    %% Re-evaluate
    
    %     new_nodes = [num_nodes + 1,num_nodes + 2];
    %     new_vessels = [num_vessels + 1,num_vessels + 2];
    
    num_nodes = numel(Node);
    num_vessels = numel(Vessel);
    
    %     if isValid
    %         figure
    %         hold on
    %         for n = 1:num_nodes
    %             plot(Node{n}.xy(1),Node{n}.xy(2),'bo')
    %         end
    %         for v = 1:num_vessels
    %             plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}.xy_End(2)],'r')
    %             %         %         pause(0.01)
    %             %     end
    %         end
    %     end
    for n = 1:num_nodes
        if Node{n}.Terminal %termination criteria > num_capillaries
            c = c + 1;
        end
    end
    
end
terminal
preterminal


figure
hold on
for n = 1:num_nodes
    plot(Node{n}.xy(1),Node{n}.xy(2),'bo')
end
for v = 1:num_vessels
    plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}.xy_End(2)],'r')
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
% P_d =[4.687488329463876;2.465387856758793];% gamma

% P_l =   [1.002871132180170;2.888884270441188e+02];% gamma
% dist_l = 'Gamma';
dist_d = 'Gamma';
P_d = [6.2078;1.5853];
for v = 1:num_vessels
    % Finds the number of generations of each vessel
    Gen(v) =  Vessel{v}.Generation;
    lll = icdf(dist_d,rand(1),P_d(1),P_d(2));
    
%     while (lll > 35) || (lll < 2.5)
%         lll = icdf(dist_d,rand(1),P_d(1),P_d(2));
%     end

    Vessel{v}.Radius = lll + 2.5;
    
    %     if Vessel{v}.Radius >
    Vessel{v}.Radius = Vessel{v}.Radius/2*10^-6;
    % Samples radii from distribution
    %     Vessel{v}.Radius = Sample_Diameter(Vessel_Type)/2*10^-6;
end
num_generations = max(Gen);

count = 0;
while count < num_vessels-1
    count = 0;
    for v = 2:num_vessels
        if Vessel{v}.Radius >= Vessel{Vessel{v}.Parent_Vessel}.Radius + 2e-6
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
Node.Growth_Angle = [];
Node.Type = "Arteriole";

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
Vessel.Type = "Arteriole";

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

% No weights
%  P_d = [2.791278438217233;6.368312347248351];% gamma
% P_l = [1.198108031965868;2.425088293496245e+02];% gamma
% P_l = [1.392950112730204;1.216845915206235e+02];
% dist_l = 'Gamma';
%     dist_d = 'Gamma';
% Vessel_Type = "Arteriole";
Node{n}.Growth_Angle = Sample_Angle;
% Node{n}.Growth_Angle

% ll = (icdf(dist_l,rand(1),P_l(1),P_l(2)));
% while (ll > 1000) || (ll < 50)
%     ll= (icdf(dist_l,rand(1),P_l(1),P_l(2)));
% end
% Length = ll;

P = [-2.002022244691608e-06,0.002000000000000];
A = P(2)*1000/2;
P = P * A^-1;
a = P(1)/2;
b = P(2);
% cdf = @(y) (-b - sqrt(b^2-4*a.*(-y)))/(2*a);
% plot(cdf(y),y)
cdf = @(y) (-b + sqrt(b^2-4*a.*(-y)))/(2*a);

% for v = 1:num_vessels
%         Vessel{v}.Length = Length_data(randi(size(Length_data,1),[1,1]))*10^-6;
    r = rand(1);
    lll = cdf(r);
    if ~isreal(lll)
        lll = (-b + sqrt(b^2+4*a.*(-r)))/(2*a);
    end
    while lll < 50
        r = rand(1);
        lll = cdf(r);
        if ~isreal(lll)
        lll = (-b + sqrt(b^2+4*a.*(-r)))/(2*a);
        end
    end
    Length = lll;


% *10^-6;
% Length = Sample_Length(Vessel_Type);

points =  -Node{Node{Node{n}.Parent_Node}.Parent_Node}.xy + Node{Node{n}.Parent_Node}.xy;
alpha = atan2d(points(2),points(1));
if (mod(n,2)== 1)
    %     theta = Node{n}.Growth_Angle + alpha;
    theta = -180 + Node{n}.Growth_Angle + alpha;
    
    Node{n}.xy = Node{Node{n}.Parent_Node}.xy + [cosd(theta)*Length,sind(theta)*Length] ;
else
    %     theta = 360 - Node{n}.Growth_Angle + alpha;
    theta = -180 + Node{n}.Growth_Angle + alpha;
    theta = -theta + 180;
    Node{n}.xy = Node{Node{n}.Parent_Node}.xy + [cosd(theta)*Length,sind(theta)*Length] ;
end

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

function [Vessel,Node] = UpdateDaughters(Vessel,Node,n,xy)
% Recursively update each daughter with updated infomation when 'shifted
% down'
Vessel{Node{n}.Parent_Vessel}.Generation =  Vessel{Node{n}.Parent_Vessel}.Generation + 1;
Node{n}.Generation  = Node{n}.Generation + 1;

[Vessel,Node] = Translate(Vessel,Node,n,xy);
if ~Node{n}.Terminal
    for d = 1:Node{n}.num_Daughter_Nodes
        [Vessel,Node] = UpdateDaughters(Vessel,Node,Node{n}.Daughter_Node(d),xy);
    end
end

end

function [Vessel,Node] = Translate(Vessel,Node,n,xy)
for d = 1:Node{n}.num_Daughter_Nodes
    Node{Node{n}.Daughter_Node(d)}.xy = Node{Node{n}.Daughter_Node(d)}.xy + xy;
    %     Vessel{Node{n}.Daughter_Vessel(d)}.xy_Start = Vessel{Node{n}.Daughter_Vessel(d)}.xy_Start+ xy;
    %     Vessel{Node{n}.Daughter_Vessel(d)}.xy_End = Vessel{Node{n}.Daughter_Vessel(d)}.xy_End + xy ;
end
end

function isValid = Check_Valid_Segment(Vessel,Node,n,isValid)

if (~Node{n}.Terminal)
    check = Check_Valid2(Vessel,Node,n);
    isValid = isValid * check;
    if ~isValid
        return
    end
    for d = 1:Node{n}.num_Daughter_Nodes
        isValid = Check_Valid_Segment(Vessel,Node,Node{n}.Daughter_Node(d),isValid);
    end
end

end

function [Vessel,Node] = Rotate_Daughters(Vessel,Node,n,sprout,angle)
% is sprout og node.
% xy = Node{n}.xy - Node{sprout}.xy;
if (~Node{n}.Terminal)
    for d = 1:Node{n}.num_Daughter_Nodes
        % Shift
        xy =  Node{Node{n}.Daughter_Node(d)}.xy - Node{sprout}.xy;
        %         Node{Node{n}.Daughter_Node(d)}.xy = Node{Node{n}.Daughter_Node(1)}.xy - xy;
        % Rotate
        Node{Node{n}.Daughter_Node(d)}.xy = Node{sprout}.xy +...
            ([cosd(angle), -sind(angle);sind(angle), cosd(angle)] * xy')';
        
        
        
        Vessel{Node{n}.Daughter_Vessel(d)}.xy_Start = Node{n}.xy;
        Vessel{Node{n}.Daughter_Vessel(d)}.xy_End = Node{Node{n}.Daughter_Node(d)}.xy;
        [Vessel,Node] = Rotate_Daughters(Vessel,Node,Node{n}.Daughter_Node(d),sprout,angle);
    end
end

end
