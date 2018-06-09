function [Vessels] = Create_Hamster_Network(isPloting)
%% Create_Hamster_Network
%{
Function Creates the vessel infomation from a given datafile.
%}

%{
Inputs:
Outputs:
    Vessels: Structure, vessel segment infomation
%}

%{
Author = Michael Zhang
Date created = 03-06-18

Superceeds readdata.m.
%}

%% Read in file
global num_vessels num_nodes

num_vessels = 2035;

fid = fopen('networkgeometry.dat','r');
datacell = textscan(fid, '%f%f%f', 'HeaderLines', 13, 'Collect', 8140);
fclose(fid);
data = datacell{1};

% Seperate all the data
for i = 0:num_vessels-1
    xyz1(i + 1,:) = data(4*i + 1,:);
    xyz2(i + 1,:) = data(4*i + 2,:);
    radius(i + 1) = data(4*i + 3,1);
    in_out(i + 1) = data(4*i + 3,2);
    num_parents(i + 1) = data(4*i + 4,1);
    parent_1(i + 1) = data(4*i + 4,2);
    parent_2(i + 1) = data(4*i + 4,3);
end

if isPloting
    % Plot the networks
    close all
    figure
    hold on
    
    for i = 1:num_vessels
        if in_out(i) == 0
            c = 'r';
        else
            c = 'b';
        end
        plot3([xyz1(i,1); xyz2(i,1)],[xyz1(i,2); xyz2(i,2)],[xyz1(i,3); xyz2(i,3)],...
            'color',c,'LineWidth',radius(i)*10000)
        axis([-0.04,0.04,-0.04,0.04,0,0.08])
    end
    
end

%% Create the strucutres
Vessels = cell(1,num_vessels);

for i = 1:num_vessels
    Vessels{i}.ID = i;
    Vessels{i}.xyz1 = xyz1(i,:);
    Vessels{i}.xyz2 = xyz2(i,:);
    Vessels{i}.Radius = radius(i);
    Vessels{i}.Area = pi*radius(i).^2;
    Vessels{i}.Length = norm([xyz2(i,:)-xyz1(i,:)]);
    Vessels{i}.In_Out = in_out(i);
    
    Vessels{i}.Parent = [];
    Vessels{i}.Daughter = [];
    Vessels{i}.Parent_Node = [];
    Vessels{i}.Daughter_Node = [];
    Vessels{i}.BC = 0;
    Vessels{i}.num_Parents = num_parents(i);
    Vessels{i}.num_Daughter = 0;
    % Assign Parent info
    if num_parents(i) > 1
        Vessels{i}.Parent = [parent_1(i), parent_2(i)];
    elseif num_parents(i) == 1
        Vessels{i}.Parent = [parent_1(i)];
    else
        Vessels{i}.Parent = 0;
        Vessels{i}.Parent_Node = 0;
    end
    
end

num_nodes = 0;
has_Daugher = zeros(num_vessels,1);
no_Parent = zeros(num_vessels,1);
% num_parent = 0;

%% Assign Daugher info
for i = 1:num_vessels
    for j = 1:num_parents(i)
        Vessels{Vessels{i}.Parent(j)}.Daughter = [Vessels{Vessels{i}.Parent(j)}.Daughter, i];
    end
end

for i = 1:num_vessels
    Vessels{i}.num_Daughter = size(Vessels{i}.Daughter,1);
end

%         if Vessels{i}.Parent(j) == 0
%             num_parent = num_parent + 1;
%         end



%         if Vessels{i}.Parent(j) >= 1
%
%             Vessels{Vessels{i}.Parent(j)}.Daughter = [Vessels{Vessels{i}.Parent(j)}.Daughter, i];
%             if has_Daugher(Vessels{Vessels{i}.Parent(j)}.Daughter_Node) == 0
%                 num_nodes = num_nodes + 1;
%                 Vessels{Vessels{i}.Parent(j)}.Daughter_Node = num_nodes;
%             end
%             has_Daugher(i) = 1;
%
%             Vessels{i}.Parent_Node = Vessels{Vessels{i}.Parent(j)}.Daughter_Node;
%
%         end

%% Check if terminal end of vessel, ie. has no daugher
terminal = 0;
for i = 1:num_vessels
    
    if Vessels{i}.num_Parents == 0
        terminal = terminal + 1;
    end
    
end


for i = 1:num_vessels
    if isempty(Vessels{i}.Daughter_Node)
        
        num_nodes = num_nodes + 1;
        Vessels{i}.Daughter_Node = num_nodes;
    end
    
end

for i = 1:num_vessels
    if isempty(Vessels{i}.Parent_Node)
        
        %             num_nodes = num_nodes + 1;
        Vessels{i}.Parent_Node =  Vessels{Vessels{i}.Parent(1)}.Daughter_Node;
    end
end
d_R = max(radius);
for i = 1:num_vessels
    Vessels{i}.n_Radius = Vessels{i}.Radius/d_R;
    Vessels{i}.n_Length = Vessels{i}.Length/d_R;
    Vessels{i}.n_Area = Vessels{i}.n_Radius.^2;
end



for i = 1:num_vessels/3
    disp(Vessels{i})
end





no_Parent = num_parents == 0;
for i = find(no_Parent)
    Vessels{i}.BC = 1;
end

no_Daughter = ones(num_nodes,1);
% aaa = ;
% % bbb = aaa;
for i = 1:num_vessels
    no_Daughter(Vessels{i}.Daughter_Node) = 0;
end

% BC = 0;
no_Daughter = no_Daughter == 1;
for i = find(no_Daughter')
    Vessels{i}.BC = -1;
end

sum(no_Parent(:))

sum(has_Daugher(:))
terminal
end

