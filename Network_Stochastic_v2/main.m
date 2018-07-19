%% Network Statistical v2
%{
The main script that runs....
%}

%{
Author = Michael Zhang
Date created = 11-06-18
%}
close all
clear
clc

%% Parameters
rng(0)

global num_vessels num_nodes num_generations num_capillaries

num_vessels = 29;
num_nodes = 30;
num_generations = 10;
num_capillaries = 4;
Vessel = cell(1,1);
Node = cell(1,1);

model = 4; 

%% Run Models

if model == 1
    %% Simple Test
    Vessel = cell(num_vessels,1);
    Node = cell(num_nodes,1);
    
    for n = 1:num_nodes
        Node{n}.xy = rand(1,2);
    end
    
    for v = 1:num_vessels
        Vessel{v}.xy_Start = Node{v}.xy;
        Vessel{v}.xy_End = Node{v+1}.xy;
        Vessel{v}.Generation = randi([1,10],1);
        Vessel{v}.Length = norm([Vessel{v}.xy_End,Vessel{v}.xy_Start]);
    end
    
    for n = 1:num_nodes
        disp(Node{n})
    end
    for v = 1:num_vessels
        disp(Vessel{v})
    end
elseif model == 2
    %% Random terminal branching
    [Vessel, Node] = RTB(Vessel,Node);
   
elseif model == 3
    %% Random segment branching
    [Vessel, Node] = RSB(Vessel,Node);
elseif model == 4
      [Vessel, Node] = CombineTree(@SSB,@SSB);
    for n = 1:num_nodes
    disp(Node{n})
end

for v = 1:num_vessels
    disp(Vessel{v})
end


end

%% Visualise
Visualise(Vessel,Node);

