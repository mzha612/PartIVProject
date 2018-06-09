function [Vessel,Node] = Solve_Network(Vessel,Node,BC_in,BC_out)
%% Solve_Network v4
%{
Function that uses the vessels and nodal infomation to solve a system of
linear equations.
%}

%{
Inputs:
    Vessel: Structure, vessel segment infomation
    Node: Structure, nodal infomation
BCs: Boundary condition pressures
Outputs:
    Vessel: Updated Structure, vessel segment infomation
    Node: Updated Structure, nodal infomation
%}

%{
Author = Michael Zhang
Date created = 09-06-18
%}

%% Parameters

global epsilon num_vessels num_nodes

%% Matrix Skeleton

Connectivity = zeros(num_nodes,num_nodes);         % Connectivity Matrix
A = Connectivity;                                  % Matrix to solve
b = zeros(num_nodes,1);                            % RHS of the matrix

%% Conectivity Matrix
% 1 if it is coming into node i, and -1 ifs its going from node i
for i = 1:num_nodes
    Connectivity(i,Node{i}.Parent_Node) = 1;
    Connectivity(i,Node{i}.Daughter_Node) = -1;
end

%% Resistance/Conductance

Resistance = zeros(num_vessels,1);
for v = 1:num_vessels
    Vessel{v}.h = 1 - (epsilon)/Vessel{v}.Radius;
%     Vessel{v}.n_Resistivity = PoiseuilleFlow(Vessel{v}.n_Radius);
    Vessel{v}.n_Resistivity = Resistance_Summets(Vessel{v}.h); % TODO if using different resistance maybe use a function handle
    Resistance(v) = Vessel{v}.n_Resistivity * Vessel{v}.n_Length;
    Vessel{v}.n_Resistance = Resistance(v);
    Vessel{v}.n_Conductance = 1/Resistance(v);
end

%% Build Matrix
% Builds up A from the connectivity matrix and the conductance values
for i = 1:num_nodes
    if Node{i}.BC ~= 0 % if the node is a boundary value
        A(i,i) = 1;
        b(i) = ((Node{i}.BC == 1) * BC_in) + ((Node{i}.BC == -1) * BC_out);
        % Assumes presure infomation is given at the boundaries
    else % If the node is not a boundary value
        for j =  1:num_nodes 
            if Connectivity(i,j) == 1 
                % for the mass conservation at node i with parent node j,
                % assign appropriate values to the matrix
                
                for k = 1:Node{j}.num_Daughter_Nodes
                    for l = 1:Node{i}.num_Parent_Nodes
                        % Checks which vessel is connected to node i,
                        if Node{j}.Daughter_Vessel(k) == Node{i}.Parent_Vessel(l)
                            index = k; 
                        end
                    end
                end
                
                A(i,j) = -Vessel{Node{j}.Daughter_Vessel(index)}.n_Conductance;
                A(i,i) = A(i,i) + Vessel{Node{j}.Daughter_Vessel(index)}.n_Conductance;
            elseif Connectivity(i,j) == -1
                % for the mass conservation at node i with daugher node j,
                % assign appropriate values to the matrix
                
                for k = 1:Node{j}.num_Parent_Nodes
                    for l = 1:Node{i}.num_Daughter_Nodes
                        % Checks which vessel is connected to node i,
                        if Node{j}.Parent_Vessel(k) == Node{i}.Daughter_Vessel(l)
                            index = k;
                        end
                    end
                end
                
                A(i,j) = -Vessel{Node{j}.Parent_Vessel(index)}.n_Conductance;
                A(i,i) = A(i,i) + Vessel{Node{j}.Parent_Vessel(index)}.n_Conductance;
            end
        end
    end
end

%% Solve

P = A\b;
% TODO maybe a different solvr is necessary for bigger networks, currently
% okey dokey.

%% Assign Flow and Pressure values

for n = 1:num_nodes
    Node{n}.n_Pressure = P(n);
end

for v = 1:num_vessels
    Vessel{v}.n_Pressure_In = Node{Vessel{v}.Parent_Node}.n_Pressure;
    Vessel{v}.n_Pressure_Out = Node{Vessel{v}.Daughter_Node}.n_Pressure;
    Vessel{v}.n_Flow = (Vessel{v}.n_Pressure_In - Vessel{v}.n_Pressure_Out) * Vessel{v}.n_Conductance;
end

end