function [Vessel,Node] = Solve_Network(Vessel,Node,BC_in,BC_out,R_type)
%% Solve_Network - Network_Model_v1
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
Date created = 17-07-18
%}

%% Parameters

global epsilon num_vessels num_nodes
global d_P_0
global cs
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
phi_f = 0.99; % Fluid phase fraction
mu_f = 10^-3; % viscosity of water, Pa.s
K = 10^10;

for v = 1:num_vessels
    Vessel{v}.h = 1 - (epsilon)/Vessel{v}.Radius;
%     v
%     chi = K*(Vessel{v}.Radius)^2/(phi_f * mu_f);
%     Vessel{v}.Resistance =  Vessel{v}.n_Length  / Vessel{v}.Radius.^2 * mu_f *Resistance_Summets(Vessel{v}.h,chi); % Pa.s.m-3.m
%     Vessel{v}.Resistance =  Vessel{v}.n_Length  / Vessel{v}.Radius.^2 * mu_f *Resistance_Summets_Quad(Vessel{v},cs); % Pa.s.m-3.m
    Vessel{v}.Resistance =  Vessel{v}.n_Length  / Vessel{v}.Radius.^2 * mu_f *Resistance_Varicose(Vessel{v},cs); % Pa.s.m-3.m
    Vessel{v}.Resistance_Poiseuille = Vessel{v}.Length*PoiseuilleFlow(Vessel{v}.Radius);
%     Vessel{v}.Resistance_Poiseuille_effective = Vessel{v}.Length*PoiseuilleFlow(Vessel{v}.Radius-epsilon);
    Vessel{v}.Resistance_Poiseuille_effective =  Vessel{v}.n_Length  / Vessel{v}.Radius.^2 * mu_f *Resistance_Summets_Quad(Vessel{v},cs);
    if (R_type == 0)
        Vessel{v}.Conductance = 1/Vessel{v}.Resistance;
    elseif (R_type == 1)
        Vessel{v}.Conductance = 1/Vessel{v}.Resistance_Poiseuille;
    elseif (R_type == 2)
        Vessel{v}.Conductance = 1/Vessel{v}.Resistance_Poiseuille_effective;
    else
        disp('error in R_type')
    end
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
                
                A(i,j) = -Vessel{Node{j}.Daughter_Vessel(index)}.Conductance;
                A(i,i) = A(i,i) + Vessel{Node{j}.Daughter_Vessel(index)}.Conductance;
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
                
                A(i,j) = -Vessel{Node{j}.Parent_Vessel(index)}.Conductance;
                A(i,i) = A(i,i) + Vessel{Node{j}.Parent_Vessel(index)}.Conductance;
            end
        end
    end
end

%% Solve

P = A\b;
% TODO maybe a different solvr is necessary for bigger networks, currently
% okie dokie.

%% Assign Pressure values

for n = 1:num_nodes
    Node{n}.n_Pressure = P(n);
end

for v = 1:num_vessels
    Vessel{v}.n_Pressure_In = Node{Vessel{v}.Parent_Node}.n_Pressure;
    Vessel{v}.n_Pressure_Out = Node{Vessel{v}.Daughter_Node}.n_Pressure;
end

%% Redimensionalise
for v = 1:num_vessels
    Vessel{v}.Pressure_In = Vessel{v}.n_Pressure_In * d_P_0;
    Vessel{v}.Pressure_Out = Vessel{v}.n_Pressure_Out * d_P_0;
    Vessel{v}.Flow = (Vessel{v}.Pressure_In - Vessel{v}.Pressure_Out) * Vessel{v}.Conductance;
%     if Vessel{v}.Flow < 0
%         disp(Vessel{v});
%     end
end

end