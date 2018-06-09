function [Vessels,Nodes] = Solve_Network(Vessels,Nodes,BC_in,BC_out)
%% Solve_Network v3
%{
Function that uses the vessels and  nodal infomation to solve a system of
linear equations.
%}

%{
Inputs:
    Vessels: Structure, vessel segment infomation
    Nodes: Structure, nodal infomation
BCs: Boundary condition pressures
Outputs:
    Q = Flows
    P = Pressures
    Vessels: Updated Structure, vessel segment infomation
%}

%{
Author = Michael Zhang
Date created = 02-06-18

Derived from Networks_v1, and v2.
%}

%% Formulate Matrix Skeleton

global epsilon num_vessels num_nodes


% mass_C = zeros(num_nodes,num_vessels + num_nodes);  % Mass conservation
% flow_C = zeros(num_vessels,num_nodes);              % Pressure driven flow
%
% % Creates mass conservation submatrix, in =  out
% for i = 1:num_vessels
%     mass_C(Vessels{i}.Parent_Node+1,i) = -1;
%     mass_C(Vessels{i}.Daughter_Node+1,i) = 1;
% end
%
% for i = 1:num_vessels
%     % zero a mass submatrix row, for the boundary condition equation
%     if (Vessels{i}.BC == -1) % Output
%         mass_C(Vessels{i}.Daughter_Node + 1,:) = 0;
%         mass_C(Vessels{i}.Daughter_Node + 1,Vessels{i}.Daughter_Node+num_vessels + 1) = 1;
%         b(Vessels{i}.Daughter_Node  + 1) = BC_out;
%     elseif (Vessels{i}.BC == 1) % Input
%         mass_C(Vessels{i}.Parent_Node + 1,:) = 0;
%         mass_C(Vessels{i}.Parent_Node + 1,Vessels{i}.Parent_Node+num_vessels+1) = 1;
%         b(Vessels{i}.Parent_Node  + 1) = BC_in;
%     end
% end
%
% % Creates pressure driven flow submatrix, parent - daughter
% for i = 1:num_vessels
%     flow_C(i,[Vessels{i}.Parent_Node]+1) = 1;
%     flow_C(i,[Vessels{i}.Daughter_Node]+1) = -1;
% end

%% Conectivity Matrix

C = zeros(num_nodes,num_nodes);
A = C;
b = zeros(num_nodes,1);               % RHS of the matrix

for i = 1:num_nodes
    C(i,Nodes{i}.Parent_Node) = 1;
    C(i,Nodes{i}.Daughter_Node) = -1;
end

%% Apply Matrix Values
% Find the radius for each vessel for the mass conservation
% for i = 1:num_vessels
%     Radius(i,1) = Vessels{i}.c_Radius;
% end

% Resistance
for i = 1:num_vessels
    h = 1 - (epsilon)/Vessels{i}.Radius;
    Vessels{i}.h = h;
    Vessels{i}.n_Resistivity = Resistance_Summets(h);
    %     TODO Where do i put the resistance thing.
    Resistance(i) = Vessels{i}.n_Resistivity * Vessels{i}.n_Length;
    Vessels{i}.n_Resistance = Resistance(i);
    Vessels{i}.n_Conductance = 1/Resistance(i);
end


for i = 1:num_nodes
    if Nodes{i}.BC ~= 0
        A(i,i) = 1;
        b(i) = (Nodes{i}.BC == 1) * BC_in + (Nodes{i}.BC == -1) * BC_out;
    else
        
        
        for j =  1:num_nodes
            if C(i,j) == 1
                for k = 1:Nodes{j}.num_Daughter_Nodes
                    for l = 1:Nodes{i}.num_Parent_Nodes
                        if Nodes{j}.Daughter_Vessel(k) == Nodes{i}.Parent_Vessel(l)
                            i1 = k;
                            %                     i2 = l;
                        end
                    end
                end
                
                A(i,j) = -Vessels{Nodes{j}.Daughter_Vessel(i1)}.n_Conductance;
                A(i,i) = A(i,i) + Vessels{Nodes{j}.Daughter_Vessel(i1)}.n_Conductance;
            elseif C(i,j) == -1
                for k = 1:Nodes{j}.num_Parent_Nodes
                    for l = 1:Nodes{i}.num_Daughter_Nodes
                        
                        if Nodes{j}.Parent_Vessel(k) == Nodes{i}.Daughter_Vessel(l)
                            i1 = k;
                            %                     i2 = l;
                        end
                    end
                end
                A(i,j) = -Vessels{Nodes{j}.Parent_Vessel(i1)}.n_Conductance;
                A(i,i) = A(i,i) + Vessels{Nodes{j}.Parent_Vessel(i1)}.n_Conductance;
            end
        end
    end
end


% Am = mass_C.*(ones(num_nodes,1)*[Radius;ones(num_nodes,1)]');
% AQ = flow_C.*((1./Resistance)*ones(num_nodes,1)');
%
% A = [Am; eye(num_vessels),-AQ];

% Solve
P = A\b;



%% Assign Flow and Pressure values

for i = 1:num_nodes
    Nodes{i}.n_Pressure = P(i);
%     disp(Nodes{i})
end

for i = 1:num_vessels
    Vessels{i}.n_Pressure_In = Nodes{Vessels{i}.Parent_Node}.n_Pressure;
    
    Vessels{i}.n_Pressure_Out = Nodes{Vessels{i}.Daughter_Node}.n_Pressure;
    
    Vessels{i}.n_Flow = (Vessels{i}.n_Pressure_In - Vessels{i}.n_Pressure_Out) * Vessels{i}.n_Conductance;
%     disp(Vessels{i})
end


% for i = 1:num_vessels
%     Vessels{i}.n_Flow = Q(i);
%     Vessels{i}.n_Pressure_In= P(Vessels{i}.Parent_Node+1);
%     Vessels{i}.n_Pressure_Out = P(Vessels{i}.Daughter_Node+1);
% end

end