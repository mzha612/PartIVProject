%% Network Main v4
%{
The main script that is able to call the different networks, and solve
given the various parameters
%}
%{
Author = Michael Zhang
Date created = 09-06-18
%}
clear
close all
clc

%% Parameters
global epsilon num_vessels isPlotting isPrinting num_nodes

isPlotting = true;             % Boolean if we are plotting values.
isPrinting = true;             % Boolean, if we are printing infomation to the console

d_P_0 = 40;                    % Inital Pressure, MPa
d_P_inf = 10;                   % Final Pressure, MPa
d_R = 40;                       % Input vessel radius, um
% epsilon = .5;                  % EGL thickness, um

% epsilon = d_epsilon/d_R;       % EGL thickness, non-dimensionalised
n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

num_bif = 9;                    % Number of bifurcatons for the binary

%% Network Model
network_type = 2;
network_type = 1;
if network_type == 1        % Create the Binary Network
    [Vessel, Node] = Create_Binary_Network(num_bif,d_R);
else
    [Vessel, Node] = Create_Real_Rat_Mesentery();
end

%% Solve
epsilon = .1;                  % EGL thickness, um
[Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);

%% Redimensionalise
for i = 1:num_vessels
    %     Vessels{i}.Radius =  Vessels{i}.n_Radius * d_R;
    %     Vessels{i}.Length = Vessels{i}.n_Length * Vessels{i}.Radius;
    %     Vessels{i}.Resistance = %TODO redimensionalise
    Vessel{i}.Pressure_In = Vessel{i}.n_Pressure_In * d_P_0;
    Vessel{i}.Pressure_Out = Vessel{i}.n_Pressure_Out * d_P_0;
end

%% Printing
% Prints all the vessel, and nodal infomation into the console.
if isPrinting
    for i = 1:num_nodes
        disp(Node{i});
    end
    for i = 1:num_vessels
        disp(Vessel{i});
    end
    if network_type == 2
        % Calculates the total flow in, then divides by the pressure
        % difference over the whole network
        % TODO not sure if this is entirely correct atm
        inflow = 0;
        for i = 1:num_nodes
            if Node{i}.BC == 1
                inflow = inflow + Vessel{Node{i}.Daughter_Vessel}.Flow;
            end
        end
        R = (d_P_0 - d_P_inf)/inflow
    else
        R = (d_P_0 - d_P_inf)/Vessel{1}.Flow % only works for binary
    end
end
