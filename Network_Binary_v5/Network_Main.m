%% Network Main v5
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

global d_P_0
d_P_0 = 13e3;                   % Inital Pressure, Pa
d_P_inf = 1000;                   % Final Pressure, Pa
d_R = 40e-6;                    % Input vessel radius, m
epsilon = 0.0001e-6;            % EGL thickness, m
% epsilon = 1e-6;                  % EGL thickness, m

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

num_bif = 9;                    % Number of bifurcatons for the binary

%% Network Model
network_type = 2;
% network_type = 1;
if network_type == 1        % Create the Binary Network
    [Vessel, Node] = Create_Binary_Network(num_bif,d_R);
else
    [Vessel, Node] = Create_Real_Rat_Mesentery();
end

%% Solve
[Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,0);

%% Printing
% Prints all the vessel, and nodal infomation into the console.
if isPrinting
    for n = 1:num_nodes
        disp(Node{n});
    end
    for v = 1:num_vessels
        disp(Vessel{v});
    end
    if network_type == 2
        % Calculates the total flow in, then divides by the pressure
        % difference over the whole network
        % TODO not sure if this is entirely correct atm
        inflow = 0;
        for n = 1:num_nodes
            if Node{n}.BC == 1
                inflow = inflow + Vessel{Node{n}.Daughter_Vessel}.Flow;
            end
        end
        R = (d_P_0 - d_P_inf)/inflow
    else
        R = (d_P_0 - d_P_inf)/Vessel{1}.Flow % only works for binary
    end
end
