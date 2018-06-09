%% Network Main v3
%{
The main script that is able to call the different networks, and solve
given the various parameters
%}
%{
Author = Michael Zhang
Date created = 02-06-18
%}
clear
close all
clc

%% Parameters
global epsilon num_vessels

d_P_0 = 40;                    % Inital Pressure, MPa
d_P_inf = 10;                   % Final Pressure, MPa
d_R = 40;                       % Input vessel radius, um
epsilon = .5;                  % EGL thickness, um

% epsilon = d_epsilon/d_R;       % EGL thickness, non-dimensionalised
n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

num_bif = 3;                    % Number of bifurcatons for the binary

%% Network Model
% network_type = 2;
network_type = 1;
if network_type == 1        % Create the Binary Network
    [Vessels, Nodes] = Create_Binary_Network(num_bif,d_R);
else                        % Create the Hamster Cheek network
    [Vessels] = Create_Hamster_Network(0);
end

%% Solve
tic
[Vessels,Nodes] = Solve_Network(Vessels,Nodes,n_P_0,n_P_inf);
toc
%% Redimensionalise
for i = 1:num_vessels
%     Vessels{i}.Radius =  Vessels{i}.n_Radius * d_R;
%     Vessels{i}.Length = Vessels{i}.n_Length * Vessels{i}.Radius;
%     Vessels{i}.Resistance = %TODO
    Vessels{i}.Flow = Vessels{i}.n_Flow * Vessels{i}.Radius;
    Vessels{i}.Pressure_In = Vessels{i}.n_Pressure_In * d_P_0;
    Vessels{i}.Pressure_Out = Vessels{i}.n_Pressure_Out * d_P_0;
    disp(Vessels{i})
end

R = (d_P_0 - d_P_inf)/Vessels{1}.Flow % only works for binary

if network_type ~= 1
    maxflow = max(u(:));
    maxpressure = max(P(:));

    figure
    hold on
    for i = 1:num_vessels
        %     if Vessels{i}.Flow > 1
        %         c = rgb(Vessels{i}.Flow);
        %     else
        %         c = rgb(Vessels{i}.Flow);
        %     end
        c = uint8(Vessels{i}.n_Pressure_In/maxpressure*255)*uint8([0,0,1]);
        plot3([Vessels{i}.xyz1(1); Vessels{i}.xyz2(1)],[Vessels{i}.xyz1(2); Vessels{i}.xyz2(2)],[Vessels{i}.xyz1(3); Vessels{i}.xyz2(3)],...
            'color',c,'LineWidth',Vessels{i}.n_Radius*2)
    end
end
 