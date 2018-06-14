%% Results - Network v4
%{
Runs the network calculation many times to obtain "results"
%}

%{
Author = Michael Zhang
Date created = 13-06-18
%}
clear
close all
clc

%% Binary Tree, Resistance vs EGL thickness.
global epsilon num_vessels num_nodes

x_step = 0.05;
EGL_Thick = 0.01:x_step:1;
d_P_0 = 50;                    % Inital Pressure, MPa
d_P_inf = 10;                   % Final Pressure, MPa
d_R = 15;                       % Input vessel radius, um

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

num_bif = 9;                    % Number of bifurcatons for the binary
for i = 1:size(EGL_Thick,2)
    Vessel = [];
    Node = [];
    disp(i)
    
    [Vessel, Node] = Create_Binary_Network(num_bif,d_R);
    epsilon = EGL_Thick(i);         % EGL thickness, um
    [Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_Bin(i) = (d_P_0 - d_P_inf)/Vessel{1}.Flow; % only works for binary
end

%% Rat Mesentery, Resistance vs EGL thickness.

for i = 1:size(EGL_Thick,2)
    Vessel = [];
    Node = [];
    disp(i)
    
    [Vessel, Node] = Create_Real_Rat_Mesentery();
    epsilon = EGL_Thick(i);          % EGL thickness, um
    [Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    
    inflow = 0;
    for n = 1:num_nodes
        if Node{n}.BC == 1
            inflow = inflow + Vessel{Node{n}.Daughter_Vessel}.Flow;
        end
    end
    Res_Rat(i) = (d_P_0 - d_P_inf)/inflow;
end

%% Plot

% TODO, normalise with posueille flow
Res_Bin = Res_Bin / max(Res_Bin);
Res_Rat = Res_Rat / max(Res_Rat);

figure
hold on
plot(EGL_Thick,Res_Bin)
plot(EGL_Thick,Res_Rat)
xlabel('EGL thickness (um)')
ylabel('Resistance')
legend('Binary','Rat')

%% Write to CSV
header = ['EGL_Thickness,','R_bin,','R_rat'];
fid = fopen('prelimresults.csv','w');
fprintf(fid,'%s\n',header)
fclose(fid);
dlmwrite('prelimresults.csv', [EGL_Thick',Res_Bin',Res_Rat'], '-append')
