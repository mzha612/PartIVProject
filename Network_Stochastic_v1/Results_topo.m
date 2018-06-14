%% Results - Topo v4
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
global num_generations num_capillaries

num_vessels = 29;
num_nodes = 30;
num_generations = 10;
num_capillaries = 180;

dth = 0.05;
EGL_Thickness = 0.01:dth:1;
d_P_0 = 50;                    % Inital Pressure, MPa
d_P_inf = 10;                   % Final Pressure, MPa

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

%% RTB
for i = 1:size(EGL_Thickness,2)
    Vessel = [];
    Node = [];
    disp(i)
    rng(0)

    [Vessel,Node] = RTB(Vessel,Node);
    epsilon = EGL_Thickness(i);                  % EGL thickness, um
    [Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_RTB(i) = (d_P_0 - d_P_inf)/Vessel{1}.Flow; % only works for binary
end

%% RSB

for i = 1:size(EGL_Thickness,2)
    Vessel = [];
    Node = [];
    disp(i)
    rng(0)
    [Vessel,Node] = RSB(Vessel,Node);
    epsilon = EGL_Thickness(i);                  % EGL thickness, um
    [Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_RSB(i) = (d_P_0 - d_P_inf)/Vessel{1}.Flow; % only works for binary
end

%% Plot

% TODO, normalise with posueille flow
Res_RTB = Res_RTB / max(Res_RTB);
Res_RSB = Res_RSB / max(Res_RSB);

figure
hold on 
plot(EGL_Thickness,Res_RTB)
plot(EGL_Thickness,Res_RSB)
xlabel('EGL thickness (um)')
ylabel('Resistance')
legend('RTB','RSB')

% Write CSV
header = ['EGL_Thickness,','R_RTB,','R_RSB'];
fid = fopen('prelimresultstopo.csv','w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite('prelimresultstopo.csv', [EGL_Thickness',Res_RTB',Res_RSB'], '-append')