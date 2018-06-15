%% Results - Topo v2
%{
Runs the network calculation many times to obtain "results"
Currently compares each network with eachother
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

% num_vessels = 29;
% num_nodes = 30;
num_generations = 10;
num_capillaries = 180;

seed = 2;

dth = 0.05;
EGL_Thickness = 0.01:dth:1;
d_P_0 = 50;                    % Inital Pressure, MPa
d_P_inf = 10;                   % Final Pressure, MPa

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

%% RTB RTB
rng(seed)
[Vessel,Node] = CombineTree(@RTB,@RTB);
for i = 1:size(EGL_Thickness,2)
    disp(i),
    epsilon = EGL_Thickness(i);                  % EGL thickness, um
    [Vessel1,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_RTBRTB(i) = (d_P_0 - d_P_inf)/Vessel1{1}.Flow; % only works for binary
end

clear Vessel1

%% RSB RSB
rng(seed)
[Vessel, Node] = CombineTree(@RSB,@RSB);
for i = 1:size(EGL_Thickness,2)
    disp(i),
    epsilon = EGL_Thickness(i);                  % EGL thickness, um
    [Vessel2,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_RSBRSB(i) = (d_P_0 - d_P_inf)/Vessel2{1}.Flow; % only works for binary
end

clear Vessel2

%% RSB RTB
rng(seed)
[Vessel, Node] = CombineTree(@RSB,@RTB);
for i = 1:size(EGL_Thickness,2)
    disp(i),
    epsilon = EGL_Thickness(i);                  % EGL thickness, um
    [Vessel3,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_RSBRTB(i) = (d_P_0 - d_P_inf)/Vessel3{1}.Flow; % only works for binary
end

clear Vessel3

%% RTB RSB
rng(seed)
[Vessel, Node] = CombineTree(@RTB,@RSB);
for i = 1:size(EGL_Thickness,2)
    disp(i),
    epsilon = EGL_Thickness(i);                  % EGL thickness, um
    [Vessel4,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf);
    Res_RTBRSB(i) = (d_P_0 - d_P_inf)/Vessel4{1}.Flow; % only works for binary
end

clear Vessel4

%% Plot

% TODO, normalise with posueille flow
Res_RTBRTB = Res_RTBRTB / (Res_RTBRTB(1));
Res_RSBRSB = Res_RSBRSB / (Res_RSBRSB(1));
Res_RSBRTB = Res_RSBRTB / (Res_RSBRTB(1));
Res_RTBRSB = Res_RTBRSB / (Res_RTBRSB(1));

figure
hold on
plot(EGL_Thickness,Res_RTBRTB)
plot(EGL_Thickness,Res_RSBRSB)
plot(EGL_Thickness,Res_RSBRTB)
plot(EGL_Thickness,Res_RTBRSB)
xlabel('EGL thickness (um)')
ylabel('Resistance')
legend('RTB-RTB','RSB-RSB','RSB-RTB','RTB-RSB','Location','southeast')

%% Write CSV
header = ['EGL_Thickness,','Res_RTBRTB,','Res_RSBRSB','Res_RSBRTB','Res_RTBRSB'];
fid = fopen('prelimresultstopo.csv','w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite('prelimresultstopo.csv', [EGL_Thickness',Res_RTBRTB',Res_RSBRSB',Res_RSBRTB',Res_RTBRSB'], '-append')