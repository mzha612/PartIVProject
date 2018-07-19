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
global d_P_0

x_step = 0.05;
EGL_Thick = [0.01:x_step:1.01]*10^-6;
d_P_0 = 14000;                    % Inital Pressure, Pa
d_P_inf = 1000;                   % Final Pressure, Pa
d_R = 15e-6;                       % Input vessel radius, m

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

num_bif = 9;                    % Number of bifurcatons for the binary

Vessel = [];
Node = [];
[Vessel, Node] = Create_Binary_Network(num_bif,d_R);
for i = 1:size(EGL_Thick,2)
    Vessel_Bin  = [];
    Node_Bin  = [];
    disp(i)
    epsilon = EGL_Thick(i);         % EGL thickness, um
    [Vessel_Bin,Node_Bin] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,0);
    Res_Bin(i) = (d_P_0 - d_P_inf)/Vessel_Bin{1}.Flow; % only works for binary

end
[Vessel_Bin_Ref,Node_Bin_Ref] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,1);
Res_Bin_Ref = (d_P_0 - d_P_inf)/Vessel_Bin_Ref{1}.Flow; % only works for binary

%% Rat Mesentery, Resistance vs EGL thickness.

Vessel = [];
Node = [];
[Vessel, Node] = Create_Real_Rat_Mesentery();
for i = 1:size(EGL_Thick,2)
    Vessel_Rat  = [];
    Node_Rat  = [];
    disp(i)
 
    epsilon = EGL_Thick(i);          % EGL thickness, um
    [Vessel_Rat,Node_Rat] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,0);
    inflow = 0;
    for n = 1:num_nodes
        if Node_Rat{n}.BC == 1
            if Vessel_Rat{Node_Rat{n}.Daughter_Vessel}.Flow > 0
                inflow = inflow + Vessel_Rat{Node_Rat{n}.Daughter_Vessel}.Flow;
            end
        elseif Node_Rat{n}.BC == -1
            if Vessel_Rat{Node_Rat{n}.Parent_Vessel}.Flow < 0
                inflow = inflow + Vessel_Rat{Node_Rat{n}.Parent_Vessel}.Flow;
            end            
        end
    end
    Res_Rat(i) = (d_P_0 - d_P_inf)/inflow;
    
end
for v = 1:num_vessels
    all_res_thick(v) = Vessel_Rat{v}.Resistance;
end


Vessel = [];
Node = [];
epsilon = EGL_Thick(1);
[Vessel, Node] = Create_Real_Rat_Mesentery();
[Vessel_Rat_Ref,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,1);
inflow = 0;
    for n = 1:num_nodes
        if Node_Rat{n}.BC == 1
            if Vessel_Rat_Ref{Node_Rat{n}.Daughter_Vessel}.Flow > 0
                inflow = inflow + Vessel_Rat_Ref{Node_Rat{n}.Daughter_Vessel}.Flow;
            end
        elseif Node_Rat{n}.BC == -1
            if Vessel_Rat_Ref{Node_Rat{n}.Parent_Vessel}.Flow < 0
                inflow = inflow + Vessel_Rat_Ref{Node_Rat{n}.Parent_Vessel}.Flow;
            end            
        end
    end
Res_Rat_Ref = (d_P_0 - d_P_inf)/inflow;
for v = 1:num_vessels
    all_res(v) = Vessel_Rat_Ref{v}.Resistance;
    disp(Vessel_Rat{v})
end

num_cap = 0;
num_art = 0;
num_ven = 0;
for v = 1:num_vessels
    num_cap = num_cap + 1*(Vessel{v}.Type == "Capillary");
    num_art = num_art + 1*(Vessel{v}.Type == "Arteriole");
    num_ven = num_ven + 1*(Vessel{v}.Type == "Venule");
end

%% Plot

Res_Bin = Res_Bin ./ Res_Bin_Ref;
Res_Rat = Res_Rat ./ Res_Rat_Ref;

figure
hold on
plot(EGL_Thick,Res_Bin)
plot(EGL_Thick,Res_Rat)
xlabel('EGL thickness (um)')
ylabel('Resistance')
legend('Binary','Rat')

%% Write to CSV
header = ['EGL_Thickness,','R_bin,','R_rat'];
fid = fopen('prelimresults3.csv','w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite('prelimresults3.csv', [EGL_Thick',Res_Bin',Res_Rat'], '-append')
