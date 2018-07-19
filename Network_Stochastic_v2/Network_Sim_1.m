%% Results - Topo v2
%{
Runs many iterations of the network calculation for each EGL thickness
%}

%{
Author = Michael Zhang
Date created = 13-06-18
%}
clear
close all
clc

%% Parameters
global epsilon num_vessels num_nodes
global num_generations num_capillaries 
global d_P_0

num_capillaries = 180;
max_generations = num_capillaries * 2 + 5;

num_iterations = 25;
num_network = 4;

seed = 3; % Random number seed
rng(seed)

% Range of egl thicknesses to run for
e_step = .05;
e_range = 0.01:e_step:1.01;
e_range = e_range *10^-6; % m

d_P_0 = 13000;                    % Inital Pressure, Pa
d_P_inf = 1000;                   % Final Pressure, Pa

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;                      % Input vessel radius, non-dimensionalised

Trees = {@RTB,@RTB;...
    @RSB,@RSB;...
    @RSB,@RTB;...
    @RTB,@RSB;
    @SSB,@SSB;};

% mmm = ['Res_RTBRTB,','Res_RSBRSB','Res_RSBRTB','Res_RTBRSB'];
% header = [['Iteration_',mmm(num_network)]];
filenames{1} = ['network_gen_dia_',num2str(num_network),'.csv'];
filenames{2} = ['network_gen_len_',num2str(num_network),'.csv'];
filenames{3} = ['network_freq_cap_',num2str(num_network),'.csv'];
fid = fopen(filenames{1},'w');
% fprintf(fid,'%s\n',header);
fclose(fid);

fid = fopen(filenames{2},'w');
% fprintf(fid,'%s\n',header);
fclose(fid);

fid = fopen(filenames{3},'w');
% fprintf(fid,'%s\n',header);
fclose(fid);


%%
R = zeros(size(e_range,2),num_iterations);
R_Ref = zeros(1,num_iterations);

for j = 1:num_iterations
    tic
    Vessel = [];
    Node = [];
    
    [Vessel,Node] = CombineTree(Trees{num_network,:}); % Create the network
    for i = 1:size(e_range,2)
        %         disp(i),
        epsilon = e_range(i);    % EGL thickness, m
        [Vessel_Temp,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,0);
        R(i,j) = (d_P_0 - d_P_inf)/Vessel_Temp{1}.Flow;
    end
    [Vessel_Ref,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,1);
    R_Ref(j) = (d_P_0 - d_P_inf)/Vessel_Ref{1}.Flow;
    
    [Vessel_Ref_thick,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,2);
    R_Ref_thick(j) = (d_P_0 - d_P_inf)/Vessel_Ref_thick{1}.Flow;
    
    obtainGenInfo(Vessel,filenames)
    clear Vessel_Ref
    toc
    
end
for v = 1:num_vessels
    disp(Vessel{v});
end
%% Plot
R = R ./ R_Ref;

err = std(R,[],2);
% /sqrt(num_iterations);

R = mean(R,2);

figure
hold on
errorbar(e_range*10^6,R,err,'color',[242,113,39]/255,'LineWidth',2)
% plot(e_range*10^6,
xlabel('EGL thickness (\mum)')
ylabel('Resistance (%)')
nnn = ([{'RTB-RTB'},{'RSB-RSB'},{'RSB-RTB'},{'RTB-RSB'},{'SSB-SSB'}]);
legend(nnn(num_network),'Location','southeast')

legend("Base Model", "Rat Mesentery Model" ,'Location','southeast')

%% Write CSV
mmm = [{'Res_RTBRTB,'},{'Res_RSBRSB,'},{'Res_RSBRTB,'},{'Res_RTBRSB,'},{'Res_SSBSSB'}];
header = ['EGL_Thickness,',mmm(num_network),'std'];
filename = ['network_resistance_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s',string(header));
fprintf(fid,'\n');
fclose(fid);
dlmwrite(filename, [e_range',R,err], '-append')

%% Read Intermediate data for generation info
real = readtable('Digitised_Data.csv');
real = table2array(real);

X_Ven_Gen	=real(:,1);
Y_Num_Cap_Ven	=real(:,2);
X_Art_Gen	=real(:,3);
Y_Num_Cap_Art	=real(:,4);
Length_Freq	=real(:,5);
Art_Len_Freq	=real(:,6);
Ven_Len_Freq	=real(:,7);
Cap_Len_Freq	=real(:,8);
Diameter_Freq =real(:,9);
Art_Dia_Freq=real(:,10);
Ven_Dia_Freq	=real(:,11);
Cap_Dia_Freq	=real(:,12);
Gen_Art	=real(:,13);
Art_Gen_Dia	=real(:,14)*10^-6;
Art_Gen_Len	=real(:,15)*10^-6;
Gen_Ven	=real(:,16);
Ven_Gen_Dia	=real(:,17)*10^-6;
Ven_Gen_Len	=real(:,18)*10^-6;
Gen_Cap	=real(:,19);
Cap_Gen_Dia	=real(:,20)*10^-6;
Cap_Gen_Len	=real(:,21)*10^-6;


dia_data = csvread(filenames{1});
len_data = csvread(filenames{2});
freq_cap = csvread(filenames{3});

n = sum(~isnan(freq_cap(1:2:end,:)) & freq_cap(1:2:end,:)~= 0);

freq_art_cap_err = std(freq_cap(1:2:end,:),'omitnan')./sqrt(sum(~isnan(freq_cap(1:2:end,:)) & freq_cap(1:2:end,:)~= 0));
freq_ven_cap_err = std(freq_cap(2:2:end,:),'omitnan')./sqrt(sum(~isnan(freq_cap(2:2:end,:)) & freq_cap(2:2:end,:)~= 0));


freq_art_cap =  nanmean(freq_cap(1:2:end,:),1);
freq_ven_cap =  nanmean(freq_cap(2:2:end,:),1);

header = ['Art,','ven,','std,'];
filename = ['network_gen_avg_freq_cap',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [freq_art_cap;freq_ven_cap;freq_art_cap_err;freq_ven_cap_err]', '-append')

figure
hold on
errorbar((1:30),freq_art_cap(1:30),freq_art_cap_err(1:30),'color',[242,113,39]/255,'LineWidth',2)
% errorbar( 1:size(freq_cap,2),freq_ven_cap,freq_ven_cap_err)
plot(X_Art_Gen	,Y_Num_Cap_Art,'-o','color',[1,64,52]/255,'LineWidth',2)
% ,X_Ven_Gen,Y_Num_Cap_Ven,'-o'
xlabel('Generation Number')
ylabel('Number of Capillaries')
% legend('Arteriole','Venule','Arteriole','Venule')
legend('Model','Data')
len_art_err = std( len_data(1:3:end,:),'omitnan')./sqrt(sum(~isnan(len_data(1:3:end,:)) & len_data(1:3:end,:)~= 0))* 10^6;

len_cap_err = std( len_data(2:3:end,:),'omitnan');
len_ven_err = std( len_data(3:3:end,:),'omitnan');
len_ven_err(40:end) = NaN;

len_art =  nanmean(len_data(1:3:end,:)) * 10^6;
len_cap =  nanmean(len_data(2:3:end,:))* 10^6;
len_ven =  nanmean(len_data(3:3:end,:))* 10^6;
len_ven(40:end) = NaN;

header = ['Art,','cap,','ven,','std,'];
filename = ['network_gen_avg_len_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [len_art;len_cap;len_ven;len_art_err;len_cap_err;len_ven_err]', '-append')

figure
hold on
errorbar((1:40),len_art(1:40),len_art_err(1:40),'color',[242,113,39]/255,'LineWidth',2)
% errorbar((1:40),len_cap(1:40),len_cap_err(1:40))
% errorbar(1:size(len_data,2),len_ven,len_ven_err)
plot(Gen_Art,Art_Gen_Len* 10^6,'-o','color',[1,64,52]/255,'LineWidth',2)
% ,Gen_Cap,Cap_Gen_Len,'-o')
% ,Gen_Ven,Ven_Gen_Len,'-o'
% Gen_Art	=real(:,13);
% Art_Gen_Dia	=real(:,14);
% Art_Gen_Len	=real(:,15);
% Gen_Ven	=real(:,16);
% Ven_Gen_Dia	=real(:,17);
% Ven_Gen_Len	=real(:,18);
% Gen_Cap	=real(:,19);
% Cap_Gen_Dia	=real(:,20);
% Cap_Gen_Len	=real(:,21);
xlabel('Generation Number')
ylabel('Mean Length (\mum)')
% legend('Arteriole','Capillary','Venule','Arteriole','Capillary','Venule')
legend('Model','Data')



dia_art_err = std(dia_data(1:3:end,:),'omitnan')./sqrt(sum(~isnan(dia_data(1:3:end,:)) & dia_data(1:3:end,:)~= 0))* 10^6;

dia_cap_err = std(dia_data(2:3:end,:),'omitnan');
dia_ven_err = std(dia_data(3:3:end,:),'omitnan');
dia_ven_err(40:end) = NaN;



dia_art =  nanmean(dia_data(1:3:end,:))* 10^6;
dia_cap =  nanmean(dia_data(2:3:end,:));
dia_ven =  nanmean(dia_data(3:3:end,:));
% dia_ven(freq_ven_cap < 1) = NaN;
dia_ven(40:end) = NaN;

header = ['Art,','cap,','ven,','std,'];
filename = ['network_gen_avg_dia_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [dia_art;dia_cap;dia_ven;dia_art_err;dia_cap_err;dia_ven_err]', '-append')

figure
hold on
% errorbar(1:size(dia_data,2),dia_art,dia_art_err)
errorbar((1:40),dia_art(1:40),dia_art_err(1:40),'color',[242,113,39]/255,'LineWidth',2)
% errorbar(1:size(dia_data,2),dia_cap,dia_cap_err)
% errorbar(1:size(dia_data,2),dia_ven,dia_ven_err)
plot(Gen_Art,Art_Gen_Dia* 10^6,'-o','color',[1,64,52]/255,'LineWidth',2)
% Gen_Cap,Cap_Gen_Dia,'-o')
% ,Gen_Ven,Ven_Gen_Dia,'-o'
xlabel('Generation Number')
ylabel('Mean Diameter (\mum)')
% legend('Arteriole','Capillary','Venule','Arteriole','Capillary','Venule')
% legend('Arteriole','Capillary','Arteriole','Capillary')
legend('Model','Data')

%%
function [] = obtainGenInfo(Vessel,filenames)
% Function that finds the diameter and (LEngths TODO) and ports against
% generation
global num_generations num_vessels

freq_cap = zeros(2,num_generations);
diameters = cell(3,num_generations);
lengths = cell(3,num_generations);
for v = 1:num_vessels
    if Vessel{v}.Type == "Arteriole"
        diameters{1,Vessel{v}.Arterial_Generation} = [diameters{1,Vessel{v}.Arterial_Generation}, Vessel{v}.Radius * 2];
        lengths{1,Vessel{v}.Arterial_Generation} = [lengths{1,Vessel{v}.Arterial_Generation}, Vessel{v}.Length];
        
    elseif Vessel{v}.Type == "Capillary"
        i = Vessel{v}.Arterial_Generation;
        j = Vessel{v}.Venous_Generation;
        
        freq_cap(1,i) = freq_cap(1,i) + 1;
        freq_cap(2,j) = freq_cap(2,j) + 1;
        
        diameters{2,i} = [diameters{2,i}, Vessel{v}.Radius * 2];
        lengths{2,i} = [lengths{2,i}, Vessel{v}.Length];
    elseif Vessel{v}.Type == "Venule"
        diameters{3,Vessel{v}.Venous_Generation} = [diameters{3,Vessel{v}.Venous_Generation}, Vessel{v}.Radius * 2];
        lengths{3,Vessel{v}.Venous_Generation} = [lengths{3,Vessel{v}.Venous_Generation}, Vessel{v}.Length];
    end
end

avg_diameters = zeros(3,num_generations);
avg_lengths = zeros(3,num_generations);
for g = 1:num_generations
    for k = 1:3
        avg_diameters(k,g) = nanmean(diameters{k,g});
        avg_lengths(k,g) = nanmean(lengths{k,g});
    end
end

dlmwrite(filenames{1}, [avg_diameters], '-append')
dlmwrite(filenames{2}, [avg_lengths], '-append')
dlmwrite(filenames{3}, [freq_cap], '-append')

end