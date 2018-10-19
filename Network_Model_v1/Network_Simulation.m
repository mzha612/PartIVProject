%% Network_Simulation - Network_Model_v1
%{
Runs many iterations of the network calculation for each EGL thickness
%}

%{
Author = Michael Zhang
Date created = 17-07-18
%}
clear
close all
clc
% 97
% 87
%% Parameters
global epsilon num_vessels num_nodes
global num_generations num_capillaries
global d_P_0

num_capillaries = 180;
max_generations = num_capillaries * 2 + 5;

num_iterations = 25;
num_network = 4;
plot_range = 40;
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
art_diameters = cell(1,num_iterations);
art_lengths = cell(1,num_iterations);

ven_diameters = cell(1,num_iterations);
ven_lengths = cell(1,num_iterations);

cap_diameters = cell(1,num_iterations);
cap_lengths = cell(1,num_iterations);
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
    
%     [Vessel_Ref_thick,~] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,2);
%     R_Ref_thick(j) = (d_P_0 - d_P_inf)/Vessel_Ref_thick{1}.Flow;
%     
    obtainGenInfo(Vessel,filenames)
    clear Vessel_Ref
    toc
    
    for v = 1:num_vessels
    if Vessel{v}.Type == "Arteriole"
        art_diameters{j} = [art_diameters{j}, Vessel{v}.Radius * 2 * 10^6];
        art_lengths{j} = [art_lengths{j}, Vessel{v}.Length * 10^6];
    
    elseif Vessel{v}.Type == "Venule"
        ven_diameters{j} = [ven_diameters{j}, Vessel{v}.Radius * 2 * 10^6];
        ven_lengths{j} = [ven_lengths{j}, Vessel{v}.Length * 10^6];
    
    elseif Vessel{v}.Type == "Capillary"
        cap_diameters{j} = [cap_diameters{j}, Vessel{v}.Radius * 2 * 10^6];
        cap_lengths{j} = [cap_lengths{j}, Vessel{v}.Length * 10^6];
    
    end
end
%     for v = 1:num_vessels
% %         if Vessel{v}.Type
%         diameters{j} = [diameters{j}, Vessel{v}.Radius * 2 * 10^6];
%         lengths{j} = [lengths{j}, Vessel{v}.Length * 10^6];
%     end
%     

end
for v = 1:num_vessels
    disp(Vessel{v});
end
%% Plot

R = R ./ R_Ref;

err = std(R,[],2)/sqrt(num_iterations);

R = mean(R,2);

% figure
subplot(3,2,1)
hold on
errorbar(e_range*10^6,R,err)
xlabel('EGL thickness (\mum)')
ylabel('Resistance (%)')
% nnn = ([{'RTB-RTB'},{'RSB-RSB'},{'RSB-RTB'},{'RTB-RSB'},{'SSB-SSB'}]);
% legend(nnn(num_network),'Location','southeast')

% legend("Base Model", "Rat Mesentery Model" ,'Location','southeast')


%% Write CSV
mmm = [{'Res_RTBRTB,'},{'Res_RSBRSB,'},{'Res_RSBRTB,'},{'Res_RTBRSB,'},{'Res_SSBSSB'}];
header = ['EGL_Thickness,',mmm(num_network),',ste'];
filename = ['network_resistance_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s',string(header));
fprintf(fid,'\n');
fclose(fid);
dlmwrite(filename, [e_range',R,err], '-append')

%% Read Intermediate data for generation info
Network_Metrics = readtable('Digitised_Data.csv');
Network_Metrics = table2array(Network_Metrics);

X_Ven_Gen	=Network_Metrics(:,1);
Y_Num_Cap_Ven	=Network_Metrics(:,2);
X_Art_Gen	=Network_Metrics(:,3);
Y_Num_Cap_Art	=Network_Metrics(:,4);
Length_Freq	=Network_Metrics(:,5);
Art_Len_Freq	=Network_Metrics(:,6);
Ven_Len_Freq	=Network_Metrics(:,7);
Cap_Len_Freq	=Network_Metrics(:,8);
Diameter_Freq =Network_Metrics(:,9);
Art_Dia_Freq=Network_Metrics(:,10);
Ven_Dia_Freq	=Network_Metrics(:,11);
Cap_Dia_Freq	=Network_Metrics(:,12);
Gen_Art	=Network_Metrics(:,13);
Art_Gen_Dia	=Network_Metrics(:,14)*10^-6;
Art_Gen_Len	=Network_Metrics(:,15)*10^-6;
Gen_Ven	=Network_Metrics(:,16);
Ven_Gen_Dia	=Network_Metrics(:,17)*10^-6;
Ven_Gen_Len	=Network_Metrics(:,18)*10^-6;
Gen_Cap	=Network_Metrics(:,19);
Cap_Gen_Dia	=Network_Metrics(:,20)*10^-6;
Cap_Gen_Len	=Network_Metrics(:,21)*10^-6;

dia_data = csvread(filenames{1});
len_data = csvread(filenames{2});
freq_cap = csvread(filenames{3});

%% Distribution of lengths and radii
d_edges = 1.25:2.5:51.25;
l_edges = 50:100:1050;

% figure
% hold on
for i = 1:num_iterations
    h = histogram(art_diameters{i},d_edges,'DisplayStyle','stairs','Normalization','probability','Visible','off');
    art_dia(i,:) = h.Values;
    h = histogram(ven_diameters{i},d_edges,'DisplayStyle','stairs','Normalization','probability','Visible','off');
    ven_dia(i,:) = h.Values;
    h = histogram(cap_diameters{i},d_edges,'DisplayStyle','stairs','Normalization','probability','Visible','off');
    cap_dia(i,:) = h.Values;
    h = histogram(art_lengths{i},l_edges,'DisplayStyle','stairs','Normalization','probability','Visible','off');
    art_len(i,:) = h.Values;
    h = histogram(ven_lengths{i},l_edges,'DisplayStyle','stairs','Normalization','probability','Visible','off');
    ven_len(i,:) = h.Values;
    h = histogram(cap_lengths{i},l_edges,'DisplayStyle','stairs','Normalization','probability','Visible','off');
    cap_len(i,:) = h.Values;
    
    
end

art_dist_dia = mean(art_dia);
ven_dist_dia = mean(ven_dia);
cap_dist_dia = mean(cap_dia);
art_dist_len = mean(art_len);
ven_dist_len = mean(ven_len);
cap_dist_len = mean(cap_len);
% 
art_dist_dia_err = std(art_dia)/sqrt(num_iterations);
ven_dist_dia_err = std(ven_dia)/sqrt(num_iterations);
cap_dist_dia_err = std(cap_dia)/sqrt(num_iterations);
art_dist_len_err = std(art_len)/sqrt(num_iterations);
ven_dist_len_err = std(ven_len)/sqrt(num_iterations);
cap_dist_len_err = std(cap_len)/sqrt(num_iterations);
% art_dist_len_err



subplot(3,2,2)
hold on
plot(Diameter_Freq,Art_Dia_Freq,'r')
plot(Diameter_Freq,Ven_Dia_Freq,'b')
plot(Diameter_Freq,Cap_Dia_Freq,'m')
errorbar(Diameter_Freq(1:20),art_dist_dia,art_dist_dia_err,'r')
errorbar(Diameter_Freq(1:20),ven_dist_dia,ven_dist_dia_err,'b')
errorbar(Diameter_Freq(1:20),cap_dist_dia,cap_dist_dia_err,'m')


subplot(3,2,3)
hold on
plot(Length_Freq,Art_Len_Freq,'r')
plot(Length_Freq,Ven_Len_Freq,'b')
plot(Length_Freq,Cap_Len_Freq,'m')
errorbar(Length_Freq(1:10),art_dist_len,art_dist_len_err,'r')
errorbar(Length_Freq(1:10),ven_dist_len,ven_dist_len_err,'b')
errorbar(Length_Freq(1:10),cap_dist_len,cap_dist_len_err,'m')



%% Frequency of capillaries per generation
freq_art_cap_err = std(freq_cap(1:2:end,:),'omitnan')./sqrt(sum(~isnan(freq_cap(1:2:end,:)) & freq_cap(1:2:end,:)~= 0));
freq_ven_cap_err = std(freq_cap(2:2:end,:),'omitnan')./sqrt(sum(~isnan(freq_cap(2:2:end,:)) & freq_cap(2:2:end,:)~= 0));

freq_art_cap =  nanmean(freq_cap(1:2:end,:),1);
freq_ven_cap =  nanmean(freq_cap(2:2:end,:),1);

% Write
header = ['Art,','ven,','std,'];
filename = ['network_gen_avg_freq_cap',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [freq_art_cap;freq_ven_cap;freq_art_cap_err;freq_ven_cap_err]', '-append')

% Plot

subplot(3,2,4)
hold on
errorbar((1:plot_range),freq_art_cap(1:plot_range),freq_art_cap_err(1:plot_range),'b')
errorbar(1:size(freq_cap,2),freq_ven_cap,freq_ven_cap_err,'r')
plot(X_Art_Gen,Y_Num_Cap_Art,'b-o',X_Ven_Gen,Y_Num_Cap_Ven,'r-o')
xlabel('Generation Number')
ylabel('Number of Capillaries')
legend('Arteriole','Venule','Arteriole','Venule')

%% Average length per generation
len_art_err = std( len_data(1:3:end,:),'omitnan')./sqrt(sum(~isnan(len_data(1:3:end,:)) & len_data(1:3:end,:)~= 0))* 10^6;
len_cap_err = std( len_data(2:3:end,:),'omitnan')./sqrt(sum(~isnan(len_data(2:3:end,:)) & len_data(2:3:end,:)~= 0))* 10^6;
len_ven_err = std( len_data(3:3:end,:),'omitnan')./sqrt(sum(~isnan(len_data(3:3:end,:)) & len_data(3:3:end,:)~= 0))* 10^6;
% len_ven_err(plot_range:end) = NaN;

len_art =  nanmean(len_data(1:3:end,:)) * 10^6;
len_cap =  nanmean(len_data(2:3:end,:))* 10^6;
len_ven =  nanmean(len_data(3:3:end,:))* 10^6;
% len_art(plot_range:end) = NaN;
% len_cap(plot_range:end) = NaN;
% len_ven(plot_range:end) = NaN;

% Write
header = ['Art,','cap,','ven,','std,'];
filename = ['network_gen_avg_len_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [len_art;len_cap;len_ven;len_art_err;len_cap_err;len_ven_err]', '-append')

% Plot

subplot(3,2,5)
hold on
errorbar((1:plot_range),len_art(1:plot_range),len_art_err(1:plot_range),'b')
errorbar((1:plot_range),len_cap(1:plot_range),len_cap_err(1:plot_range),'m')
errorbar((1:plot_range),len_ven(1:plot_range),len_ven_err(1:plot_range),'r')
plot(Gen_Art,Art_Gen_Len* 10^6,'b-o' ,Gen_Cap,Cap_Gen_Len* 10^6,'m-o',Gen_Ven,Ven_Gen_Len* 10^6,'r-o')
xlabel('Generation Number')
ylabel('Mean Length (\mum)')
legend('Arteriole','Capillary','Venule','Arteriole','Capillary','Venule')

%% Average diameter per generation
dia_art_err = std(dia_data(1:3:end,:),'omitnan')./sqrt(sum(~isnan(dia_data(1:3:end,:)) & dia_data(1:3:end,:)~= 0))* 10^6;
dia_cap_err = std(dia_data(2:3:end,:),'omitnan')./sqrt(sum(~isnan(dia_data(2:3:end,:)) & dia_data(2:3:end,:)~= 0))* 10^6;
dia_ven_err = std(dia_data(3:3:end,:),'omitnan')./sqrt(sum(~isnan(dia_data(3:3:end,:)) & dia_data(3:3:end,:)~= 0))* 10^6;
% dia_ven_err(plot_range:end) = NaN;

dia_art =  nanmean(dia_data(1:3:end,:))* 10^6;
dia_cap =  nanmean(dia_data(2:3:end,:))* 10^6;
dia_ven =  nanmean(dia_data(3:3:end,:))* 10^6;
% dia_art(plot_range:end) = NaN;
% dia_cap(plot_range:end) = NaN;
% dia_ven(plot_range:end) = NaN;

% Write
header = ['Art,','cap,','ven,','std,'];
filename = ['network_gen_avg_dia_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [dia_art;dia_cap;dia_ven;dia_art_err;dia_cap_err;dia_ven_err]', '-append')

% Plot

subplot(3,2,6)
hold on
errorbar((1:plot_range),dia_art(1:plot_range),dia_art_err(1:plot_range),'b')
errorbar((1:plot_range),dia_cap(1:plot_range),dia_cap_err(1:plot_range),'m')
errorbar((1:plot_range),dia_ven(1:plot_range),dia_ven_err(1:plot_range),'r')
plot(Gen_Art,Art_Gen_Dia* 10^6,'b-o',Gen_Cap,Cap_Gen_Dia* 10^6,'m-o',Gen_Ven,Ven_Gen_Dia* 10^6,'r-o')
xlabel('Generation Number')
ylabel('Mean Diameter (\mum)')
legend('Arteriole','Capillary','Venule','Arteriole','Capillary','Venule')


% Write
header = ['gen,','art_freq,','ven_freq,','art_freq_ste,','ven_freq_ste,',...
    'art_len,','cap_len,','ven_len,','art_len_ste,','cap_len_ste,','ven_len_ste,',...
    'art_dia,','cap_dia,','ven_dia,','art_dia_ste,','cap_dia_ste,','ven_dia_ste,'];
filename = ['network_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename, [(1:size(freq_art_cap,2));freq_art_cap;freq_ven_cap;freq_art_cap_err;freq_ven_cap_err;...
    len_art;len_cap;len_ven;len_art_err;len_cap_err;len_ven_err;...
    dia_art;dia_cap;dia_ven;dia_art_err;dia_cap_err;dia_ven_err;]', '-append')


 header = ['Diameter_Freq,','art_dia_dist,','cap_dia_dist,','ven_dia_dist,','art_dia_dist_ste,','cap_dia_dist_ste,','ven_dia_dist_ste,'];
filename = ['network__dist_dia_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename,[Diameter_Freq(1:20)';art_dist_dia;cap_dist_dia;ven_dist_dia;art_dist_dia_err;cap_dist_dia_err;ven_dist_dia_err;]','-append')


 header = ['Length_Freq,','art_len_dist,','cap_len_dist,','ven_len_dist,','art_len_dist_ste,','cap_len_dist_ste,','ven_len_dist_ste,'];
filename = ['network__dist_len_',num2str(num_network),'.csv'];
fid = fopen(filename,'w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite(filename,[Length_Freq(1:10)';art_dist_len;cap_dist_len;ven_dist_len;art_dist_len_err;cap_dist_len_err;ven_dist_len_err;]','-append')


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