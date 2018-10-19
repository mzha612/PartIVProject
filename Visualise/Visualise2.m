%% Visualise network

close all
clear all
clc

global epsilon num_vessels num_nodes
global num_generations num_capillaries
global d_P_0
global cs

% load('Network Structure');
load('Network Base Final Test 5 25');

j = 23;
% 6 7 23
num_capillaries = 180;
size_scale = 200;
Node = Network_Node{j};
Vessel = Network_Vessel{j};

e_step = 0.1;
e_range = 0.01:e_step:1.01;
e_range = e_range *10^-6; % m

d_P_0 = 13000;                    % Inital Pressure, Pa
d_P_inf = 1000;

n_P_0 = 1;                      % Inital Pressure, non-dimensionalised
n_P_inf = d_P_inf/d_P_0;        % Final Pressure, non-dimensionalised
n_R_0 = 1;

num_nodes = numel(Node);
num_vessels = numel(Vessel);

epsilon = 1*10^-6;    % EGL thickness, m
cs = 1;

[Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,0);
[Vessel_Ref,Node_Ref] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,1);

for v = 1:num_vessels
    radius(v) = Vessel{v}.Radius;
    res(v) = Vessel{v}.Resistance;
    res_ref(v) = Vessel_Ref{v}.Resistance;
    flow(v) = Vessel{v}.Flow;
    flow_ref(v) = Vessel_Ref{v}.Flow;
end

for n = 1:num_nodes
    pressure(n) = Node{n}.n_Pressure;
    pressure_ref(n) = Node_Ref{n}.n_Pressure;
end
    

radius = radius/max(radius);
linewidth = radius*7.5;

epsilon = 1*10^-6;    % EGL thickness, m
cs = 0.1;
[Vessel,Node] = Solve_Network(Vessel,Node,n_P_0,n_P_inf,0);

for v = 1:num_vessels
    res1(v) = Vessel{v}.Resistance;
    flow1(v) = Vessel{v}.Flow;
end
for n = 1:num_nodes
    pressure1(n) = Node{n}.n_Pressure;
%     pressure_ref(n) = Node_Ref{n}.n_Pressure;
end
res = (res1 - res)./res_ref; %
flow = (flow1 - flow)./flow_ref; %
pressure = (pressure1 - pressure)./pressure_ref;

max_res = max(abs(res));
max_flow = max(abs(flow));
max_pressure = max(abs(pressure));

res = res/max_res;
flow = flow/max_flow;
pressure = pressure/max_pressure;

%% Define Color
color_space = logspace(-1,0,10);
fileID = fopen('definecolor.txt','w');
for i = 1:10
    def_green = ['\definecolor{g',num2str(i-1),'}{rgb}{0,',num2str(color_space(i)),',0}'];
    def_red = ['\definecolor{r',num2str(i-1),'}{rgb}{',num2str(color_space(i)),',0,0}'];
    def_blue = ['\definecolor{b',num2str(i-1),'}{rgb}{0,0,',num2str(color_space(i)),'}'];
    fprintf(fileID,'%s\n',def_red);
    fprintf(fileID,'%s\n',def_green);
    fprintf(fileID,'%s\n',def_blue);
end
fclose(fileID);

%% Colormap
fileID = fopen('colormap.txt','w');
fprintf(fileID,'%s\n','\pgfplotsset{');
fprintf(fileID,'\t%s\n',['colormap={redgreengrad}{']);
for j = 1:10
    fprintf(fileID,'\t%s\n',['rgb =(',num2str(color_space(11-j),3) ,', 0, 0),']);
end
fprintf(fileID,'\t%s\n',['rgb =(0, 0, 0),']);
for j = 1:10
    fprintf(fileID,'\t%s\n',['rgb =(0,',num2str(color_space(j),3) ,', 0),']);
end
fprintf(fileID,'\t%s\n',['},']);

fprintf(fileID,'\t%s\n',['colormap={redbluegrad}{']);
for j = 1:10
    fprintf(fileID,'\t%s\n',['rgb =(',num2str(color_space(11-j),3) ,', 0, 0),']);
end
fprintf(fileID,'\t%s\n',['rgb =(0, 0, 0),']);
for j = 1:10
    fprintf(fileID,'\t%s\n',['rgb =(0, 0,',num2str(color_space(j),3),'),']);
end

fprintf(fileID,'\t%s\n',['},']);
fprintf(fileID,'%s\n',['}']);
fclose(fileID);

%% Drawline

for i = ["Arteriole","Venule"]
    % Resistance
    figure
    hold on
    
    %     for n = 1:num_nodes
    %         if Node{n}.Type == i
    %             %           Node{n}.Type = "Venule";
    %             plot(Node{n}.xy(1),Node{n}.xy(2),'bo')
    %         end
    %     end
    title('individual resistance')
    
    
    %          Z=1:1:100;
    %Colormap is defined as a 3 column matrix, each row being an RGB triplet
    %     map = zeros(numel(res),3);
    %     map(:,1)=0;
    %     map(:,2)=sort(res);
    %     map(:,3)=0;%./max(Z);
    %     %Set the current Colormap
    %     colormap(map);
    %     %Display Colorbar
    %     colorbar
    %     colorbar();
    % res = rand(1,360);
    
    fileID = fopen(strcat(i,"_res.txt"),'w');
    rounded_res = abs(round(res*9));
    for v = 1:num_vessels
        if Vessel{v}.Type == i
            if res(v) < 0
                plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}.xy_End(2)],'Color',[-res(v) 0 0 ] ,'LineWidth',linewidth(v))
                draw_vessel_res = ['\draw[r',num2str(rounded_res(v)),',line width = ',num2str(linewidth(v),2),']', ...
                    ' (', num2str(Vessel{v}.xy_Start(1)/size_scale,3), ',', num2str(Vessel{v}.xy_Start(2)/size_scale,3),') -- (', ...
                    num2str(Vessel{v}.xy_End(1)/size_scale,3), ',', num2str(Vessel{v}.xy_End(2)/size_scale,3),');'];
                
            else
                draw_vessel_res = ['\draw[g',num2str(rounded_res(v)),',line width = ',num2str(linewidth(v),2),']', ...
                    '(', num2str(Vessel{v}.xy_Start(1)/size_scale,3), ',', num2str(Vessel{v}.xy_Start(2)/size_scale,3),') -- (', ...
                    num2str(Vessel{v}.xy_End(1)/size_scale,3), ',', num2str(Vessel{v}.xy_End(2)/size_scale,3),');'];
                
                plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}.xy_End(2)],'Color',[0 res(v) 0] ,'LineWidth',linewidth(v))
            end
            fprintf(fileID,'%s\n',draw_vessel_res);
        end
    end
    fclose(fileID);
    
    
    % Flow
    figure
    hold on
    title('flow')
    fileID = fopen(strcat(i,"_flow.txt"),'w');
%     rounded_flow = abs(round(flow*9));
    
    for v = 1:num_vessels
        if abs(flow(v)) > 0.05
            rounded_flow(v) = 9;
        else
            rounded_flow(v) = 0;
        end 
    end
    
    for v = 1:num_vessels
        if Vessel{v}.Type == i
            if flow(v) < 0  
                %         sqrt(-flow(v))
                plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}.xy_End(2)],'Color',[-flow(v) 0 0  ] ,'LineWidth',linewidth(v))
                draw_vessel_flow = ['\draw[r',num2str(rounded_flow(v)),',line width = ',num2str(linewidth(v),2),']', ...
                    ' (', num2str(Vessel{v}.xy_Start(1)/size_scale,3), ',', num2str(Vessel{v}.xy_Start(2)/size_scale,3),') -- (', ...
                    num2str(Vessel{v}.xy_End(1)/size_scale,3), ',', num2str(Vessel{v}.xy_End(2)/size_scale,3),');'];
                
            else
                draw_vessel_flow = ['\draw[g',num2str(rounded_flow(v)),',line width = ',num2str(linewidth(v),2),']', ...
                    '(', num2str(Vessel{v}.xy_Start(1)/size_scale,3), ',', num2str(Vessel{v}.xy_Start(2)/size_scale,3),') -- (', ...
                    num2str(Vessel{v}.xy_End(1)/size_scale,3), ',', num2str(Vessel{v}.xy_End(2)/size_scale,3),');'];
                
                plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}.xy_End(2)],'Color',[0  flow(v) 0] ,'LineWidth',linewidth(v))
            end
            fprintf(fileID,'%s\n',draw_vessel_flow);
        end
    end
    fclose(fileID);
    
    
    
    % Pressure
    % Flow
    figure
    hold on
    title('pressure')
    fileID = fopen(strcat(i,"_pressure.txt"),'w');
%     rounded_pressure = abs(round(pressure*9));
    
    for n = 1:num_nodes
        if abs(pressure(n)) > 0.2
            rounded_pressure(n) = 9;
        else
            rounded_pressure(n) = 0;
        end 
    end
    
    for n = 1:num_nodes
        if Node{n}.Type == i
            if pressure(n) < 0
            scatter(Node{n}.xy(1),Node{n}.xy(2),[],[-pressure(n) 0 0 ],'filled')
            draw_pressure = ['\draw[r',num2str(rounded_pressure(n)),', fill = r',num2str(rounded_pressure(n)),'] (', num2str(Node{n}.xy(1)/size_scale,3),',',num2str(Node{n}.xy(2)/size_scale,3),') circle [radius=0.2];'];
            
            else
            scatter(Node{n}.xy(1),Node{n}.xy(2),[],[0 0 pressure(n)],'filled')
            
            draw_pressure = ['\draw[b',num2str(rounded_pressure(n)),', fill = b',num2str(rounded_pressure(n)),'] (', num2str(Node{n}.xy(1)/size_scale,3),',',num2str(Node{n}.xy(2)/size_scale,3),') circle [radius=0.2];'];
            
            end   
            fprintf(fileID,'%s\n',draw_pressure);
%           el
%             plot(Node{n}.xy(1),Node{n}.xy(2),'bo')
%             draw_vessel = ['\draw[line width = ',num2str(linewidth(v),2),']', ...
%                     ' (', num2str(Vessel{v}.xy_Start(1)/size_scale,3), ',', num2str(Vessel{v}.xy_Start(2)/size_scale,3),') -- (', ...
%                     num2str(Vessel{v}.xy_End(1)/size_scale,3), ',', num2str(Vessel{v}.xy_End(2)/size_scale,3),');'];
%             fprintf(fileID,'%s\n',draw_vessel);
        end
        
    end
    for v = 1:num_vessels
        if Vessel{v}.Type == i
                plot([Vessel{v}.xy_Start(1),Vessel{v}.xy_End(1)],[Vessel{v}.xy_Start(2),Vessel{v}. xy_End(2)],'Color',[0  0 0] ,'LineWidth',linewidth(v))
            draw_vessel_pressure = ['\draw[black, line width = ',num2str(linewidth(v),2),']', ...
                    ' (', num2str(Vessel{v}.xy_Start(1)/size_scale,3), ',', num2str(Vessel{v}.xy_Start(2)/size_scale,3),') -- (', ...
                    num2str(Vessel{v}.xy_End(1)/size_scale,3), ',', num2str(Vessel{v}.xy_End(2)/size_scale,3),');'];
            fprintf(fileID,'%s\n',draw_vessel_pressure);
        end
    end
    fclose(fileID);
end


