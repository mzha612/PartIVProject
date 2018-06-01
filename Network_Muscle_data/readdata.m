%% Part IV read txt for microcirculation model
clear 
close all
clc

num_v = 2035;

fid = fopen('networkgeometry.dat','r');
datacell = textscan(fid, '%f%f%f', 'HeaderLines', 13, 'Collect', 8140);
fclose(fid);
data = datacell{1};

for i = 0:num_v-1
    xyz1(i + 1,:) = data(4*i + 1,:);
    xyz2(i + 1,:) = data(4*i + 2,:);
    radius(i + 1) = data(4*i + 3,1);
    in_out(i + 1) = data(4*i + 3,2);
    num_parents(i + 1) = data(4*i + 4,1);
    parent_1(i + 1) = data(4*i + 4,2);
    parent_2(i + 1) = data(4*i + 4,3);
end

figure
hold on



for i = 1:num_v
    if in_out(i) == 0
        c = 'r';
    else
        c = 'b';
    end
    plot3([xyz1(i,1); xyz2(i,1)],[xyz1(i,2); xyz2(i,2)],[xyz1(i,3); xyz2(i,3)],...
        'color',c,'LineWidth',radius(i)*10000)
end
