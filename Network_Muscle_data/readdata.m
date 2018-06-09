% %% Read txt/.dat for microcirculation model
% %{
% Author = Michael Zhang
% Date created = 24-05-18
% %}
% 
% %{
% Opens the .dat and saves the data as a .mat
% %}
% 
% %{
% Inputs:
%     datafile
% BCs: Boundary condition pressures
% Outputs:
%     Vessels: Updated Structure, vessel segment infomation
% %}
% clear
% close all
% clc
% 
% num_v = 2035;
% 
% fid = fopen('networkgeometry.dat','r');
% datacell = textscan(fid, '%f%f%f', 'HeaderLines', 13, 'Collect', 8140);
% fclose(fid);
% 
% data = datacell{1};
% 
% for i = 0:num_v-1
%     xyz1(i + 1,:) = data(4*i + 1,:);
%     xyz2(i + 1,:) = data(4*i + 2,:);
%     radius(i + 1) = data(4*i + 3,1);
%     in_out(i + 1) = data(4*i + 3,2);
%     num_parents(i + 1) = data(4*i + 4,1);
%     parent_1(i + 1) = data(4*i + 4,2);
%     parent_2(i + 1) = data(4*i + 4,3);
% end
% 
% save('nw','data')
