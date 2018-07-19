%% plot_res
%{
Loads up each file, and plots
%}

%{
Author = Michael Zhang
Date created = 13-06-18
%}

clear
close all
clc

%%
m1 = readtable('network_resistance_1.csv');
m1 = table2array(m1);
egl_thick = m1(:,1);
m1 = m1(:,2);

m2 = readtable('network_resistance_2.csv');
m2 = table2array(m2);
m2 = m2(:,2);

m3 = readtable('network_resistance_3.csv');
m3 = table2array(m3);
m3 = m3(:,2);

m4 = readtable('network_resistance_4.csv');
m4 = table2array(m4);
m4 = m4(:,2);

m5 = readtable('network_resistance_5.csv');
m5 = table2array(m5);
m5 = m5(:,2);

m6 = readtable('prelimresults3.csv');
m6 = table2array(m6);
m6x = m6(:,1);
m6 = m6(:,3);

figure
plot(egl_thick,m1,egl_thick,m2,egl_thick,m3,egl_thick,m4,egl_thick,m5,m6x,m6)
xlabel('EGL thickness (um)')
ylabel('Resistance')
legend('RTB-RTB','RSB-RSB','RSB-RTB','RTB-RSB','SSB-SSB','rat','Location','southeast')

header = ['EGL_Thickness,','Res_RTBRTB,','Res_RSBRSB,','Res_RSBRTB,','Res_RTBRSB,','Res_SSBSSB,','Res_Rat,'];
fid = fopen('network_res_all.csv','w');
fprintf(fid,'%s\n',header);
fclose(fid);
dlmwrite('network_res_all.csv', [egl_thick,m1,m2,m3,m4,m5,m6], '-append')