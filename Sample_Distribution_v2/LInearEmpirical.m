close all
clear all
clc


Diameter_Data = load('Frequency_Diameter.mat');

dia = Diameter_Data.Diameter_Freq;
dia = dia - 2.5;
F_a = Diameter_Data.Art_Dia_Freq;
F_v = Diameter_Data.Ven_Dia_Freq;
F_c = Diameter_Data.Cap_Dia_Freq;

bin_a = cumsum(F_a);
bin_v = cumsum(F_v);
bin_c = cumsum(F_c);

figure
hold on
plot(dia + 2.5,bin_a)
[y, p_gamma_1] = gamma_fit(dia,bin_a,[5;3]);
plot(dia + 2.5,y)
[y, p_weibull_1] = weibull_fit(dia,bin_a,[12;2]);
plot(dia + 2.5,y)
title('dia_a')
legend('empirical','gamma','weibull')

figure
hold on
% dia_edges =0:0.5:50;
plot(dia + 2.5, F_a)
% plot(0:0.1:50,gampdf(0:0.1:50,p_gamma_3(1),p_gamma_3(2)));
plot(dia + 2.5,gampdf(dia,p_gamma_1(1),p_gamma_1(2)));
% plot(dia_edges+1.25,gampdf(dia_edges,p_gamma_2(1),p_gamma_2(2)));
% plot(dia_edges+1.25,gampdf(dia_edges,p_gamma_3(1),p_gamma_3(2)));
plot(dia + 2.5,wblpdf(dia,p_weibull_1(1),p_weibull_1(2)));
% plot(dia_edges+1.25,wblpdf(dia_edges,p_weibull_2(1),p_weibull_2(2)));
% plot(dia_edges+1.25,wblpdf(dia_edges,p_weibull_3(1),p_weibull_3(2)));
title('dia_a pdf')
legend('empirical','gamma','weibull')


figure
hold on
plot(dia + 2.5,bin_v)
[y, p_gamma_2] = gamma_fit(dia,bin_v,[3;7]);
plot(dia + 2.5,y)
[y, p_weibull_2] = weibull_fit(dia,bin_v,[20;2]);
plot(dia + 2.5,y)
title('dia_v')
legend('empirical','gamma','weibull')

figure
hold on
plot(dia + 2.5,bin_c)
[y, p_gamma_3] = gamma_fit(dia,bin_c,[6;1]);
plot(dia + 2.5,y)
[y, p_weibull_3] = weibull_fit(dia,bin_c,[8;3]);
plot(dia + 2.5,y)
title('dia_c')
legend('empirical','gamma','weibull')

% Bins for fitted
figure
hold on
plot(dia + 1.25,F_c)
dia_edges = 0:2.5:50;
for i = 1:length(dia_edges)-1
    r(i) = (gamcdf(dia_edges(i + 1),p_gamma_3(1),p_gamma_3(2)) - gamcdf(dia_edges(i),p_gamma_3(1),p_gamma_3(2)));
end
plot(dia + 1.25,r)
for i = 1:length(dia_edges)-1
    r(i) = (wblcdf(dia_edges(i + 1),p_weibull_3(1),p_weibull_3(2)) - wblcdf(dia_edges(i),p_weibull_3(1),p_weibull_3(2)));
end
plot(dia + 1.25,r)
title('pdf ida_c')
legend('empirical','gamma','weibull')

figure
hold on
dia_edges =0:0.5:50;
g_dia_a = gampdf(dia_edges,p_gamma_1(1),p_gamma_1(2));
g_dia_v = gampdf(dia_edges,p_gamma_2(1),p_gamma_2(2));
g_dia_c = gampdf(dia_edges,p_gamma_3(1),p_gamma_3(2));

% plot(0:0.1:50,gampdf(0:0.1:50,p_gamma_3(1),p_gamma_3(2)));
plot(dia_edges+2.5,gampdf(dia_edges,p_gamma_1(1),p_gamma_1(2)));
plot(dia_edges+2.5,gampdf(dia_edges,p_gamma_2(1),p_gamma_2(2)));
plot(dia_edges+2.5,gampdf(dia_edges,p_gamma_3(1),p_gamma_3(2)));
plot(dia_edges+2.5,wblpdf(dia_edges,p_weibull_1(1),p_weibull_1(2)));
plot(dia_edges+2.5,wblpdf(dia_edges,p_weibull_2(1),p_weibull_2(2)));
plot(dia_edges+2.5,wblpdf(dia_edges,p_weibull_3(1),p_weibull_3(2)));
plot(dia +2.5,F_c/2.5)
plot(dia+2.5,F_v/2.5)
plot(dia+2.5,F_a/2.5)
% plot(dia,gamcdf(dia,p_gamma_3(1),p_gamma_3(2)));
% plot(dia_edges,gamcdf(dia_edges,p_gamma_3(1),p_gamma_3(2)));
%%
% clear all
r = [];
Length_Data = load('Frequency_Length.mat');

len = [0; Length_Data.Length_Freq];
len_edges = -50:100:1050;

F_a = [.2350; Length_Data.Art_Len_Freq];
F_v = [0; Length_Data.Ven_Len_Freq];
F_c = [0; Length_Data.Cap_Len_Freq];

bin_a = cumsum(F_a);
bin_v = cumsum(F_v);
bin_c = cumsum(F_c);

figure
hold on
plot(len,bin_a-0.235)
[y, p_gamma_4] = gamma_fit(len,bin_a,[2;150]);
plot(len,y)
[y, p_weibull_4] = weibull_fit(len,bin_a,[350;1]);
plot(len,y)

P = polyfit(len,F_a,1);

% yfit = polyval(P,len);
% plot(len,yfit)

A = P(2)*1000/2;

plot(len,A^-1 *(P(1)/2*len.^2 + 0.2248*len))
title('len_a')
legend('empirical','gamma','weibull')

% Bins for fitted
figure
hold on
plot(len,F_a)
for i = 1:length(len_edges)-1
    r(i) = (gamcdf(len_edges(i + 1),p_gamma_4(1),p_gamma_4(2)) - gamcdf(len_edges(i),p_gamma_4(1),p_gamma_4(2)));
end
plot(len,r)
% len_edges = 50:100:1050;
for i = 1:length(len_edges)-1
    r(i) = (wblcdf(len_edges(i + 1),p_weibull_4(1),p_weibull_4(2)) - wblcdf(len_edges(i),p_weibull_4(1),p_weibull_4(2)));
end
plot(len,r)
% P = polyfit(len,F_a,1);
% 
% yfit = polyval(P,len);
% plot(len,yfit)
% 
% A = 0.094;
% 
% plot(len,A *(P(1)/2*len.^2 + 0.2248*len))

y = [0:0.1:1];
P = P * A^-1;
a = P(1)/2;
b = P(2);
% cdf = @(y) (-b - sqrt(b^2-4*a.*(-y)))/(2*a);
% plot(cdf(y),y)
cdf = @(y) (-b + sqrt(b^2-4*a.*(-y)))/(2*a);
% cdf = @(y) P(1) * y.^2 + P(2);
plot(cdf(y),y)
% plot(y,cdf(y))
title('pdf len_a')
legend('empirical','gamma','weibull','linear')


figure
hold on
plot(len,bin_v)
[y, p_gamma_5] = gamma_fit(len,bin_v,[2;183]);
plot(len,y)
[y, p_weibull_5] = weibull_fit(len,bin_v,[350;1]);
plot(len,y)
title('len_v')
legend('empirical','gamma','weibull')

% Bins for fitted
figure
hold on
plot(len,F_v)
for i = 1:length(len_edges)-1
    r(i) = (gamcdf(len_edges(i + 1),p_gamma_5(1),p_gamma_5(2)) - gamcdf(len_edges(i),p_gamma_5(1),p_gamma_5(2)));
end
plot(len,r)
% len_edges = 50:100:1050;
for i = 1:length(len_edges)-1
    r(i) = (wblcdf(len_edges(i + 1),p_weibull_5(1),p_weibull_5(2)) - wblcdf(len_edges(i),p_weibull_5(1),p_weibull_5(2)));
end
plot(len,r)
plot(len,100*wblpdf(len,p_weibull_5(1),p_weibull_5(2)))
plot(len,100*gampdf(len,p_gamma_5(1),p_gamma_5(2)))
title('pdf len_v')
legend('empirical','gamma','weibull')

figure
hold on
plot(len,bin_c)
[y, p_gamma_6] = gamma_fit(len,bin_c,[2;153]);
plot(len,y)
[y, p_weibull_6] = weibull_fit(len,bin_c,[400;2]);
plot(len,y)
title('len_c')
legend('empirical','gamma','weibull')

% Bins for fitted
figure
hold on
plot(len(2:end),F_c(2:end))
for i = 1:length(len_edges)-1
    r(i) = (gamcdf(len_edges(i + 1),p_gamma_6(1),p_gamma_6(2)) - gamcdf(len_edges(i),p_gamma_6(1),p_gamma_6(2)));
end
plot(len(2:end),r(2:end))
% len_edges = 50:100:1050;
for i = 1:length(len_edges)-1
    r(i) = (wblcdf(len_edges(i + 1),p_weibull_6(1),p_weibull_6(2)) - wblcdf(len_edges(i),p_weibull_6(1),p_weibull_6(2)));
end
plot(len(2:end),r(2:end))
plot(len+50,100*wblpdf(len,p_weibull_6(1),p_weibull_6(2)))
plot(len+50,100*gampdf(len,p_gamma_6(1),p_gamma_6(2)))

title('pdf len_c')
legend('empirical','gamma','weibull')


figure
hold on
len_edges_lots = 0:25:1050;

% g_len_a = gampdf(len_edges_lots,p_gamma_4(1),p_gamma_4(2));
g_len_v = gampdf(len_edges_lots,p_gamma_5(1),p_gamma_5(2));
g_len_c = gampdf(len_edges_lots,p_gamma_6(1),p_gamma_6(2));

% plot(0:0.1:50,gampdf(0:0.1:50,p_gamma_3(1),p_gamma_3(2)));
% plot(len_edges_lots+50,gampdf(len_edges_lots,p_gamma_4(1),p_gamma_4(2)));
plot(len_edges_lots + 50,100*gampdf(len_edges_lots,p_gamma_5(1),p_gamma_5(2)));
plot(len_edges_lots + 50,100*gampdf(len_edges_lots,p_gamma_6(1),p_gamma_6(2)));
plot(len,F_v)
plot(len,F_c)
axis([100,1000,0,0.25])

% % Write
% % header = ['dia_values,','len_edges_lots','g_dia_a,','g_dia_c,','g_dia_v,','g_dia_c,','g_dia_v,','g_len_c,','g_len_v,'];
% % filename = ['plotdist.csv'];
% % fid = fopen(filename,'w');
% % fprintf(fid,'%s\n',header);
% % fclose(fid);
% % dlmwrite(filename, [dia_edges;len_edges_lots;g_dia_a;g_dia_c;g_dia_v;g_len_c;g_len_v]', '-append')

