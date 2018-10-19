% Fit mean per generation
clear all
close all
clc

% Gen_Dia_Data = load('Gen_Dia.mat');
% 
% Art_Gen_Dia = Gen_Dia_Data.Art_Gen_Dia;
% Art_Gen_Dia(isnan(Art_Gen_Dia)) = [];
% 
% Cap_Gen_Dia = Gen_Dia_Data.Cap_Gen_Dia;
% Cap_Gen_Dia (isnan(Cap_Gen_Dia)) = [];
% 
% Ven_Gen_Dia = Gen_Dia_Data.Ven_Gen_Dia;
% Ven_Gen_Dia(isnan(Ven_Gen_Dia)) = [];
% 
% %% Capillary
% Cap_Mean = mean(Cap_Gen_Dia);
% var = var(Cap_Gen_Dia);
% 
% %% Arteriole
% % modelFun =  @(p,x) p(3) .* (x ./ p(1)).^(p(2)-1) .* exp(-(x ./ p(1)).^p(2));
% modelFun =  @(p,x) p(1)*x.^(p(2)) + p(3);
% 
% startingVals = [1 -.3 0];
% coefEsts = nlinfit(1:length(Art_Gen_Dia), Art_Gen_Dia', modelFun, startingVals);
% xgrid = linspace(0,length(Art_Gen_Dia),100);
% figure
% hold on
% plot(1:length(Art_Gen_Dia),Art_Gen_Dia)
% line(xgrid, modelFun(coefEsts, xgrid), 'Color','r');
% 
% %% Venule
% modelFun =  @(p,x) p(1)*x.^(p(2)) + p(3);
% % modelFun =  @(p,x) p(3) .* (x ./ p(1)).^(p(2)-1) .* exp(-(x ./ p(1)).^p(2));
% startingVals = [10 2 5];
% 
% startingVals = [1 -.3 0];
% coefEsts2 = nlinfit(1:length(Ven_Gen_Dia), Ven_Gen_Dia', modelFun, startingVals);
% xgrid = linspace(0,length(Ven_Gen_Dia),100);
% figure
% hold on
% plot(1:length(Ven_Gen_Dia),Ven_Gen_Dia)
% line(xgrid, modelFun(coefEsts2, xgrid), 'Color','r');

%% LENGTH
Gen_Len_Data = load('Gen_Len.mat');

Art_Gen_Len = Gen_Len_Data.Art_Gen_Len;
Art_Gen_Len(isnan(Art_Gen_Len)) = [];

Cap_Gen_Len = Gen_Len_Data.Cap_Gen_Len;
Cap_Gen_Len (isnan(Cap_Gen_Len)) = [];

Ven_Gen_Len = Gen_Len_Data.Ven_Gen_Len;
Ven_Gen_Len(isnan(Ven_Gen_Len)) = [];

%% Capillary
modelFun =  @(p,x) p(1)*x + (p(2));

startingVals = [-.05 500];
coefEstsc = nlinfit(1:length(Cap_Gen_Len), Cap_Gen_Len', modelFun, startingVals);
xgrid = linspace(0,length(Cap_Gen_Len),100);
figure
hold on
plot(1:length(Cap_Gen_Len),Cap_Gen_Len)
line(xgrid, modelFun(coefEstsc, xgrid), 'Color','r');

%% Arteriole

modelFun =  @(p,x) p(1)*x + (p(2));

startingVals = [-.05 500];
coefEstsa = nlinfit(1:length(Art_Gen_Len), Art_Gen_Len', modelFun, startingVals);
xgrid = linspace(0,length(Art_Gen_Len),100);
figure
hold on
plot(1:length(Art_Gen_Len),Art_Gen_Len)
line(xgrid, modelFun(coefEstsa, xgrid), 'Color','r');

%% Venule

modelFun =  @(p,x) p(1)*x + (p(2));

startingVals = [-.05 500];
coefEstsv = nlinfit(1:length(Ven_Gen_Len), Ven_Gen_Len', modelFun, startingVals);
xgrid = linspace(0,length(Ven_Gen_Len),100);
figure
hold on
plot(1:length(Ven_Gen_Len),Ven_Gen_Len)
line(xgrid, modelFun(coefEstsv, xgrid), 'Color','r');

