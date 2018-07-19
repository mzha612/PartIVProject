function [Diameter] = Sample_Diameter(Vessel_Type,Data)
% load('Frequency_Diameter.mat', 'Diameter_Freq');
Range = Data.Diameter_Freq;
% load('Frequency_Length.mat');
% clear Length_Freq
if Vessel_Type == "Arteriole"
%     load('Frequency_Diameter.mat', 'Art_Dia_Freq');
    Frequency = Data.Art_Dia_Freq;
%     clear Art_Len_Freq
elseif Vessel_Type == "Venule"
%     load('Frequency_Diameter.mat', 'Ven_Dia_Freq');
    Frequency = Data.Ven_Dia_Freq;
%     clear Ven_Len_Freq
elseif Vessel_Type == "Capillary"
%     load('Frequency_Diameter.mat', 'Cap_Dia_Freq');
    Frequency = Data.Cap_Dia_Freq;
%     clear Cap_Len_Freq
else
    disp("error in sample Diamter")
end

bin = cumsum(Frequency);
p = rand(1);
for i = 1:size(Range,1)
    if p < bin(i)
        Diameter = Range(i);
        break
    end
end

Diameter  = Diameter + 2.5*(rand(1)-0.5);
Diameter = Diameter*10^-6;



% if Vessel_Type == "Arteriole"
%     Ds = load('diameters_gamma.mat','x');
% %     Ds = load('diameters_lognorm.mat','x');
% %     Ds = load('diameters_norm.mat','x');
%     Diameter = Ds.x(randi([1,100000],[1,1]))*10^-6;
% end
end

