function [Length] = Sample_Length(Vessel_Type,data)

Range = data.Length_Freq;
% load('Frequency_Length.mat');
% clear Length_Freq
if Vessel_Type == "Arteriole"
    
%     load('Frequency_Length.mat', 'Art_Len_Freq');
    Frequency = data.Art_Len_Freq;
%     clear Art_Len_Freq
elseif Vessel_Type == "Venule"
%     load('Frequency_Length.mat', 'Ven_Len_Freq');
    Frequency = data.Ven_Len_Freq;
%     clear Ven_Len_Freq
elseif Vessel_Type == "Capillary"
%     load('Frequency_Length.mat', 'Cap_Len_Freq');
    Frequency = data.Cap_Len_Freq;
%     clear Cap_Len_Freq
else
    disp("error in sample length")
end

bin = cumsum(Frequency);
p = rand(1);
for i = 1:size(Range,1)
    if p < bin(i)
        Length = Range(i);
        break
    end
end

Length  = Length + 50*(rand(1)-0.5);
Length = Length*10^-6;
end

