function [Diameter] = Sample_Diameter(Vessel_Type)
%{
    Samples angle from gamma distribution with shape / scale form R.
%}

%{
    Michael Zhang
    21-7-18
%}
if Vessel_Type == "Arteriole"
    A = 7.102;
    B = 0.508^-1;
    Diameter = gamrnd(A,B);
elseif Vessel_Type == "Venule"
    A = 3.698;
    B = .184^-1;
    Diameter = gamrnd(A,B);
elseif Vessel_Type == "Capillary"
    A = 11.19;
    B = 1.19^-1;
    Diameter = gamrnd(A,B);
else
    disp("error in sample Diamter")
end

end

