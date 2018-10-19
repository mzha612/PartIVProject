function [Length] = Sample_Length(Vessel_Type,data)
Length = 0;
if Vessel_Type == "Arteriole"
    while(Length < 50)
    x = rand(1);
    b = -0.0021;
    a = -b/973/2;
    Length = (-b -sqrt(b^2-4*a*x))/(2*a);
    end
elseif Vessel_Type == "Venule"
    A = 1.8171;
    B = .00544^-1;
    Length = gamrnd(A,B);
elseif Vessel_Type == "Capillary"
    A = 2.49829;
    B = .00653^-1;
    Length = gamrnd(A,B);
else
    disp("error in sample Diamter")
end

end

