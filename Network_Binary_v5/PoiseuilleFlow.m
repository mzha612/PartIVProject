function Resistivity = PoiseuilleFlow(Radius)
%% PoiseuelleFlow v5
%{
Calculates the poiseuille flow resistance for a straight walled tube
Derived from Q = Gh3/12mu - wikipedia flow between two infinite parallel
plates.
h = diameter,
G = dp/dx
%}

%{
Inputs:
    Radius, m
Outputs:
   	Resistivity  dimensional resisitivity (resistance per length) % Pa.s.m^-3
%}

%{
Author = Michael Zhang
Date created = 10-06-18
%}

mu = 1e-3;
Resistivity = 3/2*mu/Radius^3; 

end