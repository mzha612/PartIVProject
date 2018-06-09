function Resistivity = PoiseuilleFlow(Radius)
%% PoiseuelleFlow v4
%{
Calculates the poiseulle flow resistance for a straight walled tube%}
%}
%{
Inputs:
   Radius 
Outputs:
   	Resistivity     Non-dimensional resisitivity (resistance per length)
%}

%{
Author = Michael Zhang
Date created = 10-06-18
%}

%% Parameters

% TODO: Be non dimensionalised as well
mu = 1;

Resistivity = 8*mu/pi/Radius^4;

end