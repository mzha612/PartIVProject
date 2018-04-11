%% Test network
% Create a small network.

clear all 
close all
clc
%% Parameters
global E
R_0 = 50;
E = 0.05; %EGL thickness
r_0 = 1;
x = 2; %number of times to branch
vessels = VesselNumbers(x);
new_radius = R_0;
h = 1 - E/new_radius;
    
for i = 1:x
    [new_radius(i+1), h(i+1)] = Bifuricate(new_radius(i));
end

for i = 1:x

    new_radius(x+i + 1) = new_radius(x-i+1);
    h(x+i +1) = h(x-i + 1);
end


[totalresistance, parallelresistance, branchresistance] = BranchResistance(vessels,h)


figure
hold on
plot(1:2*x+1,parallelresistance)

% plot(1:2*x+1,new_radius/R_0)




%%
function [new_radius, h] = Bifuricate(current_radius)
    global E
    new_radius = current_radius/(sqrt(2));
    
    h = 1 - E/new_radius;
end

function [vessels] = VesselNumbers(Branches)
   vessels = 2*ones(1,Branches*2+1);
   for i = 0:Branches
       vessels(i+1)= vessels(i+1)^i;
   end
   for i = 0:Branches-1
       vessels(Branches*2+1 - i) = vessels(i + 1);
   end
end

function [totalresistance, parallelresistance, branchresistance] = BranchResistance(vessels,h)
    branchresistance = 1./sqrt(h) ;%% TODO Julies function
    parallelresistance = vessels .* branchresistance.^-1;
       
    parallelresistance = (parallelresistance).^-1;
    totalresistance = sum(parallelresistance);
end

