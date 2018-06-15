%% Test sample length distribution
num_sample = 1000;
x = zeros(1,num_sample);

for n = 1:num_sample
%     x(n) = Sample_Length("Arteriole");
    
    x(n) = Sample_Length("Venuole");
    
%     x(n) = Sample_Length("Capillary");
end

figure 
histogram(x,'Normalization','probability','BinWidth',50)

num_sample = 1000;
x = zeros(1,num_sample);

for n = 1:num_sample
%     x(n) = Sample_Diameter("Arteriole");
    
%     x(n) = Sample_Diameter("Venuole");
    
    x(n) = Sample_Diameter("Capillary");
end


figure
histogram(x,'Normalization','probability','BinWidth',1.25)
