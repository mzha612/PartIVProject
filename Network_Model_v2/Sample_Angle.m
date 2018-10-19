function angle = Sample_Angle()%% Sample Angle
%{
    Samples angle from normal distribution with mean and std form R.
%}

%{
    Michael Zhang
    21-7-18
%}

if rand(1) < 0.025
    angle = 180;
else
    mu = 114.6;
    sigma = 43.6;  
    angle = normrnd(mu,sigma);
end
end