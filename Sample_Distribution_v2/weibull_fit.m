function [cdf_fit,p] = weibull_fit(domain,cdf_data,param)
global data x

x = domain;

data = cdf_data; 
f_fit = @weibull_fit_obj;
univ = @DirectionalGoldenSectionSearch;
n = 2;
x1 = param;
epsilon = 0.001;

[p,~] = HookeJeeves(f_fit,univ,n,x1,epsilon);

cdf_fit = wblcdf(x,p(1),p(2));
end