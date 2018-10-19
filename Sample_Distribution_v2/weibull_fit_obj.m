function val = weibull_fit_obj(param)
global data x weights

weights = ones(size(data));
% weights(3:end-3,1) = 1;
val = norm(weights.*(data - wblcdf(x,param(1),param(2))));


% val = 1;
end