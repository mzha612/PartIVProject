function val = gamma_fit_obj(param)
global data x weights

len = length(data);

weights = ones(size(data));
% weights(1:end-len/2,1) = 5;
% weights(end-len/2:end,1) = 5;
val = norm(weights.*(data - gamcdf(x,param(1),param(2))));


% val = 1;
end