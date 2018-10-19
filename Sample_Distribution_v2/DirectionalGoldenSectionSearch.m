function [ t ] = DirectionalGoldenSectionSearch(f, xi, di, eps, a, b)
% Golden Section Search from point xi in direction of vector di,
% for functions that have two or more variables:
% We minimise f(x + t*d)
%   Input
%   f = function name string of multivariable function (required input)
%   xi = point from which we are searching (required input)
%   di = search direction for multivariable function f (required input)
%   eps = length of final interval (required input)
%   a = lower interval limit for t
%   b = upper interval limit for t

% The algorithm computes
%   an = lower limit of interval of uncertainty for t
%   bn = upper limit of interval of uncertainty for t

%   Output
%   t = min (an,bn) which is at most eps from the true minimum.


large = 1; % for initial interval of uncertainty

if nargin < 6, b = large; end
if nargin < 5, a = -large; end

alpha = (sqrt(5)-1)/2;

an = a;
bn = b;

lam = an + (1-alpha)*(bn-an);
mu = an + alpha*(bn-an);

flam = f(xi + lam*di);
fmu = f(xi + mu*di);

k = 1;

while (bn-an > eps)
    
    if (flam > fmu)
        % bn unchanged
        an = lam;
        lam = mu;
        flam =fmu;
        mu = an + alpha*(bn-an);
        fmu = f(xi + mu*di);
    else
        % an unchanged
        bn = mu;
        mu = lam;
        fmu = flam;
        lam = an + (1-alpha)*(bn-an);
        flam = f(xi + lam*di);
    end
    k = k+1;
end

% final interval of uncertainty is [an,bn]

% choose the minimum of the two:
if (f(xi + an*di) < f(xi + bn*di))
    t = an;
else
    t = bn;
end

fprintf('Golden Section Search found f(xi + %f * di) = %f to be minimal, bn-an = %f\n',t,f(xi + t*di),bn-an)

% Check if this one of the initial limits of the search interval
if (t-eps < a) || (t+eps > b)
    fprintf('NOTE: The minimum was attained at the lower or upper bound for t:\n t = %f, lower = %f, upper = %f\n',t,a,b)
end


end

