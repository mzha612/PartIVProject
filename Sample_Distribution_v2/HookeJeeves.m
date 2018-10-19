function [x,val] = HookeJeeves(f,univ,n,x1,epsilon,isPlot)
%{
Hooke and Jeeves method

Input:
    @f = function handle for function f,
        input: n-dimensional vector
        output: scalar f(x)
    @univ = function handle for the univariate optimisation function
        input: n-dimensional vector x
        output: n-dimensional gradient rf(x)
    n = number of variables.
    x1 = n-dimensional starting point, x1
    epsilon = used as stopping condition
    isPlot = boolean, whether to plot iterations

Output:
    x = n-dimensional vector for the identified solution
    val = scalar objective function value.

Author: Michael Zhang
Date Created: 4-08-18 
%}
%% Initialisation
% check if isPlot is passed in
if nargin <= 5 
    isPlot = false;
end

% Inital values
y = x1;
x = [x1, x1*inf];
directions = eye(n);
k = 1; % Iteration counter
m = 1; % Step counter

% Convergence test
isConverged = @(x ,k) norm(x(:,end) - x(:,end-1)) < epsilon;

%% Hooke and Jeeves
while (~isConverged(x))
    for dim = 1:n+1
        if dim ~= n+1 % pattern search
            d = directions(:,dim);
        else % exploratory search
            k = k + 1;
            x(:,k) = y(:,m);
            
            % print output
            disp(['Iteration: ',num2str(k-1)])
            disp("Current Position:")
            disp(x(:,k))
            disp("Objective fx:")
            disp(f(x(:,k)))
            
            d = x(:,k) - x(:,k-1);
        end
        % linesearch
        alpha = univ(f,y(:,m),d,0.0001);
%         alpha = 1/k;
        % update
        y(:,m+1) = y(:,m) + alpha.*d;
        m = m + 1;
    end
    
    if k >= 10000
        disp('did not converge within 10000 iterations')
        break
    end
end

%% Plot
% to run if isPlot is passed in, and for dim 2.
if (isPlot && n == 2)
    [xminmax] = minmax(x(1,:));
    [yminmax] = minmax(x(2,:));
    
    fig = figure;
    hold on
    
    % create contour map
    X = linspace(floor(xminmax(1)),ceil(xminmax(2)), 20);
    Y = linspace(floor(yminmax(1)),ceil(yminmax(2)), 20);
%     X = linspace(-2,3, 20);
%     Y = linspace(-4,0, 20);
    [X,Y] = meshgrid(X,Y);
    for i = 1:size(X,2)
        for j = 1:size(Y,2)
            Z(i,j) = f([X(i,j),Y(i,j)]');
        end
    end
    colormap('spring')
    colorbar
    contourf(X,Y,Z,20)
    xlabel('x_1')
    ylabel('x_2')
    
    % save figure as image
    filename = @(k) ['HookeJeeves',num2str(k)];
    hgexport(fig,filename(0)) 
    for c = 1:m-1
        title(['Iteration: ', num2str(ceil(c/3))]);
        % draws lines for each step
        plot([y(1,c),y(1,c+1)],[y(2,c),y(2,c+1)],'-x','color',[80,220,100]/255,'linewidth',1.5)
        hgexport(fig,filename(c)) 
    end
end

%% Output
% assign output values
x = x(:,end);
val = f(x);
end

