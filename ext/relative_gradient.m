function grad = relative_gradient(fun,x0,grad_pct)
    if nargin==2 %Set default value for grad_pct if not specified
        grad_pct=0.01;
    end
    
    dim=length(x0); %fun has one argument for each dimension of the parspace
    grad=zeros(dim,1);
    for i=1:dim %iterate over each dimension of the parameter space
        %approximate the relative derivative as x0*(f(x+dx)-f(x-dx))/(2dx)
        %make a vector that is identical to x0, but moved slightly along
        %the i-th dimension, to get x-dx
        near=x0;
        near(i)=(1-grad_pct)*x0(i);
        %Similar for x+dx
        far=x0;
        far(i)=(1+grad_pct)*x0(i);
        % take relative partial derivative w.r.t. parameter i
        grad(i)=(fun(far)-fun(near))/(2*grad_pct);
    end
end

