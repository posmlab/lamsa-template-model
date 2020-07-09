%% linear_spring struct
% arguments in required order:
%     k - spring constant
%     m_s - mass of the spring
%       (optional)
%     F_spring_max - maximum amount of force the spring can exert
%     (optional)
% min # arguments = 1

function spring = linear_spring(k,varargin)
    % optional parameters
    varargin_param_names = {'m_s','F_spring_max'};
    varargin_default_values = {0,Inf};
    
    % check and assign optional parameters
    if (nargin < 1)
        error('Linear spring requires at least 1 argument.');
    end
    if (length(varargin)>length(varargin_param_names))
        error('Too many input parameters');
    end
    for i=1:length(varargin)
        eval([varargin_param_names{i} '=varargin{i};'])
    end
    for i=(length(varargin)+1):length(varargin_param_names)
        eval([varargin_param_names{i} '=varargin_default_values{i};'])
    end
    % model 
    if (F_spring_max == Inf)
        F_spring_max = realmax;
    end
    
    
    % only works for compressing the spring 
    spring.k = k;
    spring.Force = @(t,x) min(realmax,-k*x(1)) .* (abs(k*x(1))<F_spring_max); 
    spring.mass = m_s;
end