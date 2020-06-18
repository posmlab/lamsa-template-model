%% linear_spring struct
% arguments in required order:
%     k - spring constant
%     m_s - mass of the spring
%     F_spring_max - maximum amount of force the spring can exert
%     (optional)
% min # arguments = 2

function spring = linear_spring(k,m_s,varargin)
    % optional parameters
    varargin_param_names = {'F_spring_max'};
    varargin_default_values = {Inf};
    
    % check and assign optional parameters
    if (nargin < 2)
        error('Linear spring requires at least 2 arguments.');
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
    spring.Force = @(t,x)-k*x(1).*(abs(k*x(1))<F_spring_max);  
    spring.mass = m_s;
end

