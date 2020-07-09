%% exponential_spring struct

% arguments in required order:
%     k_0 - initial spring constant
%     characteristic_length - exponential growth constant of force
%     m_s - mass of the spring
%     (optional)
%     F_spring_max - maximum amount of force the spring can exert
%     (optional)
% min # arguments = 2

function spring = exponential_spring(k_0,characteristic_length, varargin)
    % optional parameters
    varargin_param_names = {'m_s','F_spring_max'};
    varargin_default_values = {0,Inf};
    % check and assign optional parameters
    if (nargin < 2)
        error('Exponential spring requires at least 2 arguments.');
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
    spring.Force = @(t,x)min(characteristic_length*k_0*(exp(-x(1)/characteristic_length)-1),realmax).*(abs(characteristic_length*k_0*(exp(-x(1)/characteristic_length)-1))<F_spring_max);
    spring.mass = m_s;
    spring.range= characteristic_length*log((F_spring_max/(characteristic_length*k_0))-1);
end
