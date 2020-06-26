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
    spring.Force = @(t,x)characteristic_length*k_0*(exp(-x(1)/characteristic_length)-1).*(abs(characteristic_length*k_0*(exp(-x(1)/characteristic_length)-1))<F_spring_max);
    spring.mass = m_s;
end
