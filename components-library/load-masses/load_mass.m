%% load_mass struct
% arguments:
%     m_end - mass of the load
% min # arguments = 1

function load = load_mass(m_end,varargin)
    varargin_param_names = {'m_rod','EMA'};
    varargin_default_values = {0,1};
    % check and assign optional parameters
    if (nargin < 1)
        error('load mass requires at least 1 argument');
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
    load.mass= m_end/(EMA^2) + m_rod*( (1+1/EMA)^2 + 3*(1/EMA-1)^2 ) /12;
end 