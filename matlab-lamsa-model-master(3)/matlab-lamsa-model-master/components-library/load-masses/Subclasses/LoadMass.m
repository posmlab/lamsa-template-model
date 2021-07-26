%% LoadMass class definition
% arguments:
%     mass - mass of the load
%     mass of lever arm - mass of lever arm
%     EMA - effective mechnical advantage 
% min # arguments = 1

classdef LoadMass < Mass
    
    methods (Static)
        function parameters = parameters()
            parameters = ["mass" "mass of lever arm" "EMA";
                "0.01" "0" "1";
                "0" "0" "0.001";
                "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        % constructor
        function obj = LoadMass(m_end,varargin)
            varargin_param_names = {'m_rod','EMA'};
            varargin_default_values = {0, 1};
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
            
            % model
            
            % takes y as [y0 y]

            mass = @(y) m_end/(EMA^2) + m_rod*( (1+1/EMA)^2 + 3*(1/EMA-1)^2 ) /12;
            EMA = @(y) EMA;
            
            % call parent constructor
            obj = obj@Mass(mass, EMA);
        end 
    end
end
