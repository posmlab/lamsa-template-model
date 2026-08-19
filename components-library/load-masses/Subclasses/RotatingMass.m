%% RotatingMass class definition
%
% this type of mass allows large angle lever arm rotations for variable EMA
%
% arguments:
%     m_end - mass of the load
%     m_rod - mass of lever arm
%     L1    - segment of the lever arm closest to the spring
%     L2    - segment of the lever arm furthest from the spring
%     x_rest- rest length of the spring
% min # arguments = 1

classdef RotatingMass < Mass
    
    properties
        L1, theta_0
    end
    
    methods (Static)
        function parameters = parameters()
            parameters = ["mass" "mass of lever arm" "L1" "L2" "theta initial";
                "0.01" "0" "0.001" "0.001" "0";
                "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "1.57"];
        end
    end
    
    methods
        
        % constructor
        function obj = RotatingMass(m_end,varargin)
            varargin_param_names = {'m_rod','L1','L2','theta_0'};
            varargin_default_values = {0, .001 , .001, .001, 0};
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
            EMA = L1/L2;
            mass = (m_end/(EMA^2) + m_rod*( (1+1/EMA)^2 + 3*(1/EMA-1)^2 ) /12);
            EMA = @(y) EMA;
            
            % call parent constructor
            obj = obj@Mass(mass, EMA);
            obj.L1 = L1;
            obj.theta_0 = theta_0;
        end 
    end
end
