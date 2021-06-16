%% load_mass object
% arguments:
%     mass - mass of the load
%     spring_length - rest length of the spring
%     mass of lever arm - mass of lever arm
%     EMA - effective mechnical advantage 
% min # arguments = 1

classdef NewMass < Mass
    
    methods (Static)
        function parameters = parameters()
            parameters = ["mass" "mass of lever arm" "spring_length" "L1" "L2";
                "0.01" "0" ".01" ".05" ".05";
                "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        % constructor
        function obj = NewMass(m_end,varargin)
            varargin_param_names = {'m_rod','l','L1','L2'};
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
%             l = spring_length;
%             y0 = 
%             y1 = l + 
%             y2 = 
%             d = 
%             sinalpha = sqrt(1-((y1^2-y2^2)/(2*y2*d))^2);
            mass= @(y) m_end/(EMA^2) + m_rod*( (1+1/EMA)^2 + 3*(1/EMA-1)^2 ) /12;
            EMA = @(y) L1/L2;
            obj = obj@Mass(mass, EMA);
        end 
    end
end
