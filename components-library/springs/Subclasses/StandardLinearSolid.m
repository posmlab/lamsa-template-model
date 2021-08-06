%% StandardLinearSolid class definition
% Kelvin-Voigt representation of a standard linear solid
% arguments in required order:
%     k - spring constant
%     m_s - mass of the spring
%       (optional)
%     F_spring_max - maximum amount of force the spring can exert
%     (optional)
% min # arguments = 1

classdef StandardLinearSolid < Spring
    
    properties
        F_history, k_0, k_inf, eta
    end
    
    methods(Static)
        % the necessary parameters to make a LinearSpring
        function parameters = parameters()
            parameters = ["k_0" "k_inf" "eta" "mass" "rest length";
                "2000" "2000" "1" "0" "0.01";
                "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        % constructor
        function obj = StandardLinearSolid(k_0, k_inf, eta, varargin)
            % optional parameters
            varargin_param_names = {'mass','rest_length'};
            varargin_default_values = {0,0.01};
            % check and assign optional parameters
            if (nargin < 2)
                error('StandardLinearSolid2 requires at least 3 arguments.');
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
            
            range = 1; %
            
            k1 = k_0;
            if k_0 ~= k_inf
                k2 = k_inf*k_0 / (k_0 - k_inf);
            else
                k2 = 1E20;
            end
            
            Force = @(t,y) ( -k1*k2*(y(1)) - k1*eta*y(2)) / ( k1+k2 );
            
            % call parent constructor
            obj = obj@Spring(mass, range, rest_length, Force);
            
            obj.k1 = k1;
            obj.k2 = k2;
            obj.eta = eta;
        end
        
        function obj = set.F_history(obj, new_history)
            
            obj.F_history = new_history;
            
            F_dot = 0;
            k1 = obj.k1;
            k2 = obj.k2;
            eta = obj.eta;
            hist = obj.F_history;
            
            if size(hist,1) > 1
                F_change = hist(end,2)-hist(end-1,2);
                t_change = hist(end,1)-hist(end-1,1);
                F_dot = F_change/t_change;
            end

            obj.Force = @(t,y) ( -k1*k2*(y(1)) - k1*eta*y(2) + eta*F_dot) / ( k1+k2 );
        end
    end
end