%% ExponentialElasticExtensionalSpring class definition

% arguments in required order:
%     k_0 - initial spring constant
%     characteristic_length - exponential growth constant of force
%     m_s - mass of the spring
%     (optional)
%     F_spring_max - maximum amount of force the spring can exert
%     (optional)
% min # arguments = 3

classdef ExponentialElasticExtensionalSpring < ExponentialSpring
    
    methods(Static)
        
        % the necessary parameters and default values for a LinearSpring
        function parameters = parameters()
            parameters = ["E" "A" "L" "rho" "sigma_f";
                "0.5" "0.5" "0.001" "10" "Inf";
                "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        
        function obj = ExponentialElasticExtensionalSpring(E, A, L, varargin)
            varargin_param_names = {'rho','sigma_f'};
            varargin_default_values = {0,Inf};

            % check and assign optional parameters
            if (nargin < 3)
                error('Exponential elastic extensional spring requires at least 3 argument.');
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

            k = (E*A)/L;
            m = rho*A*L;
            Fmax = sigma_f*A;
            
            obj = obj@ExponentialSpring(k, L, m, Fmax);
        end
    end
end
