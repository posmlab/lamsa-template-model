%% ExponentialElasticExtensionalSpring class definition

% arguments in required order:
%     E - modulus
%     A - cross sectional area
%     L - length
%     rho - density (optional)
%     sigma_f - failure strength in Pascals (optional)
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
        % constructor
        function obj = ExponentialElasticExtensionalSpring(E, A, L, varargin)
            varargin_param_names = {'rho','sigma_f'};
            varargin_default_values = {0,Inf};

            % check and assign optional parameters
            if (nargin < 3)
                error('Exponential elastic extensional spring requires at least 3 arguments.');
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
            k = (E*A)/L;
            m = rho*A*L;
            Fmax = sigma_f*A;
            
            % call parent constructor
            obj = obj@ExponentialSpring(k, L, m, Fmax);
        end
    end
end
