%% LinearElasticExtensionalSpring class definition
% arguments in required order:
%     E - modulus
%     A - cross sectional area
%     L - length
%     rho - density (optional)
%     sigma_f - failure strength in Pascals (optional)
% min # arguments = 1

classdef LinearElasticExtensionalSpring < LinearSpring
    
    methods(Static)
        % the necessary parameters to make a LinearElasticExtensionalSpring
        function parameters = parameters()
            parameters = ["E" "A" "L" "rho" "sigma_f";
                "1" "9.14E-4" "1" "1" "0.318";
                "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        % constructor
        function obj = LinearElasticExtensionalSpring(E, A, L, varargin)
            varargin_param_names = {'rho','sigma_f'};
            varargin_default_values = {0,Inf};

            % check and assign optional parameters
            if (nargin < 1)
                error('Linear elastic extensional spring requires at least 1 argument.');
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
            F_spring_max = sigma_f*A;
            
            % call parent constructor
            obj = obj@LinearSpring(k, m, F_spring_max, L);
        end
    end
end
