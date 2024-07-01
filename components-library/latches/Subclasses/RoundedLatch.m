%% RoundedLatch class definition
% arguments in required order:
%     R - raduis of the latch
%     m_L - mass of the latch
%     coeff_fric - coefficient of friction between the latch and load mass
%     (optional)
%     v_0 - initial velocity of the latch (optional)
%     min_latching_dist - minimum distance of loading, otherwise loading
%     fails (optional)
%     max_latching_dist - can only load until this point until unlatching
%     begins (optional)
%     runway_length - This is a number that changes the shape of the latch
%     by adding in a 'runway' portion, i.e. a portion of the latch that is
%     completely horizontal/flat before the rounded portion. The default
%     value is 0, which corresponds to a runway with a length 0 (which is
%     just the same shape as the latch) 
% min # arguments = 2

classdef RoundedLatch < Latch
    
    methods (Static)
        function parameters = parameters()
            parameters = ["radius" "mass" "μ" "v_0" "min latching dist" "max latching dist" "runway length";
                "0.0121" "0.00" "0" "1" "0" "Inf" "0";
                "0" "0" "0" "0" "0" "0" "-Inf";
                "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        % constructor
        function obj = RoundedLatch(R, m_L, varargin)
            % optional parameters
            varargin_param_names = {'coeff_fric', 'v_0','min_latching_dist','max_latching_dist','runway_length'};
            varargin_default_values = {0,0,0,Inf,0};

            % check and assign optional parameters
            if (nargin < 2)
                error('Rounded latch requires at least 2 arguments.');
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
            max_width = R + runway_length;
            mass = m_L;

            yL = @(x) (x>=runway_length)*(R*(1-sqrt(1-(x-runway_length)^2/R^2)));
            yL_prime = @(x) (x>=runway_length)*min(abs(((x-runway_length)/(R*sqrt(1-((x-runway_length)/(R))^2)))), realmax);
            yL_doubleprime = @(x) (x >= runway_length)*min(abs((((1-((x-runway_length)/(R))^2)*(R^2)+((x-runway_length)^2))/( ((1-((x-runway_length)/(R))^2)^(3/2))*(R^3)))),realmax);

            y_L = {yL, yL_prime, yL_doubleprime}; % stores yL and its derivatives
            min_latching_dist = abs(min_latching_dist);
            max_latching_dist = abs(max_latching_dist);
            runway_length = runway_length;
            
            % call parent constructor
            obj = obj@Latch(coeff_fric, v_0, max_width, mass, y_L, min_latching_dist, max_latching_dist, runway_length);
        end
    end
end

