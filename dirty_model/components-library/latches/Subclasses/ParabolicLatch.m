%% parabolic_latch object
% arguments in required order:
%     coeff_parabola - coefficient of x^2 in the parabolic shape definition
%     parabola_width - distance past the end of the runway the parabolic
%     latch extends
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
% min # arguments = 3

classdef ParabolicLatch < Latch
    
    methods (Static)
        function parameters = parameters()
            parameters = ["coeff_parabola" "parabola_width" "mass" "mu" "v_0" "min_latching_dist" "max_latching_dist" "runway_length";
                "0.005" "0.005" "0.003" "0" "0" "0" "Inf" "0"];
        end
    end
    
    methods
        % constructor
        function obj = ParabolicLatch(coeff_parabola, parabola_width, m_L, varargin)
            % optional parameters
            varargin_param_names = {'coeff_fric', 'v_0','min_latching_dist','max_latching_dist','runway_length'};
            varargin_default_values = {0,0,0,Inf,0};

            % check and assign optional parameters
            if (nargin < 3)
                error('Parabolic latch requires at least 3 arguments.');
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
            max_width = parabola_width + runway_length;
            mass = m_L;
%maybe put this back (x>=(parabola_width+runway_length))*realmax
            yL = @(x) (x>=runway_length)*(coeff_parabola*(x-runway_length)^2);%+(x>=(parabola_width+runway_length))*100000;
            yL_prime = @(x) (x>=runway_length)*2*coeff_parabola*(x-runway_length); %+(x>=(parabola_width+runway_length))*100000;
            yL_doubleprime = @(x) (x >= runway_length)*(2*coeff_parabola); %+(x>=(parabola_width+runway_length))*100000;

            y_L = {yL, yL_prime, yL_doubleprime}; % stores yL and its derivatives
            min_latching_dist = abs(min_latching_dist);
            max_latching_dist = abs(max_latching_dist);
            runway_length = runway_length;
            
            obj = obj@Latch(coeff_fric, v_0, max_width, mass, y_L, min_latching_dist, max_latching_dist, runway_length);
        end
    end
end

