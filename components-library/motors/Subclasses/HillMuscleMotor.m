%% HillMuscleMotor class definition
% arguments in required order:
%     muscle_length - rest length of the muscle
%     F_motor_max - maximum amount of force the spring can exert
%     v_motor_max - maximum velocity at which the motor can travel
%     r_activation - rate of activation of the muscle
%     L_initial - initial stretch of the muscle (optional)
%     a_L, b_L, s - phenomenological parameters fitted to describe the
%     shape of the muscle�s length�tension curve (optional)
% min # arguments = 4

classdef HillMuscleMotor < Motor
    
    methods (Static)
        % first row contains parameter names
        % second row contains default values for the loading motor
        % third row contains default values for the unlatching motor
        function parameters = parameters()
            parameters = ["muscle length" "Fmax" "Vmax" "rate of activation"...
                         "damping" "initial length" "a_L" "b_L" "s";
                "0.01" "4" "5" "200" "0.5" "0.01" "2.08" "-2.89" "-0.75";
                "4" "10" "10" "2" "0" "4" "2.08" "-2.89" "-0.75";
                "0" "0" "0" "0" "0" "0" "-Inf" "-Inf" "-Inf";
                "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        function obj = HillMuscleMotor(muscle_length,F_motor_max,v_motor_max,r_activation,motor_damping,varargin)
            % optional parameters
            varargin_param_names = {'L_initial','a_L','b_L','s'};
            varargin_default_values = {muscle_length,2.08,-2.89,-0.75};

            % check and assign optional parameters
            if (nargin < 4)
                error('Hill muscle motor requires at least 4 arguments.');
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
            if (v_motor_max == 0)
                v_motor_max = Inf;
                warning("v_max argument must be nonzero, setting to Inf")
            end
            % model
            % motor.Force = F_motor_max * F_length * F_velocity * F_activation such that:
            % F_length=@(t,x) (x(1)<=(L_initial-(0.7*muscle_length)))*(x(1)>=(L_initial-(1.3*muscle_length)))* exp(-((abs(((((L_initial-x(1))/muscle_length)^b_L)-1)/s))^a_L));
            % F_velocity=@(t,x)(1-(x(2)/v_motor_max))/(1+(x(2)/(v_motor_max/4))); 
            % F_activation=@(t,x)min(r_activation*t,1);
            % ^ obtained from Rosario et al.
            Force = @(t,x) max((x(1)<=(L_initial-(0.7*muscle_length))) * (x(1)>=(L_initial-(1.3*muscle_length))) * F_motor_max * exp(-((abs(((((L_initial-x(1))/muscle_length)^b_L)-1)/s))^a_L)) * (((1-(abs(x(2))/v_motor_max))/(1+(abs(x(2))/(v_motor_max/4))))*(x(2) >= 0) + (x(2) < 0)) .* (min(r_activation*t,1)),0);
            max_force = F_motor_max;
            range=L_initial-muscle_length+(.3*muscle_length);
        %     motor.range=muscle_length;
            velocity=v_motor_max;
            damping= @(t,x) motor_damping*x(2);
            
            % call parent constructor
            obj = obj@Motor(max_force, range, velocity, Force, muscle_length, damping);
        end  
    end
end
%% Citations
% Rosario MV, Sutton GP, Patek SN, Sawicki GS. 2016 Muscle�spring dynamics in time-limited, elastic movements.
%   Proc. R. Soc. B 283: 20161561. http://dx.doi.org/10.1098/rspb.2016.1561

