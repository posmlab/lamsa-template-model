%% DeactivatingMotor object
% arguments in required order:
%     F_motor_max      - maximum amount of force the spring can exert AT MAX
%                        VOLTAGE, i.e. voltage_fraction = 1 the F_motor_max 
%                        will get scaled down w/ the voltage fraction
%     v_motor_max      - maximum velocity at which the motor can travel AT MAX
%                        VOLTAGE, i.e. voltage_fraction = 1 the F_motor_max 
%                        will get scaled down w/ the voltage fraction 
%     range_of_motion  - how far the motor can contract
%     r_deactivation   - rate of deactivation of the muscle keeping the
%                        deactivating motor in place, as it deactivates
%     voltage_fraction - the fraction of the max voltage at which the motor
%                        is being powered - the motor's maxmax force and
%                        maxmax velocity get scaled by this number,
%                        resulting in a new Force vs Velocity 'curve' (its
%                        just linear) (optional)
%     no_braking       - this is an optional boolean flag. The default
%                        value is true. When set to true, the motor will mostly obey a
%                        linear force-velocity tradeoff, but in the special
%                        case where the motor's velocity is greater than
%                        vmax, the motor force will just be zero.
%                        If set to false, the linear motor will give a 
%                        negative force value when the
%                        motor's velocity is greater than vmax.
%                        
% min # arguments = 4

classdef DeactivatingMotor < Motor
    
    methods (Static)
        % first row contains parameter names
        % second row contains default values for the loading motor
        % third row contains default values for the unlatching motor
        function parameters = parameters()
            parameters = ["Fmax" "Vmax" "range of motion" "deactivation rate" "voltage fraction";
                "10" "10" "0.005" "1E7" "1" ;
                "0.25" "1" "0.005" "1E7" "1";
                "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        function obj = DeactivatingMotor(F_motor_max, v_motor_max, range_of_motion,r_deactivation, varargin)
            % optional parameters
            varargin_param_names = {'voltage_fraction','no_braking'};
            varargin_default_values = {1, true};

            % check and assign optional parameters
            if (nargin < 4)
                error('Deactivating motor requires at least 4 arguments.');
            end
            if (length(varargin)>length(varargin_param_names))
                error('Too many input parameters. Deactivating motor requires at least 4 paramters.');
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

            max_force = -voltage_fraction*F_motor_max;
            range = range_of_motion;
            velocity = voltage_fraction*v_motor_max;

%             if (no_braking)
%                 Force = @(t,x)-max((max_force*(1-x(2)/velocity)) .* (abs(x(1))<=range_of_motion), 0);
%             else
%                 Force = @(t,x)-(max_force*(1-x(2)/velocity)) .* (abs(x(1))<=range_of_motion);
%             end
% what they did with acitivation rate in HillMuscleMotor
%    Force = @(t,x) (x(1)<=(L_initial-(0.7*muscle_length))) * (x(1)>=(L_initial-(1.3*muscle_length))) 
%               * F_motor_max * exp(-((abs(((((L_initial-x(1))/muscle_length)^b_L)-1)/s))^a_L)) * 
%               ((1-(abs(x(2))/v_motor_max))/(1+(abs(x(2))/(v_motor_max/4)))) .* (min(r_activation*t,1));
%old force from deactivating
%            Force = @(t,x)max_force*(t==0);
            Force = @(t,x) max_force * exp(-r_deactivation*t);
            obj = obj@Motor(max_force, range, velocity, Force);
            %obj = obj@LinearMotor(-max_force, v_motor_max, range_of_motion, voltage_fraction, false);
        end 
    end
end
