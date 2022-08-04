%% DeactivatingMotor class definition
%
% Instead of pulling the latch out of the way, this deactivating motor
% pushes latch in place during loading and then deactivates to allow mass
% past during the unlatching phase.
%
% arguments in required order:
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
    
    properties
        r_deactivation
    end
    
    methods (Static)
        % first row contains parameter names
        % second row contains default values for the loading motor
        % third row contains default values for the unlatching motor
        function parameters = parameters()
            parameters = ["Vmax" "range of motion" "deactivation rate" "damping" "voltage fraction" "time delay";
                "10" "0.005" "1E7" "0.5" "1" "0";
                "1" "0.005" "1E7" "0" "1" "0";
                "0" "0" "0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        function obj = DeactivatingMotor(v_motor_max, range_of_motion,r_deactivation, motor_damping, varargin)
            % optional parameters
            varargin_param_names = {'voltage_fraction','time_delay','no_braking'};
            varargin_default_values = {1, 0, true};

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
            
            % model
            max_force = 0;
            range = range_of_motion;
            velocity = voltage_fraction*v_motor_max;
            Force = @(t,x) 0; % this is assigned in solve_lamsa
            rest_length = 0;
            damping = @(t,x) motor_damping*x(2);
            activation_delay = time_delay;

            
            % call parent constructor
            obj = obj@Motor(max_force, range, velocity, Force, rest_length,damping,activation_delay);
            
            obj.r_deactivation = r_deactivation;
        end 
    end
end
