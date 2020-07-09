%% linear_motor struct
% arguments in required order:
%     F_motor_max_max  - maximum amount of force the spring can exert AT MAX
%                        VOLTAGE, i.e. voltage_fraction = 1 the F_motor_max 
%                        will get scaled down w/ the voltage fraction
%     v_motor_max_max  - maximum velocity at which the motor can travel AT MAX
%                        VOLTAGE, i.e. voltage_fraction = 1 the F_motor_max 
%                        will get scaled down w/ the voltage fraction 
%     range_of_motion  - how far the motor can contract
%     voltage_fraction - the fraction of the max voltage at which the motor
%                        is being powered - the motor's maxmax force and
%                        maxmax velocity get scaled by this number,
%                        resulting in a new Force vs Velocity 'curve' (its
%                        just linear) (optional)
% min # arguments = 3

function motor = linear_motor(F_motor_max, v_motor_max, range_of_motion,varargin)
    % optional parameters
    varargin_param_names = {'voltage_fraction'};
    varargin_default_values = {1};

    % check and assign optional parameters
    if (nargin < 3)
        error('Linear motor requires at least 3 arguments.');
    end
    if (length(varargin)>length(varargin_param_names))
        error('Too many input parameters. Linear motor requires at least 3 paramters.');
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

    motor.max_force = voltage_fraction*F_motor_max;
    motor.range = range_of_motion;
    motor.velocity = voltage_fraction*v_motor_max;

    motor.Force = @(t,x)max((motor.max_force*(1-x(2)/motor.velocity)) .* (abs(x(1))<=range_of_motion), 0);
    % if the force expression comes out negative, sets it to 0
end 
