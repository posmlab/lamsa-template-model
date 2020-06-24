%% linear_motor struct
% arguments in required order:
%     F_motor_max - maximum amount of force the spring can exert
%     v_motor_max - maximum velocity at which the motor can travel
%     range_of_motion - how far the motor can contract
% min # arguments = 3

function motor = linear_motor(F_motor_max, v_motor_max, range_of_motion)
    if (nargin == 3)
        motor.Force = @(t,x) (F_motor_max*(1-x(2)/v_motor_max)) .* (abs(x(1))<=range_of_motion);
    else
        error('Linear motor requires 3 arguments.');
    end
    motor.max_force = F_motor_max;
    motor.range=range_of_motion;
end 
