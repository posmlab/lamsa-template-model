%% linear_voltage_motor struct
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
%                        just linear)
% min # arguments = 4


function motor = linear_voltage_motor(F_motor_max_max, v_motor_max_max, range_of_motion, voltage_fraction)

if (nargin < 4)
    error("voltage_motor requires you input a voltage_fraction, representing the fraction of the motor's max voltage")
end

motor.max_force = voltage_fraction*F_motor_max_max;
motor.range = range_of_motion;
motor.velocity = voltage_fraction*v_motor_max_max;

motor.Force = @(t,x)(motor.max_force*(1-x(2)/motor.velocity)) .* (abs(x(1))<=range_of_motion);


end