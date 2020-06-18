function motor = linear_motor(Fmax_motor, vmax_motor, range_of_motion)
if (nargin==3)
    motor.Force =@(t,x) (Fmax_motor*(1-x(2)/vmax_motor)) .* (abs(x(1))<=range_of_motion);
    motor.max_force= motor.Force(0,[0,0]);
else
    error('Linear motor requires 3 arguments');
end 