function motor = linear_motor(Fmax_motor, vmax_motor, range_of_motion); 
    motor.Force =@(t,x) (Fmax_motor*(1-x(2)/vmax_motor)) .* (abs(x(1))<=range_of_motion);
    motor.Time_independent = true;
end 
