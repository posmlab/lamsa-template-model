function motor = hill_muscle_motor(muscle_length, Fmax_motor, vmax_motor,r_activation, a_L, b_L, s,L_initial)
if (nargin==4) %(muscle_length, Fmax_motor, vmax_motor,r_activation) we assume values of a_L,b_L,s from values in Rosario paper
    a_L=2.08;
    b_L=-2.89;
    s=-.75;
    L_initial=muscle_length;
else
    error('Hill muscle motor requires 4 or 8 inputs (including a_L, b_L, s, and L_initial, or using default values');
end
motor.F_length=@(t,x) (x(1)<=(L_initial-(0.7*muscle_length)))*(x(1)>=(L_initial-(1.3*muscle_length)))* exp(-((abs(((((L_initial-x(1))/muscle_length)^b_L)-1)/s))^a_L));
motor.F_velocity=@(t,x)((1-(abs(x(2))/vmax_motor))/(1+(abs(x(2))/(vmax_motor/4))));
motor.F_activation=@(t,x)min(r_activation*t,1);
motor.Force = @(t,x) (x(1)<=(L_initial-(0.7*muscle_length)))*(x(1)>=(L_initial-(1.3*muscle_length)))*Fmax_motor * exp(-((abs(((((L_initial-x(1))/muscle_length)^b_L)-1)/s))^a_L)) * ((1-(abs(x(2))/vmax_motor))/(1+(abs(x(2))/(vmax_motor/4)))) .* (min(r_activation*t,1));
motor.max_force = Fmax_motor;
end 

