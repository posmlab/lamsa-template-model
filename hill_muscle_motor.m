function motor = hill_muscle_motor(muscle_length, Fmax_motor, vmax_motor,r_activation, a_L, b_L, s)
if (nargin==4) %(muscle_length, Fmax_motor, vmax_motor,r_activation) we assume values of a_L,b_L,s from values in Rosario paper
    motor.a_L=2.08;
    motor.b_L=2.89;
    motor.s=-.75;
elseif (nargin == 7) %assume 
    motor.a_L=a_L;
    motor.b_L=b_L;
    motor.s=s;  
else
    error('Hill muscle motor requires 4 or 7 inputs (including a_L, b_L, or s, or using default values');
end
motor.F_length=@(t,x) exp(-((abs((((x(1)/muscle_length)^motor.b_L)-1)/motor.s))^motor.a_L));
motor.F_velocity=@(t,x)(1-(x(2)/vmax_motor))/(1+(x(2)/(vmax_motor/4)));  
motor.F_activation=@(t,x)min(r_activation*t,1);
motor.Force = @(t,x) (((x(1))/muscle_length)<=1)*((x(1)/muscle_length)>=0.7)* Fmax_motor*(1-(x(2)/vmax_motor))/(1+(x(2)/(vmax_motor/4)))*(exp(-((abs(((((muscle_length-x(1))/muscle_length)^(motor.b_L))-1)/motor.s)^motor.a_L))))*min(r_activation*t,1);  
motor.max_force = Fmax_motor;
end 

