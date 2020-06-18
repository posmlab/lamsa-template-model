close all
clearvars

%motor = hill_muscle_motor(muscle_length, Fmax_motor, vmax_motor,r_activation, a_L, b_L, s);
motor = hill_muscle_motor(10, 100, 100, Inf);
%linmotor = linear_motor(Fmax_motor, vmax_motor, range_of_motion);
linmotor = linear_motor(100, 100, 3);
x=linspace(0,20,1000);
for i=1:length(x)
    f(i)=motor.Force(5,[x(i),0]);
    g(i)=linmotor.Force(5,[x(i),0]);
end
plot(x,f,'ko')
hold on
plot(x,g,'go')
