close all
clearvars

%motor = hill_muscle_motor(muscle_length, Fmax_motor, vmax_motor,r_activation, a_L, b_L, s);
motor = hill_muscle_motor(10, 100, 100, Inf);
for i=1:1000
    plot(i,motor.Force(i,[i,i]),'k')
end
