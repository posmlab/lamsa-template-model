range_of_motion=10;
F_max=10;
V_max=10;

% extra parameters for hill muscle motor
muscle_length = 10;
r_activation = Inf;

% loading motor struct initialization
 motor = linear_motor(F_max, V_max, range_of_motion);
% motor = hill_muscle_motor(muscle_length, F_max, V_max, r_activation);

m=5;
load=load_mass(m);

sol=solve_direct_actuation(motor,load);
