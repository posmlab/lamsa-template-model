close all
clearvars
tic
debug = false;

load_time_constraint=Inf;

%parameters for the loading motor
Fmax_motor = 20;
range_of_motion = 3;
vmax_motor=100.0000;
%extra parameters for hill muscle motor
muscle_length=.1;
r_activation=1E-2;
%struct initialization
loading_motor = linear_motor(Fmax_motor, vmax_motor, range_of_motion);
unlatching_motor=hill_muscle_motor(muscle_length, Fmax_motor, vmax_motor,r_activation);

%parameters for the load and struct initialization
m=10;
load = load_mass(m);

%parameters for the latch and struct initialization
R=2;
m_L= 10;
coeff_fric = 0;
v_0L=0;
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%parameters for the spring and struct initialization
k=3;
m_s=1;
F_spring_max=1E4;
% characteristic_length for exponential spring
characteristic_length = 5;
spring = linear_spring(k, m_s, F_spring_max);
% spring=exponential_spring(k, characteristic_length, m_s,F_spring_max);

%solving the model to get output 
[sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring);

%plot force outputs of motor as a function of time


