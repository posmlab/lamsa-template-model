close all
clearvars
tic
debug = false;

load_time_constraint=Inf;

%parameters for the loading motor
Fmax_motor = 10;
range_of_motion = 3;
vmax_motor=100.0000;
%extra parameters for hill muscle motor
muscle_length=.1;
r_activation=.5;
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

%plot force outputs of motor as a function of time HILL MUSCLE
sz=size(sol);
for t = 1:sz(1)
    force_array(t)=unlatching_motor.max_force*unlatching_motor.F_velocity(sol(t,1),[sol(t,2), sol(t,3)])*unlatching_motor.F_activation(sol(t,1),[sol(t,2), sol(t,3)]);
    force_activation(t)=unlatching_motor.F_activation(sol(t,1),[sol(t,2), sol(t,3)]);
    force_velocity(t)=unlatching_motor.F_velocity(sol(t,1),[sol(t,2), sol(t,3)]);
    force_length(t)=unlatching_motor.F_length(sol(t,1),[sol(t,2), sol(t,3)]);
    max_force=unlatching_motor.max_force;
end
figure
hold on
plot(sol(:,1),force_array,'r');
plot(sol(:,1),force_activation,'k');
plot(sol(:,1),force_velocity,'b');
plot(sol(:,1),force_length,'m');
plot(sol(:,1),max_force,'g');
hold off

%plot force outputs of spring as a function of time
unlatching_motor2=linear_motor(Fmax_motor, vmax_motor, range_of_motion);
spring1=linear_spring(k, m_s, F_spring_max);
[sol,transition_times]=solve_model(loading_motor,unlatching_motor2,load,latch,spring1);
sol1=sol;
sz1=size(sol1);
spring2=exponential_spring(k, characteristic_length, m_s,F_spring_max);
[sol,transition_times]=solve_model(loading_motor,unlatching_motor2,load,latch,spring2);
sol2=sol;
sz2=size(sol2);
for i = 1:sz1(1)
    linearspring(i)=spring1.Force(sol1(i,1),[sol1(i,2), sol1(i,3)]);
end
for i = 1:sz2(1)
    exponentialspring(i)=spring2.Force(sol2(i,1),[sol2(i,2), sol2(i,3)]);
end
figure
hold on
plot(sol1(:,1),linearspring,'r');
plot(sol2(:,1),exponentialspring,'k');
hold off