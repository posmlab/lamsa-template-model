close all
clearvars
tic
debug = false;
addpath(genpath(pwd));

dateString = string(datetime);
cleanDateString = regexprep(dateString, " ", "_");
cleanDateString = regexprep(cleanDateString, ":", "_");
        
load_time_constraint=Inf;


load_time_constraint = Inf;

%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 20;
loading_motor_range_of_motion = 3;
v_max_loading_motor = 10.0000;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 10;
loading_motor_r_activation = Inf;

% loading motor struct initialization
loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor2 = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
m=10;

% load mass struct initialization
load = load_mass(m);

%% latch

% latch parameters
R=2;
m_L= 100;
coeff_fric = 0;
v_0L=0;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% spring paramters
k = 6; % k or k_0 depending on linear or exponential spring
m_s=1;
F_spring_max=1E4;

% extra parameters for exponential spring
% should be a negative value
characteristic_length = -5;

% spring struct initialization
spring = linear_spring(k, m_s, F_spring_max);
%spring = exponential_spring(k, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
unlatching_motor_range_of_motion = 3;
F_max_unlatching_motor=10;
v_max_unlatching_motor=5;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 10;
unlatching_motor_r_activation = Inf;

unlatching_motor = hill_muscle_motor(unlatching_motor_muscle_length, F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_r_activation);
unlatching_motor2= linear_motor(F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_range_of_motion);

%% end editable parameters



%plot force outputs of motor as a function of time HILL MUSCLE
[sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring,cleanDateString);
sz=size(sol);
musclesol=sol;
[sol,transition_times]=solve_model(loading_motor,unlatching_motor2,load,latch,spring,cleanDateString);
linmotsz=size(sol);
for t = 1:sz(1)
    force_array(t)=unlatching_motor.max_force*unlatching_motor.F_length(musclesol(t,1),[musclesol(t,4), musclesol(t,5)])*unlatching_motor.F_velocity(musclesol(t,1),[musclesol(t,4), musclesol(t,5)])*unlatching_motor.F_activation(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    force_activation(t)=unlatching_motor.F_activation(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    force_velocity(t)=unlatching_motor.F_velocity(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    force_length(t)=unlatching_motor.F_length(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    max_force=unlatching_motor.max_force;
    
end
for t=1:linmotsz(1)
    linear_force(t)=unlatching_motor2.Force(sol(t,1),[sol(t,4), sol(t,5)]);
end
figure
hold on
plot(musclesol(:,1),force_array,'r');
plot(musclesol(:,1),force_activation,'k');
plot(musclesol(:,1),force_velocity,'b');
plot(musclesol(:,1),force_length,'m');
plot(musclesol(:,1),max_force,'g');
plot(sol(:,1),linear_force,'y');
hold off



%plot force outputs of spring as a function of time
spring1=linear_spring(k, m_s, F_spring_max);
[sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring1,cleanDateString);
sol1=sol;
sz1=size(sol1);
spring2=exponential_spring(k, characteristic_length, m_s,F_spring_max);
[sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring2,cleanDateString);
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

%comparing linear and exponential spring force vs displacement for kopt,
%k<k_opt
%values must be lower than k_opt otherwise whole range of x values isnt
%reached
k_opt=loading_motor.max_force/loading_motor_range_of_motion;
figure 
plotspot=1;
for k = [k_opt/5,k_opt/2,k_opt, k_opt*2]
    lin_spring=linear_spring(k, m_s, F_spring_max);
    expo_spring=exponential_spring(k, characteristic_length, m_s,F_spring_max);
    [sol,transition_times]=solve_model(loading_motor,unlatching_motor2,load,latch,lin_spring,cleanDateString);
    linsol=sol;
    linsz=size(linsol);
    [sol,transition_times]=solve_model(loading_motor,unlatching_motor2,load,latch,expo_spring,cleanDateString);
    exposol=sol;
    exposz=size(exposol);
    expospring=5+zeros(exposz(1));
    for i = 1:exposz(1)
        expospring(i)=expo_spring.Force(exposol(i,1),[exposol(i,2), exposol(i,3)]);
    end
    for i = 1:linsz(1)
        linspring(i)=lin_spring.Force(linsol(i,1),[linsol(i,2), linsol(i,3)]);
    end
    y=F_max_loading_motor;
    z=0;    
    expo_range=size(exposol(:,2));
    lin_range=size(linsol(:,2));
    subplot(2,2,plotspot); 
    hold on
    plot(linsol(:,2),linspring(1:lin_range(1)),'r');
    plot(exposol(:,2),expospring(1:expo_range(1)),'k');
    line([-loading_motor_range_of_motion,0],[y,y]);
    line([z,z],[0,y]);
      
    hold off
    lin_spring_work=trapz(linsol(:,2),linspring(1:lin_range(1)))
    expo_spring_work=trapz(exposol(:,2),expospring(1:expo_range(1)))
    plotspot=plotspot+1;
end  

%to add
%pdf document in paper repository, results and discussion section
%keep graphs simple to (constant v) to look at b of performance trafeoffs





    
    
    
    
    
    