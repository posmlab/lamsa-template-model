close all
clearvars
tic
debug = false;
addpath(genpath(pwd));

dateString = string(datetime);
cleanDateString = regexprep(dateString, " ", "_");
cleanDateString = regexprep(cleanDateString, ":", "_");
        
load_time_constraint=Inf;

metrics={'Pmax'};


%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 20;
loading_motor_range_of_motion = 10;
v_max_loading_motor = 10;

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
R=1;
m_L= 100;
coeff_fric = 0;
v_0L=5;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% spring paramters
k = 1; % k or k_0 depending on linear or exponential spring
m_s=1;
F_spring_max=1E4;
k_opt=loading_motor.max_force/loading_motor_range_of_motion;

% extra parameters for exponential spring
% should be a negative value
characteristic_length = 1;

% spring struct initialization
spring = linear_spring(k_opt, m_s, F_spring_max);
spring2 = exponential_spring(k_opt, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
unlatching_motor_range_of_motion = 10;
F_max_unlatching_motor=10;
v_max_unlatching_motor=10;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 5;
unlatching_motor_r_activation = .4;


unlatching_motor= linear_motor(F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_range_of_motion);
unlatching_motor2 = hill_muscle_motor(unlatching_motor_muscle_length, F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_r_activation);

%% Direct actuating motor parameters
da_range_of_motion=10;
da_F_max=20;
da_V_max=10;

% extra parameters for hill muscle motor
da_muscle_length = 10;
da_r_activation = Inf;

% loading motor struct initialization
da_motor = linear_motor(da_F_max, da_V_max, da_range_of_motion);
da_motor2 = hill_muscle_motor(da_muscle_length, da_F_max, da_V_max, da_r_activation);

%% end struct initialization


%% LaMSA Zone plotting
%zone for specific metric
%     power
%     takeoff velocity
%     kinetic energy
%     science paper figure 
    

%checking force is being calculated correctly
da_sol = solve_direct_actuation(da_motor,load);
sz=size(da_sol);
for i = 1:sz(1)
    force1(i)=da_motor.Force(da_sol(i,1),[da_sol(i,2) da_sol(i,3)]);
    force2(i)=da_sol(i,4);    
end
figure
hold on
plot(da_sol(:,1),force1,'r')
plot(da_sol(:,1),force2,'ko')
hold off

da_sol = solve_direct_actuation(da_motor,load_mass(1E3));
sz=size(da_sol);
figure
hold on
for i = 2:sz(2)
    subplot(1,3,i-1)
    plot(da_sol(:,1),da_sol(:,i))
    title("column"+i)
end
hold off



%iterating over load mass vs max power
mass=logspace(-3,3,350);
%mass=linspace(.01,100);
for m=1:length(mass)
    load_mass_arr(m)=load_mass(mass(m));
end
figure
hold on
for l = 1:length(load_mass_arr)  
    [sol,~]=solve_model(loading_motor,unlatching_motor,load_mass_arr(l),latch,spring,cleanDateString);
    LaMSAsol=sol;
    sol = solve_direct_actuation(da_motor,load_mass_arr(l));
    da_sol=sol;
    %da_max_power(l)=max(mass(l)*da_sol(:,3).*(da_sol(:,4)/mass(l)));
    da_max_power(l)=max(mass(l)*da_sol(:,3).*gradient(da_sol(:,3))./gradient(da_sol(:,1)));
    LaMSA_max_power(l)=max(mass(l)*LaMSAsol(:,3).*gradient(LaMSAsol(:,3))./gradient(LaMSAsol(:,1)));
    da_max(l)=max(da_sol(:,4).*da_sol(:,3));
end
plot(mass,da_max_power,'b')
plot(mass,LaMSA_max_power,'m')
plot(mass,da_max,'g')
set(gca, 'XScale', 'log')
title("max power vs. load mass","Interpreter","latex")
xlabel("load mass [kg]","interpreter","latex")
ylabel("max power [W]","interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
hold off



%size of solution vs load mass
mass=logspace(-3,3,350);
for m=1:length(mass)
    load_mass_arr(m)=load_mass(mass(m));
end
figure
hold on
for l = 1:length(load_mass_arr)  
    sol = solve_direct_actuation(da_motor,load_mass_arr(l));
    sz(l)=length(sol);
end
plot(mass,sz)
set(gca, 'XScale', 'log')
hold off

%iterating over load mass vs kinetic energy
% mass=logspace(-3,3,250);
% %mass=linspace(.01,100);
% for m=1:length(mass)
%     load_mass_arr(m)=load_mass(mass(m));
% end
% figure
% hold on
% for l = 1:length(load_mass_arr)  
%     [sol,~]=solve_model(loading_motor,unlatching_motor,load_mass_arr(l),latch,spring,cleanDateString);
%     LaMSAsol=sol;
%     sol = solve_direct_actuation(da_motor,load_mass_arr(l));
%     da_sol=sol;
%     da_max_ke(l)=max((1/2)*mass(l)*da_sol(:,3).^2);
%     LaMSA_ke(l)=max((1/2)*mass(l)*LaMSAsol(:,3).^2);
% end
% plot(mass,da_max_ke,'b')
% plot(mass,LaMSA_ke,'m')
% set(gca, 'XScale', 'log')
% title("max ke vs. load mass","Interpreter","latex")
% xlabel("load mass [kg]","interpreter","latex")
% ylabel("max ke [J]","interpreter","latex")
% set(gca,'TickLabelInterpreter','latex')
% hold off