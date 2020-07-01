%% initialize stuff
clearvars
addpath(genpath(pwd + "\.."))

%% initialize springs
spring1 = linear_spring(443, 0, 3.73892);
spring2 = linear_spring(2277, 0, 14.75496);
spring3 = linear_spring(4518, 0, 21.9123);

springs = [spring1, spring2, spring3];

%% initialize loading/unlatching motors based on voltages
voltsToMaxForce_motor1=@(V) (V/6)*0.0001765197;
voltsToMaxForce_motor2=@(V) (V/6)*0.0001176798;
voltsToMaxForce_motor3=@(V) (V/6)*0.0000196133;

voltsToMaxVelocity_motor1=@(V) (V/6)*0.0030000000;
voltsToMaxVelocity_motor2=@(V) (V/6)*0.0166666666;
voltsToMaxVelocity_motor3=@(V) (V/6)*0.0166666666;

voltage_motor1 = 6;
voltage_motor2 = 6;
voltage_motor3 = 6;

range_of_motion = Inf;

motor1 = linear_motor(voltsToMaxForce_motor1(voltage_motor1), voltsToMaxVelocity_motor1(voltage_motor1), range_of_motion);
motor2 = linear_motor(voltsToMaxForce_motor2(voltage_motor2), voltsToMaxVelocity_motor2(voltage_motor2), range_of_motion);
motor3 = linear_motor(voltsToMaxForce_motor3(voltage_motor3), voltsToMaxVelocity_motor3(voltage_motor3), range_of_motion);


motors = [motor1, motor2, motor3];
%% load mass
m=0.0038; % guess based on height of jumping robot from slides
load = load_mass(m);

%% latch
R = .1;
%R=0.0032; % guesses based on height of jumping robot from slides
m_L= m/15
coeff_fric = 0;
v_0L=0;
latch = rounded_latch(R,m_L);


spring_k=logspace(2,4,40);
figure
disp('starting')
max_compression = 0.00844;
for i=1:length(spring_k)
    linspring=linear_spring(spring_k(i),0,spring_k(i)*max_compression);
    %exspring=linear_spring(k_val(i),m_s,F_spring_max);
    [sol,~]=solve_model(motors(1),motors(1),load,latch,linspring);
    um_force(i)=min(sol(:,11));
    disp(['row ' num2str(i) ' of ' num2str(length(spring_k))]);
end
plot(spring_k,um_force,'.')
title("spring k vs um force","Interpreter","latex")
xlabel("spring k","interpreter","latex")
ylabel("min um force","interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
hold off



% negativeArray = zeros(3,3);
% %% call solve_model
% % make a directory for the current run
% dateString = string(datetime);wel
% cleanDateString = regexprep(dateString, " ", "_");
% cleanDateString = regexprep(cleanDateString, ":", "_");
% cleanDateString = "negative_force_tests_" + cleanDateString;
% for motorIndex=1:length(motors)
%     for springIndex=1:length(springs)
%         [sol, transition_times] = solve_model(motors(motorIndex),motors(motorIndex), load, ... 
%         latch, springs(springIndex),cleanDateString);
%         negativeArray(motorIndex, springIndex) = any(sol(:,11) < 0)
%     end
% end
