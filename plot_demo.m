close all
clearvars
tic
debug = false;
addpath(genpath(pwd));

dateString = string(datetime);
cleanDateString = regexprep(dateString, " ", "_");
cleanDateString = regexprep(cleanDateString, ":", "_");
        
load_time_constraint=Inf;


%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 20;
loading_motor_range_of_motion = 5;
v_max_loading_motor = 100.0000;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 5;
loading_motor_r_activation = Inf;

% loading motor struct initialization
loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor2 = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
m=100;
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
k = 1; % k or k_0 depending on linear or exponential spring
m_s=1;
F_spring_max=1E4;
k_opt=loading_motor.max_force/loading_motor_range_of_motion;

% extra parameters for exponential spring
% should be a negative value
characteristic_length = -5;

% spring struct initialization
spring = linear_spring(k, m_s, F_spring_max);
spring2 = exponential_spring(k, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
unlatching_motor_range_of_motion = 5;
F_max_unlatching_motor=10;
v_max_unlatching_motor=100;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 5;
unlatching_motor_r_activation = .4;

unlatching_motor = hill_muscle_motor(unlatching_motor_muscle_length, F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_r_activation);
unlatching_motor2= linear_motor(F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_range_of_motion);

%% end editable parameters

%% start time series plots 

%plot force outputs of motor as a function of time HILL MUSCLE components
[sol,~]=solve_model(loading_motor,unlatching_motor,load,latch,spring,cleanDateString);
sz=size(sol);
musclesol=sol;
for t = 1:sz(1)
    force_array(t)=unlatching_motor.max_force*unlatching_motor.F_length(musclesol(t,1),[musclesol(t,4), musclesol(t,5)])*unlatching_motor.F_velocity(musclesol(t,1),[musclesol(t,4), musclesol(t,5)])*unlatching_motor.F_activation(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    force_activation(t)=unlatching_motor.F_activation(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    force_velocity(t)=unlatching_motor.F_velocity(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    force_length(t)=unlatching_motor.F_length(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
    max_force=unlatching_motor.max_force;
    
end
figure
hold on
plot(musclesol(:,1),force_array,'r');
plot(musclesol(:,1),force_activation,'k');
plot(musclesol(:,1),force_velocity,'b');
plot(musclesol(:,1),force_length,'m');
plot(musclesol(:,1),max_force,'g');
hold off



%performance metrics as function of time for given inputs 
%compares 2 runs over every column of sol 
%code for iterating over k_vals 
k_val=[k_opt/5,k_opt*5];
for i = 1:length(k_val)
      %linspringarr(i)=linear_spring(k_val(i), m_s, F_spring_max)
      nonlinspringarr(i)=exponential_spring(k_val(i),characteristic_length, m_s, F_spring_max);
end
figure
hold on
for l = 1:length(nonlinspringarr)  
%     [sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring,cleanDateString);
%     solutionset1=sol;
%     sz1=size(solutionset1);
    [sol,~]=solve_model(loading_motor,unlatching_motor2,load,latch,nonlinspringarr(l),cleanDateString);
    solutionset2=sol;
    columntitles=["Time", "y", "ydot", "x", "xdot", "normal force on latch x", ...
        "normal force on load y", "frictional force on latch x", ...
        "frictional force on load y", "spring force", ...
        "unlatching motor force into"];
    n=1;
    col=["k","b","r","g"];
    for i = 2:3
        subplot(3,2,n)
        hold on
        plot(solutionset2(:,1),solutionset2(:,i),col(l))
%         plot(solutionset1(:,1),solutionset1(:,i),"b")
        hold off
        title(columntitles(i));
        ylabel(columntitles(i));
        xlabel(columntitles(1));
        n=n+2;
    end
    n=2;
    for i=4:5
         subplot(3,2,n);
         hold on
         plot(solutionset2(:,1),solutionset2(:,i),col(l))
%          plot(solutionset1(:,1),solutionset1(:,i),"b")
         hold off
         title(columntitles(i));
         ylabel(columntitles(i));
         xlabel(columntitles(1));
         n=n+2;
    end
    for i = 6
    %for i=6:2:10
        %force plotting on x plot 
        subplot(3,2,6)
        hold on
        plot(solutionset2(:,1),solutionset2(:,i),col(l))
    end
    title("Force Components Y");
    ylabel("Force");
    xlabel(columntitles(1));
    hold off
    for i=7
    %for i=7:2:11
        %force plotting on y plot 
%         subplot(3,2,5)
%         hold on
%         plot(solutionset2(:,1),solutionset2(:,i),col(l))
    end
    title("Force Components X");
    ylabel("Force");
    xlabel(columntitles(1));
    hold off
hold off
end
%








%% Figure with
            %power vs. time (four examples)2 of each
            %(1) 1D plot: power vs. stiffness in exponential spring model
            %   loop over k values  
            %(2) 1D plot: sweep through max muscle force
            %(3) 2D plot: combining the two
            
            
            
            
            
            
            
            
%% figure with examples of loading from changing stiffness and force   



%comparing linear and exponential spring force vs displacement for kopt,
%k<k_opt


figure 
plotspot=1;

for k = [k_opt/5,k_opt*2]
    lin_spring=linear_spring(k, m_s, F_spring_max);
    expo_spring=exponential_spring(k, characteristic_length, m_s,F_spring_max);
    [sol,~]=solve_model(loading_motor2,unlatching_motor2,load,latch,lin_spring,cleanDateString);
    linsol=sol;
    linsz=size(linsol);
    [sol,~]=solve_model(loading_motor2,unlatching_motor2,load,latch,expo_spring,cleanDateString);
    exposol=sol;
    exposz=size(exposol);
    expospring=5+zeros(exposz(1));
    for i = 1:exposz(1)
        expospring(i)=expo_spring.Force(exposol(i,1),[exposol(i,2), exposol(i,3)]);
    end
    for i = 1:linsz(1)
        linspring(i)=lin_spring.Force(linsol(i,1),[linsol(i,2), linsol(i,3)]);
    end
    expo_range=size(exposol(:,2));
    lin_range=size(linsol(:,2));
    subplot(1,2, plotspot); 
    sgtitle("$k_{opt}$="+k_opt,"Interpreter","Latex");
    hold on
    title("k="+k,"Interpreter","Latex");
    xlabel('disp.',"Interpreter","Latex");
    ylabel('F',"Interpreter","Latex");
    plot(linsol(:,2),linspring(1:lin_range(1)),'r');
    plot(exposol(:,2),expospring(1:expo_range(1)),'k');
    for i = 1:exposz(1)
        loadforce(i)=loading_motor2.Force(exposol(i,1),[exposol(i,2), exposol(i,3)]);
    end
    plot(exposol(:,2),loadforce);
    hold off
%     lin_spring_work=trapz(linsol(:,2),linspring(1:lin_range(1)))
%     expo_spring_work=trapz(exposol(:,2),expospring(1:expo_range(1)))
    plotspot=plotspot+1;
end  
%% 


    
    
    
    
    