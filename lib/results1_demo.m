close all
clearvars
tic
debug = false;
addpath(genpath(pwd));

output_directory = create_output_directory();

load_time_constraint=Inf;

metrics={'Pmax', 'vto', 'KEmax'};


%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 4;
loading_motor_range_of_motion = 5;
v_max_loading_motor = 1;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 10E-3;
loading_motor_r_activation = Inf;

% loading motor struct initialization
loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor2 = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
m=25E-3;
% load mass struct initialization
load = load_mass(m);

%% latch

% latch parameters
R=5E-3;
m_L= 3E-3;
coeff_fric = 0;
v_0L=0;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% spring paramters
k = 1; % k or k_0 depending on linear or exponential spring
m_s=2E-3;
F_spring_max=20;
k_opt=loading_motor.max_force/loading_motor_muscle_length;
k_opt=4000;

% extra parameters for exponential spring
% should be a negative value
characteristic_length = 1E-3;

% spring struct initialization
spring = linear_spring(k_opt, m_s, F_spring_max);
spring2 = exponential_spring(k_opt, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
unlatching_motor_range_of_motion = 5E-3;
F_max_unlatching_motor=.25;
v_max_unlatching_motor=1;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 5;
unlatching_motor_r_activation = .4;

unlatching_motor = hill_muscle_motor(unlatching_motor_muscle_length, F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_r_activation);
unlatching_motor2= linear_motor(F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_range_of_motion);


%% end editable parameters

%% start time series plots 

% %plot force outputs of motor as a function of time HILL MUSCLE components
% [sol,~]=solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
% sz=size(sol);
% musclesol=sol;
% for t = 1:sz(1)
%     force_array(t)=unlatching_motor.max_force*unlatching_motor.F_length(musclesol(t,1),[musclesol(t,4), musclesol(t,5)])*unlatching_motor.F_velocity(musclesol(t,1),[musclesol(t,4), musclesol(t,5)])*unlatching_motor.F_activation(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
%     force_activation(t)=unlatching_motor.F_activation(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
%     force_velocity(t)=unlatching_motor.F_velocity(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
%     force_length(t)=unlatching_motor.F_length(musclesol(t,1),[musclesol(t,4), musclesol(t,5)]);
%     max_force=unlatching_motor.max_force;
%     
% end
% figure
% hold on
% plot(musclesol(:,1),force_array,'r');
% plot(musclesol(:,1),force_activation,'k');
% plot(musclesol(:,1),force_velocity,'b');
% plot(musclesol(:,1),force_length,'m');
% plot(musclesol(:,1),max_force,'g');
% hold off


%performance metrics as function of time for given inputs 
%compares 2 runs over every column of sol 
%code for iterating over k_vals 
k_val=[k_opt/8, k_opt*2];
for i = 1:length(k_val)
      %linspringarr(i)=linear_spring(k_val(i), m_s, F_spring_max)
      nonlinspringarr(i)=exponential_spring(k_val(i),characteristic_length, m_s, F_spring_max);
end
figure
hold on
for l = 1:length(nonlinspringarr)  
%     [sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
%     solutionset1=sol;
%     sz1=size(solutionset1);
    [sol,~]=solve_model(loading_motor2,unlatching_motor2,load,latch,nonlinspringarr(l),output_directory);
    solutionset2=sol;
    columntitles=["Time", "y postition", "y velocity", "x position", "x velocity", "normal force on latch x", ...
        "normal force on load y", "frictional force on latch x", ...
        "frictional force on load y", "spring force", ...
        "unlatching motor force into"];
    %forces: x normal BLUE RED
%     y normal BLUE RED
%     y friction GREEN magenta 
%     x friction GREEN magenta 
%     spring cyan black 
%     unlatching cyan black
    n=1;
    col=["k","b","r","g"];
    units=[" [s]"," [m]"," [m/s]"];
    for i = 2:3
        subplot(3,2,n)
        box on
        set(gca,'TickLabelInterpreter','latex')
        hold on
        plot(solutionset2(:,1),solutionset2(:,i),col(l))
%         plot(solutionset1(:,1),solutionset1(:,i),"b")
        hold off
        title(columntitles(i),"Interpreter","latex");
        legend("k$<$k$_{opt}$","k$>$k$_{opt}$","Interpreter","latex")
        ylabel(columntitles(i)+units(i),"Interpreter","latex");
        xlabel(columntitles(1)+units(1),"Interpreter","latex");
        n=n+2;
    end
    n=2;
    for i=4:5
         subplot(3,2,n);
         box on
         set(gca,'TickLabelInterpreter','latex')
         hold on
         plot(solutionset2(:,1),solutionset2(:,i),col(l))
%          plot(solutionset1(:,1),solutionset1(:,i),"b")
         hold off
         title(columntitles(i),"Interpreter","latex");
         legend("k$<$k$_{opt}$","k$>$k$_{opt}$","Interpreter","latex")
         ylabel(columntitles(i)+units(i-2),"Interpreter","latex");
         xlabel(columntitles(1)+" [s]","Interpreter","latex");
         n=n+2;
    end
    col2=["Blue","Green","c","Red","m","k"];
    p=[1 4];
    for i=[6 8 11]
        %force plotting on x plot 
        subplot(3,2,6)
        box on
        set(gca,'TickLabelInterpreter','latex')
        hold on
        plot(solutionset2(:,1),solutionset2(:,i),col2(p(l)))
        p=p+1;
    end
    title("Force Components X","Interpreter","latex");
    ylabel("Force [N]","Interpreter","latex");
    legend("normal force","friction","unlatching","normal force 2","friction 2","unlatching 2","interpreter",'latex')
    xlabel(columntitles(1)+" [s]","Interpreter","latex");
    hold off
    q=[1 4];
    for i=[7 9 10]
        %force plotting on y plot 
        subplot(3,2,5)
        box on
        set(gca,'TickLabelInterpreter','latex')
        hold on
        plot(solutionset2(:,1),solutionset2(:,i),col2(q(l)));
        q=q+1;
    end
    title("Force Components Y","Interpreter","latex");
    ylabel("Force [N]","Interpreter","latex");
     legend("normal force","friction","spring","normal force 2","friction 2","spring 2","interpreter",'latex')
    xlabel(columntitles(1)+" [s]","Interpreter","latex");
    
    hold off
hold off

end
%








%% Figure with
%power vs. time (four examples)2 of each m*sol(:,3).*gradient(sol(:,3))./gradient(sol(:,1))
k_val=[k_opt/8,k_opt*2,k_opt];
for i = 1:length(k_val)
      expospringarr(i)=exponential_spring(k_val(i),characteristic_length, m_s, F_spring_max);
end
max_muscle_force = [4,8];
for i = 1:length(max_muscle_force)
    musclearr(i)=hill_muscle_motor(loading_motor_muscle_length, max_muscle_force(i), v_max_loading_motor, loading_motor_r_activation);
end
figure
hold on
for val=1:2
    subplot(2,2,val)
    set(gca,'TickLabelInterpreter','latex')
    [sol,~]=solve_model(musclearr(1),unlatching_motor2,load,latch,expospringarr(val),output_directory);
    solu1=sol;
    plot(solu1(:,1),load.mass*solu1(:,3).*gradient(solu1(:,3))./gradient(solu1(:,1)))
end
for val=1:2
    subplot(2,2,val+2)
    set(gca,'TickLabelInterpreter','latex')
    [sol,~]=solve_model(musclearr(val),unlatching_motor2,load,latch,expospringarr(3),output_directory);
    solu2=sol;
    plot(solu2(:,1),load.mass*solu2(:,3).*gradient(solu2(:,3))./gradient(solu2(:,1)))
end
subplot(2,2,1)
title("muscle max="+max_muscle_force(1)+ " [N], k="+k_val(1)+" [N/m]","Interpreter","latex");
xlabel("Time [s]","interpreter","latex")
ylabel("Power [W]","interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,2)
xlabel("Time [s]","interpreter","latex")
ylabel("Power [W]","interpreter","latex")
title("muscle max="+max_muscle_force(1)+ " [N], k="+k_val(2)+" [N/m]","Interpreter","latex");
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,3)
xlabel("Time [s]","interpreter","latex")
ylabel("Power [W]","interpreter","latex")
title("muscle max="+max_muscle_force(1)+ " [N], k="+k_val(3)+" [N/m]","Interpreter","latex");
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,4)
title("muscle max="+max_muscle_force(2)+" [N], k="+k_val(3)+" [N/m]","Interpreter","latex");
xlabel("Time [s]","interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
ylabel("Power [W]","interpreter","latex")
hold off 





%(1) 1D plot: power vs. stiffness in exponential spring model
%   loop over k values  

% k_val=linspace(k_opt/10,k_opt*20,250);
% figure
% for i=1:length(k_val)
%     exspring=exponential_spring(k_val(i),characteristic_length, m_s, F_spring_max);
%     %exspring=linear_spring(k_val(i),m_s,F_spring_max);
%     [sol,~]=solve_model(loading_motor2,unlatching_motor2,load,latch,exspring,output_directory);
%     max_power(i)=max(load.mass*sol(:,3).*gradient(sol(:,3))./gradient(sol(:,1)));
% end
% plot(k_val,max_power,'.')
% title("max power vs. stiffness","Interpreter","latex")
% xlabel("k value [N/m]","interpreter","latex")
% ylabel("max power [W]","interpreter","latex")
% set(gca,'TickLabelInterpreter','latex')
% hold off
% 
% 
% 
% %(2) 1D plot: sweep through max muscle force
% muscle_max=linspace(0,80,250);
% figure
% for i=1:length(muscle_max)
%     motor=hill_muscle_motor(loading_motor_muscle_length, muscle_max(i), v_max_loading_motor, loading_motor_r_activation);
%     [sol,~]=solve_model(motor,unlatching_motor2,load,latch,spring2,output_directory);
%     max_power1(i)=max(load.mass*sol(:,3).*gradient(sol(:,3))./gradient(sol(:,1)));
% end
% plot(muscle_max,max_power1,'.')
% title("max power vs. max muscle force","Interpreter","latex")
% xlabel("max muscle force [N]","interpreter","latex")
% ylabel("max power [W]","interpreter","latex")
% set(gca,'TickLabelInterpreter','latex')
% hold off


%(3) 2D plot: combining the two
N=100;

% setting x axis on the plot (Fmax of latch)
xname = 'k value [N/m]';
xrange = [100 10000];
k_val = linspace(xrange(1),xrange(2),N);

%setting y axis value on plot (Vmax of latch)
yname = 'muscle max force [N]';
yrange = [1E-1 10];
muscle_max_f = linspace(yrange(1),yrange(2),N);

for ii=1:length(metrics)
    outval{ii}=zeros(N);
end
for i=1:N
    for j=1:N
        motorarray(i)=hill_muscle_motor(loading_motor_muscle_length, muscle_max_f(i), v_max_loading_motor, loading_motor_r_activation);
        springarray(j)=exponential_spring(k_val(j),characteristic_length, m_s, F_spring_max);
        [sol,transition_times]=solve_model(motorarray(i),unlatching_motor2,load,latch,springarray(j),output_directory);
        met_dict=get_metrics(sol,transition_times,load ,metrics);
        for ii=1:length(metrics)
            outval{ii}(i,j)=met_dict(metrics{ii});
        end 
    end
   disp(['row ' num2str(i) ' of ' num2str(N)]);
end
%%
for ii=1:length(metrics)
    figure
    imagesc(xrange,yrange,outval{ii});
    set(gca,'YDir','normal','TickLabelInterpreter','latex');
    xlabel(xname,'Interpreter', 'Latex');
    ylabel(yname, 'Interpreter', 'Latex');
    c = colorbar;
    set(c,'TickLabelInterpreter','latex')
    c.Label.String = metrics{ii}+"[W]";
    c.Label.Interpreter="latex";
end

% power max vs muscle force DA model
max_f = linspace(1E-1, 10, N);
max_power=zeros(size(max_f));
for i = 1:length(max_f)
   motor=hill_muscle_motor(loading_motor_muscle_length, max_f(i), v_max_loading_motor, loading_motor_r_activation);
   [sol, transition_times] = solve_direct_actuation(motor,load);
   max_power(i)=max(load.mass*sol(:,3).*gradient(sol(:,3))./gradient(sol(:,1))); 
end
figure
hold on
plot(max_f,max_power)
title("Pmax vs. max muscle force", "interpreter","latex")
xlabel("max muscle force [N]", "interpreter", "latex")
ylabel("power [W]", "interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
hold off
            
            
            
%% figure with examples of loading from changing stiffness and force   



%comparing linear and exponential spring force vs displacement for kopt,
%k<k_opt


figure 
plotspot=1;

for k = [k_opt/8,k_opt, k_opt*2]
    lin_spring=linear_spring(k, m_s, F_spring_max);
    expo_spring=exponential_spring(k, characteristic_length, m_s,F_spring_max);
    [sol,~]=solve_model(loading_motor2,unlatching_motor2,load,latch,lin_spring,output_directory);
    linsol=sol;
    linsz=size(linsol);
    [sol,~]=solve_model(loading_motor2,unlatching_motor2,load,latch,expo_spring,output_directory);
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
    subplot(1,3, plotspot); 
    sgtitle("$k_{opt}$="+k_opt +" [N/m]","Interpreter","Latex");
    hold on
    title("k="+k +" [N/m]","Interpreter","Latex");
    xlabel('muscle length [m]',"Interpreter","Latex");
    ylabel('F [N]',"Interpreter","Latex");
    plot(loading_motor_muscle_length+linsol(:,2),linspring(1:lin_range(1)),'r');
    plot(loading_motor_muscle_length+exposol(:,2),expospring(1:expo_range(1)),'k');
    y=linspace(0,2*loading_motor_muscle_length,500);
    for i = 1:length(y)
        loadforce(i)=loading_motor2.Force(Inf,[loading_motor_muscle_length-y(i), 0]);
    end
    plot(y,loadforce);
    hold off
%     lin_spring_work=trapz(linsol(:,2),linspring(1:lin_range(1)))
%     expo_spring_work=trapz(exposol(:,2),expospring(1:expo_range(1)))
    plotspot=plotspot+1;
end  
%% 


    
    
    
    
    