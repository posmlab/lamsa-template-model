
function [sol,transition_times] = solve_model(loading_motor,unlatching_motor,m_eff,latch,spring)
%Solve set of differential equations for loading, unlatching, and launching
%   phases of LAMSA motion
LARGE_NUM = 1E10; % subtract a large number in fzero function to trick fzero into identifying points where motor or spring suddnely go to 0
F_s = spring.Force;
F_in = loading_motor.Force;
F_out = unlatching_motor.Force;
coeff_fric = latch.coeff_fric;


%% Loading phase: Fs vs. Fin
% Assume quasistatic loading, end position is when Fs=Fin

% finding an order-of-magnitude initial guess for fzero call from motor properties
y_list = logspace(-20,20,40); % sweep through 40 orders of magnitude
F_list = zeros(size(y_list));
for i = 1:length(y_list)
    F_list(i) = F_in(0,[y_list(i) 0]);
end
y_guess_motor = -y_list(find(F_list>0,1,'last'));

% initial guess based on initial spring stiffness
y_guess_spring = F_in(0,[0 0])/((F_s(0,eps)-F_s(0,0))/eps);

% use fzero to find when Fs=Fin
y_guess = max([y_guess_motor, y_guess_spring]);
options =  {};% optimset('Display','iter');
[y0,~,exitflag]=fzero(@(y) (F_in(0,[y 0])-F_s(0,[y 0])) - LARGE_NUM*((~F_in(0,[y 0]))||(~F_s(0,[y 0])))+LARGE_NUM*(y>0),y_guess,options);
if (exitflag<0)
    error('fzero failed');
end

%% Unlatching phase: Fs vs Flatch

%check for time independence
if spring.Time_independent == true && unlatching_motor.Time_independent == true
    F_friction = coeff_fric*F_s(0,[y0,0]);
    if F_out(0,[0 , 0]) < F_friction
        warning('latch cannot overcome friction force');
        sol = [0,0,0];
        transition_times = [0,0];
        return
    end 
else %there is time dependence
    warning('time dependence in latch or motor');
end    

[inst_check,~,~]=unlatching_end(0,[0,latch.v_0],m_eff,latch.mass,F_s,F_out,y0,latch.y_L, coeff_fric);
if inst_check>0 
    unlatch_opts=odeset('Events',@(t,y) unlatching_end(t,y,m_eff,latch.mass,F_s,F_out,y0,latch.y_L, coeff_fric),'RelTol',1E-7,'AbsTol',1E-10);
    ode=@(t,y) unlatching_ode(t,y,m_eff,latch.mass,F_s,F_out,y0,latch.y_L, coeff_fric);
    
    a_0L = F_out(0,[0 0]) / latch.mass;
    
    if (a_0L ~= 0)
        % calculate t_L_guess using quadratic formula
        % and the following kinematic equation: R = (1/2)a*t^2 + v_0*t  
        t_L_guess = (((-1*latch.v_0) + sqrt((latch.v_0)^2  + (2*a_0L*latch.max_width)))/(a_0L));
    elseif (latch.v_0 ~= 0 )
        t_L_guess = latch.max_width/latch.v_0        
    else
        warning("The latch's initial velocity and acceleration are both zero.")
        sol = [0,0,0]
        transition_times = [0,0]
        return
    end
    
    tspan=linspace(0,t_L_guess,1000);%[0,t_L_guess];
    [t_unlatch,x_unlatch]=ode45(ode,tspan,[0 latch.v_0],unlatch_opts);
    % This ODE is for the latch x-coordinate, but we want the y-coordinate, so
    % convert
    y_unlatch=zeros(size(x_unlatch));
    for i=1:length(y_unlatch)
        y_unlatch(i,1)=latch.y_L{1}(x_unlatch(i,1))+y0;
        y_unlatch(i,2)=x_unlatch(i,2)*latch.y_L{2}(x_unlatch(i,1));
    end
    if (imag(y_unlatch))
        disp(y_unlatch);
        error('y unlatch imaginary');
    end
else % instantaneous unlatching
    y_unlatch=[y0,0]; %May cause a repeated time step and give NaNs on differentiation
    t_unlatch=[0];
end
t_unlatch = real(t_unlatch);
y_unlatch = real(y_unlatch);


%% Ballistic phase:Fs only
%guess launch times by treating the spring as ideal-ish and getting the
%   frequency
stiffness=abs(F_s(0,y_unlatch(end,:))/y_unlatch(end,1)); %Here be divide by 0 errors, probably
nat_freq=sqrt(stiffness/m_eff);
t_launch_guess=pi/nat_freq;
launch_opts=odeset('Events',@(t,y) launching_end(t,y,F_s));
ode=@(t,y) launching_ode(t,y,F_s,m_eff);
tspan=linspace(0,t_launch_guess,1E3);
y0=y_unlatch(end,:)';
[t_launch,y_launch]=ode45(ode,tspan,y0,launch_opts);

%% Stitch together solutions
transition_times=[t_unlatch(end),t_unlatch(end)+t_launch(end)];
T=[t_unlatch;t_unlatch(end)+t_launch];
Y=[y_unlatch;y_launch];
sol=[T Y];
end

