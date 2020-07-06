
function [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
%Solve set of differential equations for loading, unlatching, and launching
%   phases of LAMSA motion
LARGE_NUM = 1E10; % subtract a large number in fzero function to trick fzero into identifying points where motor or spring suddnely go to 0
m_eff = load.mass + spring.mass/3;


%% Loading phase: Fs vs. Fin
% Assume quasistatic loading, end position is when Fs=Fin

% finding an order-of-magnitude initial guess for fzero call from motor properties
y_list = logspace(-20,20,40); % sweep through 40 orders of magnitude
F_list = zeros(size(y_list));
for i = 1:length(y_list)
    F_list(i) = loading_motor.Force(Inf,[y_list(i) 0]);
end
y_guess_motor = -y_list(find(F_list>0,1,'last'));

% initial guess based on initial spring stiffness
y_guess_spring = -Inf;
delta = eps;
while (y_guess_spring == -Inf)
    y_guess_spring = -loading_motor.Force(Inf,[0 0])/((spring.Force(0,-delta)-spring.Force(0,0))/(delta));
    delta = 10*delta;
    if (delta > 1000000)
        error('Unable to evaluate stiffness at 0 position')
    end
end

% use fzero to find when Fs=Fin 
y_guess = max([y_guess_motor, y_guess_spring]);
% options = optimset('Display','iter');
options = {};
[y0,~,exitflag]=fzero(@(y) (loading_motor.Force(Inf,[-y 0])-spring.Force(0,[y 0])) - LARGE_NUM*((~loading_motor.Force(Inf,[-y 0]))||(~spring.Force(0,[y 0])))+LARGE_NUM*(y>0),y_guess,options);
if (exitflag<0)
    error('fzero failed');
end


% checks latching distance conditions
if (abs(y0) < latch.min_latching_dist)
    warning('Loading failed. Does not fall within latching distance conditions.');
    sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    if (nargin >= 6)
    writeInfoToFile(m_eff, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
    end
    return
elseif (y0 > latch.max_latching_dist)
    y0 = latch.max_latching_dist;
end

%% Unlatching phase: Fs vs Flatch

if (unlatching_motor.max_force==0 && latch.v_0 == 0)
    warning('Latch gets stuck!');
     sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    return
end
try
    [inst_check,~,~]=unlatching_end(0,[0,latch.v_0],m_eff,y0,latch,spring,unlatching_motor);
catch ('Latch gets stuck!');
    warning('Latch gets stuck!')
    sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    if (nargin >= 6)
    writeInfoToFile(m_eff, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
    end
    return
end
if inst_check>0 
    unlatch_opts=odeset('Events',@(t,y) unlatching_end(t,y,m_eff,y0,latch,spring,unlatching_motor),'RelTol',1E-7,'AbsTol',1E-10);
    ode=@(t,y) unlatching_ode(t,y,m_eff,y0,latch,spring,unlatching_motor);
    
    a_0L = unlatching_motor.max_force / latch.mass;
    
    if (a_0L ~= 0)
        % calculate t_L_guess using quadratic formula
        % and the following kinematic equation: R = (1/2)a*t^2 + v_0*t  
        t_L_guess = (((-1*latch.v_0) + sqrt((latch.v_0)^2  + (2*a_0L*latch.max_width)))/(a_0L));
    elseif (latch.v_0 ~= 0 )
        t_L_guess = latch.max_width/latch.v_0;
    else
        warning("The latch's initial velocity and acceleration are both zero.")
        sol = [0,0,0]
        transition_times = [0,0]
        return
    end
    
    try
        tspan=linspace(0,t_L_guess,1000);
        [t_unlatch,x_unlatch]=ode45(ode,tspan,[0 latch.v_0], unlatch_opts);
    catch ("Latch gets stuck!");
        warning('Latch gets stuck!')
        % if the latch gets stuck, just give back
        % the initial conditions
        sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
        transition_times = [inf,inf];
        if (nargin >= 6)
        writeInfoToFile(m_eff, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
        end
        return 
    end
    
    % this while loop ensures that the system unlatches.
    % t_L_guess is a guess at the upper bound on the unlatching time.
    % Usually, integrating from t=0 to t=t_L_guess is a long enough
    % time interval to trigger the 'unlatching_end' event option. 
    % However, when it's not, we simply try longer and longer time
    % intervals until we get a long enough time interval that 
    % t_unlatch(end) ~= tspan(end), which indicates that the 
    % integration stopped early because we've activated 'unlatching_end'
    while (t_unlatch(end) == tspan(end))
        t_L_guess = 10 * t_L_guess;
        tspan = linspace(0, t_L_guess,1000);
        [t_unlatch,x_unlatch]=ode45(ode,tspan,[0 latch.v_0],unlatch_opts);
    end
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
    x_unlatch = [0,0];
end
t_unlatch = real(t_unlatch);
y_unlatch = real(y_unlatch);


%% Solving for Normal Force in the unlatching phase
for i=1:size(x_unlatch, 1)% For derivation of this equation for F_n see Overleaf doc with LaMSA derivation
    num1 = (latch.mass*spring.Force(t_unlatch(i), y_unlatch(i, :))) - ...
        (m_eff*latch.y_L{3}(x_unlatch(i,1))*(x_unlatch(i,2)^2)*latch.mass) - ...
        (unlatching_motor.Force(t_unlatch(i), x_unlatch(i,:))*m_eff*latch.y_L{2}(x_unlatch(i,1)));
    rad = 1 + ((latch.y_L{2}(x_unlatch(i,1)))^2);
    num2 = sqrt(rad);
    den1 = m_eff*latch.y_L{2}(x_unlatch(i,1))*(latch.y_L{2}(x_unlatch(i,1)) - latch.coeff_fric);
    den2 = latch.mass*(1+latch.coeff_fric*latch.y_L{2}(x_unlatch(i,1)));
    F_n(i) =(num1*num2)/(den1 + den2);%filling in the F_n vector until unlatch time
end
F_n = F_n';%switching to a column vector so we can add it to sol
%disp(F_n)

%% Components of Normal Force And Frictional Force
% Currently not working, some trig or possibly t_L issues

% Defining the geometric definitions of sine and cosine
for i=1:size(x_unlatch, 1)
    den = sqrt(1 + (latch.y_L{2}(x_unlatch(i, 1))^2));
    sin_comp(i) = (latch.y_L{2}(x_unlatch(i, 1)))/den;
    cos_comp(i) = 1/den;
end

sin_comp = sin_comp';
cos_comp = cos_comp';

%Calculating Normal and Frictional Force Components
for i=1:size(x_unlatch, 1)
    F_nx(i) = F_n(i) .* sin_comp(i);
    F_ny(i) = F_n(i) .* cos_comp(i);
    F_fx(i) = F_nx(i) * latch.coeff_fric;
    F_fy(i) = F_ny(i) * latch.coeff_fric;
end
F_nx = F_nx';
F_ny = F_ny';
F_fx = F_fx';
F_fy = F_fy';

F_comp = [F_nx F_ny F_fx F_fy];


%% Ballistic phase:Fs only
%guess launch times by treating the spring as ideal-ish and getting the
%   frequency
stiffness = abs(( spring.Force(0,y_unlatch(end,:))-spring.Force(0,y_unlatch(end,:)+(100*eps)) ) / (100*eps)); %Here be divide by 0 errors, probably
%stiffness=abs(spring.Force(0,y_unlatch(end,:))/y_unlatch(end,1)); %Here be divide by 0 errors, probably
nat_freq=sqrt(stiffness/m_eff);
t_launch_guess=pi/nat_freq;
launch_opts=odeset('Events',@(t,y) launching_end(t,y,spring));
ode=@(t,y) launching_ode(t,y,m_eff,spring);
tspan=linspace(0,t_launch_guess,1E3);
y0=y_unlatch(end,:)';
[t_launch,y_launch]=ode45(ode,tspan,y0,launch_opts);

% Solve latch dynamics during Ballistic Phase
%     Currently assuming instaneous stopping of latch

%% Stitch together solutions
T=[t_unlatch;t_unlatch(end)+t_launch];
Y=[y_unlatch;y_launch];
x_launch = repmat([x_unlatch(end, 1), 0], size(y_launch,1), 1);%Makes X the right size to fit in sol

X=[x_unlatch;x_launch];

transition_times=[t_unlatch(end),t_unlatch(end)+t_launch(end)];

for i = 1:size(T)
    fSpring(i) = spring.Force(T(i), [Y(i,1), Y(i,2)]);%fill out the fSpring vector to add to sol
    fUnlatchingMotor(i) = unlatching_motor.Force(T(i), [X(i,:)]);
end
fSpring = fSpring';
fUnlatchingMotor = fUnlatchingMotor';
% add zeros to the end of F_comp because F_comp consists of 
% forces during the unlatching phase, and the other vectors
% include forces during the ballistic phase. 
% Adding zeros makes this matrix the right size for appending 
% to the rest of the sol.

F_comp = [F_comp; zeros(size(T,1)-size(F_comp,1),4)];

% stitch together various numbers 
% for one big matrix to write to csv file 
sol=[T Y X F_comp fSpring fUnlatchingMotor];


%% Establishing Parameters for .json output
if (nargin >= 6)
writeInfoToFile(m_eff, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
end
end

