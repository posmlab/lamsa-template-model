
function [sol,transition_times] = solve_lamsa(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)

% Solve set of differential equations for loading, unlatching, and launching
% phases of LAMSA motion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading phase

N_points = 10000;

rangeBound = min([latch.max_latching_dist, spring.range, loading_motor.range, 0.1*realmax]);

% this is the range of y values where we are looking for Fspring =
% Floading_motor.
% We want to sweep over a wide range of values initially, 
% to get a rough sense of the order-of-magnitude y-value where
% the forces are equal 
y_range=-logspace(log10(eps), log10(rangeBound),N_points);

fmotor = zeros(size(y_range));
fspring = zeros(size(y_range));
for i=1:length(y_range)
    fmotor(i)=loading_motor.Force(Inf,[-y_range(i) 0]);
    fspring(i)=spring.Force(0,[y_range(i), 0]);
end
fdiff = fmotor - fspring;
index=find(fdiff<0,1,"first");
if index == 1
    y0=0;
elseif isempty(index)
    [~,ind]=max(fspring);
    if ind == length(y_range)
        y0=y_range(end);
    else
        y_range=linspace(y_range(ind),y_range(ind+1),N_points);
        fspring=zeros(size(y_range));
        for i=1:length(y_range)
            fspring(i)=spring.Force(0,[y_range(i) 0]);
        end
        [~,ind]=max(fspring);
        y0=y_range(ind);
    end
else 
    y_range=linspace(y_range(index-1),y_range(index),N_points);
    fdiff = zeros(size(y_range));
    for i=1:length(y_range)
        fdiff(i)=loading_motor.Force(Inf,[-y_range(i) 0])-spring.Force(0,[y_range(i) 0]);
    end
    index=find(fdiff<0,1,"first");
    y0=y_range(index);
end

% checks latching distance and frictional loading conditions
if (abs(y0) < latch.min_latching_dist)
    warning('Loading failed. Does not fall within latching distance conditions.');
    y0 = 0;
    sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    if (nargin >= 6)
        writeInfoToFile(m_eff(load,spring,y0,y0), transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
    end
    return
end
if isa(unlatching_motor, 'DeactivatingMotor')
    angle = atan(latch.y_L{2}(0));
    mu = latch.coeff_fric;
    numer = mu*cos(angle)-sin(angle);
    denom = cos(angle)+mu*sin(angle);
    new_max = F_eff(load, spring, y0, Inf, [y0,0])*numer/denom;
    condition = new_max < unlatching_motor.max_force; %forces are both negative
    unlatching_motor.max_force = new_max;
    unlatching_motor.Force = @(t,x) new_max * exp(-unlatching_motor.r_deactivation*t);
else
    condition = latch.coeff_fric < latch.y_L{2}(0);
end
if (condition)
    warning('Loading failed. Latch slope is too large at x=0 for the given coefficient of friction.');
    y0 = 0;
    sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    if (nargin >= 6)
        writeInfoToFile(m_eff(load,spring,y0,y0), transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
    end
    return
elseif (abs(y0) > latch.max_latching_dist)
    y0 = -latch.max_latching_dist;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unlatching phase

if (unlatching_motor.max_force==0 && latch.v_0 == 0)
    warning('Latch gets stuck!');
     sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
            latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
            unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    return
end
try
    [inst_check,~,~]=unlatching_end(0,[0,latch.v_0],load,y0,latch,spring,unlatching_motor);
catch ME
    switch ME.message
        case 'Latch gets stuck!'
            warning('Latch gets stuck!')
            sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
                    latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
                    unlatching_motor.Force(0,[0,0])];
            transition_times = [inf,inf];
            if (nargin >= 6)
                writeInfoToFile(m_eff(load,spring,y0,y0), transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
            end
            return
        otherwise
            rethrow(ME)
    end
end
if inst_check>0  % not instantaneous unlatching
    unlatch_opts=odeset('Events',@(t,y) unlatching_end(t,y,load,y0,latch,spring,unlatching_motor),...
        'RelTol',1E-6,'AbsTol',1E-9,'OutputFcn', @(t,y,flag) update_f_unlatching(t,y,flag,spring,latch,y0));
    ode=@(t,y) unlatching_ode(t,y,load,y0,latch,spring,unlatching_motor);
    
    a_0L = abs(unlatching_motor.max_force / latch.mass);
    
    if (a_0L ~= 0)
        % calculate t_L_guess using quadratic formula
        % and the following kinematic equation: R = (1/2)a*t^2 + v_0*t  
        t_L_guess = (((-1*latch.v_0) + sqrt((latch.v_0)^2  + (2*a_0L*latch.max_width)))/(a_0L));
    elseif (latch.v_0 ~= 0 )
        t_L_guess = latch.max_width/latch.v_0;
    else
        warning("The latch's initial velocity and acceleration are both zero.")
        sol = [0,0,0];
        transition_times = [0,0];
        return
    end
    
    try
        tspan=linspace(0,t_L_guess,1000);
        [t_unlatch,x_unlatch]=ode45(ode,tspan,[0 latch.v_0], unlatch_opts);
    catch ME
        switch ME.message
            case "Latch gets stuck!"
            warning('Latch gets stuck!')
            % if the latch gets stuck, just give back
            % the initial conditions
            sol = [0,y0,0,0,0,0,spring.Force(0,[y0, 0]), ...
                latch.coeff_fric*spring.Force(0,[y0, 0]),0,spring.Force(0,[y0,0]), ...
                unlatching_motor.Force(0,[0,0])];
            transition_times = [inf,inf];
            if (nargin >= 6)
                writeInfoToFile(m_eff(load,spring,y0,y0), transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
            end
            return
            otherwise
                rethrow(ME)
        end
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
        y_unlatch(i,1)=(latch.y_L{1}(x_unlatch(i,1)) - latch.y_L{1}(0)) + y0;
        y_unlatch(i,2)=x_unlatch(i,2)*latch.y_L{2}(x_unlatch(i,1));
    end
    if (imag(y_unlatch))
        disp(y_unlatch);
        error('y unlatch imaginary');
    end
else % instantaneous unlatching
    y_unlatch=[y0,0]; %May cause a repeated time step and give NaNs on differentiation
    t_unlatch= 0;
    x_unlatch = [0,0];
end
t_unlatch = real(t_unlatch);
y_unlatch = real(y_unlatch);


%% Solving for Normal Force in the unlatching phase
F_n = zeros(size(x_unlatch, 1),1);
for i=1:size(x_unlatch, 1)% For derivation of this equation for F_n see Overleaf doc with LaMSA derivation
    y=(latch.y_L{1}(x_unlatch(i,1))-latch.y_L{1}(0))+y0;
    num1 = (latch.mass*F_eff(load,spring,y0,t_unlatch(i), y_unlatch(i, :))) - ...
        (m_eff(load,spring,y0,y)*latch.y_L{3}(x_unlatch(i,1))*(x_unlatch(i,2)^2)*latch.mass) - ...
        (unlatching_motor.Force(t_unlatch(i), x_unlatch(i,:))*m_eff(load,spring,y0,y)*latch.y_L{2}(x_unlatch(i,1)));
    rad = 1 + ((latch.y_L{2}(x_unlatch(i,1)))^2);
    num2 = sqrt(rad);
    den1 = m_eff(load,spring,y0,y)*latch.y_L{2}(x_unlatch(i,1))*(latch.y_L{2}(x_unlatch(i,1)) - latch.coeff_fric);
    den2 = latch.mass*(1+latch.coeff_fric*latch.y_L{2}(x_unlatch(i,1)));
    F_n(i) =(num1*num2)/(den1 + den2);%filling in the F_n vector until unlatch time
end

%% Components of Normal Force And Frictional Force
% Currently not working, some trig or possibly t_L issues

% Defining the geometric definitions of sine and cosine
sin_comp = zeros(size(x_unlatch,1),1);
cos_comp = zeros(size(x_unlatch,1),1);
for i=1:size(x_unlatch, 1)
    den = sqrt(1 + (latch.y_L{2}(x_unlatch(i, 1))^2));
    sin_comp(i) = (latch.y_L{2}(x_unlatch(i, 1)))/den;
    cos_comp(i) = 1/den;
end

%Calculating Normal and Frictional Force Components
F_nx = zeros(size(x_unlatch, 1),1);
F_ny = zeros(size(x_unlatch, 1),1);
F_fx = zeros(size(x_unlatch, 1),1);
F_fy = zeros(size(x_unlatch, 1),1);
for i=1:size(x_unlatch, 1)
    F_nx(i) = F_n(i) .* sin_comp(i);
    F_ny(i) = F_n(i) .* cos_comp(i);
    F_fx(i) = -1*F_ny(i) * latch.coeff_fric;
    F_fy(i) = -1*F_nx(i) * latch.coeff_fric;
end
F_comp = [F_nx F_ny F_fx F_fy];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ballistic phase:Fs only
%guess launch times by treating the spring as ideal-ish and getting the
%   frequency
stiffness = abs(( spring.Force(0,y_unlatch(end,:))-spring.Force(0,y_unlatch(end,:)+(100*eps)) ) / (100*eps)); 
nat_freq=sqrt(stiffness/m_eff(load,spring,y0,y_unlatch));
t_launch_guess=pi/nat_freq;
launch_opts=odeset('Events',@(t,y) launching_end(t,y,load,spring,y0),'OutputFcn', @(t,y,flag) update_f_launching(t,y,flag,spring));
ode=@(t,y) launching_ode(t,y,load,spring,y0);
tspan=linspace(0,t_launch_guess,1E3);
y_init=y_unlatch(end,:)';
[t_launch,y_launch]=ode45(ode,tspan,y_init,launch_opts);

%% Stitch together solutions
T=[t_unlatch;t_unlatch(end)+t_launch];
Y=[y_unlatch;y_launch];
x_launch = repmat([x_unlatch(end, 1), 0], size(y_launch,1), 1);%Makes X the right size to fit in sol

X=[x_unlatch;x_launch];

transition_times=[t_unlatch(end),t_unlatch(end)+t_launch(end)];

% check to make sure that the spring didn't return to equilibrium during
% unlatching and that our solution only contains y<=0
index = find(Y(1:end-1,1)>0,1,'first');
if (~isempty(index))
    T = T(1:index);
    Y = Y(1:index,:);
    X = X(1:index,:);
    transition_times = [T(end), T(end)];
end

fSpring = zeros(size(T,1),1);
fUnlatchingMotor = zeros(size(T,1),1);
for i = 1:size(T)
    fSpring(i) = F_eff(load,spring,y0,T(i),[Y(i,1),Y(i,2)]);%fill out the fSpring vector to add to sol
    if (T(i) < t_unlatch(end))
        fUnlatchingMotor(i) = unlatching_motor.Force(T(i), X(i,:));
    else
        fUnlatchingMotor(i) = 0;
    end
    
end

% add zeros to the end of F_comp because F_comp consists of
% forces during the unlatching phase, and the other vectors
% include forces during the ballistic phase.
% Adding zeros makes this matrix the right size for appending
% to the rest of the sol.

F_comp = [F_comp; zeros(size(T,1)-size(F_comp,1),4)];
F_comp = F_comp(1:size(T,1),:);

% stitch together various numbers
% for one big matrix to write to csv file
sol=[T Y X F_comp fSpring fUnlatchingMotor];

%% Establishing Parameters for .json output
if (nargin >= 6)
writeInfoToFile(m_eff(load,spring,y0,y0), transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
end

end

function dx = launching_ode(t,y,load,spring,y0)
%ODE for the launching phase
dx(1)=y(2,1);
dx(2)=F_eff(load, spring, y0, t, y)/m_eff(load, spring, y0, y(1));
dx=dx';
end

function [value,isterminal,direction] = launching_end(t,y,load,spring,y0)
%End condition for launching
value=F_eff(load, spring, y0, t, y);
isterminal=1;
direction=0;
end

function dx = unlatching_ode(t,x,load,y0,latch,spring,unlatching_motor)
% ODE for the unlatching phase
y=(latch.y_L{1}(x(1))-latch.y_L{1}(0))+y0;
m = m_eff(load, spring, y0, y);
yL_prime = latch.y_L{2}(x(1));
yL_doubleprime = latch.y_L{3}(x(1));
Y=[y,yL_prime*x(2)];
scaling = sign(x(2));

spring_force = F_eff(load, spring, y0, t, Y);

denom=latch.mass+m*(yL_prime^2)+scaling*latch.coeff_fric*(latch.mass*yL_prime-m*yL_prime);
num=unlatching_motor.Force(t,x)+yL_prime*spring_force+scaling*latch.coeff_fric*(unlatching_motor.Force(t,x)*yL_prime-spring_force)-(x(2)^2)*(m*yL_prime*yL_doubleprime-m*scaling*latch.coeff_fric*yL_doubleprime);
dx=[x(2);num/denom];
end

function [value,isterminal,direction] = unlatching_end(t,x,load,y0,latch,spring,unlatching_motor)

%End condition for unlatching
if (x(1) > latch.max_width)
   value = 0;
   isterminal = 1;
   direction = 0;
   return
end

yL=(latch.y_L{1}(x(1))-latch.y_L{1}(0))+y0;
yL_prime = latch.y_L{2}(x(1));
yL_doubleprime = latch.y_L{3}(x(1));
y=[yL yL_prime*x(2)];
m = m_eff(load, spring, y0, y);
spring_force = F_eff(load, spring, y0, t, y);
frictionDirection = sign(x(2));
denom=latch.mass+m*(yL_prime^2)+frictionDirection*latch.coeff_fric*(latch.mass*yL_prime-m*yL_prime);
num=unlatching_motor.Force(t,x)+yL_prime*spring_force+frictionDirection*latch.coeff_fric*(unlatching_motor.Force(t,x)*yL_prime-spring_force)-(x(2)^2)*(m*yL_prime*yL_doubleprime-m*frictionDirection*latch.coeff_fric*yL_doubleprime);
xL_doubledot=num/denom;

% solving for normal force; when it's 0, unlatching is done
value = spring_force-(m*(yL_doubleprime*(x(2)^2)+(xL_doubledot*yL_prime)));


if (imag(value))
    warning('Complex value of Normal Force. ODE step size might be too large near the end of the latch');
end


stuck_threshold = 1E-9;
if ((x(2) < stuck_threshold) && (xL_doubledot < stuck_threshold))
    if (~isa(unlatching_motor, 'DeactivatingMotor') && t==0 && spring.Force(0,[y0,0])*latch.coeff_fric > unlatching_motor.max_force)
        error('Latch gets stuck!');
    elseif (unlatching_motor.Force(t, x) >= unlatching_motor.Force(t + stuck_threshold,x))
        error('Latch gets stuck!');
    else
        warning('System is moving slowly. Integration may take a long time.')
    end
end

isterminal=1;
direction=0;
end

% calculate m_eff for rotating systems
function mass = m_eff(load, spring, y0, y)

mass = load.mass;

if isa(load, 'RotatingMass')
    y = y(1);
    if y > 0
        y = 0;
    end
    L1 = load.L1;
    theta_0 = load.theta_0;
    theta = asin((-y+y0+L1*sin(theta_0))/L1);
    mass = mass/cos(theta) + 1E2*(abs(theta) > pi/2);
end

mass = mass + spring.mass/3;

end

% calculate F_eff for rotating systems
function Force = F_eff(load, spring, y0, t, y)

Force = spring.Force(t, y);

if isa(load, 'RotatingMass')
    
    y1 = spring.rest_length + abs(y0);
    y2 = spring.rest_length + abs(y(1));
    theta_0 = load.theta_0;
    L1 = load.L1;
    m = load.mass;
    
    theta = asin((-y(1)+y0+L1*sin(theta_0))/L1);
    alpha = acos((y2^2-y1^2+2*y1*L1*sin(theta_0))/(2*y2*L1));
    
    Force = real( Force*sin(alpha) - abs(m*sin(theta)*y(2)^2/(L1*cos(theta)^3))) * (abs(theta) < pi/2);
    Force = Force * (Force > 0);
end

end

% update force history for standard linear solids
function status = update_f_unlatching(t, x, ~, spring, latch, y0)

if (isa(spring, 'StandardLinearSolid2') || isa(spring, 'StandardLinearSolid')) && (size(t,2)>0)
    t = t(end);
    x = x(:,end);
    y=(latch.y_L{1}(x(1))-latch.y_L{1}(0))+y0;
    y_prime = latch.y_L{2}(x(1))*x(2);
    F = spring.Force(t,[y,y_prime]);
    spring.F_history = [spring.F_history; [t F]];
end

status = 0;

end

% update force history for standard linear solids
function status = update_f_launching(t, y, ~, spring)

if (isa(spring, 'StandardLinearSolid2') || isa(spring, 'StandardLinearSolid')) && (size(t,2)>0)
    t = t(end);
    y = y(:,end);
    F = spring.Force(t,y);
    spring.F_history = [spring.F_history; [t F]];
end

status = 0;

end
