function  [sol, transition_times] = solve_lamsa_se(tspan, loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
%SOLVE_LAMSA_SE Solves equations of motion for series elastic system
%   sol is an nx13 matrix and each column corresponds to t, arclength,
%   arcvelocity, latch displacement, latch velocity, normal forces on the
%   latch and load, frictional forces on the latch and load, spring force,
%   unlatching motor force, and the displacement and velocity of the point
%   of connection between the spring and loading motor.
% 
%   transition_times is a 1x2 vector which shows the unlatching time and
%   the time to maximum displacement respectively.


% If the latch gets stuck give back initial conditions
% Still need to consider more cases when friction is involved
if (unlatching_motor.max_force==0 && latch.v_0 == 0)
    warning('Latch gets stuck!');
    sol = [0, load.theta_0, latch.v_0,0,0,0,0,0,0,0, unlatching_motor.Force(0,[0,0])];
    transition_times = [inf,inf];
    return
end

initial_conditions = zeros(6,1);
initial_conditions(2) = load.theta_0;
initial_conditions(3) = latch.v_0;

l0 = spring.rest_length + loading_motor.rest_length;
L1 = load.lengths(1);
theta0 = load.theta_0;
theta_final = atan((l0-L1)/(L1*cos(theta0)));

% If the spring is massless, we need to do something else/
if spring.mass == 0
    [sol, transition_times] = solve_massless(tspan, loading_motor, unlatching_motor, load, latch, spring, theta_final);
else  
    options = odeset('Events', @(t,y) launching_end(t,y, theta_final), 'AbsTol', 1e-8, 'RelTol', 1e-8, ...
        'OutputFcn', @(t,y,flag) update_f(t,y,flag,load,spring,loading_motor));
    odeprob = @(t,y) se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring);
    
    [t,y,~,~,~] = ode15s(odeprob, tspan, initial_conditions, options);
    
    % Solving For Fperp and normal forces
    
    L2 = load.lengths(2);
    
    F_n = zeros(size(t));
    fUnlatchingMotor = zeros(size(t));
    F_comp = zeros(size(t,1),4);
    
    for i = 1:size(t,1)
        fUnlatchingMotor(i) = unlatching_motor.Force(t(i), [y(i,4),y(i,3)]);
        
        mu = latch.coeff_fric;
        F_n(i) =  normal_force(t(i), y(i,1), y(i,2), y(i,3), y(i,4), y(i,5), y(i,6), load, latch, spring, loading_motor, fUnlatchingMotor(i));
        phi = atan(latch.y_L{2}(y(i,4))); %angle of latch surface
        
        
        F_comp(i,1) = F_n(i)*sin(phi);    %normal force on latch
        F_comp(i,2) = F_n(i)*L2*cos(phi); %normal torque on load
        F_comp(i,3) = mu*F_n(i)*cos(phi); %frictional force on latch
        F_comp(i,4) = mu*F_n(i)*L2*sin(phi); %frictional toque on load
    end
    
    fSpring = f_perp(t, y(:,1), y(:,2), y(:,3), y(:,4), y(:,5), y(:,6), F_n, load, latch, spring, loading_motor, unlatching_motor);
    
    sol=[t y(:,2) y(:,1) y(:,4) y(:,3) F_comp fSpring fUnlatchingMotor y(:,6) y(:,5)];
    
    [~, argmaxv] = max(y(:,1));
    
    latched = y(:,4) - latch.max_width;
    transition_times = [t(find(latched > 0, 1)), t(argmaxv)];
end


if (nargin >= 7)
    writeInfoToFile(load.mass, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
end

end




function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%
%   y = [theta dot, theta, s dot, s, y1dot, y1]

ul_offset = 0.05;

dydt = zeros(6,1);

% Parameters of the system
l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L1 = load.lengths(1);
L2 = load.lengths(2);
theta0 = load.theta_0;
moI = load.mass;
mu = latch.coeff_fric;
msp = spring.mass;
mL = latch.mass;
df = latch.y_L{2}(y(4));
ddf = latch.y_L{3}(y(4));

% Values that depend on the state of the system
phi = atan(df); %angle of latch surface
beta = sqrt(2*L1^2*(1-cos(y(2)-theta0)) + l0^2 - 2*l0*L1*(sin(y(2))- sin(theta0)));
gamma = (L1^2*sin(y(2)-theta0) - l0*L1*cos(y(2)))/beta;
delta = (L1^2*cos(y(2)-theta0) - l0*L1*sin(y(2)) + gamma^2)/beta;
epsilon = mu*sin(phi) - cos(phi);
epsilonbar = -mu*cos(phi) + sin(phi);
y2 = l0 - beta;
y2dot = -gamma*y(1);
alpha = asin(L1*(cos(theta0) -  cos(y(2)))/(l0 - y2)); %Angle the spring makes with the vertical
la = L1*sin(pi/2 - alpha - y(2));

% Forces
Fsp =  spring.Force(t, [y2 - y(6), y2dot - y(5)]);
Flm = loading_motor.Force(t, [y(6), y(5)]); %Loading Motor force
Fd = loading_motor.damping(t, [y(6), y(5)]);

if t > ul_offset
    Ful = unlatching_motor.Force(t-ul_offset, [y(4),y(3)]);
else
    Ful = 0;
end
F_n =  normal_force(t, y(1), y(2), y(3), y(4), y(5), y(6), load, latch, spring, loading_motor, Ful);

% Velocities
dydt(2) = y(1);
dydt(4) = y(3);
dydt(6) = y(5);

% Different ddot theta depending on if in latched or unlatched state
if y(4) < latch.max_width && F_n >= 0 %Latched
    dydt(3) = (L2*la*epsilonbar*(-2*Flm + 6*Fsp + msp*delta*y(1)^2) - (4*moI - msp*gamma*la)*epsilonbar*ddf*y(3)^2 - 4*epsilon*L2^2*Ful )...
        /((4*moI - msp*gamma*la)*df*epsilonbar - 4*epsilon*L2^2*mL);
    dydt(1) = (ddf*y(3)^2 + df*dydt(3))/L2; 
else % Unlatched
    dydt(1) = (la*(-2*Flm + 2*Fd + 6*Fsp + msp*delta*y(1)^2))/(4*moI - msp*gamma*la); 
    dydt(3) = Ful/mL;
end

y2ddot = -gamma*dydt(1) - delta*y(1)^2;
% if Fd > Flm
%     Fd = Flm;
% end 

dydt(5) = (3/msp) * (Flm - Fd - Fsp) - y2ddot/2;

% Checking if latch gets stuck
stuck_threshold = 1E-9;
if ((y(3) < stuck_threshold) && (dydt(3) < stuck_threshold) && (t > ul_offset))
    if (~isa(unlatching_motor, 'DeactivatingMotor') && t==0 && spring.Force(0,[y2, 0])*latch.coeff_fric > unlatching_motor.max_force)
        error('Latch gets stuck!');
    elseif (unlatching_motor.Force(t, [y(4) y(3)]) >= unlatching_motor.Force(t + stuck_threshold, [y(4) y(3)]))
        error('Latch gets stuck!');
    else
        warning('System is moving slowly. Integration may take a long time.')
    end
end

end



function n = normal_force(t, thetadot, theta, sdot, s, y1dot, y1, load, latch, spring, loading_motor, Ful)
%NORMAL_FORCE is the normal force of the latch on the lever
L1 = load.lengths(1);
L2 = load.lengths(2);
df = latch.y_L{2}(s);
ddf = latch.y_L{3}(s);

% checks if latch is out of the way or if it is moving faster than the
% lever. Probably somethign wrong here.
if s < latch.max_width
    
    l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
    theta0 = load.theta_0;
    moI = load.mass;
    mu = latch.coeff_fric;
    msp = spring.mass;
    mL = latch.mass;
    
    phi = atan(df); %angle of latch surface
    beta = sqrt(2*L1^2*(1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta)- sin(theta0)));
    gamma = (L1^2*sin(theta-theta0) - l0*L1*cos(theta))/beta;
    delta = (L1^2*cos(theta-theta0) - l0*L1*sin(theta) + gamma^2)/beta;
    epsilon = mu*sin(phi) - cos(phi);
    epsilonbar = -mu*cos(phi) + sin(phi);
    y2 = l0 - beta;
    y2dot = -gamma*thetadot;
    alpha = asin(L1*(cos(theta0) -  cos(theta))/(l0 - y2)); %Angle the spring makes with the vertical
    la = L1*sin(pi/2 - alpha - theta);
    
    
    Fsp =  spring.Force(t, [y2 - y1, y2dot - y1dot]);
    Flm = loading_motor.Force(t, [y1, y1dot]); %Loading Motor force
    Fd = loading_motor.damping(t, [y1, y1dot]);
%     if Fd > Flm
%         Fd = Flm;
%     end
   
    n = ((Ful * df + mL*ddf * sdot^2)*(4*moI - msp*gamma*la) - L2*la*mL*(-2*Flm + 2*Fd + 6*Fsp + msp*delta* thetadot^2))/(4*epsilon*mL*L2^2 - epsilonbar*(4*moI - msp*gamma*la)*df);

else % If latch has been removed, no more normal force
    n = 0;
end

end


function f = f_perp(t, thetadot, theta, sdot, s, y1dot, y1, n, load, latch, spring, loading_motor, unlatching_motor)
% F_PERP takes in vectors of the solution and calculates the force 
%   perpendciular to the load applied by the spring.

l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L1 = load.lengths(1);
L2 = load.lengths(2);
theta0 = load.theta_0;
num_iter = length(t);
msp = spring.mass;
mu = latch.coeff_fric;
moI = load.mass;


beta = sqrt(2*L1^2*(1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta)- sin(theta0)));
gamma = (L1^2*sin(theta-theta0) - l0*L1*cos(theta))./beta;
delta = (L1^2*cos(theta-theta0) - l0*L1*sin(theta) + gamma.^2)./beta;
y2 = l0 - beta;
y2dot = gamma.*thetadot;
alpha = asin(L1*(cos(theta0) -  cos(theta))./(l0 - y2)); %Angle the spring makes with the vertical
    
Fsp = zeros(num_iter,1);
Flm = zeros(num_iter,1);
Fd = zeros(num_iter,1);
Ful = zeros(num_iter, 1);
phi = zeros(num_iter, 1);
df = zeros(num_iter, 1);
ddf = zeros(num_iter, 1);
for i = 1:num_iter
    Fsp(i) = spring.Force(t(i), [y2(i) - y1(i), y2dot(i) - y1dot(i)]);
    Flm(i) = loading_motor.Force(t(i), [y1(i), y1dot(i)]);
    Fd(i) = loading_motor.damping(t(i), [y1(i), y1dot(i)]);
    phi(i) = atan(latch.y_L{2}(s(i))); %angle of latch surface
    df(i) = latch.y_L{2}(s(i));
    ddf(i) = latch.y_L{3}(s(i));
    Ful(i) = unlatching_motor.Force(t(i), [s(i), sdot(i)]);
end

la = L1*sin(pi/2 - alpha - theta);
ddots = (-mu*n.*cos(phi) + n.*sin(phi) + Ful)/latch.mass;

if s < latch.max_width
    thetaddot = (ddf.*sdot.^2 + df.*ddots);
else
    thetaddot = (la.*(-2*Flm + 2*Fd + 6*Fsp + msp*delta.*thetadot.^2))./(4*moI - msp*gamma.*la);
end

% if Fd > Flm
%    Fd = Flm;
% end

f =  (1/4)*(-2* Flm + 2*Fd + 6*Fsp + msp * gamma .* thetaddot + msp * delta .* thetadot.^2 ) .* sin(pi/2 - theta - alpha); %spring force perpendicular to lever

end


function [sol, transition_times] = solve_massless(tspan, loading_motor, unlatching_motor, load, latch, spring, theta_final)

dt = 1e-5;
y(1,:) = zeros(1,6);
y(1,2) = load.theta_0;
y(1,3) = latch.v_0;
F_comp = zeros(1,4);
F_perp(1) = 0;
F_unlatching_motor(1) = 0;
ul_offset = 0.05;
NUM_ITER = tspan(2)/dt;

l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L1 = load.lengths(1);
L2 = load.lengths(2);
theta0 = load.theta_0;
mu = latch.coeff_fric;
moI = load.mass;
mL = latch.mass;

theta = y(1,2);
sdot = y(1,3);
thetadot = 0;
s = 0;
y1_old = 0;

beta = sqrt(2*L1^2*(1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta)- sin(theta0)));
gamma = (L1^2*sin(theta-theta0) - l0*L1*cos(theta))/beta;
y2 = l0 - beta;
y2dot = -gamma*thetadot;

for i = 1:NUM_ITER
    
% Solving for y1

    t = i*dt;

    % Net force on spring-muscle contact as function of y1
    y1_forces = @(x) loading_motor.Force(t, [x, (x-y1_old)/dt]) - spring.Force(t, [y2 - x, y2dot - (x-y1_old)/dt]);

    y1 = fzero(y1_forces, y1_old);
   
    
    %Forces 
    Flm = loading_motor.Force(t, [y1, (y1 - y1_old)/dt]);
    
    % New value of y1
    y1dot = (y1 - y1_old)/dt;
    y1_old = y1;
 
    
% Solving theta and s

    % Backwards Euler Method
    dydt = @(x) x - [thetadot, theta, sdot, s] - dt*se_ode_massless(t, x, theta0, l0, L1, L2, mu, moI, mL, Flm, unlatching_motor, latch, ul_offset);

    options = optimset('Display','off');
    [y_new, ~, exitflag] = fsolve(dydt, [thetadot, theta, sdot, s], options);
    
    % If fsolve fails, use forwards euler method
    if exitflag < 0
       y_new = [thetadot, theta, sdot, s] + dt*se_ode_massless(t, [thetadot,  theta, sdot, s], theta0, l0, L1, L2, mu, moI, mL, Flm, unlatching_motor, latch, ul_offset); 
    end
    
    
    % Append new values
    y(i+1, :) = [y_new y1dot y1];
    
    % Update spring Force History
    update_f(t, [y_new y1dot y1], true, load, spring, loading_motor);
    
    thetadot = y_new(1);
    theta = y_new(2);
    sdot = y_new(3);
    s = y_new(4);
    
    % F_n and Fperp values for sol
    beta = sqrt(2*L1^2*(1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta)- sin(theta0)));
    gamma = (L1^2*sin(theta-theta0) - l0*L1*cos(theta))/beta;
    y2 = l0 - beta;
    y2dot = -gamma*thetadot;
    alpha = asin(L1*(cos(theta0) -  cos(theta))/(l0 - y2));
    la = L1*sin(pi/2 - alpha - theta);
    df = latch.y_L{2}(s);
    ddf = latch.y_L{3}(s);
    phi = atan(df);
    epsilon = mu*sin(phi) - cos(phi);
    epsilonbar = -mu*cos(phi) + sin(phi);

        
    if t > ul_offset
        Ful = unlatching_motor.Force(t-ul_offset, [s, sdot]);
    else
        Ful = 0;
    end
    
    F_n = (mL*la*L2*Flm - moI*df*Ful - mL*moI*ddf*y(3)^2)/(epsilonbar*moI*df + epsilon*L2^2*mL);


    F_comp(i+1,1) = F_n*sin(phi);    %normal force on latch
    F_comp(i+1,2) = F_n*L2*cos(phi); %normal torque on load
    F_comp(i+1,3) = mu*F_n*cos(phi); %frictional force on latch
    F_comp(i+1,4) = mu*F_n*L2*sin(phi); %frictional toque on load
    
    %This code is slow. Vectorize and move it outside.
    F_perp(i+1) = Flm*sin(pi/2 - alpha - theta);
    F_unlatching_motor(i+1) = Ful;
    
    if theta >= theta_final
       break; 
    end
    
    if F_perp(i+1) < 0
        warning("Warning: Perpendicular Force is negative")
    end
    

end

t = (tspan(1):dt:t)';

sol=[t y(:,2) y(:,1) y(:,4) y(:,3) F_comp F_perp' F_unlatching_motor' y(:,6) y(:,5)];

[~, argmaxv] = max(y(:,1));

latched = y(:,4) - latch.max_width;
transition_times = [t(find(latched > 0, 1)), t(argmaxv)];  

end



function dydt = se_ode_massless(t, y, theta0, l0, L1, L2, mu, moI, mL, Flm, unlatching_motor, latch, ul_offset)

beta = sqrt(2*L1^2*(1-cos(y(2)-theta0)) + l0^2 - 2*l0*L1*(sin(y(2))- sin(theta0)));
gamma = (L1^2*sin(y(2)-theta0) - l0*L1*cos(y(2)))/beta;
y2 = l0 - beta;
y2dot = -gamma*y(1);
alpha = asin(L1*(cos(theta0) -  cos(y(2)))/(l0 - y2)); 
la = L1*sin(pi/2 - alpha - y(2));
df = latch.y_L{2}(y(4));
ddf = latch.y_L{3}(y(4));
phi = atan(df);
epsilon = mu*sin(phi) - cos(phi);
epsilonbar = -mu*cos(phi) + sin(phi);

if t > ul_offset
    Ful = unlatching_motor.Force(t-ul_offset, [y(4), y(3)]);
else
    Ful = 0;
end

F_n = (mL*la*L2*Flm - moI*df*Ful - mL*moI*ddf*y(3)^2)/(epsilonbar*moI*df - epsilon*L2^2*mL);


if y(4) < latch.max_width && F_n > 0 %Latched
    dydt(3) = (-epsilonbar*la*L2*Flm + epsilonbar*moI*ddf*y(3)^2 + epsilon*L2^2*Ful)/(epsilon*L2^2*mL - epsilonbar*moI*df);
    dydt(1) = (ddf*y(3)^2 + df*dydt(3))/L2;
else % Unlatched
    dydt(1) = (la*Flm)/moI;
    dydt(3) = Ful/mL;
end

dydt(2) = y(1);
dydt(4) = y(3);

end



% update force history for standard linear solids
function status = update_f(t, y, ~, load, spring, loading_motor)

if (isa(spring, 'StandardLinearSolid2') || isa(spring, 'StandardLinearSolid')) && (size(t,2)>0)

    l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
    L1 = load.lengths(1);
    theta0 = load.theta_0;
    
    beta = sqrt(2*L1^2*(1-cos(y(2)-theta0)) + l0^2 - 2*l0*L1*(sin(y(2))- sin(theta0)));
    gamma = (L1^2*sin(y(2)-theta0) - l0*L1*cos(y(2)))./beta;
    y2 = l0 - beta;
    y2dot = gamma*y(1);
    
    F = spring.Force(t, [y2 - y(6), y2dot - y(5)]);
    spring.F_history = [spring.F_history; [t F]];
end

status = 0;

end


function [position,isterminal,direction] = launching_end(t,y, theta_final)
%  A negative value of position ends the simulation
%  theta final is the angle at which the load is parallel to the muscle and
%  spring
position = theta_final-y(2);
isterminal = 0;
direction = 0;

end