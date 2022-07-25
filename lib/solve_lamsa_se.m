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

if spring.mass == 0
    warning("Spring is massless. Solving using direct actuation...")
    [sol,transition_times] = solve_direct_actuation(loading_motor, load);
    cols = size(sol,1);
    sol = [sol zeros(cols,8)];
    return
end

initial_conditions = zeros(6,1);
initial_conditions(2) = load.theta_0;
initial_conditions(3) = latch.v_0;
initial_conditions(5) = 1e-10;
options = odeset('Events', @(t,y) launching_end(t,y), "RelTol", 1e-8);
odeprob = @(t,y) se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring);

[t,y,~,~,~] = ode15s(odeprob, tspan, initial_conditions, options);

l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L2 = load.lengths(2);

fSpring = zeros(size(t));
fUnlatchingMotor = zeros(size(t));
F_comp = zeros(size(t,1),4);

for i = 1:size(t,1)
    fSpring(i) = f_perp(t(i), y(i,1), y(i,2), y(i,5), y(i,6), load, loading_motor, l0, spring);
    fUnlatchingMotor(i) = unlatching_motor.Force(t(i), [y(i,4),y(i,3)]);
    
    mu = latch.coeff_fric;
    n =  normal_force(y(i,1), y(i,3), y(i,4), fSpring(i), load, latch, fUnlatchingMotor(i));
    phi = atan(latch.y_L{2}(y(4))); %angle of latch surface
    
    F_comp(i,1) = n*sin(phi);    %normal force on latch
    F_comp(i,2) = n*L2*cos(phi); %normal torque on load
    F_comp(i,3) = mu*n*cos(phi); %frictional force on latch
    F_comp(i,4) = mu*L2*sin(phi); %frictional toque on load
end

sol=[t y(:,2) y(:,1) y(:,4) y(:,3) F_comp fSpring fUnlatchingMotor y(:,6) y(:,5)];

[~, argmaxv] = max(y(:,1));


latched = y(:,4) - latch.max_width;
transition_times = [t(find(latched > 0, 1)), t(argmaxv)];  

if (nargin >= 7)
    writeInfoToFile(load.mass, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory);
end
end

function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%
%   y = [theta dot, theta, s dot, s, y1dot, y1]
if imag(y(1))~= 0
   warning("Complex Numbers") 
   y = real(y);
end
dydt = zeros(6,1);

l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L1 = load.lengths(1);
L2 = load.lengths(2);
theta0 = load.theta_0;
moI = load.mass;
mu = latch.coeff_fric;
msp = spring.mass;

beta = sqrt(2*L1^2*(1-cos(y(2)-theta0)) + l0^2 - 2*l0*L1*(sin(y(2))- sin(theta0)));
gamma = (L1^2*sin(y(2)-theta0) - l0*L1*cos(y(2)))/beta;
y2 = l0 - beta;
y2dot = -gamma*y(1);
% y2 = L1*sin(y(2));
% y2dot = L1*cos(y(2))*y(1);
y2ddot = -gamma*dydt(1) - y(1)^2*(L1^2 * cos(y(2)-theta0) + l0*L1*sin(y(2)) + gamma^2)/beta;
%y2ddot = L1*(cos(y(2))-sin(y(2))*y(1)^2);
alpha = asin(L1*(cos(theta0) -  cos(y(2)))/(l0 - y2)); %Angle the spring makes with the vertical
Fsp =  spring.Force(t, [y2 - y(6), y2dot - y(5)]);
Flm = loading_motor.Force(t, [y(6), y(5)]); %Loading Motor force
Fperp = f_perp(t, y(1), y(2), y(5), y(6), load, loading_motor, l0, spring);
%Fperp = (Flm - msp/2 * ((3/msp * Flm - 3/msp * Fsp - y2ddot/2) + y2ddot)) * sin(pi/2 - y(2) - alpha); %spring force perpendicular to lever
phi = atan(latch.y_L{2}(y(4))); %angle of latch surface
Funlatch = unlatching_motor.Force(t, [y(4),y(3)]);
n =  normal_force(y(1), y(3), y(4), Fperp, load, latch, Funlatch);

dydt(1) = 1/moI * (Fperp*L1 - n*L2*cos(phi) - mu*L2*sin(phi));
dydt(2) = y(1);
dydt(3) = (-mu*n*cos(phi) + n*sin(phi) + unlatching_motor.Force(t, [y(4),y(3)]) )/latch.mass;
dydt(4) = y(3);
dydt(5) = 3/msp * Flm - 3/msp * Fsp - y2ddot/2;
dydt(6) = y(5);
%disp(dydt(5))


if Fperp >= 5
    disp("uh oh")
end 

if dydt(3) == 0 && dydt(4) == 0
   warning("Latch is Stuck")
end

end

function n = normal_force(thetadot, sdot, s, Fperp, load, latch, Ful)
%NORMAL_FORCE is the normal force of the latch on the lever
L1 = load.lengths(1);
L2 = load.lengths(2);
df = latch.y_L{2}(s);
ddf = latch.y_L{3}(s);

% checks if latch is out of the way or if it is moving faster than the
% lever
if s < latch.max_width && df*sdot <= L2*(thetadot+1e-3)
    moI = load.mass;
    mu = latch.coeff_fric;
    mL = latch.mass;
    phi = atan(df);
    
    
    n = (moI*df*Ful + mL*moI*ddf*sdot*sdot - mL * Fperp * L1)/( (moI*df*mu - mL * L2)*cos(phi) - (moI*df + mL * mu * L2)*sin(phi) );
else % If latch has been removed, no more normal force
    n = 0;
end


end

function f = f_perp(t, thetadot, theta, y1dot, y1, load, loading_motor, l0, spring)
L1 = load.lengths(1);
theta0 = load.theta_0;

beta = sqrt(2*L1^2*(1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta)- sin(theta0)));
gamma = (L1^2*sin(theta-theta0) - l0*L1*cos(theta))/beta;
y2 = l0 - beta;
y2dot = gamma*thetadot;
y2ddot = -gamma*thetadot - thetadot^2*(L1^2 * cos(theta-theta0) + l0*L1*sin(theta) + gamma^2)/beta;
msp = spring.mass;
Flm = loading_motor.Force(t, [y2 - y1, y2dot - y1dot]);
Fsp = spring.Force(t, [y2 - y1, y2dot - y1dot]);
%y2 = L1*sin(theta);
%y2dot = L1*cos(theta)*thetadot;

alpha = asin(L1*(cos(theta0) -  cos(theta))/(l0 - y2)); %Angle the spring makes with the vertical

f = 1/5 * (-Flm + 6*Fsp - msp*y2ddot) * sin(pi/2 - theta - alpha); %spring force perpendicular to lever
end


function [position,isterminal,direction] = launching_end(t,y)
position = y(1);
isterminal = 0;
direction = 0;

end