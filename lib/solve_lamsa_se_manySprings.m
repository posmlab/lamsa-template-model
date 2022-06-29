function  [sol, transition_times] = solve_lamsa_se_manySprings(tspan, loading_motor,unlatching_motor,load,latch,spring)
%SOLVE_LAMSA_SE_MANYSPRINGS Solves equations of motion for series elastic
%system, approximating the spring as a series of masses connected by
%springs with appropriately higher spring constants

NUM_STUFF = 25; % # of elements in y vector; increase by 2 to increase no. of springs by 1
initial_conditions = zeros(NUM_STUFF,1);
initial_conditions(2) = load.theta_0;
options = odeset('Events', @(t,y) launching_end(t,y)); %Stops solving when angular velocity is zero
odeprob = @(t,y) se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring);

[t,y,~,~,~] = ode15s(odeprob, tspan, initial_conditions, options);

l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L2 = load.lengths(2);
L3 = load.lengths(3);

fSpring = zeros(size(t));
fUnlatchingMotor = zeros(size(t));
F_comp = zeros(size(t,1),4);

for i = 1:size(t,1)
    fSpring(i) = f_perp(t(i), y(i,1), y(i,2), y(i,5), y(i,6), load, spring, l0);
    fUnlatchingMotor(i) = unlatching_motor.Force(t(i), [y(i,4),y(i,3)]);
    
    mu = latch.coeff_fric;
    n =  normal_force(t(i), y(i,1), y(i,2), y(i,3), y(i,4), y(i,5), y(i,6), unlatching_motor, load, latch, spring, l0);
    phi = atan(latch.y_L{2}(y(4))); %angle of latch surface
    
    F_comp(i,1) = n*sin(phi);    %normal force on latch
    F_comp(i,2) = n*L2*cos(phi); %normal torque on load
    F_comp(i,3) = mu*n*cos(phi); %frictional force on latch
    F_comp(i,4) = mu*L2*sin(phi); %frictional toque on load
end

sol=[t y(:,2).*L3 y(:,1).*L3 y(:,4) y(:,3) F_comp fSpring fUnlatchingMotor];

[~, argmaxv] = max(y(:,1));

transition_times = [t(find(flip(F_comp(:,1)), 1)), t(argmaxv)];  
end

function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%
%   y = [theta dot, theta, s dot, s, y1dot, ...]
if imag(y(1))~= 0
   warning("Complex Numbers") 
   y = real(y);
end
nsp = size(y,1) - 3;
dydt = zeros(nsp+3,1);

small_mass = 2*latch.mass/(nsp-6);
l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L1 = load.lengths(1);
L2 = load.lengths(2);
L3 = load.lengths(3);
moI = load.mass*(L3^3 + L2^3)/(3*(L3 + L2));% moment of inertia assuming uniform mass
mu = latch.coeff_fric;

fPerp = f_perp(t, y(1), y(2), y(5), y(6), load, spring, l0);
n =  normal_force(t, y(1), y(2), y(3), y(4), y(5), y(6), unlatching_motor, load, latch, spring, l0);
phi = atan(latch.y_L{2}(y(4))); %angle of latch surface
[y2, y2dot] = y_2(y(1), y(2), load, l0);

dydt(1) = 1/moI * (fPerp*L1 - n*L2*cos(phi) - mu*L2*sin(phi));
dydt(2) = y(1);
dydt(3) = (-mu*n*cos(phi) + n*sin(phi) + unlatching_motor.Force(t, [y(4),y(3)]) )/latch.mass;
dydt(4) = y(3);
dydt(5) = (spring.Force(t, [y(6) - y(8), y(5) - y(7)]) - spring.Force(t, [y2-y(6), y2dot - y(5)]))/small_mass;
dydt(6) = y(5);
for i = 7:2:nsp
   dydt(i) = (spring.Force(t, [y(i+1) - y(i+3), y(i) - y(i+2)]) - spring.Force(t, [y(i-1) - y(i+1), y(i-2) - y(i)]))/small_mass;
   dydt(i+1) = y(i);
end
dydt(nsp+2) = (loading_motor.Force(t, [y(nsp+3), y(nsp+2)]) - spring.Force(t, [y(nsp+1) - y(nsp+3), y(nsp) - y(nsp+2)]))/small_mass;
dydt(nsp+3) = y(nsp+2);
end


function f = f_perp(t, thetadot, theta, y1dot, y1, load, spring, l0)
%F_PERP is the part of the spring force perpendicular to the lever

[y2, y2dot] = y_2(thetadot, theta, load, l0);
% y1dot = y_1dot(t, thetadot, theta, y1, loading_motor, load, spring, l0);
a = alpha(theta, load, l0);

%if theta + a >= pi/2 %Why?
%    f = 0;
%else
f = spring.Force(t, [y2 - y1, y2dot - y1dot]) * sin(pi/2 - theta - a);
%end

end

function a = alpha(theta, load, l0)
%ALPHA is the angle the spring and muscle make with the vertical

[y2, ~] = y_2(0, theta, load, l0);
theta0 = load.theta_0;
L1 = load.lengths(1);

a = asin(L1*(cos(theta0) -  cos(theta))/(l0 - y2));

end

function n = normal_force(t, thetadot, theta, dsdt, s, y1dot, y1, unlatching_motor, load, latch, spring, l0)
%NORMAL_FORCE is the normal force of the latch on the lever

if s < latch.max_width
    Fperp = f_perp(t, thetadot, theta, y1dot, y1, load, spring, l0);
    df = latch.y_L{2}(s);
    ddf = latch.y_L{3}(s);
    Ful = unlatching_motor.Force(t, [s, dsdt]);
    mL = latch.mass;
    mu = latch.coeff_fric;
    L1 = load.lengths(1);
    L2 = load.lengths(2);
    L3 = load.lengths(3);
    moI = load.mass*(L3^3 + L2^3)/(3*(L3 + L2));% moment of inertia assuming uniform mass
    phi = atan(df);
    
    n = (moI*df*Ful + mL*moI*ddf*dsdt*dsdt - mL * Fperp * L1)/( (moI*df*mu - mL * L2)*cos(phi) - (moI*df + mL * mu * L2)*sin(phi) );
else % If latch has been removed, no more normal force
    n = 0;
end


end


function [y2, y2dot] = y_2(thetadot, theta, load, l0)
%Y2 is the distance the spring and muscle have contracted and its time
%derivative
L1 = load.lengths(1);
theta0 = load.theta_0;

%y2 = l0 - sqrt(L1^2*(cos(theta) - cos(theta0))^2 + (l0 - L1*(sin(theta) - sin(theta0)) )^2 );
y2 = l0 - sqrt(2*L1^2 * (1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta) - sin(theta0)));

y2dot = -thetadot*(  L1^2*sin(theta-theta0) - l0*L1*cos(theta) )   /   sqrt(2*L1^2 * (1-cos(theta-theta0)) + l0^2 - 2*l0*L1*(sin(theta) - sin(theta0)));

if y2 >= l0
    error("Spring and Muscle have zero length")
elseif (imag(y2) ~= 0) || (imag(y2dot) ~= 0)
    warning("Complex valued y2")
    y2 = real(y2);
    y2dot = real(y2dot);
end

end


function [position,isterminal,direction] = launching_end(t,y)
position = y(1)+1e-1;
isterminal = 1;
direction = 0;

end