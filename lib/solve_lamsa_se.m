function  [t,y] = solve_lamsa_se(loading_motor,unlatching_motor,load,latch,spring)
%SOLVE_LAMSA_SE Solves equations of motion for series elastic system

initial_conditions = [0; load.theta_0; 0; 0; 0];
odeprob = @(t,y) se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring);

[t,y] = ode15s(odeprob, [0, 1], initial_conditions);
end

function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%
%   y = [theta dot, theta, s dot, s, y1]
dydt = zeros(5,1);

l0 = spring.rest_length + loading_motor.rest_length;% initial length of spring + muscle
L1 = load.lengths(1);
L2 = load.lengths(2);
L3 = load.lengths(3);
moI = load.mass*(L3^3 + L2^3)/(3*(L3 + L2));% moment of inertia assuming uniform mass
mu = latch.coeff_fric;

fPerp = f_perp(t, y(1), y(2), y(5), loading_motor, load, spring, l0);
n =  normal_force(t, y(1), y(2), y(3), y(4), y(5), loading_motor, unlatching_motor, load, latch, spring, l0);
phi = atan(latch.y_L{2}(y(4))); %angle of latch surface

dydt(1) = 1/moI * (fPerp*L1 - n*L2*cos(phi) - mu*L2*sin(phi));
dydt(2) = y(1);
dydt(3) = (-mu*n*cos(phi) + n*sin(phi) + unlatching_motor.Force(t, [y(4),y(3)]) )/latch.mass;
dydt(4) = y(3);
dydt(5) = y_1dot(t, y(1), y(2), y(5), loading_motor, load, spring, l0);
end

function f = f_perp(t, thetadot, theta, y1, loading_motor, load, spring, l0)
%F_PERP is the part of the spring force perpendicular to the lever

[y2, y2dot] = y_2(thetadot, theta, load, l0);
y1dot = y_1dot(t, thetadot, theta, y1, loading_motor, load, spring, l0);
a = alpha(theta, load, l0);

if theta + a >= pi/2
    f = 0;
else
    f = spring.Force(t, [y2 - y1, y2dot-y1dot]) * sin(pi/2 - theta - a);
end

end

function a = alpha(theta, load, l0)
%ALPHA is the angle the spring and muscle make with the vertical

[y2, ~] = y_2(0, theta, load, l0);
theta0 = load.theta_0;

a = asin((cos(theta0) -  cos(theta))/(l0 - y2));

end

function n = normal_force(t, thetadot, theta, dsdt, s, y1, loading_motor, unlatching_motor, load, latch, spring, l0)
%NORMAL_FORCE is the normal force of the latch on the lever

if s < latch.max_width
    Fperp = f_perp(t, thetadot, theta, y1, loading_motor, load, spring, l0);
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
L2 = load.lengths(2);
theta0 = load.theta_0;

y2 = l0 - sqrt(L1^2*(cos(theta) - cos(theta0))^2 + (l0 - L1*(sin(theta) - sin(theta0)) )^2 );

y2dot = -thetadot*(  L1^2*(sin(theta)*cos(theta0) - cos(theta)*sin(theta0)) - l0*L1*cos(theta) )   /   sqrt(  2*L1^2*(1 + sin(theta)*sin(theta0)-cos(theta)*cos(theta0)) +  l0^2  -  2*l0*L2*(sin(theta) - sin(theta0)) );

if y2 >= l0
    error("Spring and Muscle have zero length")
end

end


function y1dot = y_1dot(t, thetadot, theta, y1, loading_motor, load, spring, l0)
%Y1DOT is the velocity of the point of intersection of the muscle and
%spring along the direction of the muscle.
% It is calculated using the force-velocity curve of the muscle

[y2, y2dot] = y_2(thetadot, theta, load, l0);

spring_strain = @(y1dot) [y1 - y2, y1dot - y2dot];

net_force = @(y1dot) spring.Force(t, spring_strain(y1dot)) - loading_motor.Force(t, [y1, y1dot]);

y1dot = fzero(net_force, 10);

if isnan(y1dot)
    error("Muscle is unable to overcome spring.")
end

end