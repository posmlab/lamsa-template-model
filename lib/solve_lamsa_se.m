function [sol,transition_times] = solve_lamsa_se(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
%SOLVE_LAMSA_SE Solves equations of motion for series elastic system

sol = se_ode(t,y,loading_motor, unlatching_motor, load, latch, spring);
transition_times = [0];
end

function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%   
%   y = [theta dot, theta, s dot, s, y1]
dydt = zeros(4,1);

n = normal_force(t, y, unlatching_motor, load, latch);

dydt(1) = 
dydt(2) = y(1);
dydt(3) = 
dydt(4) = y(3);
dydt(5) = y1dot(t, y,loading_motor, spring);
end

function f = f_perp()
%F_PERP is the part of the spring force perpendicular to the lever
a = alpha();
    if theta + a >= 90
        f = 0;
    else
        f = @(theta) spring.Force * sin(pi/2 - theta - a);
    end
end 

function a = alpha()
%ALPHA is the angle the spring and muscle make with the vertical
l0 = spring.rest_length + loading_motor.initial_length;
theta0 = load.theta_0;
    if l0 == y2
        a = 0;
    else
        a = @(theta, y2) asind((cos(theta0) -  cos(theta))/(l0 - y2));
    end 
end

function n = normal_force(t, y, unlatching_motor, load, latch, Fperp)
%NORMAL_FORCE is the normal force of the latch on the lever
I = load.moment_of_inertia;
[~, df, ddf] = latch.f_derivs;
s = y(3);
dsdt = y(4);
Ful = unlatching_motor.Force(t, [s, dsdt]);
mL = latch.mass;
mu = latch.coeff_fric;
[L1, L2, ~] = load.lengths;

phi = atan(df);

n = (I*df*Ful + mL*I*ddf*dsdt*dsdt - mL * Fperp * L1)/( (I*df*mu - mL * L2)*cos(phi) - (I*df + mL * mu * L2)*sin(phi) );

end


function y2 = y2()
%Y2 is the distance the spring and muscle have contracted

end

function y2dot = y2dot(y, load, latch)
%Y2DOT is the time derivative of y2

end

function y1dot = y1dot(t, y, loading_motor, spring)
%Y1DOT is the velocity of the point of intersection of the muscle and
%spring along the direction of the muscle.
% It is calculated using the force-velocity curve of the muscle

spring_strain = @(y1dot) [y2() - y(5), y2dot() - y1dot];

net_force = @(y1dot) spring.Force(t, spring_strain(y1dot)) - loading_motor.Force(t, [y(5), y1dot]);

y1dot = fzero(net_force, 0);

if y1dot.isNaN()
    error("Muscle is unable to overcome spring.")
end

end