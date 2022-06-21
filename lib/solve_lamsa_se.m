function [sol,transition_times] = solve_lamsa_se(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
%SOLVE_LAMSA_SE Solves equations of motion for series elastic system

sol = [0];
transition_times = [0];
end

function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%   
%   y = [theta dot, theta, s dot, s, y1]
dydt = zeros(4,1);

moI = load.mass;
n = normal_force;
L1 = load.L1;
L2 = L1/load.EMA;
phi = atan(latch.y_L);
mu = latch.coeff_fric;

dydt(1) = 1/moI * (f_perp*L1 - n*L2*cos(phi) - mu*L2*sin(phi));
dydt(2) = y(1);
dydt(3) = (-mu*n*cos(phi) + n*sin(phi) + unlatching_motor.Force)/latch.mass;
dydt(4) = y(3);
dydt(5) = 
end

function f = f_perp()
%F_PERP is the part of the spring force perpendicular to the lever
a = alpha;
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

function f = normal_force()
%NORMAL_FORCE is the normal force of the latch on the lever
end

function y2 = y2()
%Y2 is the distance the spring and muscle have contracted

end

function y2dot = y2dot()
%Y2DOT is the time derivative of y2

end

function y1dot = y1dot()
%Y1DOT is the velocity of the point of intersection of the muscle and
%spring along the direction of the muscle.
% It is calculated using the force-velocity curve of the muscle


end