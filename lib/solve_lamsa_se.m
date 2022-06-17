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

dydt(1) = 
dydt(2) = y(1);
dydt(3) = 
dydt(4) = y(3);
dydt(5) = 
end

function f = f_perp()
%F_PERP is the part of the spring force perpendicular to the lever

end

function a = alpha()
%ALPHA is the angle the spring and muscle make with the vertical

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