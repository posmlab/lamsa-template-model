function [sol,transition_times] = solve_lamsa_se(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
%SOLVE_LAMSA_SE Solves equations of motion for series elastic system
%   

% Section 1: Odeproblem

sol = [0];
transition_times = [0];
end

function dydt = se_ode(t, y, loading_motor, unlatching_motor, load, latch, spring)
%SE_ODE is the equation of motion for a series elastic system
%   
%   y = [theta dot, theta, s dot, s]
dydt = zeros(4,1);

dydt(1) = 
dydt(2) = y(1);
dydt(3) = 
dydt(4) = y(3);
end

function f_perp()

end

function alpha()

end