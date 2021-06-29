%% Latch class definition

% Latch is a superclass for latch objects

% subclass structure:
%
% Latch
%       RoundedLatch
%       ParabolicLatch

classdef Latch
    
    % common properties for all types of latches
    properties
        coeff_fric, v_0, max_width, mass, y_L, min_latching_dist, max_latching_dist, runway_length
    end
    
    methods
        % constructor
        function obj = Latch(coeff_fric, v_0, max_width, mass, y_L, min_latching_dist, max_latching_dist, runway_length)
            % making a dropoff of the form -10^6(x-a)^2+b
            % before the runway end so the latch can't
            % slip backwards

            % a = .5*y{2}(0)*10^-6;
            % b = y{1}(0)+(.5*y{2}(0))^2*10^-6;
            
            y = y_L;
            yL = @(x) (-10^6*(x-.5*y{2}(0)*10^-6)^2+(y{1}(0)+(.5*y{2}(0))^2*10^-6))*(x<0) + y{1}(x)*(x>=0);
            yL_prime = @(x) (-2*10^6*(x-.5*y{2}(0)*10^-6))*(x<0) + y{2}(x)*(x>=0);
            yL_doubleprime = @(x) (-2*10^6)*(x<0) + y{3}(x)*(x>=0);
            
            obj.coeff_fric = coeff_fric;
            obj.v_0 = v_0;
            obj.max_width = max_width;
            obj.mass = mass;
            obj.y_L = {yL, yL_prime, yL_doubleprime};
            obj.min_latching_dist = min_latching_dist;
            obj.max_latching_dist = max_latching_dist;
            obj.runway_length = runway_length;
        end
    end
end