%% Latch class definition

% Latch is a superclass for latch objects

% subclass structure:
%
% Latch
%       RoundedLatch

classdef Latch
    
    % common properties for all types of latches
    properties
        coeff_fric, v_0, max_width, mass, y_L, min_latching_dist, max_latching_dist, runway_length
    end
    
    methods
        % constructor
        function obj = Latch(coeff_fric, v_0, max_width, mass, y_L, min_latching_dist, max_latching_dist, runway_length)
            obj.coeff_fric = coeff_fric;
            obj.v_0 = v_0;
            obj.max_width = max_width;
            obj.mass = mass;
            obj.y_L = y_L;
            obj.min_latching_dist = min_latching_dist;
            obj.max_latching_dist = max_latching_dist;
            obj.runway_length = runway_length;
        end
    end
end
