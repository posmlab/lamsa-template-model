%% Motor class definition

% Motor is a superclass for motor objects

% subclass structure:
%
% Motor
%       LinearMotor
%       HillMuscleMotor
%       DeactivatingMotor

classdef Motor
    
    % common properties for all types of motors
    properties
        max_force, range, velocity, rest_length, Force, damping
    end
    
    methods
        % constructor
        function obj = Motor(max_force, range, velocity, Force, rest_length, damping)
            obj.Force  = Force;
            obj.velocity = velocity;
            obj.range  = range;
            obj.max_force = max_force;
            obj.rest_length = rest_length;
            obj.damping = damping;
        end
    end
    
end