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
        max_force, range, velocity, Force
    end
    
    methods
        % constructor
        function obj = Motor(max_force, range, velocity, Force)
            obj.Force  = Force;
            obj.velocity = velocity;
            obj.range  = range;
            obj.max_force = max_force;
         end
    end
     
end