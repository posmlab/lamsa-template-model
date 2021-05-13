%% Spring class definition

% Spring is a superclass for spring objects

% subclass structure:
%
% Spring
%       ExponentialSpring
%       LinearSpring
%               LinearElasticExtensionalSpring

classdef Spring
    
    % common properties for all types of springs
    properties
        Force, mass, range
    end
    
    methods
        % constructor
        function obj = Spring(Force, mass, range)
            obj.Force  = Force;
            obj.mass = mass;
            obj.range  = range;
         end
    end
     
end