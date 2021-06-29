%% Spring class definition

% Spring is a superclass for spring objects

% subclass structure:
%
% Spring
%       ExponentialSpring
%               ExponentialElasticExtensionalSpring
%       LinearSpring
%               LinearElasticExtensionalSpring

classdef Spring
    
    % common properties for all types of springs
    properties
        Force, mass, range, rest_length
    end
    
    methods
        % constructor
        function obj = Spring(Force, mass, range, rest_length)
            obj.Force  = Force;
            obj.mass = mass;
            obj.range  = range;
            obj.rest_length = rest_length;
         end
    end
     
end