%% Mass class definition

% Mass is a superclass for mass objects

% subclass structure:
%
% Mass
%       OneDMass
%       RotatinMass

classdef Mass
    
    % common properties shared among all mass objects
    properties
        mass, EMA
    end
    
    methods
        % constructor
        function obj = Mass(mass, EMA)
            obj.mass = mass;
            obj.EMA = EMA;
        end
    end
end

