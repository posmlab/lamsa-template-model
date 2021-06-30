%% Mass class definition

% Mass is a superclass for mass objects

% subclass structure:
%
% Mass
%       LoadMass

classdef Mass
    
    % common properties shared among all mass objects
    properties
        mass, real_mass, EMA
    end
    
    methods
        % constructor
        function obj = Mass(mass, real_mass, EMA)
            obj.mass = mass;
            obj.real_mass = real_mass;
            obj.EMA = EMA;
        end
    end
end

