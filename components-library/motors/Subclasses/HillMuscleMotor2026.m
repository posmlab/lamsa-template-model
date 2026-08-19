%% HillMuscleMotor2026 class definition
% arguments in required order:
%     muscle_length - rest length of the muscle
%     muscle_area - cross-sectional area of the muscle
%     pennation_angle - angle at which muscle fibers attach to apodeme/tendon
%     muscle_vmax - maximum velocity at which the muscle can travel (in lengths/sec)
%     fv_curvature - curvature constant for force-velocity curve
%     min_length - smallest length at which muscle can produce force (normalized by rest length)
%     fl_curvature - curvature constant for force-length curve
%     specific_tension - maximum isometric stress of the muscle
%     r_activation - rate of activation of the muscle
%     L_initial - initial stretch of the muscle (normalized by rest length)
% min # arguments = 3

classdef HillMuscleMotor2026 < Motor
    
    methods (Static)
        % first row contains parameter names
        % second row contains default values for the loading motor
        % third row contains default values for the unlatching motor
        function parameters = parameters()
            parameters = ["muscle_length" "muscle_area" "pennation_angle" "muscle_vmax" "fv_curvature" "min_length" "fl_curvature" "specific_tension" "r_activation" "l_initial";
                "0" "0" "0" "5" "2.8" "0.5" "-8" "100" "200" "1";
                "0" "0" "0" "5" "2.8" "0.5" "-8" "100" "200" "1";
                "-Inf" "-Inf" "-Inf" "-Inf" "-Inf" "-Inf" "-Inf" "-Inf" "-Inf" "-Inf";
                "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        function obj = HillMuscleMotor2026(
            muscle_length,
            muscle_area,
            pennation_angle,
            muscle_vmax,
            fv_curvature,
            min_length,
            fl_curvature,
            specific_tension,
            r_activation,
            L_initial)
            
            % check and assign optional parameters
            if (nargin < 3)
                error('Hill muscle motor requires at least 3 arguments.');
            end
            if (length(varargin)>length(varargin_param_names))
                error('Too many input parameters');
            end
            for i=1:length(varargin)
                eval([varargin_param_names{i} '=varargin{i};'])
            end
            for i=(length(varargin)+1):length(varargin_param_names)
                eval([varargin_param_names{i} '=varargin_default_values{i};'])
            end

            % parameter range verification
            if (muscle_vmax <= 0)
                v_motor_max = 5;
                warning("muscle_vmax must be strictly positive, setting to 5 lengths/sec")
            end
            if (fv_curvature < 0)
                fv_curvature = 2.8;
                warning("fv_curvature must be positive, setting to 2.8")
            end
            if (min_length < 0) | (min_length > 1)
                min_length = 0.5;
                warning("min_length must be between 0 and 1 lengths, setting to 0.5 lengths")
            end
            if (fl_curvature < (min_length-1)^-3) | (min_length > 0)
                fl_curvature = (min_length-1)^-3;
                warning("fl_curvature must be between (min_length - 1)^-3 and 0, setting to minimum")
            end
            if (specific_tension <= 0)
                specific_tension = 100;
                warning("specific_tension must be strictly positive, setting to 100 kPa")
            end
            if (r_activation <= 0)
                r_activation = 200;
                warning("r_activation must be strictly positive, setting to 200 Hz")
            end
            if (L_initial < min_length) | (L_initial > 1)
                L_initial = 1;
                warning("min_length must be between min_length and 1 lengths, setting to 1 length")
            end

            %%% model
            % see here for equations: https://www.desmos.com/3d/4rgk6al4k1
            % where alpha is in degrees, t is in seconds, l is in relative lengths, and v is in relative lengths/sec
            % see here for full literature referenced:

            % normalized length and velocity, angle states
            l = l_initial - x(1) / muscle_length; % convert to relative length
            v = x(2) / muscle_length; % convert to relative velocity 
            alpha = atand(sind(pennation_angle) / (cosd(pennation_angle)-1+l));

            % repeated values
            fl_frac = 1/(min_length-1)^2

            Force = @(t,x) (l <= 1) % no over extension
                        * (l >= min_length) % zero negative range of fl curve
                        * (v >= 0) % no modeling of negative velocity
                        * (v <= muscle_vmax) % zero negative range of fv curve
                        * muscle_area * specific_tension % maximum isometric force
                        * cosd(alpha) % angular correction
                        * min(r_activation*t, 1) % activation curve F_time(t)
                        * (
                            (min_length*fl_curvature + 1 - fl_frac)
                            + (-(2*min_length+1)*fl_curvature + 2*fl_frac) * l
                            + ((min_length_2)*fl_curvature - fl_frac) * l^2
                            + (-fl_curvature) * l^3
                            ) % third order polynomial fl curve
                        * (1 - v/(vmax*cosd(alpha))) / (1 + fv_curvature*v/(vmax*cosd(alpha))); % fv curve

            max_force = muscle_area * specific_tension * cosd(pennation_angle);
            range= (l_initial - min_length) * muscle_length;
            velocity = muscle_vmax * muscle_length * cosd(pennation_angle);
            
            % call parent constructor
            obj = obj@Motor(max_force, range, velocity, Force, muscle_length);
        end  
    end
end

%% Citations
% Rosario MV, Sutton GP, Patek SN, Sawicki GS. 2016 Muscle�spring dynamics in time-limited, elastic movements.
%   Proc. R. Soc. B 283: 20161561. http://dx.doi.org/10.1098/rspb.2016.1561
% this is out of date!!!