%% hill_muscle_motor struct
% arguments in required order:
%     muscle_length - rest length of the muscle
%     F_motor_max - maximum amount of force the spring can exert
%     v_motor_max - maximum velocity at which the motor can travel
%     r_activation - rate of activation of the muscle
%     L_initial - initial stretch of the muscle (optional)
%     a_L, b_L, s - phenomenological parameters fitted to describe the
%     shape of the muscle’s length–tension curve (optional)
% min # arguments = 4

function motor = hill_muscle_motor(muscle_length,F_motor_max,v_motor_max,r_activation,varargin)
    % optional parameters
    varargin_param_names = {'L_initial','a_L','b_L','s'};
    varargin_default_values = {muscle_length,2.08,-2.89,-0.75};
    
    % check and assign optional parameters
    if (length(varargin)>length(varargin_param_names))
        error('Too many input parameters');
    end
    for i=1:length(varargin)
        eval([varargin_param_names{i} '=varargin{i};'])
    end
    for i=(length(varargin)+1):length(varargin_param_names)
        eval([varargin_param_names{i} '=varargin_default_values{i};'])
    end
    
    % model
    % motor.Force = F_motor_max * F_length * F_velocity * F_activation such that:
    % F_length=@(t,x) exp(-((abs((((x(1)/muscle_length)^motor.b_L)-1)/motor.s))^motor.a_L));
    % F_velocity=@(t,x)(1-(x(2)/vmax_motor))/(1+(x(2)/(vmax_motor/4))); 
    % F_activation=@(t,x)min(r_activation*t,1);
    % ^ obtained from Rosario et al.
    motor.Force = @(t,x) (x(1)<=(L_initial-(0.7*muscle_length))) * (x(1)>=(L_initial-(1.3*muscle_length))) * F_motor_max * exp(-((abs(((((L_initial-x(1))/muscle_length)^b_L)-1)/s))^a_L)) * ((1-(abs(x(2))/v_motor_max))/(1+(abs(x(2))/(v_motor_max/4)))) .* (min(r_activation*t,1));
    motor.max_force = F_motor_max;
end  

%% Citations
% Rosario MV, Sutton GP, Patek SN, Sawicki GS. 2016 Muscle–spring dynamics in time-limited, elastic movements.
%   Proc. R. Soc. B 283: 20161561. http://dx.doi.org/10.1098/rspb.2016.1561

