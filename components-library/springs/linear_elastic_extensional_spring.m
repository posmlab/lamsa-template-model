%% linear_elastic_extensional_spring struct
% arguments in required order:
%     E - modulus
%     A - cross sectional area
%     L - length
%     rho - density (optional)
%     sigma_f - failure strength in Pascals (optional)
% min # arguments = 0

function spring = linear_elastic_extensional_spring(varargin)
    %% linear_elastic_extensional_spring input parameters
    spring.param_names = {'E','A','L','rho','sigma_f'};
    spring.param_default_values = {1,1,1,0,Inf};
    
    %% General code for spring components
    spring.param_values = containers.Map(spring.param_names,spring.param_default_values,'UniformValues',false);
    
    % check and assign optional parameters
    if (length(varargin)>length(spring.param_names))
        error('Too many input parameters');
    end
    
    % loop over input arguments and assign them to spring.param_values container Map
    for i=1:length(varargin)
        spring.param_values(spring.param_names{i}) = varargin{i};
    end
    
    % loop over parameter names and store them outside the container map
    % (making a local copy to make their use in any equations more human
    % readable)
    for i=1:length(spring.param_names)
        eval([spring.param_names{i} ' = spring.param_values(''' spring.param_names{i} ''')';]);
    end    

    %% Defining the linear_elastic_extensional_spring Force, mass, and range
    k = (E*A)/L;
    m = rho*A*L;
    F_spring_max = sigma_f*A;
    lin_spring = linear_spring(k,m,F_spring_max);
    
    spring.Force = lin_spring.Force;  
    spring.mass = m;
    spring.range=F_spring_max/k;
end