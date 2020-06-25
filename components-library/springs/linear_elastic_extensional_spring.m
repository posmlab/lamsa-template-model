%% linear_elastic_extensional_spring struct
% arguments in required order:
%     E - modulus
%     A - cross sectional area
%     L - length
%     rho - density (optional)
%     sigma_f - failure strength in Pascals (optional)
% min # arguments = 1

function spring = linear_elastic_extensional_spring(E,A,L,varargin)
    varargin_param_names = {'rho','sigma_f'};
    varargin_default_values = {0,Inf};
    
    % check and assign optional parameters
    if (nargin < 1)
        error('Linear elastic extensional spring requires at least 1 argument.');
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

    k = (E*A)/L;
    m = rho*A*L;
    F_spring_max = sigma_f*A;
    lin_spring = linear_spring(k,m,F_spring_max);
    
    spring.Force = lin_spring.Force;  
    spring.mass = m;
end