%% rounded_latch struct
% arguments in required order:
%     R - raduis of the latch
%     m_L - mass of the latch
%     coeff_fric - coefficient of friction between the latch and load mass
%     (optional)
%     v_0 - initial velocity of the latch (optional)
% min # arguments = 2

function latch = rounded_latch(R, m_L, varargin)
    % optional parameters
    varargin_param_names = {'latch.coeff_fric', 'latch.v_0'};
    varargin_default_values = {0,0};
    
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
    latch.max_width = R;
    latch.mass = m_L;
    
    yL = @(x) latch.max_width*(1-sqrt(1-x^2/latch.max_width^2)); % function for the shape of the latch
    syms x;
    yL_prime = diff(yL(x));
    yL_doubleprime = diff(yL(x),2);
    latch.y_L = {yL, matlabFunction(yL_prime), matlabFunction(yL_doubleprime)}; % stores yL and its derivatives
end 