%% rounded_latch struct
% arguments:
%     R - raduis of the latch
%     m_L - mass of the latch
%     coeff_fric - coefficient of friction between the latch and load mass
%     v_0 - initial velocity of the latch
% min # arguments = 2

function latch = rounded_latch(R, m_L, coeff_fric, v_0)
    if (nargin == 2) % only R and m_L passed in
        latch.coeff_fric = 0; % assume no friction
        latch.v_0 = 0; % assume 0 initial velocity
    elseif (nargin == 3) % only R, m_L , and coeff_fric passed in
        latch.coeff_fric = coeff_fric;
        latch.v_0 = 0; % assume 0 initial velocity 
    elseif (nargin == 4) % R, m_L , coeff_fric, and v_0 passed in
        latch.coeff_fric = coeff_fric;
        latch.v_0 = v_0;
    else 
        error('Rounded latch requires 2, 3, or 4 arguments. See comments in rounded_latch.m for more information');
    end
    latch.max_width = R;
    latch.mass = m_L;
    
    yL = @(x) latch.max_width*(1-sqrt(1-x^2/latch.max_width^2)); % function for the shape of the latch
    syms x;
    yL_prime = diff(yL(x));
    yL_doubleprime = diff(yL(x),2);
    latch.y_L = {yL, matlabFunction(yL_prime), matlabFunction(yL_doubleprime)}; % stores yL and its derivatives
end 