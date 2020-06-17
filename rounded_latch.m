function latch = rounded_latch(R, m_L, coeff_fric, v_0)
if (nargin==2) %assume only input (R, m_L), assume v_0,coeff_fric=0
    latch.v_0=0;
    latch.coeff_fric = 0;
elseif (nargin==3)%assume only (R,m_L,coeff_fric), assume v_0 = 0
    latch.coeff_fric = coeff_fric;
    latch.v_0=0;
elseif (nargin==4) %if all arguments are inputted
    latch.coeff_fric = coeff_fric;
    latch.v_0=v_0;
else 
    error('Rounded latch requires 2-4 arguments');
end

%initializing latch parameters within struct that are unchanged based on
%nargin
latch.max_width=R;
latch.mass=m_L;
yL=@(x) latch.max_width*(1-sqrt(1-x^2/latch.max_width^2));
syms x;
yL_prime = diff(yL(x));
yL_doubleprime = diff(yL(x),2);
latch.y_L = {yL, matlabFunction(yL_prime), matlabFunction(yL_doubleprime)};

end 