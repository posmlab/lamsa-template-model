%% linear_spring struct
% arguments:
%     k - spring constant
%     m_s - mass of the spring
%     F_spring_max - maximum amount of force the spring can exert
% min # arguments = 2

function spring = linear_spring(k,m_s,F_spring_max)
    if (nargin == 2) % only k and m_s passed in
        spring.Force = @(t,x)-k*x(1);
    elseif (nargin == 3) % k_0, m_s, and F_spring_max passed in
        spring.Force = @(t,x)-k*x(1).*(abs(k*x(1))<F_spring_max);  
    else
        error('Linear spring requires 2 or 3 inputs. See comments in linear_spring.m for more information.');
    end
    spring.mass = m_s;
end

