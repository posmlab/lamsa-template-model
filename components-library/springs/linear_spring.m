function spring = linear_spring(k,m_s,F_spring_max,t,x)

if (nargin==2) % only k, m_s passed in, assume F_spring_max = inf, return function handle
    spring.Force = @(t,x)-k*x(1);
elseif (nargin==3) % if only (k, m_s,F_spring_max) passsed in, then return function handle
    spring.Force = @(t,x)-k*x(1).*(abs(k*x(1))<F_spring_max);
elseif (nargin==4) % only (k,m_s,t,x) passed in, assume F_spring_max = inf, return force value
    spring.Force = -k*x(1);
elseif (nargin==5) % else return value of force for a specicifc  
    spring.Force = -k*x(1).*(abs(k*x(1))<F_spring_max);  
else
    error('Linear spring needs between one to four arguments');
end
spring.mass = m_s;
end % end linear_spring

