function spring = linear_spring(k,F_spring_max,t,x)

if (nargin==1) % only k passed in, assume F_spring_max = inf, return function handle
    spring.Force = @(t,x)-k*x(1);
elseif (nargin==2) % if only (k, F_spring_max) passsed in, then return function handle
    spring.Force = @(t,x)-k*x(1).*(abs(k*x(1))<F_spring_max);
elseif (nargin==3) % only (k,t,x) passed in, assume F_spring_max = inf, return force value
    spring.Force = -k*x(1);
elseif (nargin==4) % else return value of force for a specicifc  
    spring.Force = -k*x(1).*(abs(k*x(1))<F_spring_max);  
else
    error('Linear spring needs between one to four arguments');
end
spring.Time_independent = true;
end % end linear_spring

