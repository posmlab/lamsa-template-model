function spring = exponential_spring(k_0, characteristic_length, m_s, F_spring_max)
   spring.mass=m_s;
   if (nargin==3) % only k_0, m_s, characteristic_length passed in, assume F_spring_max = inf, return function handle
    spring.Force=@(t,x)-characteristic_length*k_0*(exp(x(1)/characteristic_length)-1);
   elseif (nargin==4) %return function handle
    spring.Force=@(t,x)-characteristic_length*k_0*(exp(x(1)/characteristic_length)-1).*(abs(characteristic_length*k_0*(exp(x(1)/characteristic_length)-1))<F_spring_max);
   else
    error('Exponential spring needs three or four arguments');
   end
end


   