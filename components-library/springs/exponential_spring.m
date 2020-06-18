%% exponential_spring struct
% arguments:
%     k_0 - initial spring constant
%     m_s - mass of the spring
%     characteristic_length - length of the spring at which force = 0
%     F_spring_max - maximum amount of force the spring can exert
% min # arguments = 3

function spring = exponential_spring(k_0, m_s, characteristic_length, F_spring_max)
   if (nargin == 3) % only k_0, m_s, and characteristic_length passed in
        spring.Force = @(t,x)-characteristic_length*k_0*(exp(x(1)/characteristic_length)-1);
   elseif (nargin == 4) % k_0, m_s, characteristic_length, and F_spring_max passed in
        spring.Force = @(t,x)-characteristic_length*k_0*(exp(x(1)/characteristic_length)-1).*(abs(characteristic_length*k_0*(exp(x(1)/characteristic_length)-1))<F_spring_max);
   else
        error('Exponential spring requires 3 or 4 inputs. See comments in exponential_spring.m for more information.');
   end
   spring.mass=m_s;
end


   