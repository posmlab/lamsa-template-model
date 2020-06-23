function [value,isterminal,direction] = unlatching_end(t,x,m_eff,y0,latch,spring,unlatching_motor)
%End condition for unlatching
%First we find the second time derivative of x

if (x(1) > latch.max_width)
   value = 0;
   isterminal = 1;
   direction = 0;
   return
end

yL=latch.y_L{1}(x(1))+y0;
yL_prime = latch.y_L{2}(x(1));
yL_doubleprime = latch.y_L{3}(x(1));
y=[yL yL_prime*x(2)];
frictionDirection = sign(x(2));
denom=latch.mass+m_eff*(yL_prime^2)+frictionDirection*latch.coeff_fric*(latch.mass*yL_prime-m_eff*yL_prime);
num=unlatching_motor.Force(t,x)+yL_prime*spring.Force(t,y)+frictionDirection*latch.coeff_fric*(unlatching_motor.Force(t,x)*yL_prime-spring.Force(t,y))-(x(2)^2)*(m_eff*yL_prime*yL_doubleprime-m_eff*frictionDirection*latch.coeff_fric*yL_doubleprime);
xL_doubledot=num/denom;



value = spring.Force(t,y)-(m_eff*(yL_doubleprime*(x(2)^2)+(xL_doubledot*yL_prime)));

if (imag(value))
    warning('Complex value of Normal Force. ODE step size might be too large near the end of the latch');
end

if (( x(2) < 1E-3) && (xL_doubledot < 1E-3))
    error('Latch gets stuck!')
end

isterminal=1;
direction=0;
end

