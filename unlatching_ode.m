function dx = unlatching_ode(t,x,m_eff,y0,latch,spring,unlatching_motor)
%function dx = unlatching_ode(t,x,m_eff,m_L,F_s,F_out,y0,y_L, coeff_fric)
% ODE for the unlatching phase
y=latch.y_L{1}(x(1))+y0;
yL_prime = latch.y_L{2}(x(1));
yL_doubleprime = latch.y_L{3}(x(1));
Y=[y,yL_prime*x(2)];
scaling = sign(x(2));
denom=latch.mass+m_eff*(yL_prime^2)+scaling*latch.coeff_fric*(latch.mass*yL_prime-m_eff*yL_prime);
num=unlatching_motor.Force(t,x)+yL_prime*spring.Force(t,y)+scaling*latch.coeff_fric*(unlatching_motor.Force(t,x)*yL_prime-spring.Force(t,y))-(x(2)^2)*(m_eff*yL_prime*yL_doubleprime-m_eff*scaling*latch.coeff_fric*yL_doubleprime);
dx=[x(2);num/denom];
end