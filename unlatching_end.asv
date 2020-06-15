function [value,isterminal,direction] = unlatching_end(t,x,m_eff,m_L,F_s,F_out,y0,y_L,coeff_fric)
%End condition for unlatching
%First we find the second time derivative of x

yL=y_L{1}(x(1))+y0;
yL_prime = y_L{2}(x(1));
yL_doubleprime = y_L{3}(x(1));
y=[yL yL_prime*x(2)];
scaling = sign(x(2));
denom=m_L+m_eff*(yL_prime^2)+scaling*coeff_fric*(m_L*yL_prime-m_eff*yL_prime);
num=F_out(t,x)+yL_prime*F_s(t,y)+scaling*coeff_fric*(F_out(t,x)*yL_prime-F_s(t,y))-(x(2)^2)*(m_eff*yL_prime*yL_doubleprime-m_eff*scaling*coeff_fric*yL_doubleprime);
xL_doubledot=num/denom;



value = F_s(t,y)-m_eff*(yL_doubleprime*x(2)^2+xL_doubledot*yL_prime);
if (imag(value))
    warning('Complex value of Normal Force. ODE step size might be too large near the end of the latch');
end

isterminal=1;
direction=0;
end

