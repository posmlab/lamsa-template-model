function dx = unlatching_ode(t,x,m_eff,m_L,F_s,F_out,y0,y_L)
% ODE for the unlatching phase
y=y_L{1}(x(1))+y0;

yL_prime = y_L{2}(x(1));
yL_doubleprime = y_L{3}(x(1));
Y=[y,yL_prime*x(2)];
denom=m_L+m_eff*yL_prime^2;
num=F_out(t,x)+yL_prime*F_s(t,Y)-m_eff*yL_prime*yL_doubleprime*x(2)^2;
dx=[x(2);num/denom];
end