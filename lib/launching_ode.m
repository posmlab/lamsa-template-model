function dx = launching_ode(t,x,m_eff,spring)
%ODE for the launching phase
dx(1)=x(2,1);
dx(2)=spring.Force(t,x)/m_eff;
dx=dx';
end

