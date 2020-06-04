function dx = launching_ode(t,x,Fs,m_eff)
%ODE for the launching phase
dx(1)=x(2,1);
dx(2)=Fs(t,x)/m_eff;
dx=dx';
end

