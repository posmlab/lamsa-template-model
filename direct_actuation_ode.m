function dx=direct_actuation_ode(t,x,load,motor)
%ODE for direct actuation: a=F/m
dx=zeros(2,1);
dx(1)=x(2,1);
dx(2)=motor.Force(t,x)/load.mass;
dx=dx';
end
