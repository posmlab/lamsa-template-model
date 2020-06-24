function dy=direct_actuation_ode(t,y,load,motor)
%ODE for direct actuation: a=F/m
dy=zeros(2,1);
dy(1)=y(2);
dy(2)=motor.Force(t,y)/load.mass .*(motor.Force(t,y)>=0);
end
