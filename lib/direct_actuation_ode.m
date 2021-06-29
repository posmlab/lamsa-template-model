function dy=direct_actuation_ode(t,y,load,motor)
%ODE for direct actuation: a=F/m
dy=zeros(2,1);
dy(1)=y(2);

y1 = motor.rest_length;
y2 = motor.rest_length - y(1);

dy(2)=motor.Force(t,y)/load.mass([y1,y2]);

end

