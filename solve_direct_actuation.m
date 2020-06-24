function sol = solve_direct_actuation(motor,load)
%Solve differential equation launching phase of direction of direct
%actuated motion usinf F=ma
%motor struct
    %motor.Force
    %motor.max_force
%load_mass struct
    %load.mass
    
%is there a drag force/gravity?


%% Ballistic phase:F motor only 
launch_opts=odeset('Events',@(t,y) direct_actuation_end(t,y,motor));
ode=@(t,y) direct_actuation_ode(t,y,load,motor);
tspan=linspace(0,10,1000);
y0=[0,0];
[t,y]=ode45(ode,tspan,y0,launch_opts);
T=[t];
Y=[y];
for i = 1:size(T)
    fMotor(i) = motor.Force(T(i), [Y(i,:)]);
end
fMotor = fMotor';
index=find(fMotor==0.0000,1,"first")
fMotor = fMotor(1:index,1);
T=T(1:index,1);
Y=Y(1:index,:);

sol=[T,Y,fMotor]