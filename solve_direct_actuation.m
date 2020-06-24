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
launch_opts=odeset('Events',@(t,y) direct_actuation_end(t,y,motor),'RelTol',1E-7,'AbsTol',1E-10);
ode=@(t,y) direct_actuation_ode(t,y,load,motor);
t_guess=sqrt((2*motor.range*load.mass)/motor.max_force);
tspan=linspace(0,t_guess*100,10000);
y0=[0,0];
sol=ode45(ode,tspan,y0,launch_opts);

T=linspace(0,sol.xe,1000);
Y=deval(sol,T);
Y=Y';
for i = 1:size(T)
    fMotor(i) = motor.Force(T(i), [Y(i,:)]);
end
fMotor = fMotor';

sol=[T,Y,fMotor];