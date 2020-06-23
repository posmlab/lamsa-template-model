function [sol] = solve_direct_actuation(motor,load)
%Solve differential equation launching phase of direction of direct
%actuated motion usinf F=ma
%motor struct
    %motor.Force
    %motor.max_force
%load_mass struct
    %load.mass
    
%is there a drag force/gravity

%% Ballistic phase:F motor only 
launch_opts=odeset('Events',@(t,x) direct_actutation_end(t,x,motor));
ode=@(t,x) direct_actuation_ode(t,x,load,motor);
tspan=linspace(0,10,1000);
x0=0;
[t,x]=ode45(ode,tspan,x0,launch_opts);

% Solve latch dynamics during Ballistic Phase

end