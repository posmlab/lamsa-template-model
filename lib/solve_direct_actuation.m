function [sol, transition_times] = solve_direct_actuation(motor,load)
%Solve differential equation launching phase of direction of direct
%actuated motion usinf F=ma
%motor struct
    %motor.Force
    %motor.max_force
    %motor.r
%load_mass struct
    %load.mass


%% Ballistic phase:F motor only 
launch_opts=odeset('Events',@(t,y) direct_actuation_end(t,y,motor),'RelTol',1E-7,'AbsTol',1E-10);
ode=@(t,y) direct_actuation_ode(t,y,load,motor);

t_guess_v=(motor.velocity*load.mass)/motor.max_force;
t_guess_pos=sqrt((2*motor.range*load.mass)/motor.max_force);

t_guess=max(t_guess_v,t_guess_pos);
tspan=linspace(0,t_guess,10000);
y0=[0,0];
[t,y]=ode45(ode,tspan,y0,launch_opts);

% run ode45 until the projectile launches
while (t(end) == tspan(end))
    t_guess = 10 *t_guess;
    tspan = linspace(0, t_guess,10000);
    [t,y]=ode45(ode,tspan,y0,launch_opts);
end

T=[t];
Y=[y];

for i = 1:size(T)
    fMotor(i) = motor.Force(T(i), [Y(i,:)]);
end
fMotor = fMotor';

transition_times = [T(1), T(end)];
sol=[T,Y,fMotor];
end

function dy=direct_actuation_ode(t,y,load,motor)
%ODE for direct actuation: a=F/m
dy=zeros(2,1);
dy(1)=y(2);

Force = F_eff(load, motor, t, y);
mass = m_eff(load, y);

dy(2)=Force/mass;

end

function [value,isterminal,direction]=direct_actuation_end(t,y,motor)
%End condition for loading
value=motor.Force(t,y);
isterminal=1;
direction=0;
end

function mass = m_eff(load, y)

mass = load.mass;

if isa(load, 'RotatingMass')
    L1 = load.L1;
    theta_0 = load.theta_0;
    theta = asin((-y(1)+L1*sin(theta_0))/L1);
    mass = mass/cos(theta) + 1E2*(abs(theta) > pi/2);
end

if isa(load, 'RotatingMassSE')
    L1 = load.lengths(1);
    moI = load.mass;
    mass = moI/L1;
end

end

function Force = F_eff(load, motor, t, y)

Force = motor.Force(t, y);

if isa(load, 'RotatingMass')
    
    y1 = motor.rest_length;
    y2 = motor.rest_length - y(1);
    theta_0 = load.theta_0;
    L1 = load.L1;
    m = load.mass;
    
    theta = asin((-y(1)+L1*sin(theta_0))/L1);
    alpha = acos((y2^2-y1^2+2*y1*L1*sin(theta_0))/(2*y2*L1));
    
    Force = real( Force*sin(alpha) - abs(m*sin(theta)*y(2)^2/(cos(theta)^3))) * (abs(theta) < pi/2);
end


if isa(load, 'RotatingMassSE')
    l0 = motor.rest_length;
    L1 = load.lengths(1);
    theta0 = load.theta_0;
    beta = sqrt(2*L1^2*(1-cos(y(1)-theta0)) + l0^2 - 2*l0*L1*(sin(y(1))- sin(theta0)));
    gamma = (L1^2*sin(y(1)-theta0) - l0*L1*cos(y(1)))/beta;
    y2 = l0 - beta;
    y2dot = -gamma*y(2);
    alpha = asin(L1*(cos(theta0) -  cos(y(2)))/(l0 - y2)); %Angle the spring makes with the vertical

    Force = motor.Force(t, [y2 y2dot])*sin(pi/2 - alpha - y(1));
end

end