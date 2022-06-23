metrics = {'tto', 'vto', 'Pmax', 'ymax', 'tL', 'KEmax', 'yunlatch', 'amax'};
loading_motor =LinearMotor(20, 20, 0.25);
load = RotatingMassSE(0, 0.01, 0.2, 0.2, 1, -pi/4);
unlatching_motor = LinearMotor(1, 1, .05);
latch = ParabolicLatch(150, .05, .03, 0,1000);
spring = LinearSpring(2, 0.1, Inf, 0.5);

tspan = [0, 1];

[t,y] = solve_lamsa_se(tspan, loading_motor,unlatching_motor,load,latch,spring);

plot(t,real(y(:,2)))