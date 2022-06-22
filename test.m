metrics = {'tto', 'vto', 'Pmax', 'ymax', 'tL', 'KEmax', 'yunlatch', 'amax'};
loading_motor =LinearMotor(5, 5, .1);
load = RotatingMassSE(0, 0.01, 0.001, 0.001, 0.01, -pi/4);
unlatching_motor = LinearMotor(10, 10, .005);
latch = ParabolicLatch(150, .0005, .0003, 0, 0);
spring = ExponentialSpring(200, 1);

[t,y] = solve_lamsa_se(loading_motor,unlatching_motor,load,latch,spring);