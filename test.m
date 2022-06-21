metrics = {'tto', 'vto', 'Pmax', 'ymax', 'tL', 'KEmax', 'yunlatch', 'amax'};
loading_motor = HillMuscleMotor(.01, 20, 5, 200);
load = RotatingMassSE(0);
unlatching_motor = LinearMotor(10, 10, .005);
latch = ParabolicLatch(150, .0005, .0003, 0, 0);
spring = ExponentialSpring(2000, 1e-4);
[t,y] = solve_lamsa_se(loading_motor,unlatching_motor,load,latch,spring);