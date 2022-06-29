metrics = {'tto', 'vto', 'Pmax', 'ymax', 'tL', 'KEmax', 'yunlatch', 'amax'};
loading_motor = HillMuscleMotor(.03, 10, 5, 100);
load = RotatingMassSE(1e-2, 3e-2 ,1e-3, 1e-3, 5e-2, -pi/4);
load2 = RotatingMass(0, 3e-2, 1e-3, 5e-2, -pi/4);
unlatching_motor = LinearMotor(1, 1, .05);
latch = ParabolicLatch(150, 1e-4, 3, 0, 1000);
spring = LinearSpring(1e8, 0.05, Inf, 5e-3);

tspan = [0, 1];

[sol, t_times] = solve_lamsa_se(tspan, loading_motor,unlatching_motor,load,latch,spring);
[sol2, t_times2] = solve_direct_actuation(loading_motor, load2);

figure();

plot(sol(:,1), sol(:,2))
hold on;
xline(t_times(1))
title("theta")
%plot(sol2(:,1), sol2(:,2))
hold off;

figure();

plot(sol(:,1), sol(:,3))
title("theta dot")
hold on;
%plot(sol2(:,1), sol2(:,3))
hold off;

figure();
 
plot(sol(:,1), sol(:,10))
title("F perp")

figure();
 
plot(sol(:,1), sol(:,12))
title("y1")

figure();
 
plot(sol(:,1), sol(:,13))
title("y1dot")