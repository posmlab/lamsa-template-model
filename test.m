metrics = {'tto', 'vto', 'Pmax', 'ymax', 'tL', 'KEmax', 'yunlatch', 'amax'};
loading_motor =LinearMotor(20, 20, 2e-2);
load = RotatingMassSE(0, 3e-2 ,1e-3, 1e-3, 5e-2, -pi/4);
unlatching_motor = LinearMotor(1, 1, .05);
latch = ParabolicLatch(150, 1e-4, 3e-3, 10);
spring = LinearSpring(2000, 0.1, Inf, 5e-3);
springms = LinearSpring(200, 0.1, Inf, 5e-3);

tspan = [0, 1];

[sol, t_times] = solve_lamsa_se(tspan, loading_motor,unlatching_motor,load,latch,spring);

tspan = [0, t_times(2)];

[solms, ~] = solve_lamsa_se_manySprings(tspan, loading_motor,unlatching_motor,load,latch,springms);

plot(sol(:,1), sol(:,2))
hold on;
xline(t_times(1))
plot(solms(:,1), solms(:,2))
title("theta")
hold off;

figure();

plot(sol(:,1), sol(:,3))
hold on;
plot(solms(:,1), solms(:,3))
title("theta dot")
hold off;

figure();
 
plot(sol(:,1), sol(:,10))
hold on;
plot(solms(:,1), solms(:,10))
title("F perp")
hold off;

figure();
 
plot(sol(:,1), sol(:,12))
hold on;
plot(solms(:,1), solms(:,12))
title("y1")
hold off;

figure();
 
plot(sol(:,1), sol(:,13))
hold on;
plot(solms(:,1), solms(:,13))
title("y1dot")
hold off;

