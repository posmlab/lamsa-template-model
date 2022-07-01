warning('off','all')

metrics = {'tto', 'vto', 'Pmax', 'ymax', 'tL', 'KEmax', 'yunlatch', 'amax'};
% "muscle length" "Fmax" "Vmax" "rate of activation" "initial length" "a_L" "b_L" "s"
loading_motor = HillMuscleMotor(.001, 20, 0.1, 200,.001, 2.08, -2.89, -0.75);
% "mass" "mass of lever arm" "L1" "L2" "L3" "theta initial"
load = RotatingMassSE(0, 3e-2 ,2e-2, 2e-2, 5e-2, -pi/4);
%"mass" "mass of lever arm" "L1" "L2" "theta initial"
load2 = RotatingMass(1e-2, 3e-2, 1e-3, 5e-2, -pi/4);
% "Fmax" "Vmax" "range of motion" "voltage fraction" "no braking" "muscle length"
unlatching_motor = LinearMotor(1e-4, 1e-4, .5, 1, .05);
% "coeff_parabola" "parabola_width" "mass" "μ" "v_0" "min_latching_dist" "max_latching_dist" "runway_length"
latch = ParabolicLatch(150, 1e-3, 1e-2, 0, 1e-3, 0, Inf, 0);
% "k" "m_s" "F_spring_max" "rest length"
spring = LinearSpring(2000, 5, Inf, 0.02);

tspan = [0, 1];

[sol, t_times] = solve_lamsa_se(tspan, loading_motor,unlatching_motor,load,latch,spring);
%[sol2, t_times2] = solve_direct_actuation(loading_motor, load2);

figure();

plot(sol(:,1), sol(:,2))
hold on;
xline(t_times(1))
title("theta")
hold off;

figure();

plot(sol(:,1), sol(:,3))
title("theta dot")
hold on;
%plot(sol2(:,1), sol2(:,3))
hold off;


figure();

plot(sol(:,1), sol(:,4))
hold on;
xline(t_times(1))
yline(latch.max_width)
title("s")
hold off;

figure();

plot(sol(:,1), sol(:,5))
hold on;
xline(t_times(1))
title("s dot")
hold off;

figure();
 
plot(sol(:,1), sol(:,10))
title("F perp")

figure();
 
plot(sol(:,1), sol(:,7))
title("Normal Torque on Load")

figure();
 
plot(sol(:,1), sol(:,12))
title("y1")

figure();
 
plot(sol(:,1), sol(:,13))
title("y1dot")