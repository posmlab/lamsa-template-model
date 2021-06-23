function x = exit_lamsa_zone()

    x0 = [0.01 0 1 0.005 0.003 0 0 0 1000 0 2000 .001 .002 20 .01 30 1 200 .01 2.08 -2.89 -.75 .25 1 .005 1];
    dx = .001*x0';
    metric = 'amax';
    
    x = x0';
    while true
        load_f = @(x0)LoadMass(x0(1),x0(2),x0(3));
        latch_f = @(x0)RoundedLatch(x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10));
        spring_f = @(x0)ExponentialSpring(x0(11),x0(12),x0(13),x0(14));
        lm_f = @(x0)HillMuscleMotor(x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22));
        um_f = @(x0)LinearMotor(x0(23),x0(24),x0(25),x0(26));

        wrapper_func = @(x0)init_solve(lm_f,um_f,load_f,latch_f,spring_f,x0,metric);
        grad = relative_gradient(wrapper_func,x0,0.01)*1.0e-3;
        x = x - dx.*grad;
        
        load = load_f(x);
        latch = latch_f(x);
        spring = spring_f(x);
        lm = lm_f(x);
        um = um_f(x);
        
        [sol,transition_times] = solve_model(lm,um,load,latch,spring);
        [direct_sol, direct_transition_times] = solve_direct_actuation(lm,load);
        met_dict=get_metrics(sol,transition_times,load,metric);
        direct_met_dict=get_metrics(direct_sol,direct_transition_times,load,metric);
        ratio = met_dict('amax')/direct_met_dict('amax')
        if ratio <= 1
            return
        end
    end

end