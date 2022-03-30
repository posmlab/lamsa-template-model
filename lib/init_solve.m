%%init_solve
% initializes the structs using inputs from the GUI, and calls solve_lamsa
% inputs:
%   lm, um, ld, lt, sp are function handles that take in x0
% outputs:
%   value of a metric from solve_lamsa

function output = init_solve(lm, um, ld, lt, sp, x0, metric, DA_flag)
    if (nargin<7)
        metric = 'Pmax';
    end
    if (nargin<8)
        DA_flag = false;
    end
    
    loading_motor = lm(x0);
    unlatching_motor = um(x0);
    load = ld(x0);
    latch = lt(x0);
    spring = sp(x0);
    
    if (~DA_flag)
        [sol, trans_times] = solve_lamsa(loading_motor, unlatching_motor, load, latch, spring);
    else
        [sol, trans_times] = solve_direct_actuation(loading_motor,load);
    end
    
    cont = get_metrics(sol, trans_times, load, metric);
    output = cont(metric);
end