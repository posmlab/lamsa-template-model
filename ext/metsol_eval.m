%% solve model, evaluate metric
% inputs:
%       loading motor
%       unlatching motor
%       load
%       latch
%       spring
%       metric

function output = metsol_eval(lm, um, ld, lt, sp, metric)
    if nargin==5
        metric=='Pmax';
    end
    [sol, trans_times] = solve_model(lm, um, ld, lt, sp);
    cont = get_metrics(sol, trans_times, ld, metric);
    output = cont(metric);
end