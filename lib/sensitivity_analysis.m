%%sensitivity_analysis
% creates wrapper function to put into sensitive_axes
% finds ordered list of sensitive axes
% inputs:
%   function handles from GUI, requested metric, and labels (array of variables chosen)
% outputs:
%   matrix of most sensitive axes (directions)
%   list of variables in order of sensitivity

function [combos, var_list] = sensitivity_analysis(lm,um,ld,lt,sp,x0,metric,labels,DA_flag,first_order_flag)
    if (nargin<9)
        DA_flag = false;
    end
    metrics = {'tto','vto','Pmax','ymax','tL','KEmax','yunlatch','amax'};
    
    start = x0;
    
    %%Function
    wrapper_func = @(x0)init_solve(lm,um,ld,lt,sp,x0,metric,DA_flag);
    
    %%Data Dynamics
       
    [ax,pac] = sensitive_axes(wrapper_func,x0,0.05,0.05,0.05,8,first_order_flag);
    first_axis = ax(:,1).^2*100+100;
    for i=1:length(start)
        if first_axis(i) == 100
            first_axis(i) = pac(i);
        end
    end
    
    [~,R]=sort(first_axis,'descend');
    R=R';
    ordered_vars = strings(size(R));
    for i = 1:length(R)
        ordered_vars(i) = labels(R(i));
    end
    
    var_list = ordered_vars; % ordered from most to least sensitive
    combos = ax;
end
