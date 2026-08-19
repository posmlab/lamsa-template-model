function metrics = get_metrics(sol, transition_times, load, met_names)
    % inputs:
    %   solution matrix with columns [time, position, velocity]
    %   and the effective mass (EMA)
    % outputs: 
    %   metrics specified in met_names

    m = load.mass;
    L1 = load.L1
    L2 = L1 ./ load.EMA([0 0]);
    theta_0 = load.theta_0;
    theta_f = load.theta_f;
    y0 = sol(1,2)

    %%% convert to angular space
    sol(:,2) = asin((-sol(:,2).+y0.+L1.*sin(theta_0))./L1);
    sol(:,3) = sol(:,3)./L1;

    %%% determine metrics and calculate each one
    metrics=containers.Map(met_names,zeros(length(met_names),1),'UniformValues',false);

    if isKey(metrics, 't_unlatch')
        metrics('t_unlatch') = transition_times(1);
    end
    if isKey(metircs, 'theta_unlatch')
        index = find(sol(:,1)<=transition_times(1), 1, 'last');
        metrics('theta_unlatch') = sol(index, 2);
    end
    if isKey(metrics, 'omega_unlatch')
        index = find(sol(:,1)<=transition_times(1), 1, 'last');
        metrics('omega_unlatch') = sol(index, 3);
    end
    if isKey(metrics, 'v_unlatch')
        index = find(sol(:,1)<=transition_times(1), 1, 'last');
        metrics('v_unlatch') = sol(index, 3) .* L2;
    end
    if isKey(metrics, 'KE_unlatch')
        index = find(sol(:,1)<=transition_times(1), 1, 'last');
        metrics('v_unlatch') = 0.5 .* m .* L1.^2 .* sol(index, 3).^2;
    end

    if isKey(metrics, 't_takeoff')
        metrics('t_takeoff') = sol(end, 1);
    end
    if isKey(metircs, 'theta_takeoff')
        metrics('theta_takeoff') = sol(end, 2);
    end
    if isKey(metrics, 'omega_takeoff')
        metrics('omega_takeoff') = sol(end, 3);
    end
    if isKey(metrics, 'v_takeoff')
        metrics('v_takeoff') = sol(end, 3) .* L2;
    end
    if isKey(metrics, 'KE_takeoff')
        metrics('v_takeoff') = 0.5 .* m .* L1.^2 .* sol(end, 3).^2;
    end

    if isKey(metrics, 't_close')
        index = find(sol(:,2)<=theta_f, 1, 'first')
        if isempty(index)
            metrics('t_close') = sol(end, 1) + (sol(end, 2) - theta_f) / sol(end, 3);
        else
            metrics('t_close') = sol(index, 1);
        end
    end
