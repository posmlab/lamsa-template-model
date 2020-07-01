function metrics = get_metrics(sol,transition_times,load,met_names)
%Takes solution matrix (columns time, position, velocity)and the effective 
%    mass and returns the metrics specified in met_names
    m = load.mass;
    kinetic_energy=0.5*m*sol(:,3).^2;
    metrics=containers.Map(met_names,zeros(length(met_names),1),'UniformValues',false);
    if isKey(metrics,'vto')
        metrics('vto')=sol(end,3);
    end
    if isKey(metrics,'vmax')
        metrics('vmax')=max(sol(:,3));
    end
    if isKey(metrics,'tto')
        metrics('tto')=sol(end,1);
    end
    if isKey(metrics,'KEmax')
        metrics('KEmax')=max(0.5*m*sol(:,3).^2);
    end
    if isKey(metrics,'Pmax')
        metrics('Pmax')=max(gradient(kinetic_energy)./gradient(sol(:,1)));
    end
    if isKey(metrics,'ymax')
        metrics('ymax')=sol(1,2);
    end
    if isKey(metrics,'tL')
        metrics('tL')=transition_times(1);
    end
    if isKey(metrics,'yunlatch')
        index=find(sol(:,1)<=transition_times(1),1,"last");
        metrics('yunlatch')= sol(index, 2);
    end
    if isKey(metrics,'minumforce')
        metrics('minumforce')= min(sol(:,11));
    end
end