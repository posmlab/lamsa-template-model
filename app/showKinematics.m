function showKinematics(axesHandle,evnt, xAxisVariable, yAxisVariable, ...
                                         xLogspaceFlag, yLogspaceFlag, ...
                                         app)

% ptH = getappdata(axesHandle,'CurrentPoint')
% pt = get(gcf,'CurrentPoint') %Getting click position

% properties(axesHandle)
x = evnt.IntersectionPoint(1);
y = evnt.IntersectionPoint(2);
varnames = {'load_mass_mass',...
                'latch_mass','latch_coeff_fric','latch_radius','latch_v_0','min_latching_dist','max_latching_dist',...
                'linear_spring_k','linear_spring_mass','linear_spring_Fmax',...
                'exp_spring_k','exp_spring_char_len','exp_spring_Fmax','exp_spring_mass',...
                'lee_spring_E','lee_spring_A','lee_spring_L','lee_spring_rho','lee_spring_sigma_f',...
                'lm_linear_motor_Fmax','lm_linear_motor_Vmax','lm_linear_motor_range_of_motion','lm_linear_motor_voltage_frac'...
                'lm_hill_motor_Fmax','lm_hill_motor_Vmax','lm_hill_motor_muscle_length','lm_hill_motor_rate_of_activation','lm_hill_motor_L_i','lm_hill_motor_a_L','lm_hill_motor_b_L','lm_hill_motor_s',...
                'um_linear_motor_Fmax','um_linear_motor_Vmax','um_linear_motor_range_of_motion','um_linear_motor_voltage_frac'...
                'um_hill_motor_Fmax','um_hill_motor_Vmax','um_hill_motor_muscle_length','um_hill_motor_rate_of_activation','um_hill_motor_L_i','um_hill_motor_a_L','um_hill_motor_b_L','um_hill_motor_s'};
            latexlabels = {'load mass',...
                'latch mass','latch $\mu$','latch radius','latch $v_0$','min latching distance','max latching distance',...
                'linear spring k','linear spring mass','linear spring Fmax',...
                'exponential spring $k_0$','exponential spring characteristic length','exponential spring Fmax','exponential spring mass',...
                'spring E','spring A','spring L','spring $\rho$','spring $\sigma_f$',...
                'loading motor (linear) Fmax','loading motor (linear) motor Vmax','loading motor (linear) range of motion','loading motor (linear) voltage fraction'...
                'loading motor (Hill) Fmax','loading motor (Hill) Vmax','loading motor (Hill) muscle length','loading motor (Hill) rate of activation','loading motor (Hill) motor $L_i$','loading motor (Hill) $a_L$','loading motor (Hill) $b_L$','loading motor (Hill) s',...
                'unlatching motor (linear) Fmax','unlatching motor (linear) motor Vmax','unlatching motor (linear) range of motion','unlatching motor (linear) voltage fraction'...
                'unlatching motor (Hill) Fmax','unlatching motor (Hill) Vmax','unlatching motor (Hill) muscle length','unlatching motor (Hill) rate of activation','unlatching motor (Hill) motor $L_i$','unlatching motor (Hill) $a_L$','unlatching motor (Hill) $b_L$','unlatching motor (Hill) s'};
            nonlatexlabels = {'load mass',...
                'latch mass','latch mu','latch radius','latch v_0','min latching distance','max latching distance',...
                'linear spring k','linear spring mass','linear spring Fmax',...
                'exponential spring k_0','exponential spring characteristic length','exponential spring Fmax','exponential spring mass',...
                'spring E','spring A','spring L','spring rho','spring sigma_f',...
                'loading motor (linear) Fmax','loading motor (linear) motor Vmax','loading motor (linear) range of motion', 'loading motor (linear) voltage fraction'...
                'loading motor (Hill) Fmax','loading motor (Hill) Vmax','loading motor (Hill) muscle length','loading motor (Hill) rate of activation','loading motor (Hill) motor L_i','loading motor (Hill) a_L','loading motor (Hill) b_L','loading motor (Hill) s',...
                'unlatching motor (linear) Fmax','unlatching motor (linear) motor Vmax','unlatching motor (linear) range of motion','unlatching motor (linear) voltage fraction'...
                'unlatching motor (Hill) Fmax','unlatching motor (Hill) Vmax','unlatching motor (Hill) muscle length','unlatching motor (Hill) rate of activation','unlatching motor (Hill) motor $L_i$','unlatching motor (Hill) $a_L$','unlatching motor (Hill) $b_L$','unlatching motor (Hill) s'};
            
axis_labels_dict = containers.Map(varnames,latexlabels);
dropdown_items_dict = containers.Map(varnames,nonlatexlabels);
dropdown_items_opposite_dict = containers.Map(nonlatexlabels,varnames);

if (xLogspaceFlag)
    eval(['app.' xAxisVariable '.Value = ' num2str(10^(x)) ';']);
else
    eval(['app.' xAxisVariable '.Value = ' num2str(x) ';']);
end
if (yLogspaceFlag)
    eval(['app.' yAxisVariable '.Value = ' num2str(10^(y)) ';']);
else
    eval(['app.' yAxisVariable '.Value = ' num2str(y) ';']);
end
%% initializing LaMSA component structs

% load mass struct initialization
load = load_mass(app.load_mass_mass.Value,app.load_m_rod.Value,app.load_EMA.Value);

% latch struct initialization
latch = rounded_latch(app.latch_radius.Value,app.latch_mass.Value,app.latch_coeff_fric.Value, app.latch_v_0.Value, app.min_latching_dist.Value, app.max_latching_dist.Value, app.runway_length.Value);

% spring struct initialization
if (app.spring.SelectedTab == app.linear_spring)
    spring = linear_spring(app.linear_spring_k.Value,app.linear_spring_mass.Value,app.linear_spring_Fmax.Value);
elseif (app.spring.SelectedTab == app.exponential_spring)
    spring = exponential_spring(app.exp_spring_k.Value,app.exp_spring_char_len.Value,app.exp_spring_mass.Value,app.exp_spring_Fmax.Value);
elseif (app.spring.SelectedTab == app.linear_elastic_extensional_spring)
    spring = linear_elastic_extensional_spring(app.lee_spring_E.Value,app.lee_spring_A.Value,app.lee_spring_L.Value,app.lee_spring_rho.Value,app.lee_spring_sigma_f.Value);
end

% loading motor struct initialization
if (app.loading_motor.SelectedTab == app.lm_linear_motor)
    loading_motor = linear_motor(app.lm_linear_motor_Fmax.Value,app.lm_linear_motor_Vmax.Value,app.lm_linear_motor_range_of_motion.Value,app.lm_linear_motor_voltage_frac.Value);
elseif (app.loading_motor.SelectedTab == app.lm_hill_muscle_motor)
    loading_motor = hill_muscle_motor(app.lm_hill_motor_muscle_length.Value,app.lm_hill_motor_Fmax.Value,app.lm_hill_motor_Vmax.Value,app.lm_hill_motor_rate_of_activation.Value,app.lm_hill_motor_L_i.Value,app.lm_hill_motor_a_L.Value,app.lm_hill_motor_b_L.Value,app.lm_hill_motor_s.Value);
end

% unlatching motor struct initialization
if (app.unlatching_motor.SelectedTab == app.um_linear_motor)
    unlatching_motor = linear_motor(app.um_linear_motor_Fmax.Value,app.um_linear_motor_Vmax.Value,app.um_linear_motor_range_of_motion.Value,app.um_linear_motor_voltage_frac.Value);
elseif (app.unlatching_motor.SelectedTab == app.um_hill_muscle_motor)
    unlatching_motor = hill_muscle_motor(app.um_hill_motor_muscle_length.Value,app.um_hill_motor_Fmax.Value,app.um_hill_motor_Vmax.Value,app.um_hill_motor_rate_of_activation.Value,app.um_hill_motor_L_i.Value,app.um_hill_motor_a_L.Value,app.um_hill_motor_b_L.Value,app.um_hill_motor_s.Value);
end

output_directory = create_output_directory();

% calling solve model
[sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);

plotNames = {'Time vs. y coordinate of load', 'Time vs. y velocity of load', 'Time vs. y-forces on load',...
             'Time vs. x coordinate of latch','Time vs. x velocity of latch','Time vs. x-forces on latch'};

[indexs, booleanValue] = listdlg('ListString', plotNames, ...
                             'Name','Plot Selection',...
                             'PromptString','Which kinematic plot(s) would you like to see?',... 
                             'ListSize', [800, 500]);

if (booleanValue)
figure;

if ((ismember(1,indexs) || ismember(2,indexs) || ismember(3,indexs)) && ...
    (ismember(4,indexs) || ismember(5,indexs) || ismember(6,indexs)))
    numColumns = 2;  
else
    numColumns = 1;
end

numYPlots = ismember(1, indexs) + ismember(2, indexs) + ismember(3, indexs);
numXPlots = ismember(4, indexs) + ismember(5, indexs) + ismember(6, indexs);
numRows = max(numXPlots, numYPlots);


% I'm sorry. This is awful. 
% This is textbook how-NOT-to do if statements.
% but this is what seems readable and how my brain structured it in my
% head.

if (numColumns == 1)
    if (ismember(1,indexs) || ismember(2, indexs) || ismember(3, indexs)) % if plotting y column
        if (ismember(1, indexs))
            subplot(numRows, numColumns, 1)
            plot(sol(:,1),sol(:,2))
            xlabel("time",'Interpreter','latex')
            ylabel("y coordinate of load",'Interpreter','latex')
            title("y coordinate of load vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(2, indexs))
            if (ismember(1,indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            plot(sol(:,1),sol(:,3))
            xlabel("time",'Interpreter','latex')
            ylabel("y velocity of load",'Interpreter','latex')
            title("y velocity of load vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(3, indexs))
            if (ismember(1, indexs) && ismember(2, indexs))
                subplot(numRows, numColumns, 3)
            elseif (ismember(1, indexs) || ismember(2, indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            y_forces = sol(:,10)+sol(:,9)-sol(:,7);
            plot(sol(:,1),y_forces)
            xlabel("time",'Interpreter','latex')
            ylabel("y forces on load",'Interpreter','latex')
            title("y forces on load vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        
    else % else: plotting x column
        if (ismember(4, indexs))
            subplot(numRows, numColumns, 1)
            plot(sol(:,1),sol(:,4))
            xlabel("time",'Interpreter','latex')
            ylabel("x coordinate of latch",'Interpreter','latex')
            title("x coordinate of latch vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(5,indexs))
            if (ismember(4, indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            plot(sol(:,1),sol(:,5))
            xlabel("time",'Interpreter','latex')
            ylabel("x velocity of latch",'Interpreter','latex')
            title("x velocity of latch vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(6, indexs))
            if (ismember(4, indexs) && ismember(5, indexs))
                subplot(numRows, numColumns, 3)
            elseif (ismember(4, indexs) || ismember(5, indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            x_forces = sol(:,11)+sol(:,8)+sol(:,6);
            plot(sol(:,1),x_forces)
            xlabel("time",'Interpreter','latex')
            ylabel("x forces on latch",'Interpreter','latex')
            title("x forces on latch vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
    end
end

if (numColumns == 2)
    if (ismember(1, indexs)) % t vs y
        subplot(numRows, numColumns,1)
        plot(sol(:,1),sol(:,2))
        xlabel("time",'Interpreter','latex')
        ylabel("y coordinate of load",'Interpreter','latex')
        title("y coordinate of load vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
    if (ismember(4, indexs)) % t vs x
        subplot(numRows, numColumns,2)
        plot(sol(:,1),sol(:,4))
        xlabel("time",'Interpreter','latex')
        ylabel("x coordinate of latch",'Interpreter','latex')
        title("x coordinate of latch vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
    if (ismember(2, indexs)) % t vs ydot
        if (ismember(1, indexs))
            subplot(numRows, numColumns, 3)
        else
            subplot(numRows, numColumns, 1)
        end
        plot(sol(:,1),sol(:,3))
        xlabel("time",'Interpreter','latex')
        ylabel("y velocity of load",'Interpreter','latex')
        title("y velocity of load vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
    if (ismember(5, indexs)) % t vs xdot
        if (ismember(4, indexs))
            subplot(numRows, numColumns, 4)
        else
            subplot(numRows, numColumns, 2)
        end
        plot(sol(:,1),sol(:,5))
        xlabel("time",'Interpreter','latex')
        ylabel("x velocity of latch",'Interpreter','latex')
        title("x velocity of latch vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end

    if (ismember(3, indexs)) % t vs yforces
        if ((~ismember(1, indexs)) && (~ismember(2, indexs)))
            subplot(numRows, numColumns, 1)
        elseif xor(ismember(1, indexs),ismember(2, indexs))
            subplot(numRows, numColumns, 3)
        else
            subplot(numRows, numColumns, 5)
        end
        y_forces = sol(:,10)+sol(:,9)-sol(:,7);
        plot(sol(:,1), y_forces)
        xlabel("time",'Interpreter','latex')
        ylabel("y forces on load",'Interpreter','latex')
        title("y forces on load vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end

    if (ismember(6, indexs)) % t vs xforces
        if ((~ismember(4, indexs)) && (~ismember(5, indexs)))
            subplot(numRows, numColumns, 2)
        elseif xor(ismember(4, indexs),ismember(5, indexs))
            subplot(numRows, numColumns, 4)
        else
            subplot(numRows, numColumns, 6)
        end
        x_forces = sol(:,11)+sol(:,8)+sol(:,6);
        plot(sol(:,1), x_forces)
        xlabel("time",'Interpreter','latex')
        ylabel("x forces on latch",'Interpreter','latex')
        title("x forces on latch vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
end

end

end