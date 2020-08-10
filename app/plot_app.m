classdef plot_app < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ParametersLabel                 matlab.ui.control.Label
        spring                          matlab.ui.container.TabGroup
        linear_spring                   matlab.ui.container.Tab
        kLabel                          matlab.ui.control.Label
        linear_spring_k                 matlab.ui.control.NumericEditField
        massEditFieldLabel_4            matlab.ui.control.Label
        linear_spring_mass              matlab.ui.control.NumericEditField
        FmaxLabel_3                     matlab.ui.control.Label
        linear_spring_Fmax              matlab.ui.control.NumericEditField
        linear_elastic_extensional_spring  matlab.ui.container.Tab
        ELabel                          matlab.ui.control.Label
        lee_spring_E                    matlab.ui.control.NumericEditField
        rhoLabel                        matlab.ui.control.Label
        lee_spring_rho                  matlab.ui.control.NumericEditField
        ALabel                          matlab.ui.control.Label
        lee_spring_A                    matlab.ui.control.NumericEditField
        LLabel                          matlab.ui.control.Label
        lee_spring_L                    matlab.ui.control.NumericEditField
        sigma_fLabel                    matlab.ui.control.Label
        lee_spring_sigma_f              matlab.ui.control.NumericEditField
        exponential_spring              matlab.ui.container.Tab
        k_0Label                        matlab.ui.control.Label
        exp_spring_k                    matlab.ui.control.NumericEditField
        characteristiclengthLabel       matlab.ui.control.Label
        exp_spring_char_len             matlab.ui.control.NumericEditField
        FmaxLabel_4                     matlab.ui.control.Label
        exp_spring_Fmax                 matlab.ui.control.NumericEditField
        massEditFieldLabel_5            matlab.ui.control.Label
        exp_spring_mass                 matlab.ui.control.NumericEditField
        loading_motor                   matlab.ui.container.TabGroup
        lm_linear_motor                 matlab.ui.container.Tab
        FmaxLabel_2                     matlab.ui.control.Label
        lm_linear_motor_Fmax            matlab.ui.control.NumericEditField
        VmaxLabel_2                     matlab.ui.control.Label
        lm_linear_motor_Vmax            matlab.ui.control.NumericEditField
        rangeofmotionLabel              matlab.ui.control.Label
        lm_linear_motor_range_of_motion  matlab.ui.control.NumericEditField
        voltagefracLabel                matlab.ui.control.Label
        lm_linear_motor_voltage_frac    matlab.ui.control.NumericEditField
        lm_hill_muscle_motor            matlab.ui.container.Tab
        FmaxLabel                       matlab.ui.control.Label
        lm_hill_motor_Fmax              matlab.ui.control.NumericEditField
        VmaxLabel                       matlab.ui.control.Label
        lm_hill_motor_Vmax              matlab.ui.control.NumericEditField
        musclelengthLabel               matlab.ui.control.Label
        lm_hill_motor_muscle_length     matlab.ui.control.NumericEditField
        rateofactivationLabel           matlab.ui.control.Label
        lm_hill_motor_rate_of_activation  matlab.ui.control.NumericEditField
        initiallengthLabel              matlab.ui.control.Label
        lm_hill_motor_L_i               matlab.ui.control.NumericEditField
        a_LLabel                        matlab.ui.control.Label
        lm_hill_motor_a_L               matlab.ui.control.NumericEditField
        b_LLabel                        matlab.ui.control.Label
        lm_hill_motor_b_L               matlab.ui.control.NumericEditField
        sLabel                          matlab.ui.control.Label
        lm_hill_motor_s                 matlab.ui.control.NumericEditField
        latch                           matlab.ui.container.TabGroup
        rounded_latch                   matlab.ui.container.Tab
        massLabel                       matlab.ui.control.Label
        latch_mass                      matlab.ui.control.NumericEditField
        Label                           matlab.ui.control.Label
        latch_coeff_fric                matlab.ui.control.NumericEditField
        radiusLabel                     matlab.ui.control.Label
        latch_radius                    matlab.ui.control.NumericEditField
        v_0Label                        matlab.ui.control.Label
        latch_v_0                       matlab.ui.control.NumericEditField
        minlatchingdistLabel            matlab.ui.control.Label
        min_latching_dist               matlab.ui.control.NumericEditField
        maxlatchingdistLabel            matlab.ui.control.Label
        max_latching_dist               matlab.ui.control.NumericEditField
        runwaylengthLabel               matlab.ui.control.Label
        runway_length                   matlab.ui.control.NumericEditField
        load                            matlab.ui.container.TabGroup
        load_mass                       matlab.ui.container.Tab
        massEditFieldLabel              matlab.ui.control.Label
        load_mass_mass                  matlab.ui.control.NumericEditField
        massofleverarmLabel             matlab.ui.control.Label
        load_m_rod                      matlab.ui.control.NumericEditField
        EMALabel                        matlab.ui.control.Label
        load_EMA                        matlab.ui.control.NumericEditField
        LoadingMotorLabel               matlab.ui.control.Label
        unlatching_motor                matlab.ui.container.TabGroup
        um_linear_motor                 matlab.ui.container.Tab
        FmaxLabel_5                     matlab.ui.control.Label
        um_linear_motor_Fmax            matlab.ui.control.NumericEditField
        VmaxLabel_3                     matlab.ui.control.Label
        um_linear_motor_Vmax            matlab.ui.control.NumericEditField
        rangeofmotionLabel_2            matlab.ui.control.Label
        um_linear_motor_range_of_motion  matlab.ui.control.NumericEditField
        voltagefracLabel_2              matlab.ui.control.Label
        um_linear_motor_voltage_frac    matlab.ui.control.NumericEditField
        um_hill_muscle_motor            matlab.ui.container.Tab
        FmaxLabel_6                     matlab.ui.control.Label
        um_hill_motor_Fmax              matlab.ui.control.NumericEditField
        VmaxLabel_4                     matlab.ui.control.Label
        um_hill_motor_Vmax              matlab.ui.control.NumericEditField
        musclelengthLabel_2             matlab.ui.control.Label
        um_hill_motor_muscle_length     matlab.ui.control.NumericEditField
        rateofactivationLabel_2         matlab.ui.control.Label
        um_hill_motor_rate_of_activation  matlab.ui.control.NumericEditField
        initiallengthLabel_2            matlab.ui.control.Label
        um_hill_motor_L_i               matlab.ui.control.NumericEditField
        a_LLabel_2                      matlab.ui.control.Label
        um_hill_motor_a_L               matlab.ui.control.NumericEditField
        b_LLabel_2                      matlab.ui.control.Label
        um_hill_motor_b_L               matlab.ui.control.NumericEditField
        sLabel_2                        matlab.ui.control.Label
        um_hill_motor_s                 matlab.ui.control.NumericEditField
        UnlatchingMotorLabel            matlab.ui.control.Label
        GraphingCornerLabel             matlab.ui.control.Label
        graphing_corner                 matlab.ui.container.TabGroup
        graphing_corner_kinematics      matlab.ui.container.Tab
        KinematicsOutputOptionsLabel    matlab.ui.control.Label
        forcedisp                       matlab.ui.control.CheckBox
        latchkinematicsCheckBox         matlab.ui.control.CheckBox
        loadkinematicsCheckBox          matlab.ui.control.CheckBox
        LoadingphaseLabel               matlab.ui.control.Label
        UnlatchingandlaunchingphaseLabel  matlab.ui.control.Label
        graphing_corner_one_D           matlab.ui.container.Tab
        OD_minunlatchingmotorforceCheckBox  matlab.ui.control.CheckBox
        xaxisLabel_2                    matlab.ui.control.Label
        OD_y_maxCheckBox                matlab.ui.control.CheckBox
        OD_y_unlatchCheckBox            matlab.ui.control.CheckBox
        OD_t_LCheckBox                  matlab.ui.control.CheckBox
        OD_v_toCheckBox                 matlab.ui.control.CheckBox
        OD_P_maxCheckBox                matlab.ui.control.CheckBox
        OD_t_toCheckBox                 matlab.ui.control.CheckBox
        OD_KE_maxCheckBox               matlab.ui.control.CheckBox
        OD_IV1DropDown                  matlab.ui.control.DropDown
        yaxisOutputOptionsLabel         matlab.ui.control.Label
        pixelsofresolutionLabel_2       matlab.ui.control.Label
        OD_n                            matlab.ui.control.NumericEditField
        xminEditFieldLabel_2            matlab.ui.control.Label
        OD_xmin                         matlab.ui.control.NumericEditField
        xmaxEditFieldLabel_2            matlab.ui.control.Label
        OD_xmax                         matlab.ui.control.NumericEditField
        OD_x_log_space                  matlab.ui.control.Switch
        graphing_corner_heatmap         matlab.ui.container.Tab
        minunlatchingmotorforceCheckBox  matlab.ui.control.CheckBox
        xaxisLabel                      matlab.ui.control.Label
        yaxisLabel                      matlab.ui.control.Label
        y_maxCheckBox                   matlab.ui.control.CheckBox
        y_unlatchCheckBox               matlab.ui.control.CheckBox
        t_LCheckBox                     matlab.ui.control.CheckBox
        v_toCheckBox                    matlab.ui.control.CheckBox
        P_maxCheckBox                   matlab.ui.control.CheckBox
        t_toCheckBox                    matlab.ui.control.CheckBox
        KE_maxCheckBox                  matlab.ui.control.CheckBox
        IV1DropDown                     matlab.ui.control.DropDown
        IV2DropDown                     matlab.ui.control.DropDown
        xminEditFieldLabel              matlab.ui.control.Label
        xmin                            matlab.ui.control.NumericEditField
        xmaxEditFieldLabel              matlab.ui.control.Label
        xmax                            matlab.ui.control.NumericEditField
        yminLabel                       matlab.ui.control.Label
        ymin                            matlab.ui.control.NumericEditField
        ymaxLabel                       matlab.ui.control.Label
        ymax                            matlab.ui.control.NumericEditField
        HeatmapOutputOptionsLabel       matlab.ui.control.Label
        pixelsofresolutionLabel         matlab.ui.control.Label
        n                               matlab.ui.control.NumericEditField
        x_log_space                     matlab.ui.control.Switch
        y_log_space                     matlab.ui.control.Switch
        graphing_corner_sensitivity     matlab.ui.container.Tab
        VariablesListBoxLabel           matlab.ui.control.Label
        sensitivity_vars                matlab.ui.control.ListBox
        pressCTRLtoselectmultipleLabel  matlab.ui.control.Label
        MetricButtonGroup               matlab.ui.container.ButtonGroup
        ymaxButton                      matlab.ui.control.RadioButton
        yunlatchButton                  matlab.ui.control.RadioButton
        tLButton                        matlab.ui.control.RadioButton
        vtoButton                       matlab.ui.control.RadioButton
        PmaxButton                      matlab.ui.control.RadioButton
        ttoButton                       matlab.ui.control.RadioButton
        KEmaxButton                     matlab.ui.control.RadioButton
        go                              matlab.ui.control.Button
        ShowModelSchematicButton        matlab.ui.control.Button
        Image                           matlab.ui.control.Image
        savesolutionCheckBox            matlab.ui.control.CheckBox
    end

    properties (Access = private)
        load_vars;
        latch_vars;
        spring_vars;
        lm_vars;
        um_vars;
        dd_vars;
        dd_vars_labels;
        axis_labels_dict;
        dropdown_items_dict;
        dropdown_items_opposite_dict;
    end
    
    methods (Access = private)
        function update_dd_vars(app)
            app.dd_vars = [app.load_vars app.latch_vars app.spring_vars app.lm_vars app.um_vars];
            for i=1:length(app.dd_vars)
                app.dd_vars_labels{i} = app.dropdown_items_dict(app.dd_vars{i});
            end
            app.IV1DropDown.Items = app.dd_vars_labels;
            app.IV2DropDown.Items = app.dd_vars_labels;
            app.OD_IV1DropDown.Items = app.dd_vars_labels;
            app.sensitivity_vars.Items = app.dd_vars_labels;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function heatmap(app)
            
            if strcmp(app.x_log_space.Value,'log')
                if ((app.xmin.Value <= 0) ||(app.xmax.Value <= 0))
                    warndlg("One or both of the values you've entered for the xmin/xmax " + ...
                        "looping range is non-positive, " + ...
                        "AND you have chosen a logspace x-axis. " + newline +...
                        "Either change to linear spacing or use a positive bound.", 'Warning')
                    return
                end
            end
            
            if strcmp(app.y_log_space.Value,'log')
                if ((app.ymin.Value <= 0) || (app.ymax.Value <= 0))
                    warndlg("One or both of the values you've entered for the ymin/ymax " + ...
                        "looping range is non-positive,..." + ...
                        "AND you have chosen a logspace y-axis. " + newline +...
                        "Either change to linear spacing or use a positive bound.", 'Warning')
                    return
                end
            end
            
            if ~(app.y_maxCheckBox.Value || app.y_unlatchCheckBox.Value || app.t_LCheckBox.Value || ...
                    app.v_toCheckBox.Value || app.P_maxCheckBox.Value || app.t_toCheckBox.Value || app.KE_maxCheckBox.Value)
                warndlg('You must select at least one output option', 'Warning');
                return
            end
            
            % determines resolution of heatplots
            N=app.n.Value;
            
            % output directory initialization
            if app.savesolutionCheckBox.Value
                output_directory = create_output_directory();
            end
            
            %% plot parameters
            
            looping_param_x = app.dropdown_items_opposite_dict(app.IV1DropDown.Value);
            looping_param_y = app.dropdown_items_opposite_dict(app.IV2DropDown.Value);
            eval(['looping_param_x_value = app.' app.dropdown_items_opposite_dict(app.IV1DropDown.Value) '.Value;']);
            eval(['looping_param_y_value = app.' app.dropdown_items_opposite_dict(app.IV2DropDown.Value) '.Value;']);
            
            eval(['looping_param_x_limits = app.' app.dropdown_items_opposite_dict(app.IV1DropDown.Value) '.Limits;']);
            eval(['looping_param_y_limits = app.' app.dropdown_items_opposite_dict(app.IV1DropDown.Value) '.Limits;']);
            
            if ((app.xmin.Value < looping_param_x_limits(1)) | ...
                (app.xmin.Value > looping_param_x_limits(2)) | ...
                (app.xmax.Value < looping_param_x_limits(1)) | ...
                (app.xmax.Value > looping_param_x_limits(2)))     
                warndlg("You have chosen an xmin or xmax value that is out of the bounds for the looping quantity you have selected." + newline + ...
                    "i.e. you have chosen to loop through negative masses, or negative coefficients of friction, or something similarly bad.", ...
                    'Warning')
                return
            end
            if ((app.ymin.Value < looping_param_y_limits(1)) | ...
                (app.ymin.Value > looping_param_y_limits(2)) | ...
                (app.ymax.Value < looping_param_y_limits(1)) | ...
                (app.ymax.Value > looping_param_y_limits(2)))     
                warndlg("You have chosen an ymin or ymax value that is out of the bounds for the looping quantity you have selected." + newline + ...
                    "i.e. you have chosen to loop through negative masses, or negative coefficients of friction, or something similarly bad.", ...
                    'Warning')
                return
            end
            
            % initializing waitbar
            f = waitbar(0,'Please wait...');
            load_bar_value = 1/N;
            load_bar_increment = 1/N;
            
            % setting x axis on the plot
            xname = app.dropdown_items_opposite_dict(app.IV1DropDown.Value);
            if strcmp(app.x_log_space.Value,'log')
                xrange = [log10(app.xmin.Value) log10(app.xmax.Value)];
                looping_value_x = logspace(xrange(1),xrange(2),N);
            else
                xrange = [app.xmin.Value app.xmax.Value];
                looping_value_x = linspace(xrange(1),xrange(2),N);
            end
            
            %setting y axis value on plot
            yname = app.dropdown_items_opposite_dict(app.IV2DropDown.Value);
            if strcmp(app.y_log_space.Value,'log')
                yrange = [log10(app.ymin.Value) log10(app.ymax.Value)];
                looping_value_y = logspace(yrange(1),yrange(2),N);
            else
                yrange = [app.ymin.Value app.ymax.Value];
                looping_value_y = linspace(yrange(1),yrange(2),N);
            end
            
            % add things to metrics
            metrics = {};
            metrics_names = {};
            metrics_labels = {};
            if app.y_maxCheckBox.Value
                metrics{end+1} = 'ymax';
                metrics_names{end+1} = 'ymax';
                metrics_labels{end+1} = '$y_{\textrm{max}}$';
            end
            if app.y_unlatchCheckBox.Value
                metrics{end+1} = 'yunlatch';
                metrics_names{end+1} = 'yunlatch';
                metrics_labels{end+1} = '$y_{\textrm{unlatch}}$';
            end
            if app.t_LCheckBox.Value
                metrics{end+1} = 'tL';
                metrics_names{end+1} = 'tL';
                metrics_labels{end+1} = '$t_L$';
            end
            if app.v_toCheckBox.Value
                metrics{end+1} = 'vto';
                metrics_names{end+1} = 'vto';
                metrics_labels{end+1} = '$v_{\textrm{to}}$';
            end
            if app.P_maxCheckBox.Value
                metrics{end+1} = 'Pmax';
                metrics_names{end+1} = 'Pmax';
                metrics_labels{end+1} = '$P_{\textrm{max}}$';
            end
            if app.t_toCheckBox.Value
                metrics{end+1} = 'tto';
                metrics_names{end+1} = 'tto';
                metrics_labels{end+1} = '$t_{\textrm{to}}$';
            end
            if app.minunlatchingmotorforceCheckBox.Value
                metrics{end+1} = 'minumforce';
                metrics_names{end+1} = 'minumforce';
                metrics_labels{end+1} = 'min unlatching force';
            end
            if app.KE_maxCheckBox.Value
                metrics{end+1} = 'KEmax';
                metrics_names{end+1} = 'KEmax';
                metrics_labels{end+1} = '$KE_{\textrm{max}}$';
            end
%             if app.KERatioCheckBox.Value
%                 metrics{end+1} = 'KE_Ratio';
%                 metrics_names{end+1} = 'KE_Ratio';
%                 metrics_labels{end+1} = 'KE Ratio';
%             end
%             if app.unlatchingmotorworkdoneCheckBox.Value
%                 metrics{end+1} = 'unlatching_motor_work_done';
%                 metrics_names{end+1} = 'unlatching_motor_work_done';
%                 metrics_labels{end+1} = 'unlatching motor work done';
%             end
            
            
            if isempty(metrics)
                warndlg('You must pick at least one output option','Error');
                return
            end
            
            
            metrics_dict = containers.Map(metrics_names,metrics_labels);
            
            
%             looping_param_x = app.dropdown_items_opposite_dict(app.IV1DropDown.Value);
%             looping_param_y = app.dropdown_items_opposite_dict(app.IV2DropDown.Value);
%             eval(['looping_param_x_value = app.' app.dropdown_items_opposite_dict(app.IV1DropDown.Value) '.Value;']);
%             eval(['looping_param_y_value = app.' app.dropdown_items_opposite_dict(app.IV2DropDown.Value) '.Value;']);
            
            % ensures that both independent variables that we are varying over are not the same
            if (strcmp(looping_param_y,looping_param_x))
                errordlg('The two independent variables must be different. Please change them and retry.','Error');
                return
            end
            
            for ii=1:length(metrics)
                outval{ii}=zeros(N);
            end
            for i=1:N %iterate over y-axis-variable of plot
                for j=1:N %iterate over x-axis-variable of plot       
                    
                    eval(['app.' looping_param_x '.Value = ' num2str(looping_value_x(j)) ';'])
                    eval(['app.' looping_param_y '.Value = ' num2str(looping_value_y(i)) ';'])
                    %% initializing LaMSA component structs
                    
                    % load mass struct initialization
                    load = load_mass(app.load_mass_mass.Value,app.load_m_rod.Value,app.load_EMA.Value);
                    
                    % latch struct initialization
                    latch = rounded_latch(app.latch_radius.Value,app.latch_mass.Value,app.latch_coeff_fric.Value, app.latch_v_0.Value, app.min_latching_dist.Value, app.max_latching_dist.Value, app.runway_length.Value);
                    
                    % spring struct initialization
                    if (app.spring.SelectedTab == app.linear_spring)
                        spring = linear_spring(app.linear_spring_k.Value,app.linear_spring_mass.Value,app.linear_spring_Fmax.Value);
                    elseif (app.spring.SelectedTab == app.exponential_spring)
%                         disp("THIS IS THE APP CODE")
%                         app.exp_spring_k.Value
%                         app.exp_spring_char_len.Value
%                         app.exp_spring_mass.Value
%                         app.exp_spring_Fmax.Value
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
                    
                    % calling solve model
                    if app.savesolutionCheckBox.Value
                        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
                    else
                        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);
                    end
                    
                    % warning messages
                    [warnMsg, warnId] = lastwarn;
                    if (strcmp(warnMsg,'Loading failed. Does not fall within latching distance conditions.'))
                        if ~exist('loading_failed','var')
                            loading_failed = warndlg('Loading failed. Does not fall within latching distance conditions.','Warning');
                        end
                    end
                    if (strcmp(warnMsg,'Latch gets stuck!'))
                        if ~exist('latch_stuck','var')
                            latch_stuck = warndlg('Latch gets stuck!','Warning');
                        end
                    end
                    if (strcmp(warnMsg,"The latch's initial velocity and acceleration are both zero."))
                        if ~exist('latch_vi_a_zero','var')
                            latch_vi_a_zero = warndlg("The latch's initial velocity and acceleration are both zero.",'Warning');
                        end
                    end

                    
                    met_dict=get_metrics(sol,transition_times,load,metrics);
                    for ii=1:length(metrics)
                        % the KE_Ratio gets special treatment because it's
                        % weird. This is the ratio of the load mass's final
                        % kinetic energy in braking vs. non-braking
                        % scenarios. The difference between these scenarios
                        % is that in one scenario, we allow a linear
                        % unlatching motor to behave as normal (exerts a
                        % negative force on the latch if the motor moves
                        % faster than vmax), whereas in the other, we
                        % specify that if the unlatching motor is moving so
                        % fast that it would normally brake, we
                        % artificially set the force to 0 instead. 
                        % Comparing these two physical scenarios requires a
                        % second call to solve_model, hence the special
                        % treatment.
%                         if (strcmp(metrics{ii},'KE_Ratio') && app.unlatching_motor.SelectedTab == app.um_linear_motor)
% 
%                             unlatching_motor_no_braking = linear_motor(app.um_linear_motor_Fmax.Value,app.um_linear_motor_Vmax.Value,app.um_linear_motor_range_of_motion.Value,app.um_linear_motor_voltage_frac.Value, true);
%                             [sol_no_braking, tt_no_braking] = solve_model(loading_motor,unlatching_motor_no_braking,load,latch,spring,output_directory);
%                             
%                             KE_no_braking = (0.5*load.mass*(sol_no_braking(end,3)^2));
%                             KE_braking = (0.5*load.mass*(sol(end,3)^2));
%                             ratio = KE_no_braking/KE_braking;
%                             outval{ii}(i,j)=ratio;
%                         elseif (strcmp(metrics{ii},'KE_Ratio') && app.unlatching_motor.SelectedTab ~= app.um_linear_motor)
%                             error("The KE_Ratio metric is only available for a linear unlatching motor (for now!)")
%                         else
%                             outval{ii}(i,j)=met_dict(metrics{ii});
%                         end
                        outval{ii}(i,j)=met_dict(metrics{ii});
                        metrics{ii};
                    end
                end
                disp(['row ' num2str(i) ' of ' num2str(N)]);
                load_bar_value = load_bar_value + load_bar_increment;
                waitbar(load_bar_value,f,'Processing...');
            end
            
            close(f)
            
            eval(['app.' app.dropdown_items_opposite_dict(app.IV1DropDown.Value) '.Value = looping_param_x_value;']);
            eval(['app.' app.dropdown_items_opposite_dict(app.IV2DropDown.Value) '.Value = looping_param_y_value;']);
            
            % plot output
            fh = figure('Name','Heatmaps');
            fh.WindowState = 'maximized';
            subplot_rows = floor(sqrt(length(metrics)));
            subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
            for ii=1:length(metrics)
                subplot(subplot_rows,subplot_cols,ii);
                
                % so that you can click on the heatplot and have a
                % kinematics
                % plot show up
                imageHandle = imagesc(xrange,yrange,outval{ii});
                set(gca,'YDir','normal');
                xlabel(app.axis_labels_dict(xname),'Interpreter', 'Latex');
                ylabel(app.axis_labels_dict(yname), 'Interpreter', 'Latex');
                
                c = colorbar;
                colormap(linspecer);
                c.Label.Interpreter = 'Latex';
                set(c,'TickLabelInterpreter','Latex');
                c.Label.String = metrics_dict(metrics_names{ii});
                
                
                % makes the x and y axes in log scale
                set(gca,'TickLabelInterpreter','Latex');
                
                xtick = get(gca,'XTick');
                ytick = get(gca,'YTick');
                xticklabel = get(gca,'XTickLabel');
                yticklabel = get(gca,'YTickLabel');
                if strcmp(app.x_log_space.Value,'log')
                    for i=1:length(xticks)
                        xticklabel{i} = strcat('$10^{',num2str(xtick(i)),'}$');
                    end
                end
                set(gca,'XTickLabel',xticklabel);
                set(gca,'XTick',xtick);
                if strcmp(app.y_log_space.Value,'log')
                    for j=1:length(yticks)
                        yticklabel{j} = strcat('$10^{',num2str(ytick(j)),'}$');
                    end
                end
                set(gca,'YTickLabel',yticklabel);
                set(gca,'YTick',ytick);
                set(imageHandle,'ButtonDownFcn',{@showKinematics, looping_param_x, looping_param_y, app});
                
            end
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function one_D_plot(app)
            if ~(app.OD_y_maxCheckBox.Value || app.OD_y_unlatchCheckBox.Value || app.OD_t_LCheckBox.Value || ...
                    app.OD_v_toCheckBox.Value || app.OD_P_maxCheckBox.Value || app.OD_t_toCheckBox.Value || app.OD_KE_maxCheckBox.Value)
                warndlg('You must select at least one output option', 'Warning');
                return
            end
            if strcmp(app.OD_x_log_space.Value,'log')
                if ((app.OD_xmin.Value <= 0) |(app.OD_xmax.Value <= 0))
                    warndlg("One or both of the values you've entered for the xmin/xmax " + ...
                        "looping range is non-positive, " + ...
                        "AND you have chosen a logspace x-axis. " + newline +...
                        "Either change to linear spacing or use a positive bound.", 'Warning')
                    return
                end
            end
            
            % determines resolution of heatplots
            N=app.OD_n.Value; 
            
            % output directory initialization
            if app.savesolutionCheckBox.Value
                output_directory = create_output_directory();
            end
            
            %% plot parameters
            looping_param_x = app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value);
            eval(['looping_param_x_value = app.' app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value) '.Value;']);
            eval(['looping_param_x_limits = app.' app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value) '.Limits;']);
            
            if ((app.OD_xmin.Value < looping_param_x_limits(1)) | ...
                (app.OD_xmin.Value > looping_param_x_limits(2)) | ...
                (app.OD_xmax.Value < looping_param_x_limits(1)) | ...
                (app.OD_xmax.Value > looping_param_x_limits(2)))     
                warndlg("You have chosen an xmin or xmax value that is out of the bounds for the looping quantity you have selected." + newline + ...
                    "i.e. you have chosen to loop through negative masses, or negative coefficients of friction, or something similarly bad.", ...
                    'Warning')
                return
            end
            
            % initializing waitbar
            f = waitbar(0,'Please wait...');
            load_bar_value = 1/N;
            load_bar_increment = 1/N;
            
            % setting x axis on the plot
            xname = app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value);
            if strcmp(app.OD_x_log_space.Value,'log')
                xrange = [log10(app.OD_xmin.Value) log10(app.OD_xmax.Value)];
                looping_value_x = logspace(xrange(1),xrange(2),N);
            else
                xrange = [app.OD_xmin.Value app.OD_xmax.Value];
                looping_value_x = linspace(xrange(1),xrange(2),N);
            end
            
            % add things to metrics
            metrics = {};
            
            if app.OD_y_maxCheckBox.Value
                metrics{end+1} = 'ymax';
            end
            if app.OD_y_unlatchCheckBox.Value
                metrics{end+1} = 'yunlatch';
            end
            if app.OD_t_LCheckBox.Value
                metrics{end+1} = 'tL';
            end
            if app.OD_v_toCheckBox.Value
                metrics{end+1} = 'vto';
            end
            if app.OD_P_maxCheckBox.Value
                metrics{end+1} = 'Pmax';
            end
            if app.OD_t_toCheckBox.Value
                metrics{end+1} = 'tto';
            end
            if app.OD_minunlatchingmotorforceCheckBox.Value
                metrics{end+1} = 'minumforce';
            end
            if app.OD_KE_maxCheckBox.Value
                metrics{end+1} = 'KEmax';
            end
%             if app.OD_KE_ratioCheckBox.Value
%                 metrics{end+1} = 'KE_Ratio';
%             end
%             if app.OD_unlatchingmotorworkdoneCheckBox.Value
%                 metrics{end+1} = 'unlatching_motor_work_done';
%             end
            
            if isempty(metrics)
                warndlg('You must pick at least one output option','Error');
                return
            end
            
            metrics_names = {'ymax','yunlatch','tL','vto','Pmax','tto','minumforce','KEmax', 'KE_Ratio', 'unlatching_motor_work_done'};
            metrics_labels = {'$y_{\textrm{max}}$','$y_{\textrm{unlatch}}$','$t_L$','$v_{\textrm{to}}$','$P_{\textrm{max}}$','$t_{\textrm{to}}$','min unlatching force','$KE_{\textrm{max}}$', 'KE Ratio', 'unlatching motor work done'};
            metrics_dict = containers.Map(metrics_names,metrics_labels);
            
            for ii=1:length(metrics)
                outval{ii} = zeros(size(looping_value_x));
            end
            for i=1:N    
                eval(['app.' looping_param_x '.Value = ' num2str(looping_value_x(i)) ';']);
                
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
                
                % calling solve model
                if app.savesolutionCheckBox.Value
                    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
                else
                    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);
                end
                
                % warning messages
                [warnMsg, warnId] = lastwarn;
                if (strcmp(warnMsg,'Loading failed. Does not fall within latching distance conditions.'))
                    if ~exist('loading_failed','var')
                        loading_failed = warndlg('Loading failed. Does not fall within latching distance conditions.','Warning');
                    end
                end
                if (strcmp(warnMsg,'Latch gets stuck!'))
                    if ~exist('latch_stuck','var')
                        latch_stuck = warndlg('Latch gets stuck!','Warning');
                    end
                end
                if (strcmp(warnMsg,"The latch's initial velocity and acceleration are both zero."))
                    if ~exist('latch_vi_a_zero','var')
                        latch_vi_a_zero = warndlg("The latch's initial velocity and acceleration are both zero.",'Warning');
                    end
                end
                
                met_dict=get_metrics(sol,transition_times,load,metrics);
                for ii=1:length(metrics)
                    % the KE_Ratio gets special treatment because it's
                    % weird. This is the ratio of the load mass's final
                    % kinetic energy in braking vs. non-braking
                    % scenarios. The difference between these scenarios
                    % is that in one scenario, we allow a linear
                    % unlatching motor to behave as normal (exerts a
                    % negative force on the latch if the motor moves
                    % faster than vmax), whereas in the other, we
                    % specify that if the unlatching motor is moving so
                    % fast that it would normally brake, we
                    % artificially set the force to 0 instead. 
                    % Comparing these two physical scenarios requires a
                    % second call to solve_model, hence the special
                    % treatment.
                    if (strcmp(metrics{ii},'KE_Ratio') && app.unlatching_motor.SelectedTab == app.um_linear_motor)
                        
                        unlatching_motor_no_braking = linear_motor(app.um_linear_motor_Fmax.Value,app.um_linear_motor_Vmax.Value,app.um_linear_motor_range_of_motion.Value,app.um_linear_motor_voltage_frac.Value, true);
                        [sol_no_braking, tt_no_braking] = solve_model(loading_motor,unlatching_motor_no_braking,load,latch,spring,output_directory);
                        KE_no_braking = (0.5*load.mass*(sol_no_braking(end,3)^2));
                        
                        KE_braking = (0.5*load.mass*(sol(end,3)^2));
                        ratio = KE_no_braking/KE_braking;
                        outval{ii}(i)=ratio;
                    elseif (strcmp(metrics{ii},'KE_Ratio') && app.unlatching_motor.SelectedTab ~= app.um_linear_motor)
                        error("The KE_Ratio metric is only available for a linear unlatching motor (for now!)")
                    else
                        outval{ii}(i)=met_dict(metrics{ii});
                    end
                end
                
                disp(['point ' num2str(i) ' of ' num2str(N)]);
                load_bar_value = load_bar_value + load_bar_increment;
                waitbar(load_bar_value,f,'Processing...');
            end

            close(f)
            eval(['app.' app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value) '.Value = looping_param_x_value;']);
            
            % plot output
            fh = figure('Name','1D Plot');
            fh.WindowState = 'maximized';
            subplot_rows = floor(sqrt(length(metrics)));
            subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
            for ii=1:length(metrics)
                subplot(subplot_rows,subplot_cols,ii);
                plot(looping_value_x,outval{ii},'.'); 
                xlabel(app.axis_labels_dict(xname),'Interpreter', 'Latex');
                ylabel(metrics_dict(metrics{ii}), 'Interpreter', 'Latex');
                ylim_range = max(outval{ii})-min(outval{ii});
                if (ylim_range == 0)
                    ylim_range = abs(max(outval{ii}));
                    if ~(ylim_range == 0)
                        ylim([min([0, min(outval{ii})-0.05*ylim_range]) max([0,max(outval{ii})+0.05*ylim_range]) ])
                    end
                else
                    ylim([min(outval{ii})-0.05*ylim_range max(outval{ii})+0.05*ylim_range])
                end
                set(gca,'TickLabelInterpreter','latex')
                % makes the x in log scale if needed
                if strcmp(app.OD_x_log_space.Value,'log')
                    set(gca,'XScale','log')
                end
               
            end
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function kinematics(app)
            if ~(app.loadkinematicsCheckBox.Value || app.latchkinematicsCheckBox.Value || app.forcedisp.Value)
                warndlg('You must select at least one kinematics option', 'Warning');
                return
            end
            
            f = waitbar(0,'Please wait...');
            load_bar_value = 1/4;
            load_bar_increment = 1/4;
            
            % initializing output directory
            if app.savesolutionCheckBox.Value
                output_directory = create_output_directory();
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
            
            load_bar_value = load_bar_value + 2*load_bar_increment;
            waitbar(load_bar_value,f,'Processing...');
            
            % loading motor struct initialization
            if (app.loading_motor.SelectedTab == app.lm_linear_motor)
                loading_motor = linear_motor(app.lm_linear_motor_Fmax.Value,app.lm_linear_motor_Vmax.Value,app.lm_linear_motor_range_of_motion.Value,app.lm_linear_motor_voltage_frac.Value);
            elseif (app.loading_motor.SelectedTab == app.lm_hill_muscle_motor)
                loading_motor = hill_muscle_motor(app.lm_hill_motor_muscle_length.Value,app.lm_hill_motor_Fmax.Value,app.lm_hill_motor_Vmax.Value,app.lm_hill_motor_rate_of_activation.Value,app.lm_hill_motor_L_i.Value,app.lm_hill_motor_a_L.Value,app.lm_hill_motor_b_L.Value,app.lm_hill_motor_s.Value);
            end
            
            load_bar_value = load_bar_value + load_bar_increment;
            waitbar(load_bar_value,f,'Processing...');
            
            % unlatching motor struct initialization
            if (app.unlatching_motor.SelectedTab == app.um_linear_motor)
                unlatching_motor = linear_motor(app.um_linear_motor_Fmax.Value,app.um_linear_motor_Vmax.Value,app.um_linear_motor_range_of_motion.Value,app.um_linear_motor_voltage_frac.Value);
            elseif (app.unlatching_motor.SelectedTab == app.um_hill_muscle_motor)
                unlatching_motor = hill_muscle_motor(app.um_hill_motor_muscle_length.Value,app.um_hill_motor_Fmax.Value,app.um_hill_motor_Vmax.Value,app.um_hill_motor_rate_of_activation.Value,app.um_hill_motor_L_i.Value,app.um_hill_motor_a_L.Value,app.um_hill_motor_b_L.Value,app.um_hill_motor_s.Value);
            end
            
            if app.savesolutionCheckBox.Value
                    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
                else
                    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);
            end
            
            columntitles=["Time", "y postition", "y velocity", "x position", "x velocity", "normal force on latch x", ...
                "normal force on load y", "frictional force on latch x", ...
                "frictional force on load y", "spring force", ...
                "unlatching motor force into"];
            
            close(f)
            
            % latch kinematics
            if (app.latchkinematicsCheckBox.Value)    
                figure('Name','Latch Kinematics')
                for i = 4:5
                    subplot(1,3,i-3)
                    box on
                    set(gca,'TickLabelInterpreter','latex')
                    hold on
                    plot(sol(:,1),sol(:,i));
                    hold off
                    title(columntitles(i),"Interpreter","latex");
                    ylabel(columntitles(i),"Interpreter","latex");
                    xlabel(columntitles(1),"Interpreter","latex");
                end
                
                for i = [6 8 11]
                    subplot(1,3,3)
                    box on
                    set(gca,'TickLabelInterpreter','latex')
                    hold on
                    plot(sol(:,1),sol(:,i),"DisplayName",columntitles(i))
                end
                title("Latch Force Components","Interpreter","latex");
                ylabel("Force","Interpreter","latex");
                xlabel(columntitles(1),"Interpreter","latex");
                legend("show")
                hold off
            end
            
            % load kinematics
            if (app.loadkinematicsCheckBox.Value)
                figure('Name','Load Kinematics')
                for i = 2:3
                    subplot(1,3,i-1)
                    box on
                    set(gca,'TickLabelInterpreter','latex')
                    hold on
                    plot(sol(:,1),sol(:,i));
                    hold off
                    title(columntitles(i),"Interpreter","latex");
                    ylabel(columntitles(i),"Interpreter","latex");
                    xlabel(columntitles(1),"Interpreter","latex");
                end
                
                for i = [7 9 10]
                    subplot(1,3,3)
                    box on
                    set(gca,'TickLabelInterpreter','latex')
                    hold on
                    plot(sol(:,1),sol(:,i),"DisplayName",columntitles(i))
                end
                title("Load Force Components","Interpreter","latex");
                ylabel("Force","Interpreter","latex");
                xlabel(columntitles(1),"Interpreter","latex");
                legend("show")
                hold off
            end    
            
            % force displacement curves for motor and spring
            if (app.forcedisp.Value)    
                figure('Name','Motor and Spring Force Displacement Curves')
                if (loading_motor.range == Inf)
                    y_max = 5*sol(1,2);
                else
                    y_max = -loading_motor.range;
                end
                y_range = linspace(0,y_max,100);
                motor_force_array = [];
                spring_force_array = [];
                for i = 1:length(y_range)
                    motor_force_array(end+1) = loading_motor.Force(Inf,[-y_range(i),0]);
                    spring_force_array(end+1) = spring.Force(Inf,[y_range(i),0]);
                end
                force_arrays = {spring_force_array motor_force_array};
                for i = 1:2
                    box on
                    set(gca,'TickLabelInterpreter','latex')
                    hold on
                    plot(abs(y_range),force_arrays{i});
                end
                plot([-sol(1,2),-sol(1,2)],[0,spring.Force(Inf,[sol(1,2),0])],"--k")
                plot(-sol(1,2),spring.Force(Inf,[sol(1,2),0]),'ok')
                title('Force Displacement Curve',"Interpreter","latex");
                ylabel('Force',"Interpreter","latex");
                xlabel('Displacement',"Interpreter","latex");
                legend('spring','motor','loaded displacement')
            end
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sensitivity(app)
            f = waitbar(1/4,'Please wait...');
            load_bar_value = 1/4;
            load_bar_increment = 1/4;
            
            labels = strings([1 length(app.sensitivity_vars.Value)]);
            x0 = zeros([1 length(labels)]);
            for i=1:length(labels)
                labels(i) = app.sensitivity_vars.Value(i);
                x0(i) = eval(['app.' app.dropdown_items_opposite_dict(labels(i)) '.Value;']);
            end
            load_str = "mass = @(x0)load_mass(app.load_mass_mass.Value,app.load_m_rod.Value,app.load_EMA.Value);";
            
            latch_str = "latch = @(x0)rounded_latch(app.latch_radius.Value,app.latch_mass.Value,app.latch_coeff_fric.Value, app.latch_v_0.Value, app.min_latching_dist.Value, app.max_latching_dist.Value, app.runway_length.Value);";
            
            if (app.spring.SelectedTab == app.linear_spring)
                spring_str = "spring_fun = @(x0)linear_spring(app.linear_spring_k.Value,app.linear_spring_mass.Value,app.linear_spring_Fmax.Value);";
            elseif (app.spring.SelectedTab == app.exponential_spring)
                spring_str = "spring_fun = @(x0)exponential_spring(app.exp_spring_k.Value,app.exp_spring_char_len.Value,app.exp_spring_mass.Value,app.exp_spring_Fmax.Value);";
            elseif (app.spring.SelectedTab == app.linear_elastic_extensional_spring)
                spring_str = "spring_fun = @(x0)linear_elastic_extensional_spring(app.lee_spring_E.Value,app.lee_spring_A.Value,app.lee_spring_L.Value,app.lee_spring_rho.Value,app.lee_spring_sigma_f.Value);";
            end

            if (app.loading_motor.SelectedTab == app.lm_linear_motor)
                loading_motor_str = "loading_motor = @(x0)linear_motor(app.lm_linear_motor_Fmax.Value,app.lm_linear_motor_Vmax.Value,app.lm_linear_motor_range_of_motion.Value,app.lm_linear_motor_voltage_frac.Value);";
            elseif (app.loading_motor.SelectedTab == app.lm_hill_muscle_motor)
                loading_motor_str = "loading_motor = @(x0)hill_muscle_motor(app.lm_hill_motor_muscle_length.Value,app.lm_hill_motor_Fmax.Value,app.lm_hill_motor_Vmax.Value,app.lm_hill_motor_rate_of_activation.Value,app.lm_hill_motor_L_i.Value,app.lm_hill_motor_a_L.Value,app.lm_hill_motor_b_L.Value,app.lm_hill_motor_s.Value);";
            end

            if (app.unlatching_motor.SelectedTab == app.um_linear_motor)
                unlatching_motor_str = "unlatching_motor = @(x0)linear_motor(app.um_linear_motor_Fmax.Value,app.um_linear_motor_Vmax.Value,app.um_linear_motor_range_of_motion.Value,app.um_linear_motor_voltage_frac.Value);";
            elseif (app.unlatching_motor.SelectedTab == app.um_hill_muscle_motor)
                unlatching_motor_str = "unlatchin_motor = @(x0)hill_muscle_motor(app.um_hill_motor_muscle_length.Value,app.um_hill_motor_Fmax.Value,app.um_hill_motor_Vmax.Value,app.um_hill_motor_rate_of_activation.Value,app.um_hill_motor_L_i.Value,app.um_hill_motor_a_L.Value,app.um_hill_motor_b_L.Value,app.um_hill_motor_s.Value);";
            end
            
            load_bar_value = load_bar_value + load_bar_increment;
            waitbar(load_bar_value,f,'Processing...');
            
            for i=1:length(labels)
                load_str = strrep(load_str,strcat("app.",app.dropdown_items_opposite_dict(labels(i)),".Value"),strcat("x0(",num2str(i),")"));
                latch_str = strrep(latch_str,strcat("app.",app.dropdown_items_opposite_dict(labels(i)),".Value"),strcat("x0(",num2str(i),")"));
                spring_str = strrep(spring_str,strcat("app.",app.dropdown_items_opposite_dict(labels(i)),".Value"),strcat("x0(",num2str(i),")"));
                loading_motor_str = strrep(loading_motor_str,strcat("app.",app.dropdown_items_opposite_dict(labels(i)),".Value"),strcat("x0(",num2str(i),")"));
                unlatching_motor_str = strrep(unlatching_motor_str,strcat("app.",app.dropdown_items_opposite_dict(labels(i)),".Value"),strcat("x0(",num2str(i),")"));
            end
            eval(load_str);
            eval(latch_str);
            eval(spring_str);
            eval(loading_motor_str);
            eval(unlatching_motor_str);
            
            load_bar_value = load_bar_value + load_bar_increment;
            waitbar(load_bar_value,f,'Processing...');
            
            [combos,var_list] = sensitivity_analysis(loading_motor,unlatching_motor,mass,latch,spring_fun,x0,app.MetricButtonGroup.SelectedObject.Text,labels);
            commandwindow;
            disp(' ');
            disp('----------------------');
            disp('Sensitivity Analysis');
            disp('----------------------');
            disp(' ');
            column_labels = strings([1, size(combos,2)]);
            for i=1:length(column_labels)
                column_labels(i) = "combo #"+num2str(i);
            end
            combo_table = array2table(combos,'RowNames',labels,'VariableNames',column_labels);
            disp('Most Sensitive Combinations of Parameters')
            disp(combo_table);
            disp('Most Sensitive Parameter Axes')
            for i=1:length(var_list)
                disp(var_list(i));
            end
            close(f)
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            drawnow;
            app.UIFigure.WindowState = 'maximized';
            % makes comparing changes in the app in github easier
            pathname = char(which("plot_app.mlapp"));
            lastIndex = strlength(pathname);
            nameLength = strlength("plot_app.mlapp");
            correctIndex = lastIndex-nameLength;
            pathname = pathname(1:correctIndex);
            pathname = pathname +"/..";
            addpath(genpath(fullfile(pathname)));
            mlapp2classdef("plot_app.mlapp");
            
            varnames = {'load_mass_mass','load_m_rod','load_EMA',...
                'latch_mass','latch_coeff_fric','latch_radius','latch_v_0','min_latching_dist','max_latching_dist','runway_length'...
                'linear_spring_k','linear_spring_mass','linear_spring_Fmax',...
                'exp_spring_k','exp_spring_char_len','exp_spring_Fmax','exp_spring_mass',...
                'lee_spring_E','lee_spring_A','lee_spring_L','lee_spring_rho','lee_spring_sigma_f',...
                'lm_linear_motor_Fmax','lm_linear_motor_Vmax','lm_linear_motor_range_of_motion','lm_linear_motor_voltage_frac'...
                'lm_hill_motor_Fmax','lm_hill_motor_Vmax','lm_hill_motor_muscle_length','lm_hill_motor_rate_of_activation','lm_hill_motor_L_i','lm_hill_motor_a_L','lm_hill_motor_b_L','lm_hill_motor_s',...
                'um_linear_motor_Fmax','um_linear_motor_Vmax','um_linear_motor_range_of_motion','um_linear_motor_voltage_frac'...
                'um_hill_motor_Fmax','um_hill_motor_Vmax','um_hill_motor_muscle_length','um_hill_motor_rate_of_activation','um_hill_motor_L_i','um_hill_motor_a_L','um_hill_motor_b_L','um_hill_motor_s'};
            latexlabels = {'load mass','lever arm mass','EMA'...
                'latch mass','latch $\mu$','latch radius','latch $v_0$','min latching distance','max latching distance','runway length'...
                'linear spring k','linear spring mass','linear spring Fmax',...
                'exponential spring $k_0$','exponential spring characteristic length','exponential spring Fmax','exponential spring mass',...
                'spring E','spring A','spring L','spring $\rho$','spring $\sigma_f$',...
                'loading motor (linear) Fmax','loading motor (linear) motor Vmax','loading motor (linear) range of motion','loading motor (linear) voltage fraction'...
                'loading motor (Hill) Fmax','loading motor (Hill) Vmax','loading motor (Hill) muscle length','loading motor (Hill) rate of activation','loading motor (Hill) motor $L_i$','loading motor (Hill) $a_L$','loading motor (Hill) $b_L$','loading motor (Hill) s',...
                'unlatching motor (linear) Fmax','unlatching motor (linear) motor Vmax','unlatching motor (linear) range of motion','unlatching motor (linear) voltage fraction'...
                'unlatching motor (Hill) Fmax','unlatching motor (Hill) Vmax','unlatching motor (Hill) muscle length','unlatching motor (Hill) rate of activation','unlatching motor (Hill) motor $L_i$','unlatching motor (Hill) $a_L$','unlatching motor (Hill) $b_L$','unlatching motor (Hill) s'};
            nonlatexlabels = {'load mass','lever arm mass','EMA'...
                'latch mass','latch mu','latch radius','latch v_0','min latching distance','max latching distance','runway length'...
                'linear spring k','linear spring mass','linear spring Fmax',...
                'exponential spring k_0','exponential spring characteristic length','exponential spring Fmax','exponential spring mass',...
                'spring E','spring A','spring L','spring rho','spring sigma_f',...
                'loading motor (linear) Fmax','loading motor (linear) motor Vmax','loading motor (linear) range of motion', 'loading motor (linear) voltage fraction'...
                'loading motor (Hill) Fmax','loading motor (Hill) Vmax','loading motor (Hill) muscle length','loading motor (Hill) rate of activation','loading motor (Hill) motor L_i','loading motor (Hill) a_L','loading motor (Hill) b_L','loading motor (Hill) s',...
                'unlatching motor (linear) Fmax','unlatching motor (linear) motor Vmax','unlatching motor (linear) range of motion','unlatching motor (linear) voltage fraction'...
                'unlatching motor (Hill) Fmax','unlatching motor (Hill) Vmax','unlatching motor (Hill) muscle length','unlatching motor (Hill) rate of activation','unlatching motor (Hill) motor $L_i$','unlatching motor (Hill) $a_L$','unlatching motor (Hill) $b_L$','unlatching motor (Hill) s'};
            
            app.axis_labels_dict = containers.Map(varnames,latexlabels);
            app.dropdown_items_dict = containers.Map(varnames,nonlatexlabels);
            app.dropdown_items_opposite_dict = containers.Map(nonlatexlabels,varnames);

            app.sensitivity_vars.Multiselect = 'on';
            
            app.load_vars = {'load_mass_mass','load_m_rod','load_EMA'};
            app.latch_vars = {'latch_mass','latch_coeff_fric','latch_radius','latch_v_0','min_latching_dist','max_latching_dist','runway_length'};
            app.spring_vars = {'exp_spring_k','exp_spring_char_len','exp_spring_Fmax','exp_spring_mass'};
            %app.spring_vars = {'linear_spring_k','linear_spring_mass','linear_spring_Fmax'};
            %app.lm_vars = {'lm_linear_motor_Fmax','lm_linear_motor_Vmax','lm_linear_motor_range_of_motion','lm_linear_motor_voltage_frac'};
            app.lm_vars = {'lm_hill_motor_Fmax','lm_hill_motor_Vmax','lm_hill_motor_muscle_length','lm_hill_motor_rate_of_activation','lm_hill_motor_L_i','lm_hill_motor_a_L','lm_hill_motor_b_L','lm_hill_motor_s'};
            app.um_vars = {'um_linear_motor_Fmax','um_linear_motor_Vmax','um_linear_motor_range_of_motion','um_linear_motor_voltage_frac'};
            update_dd_vars(app);
            
            app.IV2DropDown.Items(ismember(app.IV2DropDown.Items,app.IV1DropDown.Value)) = [];
            
            % parameter defaults
            app.spring.SelectedTab = app.exponential_spring;
            app.loading_motor.SelectedTab = app.lm_hill_muscle_motor;
            
            % 1D plots defaults
            app.OD_IV1DropDown.Value = 'loading motor (Hill) Fmax';
            
            % 2D plots defaults
            app.IV1DropDown.Value = 'exponential spring k_0';
            app.IV2DropDown.Value = 'loading motor (Hill) Fmax';
        end

        % Button pushed function: go
        function goButtonPushed(app, event)
            lastwarn('');
            
            if (app.graphing_corner.SelectedTab == app.graphing_corner_heatmap)
                heatmap(app);
            elseif (app.graphing_corner.SelectedTab == app.graphing_corner_one_D)
                one_D_plot(app);
            elseif (app.graphing_corner.SelectedTab == app.graphing_corner_kinematics)
                kinematics(app);
            elseif (app.graphing_corner.SelectedTab == app.graphing_corner_sensitivity)
                sensitivity(app);
            end
        end

        % Value changed function: lm_hill_motor_muscle_length
        function lm_hill_motor_muscle_lengthValueChanged(app, event)
            % automatically sets lm_hill_motor_L_i to the
            % lm_hill_motor_muscle_length value
            value = app.lm_hill_motor_muscle_length.Value;
            app.lm_hill_motor_L_i.Value = value;
        end

        % Value changed function: um_hill_motor_muscle_length
        function um_hill_motor_muscle_lengthValueChanged(app, event)
            % automatically sets um_hill_motor_L_i to the
            % um_hill_motor_muscle_length value
            value = app.um_hill_motor_muscle_length.Value;
            app.um_hill_motor_L_i.Value = value;
        end

        % Value changed function: IV1DropDown
        function IV1DropDownValueChanged(app, event)
            value = app.IV1DropDown.Value;
            default = app.dd_vars_labels;
            app.IV2DropDown.Items = default;
            app.IV2DropDown.Items(ismember(app.IV2DropDown.Items,value)) = [];
        end

        % Value changed function: IV2DropDown
        function IV2DropDownValueChanged(app, event)
            value = app.IV2DropDown.Value;
            default = app.dd_vars_labels;
            app.IV1DropDown.Items = default;
            app.IV1DropDown.Items(ismember(app.IV1DropDown.Items,value)) = [];
        end

        % Selection change function: spring
        function springSelectionChanged(app, event)
            selectedTab = app.spring.SelectedTab;
            if (selectedTab == app.linear_spring)
                app.spring_vars = {'linear_spring_k','linear_spring_mass','linear_spring_Fmax'};
            elseif (selectedTab == app.exponential_spring)
                app.spring_vars = {'exp_spring_k','exp_spring_char_len','exp_spring_Fmax','exp_spring_mass'};
            elseif (selectedTab == app.linear_elastic_extensional_spring)
                app.spring_vars = {'lee_spring_E','lee_spring_A','lee_spring_L','lee_spring_rho','lee_spring_sigma_f'};
            end
            update_dd_vars(app);
        end

        % Selection change function: loading_motor
        function loading_motorSelectionChanged(app, event)
            selectedTab = app.loading_motor.SelectedTab;
            if (selectedTab == app.lm_linear_motor)
                app.lm_vars = {'lm_linear_motor_Fmax','lm_linear_motor_Vmax','lm_linear_motor_range_of_motion','lm_linear_motor_voltage_frac'};
            elseif (selectedTab == app.lm_hill_muscle_motor)
                app.lm_vars = {'lm_hill_motor_Fmax','lm_hill_motor_Vmax','lm_hill_motor_muscle_length','lm_hill_motor_rate_of_activation','lm_hill_motor_L_i','lm_hill_motor_a_L','lm_hill_motor_b_L','lm_hill_motor_s'};
            end
            update_dd_vars(app);
        end

        % Selection change function: unlatching_motor
        function unlatching_motorSelectionChanged(app, event)
            selectedTab = app.unlatching_motor.SelectedTab;
            if (selectedTab == app.um_linear_motor)
                app.um_vars = {'um_linear_motor_Fmax','um_linear_motor_Vmax','um_linear_motor_range_of_motion','um_linear_motor_voltage_frac'};
            elseif (selectedTab == app.um_hill_muscle_motor)
                app.um_vars = {'um_hill_motor_Fmax','um_hill_motor_Vmax','um_hill_motor_muscle_length','um_hill_motor_rate_of_activation','um_hill_motor_L_i','um_hill_motor_a_L','um_hill_motor_b_L','um_hill_motor_s'};
            end
            update_dd_vars(app);
        end

        % Button pushed function: ShowModelSchematicButton
        function ShowModelSchematicButtonPushed(app, event)
%             image = imread('model.pdf');
%             figure
%             imshow(image,'Border','tight','InitialMagnification','fit')
              winopen("model.pdf");
        end

        % Selection change function: graphing_corner
        function graphing_cornerSelectionChanged(app, event)
            selectedTab = app.graphing_corner.SelectedTab;
            if (selectedTab == app.graphing_corner_sensitivity)
                set(app.savesolutionCheckBox,'visible','off');
                app.go.Text = 'Compute!';
            else
                set(app.savesolutionCheckBox,'visible','on');
                app.go.Text = 'Graph!';
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.9373 0.9882 0.8627];
            app.UIFigure.Position = [60 -10 1132 677];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Scrollable = 'on';

            % Create ParametersLabel
            app.ParametersLabel = uilabel(app.UIFigure);
            app.ParametersLabel.FontSize = 22;
            app.ParametersLabel.FontWeight = 'bold';
            app.ParametersLabel.Position = [29 628 125 27];
            app.ParametersLabel.Text = 'Parameters';

            % Create spring
            app.spring = uitabgroup(app.UIFigure);
            app.spring.AutoResizeChildren = 'off';
            app.spring.SelectionChangedFcn = createCallbackFcn(app, @springSelectionChanged, true);
            app.spring.Position = [29 307 712 81];

            % Create linear_spring
            app.linear_spring = uitab(app.spring);
            app.linear_spring.AutoResizeChildren = 'off';
            app.linear_spring.Title = 'Linear Spring';
            app.linear_spring.BackgroundColor = [0.7686 0.9608 0.9608];

            % Create kLabel
            app.kLabel = uilabel(app.linear_spring);
            app.kLabel.HorizontalAlignment = 'right';
            app.kLabel.Position = [122 21 44 22];
            app.kLabel.Text = 'k =';

            % Create linear_spring_k
            app.linear_spring_k = uieditfield(app.linear_spring, 'numeric');
            app.linear_spring_k.Limits = [0 Inf];
            app.linear_spring_k.Tooltip = {'stiffness of the spring'};
            app.linear_spring_k.Position = [173 21 49 22];
            app.linear_spring_k.Value = 0.5;

            % Create massEditFieldLabel_4
            app.massEditFieldLabel_4 = uilabel(app.linear_spring);
            app.massEditFieldLabel_4.HorizontalAlignment = 'right';
            app.massEditFieldLabel_4.Position = [14 21 44 22];
            app.massEditFieldLabel_4.Text = 'mass =';

            % Create linear_spring_mass
            app.linear_spring_mass = uieditfield(app.linear_spring, 'numeric');
            app.linear_spring_mass.Limits = [0 Inf];
            app.linear_spring_mass.Tooltip = {'mass of the spring'};
            app.linear_spring_mass.Position = [65 21 49 22];

            % Create FmaxLabel_3
            app.FmaxLabel_3 = uilabel(app.linear_spring);
            app.FmaxLabel_3.HorizontalAlignment = 'right';
            app.FmaxLabel_3.Position = [252 21 46 22];
            app.FmaxLabel_3.Text = 'Fmax =';

            % Create linear_spring_Fmax
            app.linear_spring_Fmax = uieditfield(app.linear_spring, 'numeric');
            app.linear_spring_Fmax.Limits = [0 Inf];
            app.linear_spring_Fmax.Tooltip = {'max force of spring'};
            app.linear_spring_Fmax.Position = [307 21 49 22];
            app.linear_spring_Fmax.Value = Inf;

            % Create linear_elastic_extensional_spring
            app.linear_elastic_extensional_spring = uitab(app.spring);
            app.linear_elastic_extensional_spring.Title = 'Linear Elastic Extensional Spring';
            app.linear_elastic_extensional_spring.BackgroundColor = [0.4941 0.8784 0.8784];

            % Create ELabel
            app.ELabel = uilabel(app.linear_elastic_extensional_spring);
            app.ELabel.HorizontalAlignment = 'right';
            app.ELabel.Position = [14 21 34 22];
            app.ELabel.Text = 'E =';

            % Create lee_spring_E
            app.lee_spring_E = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_E.Limits = [0 Inf];
            app.lee_spring_E.Tooltip = {'modulus'};
            app.lee_spring_E.Position = [57 21 51 22];
            app.lee_spring_E.Value = 0.5;

            % Create rhoLabel
            app.rhoLabel = uilabel(app.linear_elastic_extensional_spring);
            app.rhoLabel.HorizontalAlignment = 'right';
            app.rhoLabel.Position = [351 22 33 22];
            app.rhoLabel.Text = 'rho =';

            % Create lee_spring_rho
            app.lee_spring_rho = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_rho.Limits = [0 Inf];
            app.lee_spring_rho.Tooltip = {'density'};
            app.lee_spring_rho.Position = [395 22 49 22];
            app.lee_spring_rho.Value = 10;

            % Create ALabel
            app.ALabel = uilabel(app.linear_elastic_extensional_spring);
            app.ALabel.HorizontalAlignment = 'right';
            app.ALabel.Position = [132 21 33 22];
            app.ALabel.Text = 'A =';

            % Create lee_spring_A
            app.lee_spring_A = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_A.Limits = [0 Inf];
            app.lee_spring_A.Tooltip = {'cross sectional area'};
            app.lee_spring_A.Position = [176 21 49 22];
            app.lee_spring_A.Value = 0.5;

            % Create LLabel
            app.LLabel = uilabel(app.linear_elastic_extensional_spring);
            app.LLabel.HorizontalAlignment = 'right';
            app.LLabel.Position = [236 22 33 22];
            app.LLabel.Text = 'L =';

            % Create lee_spring_L
            app.lee_spring_L = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_L.Limits = [0 Inf];
            app.lee_spring_L.Tooltip = {'length'};
            app.lee_spring_L.Position = [280 22 49 22];
            app.lee_spring_L.Value = 0.001;

            % Create sigma_fLabel
            app.sigma_fLabel = uilabel(app.linear_elastic_extensional_spring);
            app.sigma_fLabel.HorizontalAlignment = 'right';
            app.sigma_fLabel.Position = [456 22 71 22];
            app.sigma_fLabel.Text = 'sigma_f = ';

            % Create lee_spring_sigma_f
            app.lee_spring_sigma_f = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_sigma_f.Limits = [0 Inf];
            app.lee_spring_sigma_f.Tooltip = {'failure strength in Pa'};
            app.lee_spring_sigma_f.Position = [533 22 50 22];
            app.lee_spring_sigma_f.Value = Inf;

            % Create exponential_spring
            app.exponential_spring = uitab(app.spring);
            app.exponential_spring.AutoResizeChildren = 'off';
            app.exponential_spring.Title = 'Exponential Spring';
            app.exponential_spring.BackgroundColor = [0.3765 0.749 0.749];

            % Create k_0Label
            app.k_0Label = uilabel(app.exponential_spring);
            app.k_0Label.HorizontalAlignment = 'right';
            app.k_0Label.Position = [18 21 44 22];
            app.k_0Label.Text = 'k_0 =';

            % Create exp_spring_k
            app.exp_spring_k = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_k.Limits = [0 Inf];
            app.exp_spring_k.Tooltip = {'spring stiffness'};
            app.exp_spring_k.Position = [72 21 51 22];
            app.exp_spring_k.Value = 2000;

            % Create characteristiclengthLabel
            app.characteristiclengthLabel = uilabel(app.exponential_spring);
            app.characteristiclengthLabel.HorizontalAlignment = 'right';
            app.characteristiclengthLabel.Position = [277 22 123 22];
            app.characteristiclengthLabel.Text = 'characteristic length =';

            % Create exp_spring_char_len
            app.exp_spring_char_len = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_char_len.Limits = [0 Inf];
            app.exp_spring_char_len.Tooltip = {'resting length where force is 0'};
            app.exp_spring_char_len.Position = [409 22 51 22];
            app.exp_spring_char_len.Value = 0.001;

            % Create FmaxLabel_4
            app.FmaxLabel_4 = uilabel(app.exponential_spring);
            app.FmaxLabel_4.HorizontalAlignment = 'right';
            app.FmaxLabel_4.Position = [487 22 46 22];
            app.FmaxLabel_4.Text = 'Fmax =';

            % Create exp_spring_Fmax
            app.exp_spring_Fmax = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_Fmax.Limits = [0 Inf];
            app.exp_spring_Fmax.Tooltip = {'max force of the spring'};
            app.exp_spring_Fmax.Position = [540 22 51 22];
            app.exp_spring_Fmax.Value = 20;

            % Create massEditFieldLabel_5
            app.massEditFieldLabel_5 = uilabel(app.exponential_spring);
            app.massEditFieldLabel_5.HorizontalAlignment = 'right';
            app.massEditFieldLabel_5.Position = [147 21 44 22];
            app.massEditFieldLabel_5.Text = 'mass =';

            % Create exp_spring_mass
            app.exp_spring_mass = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_mass.Limits = [0 Inf];
            app.exp_spring_mass.Tooltip = {'mass of the spring'};
            app.exp_spring_mass.Position = [204 21 51 22];
            app.exp_spring_mass.Value = 0.002;

            % Create loading_motor
            app.loading_motor = uitabgroup(app.UIFigure);
            app.loading_motor.AutoResizeChildren = 'off';
            app.loading_motor.SelectionChangedFcn = createCallbackFcn(app, @loading_motorSelectionChanged, true);
            app.loading_motor.Position = [30 168 711 110];

            % Create lm_linear_motor
            app.lm_linear_motor = uitab(app.loading_motor);
            app.lm_linear_motor.AutoResizeChildren = 'off';
            app.lm_linear_motor.Title = 'Linear Motor';
            app.lm_linear_motor.BackgroundColor = [1 0.7686 0.7686];

            % Create FmaxLabel_2
            app.FmaxLabel_2 = uilabel(app.lm_linear_motor);
            app.FmaxLabel_2.HorizontalAlignment = 'right';
            app.FmaxLabel_2.Position = [12 52 46 22];
            app.FmaxLabel_2.Text = 'Fmax =';

            % Create lm_linear_motor_Fmax
            app.lm_linear_motor_Fmax = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_Fmax.Limits = [0 Inf];
            app.lm_linear_motor_Fmax.Tooltip = {'max force of the motor'};
            app.lm_linear_motor_Fmax.Position = [68 52 49 22];
            app.lm_linear_motor_Fmax.Value = 10;

            % Create VmaxLabel_2
            app.VmaxLabel_2 = uilabel(app.lm_linear_motor);
            app.VmaxLabel_2.HorizontalAlignment = 'right';
            app.VmaxLabel_2.Position = [137 52 47 22];
            app.VmaxLabel_2.Text = 'Vmax =';

            % Create lm_linear_motor_Vmax
            app.lm_linear_motor_Vmax = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_Vmax.Limits = [0 Inf];
            app.lm_linear_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.lm_linear_motor_Vmax.Position = [195 52 49 22];
            app.lm_linear_motor_Vmax.Value = 10;

            % Create rangeofmotionLabel
            app.rangeofmotionLabel = uilabel(app.lm_linear_motor);
            app.rangeofmotionLabel.HorizontalAlignment = 'right';
            app.rangeofmotionLabel.Position = [433 52 99 22];
            app.rangeofmotionLabel.Text = 'range of motion =';

            % Create lm_linear_motor_range_of_motion
            app.lm_linear_motor_range_of_motion = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_range_of_motion.Limits = [0 Inf];
            app.lm_linear_motor_range_of_motion.Tooltip = {'total possible range of motion of the motor'};
            app.lm_linear_motor_range_of_motion.Position = [545 52 48 22];
            app.lm_linear_motor_range_of_motion.Value = 0.005;

            % Create voltagefracLabel
            app.voltagefracLabel = uilabel(app.lm_linear_motor);
            app.voltagefracLabel.HorizontalAlignment = 'right';
            app.voltagefracLabel.Position = [270 52 78 22];
            app.voltagefracLabel.Text = 'voltage frac =';

            % Create lm_linear_motor_voltage_frac
            app.lm_linear_motor_voltage_frac = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_voltage_frac.Limits = [0 Inf];
            app.lm_linear_motor_voltage_frac.Tooltip = {'fraction of max voltage of motor'};
            app.lm_linear_motor_voltage_frac.Position = [359 52 48 22];
            app.lm_linear_motor_voltage_frac.Value = 1;

            % Create lm_hill_muscle_motor
            app.lm_hill_muscle_motor = uitab(app.loading_motor);
            app.lm_hill_muscle_motor.AutoResizeChildren = 'off';
            app.lm_hill_muscle_motor.Title = 'Hill Muscle Motor';
            app.lm_hill_muscle_motor.BackgroundColor = [1 0.5686 0.5686];

            % Create FmaxLabel
            app.FmaxLabel = uilabel(app.lm_hill_muscle_motor);
            app.FmaxLabel.HorizontalAlignment = 'right';
            app.FmaxLabel.Position = [43 52 46 22];
            app.FmaxLabel.Text = 'Fmax =';

            % Create lm_hill_motor_Fmax
            app.lm_hill_motor_Fmax = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_Fmax.Limits = [0 Inf];
            app.lm_hill_motor_Fmax.Tooltip = {'max force of the motor'};
            app.lm_hill_motor_Fmax.Position = [101 52 49 22];
            app.lm_hill_motor_Fmax.Value = 4;

            % Create VmaxLabel
            app.VmaxLabel = uilabel(app.lm_hill_muscle_motor);
            app.VmaxLabel.HorizontalAlignment = 'right';
            app.VmaxLabel.Position = [177 52 47 22];
            app.VmaxLabel.Text = 'Vmax =';

            % Create lm_hill_motor_Vmax
            app.lm_hill_motor_Vmax = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_Vmax.Limits = [0 Inf];
            app.lm_hill_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.lm_hill_motor_Vmax.Position = [236 52 49 22];
            app.lm_hill_motor_Vmax.Value = 1;

            % Create musclelengthLabel
            app.musclelengthLabel = uilabel(app.lm_hill_muscle_motor);
            app.musclelengthLabel.HorizontalAlignment = 'right';
            app.musclelengthLabel.Position = [311 52 90 22];
            app.musclelengthLabel.Text = 'muscle length =';

            % Create lm_hill_motor_muscle_length
            app.lm_hill_motor_muscle_length = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_muscle_length.Limits = [0 Inf];
            app.lm_hill_motor_muscle_length.ValueChangedFcn = createCallbackFcn(app, @lm_hill_motor_muscle_lengthValueChanged, true);
            app.lm_hill_motor_muscle_length.Tooltip = {''};
            app.lm_hill_motor_muscle_length.Position = [416 52 47 22];
            app.lm_hill_motor_muscle_length.Value = 0.01;

            % Create rateofactivationLabel
            app.rateofactivationLabel = uilabel(app.lm_hill_muscle_motor);
            app.rateofactivationLabel.HorizontalAlignment = 'right';
            app.rateofactivationLabel.Position = [486 51 104 22];
            app.rateofactivationLabel.Text = 'rate of activation =';

            % Create lm_hill_motor_rate_of_activation
            app.lm_hill_motor_rate_of_activation = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_rate_of_activation.Limits = [0 Inf];
            app.lm_hill_motor_rate_of_activation.Tooltip = {''};
            app.lm_hill_motor_rate_of_activation.Position = [605 51 48 22];
            app.lm_hill_motor_rate_of_activation.Value = Inf;

            % Create initiallengthLabel
            app.initiallengthLabel = uilabel(app.lm_hill_muscle_motor);
            app.initiallengthLabel.HorizontalAlignment = 'right';
            app.initiallengthLabel.Position = [10 18 79 22];
            app.initiallengthLabel.Text = 'initial length =';

            % Create lm_hill_motor_L_i
            app.lm_hill_motor_L_i = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_L_i.Limits = [0 Inf];
            app.lm_hill_motor_L_i.Tooltip = {''};
            app.lm_hill_motor_L_i.Position = [101 18 49 22];
            app.lm_hill_motor_L_i.Value = 0.01;

            % Create a_LLabel
            app.a_LLabel = uilabel(app.lm_hill_muscle_motor);
            app.a_LLabel.HorizontalAlignment = 'right';
            app.a_LLabel.Position = [188 18 36 22];
            app.a_LLabel.Text = 'a_L =';

            % Create lm_hill_motor_a_L
            app.lm_hill_motor_a_L = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_a_L.Tooltip = {'determines shape of length-tension relationship'};
            app.lm_hill_motor_a_L.Position = [239 18 46 22];
            app.lm_hill_motor_a_L.Value = 2.08;

            % Create b_LLabel
            app.b_LLabel = uilabel(app.lm_hill_muscle_motor);
            app.b_LLabel.HorizontalAlignment = 'right';
            app.b_LLabel.Position = [365 17 36 22];
            app.b_LLabel.Text = 'b_L =';

            % Create lm_hill_motor_b_L
            app.lm_hill_motor_b_L = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_b_L.Tooltip = {'determines shape of length-tension relationship'};
            app.lm_hill_motor_b_L.Position = [415 17 49 22];
            app.lm_hill_motor_b_L.Value = -2.89;

            % Create sLabel
            app.sLabel = uilabel(app.lm_hill_muscle_motor);
            app.sLabel.HorizontalAlignment = 'right';
            app.sLabel.Position = [566 17 25 22];
            app.sLabel.Text = 's =';

            % Create lm_hill_motor_s
            app.lm_hill_motor_s = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_s.Tooltip = {'determines shape of length-tension relationship'};
            app.lm_hill_motor_s.Position = [605 17 48 22];
            app.lm_hill_motor_s.Value = -0.75;

            % Create latch
            app.latch = uitabgroup(app.UIFigure);
            app.latch.AutoResizeChildren = 'off';
            app.latch.Position = [29 404 712 114];

            % Create rounded_latch
            app.rounded_latch = uitab(app.latch);
            app.rounded_latch.AutoResizeChildren = 'off';
            app.rounded_latch.Title = 'Rounded Latch';
            app.rounded_latch.BackgroundColor = [0.7137 0.9804 0.7137];

            % Create massLabel
            app.massLabel = uilabel(app.rounded_latch);
            app.massLabel.HorizontalAlignment = 'right';
            app.massLabel.Position = [59 53 54 22];
            app.massLabel.Text = '   mass =';

            % Create latch_mass
            app.latch_mass = uieditfield(app.rounded_latch, 'numeric');
            app.latch_mass.Limits = [0 Inf];
            app.latch_mass.Tooltip = {'mass of the latch'};
            app.latch_mass.Position = [122 53 51 22];
            app.latch_mass.Value = 0.003;

            % Create Label
            app.Label = uilabel(app.rounded_latch);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [441 55 25 22];
            app.Label.Text = 'μ =';

            % Create latch_coeff_fric
            app.latch_coeff_fric = uieditfield(app.rounded_latch, 'numeric');
            app.latch_coeff_fric.Limits = [0 Inf];
            app.latch_coeff_fric.Tooltip = {'coefficient of friction between the latch and load mass'};
            app.latch_coeff_fric.Position = [476 55 51 22];

            % Create radiusLabel
            app.radiusLabel = uilabel(app.rounded_latch);
            app.radiusLabel.HorizontalAlignment = 'right';
            app.radiusLabel.Position = [249 53 48 22];
            app.radiusLabel.Text = 'radius =';

            % Create latch_radius
            app.latch_radius = uieditfield(app.rounded_latch, 'numeric');
            app.latch_radius.Limits = [0 Inf];
            app.latch_radius.Tooltip = {'latch radius'};
            app.latch_radius.Position = [309 53 48 22];
            app.latch_radius.Value = 0.005;

            % Create v_0Label
            app.v_0Label = uilabel(app.rounded_latch);
            app.v_0Label.HorizontalAlignment = 'right';
            app.v_0Label.Position = [592 53 38 22];
            app.v_0Label.Text = 'v_0 =';

            % Create latch_v_0
            app.latch_v_0 = uieditfield(app.rounded_latch, 'numeric');
            app.latch_v_0.Limits = [0 Inf];
            app.latch_v_0.Tooltip = {'initial velocity of the latch'};
            app.latch_v_0.Position = [642 53 48 22];

            % Create minlatchingdistLabel
            app.minlatchingdistLabel = uilabel(app.rounded_latch);
            app.minlatchingdistLabel.HorizontalAlignment = 'right';
            app.minlatchingdistLabel.Position = [11 19 102 22];
            app.minlatchingdistLabel.Text = 'min latching dist =';

            % Create min_latching_dist
            app.min_latching_dist = uieditfield(app.rounded_latch, 'numeric');
            app.min_latching_dist.Limits = [0 Inf];
            app.min_latching_dist.Tooltip = {'lowest position the latch can engage'};
            app.min_latching_dist.Position = [122 19 51 22];

            % Create maxlatchingdistLabel
            app.maxlatchingdistLabel = uilabel(app.rounded_latch);
            app.maxlatchingdistLabel.HorizontalAlignment = 'right';
            app.maxlatchingdistLabel.Position = [192 19 105 22];
            app.maxlatchingdistLabel.Text = 'max latching dist =';

            % Create max_latching_dist
            app.max_latching_dist = uieditfield(app.rounded_latch, 'numeric');
            app.max_latching_dist.Limits = [0 Inf];
            app.max_latching_dist.Tooltip = {'highest position the latch can engage'};
            app.max_latching_dist.Position = [309 19 48 22];
            app.max_latching_dist.Value = Inf;

            % Create runwaylengthLabel
            app.runwaylengthLabel = uilabel(app.rounded_latch);
            app.runwaylengthLabel.HorizontalAlignment = 'right';
            app.runwaylengthLabel.Position = [374 19 91 22];
            app.runwaylengthLabel.Text = 'runway length =';

            % Create runway_length
            app.runway_length = uieditfield(app.rounded_latch, 'numeric');
            app.runway_length.Limits = [0 Inf];
            app.runway_length.Tooltip = {'coefficient of friction between the latch and load mass'};
            app.runway_length.Position = [475 19 52 22];

            % Create load
            app.load = uitabgroup(app.UIFigure);
            app.load.AutoResizeChildren = 'off';
            app.load.Position = [30 533 711 79];

            % Create load_mass
            app.load_mass = uitab(app.load);
            app.load_mass.AutoResizeChildren = 'off';
            app.load_mass.Title = 'Load Mass';
            app.load_mass.BackgroundColor = [0.9137 0.8039 0.9804];

            % Create massEditFieldLabel
            app.massEditFieldLabel = uilabel(app.load_mass);
            app.massEditFieldLabel.HorizontalAlignment = 'right';
            app.massEditFieldLabel.Position = [13 18 44 22];
            app.massEditFieldLabel.Text = 'mass =';

            % Create load_mass_mass
            app.load_mass_mass = uieditfield(app.load_mass, 'numeric');
            app.load_mass_mass.Limits = [0 Inf];
            app.load_mass_mass.Tooltip = {'mass of the load'};
            app.load_mass_mass.Position = [72 18 52 22];
            app.load_mass_mass.Value = 0.01;

            % Create massofleverarmLabel
            app.massofleverarmLabel = uilabel(app.load_mass);
            app.massofleverarmLabel.HorizontalAlignment = 'right';
            app.massofleverarmLabel.Position = [132 18 111 22];
            app.massofleverarmLabel.Text = 'mass of lever arm =';

            % Create load_m_rod
            app.load_m_rod = uieditfield(app.load_mass, 'numeric');
            app.load_m_rod.Limits = [0 Inf];
            app.load_m_rod.Tooltip = {''};
            app.load_m_rod.Position = [258 18 46 22];

            % Create EMALabel
            app.EMALabel = uilabel(app.load_mass);
            app.EMALabel.HorizontalAlignment = 'right';
            app.EMALabel.Position = [320 18 42 22];
            app.EMALabel.Text = 'EMA =';

            % Create load_EMA
            app.load_EMA = uieditfield(app.load_mass, 'numeric');
            app.load_EMA.Limits = [0.001 Inf];
            app.load_EMA.Tooltip = {'effective mechanical advantage'};
            app.load_EMA.Position = [374 18 48 22];
            app.load_EMA.Value = 1;

            % Create LoadingMotorLabel
            app.LoadingMotorLabel = uilabel(app.UIFigure);
            app.LoadingMotorLabel.FontWeight = 'bold';
            app.LoadingMotorLabel.Position = [32 277 92 22];
            app.LoadingMotorLabel.Text = 'Loading Motor:';

            % Create unlatching_motor
            app.unlatching_motor = uitabgroup(app.UIFigure);
            app.unlatching_motor.AutoResizeChildren = 'off';
            app.unlatching_motor.SelectionChangedFcn = createCallbackFcn(app, @unlatching_motorSelectionChanged, true);
            app.unlatching_motor.Position = [29 27 712 110];

            % Create um_linear_motor
            app.um_linear_motor = uitab(app.unlatching_motor);
            app.um_linear_motor.AutoResizeChildren = 'off';
            app.um_linear_motor.Title = 'Linear Motor';
            app.um_linear_motor.BackgroundColor = [1 0.7686 0.7686];

            % Create FmaxLabel_5
            app.FmaxLabel_5 = uilabel(app.um_linear_motor);
            app.FmaxLabel_5.HorizontalAlignment = 'right';
            app.FmaxLabel_5.Position = [14 51 46 22];
            app.FmaxLabel_5.Text = 'Fmax =';

            % Create um_linear_motor_Fmax
            app.um_linear_motor_Fmax = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_Fmax.Limits = [0 Inf];
            app.um_linear_motor_Fmax.Tooltip = {'max force of the motor'};
            app.um_linear_motor_Fmax.Position = [68 51 50 22];
            app.um_linear_motor_Fmax.Value = 0.25;

            % Create VmaxLabel_3
            app.VmaxLabel_3 = uilabel(app.um_linear_motor);
            app.VmaxLabel_3.HorizontalAlignment = 'right';
            app.VmaxLabel_3.Position = [140 51 47 22];
            app.VmaxLabel_3.Text = 'Vmax =';

            % Create um_linear_motor_Vmax
            app.um_linear_motor_Vmax = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_Vmax.Limits = [0 Inf];
            app.um_linear_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.um_linear_motor_Vmax.Position = [196 51 49 22];
            app.um_linear_motor_Vmax.Value = 1;

            % Create rangeofmotionLabel_2
            app.rangeofmotionLabel_2 = uilabel(app.um_linear_motor);
            app.rangeofmotionLabel_2.HorizontalAlignment = 'right';
            app.rangeofmotionLabel_2.Position = [433 52 99 22];
            app.rangeofmotionLabel_2.Text = 'range of motion =';

            % Create um_linear_motor_range_of_motion
            app.um_linear_motor_range_of_motion = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_range_of_motion.Limits = [0 Inf];
            app.um_linear_motor_range_of_motion.Tooltip = {'total possible range of motion of the motor'};
            app.um_linear_motor_range_of_motion.Position = [547 52 47 22];
            app.um_linear_motor_range_of_motion.Value = 0.005;

            % Create voltagefracLabel_2
            app.voltagefracLabel_2 = uilabel(app.um_linear_motor);
            app.voltagefracLabel_2.HorizontalAlignment = 'right';
            app.voltagefracLabel_2.Position = [272 51 78 22];
            app.voltagefracLabel_2.Text = 'voltage frac =';

            % Create um_linear_motor_voltage_frac
            app.um_linear_motor_voltage_frac = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_voltage_frac.Limits = [0 Inf];
            app.um_linear_motor_voltage_frac.Tooltip = {'fraction of max voltage of motor'};
            app.um_linear_motor_voltage_frac.Position = [359 51 49 22];
            app.um_linear_motor_voltage_frac.Value = 1;

            % Create um_hill_muscle_motor
            app.um_hill_muscle_motor = uitab(app.unlatching_motor);
            app.um_hill_muscle_motor.AutoResizeChildren = 'off';
            app.um_hill_muscle_motor.Title = 'Hill Muscle Motor';
            app.um_hill_muscle_motor.BackgroundColor = [1 0.5686 0.5686];

            % Create FmaxLabel_6
            app.FmaxLabel_6 = uilabel(app.um_hill_muscle_motor);
            app.FmaxLabel_6.HorizontalAlignment = 'right';
            app.FmaxLabel_6.Position = [41 52 46 22];
            app.FmaxLabel_6.Text = 'Fmax =';

            % Create um_hill_motor_Fmax
            app.um_hill_motor_Fmax = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_Fmax.Limits = [0 Inf];
            app.um_hill_motor_Fmax.Tooltip = {'max force of the motor'};
            app.um_hill_motor_Fmax.Position = [99 52 49 22];
            app.um_hill_motor_Fmax.Value = 10;

            % Create VmaxLabel_4
            app.VmaxLabel_4 = uilabel(app.um_hill_muscle_motor);
            app.VmaxLabel_4.HorizontalAlignment = 'right';
            app.VmaxLabel_4.Position = [180 52 47 22];
            app.VmaxLabel_4.Text = 'Vmax =';

            % Create um_hill_motor_Vmax
            app.um_hill_motor_Vmax = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_Vmax.Limits = [0 Inf];
            app.um_hill_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.um_hill_motor_Vmax.Position = [238 52 48 22];
            app.um_hill_motor_Vmax.Value = 10;

            % Create musclelengthLabel_2
            app.musclelengthLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.musclelengthLabel_2.HorizontalAlignment = 'right';
            app.musclelengthLabel_2.Position = [309 52 90 22];
            app.musclelengthLabel_2.Text = 'muscle length =';

            % Create um_hill_motor_muscle_length
            app.um_hill_motor_muscle_length = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_muscle_length.Limits = [0 Inf];
            app.um_hill_motor_muscle_length.ValueChangedFcn = createCallbackFcn(app, @um_hill_motor_muscle_lengthValueChanged, true);
            app.um_hill_motor_muscle_length.Tooltip = {''};
            app.um_hill_motor_muscle_length.Position = [409 52 50 22];
            app.um_hill_motor_muscle_length.Value = 4;

            % Create rateofactivationLabel_2
            app.rateofactivationLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.rateofactivationLabel_2.HorizontalAlignment = 'right';
            app.rateofactivationLabel_2.Position = [489 52 104 22];
            app.rateofactivationLabel_2.Text = 'rate of activation =';

            % Create um_hill_motor_rate_of_activation
            app.um_hill_motor_rate_of_activation = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_rate_of_activation.Limits = [0 Inf];
            app.um_hill_motor_rate_of_activation.Tooltip = {''};
            app.um_hill_motor_rate_of_activation.Position = [606 52 48 22];
            app.um_hill_motor_rate_of_activation.Value = 2;

            % Create initiallengthLabel_2
            app.initiallengthLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.initiallengthLabel_2.HorizontalAlignment = 'right';
            app.initiallengthLabel_2.Position = [9 19 79 22];
            app.initiallengthLabel_2.Text = 'initial length =';

            % Create um_hill_motor_L_i
            app.um_hill_motor_L_i = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_L_i.Limits = [0 Inf];
            app.um_hill_motor_L_i.Tooltip = {''};
            app.um_hill_motor_L_i.Position = [99 19 49 22];
            app.um_hill_motor_L_i.Value = 4;

            % Create a_LLabel_2
            app.a_LLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.a_LLabel_2.HorizontalAlignment = 'right';
            app.a_LLabel_2.Position = [191 19 36 22];
            app.a_LLabel_2.Text = 'a_L =';

            % Create um_hill_motor_a_L
            app.um_hill_motor_a_L = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_a_L.Tooltip = {'determines shape of length-tension relationship'};
            app.um_hill_motor_a_L.Position = [238 19 49 22];
            app.um_hill_motor_a_L.Value = 2.08;

            % Create b_LLabel_2
            app.b_LLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.b_LLabel_2.HorizontalAlignment = 'right';
            app.b_LLabel_2.Position = [363 19 36 22];
            app.b_LLabel_2.Text = 'b_L =';

            % Create um_hill_motor_b_L
            app.um_hill_motor_b_L = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_b_L.Tooltip = {'determines shape of length-tension relationship'};
            app.um_hill_motor_b_L.Position = [409 19 50 22];
            app.um_hill_motor_b_L.Value = -2.89;

            % Create sLabel_2
            app.sLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.sLabel_2.HorizontalAlignment = 'right';
            app.sLabel_2.Position = [568 19 25 22];
            app.sLabel_2.Text = 's =';

            % Create um_hill_motor_s
            app.um_hill_motor_s = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_s.Tooltip = {'determines shape of length-tension relationship'};
            app.um_hill_motor_s.Position = [606 19 48 22];
            app.um_hill_motor_s.Value = -0.75;

            % Create UnlatchingMotorLabel
            app.UnlatchingMotorLabel = uilabel(app.UIFigure);
            app.UnlatchingMotorLabel.FontWeight = 'bold';
            app.UnlatchingMotorLabel.Position = [31 138 108 22];
            app.UnlatchingMotorLabel.Text = 'Unlatching Motor:';

            % Create GraphingCornerLabel
            app.GraphingCornerLabel = uilabel(app.UIFigure);
            app.GraphingCornerLabel.FontSize = 22;
            app.GraphingCornerLabel.FontWeight = 'bold';
            app.GraphingCornerLabel.Position = [854 627 181 27];
            app.GraphingCornerLabel.Text = 'Graphing Corner';

            % Create graphing_corner
            app.graphing_corner = uitabgroup(app.UIFigure);
            app.graphing_corner.AutoResizeChildren = 'off';
            app.graphing_corner.SelectionChangedFcn = createCallbackFcn(app, @graphing_cornerSelectionChanged, true);
            app.graphing_corner.Position = [783 218 322 390];

            % Create graphing_corner_kinematics
            app.graphing_corner_kinematics = uitab(app.graphing_corner);
            app.graphing_corner_kinematics.Title = 'Kinematics';
            app.graphing_corner_kinematics.BackgroundColor = [0.949 0.8471 0.6078];

            % Create KinematicsOutputOptionsLabel
            app.KinematicsOutputOptionsLabel = uilabel(app.graphing_corner_kinematics);
            app.KinematicsOutputOptionsLabel.FontSize = 22;
            app.KinematicsOutputOptionsLabel.FontWeight = 'bold';
            app.KinematicsOutputOptionsLabel.Position = [16 324 294 28];
            app.KinematicsOutputOptionsLabel.Text = 'Kinematics Output Options';

            % Create forcedisp
            app.forcedisp = uicheckbox(app.graphing_corner_kinematics);
            app.forcedisp.Text = 'motor and spring force displacement curves';
            app.forcedisp.Position = [54 265 256 22];
            app.forcedisp.Value = true;

            % Create latchkinematicsCheckBox
            app.latchkinematicsCheckBox = uicheckbox(app.graphing_corner_kinematics);
            app.latchkinematicsCheckBox.Text = 'latch kinematics';
            app.latchkinematicsCheckBox.Position = [54 210 107 22];
            app.latchkinematicsCheckBox.Value = true;

            % Create loadkinematicsCheckBox
            app.loadkinematicsCheckBox = uicheckbox(app.graphing_corner_kinematics);
            app.loadkinematicsCheckBox.Text = 'load kinematics';
            app.loadkinematicsCheckBox.Position = [54 183 105 22];
            app.loadkinematicsCheckBox.Value = true;

            % Create LoadingphaseLabel
            app.LoadingphaseLabel = uilabel(app.graphing_corner_kinematics);
            app.LoadingphaseLabel.Position = [34 290 111 22];
            app.LoadingphaseLabel.Text = 'Loading phase:';

            % Create UnlatchingandlaunchingphaseLabel
            app.UnlatchingandlaunchingphaseLabel = uilabel(app.graphing_corner_kinematics);
            app.UnlatchingandlaunchingphaseLabel.Position = [33 238 180 22];
            app.UnlatchingandlaunchingphaseLabel.Text = 'Unlatching and launching phase:';

            % Create graphing_corner_one_D
            app.graphing_corner_one_D = uitab(app.graphing_corner);
            app.graphing_corner_one_D.AutoResizeChildren = 'off';
            app.graphing_corner_one_D.Title = '1D Plot';
            app.graphing_corner_one_D.BackgroundColor = [0.9608 0.7804 0.3647];

            % Create OD_minunlatchingmotorforceCheckBox
            app.OD_minunlatchingmotorforceCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_minunlatchingmotorforceCheckBox.Text = 'min unlatching motor force';
            app.OD_minunlatchingmotorforceCheckBox.Position = [127 173 164 22];

            % Create xaxisLabel_2
            app.xaxisLabel_2 = uilabel(app.graphing_corner_one_D);
            app.xaxisLabel_2.FontSize = 22;
            app.xaxisLabel_2.FontWeight = 'bold';
            app.xaxisLabel_2.Position = [50 329 67 28];
            app.xaxisLabel_2.Text = 'x-axis';

            % Create OD_y_maxCheckBox
            app.OD_y_maxCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_y_maxCheckBox.Text = 'y_max';
            app.OD_y_maxCheckBox.Position = [32 221 58 22];

            % Create OD_y_unlatchCheckBox
            app.OD_y_unlatchCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_y_unlatchCheckBox.Text = 'y_unlatch';
            app.OD_y_unlatchCheckBox.Position = [32 197 74 22];

            % Create OD_t_LCheckBox
            app.OD_t_LCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_t_LCheckBox.Text = 't_L';
            app.OD_t_LCheckBox.Position = [32 173 39 22];

            % Create OD_v_toCheckBox
            app.OD_v_toCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_v_toCheckBox.Text = 'v_to';
            app.OD_v_toCheckBox.Position = [32 149 45 22];

            % Create OD_P_maxCheckBox
            app.OD_P_maxCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_P_maxCheckBox.Text = 'P_max';
            app.OD_P_maxCheckBox.Position = [127 221 60 22];
            app.OD_P_maxCheckBox.Value = true;

            % Create OD_t_toCheckBox
            app.OD_t_toCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_t_toCheckBox.Text = 't_to';
            app.OD_t_toCheckBox.Position = [127 197 42 22];

            % Create OD_KE_maxCheckBox
            app.OD_KE_maxCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_KE_maxCheckBox.Text = 'KE_max';
            app.OD_KE_maxCheckBox.Position = [127 149 68 22];

            % Create OD_IV1DropDown
            app.OD_IV1DropDown = uidropdown(app.graphing_corner_one_D);
            app.OD_IV1DropDown.Items = {'--no selection--'};
            app.OD_IV1DropDown.Tooltip = {'Independent Variable 1'};
            app.OD_IV1DropDown.FontSize = 16;
            app.OD_IV1DropDown.Position = [128 329 150 22];
            app.OD_IV1DropDown.Value = '--no selection--';

            % Create yaxisOutputOptionsLabel
            app.yaxisOutputOptionsLabel = uilabel(app.graphing_corner_one_D);
            app.yaxisOutputOptionsLabel.FontSize = 22;
            app.yaxisOutputOptionsLabel.FontWeight = 'bold';
            app.yaxisOutputOptionsLabel.Position = [44 257 237 28];
            app.yaxisOutputOptionsLabel.Text = 'y-axis Output Options';

            % Create pixelsofresolutionLabel_2
            app.pixelsofresolutionLabel_2 = uilabel(app.graphing_corner_one_D);
            app.pixelsofresolutionLabel_2.HorizontalAlignment = 'right';
            app.pixelsofresolutionLabel_2.FontSize = 16;
            app.pixelsofresolutionLabel_2.Position = [64 75 151 22];
            app.pixelsofresolutionLabel_2.Text = 'pixels of resolution =';

            % Create OD_n
            app.OD_n = uieditfield(app.graphing_corner_one_D, 'numeric');
            app.OD_n.Limits = [1 1024];
            app.OD_n.RoundFractionalValues = 'on';
            app.OD_n.ValueDisplayFormat = '%.0f';
            app.OD_n.Position = [225 77 39 22];
            app.OD_n.Value = 100;

            % Create xminEditFieldLabel_2
            app.xminEditFieldLabel_2 = uilabel(app.graphing_corner_one_D);
            app.xminEditFieldLabel_2.HorizontalAlignment = 'right';
            app.xminEditFieldLabel_2.Position = [19 298 34 22];
            app.xminEditFieldLabel_2.Text = 'x min';

            % Create OD_xmin
            app.OD_xmin = uieditfield(app.graphing_corner_one_D, 'numeric');
            app.OD_xmin.Position = [68 298 40 22];

            % Create xmaxEditFieldLabel_2
            app.xmaxEditFieldLabel_2 = uilabel(app.graphing_corner_one_D);
            app.xmaxEditFieldLabel_2.HorizontalAlignment = 'right';
            app.xmaxEditFieldLabel_2.Position = [116 298 38 22];
            app.xmaxEditFieldLabel_2.Text = 'x max';

            % Create OD_xmax
            app.OD_xmax = uieditfield(app.graphing_corner_one_D, 'numeric');
            app.OD_xmax.Position = [169 298 40 22];
            app.OD_xmax.Value = 42;

            % Create OD_x_log_space
            app.OD_x_log_space = uiswitch(app.graphing_corner_one_D, 'slider');
            app.OD_x_log_space.Items = {'lin', 'log'};
            app.OD_x_log_space.Position = [251 303 27 12];
            app.OD_x_log_space.Value = 'lin';

            % Create graphing_corner_heatmap
            app.graphing_corner_heatmap = uitab(app.graphing_corner);
            app.graphing_corner_heatmap.AutoResizeChildren = 'off';
            app.graphing_corner_heatmap.Title = '2D Plot';
            app.graphing_corner_heatmap.BackgroundColor = [0.9804 0.6941 0.0275];

            % Create minunlatchingmotorforceCheckBox
            app.minunlatchingmotorforceCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.minunlatchingmotorforceCheckBox.Text = 'min unlatching motor force';
            app.minunlatchingmotorforceCheckBox.Position = [127 100 164 22];

            % Create xaxisLabel
            app.xaxisLabel = uilabel(app.graphing_corner_heatmap);
            app.xaxisLabel.FontSize = 22;
            app.xaxisLabel.FontWeight = 'bold';
            app.xaxisLabel.Position = [21 329 67 28];
            app.xaxisLabel.Text = 'x-axis';

            % Create yaxisLabel
            app.yaxisLabel = uilabel(app.graphing_corner_heatmap);
            app.yaxisLabel.FontSize = 22;
            app.yaxisLabel.FontWeight = 'bold';
            app.yaxisLabel.Position = [21 252 67 28];
            app.yaxisLabel.Text = 'y-axis';

            % Create y_maxCheckBox
            app.y_maxCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.y_maxCheckBox.Text = 'y_max';
            app.y_maxCheckBox.Position = [32 148 58 22];

            % Create y_unlatchCheckBox
            app.y_unlatchCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.y_unlatchCheckBox.Text = 'y_unlatch';
            app.y_unlatchCheckBox.Position = [32 124 74 22];

            % Create t_LCheckBox
            app.t_LCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.t_LCheckBox.Text = 't_L';
            app.t_LCheckBox.Position = [32 100 39 22];

            % Create v_toCheckBox
            app.v_toCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.v_toCheckBox.Text = 'v_to';
            app.v_toCheckBox.Position = [32 76 45 22];

            % Create P_maxCheckBox
            app.P_maxCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.P_maxCheckBox.Text = 'P_max';
            app.P_maxCheckBox.Position = [127 148 60 22];
            app.P_maxCheckBox.Value = true;

            % Create t_toCheckBox
            app.t_toCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.t_toCheckBox.Text = 't_to';
            app.t_toCheckBox.Position = [127 124 42 22];

            % Create KE_maxCheckBox
            app.KE_maxCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.KE_maxCheckBox.Text = 'KE_max';
            app.KE_maxCheckBox.Position = [127 76 68 22];

            % Create IV1DropDown
            app.IV1DropDown = uidropdown(app.graphing_corner_heatmap);
            app.IV1DropDown.Items = {'--no selection--'};
            app.IV1DropDown.ValueChangedFcn = createCallbackFcn(app, @IV1DropDownValueChanged, true);
            app.IV1DropDown.Tooltip = {'Independent Variable 1'};
            app.IV1DropDown.FontSize = 16;
            app.IV1DropDown.Position = [110 329 150 22];
            app.IV1DropDown.Value = '--no selection--';

            % Create IV2DropDown
            app.IV2DropDown = uidropdown(app.graphing_corner_heatmap);
            app.IV2DropDown.Items = {'--no selection--'};
            app.IV2DropDown.ValueChangedFcn = createCallbackFcn(app, @IV2DropDownValueChanged, true);
            app.IV2DropDown.FontSize = 16;
            app.IV2DropDown.Position = [110 252 150 22];
            app.IV2DropDown.Value = '--no selection--';

            % Create xminEditFieldLabel
            app.xminEditFieldLabel = uilabel(app.graphing_corner_heatmap);
            app.xminEditFieldLabel.HorizontalAlignment = 'right';
            app.xminEditFieldLabel.Position = [19 294 34 22];
            app.xminEditFieldLabel.Text = 'x min';

            % Create xmin
            app.xmin = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.xmin.Position = [68 294 43 22];
            app.xmin.Value = 500;

            % Create xmaxEditFieldLabel
            app.xmaxEditFieldLabel = uilabel(app.graphing_corner_heatmap);
            app.xmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.xmaxEditFieldLabel.Position = [123 294 38 22];
            app.xmaxEditFieldLabel.Text = 'x max';

            % Create xmax
            app.xmax = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.xmax.Position = [174 294 42 22];
            app.xmax.Value = 8000;

            % Create yminLabel
            app.yminLabel = uilabel(app.graphing_corner_heatmap);
            app.yminLabel.HorizontalAlignment = 'right';
            app.yminLabel.Position = [19 217 34 22];
            app.yminLabel.Text = 'y min';

            % Create ymin
            app.ymin = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.ymin.Position = [68 217 43 22];
            app.ymin.Value = 0.5;

            % Create ymaxLabel
            app.ymaxLabel = uilabel(app.graphing_corner_heatmap);
            app.ymaxLabel.HorizontalAlignment = 'right';
            app.ymaxLabel.Position = [121 217 38 22];
            app.ymaxLabel.Text = 'y max';

            % Create ymax
            app.ymax = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.ymax.Position = [174 217 42 22];
            app.ymax.Value = 15;

            % Create HeatmapOutputOptionsLabel
            app.HeatmapOutputOptionsLabel = uilabel(app.graphing_corner_heatmap);
            app.HeatmapOutputOptionsLabel.FontSize = 22;
            app.HeatmapOutputOptionsLabel.FontWeight = 'bold';
            app.HeatmapOutputOptionsLabel.Position = [28 178 269 28];
            app.HeatmapOutputOptionsLabel.Text = 'Heatmap Output Options';

            % Create pixelsofresolutionLabel
            app.pixelsofresolutionLabel = uilabel(app.graphing_corner_heatmap);
            app.pixelsofresolutionLabel.HorizontalAlignment = 'right';
            app.pixelsofresolutionLabel.FontSize = 16;
            app.pixelsofresolutionLabel.Position = [58 34 151 22];
            app.pixelsofresolutionLabel.Text = 'pixels of resolution =';

            % Create n
            app.n = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.n.Limits = [1 1024];
            app.n.RoundFractionalValues = 'on';
            app.n.ValueDisplayFormat = '%.0f';
            app.n.Position = [224 34 39 22];
            app.n.Value = 30;

            % Create x_log_space
            app.x_log_space = uiswitch(app.graphing_corner_heatmap, 'slider');
            app.x_log_space.Items = {'lin', 'log'};
            app.x_log_space.Position = [251 299 27 12];
            app.x_log_space.Value = 'lin';

            % Create y_log_space
            app.y_log_space = uiswitch(app.graphing_corner_heatmap, 'slider');
            app.y_log_space.Items = {'lin', 'log'};
            app.y_log_space.Position = [251 222 27 12];
            app.y_log_space.Value = 'lin';

            % Create graphing_corner_sensitivity
            app.graphing_corner_sensitivity = uitab(app.graphing_corner);
            app.graphing_corner_sensitivity.Title = 'Sensitivity';
            app.graphing_corner_sensitivity.BackgroundColor = [1 1 0.6392];

            % Create VariablesListBoxLabel
            app.VariablesListBoxLabel = uilabel(app.graphing_corner_sensitivity);
            app.VariablesListBoxLabel.HorizontalAlignment = 'right';
            app.VariablesListBoxLabel.Position = [16 314 54 22];
            app.VariablesListBoxLabel.Text = 'Variables';

            % Create sensitivity_vars
            app.sensitivity_vars = uilistbox(app.graphing_corner_sensitivity);
            app.sensitivity_vars.Position = [85 124 225 214];

            % Create pressCTRLtoselectmultipleLabel
            app.pressCTRLtoselectmultipleLabel = uilabel(app.graphing_corner_sensitivity);
            app.pressCTRLtoselectmultipleLabel.Position = [26 238 51 73];
            app.pressCTRLtoselectmultipleLabel.Text = {'(press '; 'CTRL to '; 'select '; 'multiple)'};

            % Create MetricButtonGroup
            app.MetricButtonGroup = uibuttongroup(app.graphing_corner_sensitivity);
            app.MetricButtonGroup.Title = 'Metric';
            app.MetricButtonGroup.Position = [10 13 303 98];

            % Create ymaxButton
            app.ymaxButton = uiradiobutton(app.MetricButtonGroup);
            app.ymaxButton.Text = 'ymax';
            app.ymaxButton.Position = [11 51 58 22];

            % Create yunlatchButton
            app.yunlatchButton = uiradiobutton(app.MetricButtonGroup);
            app.yunlatchButton.Text = 'yunlatch';
            app.yunlatchButton.Position = [11 29 74 22];

            % Create tLButton
            app.tLButton = uiradiobutton(app.MetricButtonGroup);
            app.tLButton.Text = 'tL';
            app.tLButton.Position = [11 7 65 22];

            % Create vtoButton
            app.vtoButton = uiradiobutton(app.MetricButtonGroup);
            app.vtoButton.Text = 'vto';
            app.vtoButton.Position = [118 51 58 22];

            % Create PmaxButton
            app.PmaxButton = uiradiobutton(app.MetricButtonGroup);
            app.PmaxButton.Text = 'Pmax';
            app.PmaxButton.Position = [118 29 65 22];
            app.PmaxButton.Value = true;

            % Create ttoButton
            app.ttoButton = uiradiobutton(app.MetricButtonGroup);
            app.ttoButton.Text = 'tto';
            app.ttoButton.Position = [118 7 65 22];

            % Create KEmaxButton
            app.KEmaxButton = uiradiobutton(app.MetricButtonGroup);
            app.KEmaxButton.Text = 'KEmax';
            app.KEmaxButton.Position = [221 51 68 22];

            % Create go
            app.go = uibutton(app.UIFigure, 'push');
            app.go.ButtonPushedFcn = createCallbackFcn(app, @goButtonPushed, true);
            app.go.IconAlignment = 'bottom';
            app.go.BackgroundColor = [1 1 1];
            app.go.FontSize = 22;
            app.go.FontWeight = 'bold';
            app.go.Position = [799 167 138 38];
            app.go.Text = 'Graph!';

            % Create ShowModelSchematicButton
            app.ShowModelSchematicButton = uibutton(app.UIFigure, 'push');
            app.ShowModelSchematicButton.ButtonPushedFcn = createCallbackFcn(app, @ShowModelSchematicButtonPushed, true);
            app.ShowModelSchematicButton.BackgroundColor = [1 1 1];
            app.ShowModelSchematicButton.Position = [175 626 140 28];
            app.ShowModelSchematicButton.Text = 'Show Model Schematic';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [752 28 371 128];
            app.Image.ImageSource = 'lamsa.PNG';

            % Create savesolutionCheckBox
            app.savesolutionCheckBox = uicheckbox(app.UIFigure);
            app.savesolutionCheckBox.Text = 'save solution?';
            app.savesolutionCheckBox.Position = [952 175 99 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = plot_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

