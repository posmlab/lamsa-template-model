classdef plot_app < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        ParametersLabel                matlab.ui.control.Label
        spring                         matlab.ui.container.TabGroup
        linear_spring                  matlab.ui.container.Tab
        kLabel                         matlab.ui.control.Label
        linear_spring_k                matlab.ui.control.NumericEditField
        massEditFieldLabel_4           matlab.ui.control.Label
        linear_spring_mass             matlab.ui.control.NumericEditField
        FmaxLabel_3                    matlab.ui.control.Label
        linear_spring_Fmax             matlab.ui.control.NumericEditField
        RequiredFieldsLabel_2          matlab.ui.control.Label
        OptionalFieldsLabel_2          matlab.ui.control.Label
        defaultvaluesfilledinLabel_2   matlab.ui.control.Label
        exponential_spring             matlab.ui.container.Tab
        k_0Label                       matlab.ui.control.Label
        exp_spring_k                   matlab.ui.control.NumericEditField
        characteristiclengthLabel      matlab.ui.control.Label
        exp_spring_char_len            matlab.ui.control.NumericEditField
        FmaxLabel_4                    matlab.ui.control.Label
        exp_spring_Fmax                matlab.ui.control.NumericEditField
        massEditFieldLabel_5           matlab.ui.control.Label
        exp_spring_mass                matlab.ui.control.NumericEditField
        RequiredFieldsLabel            matlab.ui.control.Label
        OptionalFieldsLabel            matlab.ui.control.Label
        defaultvaluesfilledinLabel     matlab.ui.control.Label
        linear_elastic_extensional_spring  matlab.ui.container.Tab
        RequiredFieldsLabel_9          matlab.ui.control.Label
        OptionalFieldsLabel_7          matlab.ui.control.Label
        defaultvaluesfilledinLabel_7   matlab.ui.control.Label
        ELabel                         matlab.ui.control.Label
        lee_spring_E                   matlab.ui.control.NumericEditField
        rhoLabel                       matlab.ui.control.Label
        lee_spring_rho                 matlab.ui.control.NumericEditField
        ALabel                         matlab.ui.control.Label
        lee_spring_A                   matlab.ui.control.NumericEditField
        LLabel                         matlab.ui.control.Label
        lee_spring_L                   matlab.ui.control.NumericEditField
        sigma_fLabel                   matlab.ui.control.Label
        lee_spring_sigma_f             matlab.ui.control.NumericEditField
        loading_motor                  matlab.ui.container.TabGroup
        lm_linear_motor                matlab.ui.container.Tab
        FmaxLabel_2                    matlab.ui.control.Label
        lm_linear_motor_Fmax           matlab.ui.control.NumericEditField
        VmaxLabel_2                    matlab.ui.control.Label
        lm_linear_motor_Vmax           matlab.ui.control.NumericEditField
        RequiredFieldsLabel_3          matlab.ui.control.Label
        rangeofmotionLabel             matlab.ui.control.Label
        lm_linear_motor_range_of_motion  matlab.ui.control.NumericEditField
        OptionalFieldsLabel_8          matlab.ui.control.Label
        defaultvaluesfilledinLabel_8   matlab.ui.control.Label
        voltagefracLabel               matlab.ui.control.Label
        lm_linear_motor_voltage_frac   matlab.ui.control.NumericEditField
        lm_hill_muscle_motor           matlab.ui.container.Tab
        FmaxLabel                      matlab.ui.control.Label
        lm_hill_motor_Fmax             matlab.ui.control.NumericEditField
        VmaxLabel                      matlab.ui.control.Label
        lm_hill_motor_Vmax             matlab.ui.control.NumericEditField
        musclelengthLabel              matlab.ui.control.Label
        lm_hill_motor_muscle_length    matlab.ui.control.NumericEditField
        RequiredFieldsLabel_4          matlab.ui.control.Label
        OptionalFieldsLabel_4          matlab.ui.control.Label
        defaultvaluesfilledinLabel_4   matlab.ui.control.Label
        rateofactivationLabel          matlab.ui.control.Label
        lm_hill_motor_rate_of_activation  matlab.ui.control.NumericEditField
        initiallengthLabel             matlab.ui.control.Label
        lm_hill_motor_L_i              matlab.ui.control.NumericEditField
        a_LLabel                       matlab.ui.control.Label
        lm_hill_motor_a_L              matlab.ui.control.NumericEditField
        b_LLabel                       matlab.ui.control.Label
        lm_hill_motor_b_L              matlab.ui.control.NumericEditField
        sLabel                         matlab.ui.control.Label
        lm_hill_motor_s                matlab.ui.control.NumericEditField
        latch                          matlab.ui.container.TabGroup
        rounded_latch                  matlab.ui.container.Tab
        RequiredFieldsLabel_5          matlab.ui.control.Label
        OptionalFieldsLabel_5          matlab.ui.control.Label
        defaultvaluesfilledinLabel_5   matlab.ui.control.Label
        massEditFieldLabel_2           matlab.ui.control.Label
        latch_mass                     matlab.ui.control.NumericEditField
        Label                          matlab.ui.control.Label
        latch_coeff_fric               matlab.ui.control.NumericEditField
        radiusLabel                    matlab.ui.control.Label
        latch_radius                   matlab.ui.control.NumericEditField
        v_0Label                       matlab.ui.control.Label
        latch_v_0                      matlab.ui.control.NumericEditField
        minlatchingdistLabel           matlab.ui.control.Label
        min_latching_dist              matlab.ui.control.NumericEditField
        maxlatchingdistLabel           matlab.ui.control.Label
        max_latching_dist              matlab.ui.control.NumericEditField
        load                           matlab.ui.container.TabGroup
        load_mass                      matlab.ui.container.Tab
        RequiredFieldsLabel_6          matlab.ui.control.Label
        massEditFieldLabel             matlab.ui.control.Label
        load_mass_mass                 matlab.ui.control.NumericEditField
        OptionalFieldsLabel_10         matlab.ui.control.Label
        defaultvaluesfilledinLabel_10  matlab.ui.control.Label
        massofleverarmLabel            matlab.ui.control.Label
        load_m_rod                     matlab.ui.control.NumericEditField
        EMALabel                       matlab.ui.control.Label
        load_EMA                       matlab.ui.control.NumericEditField
        LoadingMotorLabel              matlab.ui.control.Label
        unlatching_motor               matlab.ui.container.TabGroup
        um_linear_motor                matlab.ui.container.Tab
        FmaxLabel_5                    matlab.ui.control.Label
        um_linear_motor_Fmax           matlab.ui.control.NumericEditField
        VmaxLabel_3                    matlab.ui.control.Label
        um_linear_motor_Vmax           matlab.ui.control.NumericEditField
        RequiredFieldsLabel_7          matlab.ui.control.Label
        rangeofmotionLabel_2           matlab.ui.control.Label
        um_linear_motor_range_of_motion  matlab.ui.control.NumericEditField
        voltagefracLabel_2             matlab.ui.control.Label
        um_linear_motor_voltage_frac   matlab.ui.control.NumericEditField
        OptionalFieldsLabel_9          matlab.ui.control.Label
        defaultvaluesfilledinLabel_9   matlab.ui.control.Label
        um_hill_muscle_motor           matlab.ui.container.Tab
        FmaxLabel_6                    matlab.ui.control.Label
        um_hill_motor_Fmax             matlab.ui.control.NumericEditField
        VmaxLabel_4                    matlab.ui.control.Label
        um_hill_motor_Vmax             matlab.ui.control.NumericEditField
        musclelengthLabel_2            matlab.ui.control.Label
        um_hill_motor_muscle_length    matlab.ui.control.NumericEditField
        RequiredFieldsLabel_8          matlab.ui.control.Label
        OptionalFieldsLabel_6          matlab.ui.control.Label
        defaultvaluesfilledinLabel_6   matlab.ui.control.Label
        rateofactivationLabel_2        matlab.ui.control.Label
        um_hill_motor_rate_of_activation  matlab.ui.control.NumericEditField
        initiallengthLabel_2           matlab.ui.control.Label
        um_hill_motor_L_i              matlab.ui.control.NumericEditField
        a_LLabel_2                     matlab.ui.control.Label
        um_hill_motor_a_L              matlab.ui.control.NumericEditField
        b_LLabel_2                     matlab.ui.control.Label
        um_hill_motor_b_L              matlab.ui.control.NumericEditField
        sLabel_2                       matlab.ui.control.Label
        um_hill_motor_s                matlab.ui.control.NumericEditField
        UnlatchingMotorLabel           matlab.ui.control.Label
        GraphingCornerLabel            matlab.ui.control.Label
        graphing_corner                matlab.ui.container.TabGroup
        graphing_corner_heatmap        matlab.ui.container.Tab
        minunlatchingmotorforceCheckBox  matlab.ui.control.CheckBox
        xaxisLabel                     matlab.ui.control.Label
        yaxisLabel                     matlab.ui.control.Label
        x_log_space                    matlab.ui.control.CheckBox
        y_log_space                    matlab.ui.control.CheckBox
        y_maxCheckBox                  matlab.ui.control.CheckBox
        y_unlatchCheckBox              matlab.ui.control.CheckBox
        t_LCheckBox                    matlab.ui.control.CheckBox
        v_toCheckBox                   matlab.ui.control.CheckBox
        P_maxCheckBox                  matlab.ui.control.CheckBox
        t_toCheckBox                   matlab.ui.control.CheckBox
        KE_maxCheckBox                 matlab.ui.control.CheckBox
        IV1DropDownLabel               matlab.ui.control.Label
        IV1DropDown                    matlab.ui.control.DropDown
        IV2DropDownLabel               matlab.ui.control.Label
        IV2DropDown                    matlab.ui.control.DropDown
        xminEditFieldLabel             matlab.ui.control.Label
        xmin                           matlab.ui.control.NumericEditField
        xmaxEditFieldLabel             matlab.ui.control.Label
        xmax                           matlab.ui.control.NumericEditField
        yminLabel                      matlab.ui.control.Label
        ymin                           matlab.ui.control.NumericEditField
        ymaxLabel                      matlab.ui.control.Label
        ymax                           matlab.ui.control.NumericEditField
        heatmapoutputoptionsLabel      matlab.ui.control.Label
        pixelsofresolutionLabel        matlab.ui.control.Label
        n                              matlab.ui.control.NumericEditField
        graphing_corner_one_D          matlab.ui.container.Tab
        OD_minunlatchingmotorforceCheckBox  matlab.ui.control.CheckBox
        xaxisLabel_2                   matlab.ui.control.Label
        OD_x_log_space                 matlab.ui.control.CheckBox
        OD_y_maxCheckBox               matlab.ui.control.CheckBox
        OD_y_unlatchCheckBox           matlab.ui.control.CheckBox
        OD_t_LCheckBox                 matlab.ui.control.CheckBox
        OD_v_toCheckBox                matlab.ui.control.CheckBox
        OD_P_maxCheckBox               matlab.ui.control.CheckBox
        OD_t_toCheckBox                matlab.ui.control.CheckBox
        OD_KE_maxCheckBox              matlab.ui.control.CheckBox
        IV1DropDownLabel_2             matlab.ui.control.Label
        OD_IV1DropDown                 matlab.ui.control.DropDown
        yaxisoutputoptionsLabel        matlab.ui.control.Label
        pixelsofresolutionLabel_2      matlab.ui.control.Label
        OD_n                           matlab.ui.control.NumericEditField
        xminEditFieldLabel_2           matlab.ui.control.Label
        OD_xmin                        matlab.ui.control.NumericEditField
        xmaxEditFieldLabel_2           matlab.ui.control.Label
        OD_xmax                        matlab.ui.control.NumericEditField
        graphing_corner_kinematics     matlab.ui.container.Tab
        go                             matlab.ui.control.Button
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
        end
        
        function heatmap(app)
            % determines resolution of heatplots
            N=app.n.Value; 
            
            % initializing waitbar
            f = waitbar(0,'Please wait...');
            load_bar_value = 1/N;
            load_bar_increment = 1/N;
            
            % output directory initialization
            output_directory = create_output_directory();
            
            %% plot parameters

            % setting x axis on the plot
            xname = app.dropdown_items_opposite_dict(app.IV1DropDown.Value);
            xrange = [app.xmin.Value app.xmax.Value];
            if app.x_log_space.Value
                looping_value_x = logspace(xrange(1),xrange(2),N);
            else
                looping_value_x = linspace(xrange(1),xrange(2),N);
            end
            
            %setting y axis value on plot
            yname = app.dropdown_items_opposite_dict(app.IV2DropDown.Value);
            yrange = [app.ymin.Value app.ymax.Value];
            if app.y_log_space.Value
                looping_value_y = logspace(yrange(1),yrange(2),N);
            else
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
            
            if isempty(metrics)
                errordlg('You must pick at least one output option','Error');
                error('You must pick at least one output')
            end
            
            metrics_dict = containers.Map(metrics_names,metrics_labels);
            
            looping_param_x = app.dropdown_items_opposite_dict(app.IV1DropDown.Value);
            looping_param_y = app.dropdown_items_opposite_dict(app.IV2DropDown.Value);
            
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
                    
                    eval(['app.' looping_param_x '.Value = ' num2str(looping_value_x(j)) ';']);
                    eval(['app.' looping_param_y '.Value = ' num2str(looping_value_y(i)) ';']);
                    
                    %% initializing LaMSA component structs
            
                    % load mass struct initialization
                    load = load_mass(app.load_mass_mass.Value,app.load_m_rod.Value,app.load_EMA.Value);
                    
                    % latch struct initialization
                    latch = rounded_latch(app.latch_radius.Value,app.latch_mass.Value,app.latch_coeff_fric.Value, app.latch_v_0.Value, app.min_latching_dist.Value, app.max_latching_dist.Value);
                    
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
                    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
                    
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
                        outval{ii}(i,j)=met_dict(metrics{ii});
                    end
                end
                disp(['row ' num2str(i) ' of ' num2str(N)]);
                load_bar_value = load_bar_value + load_bar_increment;
                waitbar(load_bar_value,f,'Processing...');
            end
            
            close(f)
            
            % plot output
            fh = figure();
            fh.WindowState = 'maximized';
            subplot_rows = floor(sqrt(length(metrics)));
            subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
            for ii=1:length(metrics)
                subplot(subplot_rows,subplot_cols,ii);
                
                imagesc(xrange,yrange,outval{ii});
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
                if app.x_log_space.Value
                    for i=1:length(xticks)
                        xticklabel{i} = strcat('$10^{',num2str(xtick(i)),'}$');
                    end
                end
                set(gca,'XTickLabel',xticklabel);
                set(gca,'XTick',xtick);
                if app.y_log_space.Value
                    for j=1:length(yticks)
                        yticklabel{j} = strcat('$10^{',num2str(ytick(j)),'}$');
                    end
                end
                set(gca,'YTickLabel',yticklabel);
                set(gca,'YTick',ytick);
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function one_D_plot(app)
            N=app.OD_n.Value; 
            
            % initializing waitbar
            f = waitbar(0,'Please wait...');
            load_bar_value = 1/N;
            load_bar_increment = 1/N;
            
            % output directory initialization
            output_directory = create_output_directory();
            
            %% plot parameters

            % setting x axis on the plot
            xname = app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value);
            xrange = [app.OD_xmin.Value app.OD_xmax.Value];
            if app.OD_x_log_space.Value
                looping_value_x = logspace(xrange(1),xrange(2),N);
            else
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
            
            if isempty(metrics)
                errordlg('You must pick at least one output option','Error');
                error('You must pick at least one output')
            end
            
            metrics_names = {'ymax','yunlatch','tL','vto','Pmax','tto','minumforce','KEmax'};
            metrics_labels = {'$y_{\textrm{max}}$','$y_{\textrm{unlatch}}$','$t_L$','$v_{\textrm{to}}$','$P_{\textrm{max}}$','$t_{\textrm{to}}$','min unlatching force','$KE_{\textrm{max}}$'};
            metrics_dict = containers.Map(metrics_names,metrics_labels);
            
            looping_param_x = app.dropdown_items_opposite_dict(app.OD_IV1DropDown.Value);
            
            for ii=1:length(metrics)
                outval{ii} = zeros(size(looping_value_x));
            end
            for i=1:N    
                eval(['app.' looping_param_x '.Value = ' num2str(looping_value_x(i)) ';']);
                
                %% initializing LaMSA component structs
        
                % load mass struct initialization
                load = load_mass(app.load_mass_mass.Value,app.load_m_rod.Value,app.load_EMA.Value);
                
                % latch struct initialization
                latch = rounded_latch(app.latch_radius.Value,app.latch_mass.Value,app.latch_coeff_fric.Value, app.latch_v_0.Value, app.min_latching_dist.Value, app.max_latching_dist.Value);
                
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
                [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
                
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
                    outval{ii}(i)=met_dict(metrics{ii});
                end
                
                disp(['row ' num2str(i) ' of ' num2str(N)]);
                load_bar_value = load_bar_value + load_bar_increment;
                waitbar(load_bar_value,f,'Processing...');
            end

            close(f)
            
            % plot output
            fh = figure();
            fh.WindowState = 'maximized';
            subplot_rows = floor(sqrt(length(metrics)));
            subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
            for ii=1:length(metrics)
                subplot(subplot_rows,subplot_cols,ii);
                size(outval{ii})
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
                
                % makes the x in log scale if needed
                if app.OD_x_log_space.Value
                    set(gca,'XScale','log')
                end
               
            end
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function kinematics(app)
            % initializing output directory
            output_directory = create_output_directory();
            
            %% initializing LaMSA component structs
            
            % load mass struct initialization
            load = load_mass(app.load_mass_mass.Value,app.load_m_rod.Value,app.load_EMA.Value);
            
            % latch struct initialization
            latch = rounded_latch(app.latch_radius.Value,app.latch_mass.Value,app.latch_coeff_fric.Value, app.latch_v_0.Value, app.min_latching_dist.Value, app.max_latching_dist.Value);
            
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
        
            [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring,output_directory);
            columntitles=["Time", "y postition", "y velocity", "x position", "x velocity", "normal force on latch x", ...
                "normal force on load y", "frictional force on latch x", ...
                "frictional force on load y", "spring force", ...
                "unlatching motor force into"];
            units=[" [s]"," [m]"," [m/s]"];
            
            % latch kinematics
            figure
            for i = 4:5
                subplot(3,1,i-3)
                box on
                set(gca,'TickLabelInterpreter','latex')
                hold on
                plot(sol(:,1),sol(:,i));
                hold off
                title(columntitles(i),"Interpreter","latex");
                ylabel(columntitles(i)+units(i-2),"Interpreter","latex");
                xlabel(columntitles(1)+units(1),"Interpreter","latex");
            end
            
            for i = [6 8 11]
                subplot(3,1,3)
                box on
                set(gca,'TickLabelInterpreter','latex')
                hold on
                plot(sol(:,1),sol(:,i),"DisplayName",columntitles(i))
            end
            title("Latch Force Components","Interpreter","latex");
            ylabel("Force [N]","Interpreter","latex");
            xlabel(columntitles(1)+" [s]","Interpreter","latex");
            legend("show")
            hold off
            
            % load kinematics
            figure
            for i = 2:3
                subplot(3,1,i-1)
                box on
                set(gca,'TickLabelInterpreter','latex')
                hold on
                plot(sol(:,1),sol(:,i));
                hold off
                title(columntitles(i),"Interpreter","latex");
                ylabel(columntitles(i)+units(i),"Interpreter","latex");
                xlabel(columntitles(1)+units(1),"Interpreter","latex");
            end
            
            for i = [7 9 10]
                subplot(3,1,3)
                box on
                set(gca,'TickLabelInterpreter','latex')
                hold on
                plot(sol(:,1),sol(:,i),"DisplayName",columntitles(i))
            end
            title("Load Force Components","Interpreter","latex");
            ylabel("Force [N]","Interpreter","latex");
            xlabel(columntitles(1)+" [s]","Interpreter","latex");
            legend("show")
            hold off
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            
            % makes comparing changes in the app in github easier
            pathname = char(which("plot_app.mlapp"));
            lastIndex = strlength(pathname);
            nameLength = strlength("plot_app.mlapp");
            correctIndex = lastIndex-nameLength;
            pathname = pathname(1:correctIndex);
            pathname = pathname +"/../ext";
            addpath(fullfile(pathname))
            mlapp2classdef("plot_app.mlapp")
            
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
            
            app.axis_labels_dict = containers.Map(varnames,latexlabels);
            app.dropdown_items_dict = containers.Map(varnames,nonlatexlabels);
            app.dropdown_items_opposite_dict = containers.Map(nonlatexlabels,varnames);

            
            app.load_vars = {'load_mass_mass'};
            app.latch_vars = {'latch_mass','latch_coeff_fric','latch_radius','latch_v_0','min_latching_dist','max_latching_dist'};
            app.spring_vars = {'linear_spring_k','linear_spring_mass','linear_spring_Fmax'};
            app.lm_vars = {'lm_linear_motor_Fmax','lm_linear_motor_Vmax','lm_linear_motor_range_of_motion','lm_linear_motor_voltage_frac'};
            app.um_vars = {'um_linear_motor_Fmax','um_linear_motor_Vmax','um_linear_motor_range_of_motion','um_linear_motor_voltage_frac'};
            update_dd_vars(app);
            
            app.IV2DropDown.Items(ismember(app.IV2DropDown.Items,app.IV1DropDown.Value)) = [];
        end

        % Button pushed function: go
        function goButtonPushed(app, event)
            addpath(genpath(fullfile(pwd,'..')));
            lastwarn('');
            
            if (app.graphing_corner.SelectedTab == app.graphing_corner_heatmap)
                heatmap(app);
            elseif (app.graphing_corner.SelectedTab == app.graphing_corner_one_D)
                one_D_plot(app);
            elseif (app.graphing_corner.SelectedTab == app.graphing_corner_kinematics)
                kinematics(app);
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
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.902 0.902 0.902];
            app.UIFigure.Position = [60 -10 1207 803];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Scrollable = 'on';

            % Create ParametersLabel
            app.ParametersLabel = uilabel(app.UIFigure);
            app.ParametersLabel.FontSize = 22;
            app.ParametersLabel.FontWeight = 'bold';
            app.ParametersLabel.Position = [29 754 125 27];
            app.ParametersLabel.Text = 'Parameters';

            % Create spring
            app.spring = uitabgroup(app.UIFigure);
            app.spring.AutoResizeChildren = 'off';
            app.spring.SelectionChangedFcn = createCallbackFcn(app, @springSelectionChanged, true);
            app.spring.Position = [32 338 795 115];

            % Create linear_spring
            app.linear_spring = uitab(app.spring);
            app.linear_spring.AutoResizeChildren = 'off';
            app.linear_spring.Title = 'Linear Spring';

            % Create kLabel
            app.kLabel = uilabel(app.linear_spring);
            app.kLabel.HorizontalAlignment = 'right';
            app.kLabel.Position = [162 56 44 22];
            app.kLabel.Text = 'k =';

            % Create linear_spring_k
            app.linear_spring_k = uieditfield(app.linear_spring, 'numeric');
            app.linear_spring_k.Limits = [0 Inf];
            app.linear_spring_k.Tooltip = {'stiffness of the spring'};
            app.linear_spring_k.Position = [221 56 41 22];
            app.linear_spring_k.Value = 0.5;

            % Create massEditFieldLabel_4
            app.massEditFieldLabel_4 = uilabel(app.linear_spring);
            app.massEditFieldLabel_4.HorizontalAlignment = 'right';
            app.massEditFieldLabel_4.Position = [162 17 44 22];
            app.massEditFieldLabel_4.Text = 'mass =';

            % Create linear_spring_mass
            app.linear_spring_mass = uieditfield(app.linear_spring, 'numeric');
            app.linear_spring_mass.Limits = [0 Inf];
            app.linear_spring_mass.Tooltip = {'mass of the spring'};
            app.linear_spring_mass.Position = [221 17 41 22];

            % Create FmaxLabel_3
            app.FmaxLabel_3 = uilabel(app.linear_spring);
            app.FmaxLabel_3.HorizontalAlignment = 'right';
            app.FmaxLabel_3.Position = [325 17 46 22];
            app.FmaxLabel_3.Text = 'Fmax =';

            % Create linear_spring_Fmax
            app.linear_spring_Fmax = uieditfield(app.linear_spring, 'numeric');
            app.linear_spring_Fmax.Limits = [0 Inf];
            app.linear_spring_Fmax.Tooltip = {'max force of spring'};
            app.linear_spring_Fmax.Position = [386 17 43 22];
            app.linear_spring_Fmax.Value = Inf;

            % Create RequiredFieldsLabel_2
            app.RequiredFieldsLabel_2 = uilabel(app.linear_spring);
            app.RequiredFieldsLabel_2.Position = [18 56 92 22];
            app.RequiredFieldsLabel_2.Text = 'Required Fields:';

            % Create OptionalFieldsLabel_2
            app.OptionalFieldsLabel_2 = uilabel(app.linear_spring);
            app.OptionalFieldsLabel_2.Position = [18 17 89 22];
            app.OptionalFieldsLabel_2.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_2
            app.defaultvaluesfilledinLabel_2 = uilabel(app.linear_spring);
            app.defaultvaluesfilledinLabel_2.FontSize = 9;
            app.defaultvaluesfilledinLabel_2.Position = [14 6 98 13];
            app.defaultvaluesfilledinLabel_2.Text = '(default values filled in)';

            % Create exponential_spring
            app.exponential_spring = uitab(app.spring);
            app.exponential_spring.AutoResizeChildren = 'off';
            app.exponential_spring.Title = 'Exponential Spring';

            % Create k_0Label
            app.k_0Label = uilabel(app.exponential_spring);
            app.k_0Label.HorizontalAlignment = 'right';
            app.k_0Label.Position = [162 56 44 22];
            app.k_0Label.Text = 'k_0 =';

            % Create exp_spring_k
            app.exp_spring_k = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_k.Limits = [0 Inf];
            app.exp_spring_k.Tooltip = {'spring stiffness'};
            app.exp_spring_k.Position = [221 56 46 22];
            app.exp_spring_k.Value = 2000;

            % Create characteristiclengthLabel
            app.characteristiclengthLabel = uilabel(app.exponential_spring);
            app.characteristiclengthLabel.HorizontalAlignment = 'right';
            app.characteristiclengthLabel.Position = [332 56 123 22];
            app.characteristiclengthLabel.Text = 'characteristic length =';

            % Create exp_spring_char_len
            app.exp_spring_char_len = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_char_len.Limits = [0 Inf];
            app.exp_spring_char_len.Tooltip = {'resting length where force is 0'};
            app.exp_spring_char_len.Position = [470 56 38 22];
            app.exp_spring_char_len.Value = 0.001;

            % Create FmaxLabel_4
            app.FmaxLabel_4 = uilabel(app.exponential_spring);
            app.FmaxLabel_4.HorizontalAlignment = 'right';
            app.FmaxLabel_4.Position = [411 17 46 22];
            app.FmaxLabel_4.Text = 'Fmax =';

            % Create exp_spring_Fmax
            app.exp_spring_Fmax = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_Fmax.Limits = [0 Inf];
            app.exp_spring_Fmax.Tooltip = {'max force of the spring'};
            app.exp_spring_Fmax.Position = [472 17 36 22];
            app.exp_spring_Fmax.Value = Inf;

            % Create massEditFieldLabel_5
            app.massEditFieldLabel_5 = uilabel(app.exponential_spring);
            app.massEditFieldLabel_5.HorizontalAlignment = 'right';
            app.massEditFieldLabel_5.Position = [162 17 44 22];
            app.massEditFieldLabel_5.Text = 'mass =';

            % Create exp_spring_mass
            app.exp_spring_mass = uieditfield(app.exponential_spring, 'numeric');
            app.exp_spring_mass.Limits = [0 Inf];
            app.exp_spring_mass.Tooltip = {'mass of the spring'};
            app.exp_spring_mass.Position = [221 17 46 22];
            app.exp_spring_mass.Value = 0.01;

            % Create RequiredFieldsLabel
            app.RequiredFieldsLabel = uilabel(app.exponential_spring);
            app.RequiredFieldsLabel.Position = [18 56 92 22];
            app.RequiredFieldsLabel.Text = 'Required Fields:';

            % Create OptionalFieldsLabel
            app.OptionalFieldsLabel = uilabel(app.exponential_spring);
            app.OptionalFieldsLabel.Position = [18 17 89 22];
            app.OptionalFieldsLabel.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel
            app.defaultvaluesfilledinLabel = uilabel(app.exponential_spring);
            app.defaultvaluesfilledinLabel.FontSize = 9;
            app.defaultvaluesfilledinLabel.Position = [14 6 98 13];
            app.defaultvaluesfilledinLabel.Text = '(default values filled in)';

            % Create linear_elastic_extensional_spring
            app.linear_elastic_extensional_spring = uitab(app.spring);
            app.linear_elastic_extensional_spring.Title = 'Linear Elastic Extensional Spring';

            % Create RequiredFieldsLabel_9
            app.RequiredFieldsLabel_9 = uilabel(app.linear_elastic_extensional_spring);
            app.RequiredFieldsLabel_9.Position = [18 56 92 22];
            app.RequiredFieldsLabel_9.Text = 'Required Fields:';

            % Create OptionalFieldsLabel_7
            app.OptionalFieldsLabel_7 = uilabel(app.linear_elastic_extensional_spring);
            app.OptionalFieldsLabel_7.Position = [18 17 89 22];
            app.OptionalFieldsLabel_7.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_7
            app.defaultvaluesfilledinLabel_7 = uilabel(app.linear_elastic_extensional_spring);
            app.defaultvaluesfilledinLabel_7.FontSize = 9;
            app.defaultvaluesfilledinLabel_7.Position = [14 6 98 13];
            app.defaultvaluesfilledinLabel_7.Text = '(default values filled in)';

            % Create ELabel
            app.ELabel = uilabel(app.linear_elastic_extensional_spring);
            app.ELabel.HorizontalAlignment = 'right';
            app.ELabel.Position = [162 56 44 22];
            app.ELabel.Text = 'E =';

            % Create lee_spring_E
            app.lee_spring_E = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_E.Limits = [0 Inf];
            app.lee_spring_E.Tooltip = {'modulus'};
            app.lee_spring_E.Position = [221 56 34 22];
            app.lee_spring_E.Value = 0.5;

            % Create rhoLabel
            app.rhoLabel = uilabel(app.linear_elastic_extensional_spring);
            app.rhoLabel.HorizontalAlignment = 'right';
            app.rhoLabel.Position = [162 17 44 22];
            app.rhoLabel.Text = 'rho =';

            % Create lee_spring_rho
            app.lee_spring_rho = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_rho.Limits = [0 Inf];
            app.lee_spring_rho.Tooltip = {'density'};
            app.lee_spring_rho.Position = [221 17 34 22];
            app.lee_spring_rho.Value = 10;

            % Create ALabel
            app.ALabel = uilabel(app.linear_elastic_extensional_spring);
            app.ALabel.HorizontalAlignment = 'right';
            app.ALabel.Position = [307 56 44 22];
            app.ALabel.Text = 'A =';

            % Create lee_spring_A
            app.lee_spring_A = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_A.Limits = [0 Inf];
            app.lee_spring_A.Tooltip = {'cross sectional area'};
            app.lee_spring_A.Position = [366 56 34 22];
            app.lee_spring_A.Value = 0.5;

            % Create LLabel
            app.LLabel = uilabel(app.linear_elastic_extensional_spring);
            app.LLabel.HorizontalAlignment = 'right';
            app.LLabel.Position = [445 56 44 22];
            app.LLabel.Text = 'L =';

            % Create lee_spring_L
            app.lee_spring_L = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_L.Limits = [0 Inf];
            app.lee_spring_L.Tooltip = {'length'};
            app.lee_spring_L.Position = [504 56 34 22];
            app.lee_spring_L.Value = 0.001;

            % Create sigma_fLabel
            app.sigma_fLabel = uilabel(app.linear_elastic_extensional_spring);
            app.sigma_fLabel.HorizontalAlignment = 'right';
            app.sigma_fLabel.Position = [293 17 58 22];
            app.sigma_fLabel.Text = 'sigma_f =';

            % Create lee_spring_sigma_f
            app.lee_spring_sigma_f = uieditfield(app.linear_elastic_extensional_spring, 'numeric');
            app.lee_spring_sigma_f.Limits = [0 Inf];
            app.lee_spring_sigma_f.Tooltip = {'failure strength in Pa'};
            app.lee_spring_sigma_f.Position = [366 17 34 22];
            app.lee_spring_sigma_f.Value = Inf;

            % Create loading_motor
            app.loading_motor = uitabgroup(app.UIFigure);
            app.loading_motor.AutoResizeChildren = 'off';
            app.loading_motor.SelectionChangedFcn = createCallbackFcn(app, @loading_motorSelectionChanged, true);
            app.loading_motor.Position = [31 185 795 121];

            % Create lm_linear_motor
            app.lm_linear_motor = uitab(app.loading_motor);
            app.lm_linear_motor.AutoResizeChildren = 'off';
            app.lm_linear_motor.Title = 'Linear Motor';

            % Create FmaxLabel_2
            app.FmaxLabel_2 = uilabel(app.lm_linear_motor);
            app.FmaxLabel_2.HorizontalAlignment = 'right';
            app.FmaxLabel_2.Position = [157 62 46 22];
            app.FmaxLabel_2.Text = 'Fmax =';

            % Create lm_linear_motor_Fmax
            app.lm_linear_motor_Fmax = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_Fmax.Limits = [0 Inf];
            app.lm_linear_motor_Fmax.Tooltip = {'max force of the motor'};
            app.lm_linear_motor_Fmax.Position = [218 62 40 22];
            app.lm_linear_motor_Fmax.Value = 10;

            % Create VmaxLabel_2
            app.VmaxLabel_2 = uilabel(app.lm_linear_motor);
            app.VmaxLabel_2.HorizontalAlignment = 'right';
            app.VmaxLabel_2.Position = [293 62 47 22];
            app.VmaxLabel_2.Text = 'Vmax =';

            % Create lm_linear_motor_Vmax
            app.lm_linear_motor_Vmax = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_Vmax.Limits = [0 Inf];
            app.lm_linear_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.lm_linear_motor_Vmax.Position = [355 62 37 22];
            app.lm_linear_motor_Vmax.Value = 10;

            % Create RequiredFieldsLabel_3
            app.RequiredFieldsLabel_3 = uilabel(app.lm_linear_motor);
            app.RequiredFieldsLabel_3.Position = [16 62 92 22];
            app.RequiredFieldsLabel_3.Text = 'Required Fields:';

            % Create rangeofmotionLabel
            app.rangeofmotionLabel = uilabel(app.lm_linear_motor);
            app.rangeofmotionLabel.HorizontalAlignment = 'right';
            app.rangeofmotionLabel.Position = [430 62 99 22];
            app.rangeofmotionLabel.Text = 'range of motion =';

            % Create lm_linear_motor_range_of_motion
            app.lm_linear_motor_range_of_motion = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_range_of_motion.Limits = [0 Inf];
            app.lm_linear_motor_range_of_motion.Tooltip = {'total possible range of motion of the motor'};
            app.lm_linear_motor_range_of_motion.Position = [544 62 37 22];
            app.lm_linear_motor_range_of_motion.Value = 0.005;

            % Create OptionalFieldsLabel_8
            app.OptionalFieldsLabel_8 = uilabel(app.lm_linear_motor);
            app.OptionalFieldsLabel_8.Position = [18 24 89 22];
            app.OptionalFieldsLabel_8.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_8
            app.defaultvaluesfilledinLabel_8 = uilabel(app.lm_linear_motor);
            app.defaultvaluesfilledinLabel_8.FontSize = 9;
            app.defaultvaluesfilledinLabel_8.Position = [14 13 98 13];
            app.defaultvaluesfilledinLabel_8.Text = '(default values filled in)';

            % Create voltagefracLabel
            app.voltagefracLabel = uilabel(app.lm_linear_motor);
            app.voltagefracLabel.HorizontalAlignment = 'right';
            app.voltagefracLabel.Position = [127 23 78 22];
            app.voltagefracLabel.Text = 'voltage frac =';

            % Create lm_linear_motor_voltage_frac
            app.lm_linear_motor_voltage_frac = uieditfield(app.lm_linear_motor, 'numeric');
            app.lm_linear_motor_voltage_frac.Limits = [0 Inf];
            app.lm_linear_motor_voltage_frac.Tooltip = {'fraction of max voltage of motor'};
            app.lm_linear_motor_voltage_frac.Position = [220 23 40 22];
            app.lm_linear_motor_voltage_frac.Value = 1;

            % Create lm_hill_muscle_motor
            app.lm_hill_muscle_motor = uitab(app.loading_motor);
            app.lm_hill_muscle_motor.AutoResizeChildren = 'off';
            app.lm_hill_muscle_motor.Title = 'Hill Muscle Motor';

            % Create FmaxLabel
            app.FmaxLabel = uilabel(app.lm_hill_muscle_motor);
            app.FmaxLabel.HorizontalAlignment = 'right';
            app.FmaxLabel.Position = [161 63 46 22];
            app.FmaxLabel.Text = 'Fmax =';

            % Create lm_hill_motor_Fmax
            app.lm_hill_motor_Fmax = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_Fmax.Limits = [0 Inf];
            app.lm_hill_motor_Fmax.Tooltip = {'max force of the motor'};
            app.lm_hill_motor_Fmax.Position = [222 63 34 22];
            app.lm_hill_motor_Fmax.Value = 10;

            % Create VmaxLabel
            app.VmaxLabel = uilabel(app.lm_hill_muscle_motor);
            app.VmaxLabel.HorizontalAlignment = 'right';
            app.VmaxLabel.Position = [292 62 47 22];
            app.VmaxLabel.Text = 'Vmax =';

            % Create lm_hill_motor_Vmax
            app.lm_hill_motor_Vmax = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_Vmax.Limits = [0 Inf];
            app.lm_hill_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.lm_hill_motor_Vmax.Position = [354 62 37 22];
            app.lm_hill_motor_Vmax.Value = 10;

            % Create musclelengthLabel
            app.musclelengthLabel = uilabel(app.lm_hill_muscle_motor);
            app.musclelengthLabel.HorizontalAlignment = 'right';
            app.musclelengthLabel.Position = [420 62 90 22];
            app.musclelengthLabel.Text = 'muscle length =';

            % Create lm_hill_motor_muscle_length
            app.lm_hill_motor_muscle_length = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_muscle_length.Limits = [0 Inf];
            app.lm_hill_motor_muscle_length.ValueChangedFcn = createCallbackFcn(app, @lm_hill_motor_muscle_lengthValueChanged, true);
            app.lm_hill_motor_muscle_length.Tooltip = {''};
            app.lm_hill_motor_muscle_length.Position = [525 62 36 22];
            app.lm_hill_motor_muscle_length.Value = 0.01;

            % Create RequiredFieldsLabel_4
            app.RequiredFieldsLabel_4 = uilabel(app.lm_hill_muscle_motor);
            app.RequiredFieldsLabel_4.Position = [16 62 92 22];
            app.RequiredFieldsLabel_4.Text = 'Required Fields:';

            % Create OptionalFieldsLabel_4
            app.OptionalFieldsLabel_4 = uilabel(app.lm_hill_muscle_motor);
            app.OptionalFieldsLabel_4.Position = [17 23 89 22];
            app.OptionalFieldsLabel_4.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_4
            app.defaultvaluesfilledinLabel_4 = uilabel(app.lm_hill_muscle_motor);
            app.defaultvaluesfilledinLabel_4.FontSize = 9;
            app.defaultvaluesfilledinLabel_4.Position = [13 12 98 13];
            app.defaultvaluesfilledinLabel_4.Text = '(default values filled in)';

            % Create rateofactivationLabel
            app.rateofactivationLabel = uilabel(app.lm_hill_muscle_motor);
            app.rateofactivationLabel.HorizontalAlignment = 'right';
            app.rateofactivationLabel.Position = [589 63 104 22];
            app.rateofactivationLabel.Text = 'rate of activation =';

            % Create lm_hill_motor_rate_of_activation
            app.lm_hill_motor_rate_of_activation = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_rate_of_activation.Limits = [0 Inf];
            app.lm_hill_motor_rate_of_activation.Tooltip = {''};
            app.lm_hill_motor_rate_of_activation.Position = [708 63 37 22];
            app.lm_hill_motor_rate_of_activation.Value = Inf;

            % Create initiallengthLabel
            app.initiallengthLabel = uilabel(app.lm_hill_muscle_motor);
            app.initiallengthLabel.HorizontalAlignment = 'right';
            app.initiallengthLabel.Position = [128 23 79 22];
            app.initiallengthLabel.Text = 'initial length =';

            % Create lm_hill_motor_L_i
            app.lm_hill_motor_L_i = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_L_i.Limits = [0 Inf];
            app.lm_hill_motor_L_i.Tooltip = {''};
            app.lm_hill_motor_L_i.Position = [222 23 33 22];
            app.lm_hill_motor_L_i.Value = 0.01;

            % Create a_LLabel
            app.a_LLabel = uilabel(app.lm_hill_muscle_motor);
            app.a_LLabel.HorizontalAlignment = 'right';
            app.a_LLabel.Position = [303 23 36 22];
            app.a_LLabel.Text = 'a_L =';

            % Create lm_hill_motor_a_L
            app.lm_hill_motor_a_L = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_a_L.Tooltip = {'determines shape of length-tension relationship'};
            app.lm_hill_motor_a_L.Position = [354 23 37 22];
            app.lm_hill_motor_a_L.Value = 2.08;

            % Create b_LLabel
            app.b_LLabel = uilabel(app.lm_hill_muscle_motor);
            app.b_LLabel.HorizontalAlignment = 'right';
            app.b_LLabel.Position = [472 24 36 22];
            app.b_LLabel.Text = 'b_L =';

            % Create lm_hill_motor_b_L
            app.lm_hill_motor_b_L = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_b_L.Tooltip = {'determines shape of length-tension relationship'};
            app.lm_hill_motor_b_L.Position = [523 24 39 22];
            app.lm_hill_motor_b_L.Value = -2.89;

            % Create sLabel
            app.sLabel = uilabel(app.lm_hill_muscle_motor);
            app.sLabel.HorizontalAlignment = 'right';
            app.sLabel.Position = [668 23 25 22];
            app.sLabel.Text = 's =';

            % Create lm_hill_motor_s
            app.lm_hill_motor_s = uieditfield(app.lm_hill_muscle_motor, 'numeric');
            app.lm_hill_motor_s.Tooltip = {'determines shape of length-tension relationship'};
            app.lm_hill_motor_s.Position = [708 23 37 22];
            app.lm_hill_motor_s.Value = -0.75;

            % Create latch
            app.latch = uitabgroup(app.UIFigure);
            app.latch.AutoResizeChildren = 'off';
            app.latch.Position = [32 472 795 116];

            % Create rounded_latch
            app.rounded_latch = uitab(app.latch);
            app.rounded_latch.AutoResizeChildren = 'off';
            app.rounded_latch.Title = 'Rounded Latch';

            % Create RequiredFieldsLabel_5
            app.RequiredFieldsLabel_5 = uilabel(app.rounded_latch);
            app.RequiredFieldsLabel_5.Position = [18 57 92 22];
            app.RequiredFieldsLabel_5.Text = 'Required Fields:';

            % Create OptionalFieldsLabel_5
            app.OptionalFieldsLabel_5 = uilabel(app.rounded_latch);
            app.OptionalFieldsLabel_5.Position = [18 18 89 22];
            app.OptionalFieldsLabel_5.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_5
            app.defaultvaluesfilledinLabel_5 = uilabel(app.rounded_latch);
            app.defaultvaluesfilledinLabel_5.FontSize = 9;
            app.defaultvaluesfilledinLabel_5.Position = [14 7 98 13];
            app.defaultvaluesfilledinLabel_5.Text = '(default values filled in)';

            % Create massEditFieldLabel_2
            app.massEditFieldLabel_2 = uilabel(app.rounded_latch);
            app.massEditFieldLabel_2.HorizontalAlignment = 'right';
            app.massEditFieldLabel_2.Position = [162 57 44 22];
            app.massEditFieldLabel_2.Text = 'mass =';

            % Create latch_mass
            app.latch_mass = uieditfield(app.rounded_latch, 'numeric');
            app.latch_mass.Limits = [0 Inf];
            app.latch_mass.Tooltip = {'mass of the latch'};
            app.latch_mass.Position = [221 57 40 22];
            app.latch_mass.Value = 100;

            % Create Label
            app.Label = uilabel(app.rounded_latch);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [181 18 25 22];
            app.Label.Text = 'μ =';

            % Create latch_coeff_fric
            app.latch_coeff_fric = uieditfield(app.rounded_latch, 'numeric');
            app.latch_coeff_fric.Limits = [0 Inf];
            app.latch_coeff_fric.Tooltip = {'coefficient of friction between the latch and load mass'};
            app.latch_coeff_fric.Position = [221 18 40 22];

            % Create radiusLabel
            app.radiusLabel = uilabel(app.rounded_latch);
            app.radiusLabel.HorizontalAlignment = 'right';
            app.radiusLabel.Position = [288 57 48 22];
            app.radiusLabel.Text = 'radius =';

            % Create latch_radius
            app.latch_radius = uieditfield(app.rounded_latch, 'numeric');
            app.latch_radius.Limits = [0 Inf];
            app.latch_radius.Tooltip = {'latch radius'};
            app.latch_radius.Position = [351 57 48 22];
            app.latch_radius.Value = 0.005;

            % Create v_0Label
            app.v_0Label = uilabel(app.rounded_latch);
            app.v_0Label.HorizontalAlignment = 'right';
            app.v_0Label.Position = [301 18 38 22];
            app.v_0Label.Text = 'v_0 =';

            % Create latch_v_0
            app.latch_v_0 = uieditfield(app.rounded_latch, 'numeric');
            app.latch_v_0.Limits = [0 Inf];
            app.latch_v_0.Tooltip = {'initial velocity of the latch'};
            app.latch_v_0.Position = [351 18 48 22];
            app.latch_v_0.Value = 100;

            % Create minlatchingdistLabel
            app.minlatchingdistLabel = uilabel(app.rounded_latch);
            app.minlatchingdistLabel.HorizontalAlignment = 'right';
            app.minlatchingdistLabel.Position = [419 18 102 22];
            app.minlatchingdistLabel.Text = 'min latching dist =';

            % Create min_latching_dist
            app.min_latching_dist = uieditfield(app.rounded_latch, 'numeric');
            app.min_latching_dist.Limits = [0 Inf];
            app.min_latching_dist.Tooltip = {'lowest position the latch can engage'};
            app.min_latching_dist.Position = [533 18 48 22];

            % Create maxlatchingdistLabel
            app.maxlatchingdistLabel = uilabel(app.rounded_latch);
            app.maxlatchingdistLabel.HorizontalAlignment = 'right';
            app.maxlatchingdistLabel.Position = [590 18 105 22];
            app.maxlatchingdistLabel.Text = 'max latching dist =';

            % Create max_latching_dist
            app.max_latching_dist = uieditfield(app.rounded_latch, 'numeric');
            app.max_latching_dist.Limits = [0 Inf];
            app.max_latching_dist.Tooltip = {'highest position the latch can engage'};
            app.max_latching_dist.Position = [707 18 48 22];
            app.max_latching_dist.Value = Inf;

            % Create load
            app.load = uitabgroup(app.UIFigure);
            app.load.AutoResizeChildren = 'off';
            app.load.Position = [30 611 795 123];

            % Create load_mass
            app.load_mass = uitab(app.load);
            app.load_mass.AutoResizeChildren = 'off';
            app.load_mass.Title = 'Load Mass';

            % Create RequiredFieldsLabel_6
            app.RequiredFieldsLabel_6 = uilabel(app.load_mass);
            app.RequiredFieldsLabel_6.Position = [18 65 92 22];
            app.RequiredFieldsLabel_6.Text = 'Required Fields:';

            % Create massEditFieldLabel
            app.massEditFieldLabel = uilabel(app.load_mass);
            app.massEditFieldLabel.HorizontalAlignment = 'right';
            app.massEditFieldLabel.Position = [199 65 44 22];
            app.massEditFieldLabel.Text = 'mass =';

            % Create load_mass_mass
            app.load_mass_mass = uieditfield(app.load_mass, 'numeric');
            app.load_mass_mass.Limits = [0 Inf];
            app.load_mass_mass.Tooltip = {'mass of the load'};
            app.load_mass_mass.Position = [258 65 40 22];
            app.load_mass_mass.Value = 0.001;

            % Create OptionalFieldsLabel_10
            app.OptionalFieldsLabel_10 = uilabel(app.load_mass);
            app.OptionalFieldsLabel_10.Position = [19 24 89 22];
            app.OptionalFieldsLabel_10.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_10
            app.defaultvaluesfilledinLabel_10 = uilabel(app.load_mass);
            app.defaultvaluesfilledinLabel_10.FontSize = 9;
            app.defaultvaluesfilledinLabel_10.Position = [15 13 98 13];
            app.defaultvaluesfilledinLabel_10.Text = '(default values filled in)';

            % Create massofleverarmLabel
            app.massofleverarmLabel = uilabel(app.load_mass);
            app.massofleverarmLabel.HorizontalAlignment = 'right';
            app.massofleverarmLabel.Position = [133 24 111 22];
            app.massofleverarmLabel.Text = 'mass of lever arm =';

            % Create load_m_rod
            app.load_m_rod = uieditfield(app.load_mass, 'numeric');
            app.load_m_rod.Limits = [0 Inf];
            app.load_m_rod.Tooltip = {''};
            app.load_m_rod.Position = [259 24 40 22];

            % Create EMALabel
            app.EMALabel = uilabel(app.load_mass);
            app.EMALabel.HorizontalAlignment = 'right';
            app.EMALabel.Position = [353 24 42 22];
            app.EMALabel.Text = 'EMA =';

            % Create load_EMA
            app.load_EMA = uieditfield(app.load_mass, 'numeric');
            app.load_EMA.Limits = [0 Inf];
            app.load_EMA.Tooltip = {'effective mechanical advantage'};
            app.load_EMA.Position = [407 24 48 22];
            app.load_EMA.Value = 1;

            % Create LoadingMotorLabel
            app.LoadingMotorLabel = uilabel(app.UIFigure);
            app.LoadingMotorLabel.Position = [34 306 85 22];
            app.LoadingMotorLabel.Text = 'Loading Motor:';

            % Create unlatching_motor
            app.unlatching_motor = uitabgroup(app.UIFigure);
            app.unlatching_motor.AutoResizeChildren = 'off';
            app.unlatching_motor.SelectionChangedFcn = createCallbackFcn(app, @unlatching_motorSelectionChanged, true);
            app.unlatching_motor.Position = [33 30 795 121];

            % Create um_linear_motor
            app.um_linear_motor = uitab(app.unlatching_motor);
            app.um_linear_motor.AutoResizeChildren = 'off';
            app.um_linear_motor.Title = 'Linear Motor';

            % Create FmaxLabel_5
            app.FmaxLabel_5 = uilabel(app.um_linear_motor);
            app.FmaxLabel_5.HorizontalAlignment = 'right';
            app.FmaxLabel_5.Position = [166 62 46 22];
            app.FmaxLabel_5.Text = 'Fmax =';

            % Create um_linear_motor_Fmax
            app.um_linear_motor_Fmax = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_Fmax.Limits = [0 Inf];
            app.um_linear_motor_Fmax.Tooltip = {'max force of the motor'};
            app.um_linear_motor_Fmax.Position = [227 62 39 22];
            app.um_linear_motor_Fmax.Value = 1;

            % Create VmaxLabel_3
            app.VmaxLabel_3 = uilabel(app.um_linear_motor);
            app.VmaxLabel_3.HorizontalAlignment = 'right';
            app.VmaxLabel_3.Position = [302 62 47 22];
            app.VmaxLabel_3.Text = 'Vmax =';

            % Create um_linear_motor_Vmax
            app.um_linear_motor_Vmax = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_Vmax.Limits = [0 Inf];
            app.um_linear_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.um_linear_motor_Vmax.Position = [364 62 37 22];
            app.um_linear_motor_Vmax.Value = 1;

            % Create RequiredFieldsLabel_7
            app.RequiredFieldsLabel_7 = uilabel(app.um_linear_motor);
            app.RequiredFieldsLabel_7.Position = [16 62 92 22];
            app.RequiredFieldsLabel_7.Text = 'Required Fields:';

            % Create rangeofmotionLabel_2
            app.rangeofmotionLabel_2 = uilabel(app.um_linear_motor);
            app.rangeofmotionLabel_2.HorizontalAlignment = 'right';
            app.rangeofmotionLabel_2.Position = [439 62 99 22];
            app.rangeofmotionLabel_2.Text = 'range of motion =';

            % Create um_linear_motor_range_of_motion
            app.um_linear_motor_range_of_motion = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_range_of_motion.Limits = [0 Inf];
            app.um_linear_motor_range_of_motion.Tooltip = {'total possible range of motion of the motor'};
            app.um_linear_motor_range_of_motion.Position = [553 62 47 22];
            app.um_linear_motor_range_of_motion.Value = 0.005;

            % Create voltagefracLabel_2
            app.voltagefracLabel_2 = uilabel(app.um_linear_motor);
            app.voltagefracLabel_2.HorizontalAlignment = 'right';
            app.voltagefracLabel_2.Position = [133 23 78 22];
            app.voltagefracLabel_2.Text = 'voltage frac =';

            % Create um_linear_motor_voltage_frac
            app.um_linear_motor_voltage_frac = uieditfield(app.um_linear_motor, 'numeric');
            app.um_linear_motor_voltage_frac.Limits = [0 Inf];
            app.um_linear_motor_voltage_frac.Tooltip = {'fraction of max voltage of motor'};
            app.um_linear_motor_voltage_frac.Position = [226 23 40 22];
            app.um_linear_motor_voltage_frac.Value = 1;

            % Create OptionalFieldsLabel_9
            app.OptionalFieldsLabel_9 = uilabel(app.um_linear_motor);
            app.OptionalFieldsLabel_9.Position = [15 26 89 22];
            app.OptionalFieldsLabel_9.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_9
            app.defaultvaluesfilledinLabel_9 = uilabel(app.um_linear_motor);
            app.defaultvaluesfilledinLabel_9.FontSize = 9;
            app.defaultvaluesfilledinLabel_9.Position = [11 15 98 13];
            app.defaultvaluesfilledinLabel_9.Text = '(default values filled in)';

            % Create um_hill_muscle_motor
            app.um_hill_muscle_motor = uitab(app.unlatching_motor);
            app.um_hill_muscle_motor.AutoResizeChildren = 'off';
            app.um_hill_muscle_motor.Title = 'Hill Muscle Motor';

            % Create FmaxLabel_6
            app.FmaxLabel_6 = uilabel(app.um_hill_muscle_motor);
            app.FmaxLabel_6.HorizontalAlignment = 'right';
            app.FmaxLabel_6.Position = [161 63 46 22];
            app.FmaxLabel_6.Text = 'Fmax =';

            % Create um_hill_motor_Fmax
            app.um_hill_motor_Fmax = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_Fmax.Limits = [0 Inf];
            app.um_hill_motor_Fmax.Tooltip = {'max force of the motor'};
            app.um_hill_motor_Fmax.Position = [222 63 33 22];
            app.um_hill_motor_Fmax.Value = 10;

            % Create VmaxLabel_4
            app.VmaxLabel_4 = uilabel(app.um_hill_muscle_motor);
            app.VmaxLabel_4.HorizontalAlignment = 'right';
            app.VmaxLabel_4.Position = [292 62 47 22];
            app.VmaxLabel_4.Text = 'Vmax =';

            % Create um_hill_motor_Vmax
            app.um_hill_motor_Vmax = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_Vmax.Limits = [0 Inf];
            app.um_hill_motor_Vmax.Tooltip = {'max velocity of the motor'};
            app.um_hill_motor_Vmax.Position = [354 62 37 22];
            app.um_hill_motor_Vmax.Value = 10;

            % Create musclelengthLabel_2
            app.musclelengthLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.musclelengthLabel_2.HorizontalAlignment = 'right';
            app.musclelengthLabel_2.Position = [420 62 90 22];
            app.musclelengthLabel_2.Text = 'muscle length =';

            % Create um_hill_motor_muscle_length
            app.um_hill_motor_muscle_length = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_muscle_length.Limits = [0 Inf];
            app.um_hill_motor_muscle_length.ValueChangedFcn = createCallbackFcn(app, @um_hill_motor_muscle_lengthValueChanged, true);
            app.um_hill_motor_muscle_length.Tooltip = {''};
            app.um_hill_motor_muscle_length.Position = [525 62 36 22];
            app.um_hill_motor_muscle_length.Value = 4;

            % Create RequiredFieldsLabel_8
            app.RequiredFieldsLabel_8 = uilabel(app.um_hill_muscle_motor);
            app.RequiredFieldsLabel_8.Position = [16 62 92 22];
            app.RequiredFieldsLabel_8.Text = 'Required Fields:';

            % Create OptionalFieldsLabel_6
            app.OptionalFieldsLabel_6 = uilabel(app.um_hill_muscle_motor);
            app.OptionalFieldsLabel_6.Position = [17 23 89 22];
            app.OptionalFieldsLabel_6.Text = 'Optional Fields:';

            % Create defaultvaluesfilledinLabel_6
            app.defaultvaluesfilledinLabel_6 = uilabel(app.um_hill_muscle_motor);
            app.defaultvaluesfilledinLabel_6.FontSize = 9;
            app.defaultvaluesfilledinLabel_6.Position = [13 12 98 13];
            app.defaultvaluesfilledinLabel_6.Text = '(default values filled in)';

            % Create rateofactivationLabel_2
            app.rateofactivationLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.rateofactivationLabel_2.HorizontalAlignment = 'right';
            app.rateofactivationLabel_2.Position = [589 63 104 22];
            app.rateofactivationLabel_2.Text = 'rate of activation =';

            % Create um_hill_motor_rate_of_activation
            app.um_hill_motor_rate_of_activation = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_rate_of_activation.Limits = [0 Inf];
            app.um_hill_motor_rate_of_activation.Tooltip = {''};
            app.um_hill_motor_rate_of_activation.Position = [708 63 36 22];
            app.um_hill_motor_rate_of_activation.Value = 2;

            % Create initiallengthLabel_2
            app.initiallengthLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.initiallengthLabel_2.HorizontalAlignment = 'right';
            app.initiallengthLabel_2.Position = [128 23 79 22];
            app.initiallengthLabel_2.Text = 'initial length =';

            % Create um_hill_motor_L_i
            app.um_hill_motor_L_i = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_L_i.Limits = [0 Inf];
            app.um_hill_motor_L_i.Tooltip = {''};
            app.um_hill_motor_L_i.Position = [222 23 34 22];
            app.um_hill_motor_L_i.Value = 4;

            % Create a_LLabel_2
            app.a_LLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.a_LLabel_2.HorizontalAlignment = 'right';
            app.a_LLabel_2.Position = [303 23 36 22];
            app.a_LLabel_2.Text = 'a_L =';

            % Create um_hill_motor_a_L
            app.um_hill_motor_a_L = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_a_L.Tooltip = {'determines shape of length-tension relationship'};
            app.um_hill_motor_a_L.Position = [354 23 37 22];
            app.um_hill_motor_a_L.Value = 2.08;

            % Create b_LLabel_2
            app.b_LLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.b_LLabel_2.HorizontalAlignment = 'right';
            app.b_LLabel_2.Position = [472 24 36 22];
            app.b_LLabel_2.Text = 'b_L =';

            % Create um_hill_motor_b_L
            app.um_hill_motor_b_L = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_b_L.Tooltip = {'determines shape of length-tension relationship'};
            app.um_hill_motor_b_L.Position = [523 24 39 22];
            app.um_hill_motor_b_L.Value = -2.89;

            % Create sLabel_2
            app.sLabel_2 = uilabel(app.um_hill_muscle_motor);
            app.sLabel_2.HorizontalAlignment = 'right';
            app.sLabel_2.Position = [668 23 25 22];
            app.sLabel_2.Text = 's =';

            % Create um_hill_motor_s
            app.um_hill_motor_s = uieditfield(app.um_hill_muscle_motor, 'numeric');
            app.um_hill_motor_s.Tooltip = {'determines shape of length-tension relationship'};
            app.um_hill_motor_s.Position = [708 23 37 22];
            app.um_hill_motor_s.Value = -0.75;

            % Create UnlatchingMotorLabel
            app.UnlatchingMotorLabel = uilabel(app.UIFigure);
            app.UnlatchingMotorLabel.Position = [38 151 99 22];
            app.UnlatchingMotorLabel.Text = 'Unlatching Motor:';

            % Create GraphingCornerLabel
            app.GraphingCornerLabel = uilabel(app.UIFigure);
            app.GraphingCornerLabel.FontSize = 22;
            app.GraphingCornerLabel.FontWeight = 'bold';
            app.GraphingCornerLabel.Position = [924 754 181 27];
            app.GraphingCornerLabel.Text = 'Graphing Corner';

            % Create graphing_corner
            app.graphing_corner = uitabgroup(app.UIFigure);
            app.graphing_corner.AutoResizeChildren = 'off';
            app.graphing_corner.Position = [856 238 322 496];

            % Create graphing_corner_heatmap
            app.graphing_corner_heatmap = uitab(app.graphing_corner);
            app.graphing_corner_heatmap.AutoResizeChildren = 'off';
            app.graphing_corner_heatmap.Title = '2D Plot';

            % Create minunlatchingmotorforceCheckBox
            app.minunlatchingmotorforceCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.minunlatchingmotorforceCheckBox.Text = 'min unlatching motor force';
            app.minunlatchingmotorforceCheckBox.Position = [127 97 164 22];

            % Create xaxisLabel
            app.xaxisLabel = uilabel(app.graphing_corner_heatmap);
            app.xaxisLabel.FontSize = 22;
            app.xaxisLabel.Position = [128 436 63 27];
            app.xaxisLabel.Text = 'x-axis';

            % Create yaxisLabel
            app.yaxisLabel = uilabel(app.graphing_corner_heatmap);
            app.yaxisLabel.FontSize = 22;
            app.yaxisLabel.Position = [128 313 63 27];
            app.yaxisLabel.Text = 'y-axis';

            % Create x_log_space
            app.x_log_space = uicheckbox(app.graphing_corner_heatmap);
            app.x_log_space.Tooltip = {'runs in linspace unless this box is checked'};
            app.x_log_space.Text = 'log space';
            app.x_log_space.Position = [235 364 74 22];
            app.x_log_space.Value = true;

            % Create y_log_space
            app.y_log_space = uicheckbox(app.graphing_corner_heatmap);
            app.y_log_space.Tooltip = {'runs in linspace unless this box is checked'};
            app.y_log_space.Text = 'log space';
            app.y_log_space.Position = [235 235 74 22];
            app.y_log_space.Value = true;

            % Create y_maxCheckBox
            app.y_maxCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.y_maxCheckBox.Text = 'y_max';
            app.y_maxCheckBox.Position = [32 145 58 22];
            app.y_maxCheckBox.Value = true;

            % Create y_unlatchCheckBox
            app.y_unlatchCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.y_unlatchCheckBox.Text = 'y_unlatch';
            app.y_unlatchCheckBox.Position = [32 121 74 22];

            % Create t_LCheckBox
            app.t_LCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.t_LCheckBox.Text = 't_L';
            app.t_LCheckBox.Position = [32 97 39 22];

            % Create v_toCheckBox
            app.v_toCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.v_toCheckBox.Text = 'v_to';
            app.v_toCheckBox.Position = [32 73 45 22];

            % Create P_maxCheckBox
            app.P_maxCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.P_maxCheckBox.Text = 'P_max';
            app.P_maxCheckBox.Position = [127 145 60 22];

            % Create t_toCheckBox
            app.t_toCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.t_toCheckBox.Text = 't_to';
            app.t_toCheckBox.Position = [127 121 42 22];

            % Create KE_maxCheckBox
            app.KE_maxCheckBox = uicheckbox(app.graphing_corner_heatmap);
            app.KE_maxCheckBox.Text = 'KE_max';
            app.KE_maxCheckBox.Position = [127 73 68 22];

            % Create IV1DropDownLabel
            app.IV1DropDownLabel = uilabel(app.graphing_corner_heatmap);
            app.IV1DropDownLabel.HorizontalAlignment = 'right';
            app.IV1DropDownLabel.FontSize = 16;
            app.IV1DropDownLabel.Position = [61 397 34 22];
            app.IV1DropDownLabel.Text = 'IV 1';

            % Create IV1DropDown
            app.IV1DropDown = uidropdown(app.graphing_corner_heatmap);
            app.IV1DropDown.Items = {'--no selection--'};
            app.IV1DropDown.ValueChangedFcn = createCallbackFcn(app, @IV1DropDownValueChanged, true);
            app.IV1DropDown.Tooltip = {'Independent Variable 1'};
            app.IV1DropDown.FontSize = 16;
            app.IV1DropDown.Position = [110 397 150 22];
            app.IV1DropDown.Value = '--no selection--';

            % Create IV2DropDownLabel
            app.IV2DropDownLabel = uilabel(app.graphing_corner_heatmap);
            app.IV2DropDownLabel.HorizontalAlignment = 'right';
            app.IV2DropDownLabel.FontSize = 16;
            app.IV2DropDownLabel.Position = [61 274 34 22];
            app.IV2DropDownLabel.Text = 'IV 2';

            % Create IV2DropDown
            app.IV2DropDown = uidropdown(app.graphing_corner_heatmap);
            app.IV2DropDown.Items = {'--no selection--'};
            app.IV2DropDown.ValueChangedFcn = createCallbackFcn(app, @IV2DropDownValueChanged, true);
            app.IV2DropDown.FontSize = 16;
            app.IV2DropDown.Position = [110 274 150 22];
            app.IV2DropDown.Value = '--no selection--';

            % Create xminEditFieldLabel
            app.xminEditFieldLabel = uilabel(app.graphing_corner_heatmap);
            app.xminEditFieldLabel.HorizontalAlignment = 'right';
            app.xminEditFieldLabel.Position = [19 364 34 22];
            app.xminEditFieldLabel.Text = 'x min';

            % Create xmin
            app.xmin = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.xmin.Position = [68 364 40 22];
            app.xmin.Value = 4;

            % Create xmaxEditFieldLabel
            app.xmaxEditFieldLabel = uilabel(app.graphing_corner_heatmap);
            app.xmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.xmaxEditFieldLabel.Position = [116 364 38 22];
            app.xmaxEditFieldLabel.Text = 'x max';

            % Create xmax
            app.xmax = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.xmax.Position = [169 364 40 22];
            app.xmax.Value = 9;

            % Create yminLabel
            app.yminLabel = uilabel(app.graphing_corner_heatmap);
            app.yminLabel.HorizontalAlignment = 'right';
            app.yminLabel.Position = [19 235 34 22];
            app.yminLabel.Text = 'y min';

            % Create ymin
            app.ymin = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.ymin.Position = [68 235 40 22];
            app.ymin.Value = -5;

            % Create ymaxLabel
            app.ymaxLabel = uilabel(app.graphing_corner_heatmap);
            app.ymaxLabel.HorizontalAlignment = 'right';
            app.ymaxLabel.Position = [121 235 38 22];
            app.ymaxLabel.Text = 'y max';

            % Create ymax
            app.ymax = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.ymax.Position = [174 235 40 22];
            app.ymax.Value = -2;

            % Create heatmapoutputoptionsLabel
            app.heatmapoutputoptionsLabel = uilabel(app.graphing_corner_heatmap);
            app.heatmapoutputoptionsLabel.FontSize = 22;
            app.heatmapoutputoptionsLabel.Position = [42 177 237 27];
            app.heatmapoutputoptionsLabel.Text = 'heatmap output options';

            % Create pixelsofresolutionLabel
            app.pixelsofresolutionLabel = uilabel(app.graphing_corner_heatmap);
            app.pixelsofresolutionLabel.HorizontalAlignment = 'right';
            app.pixelsofresolutionLabel.FontSize = 16;
            app.pixelsofresolutionLabel.Position = [58 30 151 22];
            app.pixelsofresolutionLabel.Text = 'pixels of resolution =';

            % Create n
            app.n = uieditfield(app.graphing_corner_heatmap, 'numeric');
            app.n.Limits = [1 1024];
            app.n.RoundFractionalValues = 'on';
            app.n.ValueDisplayFormat = '%.0f';
            app.n.Position = [224 30 39 22];
            app.n.Value = 10;

            % Create graphing_corner_one_D
            app.graphing_corner_one_D = uitab(app.graphing_corner);
            app.graphing_corner_one_D.AutoResizeChildren = 'off';
            app.graphing_corner_one_D.Title = '1D Plot';

            % Create OD_minunlatchingmotorforceCheckBox
            app.OD_minunlatchingmotorforceCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_minunlatchingmotorforceCheckBox.Text = 'min unlatching motor force';
            app.OD_minunlatchingmotorforceCheckBox.Position = [127 228 164 22];

            % Create xaxisLabel_2
            app.xaxisLabel_2 = uilabel(app.graphing_corner_one_D);
            app.xaxisLabel_2.FontSize = 22;
            app.xaxisLabel_2.Position = [128 436 63 27];
            app.xaxisLabel_2.Text = 'x-axis';

            % Create OD_x_log_space
            app.OD_x_log_space = uicheckbox(app.graphing_corner_one_D);
            app.OD_x_log_space.Tooltip = {'runs in linspace unless this box is checked'};
            app.OD_x_log_space.Text = 'log space';
            app.OD_x_log_space.Position = [235 364 74 22];
            app.OD_x_log_space.Value = true;

            % Create OD_y_maxCheckBox
            app.OD_y_maxCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_y_maxCheckBox.Text = 'y_max';
            app.OD_y_maxCheckBox.Position = [32 276 58 22];
            app.OD_y_maxCheckBox.Value = true;

            % Create OD_y_unlatchCheckBox
            app.OD_y_unlatchCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_y_unlatchCheckBox.Text = 'y_unlatch';
            app.OD_y_unlatchCheckBox.Position = [32 252 74 22];

            % Create OD_t_LCheckBox
            app.OD_t_LCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_t_LCheckBox.Text = 't_L';
            app.OD_t_LCheckBox.Position = [32 228 39 22];

            % Create OD_v_toCheckBox
            app.OD_v_toCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_v_toCheckBox.Text = 'v_to';
            app.OD_v_toCheckBox.Position = [32 204 45 22];

            % Create OD_P_maxCheckBox
            app.OD_P_maxCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_P_maxCheckBox.Text = 'P_max';
            app.OD_P_maxCheckBox.Position = [127 276 60 22];

            % Create OD_t_toCheckBox
            app.OD_t_toCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_t_toCheckBox.Text = 't_to';
            app.OD_t_toCheckBox.Position = [127 252 42 22];

            % Create OD_KE_maxCheckBox
            app.OD_KE_maxCheckBox = uicheckbox(app.graphing_corner_one_D);
            app.OD_KE_maxCheckBox.Text = 'KE_max';
            app.OD_KE_maxCheckBox.Position = [127 204 68 22];

            % Create IV1DropDownLabel_2
            app.IV1DropDownLabel_2 = uilabel(app.graphing_corner_one_D);
            app.IV1DropDownLabel_2.HorizontalAlignment = 'right';
            app.IV1DropDownLabel_2.FontSize = 16;
            app.IV1DropDownLabel_2.Position = [61 397 34 22];
            app.IV1DropDownLabel_2.Text = 'IV 1';

            % Create OD_IV1DropDown
            app.OD_IV1DropDown = uidropdown(app.graphing_corner_one_D);
            app.OD_IV1DropDown.Items = {'--no selection--'};
            app.OD_IV1DropDown.Tooltip = {'Independent Variable 1'};
            app.OD_IV1DropDown.FontSize = 16;
            app.OD_IV1DropDown.Position = [110 397 150 22];
            app.OD_IV1DropDown.Value = '--no selection--';

            % Create yaxisoutputoptionsLabel
            app.yaxisoutputoptionsLabel = uilabel(app.graphing_corner_one_D);
            app.yaxisoutputoptionsLabel.FontSize = 22;
            app.yaxisoutputoptionsLabel.Position = [60 313 208 27];
            app.yaxisoutputoptionsLabel.Text = 'y-axis output options';

            % Create pixelsofresolutionLabel_2
            app.pixelsofresolutionLabel_2 = uilabel(app.graphing_corner_one_D);
            app.pixelsofresolutionLabel_2.HorizontalAlignment = 'right';
            app.pixelsofresolutionLabel_2.FontSize = 16;
            app.pixelsofresolutionLabel_2.Position = [58 161 151 22];
            app.pixelsofresolutionLabel_2.Text = 'pixels of resolution =';

            % Create OD_n
            app.OD_n = uieditfield(app.graphing_corner_one_D, 'numeric');
            app.OD_n.Limits = [1 1024];
            app.OD_n.RoundFractionalValues = 'on';
            app.OD_n.ValueDisplayFormat = '%.0f';
            app.OD_n.Position = [224 161 39 22];
            app.OD_n.Value = 10;

            % Create xminEditFieldLabel_2
            app.xminEditFieldLabel_2 = uilabel(app.graphing_corner_one_D);
            app.xminEditFieldLabel_2.HorizontalAlignment = 'right';
            app.xminEditFieldLabel_2.Position = [19 364 34 22];
            app.xminEditFieldLabel_2.Text = 'x min';

            % Create OD_xmin
            app.OD_xmin = uieditfield(app.graphing_corner_one_D, 'numeric');
            app.OD_xmin.Position = [68 364 40 22];
            app.OD_xmin.Value = 4;

            % Create xmaxEditFieldLabel_2
            app.xmaxEditFieldLabel_2 = uilabel(app.graphing_corner_one_D);
            app.xmaxEditFieldLabel_2.HorizontalAlignment = 'right';
            app.xmaxEditFieldLabel_2.Position = [116 364 38 22];
            app.xmaxEditFieldLabel_2.Text = 'x max';

            % Create OD_xmax
            app.OD_xmax = uieditfield(app.graphing_corner_one_D, 'numeric');
            app.OD_xmax.Position = [169 364 40 22];
            app.OD_xmax.Value = 9;

            % Create graphing_corner_kinematics
            app.graphing_corner_kinematics = uitab(app.graphing_corner);
            app.graphing_corner_kinematics.Title = 'Kinematics';

            % Create go
            app.go = uibutton(app.UIFigure, 'push');
            app.go.ButtonPushedFcn = createCallbackFcn(app, @goButtonPushed, true);
            app.go.FontSize = 22;
            app.go.Position = [945 133 150 76];
            app.go.Text = 'Graph!';

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

