addpath(genpath(fullfile(pwd,'..')));

% determines resolution of heatplots
N=1000; 

% setting x axis on the plot
xrange = [1 100000];
looping_value_x = linspace(xrange(1),xrange(2),N);

% add things to metrics
metrics = {'Pmax'};

% load mass struct initialization
load = load_mass(0.01,0,1);

% latch struct initialization
%latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);


% spring struct initialization
%spring = exponential_spring(app.exp_spring_k.Value,app.exp_spring_char_len.Value,app.exp_spring_mass.Value,app.exp_spring_Fmax.Value);

% loading motor struct initialization
loading_motor = linear_motor(10,10,0.005,1);
%loading_motor = hill_muscle_motor(app.lm_hill_motor_muscle_length.Value,app.lm_hill_motor_Fmax.Value,app.lm_hill_motor_Vmax.Value,app.lm_hill_motor_rate_of_activation.Value,app.lm_hill_motor_L_i.Value,app.lm_hill_motor_a_L.Value,app.lm_hill_motor_b_L.Value,app.lm_hill_motor_s.Value);

% unlatching motor struct initialization
unlatching_motor = linear_motor(0.25,1,0.005,1);
%unlatching_motor = hill_muscle_motor(app.um_hill_motor_muscle_length.Value,app.um_hill_motor_Fmax.Value,app.um_hill_motor_Vmax.Value,app.um_hill_motor_rate_of_activation.Value,app.um_hill_motor_L_i.Value,app.um_hill_motor_a_L.Value,app.um_hill_motor_b_L.Value,app.um_hill_motor_s.Value);

k_index = 1;
M=20;
kopt_values = zeros(1,M);
range = [0.0002 0.0048];
linsp = linspace(range(1), range(2), M);
for mld = 1:length(linsp)
    latch = rounded_latch(0.005,0.003,0,0,0,linsp(mld),0);

    for ii=1:length(metrics)
        outval{ii} = zeros(size(looping_value_x));
    end
    for i=1:N    
        spring = linear_spring(looping_value_x(i),0,Inf);

        % calling solve model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,metrics);
        for ii=1:length(metrics)
          outval{ii}(i)=met_dict(metrics{ii});
        end

        %disp(['A point ' num2str(i) ' of ' num2str(N)]);
    end

    % plot output
%     fh = figure('Name','1D Plot');
%     fh.WindowState = 'maximized';
%     subplot_rows = floor(sqrt(length(metrics)));
%     subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
%     for ii=1:length(metrics)
%         subplot(subplot_rows,subplot_cols,ii);
%         plot(looping_value_x,outval{ii},'.'); 
%         xlabel('k','Interpreter', 'Latex');
%         ylabel('Pmax', 'Interpreter', 'Latex');
%         ylim_range = max(outval{ii})-min(outval{ii});
%         if (ylim_range == 0)
%             ylim_range = abs(max(outval{ii}));
%             if ~(ylim_range == 0)
%                 ylim([min([0, min(outval{ii})-0.05*ylim_range]) max([0,max(outval{ii})+0.05*ylim_range]) ])
%             end
%         else
%             ylim([min(outval{ii})-0.05*ylim_range max(outval{ii})+0.05*ylim_range])
%         end
%         set(gca,'TickLabelInterpreter','latex')
%         % makes the x in log scale if needed
% %         if app.OD_x_log_space.Value
% %             set(gca,'XScale','log')
% %         end
% 
%     end
%  
%      hold on

    Pmax = outval{1};
    k = looping_value_x;
    Pmax_sigma5 = zeros(size(Pmax));
    for i = 1:length(k)
        mu = looping_value_x(i);
        sigma = 0.1*mu;
        gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
        Pmax_sigma5(i) = trapz(k,gaus(k).*Pmax);
    end

    %plot(k,Pmax_sigma5,'ok')
    
    [maxval, index] = max(Pmax_sigma5);
    kopt_values(k_index) = looping_value_x(index);
    k_index = k_index + 1;
    
    disp(['point ' num2str(mld) ' of ' num2str(M)]);
end
figure
plot(linsp,kopt_values,'ok');
xlabel('max latching distance');
ylabel('k_optimal');
hold on

%% second sigma
kopt_values2 = zeros(1,M);
k_index = 1;
for mld = 1:length(linsp)
    latch = rounded_latch(0.005,0.003,0,0,0,linsp(mld),0);

    for ii=1:length(metrics)
        outval{ii} = zeros(size(looping_value_x));
    end
    for i=1:N    
        spring = linear_spring(looping_value_x(i),0,Inf);

        % calling solve model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,metrics);
        for ii=1:length(metrics)
          outval{ii}(i)=met_dict(metrics{ii});
        end

        %disp(['A point ' num2str(i) ' of ' num2str(N)]);
    end

    % plot output
%     fh = figure('Name','1D Plot');
%     fh.WindowState = 'maximized';
%     subplot_rows = floor(sqrt(length(metrics)));
%     subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
%     for ii=1:length(metrics)
%         subplot(subplot_rows,subplot_cols,ii);
%         plot(looping_value_x,outval{ii},'.'); 
%         xlabel('k','Interpreter', 'Latex');
%         ylabel('Pmax', 'Interpreter', 'Latex');
%         ylim_range = max(outval{ii})-min(outval{ii});
%         if (ylim_range == 0)
%             ylim_range = abs(max(outval{ii}));
%             if ~(ylim_range == 0)
%                 ylim([min([0, min(outval{ii})-0.05*ylim_range]) max([0,max(outval{ii})+0.05*ylim_range]) ])
%             end
%         else
%             ylim([min(outval{ii})-0.05*ylim_range max(outval{ii})+0.05*ylim_range])
%         end
%         set(gca,'TickLabelInterpreter','latex')
%         % makes the x in log scale if needed
% %         if app.OD_x_log_space.Value
% %             set(gca,'XScale','log')
% %         end
% 
%     end
%  
%      hold on

    Pmax = outval{1};
    k = looping_value_x;
    Pmax_sigma5 = zeros(size(Pmax));
    for i = 1:length(k)
        mu = looping_value_x(i);
        sigma = 0.05*mu;
        gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
        Pmax_sigma5(i) = trapz(k,gaus(k).*Pmax);
    end

    %plot(k,Pmax_sigma5,'ok')
    
    [maxval, index] = max(Pmax_sigma5);
    kopt_values2(k_index) = looping_value_x(index);
    k_index = k_index + 1;
    
    disp(['point ' num2str(mld) ' of ' num2str(M)]);
end

plot(linsp,kopt_values2,'or');

%% third sigma

kopt_values3 = zeros(1,M);
k_index = 1;
for mld = 1:length(linsp)
    latch = rounded_latch(0.005,0.003,0,0,0,linsp(mld),0);

    for ii=1:length(metrics)
        outval{ii} = zeros(size(looping_value_x));
    end
    for i=1:N    
        spring = linear_spring(looping_value_x(i),0,Inf);

        % calling solve model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,metrics);
        for ii=1:length(metrics)
          outval{ii}(i)=met_dict(metrics{ii});
        end

        %disp(['A point ' num2str(i) ' of ' num2str(N)]);
    end

    % plot output
%     fh = figure('Name','1D Plot');
%     fh.WindowState = 'maximized';
%     subplot_rows = floor(sqrt(length(metrics)));
%     subplot_cols = ceil(length(metrics)/floor(sqrt(length(metrics))));
%     for ii=1:length(metrics)
%         subplot(subplot_rows,subplot_cols,ii);
%         plot(looping_value_x,outval{ii},'.'); 
%         xlabel('k','Interpreter', 'Latex');
%         ylabel('Pmax', 'Interpreter', 'Latex');
%         ylim_range = max(outval{ii})-min(outval{ii});
%         if (ylim_range == 0)
%             ylim_range = abs(max(outval{ii}));
%             if ~(ylim_range == 0)
%                 ylim([min([0, min(outval{ii})-0.05*ylim_range]) max([0,max(outval{ii})+0.05*ylim_range]) ])
%             end
%         else
%             ylim([min(outval{ii})-0.05*ylim_range max(outval{ii})+0.05*ylim_range])
%         end
%         set(gca,'TickLabelInterpreter','latex')
%         % makes the x in log scale if needed
% %         if app.OD_x_log_space.Value
% %             set(gca,'XScale','log')
% %         end
% 
%     end
%  
%      hold on

    Pmax = outval{1};
    k = looping_value_x;
    Pmax_sigma5 = zeros(size(Pmax));
    for i = 1:length(k)
        mu = looping_value_x(i);
        sigma = 0.2*mu;
        gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
        Pmax_sigma5(i) = trapz(k,gaus(k).*Pmax);
    end

    %plot(k,Pmax_sigma5,'ok')
    
    [maxval, index] = max(Pmax_sigma5);
    kopt_values3(k_index) = looping_value_x(index);
    k_index = k_index + 1;
    
    disp(['point ' num2str(mld) ' of ' num2str(M)]);
end

plot(linsp,kopt_values3,'ob');