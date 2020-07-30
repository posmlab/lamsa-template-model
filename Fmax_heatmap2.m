clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines resolution of heatplots
N=10; 

% setting x axis on the plot 
xname = 'loading motor Fmax';
xrange = [0 40];
looping_value_x = linspace(xrange(1),xrange(2),N);

%setting y axis value on plot 
yname = 'standard deviation';
yrange = [0 1];
looping_value_y = linspace(yrange(1),yrange(2),N);

krange = [0 1E4];
k_values = linspace(krange(1),krange(2),N);

% add things to metrics
metrics = {'Pmax'};

% load mass struct initialization
load = load_mass(0.01,0,1);

% latch struct initialization
latch = rounded_latch(0.005,0.003,0,0,0.003,0.0035,0);
%latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);

% spring struct initialization

%spring = exponential_spring(app.exp_spring_k.Value,app.exp_spring_char_len.Value,app.exp_spring_mass.Value,app.exp_spring_Fmax.Value);

% loading motor struct initialization
%loading_motor = hill_muscle_motor(app.lm_hill_motor_muscle_length.Value,app.lm_hill_motor_Fmax.Value,app.lm_hill_motor_Vmax.Value,app.lm_hill_motor_rate_of_activation.Value,app.lm_hill_motor_L_i.Value,app.lm_hill_motor_a_L.Value,app.lm_hill_motor_b_L.Value,app.lm_hill_motor_s.Value);

% unlatching motor struct initialization
unlatching_motor = linear_motor(0.25,1,0.005,1);
%unlatching_motor = hill_muscle_motor(app.um_hill_motor_muscle_length.Value,app.um_hill_motor_Fmax.Value,app.um_hill_motor_Vmax.Value,app.um_hill_motor_rate_of_activation.Value,app.um_hill_motor_L_i.Value,app.um_hill_motor_a_L.Value,app.um_hill_motor_b_L.Value,app.um_hill_motor_s.Value);

disp('starting')
highest_P = zeros(N);
best_k = zeros(1,N);
for jj = 1:length(k_values)
    spring = linear_spring(k_values(jj),0,Inf);
    
    for ii=1:length(metrics)
        outval{ii}=zeros(N);
    end
    for ii=1:length(metrics)
        outval1{ii} = zeros(size(looping_value_x));
    end
    for j = 1:N %x var
        loading_motor = linear_motor(looping_value_x(j),10,0.005,1);

        %calling solve_model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,metrics);
        for ii=1:length(metrics)
          outval1{ii}(j)=met_dict(metrics{ii});
        end
    end
    Pmax = outval1{1};
    Fmax = looping_value_x;
    outval{1}(1,:) = Pmax;
    for i = 2:N %y var
        Pmax_sigma = zeros(size(Pmax));
        for j = 1:N %x var        
            mu = looping_value_x(j);
            sigma = looping_value_y(i)*mu;
            if (sigma < (Fmax(2)-Fmax(1)))
                Pmax_sigma(j) = Pmax(j);
            else
                gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
                Pmax_sigma(j) = trapz(Fmax,gaus(Fmax).*Pmax);
            end
            [maxval, index] = max(Pmax_sigma);
            outval{ii}(i,j) = maxval;
        end
    end

    for x = 1:N
        if (outval{1}(1,x) > highest_P(1,x))
            highest_P(1,x) = outval{1}(1,x);
            best_k(x) = k_values(jj);
        end
    end
    
    disp(['row ' num2str(jj) ' of ' num2str(N)]);
end
%% calculating new max power
for ii=1:length(metrics)
    outval{ii}=zeros(N);
end
for ii=1:length(metrics)
    outval1{ii} = zeros(size(looping_value_x));
end
for j = 1:N %x var
    spring = linear_spring(best_k(j),0,Inf);
    loading_motor = linear_motor(looping_value_x(j),10,0.005,1);

    %calling solve_model
    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

    met_dict=get_metrics(sol,transition_times,load,metrics);
    for ii=1:length(metrics)
      outval1{ii}(j)=met_dict(metrics{ii});
    end
end
Pmax = outval1{1};
Fmax = looping_value_x;
outval{1}(1,:) = Pmax;
for i = 2:N %y var
    Pmax_sigma = zeros(size(Pmax));
    for j = 1:N %x var        
        mu = looping_value_x(j);
        sigma = looping_value_y(i)*mu;
        if (sigma < (Fmax(2)-Fmax(1)))
            Pmax_sigma(j) = Pmax(j);
        else
            gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
            Pmax_sigma(j) = trapz(Fmax,gaus(Fmax).*Pmax);
        end
        [maxval, index] = max(Pmax_sigma);
        outval{ii}(i,j) = maxval;
    end
    disp(['row ' num2str(j) ' of ' num2str(N)]);
end
    
%% plotting
hmap = figure;
n=1;
for ii=1:length(metrics)
%Establishing axes and labels
    imagesc(xrange,yrange,repmat(best_k,[N,1]));
    set(gca,'YDir','normal');
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xname,'Interpreter', 'Latex');
    ylabel(yname, 'Interpreter', 'Latex');
%Colorbar settings
    c = colorbar;
    c.Label.String = 'k';
    set(c,'TickLabelInterpreter','latex')
    c.Label.Interpreter="latex";
    n=n+1;
end

hmap2 = figure;
n=1;
for ii=1:length(metrics)
%Establishing axes and labels
    imagesc(xrange,yrange,outval{1});
    set(gca,'YDir','normal');
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xname,'Interpreter', 'Latex');
    ylabel(yname, 'Interpreter', 'Latex');
%Colorbar settings
    c = colorbar;
    c.Label.String = 'highest P with best k';
    set(c,'TickLabelInterpreter','latex')
    c.Label.Interpreter="latex";
    n=n+1;
end