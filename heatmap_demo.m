%% run this file to generate heatplots for the 6 performance metrics

%% don't touch these
close all
clearvars
tic
debug = false;
addpath(genpath(pwd)); % add all subdirectories to path to access the files in components-library

%% plot parameters 
N=25; % determines resolution of heatplots

% setting x axis on the plot (Fmax of latch)
xname = 'Fmax';
xrange = [-1 3];
Fmaxs = logspace(xrange(1),xrange(2),N);

%setting y axis value on plot (Vmax of latch)
yname = 'vmax';
yrange = [-2 2];
v_maxs = logspace(yrange(1),yrange(2),N);

metrics = {'tto','vto','Pmax','ymax','tL','KEmax','yunlatch'};

load_time_constraint = Inf;

%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 20;
loading_motor_range_of_motion = 3;
v_max_loading_motor = 10.0000;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 10;
loading_motor_r_activation = Inf;

% loading motor struct initialization
%loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
m=100;

% load mass struct initialization
load = load_mass(m);

%% latch

% latch parameters
R=2E-1;
m_L= 100;
coeff_fric = 0;
v_0L=0;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% spring paramters
k = 6; % k or k_0 depending on linear or exponential spring
m_s=1;
F_spring_max=1E4;

% extra parameters for exponential spring
characteristic_length = 5;

% spring struct initialization
spring = linear_spring(k, m_s, F_spring_max);
%spring = exponential_spring(k, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
unlatchinging_motor_range_of_motion = 3;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 10;
unlatching_motor_r_activation = Inf;

% unlatching motor struct initialization happens in next section

%% initializing an output value matrix for each metric

for ii=1:length(metrics)
    outval{ii}=zeros(N);
end
if (debug)
    h1 = figure()
    h2 = figure()
end
for i = 1:N % iterate over y-axis-variable of plot
    for j = 1:N % iterate over x-axis-variable of plot
        
        % unlatching motor struct initialization
        unlatching_motor = hill_muscle_motor(unlatching_motor_muscle_length, Fmaxs(j), v_maxs(i), unlatching_motor_r_activation);
        %unlatching_motor = linear_motor(Fmaxs(j),v_maxs(i), unlatching_motor_range_of_motion);
        
        % input structs for each component of LaMSA system into solve_model
        [sol,transition_times]=solve_model(loading_motor, unlatching_motor, load, latch, spring);

        if (debug)
            figure(h1)
            plot(sol(:,1),sol(:,2),'.');
            hold on;
            figure(h2)
            plot(sol(:,1),sol(:,3),'.');
            hold on;
            ginput(1)
        end
        met_dict=get_metrics(sol,transition_times,load, spring ,metrics);
        for ii=1:length(metrics)
            outval{ii}(i,j)=met_dict(metrics{ii});
        end
         
    end
   disp(['row ' num2str(i) ' of ' num2str(N)]);
end
toc

%% Plot the output data

for ii=1:length(metrics)
    figure();
    imagesc(xrange,yrange,outval{ii});
    set(gca,'YDir','normal')
    xlabel(xname)
    ylabel(yname);
    c = colorbar;
    c.Label.String = metrics{ii};
end