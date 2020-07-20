%% run this file to generate heatplots for the 6 performance metrics

%% don't touch these
close all
clearvars
tic
debug = false;
addpath(genpath(fullfile(pwd,'..'))); % add all subdirectories to path to access the files in components-library

%% edit the following parameters

%% plot parameters 
N=256; % determines resolution of heatplots

% setting x axis on the plot (Fmax of latch)
xname = 'k value [N/m]';
xrange = [100 10000];
k_val = linspace(xrange(1),xrange(2),N);

%setting y axis value on plot (Vmax of latch)
yname = 'muscle max force [N]';
yrange = [1E-1 10];
muscle_max_f = linspace(yrange(1),yrange(2),N);

metrics = {'tto','vto','Pmax','KEmax'};


load_time_constraint = Inf;

%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 4;
loading_motor_range_of_motion = 5;
v_max_loading_motor = 1;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 10E-3;
loading_motor_r_activation = Inf;

% loading motor struct initialization
%loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
m=25E-3;

%trap jaw ant load mass
% L = 1E-2;
% rho = 10;
% A = 1E-2;
% E = 10;
% sigma_f = 10E6;
% where m_s = L*rho*A, k=(E*A)/L


EMA = 1;
m_rod = .1;
m_end = .1;

% load mass struct initialization
% load = load_mass(m_end,m_rod,EMA);
load = load_mass(m);


%% latch

% latch parameters
R=5E-3;
m_L= 3E-3;

coeff_fric = 0;
v_0L=0;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% spring paramters
k = 1; % k or k_0 depending on linear or exponential spring
m_s=2E-3;
F_spring_max=20;

% extra parameters for exponential spring
% should be a negative value
characteristic_length = 1E-3;

% spring struct initialization
spring = linear_spring(k, m_s, F_spring_max);
%spring = exponential_spring(k, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
F_max_unlatching_motor = .25;
unlatching_motor_range_of_motion = 5E-3;
v_max_unlatching_motor=1;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 10;
unlatching_motor_r_activation = Inf;
unlatching_motor= linear_motor(F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_range_of_motion);
% unlatching motor struct initialization happens in next section
%% end editable parameters

% make a directory for every run
output_directory = create_output_directory();
        


%% initializing an output value matrix for each metric
% Establishing a cell array of the output matrices with resolution N
for ii=1:length(metrics)
    outval{ii}=zeros(N);
end
if (debug)
    h1 = figure()
    h2 = figure()
end
for i=1:N %iterate over y-axis-variable of plot
    for j=1:N %iterate over x-axis-variable of plot
        % unlatching motor struct initialization
        %unlatching_motor = hill_muscle_motor(unlatching_motor_muscle_length, Fmaxs(j), v_maxs(i),unlatching_motor_r_activation);
        %input structs for each component
%         k = Es(j)*As(i)/L;
%         F_spring_max= sigma_f*As(i);
%         m_s=As(i)*L*rho;
        loading_motor=hill_muscle_motor(loading_motor_muscle_length, muscle_max_f(i), v_max_loading_motor, loading_motor_r_activation);
        spring=exponential_spring(k_val(j),characteristic_length, m_s, F_spring_max);
        % input structs for each component of LaMSA system into solve_model
        [sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring);%Add 6th input to write csv files with data: output_directory
        [solDA, ttDA] = solve_direct_actuation(loading_motor, load);
        if (debug)
            figure(h1)
            plot(sol(:,1),sol(:,2),'.');
            hold on;
            figure(h2)
            plot(sol(:,1),sol(:,3),'.');
            hold on;
            ginput(1)
        end
        met_dict=get_metrics(sol,transition_times,load ,metrics);%getting metrics for LaMSA system
        met_dict_DA = get_metrics(solDA, ttDA, load, metrics);%getting metrics for direct actuation system
        for ii=1:length(metrics)
            outval{ii}(i,j)=(met_dict(metrics{ii}))/(met_dict_DA(metrics{ii}));%assigning the values in the output matrix as the ratio between the systems
        end
        
    end
   disp(['row ' num2str(i) ' of ' num2str(N)]);
end
toc

%% Plot the output data
lPlot = figure;
n=1;
for ii=1:length(metrics)
%Establishing axes and labels
    name = {'ax1' 'ax2' 'ax3' 'ax4'};
    name{ii} = subplot(2,2,n);
    imagesc(xrange,yrange,outval{ii});
    set(gca,'YDir','normal');
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xname,'Interpreter', 'Latex');
    ylabel(yname, 'Interpreter', 'Latex');
%Colorbar settings
    c = colorbar;
    map = linspecer(N);
    colormap(name{ii}, map);
    c.Label.String = metrics{ii};
    set(c,'TickLabelInterpreter','latex')
    c.Label.Interpreter="latex";
    n=n+1;
%Finding the LaMSA Zone Boundary
    hold on
    bndry = zeros(size(outval{ii}) + [2,2]);
    bndry(2:end-1, 2:end-1) = outval{ii};
    bndry(1,:) = bndry(2,:);
    bndry(end,:) = bndry(end-1,:);
    bndry(:,1) = bndry(:,2);
    bndry(:,end) = bndry(:,end-1);
    LaMSA = zeros(size(bndry));
    index1 = find(bndry > 1);
    LaMSA(index1) = 1;
    B = bwboundaries(LaMSA);
    for k = 1:length(B)
        boundary = B{k}-1;
        out_of_range = boundary<1|boundary>N;
        out_of_range = out_of_range(:,1)|out_of_range(:,2);
        boundary(out_of_range,:)=[];
%Make sure to change the plot variables to match the x and y variable
%spaces in the next line
        plot(k_val(boundary(:,2)), muscle_max_f(boundary(:,1)), '.k', 'LineWidth', 3)
        hold on
    end
end

prompt = sprintf('Please enter a name for the LaMSA Zone plot:\n');
str = input(prompt, 's');
saveas(lPlot, strcat(str, '.fig'));

