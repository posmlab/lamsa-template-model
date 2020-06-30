%% run this file to generate heatplots for the 6 performance metrics

%% don't touch these
close all
clearvars
tic
debug = false;
addpath(genpath(fullfile(pwd,'..'))); % add all subdirectories to path to access the files in components-library

%% edit the following parameters

%% plot parameters 
N=20; % determines resolution of heatplots

% setting x axis on the plot (Fmax of latch)
xname = 'k value [N/m]';
xrange = [-1 3];
k_val = logspace(xrange(1),xrange(2),N);

%setting y axis value on plot (Vmax of latch)
yname = 'muscle max force [N]';
yrange = [-1 2];
muscle_max_f = logspace(yrange(1),yrange(2),N);

metrics = {'tto','vto','Pmax','KEmax'};
% hold on
% close all
% clearvars
% tic
% debug = false;
% N=50;


% xname = 'Young''s Modulus, $E$[Pa]';
% xrange = [2 9];
% Es=logspace(xrange(1),xrange(2),N);
% yname = 'Cross-Sectional Area, $A$[m$^2$]';
% yrange = [-10 -0];
% As=logspace(yrange(1),yrange(2),N);
% metrics = {'KEmax'};
% L = 10E-3;
% rho = 10;
% sigma_f = 10E6;

load_time_constraint = Inf;

%% loading motor

% loading motor parameters for linear motor
F_max_loading_motor = 20;
loading_motor_range_of_motion = 5;
v_max_loading_motor = 100.0000;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 10;
loading_motor_r_activation = Inf;

% loading motor struct initialization
%loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
m=.001;

%trap jaw ant load mass
% L = 1E-2;
% rho = 10;
% A = 1E-2;
% E = 10;
% where m_s = L*rho*A, k=(E*A)/L


EMA = 1;
m_rod = .1;
m_end = .1;

% load mass struct initialization
load = load_mass(m_end,m_rod,EMA);


%% latch

% latch parameters
R=.01;
m_L= 1000;

coeff_fric = 0;
v_0L=0;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% spring paramters
k = 1; % k or k_0 depending on linear or exponential spring
m_s=0;
F_spring_max=1E4;

% extra parameters for exponential spring
% should be a negative value
characteristic_length = 5;

% spring struct initialization
spring = linear_spring(k, m_s, F_spring_max);
%spring = exponential_spring(k, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
F_max_unlatching_motor = 100;
unlatching_motor_range_of_motion = 50;
v_max_unlatching_motor=50;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 10;
unlatching_motor_r_activation = Inf;
unlatching_motor= linear_motor(F_max_unlatching_motor, v_max_unlatching_motor, unlatching_motor_range_of_motion);
% unlatching motor struct initialization happens in next section
%% end editable parameters

% make a directory for every run
output_directory = create_output_directory();
        


%% initializing an output value matrix for each metric

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
        [sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring, output_directory);
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
        met_dict=get_metrics(sol,transition_times,load ,metrics);
        met_dict_DA = get_metrics(solDA, ttDA, load, metrics);
        for ii=1:length(metrics)
            outval{ii}(i,j)=(met_dict(metrics{ii}))/(met_dict_DA(metrics{ii}));
        end
        
    end
   disp(['row ' num2str(i) ' of ' num2str(N)]);
end
toc

%% Plot the output data
figure
n=1;
for ii=1:length(metrics)
    name = {'ax1' 'ax2' 'ax3' 'ax4'};
    name{ii} = subplot(2,2,n);
    imagesc(xrange,yrange,outval{ii});
    set(gca,'YDir','normal');
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xname,'Interpreter', 'Latex');
    ylabel(yname, 'Interpreter', 'Latex');
    c = colorbar;
    map = linspecer(N);
    colormap(name{ii}, map);
    c.Label.String = metrics{ii};
    set(c,'TickLabelInterpreter','latex')
    c.Label.Interpreter="latex";
    n=n+1;
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
        plot(log10(k_val(boundary(:,2))), log10(muscle_max_f(boundary(:,1))), '.k', 'LineWidth', 3)
    end
end

%%Comparison
%big_Diff=max(max(outval{1}-old))