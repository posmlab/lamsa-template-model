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

% setting x axis on the plot 
xname = 'mandible mass';
xrange = [-3 -1];
%Fmaxs = logspace(xrange(1),xrange(2),N);
m_m = logspace(xrange(1),xrange(2),N);


%setting y axis value on plot 
yname = 'ema';
yrange = [-2 0];
ema = logspace(yrange(1),yrange(2),N);

metrics = {'tto','vto','Pmax','ymax','tL','KEmax','yunlatch','amax'};
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
F_max_loading_motor = 4.72E-6;
loading_motor_range_of_motion = 3;
v_max_loading_motor = 10.0000;

% extra parameters for hill muscle motor
loading_motor_muscle_length = 1.5*1;
loading_motor_r_activation = Inf;

% loading motor struct initialization
%loading_motor = linear_motor(F_max_loading_motor, v_max_loading_motor, loading_motor_range_of_motion);
loading_motor = hill_muscle_motor(loading_motor_muscle_length, F_max_loading_motor, v_max_loading_motor, loading_motor_r_activation);

%% load mass

% load mass parameters
%m=1;

%trap jaw ant load mass
EMA = 1.16E-1;
m_rod = 1.56E-2;
m_end = 0;

% load mass struct initialization
load = load_mass(m_end,m_rod,EMA);


%% latch

% latch parameters
R=3.96E-1;
m_L= 1;

coeff_fric = 0;
v_0L=1;

% latch struct initialization
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%% spring

% extra parameters for exponential spring
% should be a negative value
characteristic_length = 1;

%trap jaw ant spring parameters 
 L = 1;
 rho = 1.59E3;
 A = 2.56E-3;
 E = 1;
 sigma_f = 3.18E-1;
 m_s = L*rho*A;
 k=(E*A)/L;
 F_spring_max= sigma_f*A;
 

% spring struct initialization
spring = linear_spring(k, m_s, F_spring_max);
%spring = exponential_spring(k, characteristic_length, m_s, F_spring_max);

%% unlatching motor

% unlatching motor paramters for linear motor
F_max_unlatching_motor = 0;
unlatchinging_motor_range_of_motion = Inf;

% extra parameters for hill muscle motor
unlatching_motor_muscle_length = 0;
unlatching_motor_r_activation = Inf;

% unlatching motor struct initialization happens in next section
unlatching_motor = linear_motor(0,0,0);
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

        load = load_mass(0,m_m(j),ema(i));
        [sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring);
        if (debug)
            figure(h1)
            plot(sol(:,1),sol(:,2),'.');
            hold on;
            figure(h2)
            plot(sol(:,1),sol(:,3),'.');
            hold on;
            ginput(1)
        end
        met_dict=get_metrics(sol,transition_times,load,metrics);
        for ii=1:length(metrics)
            outval{ii}(i,j)=met_dict(metrics{ii});
        end
         
    end
   disp(['row ' num2str(i) ' of ' num2str(N)]);
end
toc

%% Plot the output data
hmap = figure;
n=1;
for ii=1:length(metrics)
    subplot(2,4,n);
    imagesc(xrange,yrange,outval{ii});
    set(gca,'YDir','normal');
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xname,'Interpreter', 'Latex');
    ylabel(yname, 'Interpreter', 'Latex');
    c = colorbar;
    c.Label.String = metrics{ii};
    set(c,'TickLabelInterpreter','latex')
    c.Label.Interpreter="latex";
    n=n+1;
end
prompt = sprintf('Please enter a name for the heatmap figure:\n');
str = input(prompt, 's');
saveas(hmap, strcat(str, '.fig'));

%%Comparison
%big_Diff=max(max(outval{1}-old))