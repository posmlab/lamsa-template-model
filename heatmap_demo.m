%run this file to generate heatplots for the 6 performance metrics
%vary parameters below 
close all
clearvars
tic
debug = false;
N=50;



%setting x axis on the plot (Fmax of latch)
xname = 'Fmax';
xrange = [-1 3];
Fmaxs=logspace(xrange(1),xrange(2),N);

%setting y axis value on plot (Vmax of latch)
yname = 'vmax';
yrange = [-2 2];
v_maxs=logspace(yrange(1),yrange(2),N);

metrics = {'tto','vto','Pmax','ymax','tL','KEmax','yunlatch'};

load_time_constraint=Inf;

%parameters for the loading motor
Fmax_motor = 20;
range_of_motion = 3;
vmax_motor=10.0000;
%extra parameters for hill muscle motor
muscle_length=.1;
r_activation=1E-5;
%struct initialization
loading_motor = linear_motor(Fmax_motor, vmax_motor, range_of_motion);

%parameters for the load and struct initialization
m=100;
load = load_mass(m);

%parameters for the latch and struct initialization
R=2E-1;
m_L= 1E5;
coeff_fric = .1;
v_0L=0;
latch = rounded_latch(R, m_L, coeff_fric, v_0L);

%parameters for the spring and struct initialization
k=.3;
m_s=1E-4;
F_spring_max=1E4;
%characteristic_length for exponential spring
characteristic_length = 5;
spring = linear_spring(k, m_s, F_spring_max);
%spring=exponential_spring(k, characteristic_length, m_s,F_spring_max);





% initialize an output value matrix for each metric
for ii=1:length(metrics)
    outval{ii}=zeros(N);
end
if (debug)
    h1 = figure()
    h2 = figure()
end
for i=1:N %iterate over y-axis-variable of plot
    for j=1:N %iterate over x-axis-variable of plot
        %unlatching_motor = hill_muscle_motor(muscle_length, Fmaxs(j), v_maxs(i),r_activation);
        unlatching_motor = linear_motor(Fmaxs(j),v_maxs(i), range_of_motion);
        %input structs for each component
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