%% loading motor Fmax vs. Pmax 1D plot

% this script plots maximum power as a function of the loading motor's Fmax

% x-axis: Fmax
% y-axis: Pmax

% two sets of points are plotted:
% 1. Fmax vs. Pmax for a regular motor
% 2. Fmax vs. Pmax for a stochastic stochastic motor

% With the stochastic motor, we are considering: what happens to 
% Pmax if we use a motor where the motor force is slightly different
% every time the motor is used? To model this, we generate a
% gaussian where the each Fmax value is the mean and the standard deviation
% is some fraction of that mean, and integrate over the gaussian.

% Make sure to ignore the stochastic motor data as soon as it starts
% curving down, because the downward curve is due to a lack of data to the 
% right of the graph, so it is misleading. Basically, if your x range is
% 0-40, then only consider 0-20.

%% initializing variables

clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines number of points to plot
N=100; 

% determines range of x values
xrange = [0 40];
Fmax_values = linspace(xrange(1),xrange(2),N);

load = load_mass(0.01,0,1);
spring = linear_spring(2000,0,Inf);
unlatching_motor = linear_motor(0.25,1,0.005,1);

% play with the min and max latching distances!
latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);
% latch = rounded_latch(0.005,0.003,0,0,0.003,0.0035,0);

%% calculating Pmax as a function of Fmax

% initializing output array
Pmax_values = zeros(size(Fmax_values));

for i=1:N
    % looping over Fmax values
    loading_motor = linear_motor(Fmax_values(i),10,0.005,1);

    % calling solve model
    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);
    
    % populating output array
    met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
    Pmax_values(i)=met_dict('Pmax');

    disp(['point ' num2str(i) ' of ' num2str(N)]);
end

%% plotting Fmax vs. Pmax for non-stochastic motor (sigma = 0)

% plot output
figure;

% plotting non-stochastic spring values
plot(Fmax_values,Pmax_values,'.'); 

% setting axis labels
xlabel('Fmax','Interpreter', 'Latex');
ylabel('Pmax', 'Interpreter', 'Latex');

% sets y-axis to something reasonable
ylim_range = max(Pmax_values)-min(Pmax_values);
if (ylim_range == 0)
    ylim_range = abs(max(Pmax_values));
    if ~(ylim_range == 0)
        ylim([min([0, min(Pmax_values)-0.05*ylim_range]) max([0,max(Pmax_values)+0.05*ylim_range]) ])
    end
else
    ylim([min(Pmax_values)-0.05*ylim_range max(Pmax_values)+0.05*ylim_range])
end
set(gca,'TickLabelInterpreter','latex')

hold on

%% calculating Fmax vs. Pmax for a stochastic motor

% initializing output matrix
Pmax_sigma = zeros(size(Pmax_values));

for i = 1:length(Fmax_values)
    % sets each Fmax values as the mean of the gaussian
    mu = Fmax_values(i);
    
    % sets the standard deviation as a percent of the mean
    sigma = 0.1*mu;
    
    % prevents numerical integration issue
    if (sigma < (Fmax_values(2)-Fmax_values(1)))
        Pmax_sigma(i) = Pmax_values(i);
    else
        %function for the gaussian
        gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
        % integrate over the gaussian
        Pmax_sigma(i) = trapz(Fmax_values,gaus(Fmax_values).*Pmax_values);
    end
end

%% plotting Fmax vs. Pmax for stochastic motor (sigma ~= 0)

plot(Fmax_values,Pmax_sigma,'ok')
