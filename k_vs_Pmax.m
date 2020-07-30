%% k vs. Pmax 1D plot

% this script plots maximum power as a function of spring constant

% x-axis: k
% y-axis: Pmax

% two sets of points are plotted:
% 1. k vs. Pmax for a regular spring
% 2. k vs. Pmax for a stochastic spring

% With the stochastic spring, we are considering: what happens to 
% Pmax if we use a spring where the spring constant is slightly different
% every time the spring is used? To model this, we generate a
% gaussian where the each k value is the mean and the standard deviation
% is some fraction of that mean, and integrate over the gaussian.

% Make sure to ignore the stochastic spring data as soon as it starts
% curving down (this will only happen if Pmax doesnt go to almost 0), 
% because the downward curve is due to a lack of data to the 
% right of the graph, so it is misleading. If Pmax is going almost to
% 0 though, it is not a problem.

%% initializing variables

clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines number of points to plot
N=500; 

% determines range of x values
xrange = [1 5E4];
k_values = linspace(xrange(1),xrange(2),N);

% initializing component structs
load = load_mass(0.01,0,1);
loading_motor = linear_motor(10,10,0.005,1);
unlatching_motor = linear_motor(0.25,1,0.005,1);

% play with the min and max latching distance!
latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);
% latch = rounded_latch(0.005,0.003,0,0,0.002,0.003,0);

%% calculating Pmax as a function of k

% initializing Pmax output array
Pmax_values = zeros(size(k_values));

for i=1:N
    % looping through different spring constants
    spring = linear_spring(k_values(i),0,Inf);

    % calling solve model
    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

    % populating output array
    met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
    Pmax_values(i)=met_dict('Pmax');

    disp(['point ' num2str(i) ' of ' num2str(N)]);
end

%% plotting k vs. Pmax for non-stochastic spring (sigma = 0)

fh = figure;

% plot non-stochastic spring points
plot(k_values,Pmax_values,'.');

% setting axis labels
xlabel('k','Interpreter', 'Latex');
ylabel('Pmax', 'Interpreter', 'Latex');

% sets the y-axis to something reasonable
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

%% calculating k vs. Pmax for a stochastic spring 

% initializing output matrix
Pmax_sigma = zeros(size(Pmax_values));

% looping through all k values
for i = 1:length(k_values)
    % sets each k value as the mean of a gaussian
    mu = k_values(i);
    
    % sets the standard deviation as a percent of the mean
    sigma = 0.1*mu;
    
    % function for the gaussian
    gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
    
    % integrate over the gaussian
    Pmax_sigma(i) = trapz(k_values,gaus(k_values).*Pmax_values);
end

%% plotting k vs. Pmax for stochastic spring (sigma ~=0)

plot(k_values,Pmax_sigma,'ok')
