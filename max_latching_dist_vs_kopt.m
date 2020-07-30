%% max latching distance vs. k_optimal

% this script plots k_optimal as a function of max latching distance with a 
% stochastic spring

% x-axis: max latching distance
% y-axis: k_optimal

% How this works:
% 1. For each max latching distance point that is sampled, a k vs. Pmax 1D
% plot is generated (see k_vs_Pmax_1Dplot). 
% 2. The k value that produces the highest power is called "k_optimal"

% Change the sigma percentage on line 82 to see how the degree of variation
% in the stochastic spring to see how it affects the optimal k value
% spoiler alert: a stiffer spring is more optimal for a larger standard
% deviation

%% initializing variables

clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines resolution of 1D plots used to find k_optimal
N = 500; 

% determines range of x values used in 1D plot
xrange = [1 100000];
k_values = linspace(xrange(1),xrange(2),N);

% initializing component structs
load = load_mass(0.01,0,1);
loading_motor = linear_motor(10,10,0.005,1);
unlatching_motor = linear_motor(0.25,1,0.005,1);

% determines number of points sampled in max latching distance range
M = 20;

% determines range of max latching distances 
mld_range = [0.001 0.006];
mld_range = linspace(mld_range(1), mld_range(2), M);

% initializing max latching distance vs. k_optimal output array
kopt_values = zeros(1,M);

%% calculating

for mld = 1:M
    % looping through max latching distances
    latch = rounded_latch(0.005,0.003,0,0,0,mld_range(mld),0);
    
    % initializing k vs. Pmax output array
    Pmax_values = zeros(size(k_values));

    % looping over k values to find k_optimal at that max latching distance
    % (see k_vs_Pmax_1Dplot)
    for i=1:N    
        % looping over k values
        spring = linear_spring(k_values(i),0,Inf);

        % calling solve model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        % populating output array
        met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
        Pmax_values(i)=met_dict('Pmax');
    end
    % making the spring stochastic
    Pmax = Pmax_values;
    k = k_values;
    Pmax_sigma = zeros(size(Pmax));
    for i = 1:length(k)
        % mean of the gaussian
        mu = k_values(i);
        % standard deviation as a percent of the mean
        % play around with this!
        sigma = 0.1*mu;
        % gaussian function
        gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
        % integrate over the gaussian
        Pmax_sigma(i) = trapz(k,gaus(k).*Pmax);
    end
    
    % finds the max power in the 1D plot
    [maxval, index] = max(Pmax_sigma);
    
    % saves the k value that produced that Pmax into the output array
    kopt_values(mld) = k_values(index);
    
    disp(['point ' num2str(mld) ' of ' num2str(M)]);
end

%% plotting

figure

plot(mld_range,kopt_values,'ok');

xlabel('max latching distance','Interpreter','latex');
ylabel('kopt','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')