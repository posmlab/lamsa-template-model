%% Fmax vs. sigma vs k_optimal/highest power 2

% this script does the same thing as Fmax vs. sigma vs k_optimal/highest
% power, except instead of finding a k_optimal for each Fmax vs. sigma, we
% find the k_optimal for Fmax vs. sigma = 0, and set that same k_optimal
% for sigma at that Fmax.

%% initializing variables

clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines resolution of heatplots
N=20; 

% setting x axis on the plot 
xname = 'loading motor Fmax';
xrange = [0 40];
Fmax_values = linspace(xrange(1),xrange(2),N);

%setting y axis value on plot 
yname = 'standard deviation';
yrange = [0 1];
sigma_values = linspace(yrange(1),yrange(2),N);

% range of k values to loop through to find the best one
krange = [0 1E4];
k_values = linspace(krange(1),krange(2),N);

% load mass struct initialization
load = load_mass(0.01,0,1);
unlatching_motor = linear_motor(0.25,1,0.005,1);

% play with the min and max latching distance!
latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);
% latch = rounded_latch(0.005,0.003,0,0,0.003,0.0035,0);

%% calculating

% initializing output matrices
Pmax_values_nonstochastic = zeros(N);

for jj = 1:length(Fmax_values)
    % looping through k values to find the best one
    loading_motor = linear_motor(Fmax_values(jj),10,0.005,1);
    
    for j = 1:length(k_values) %x var
        spring = linear_spring(k_values(j),0,Inf);

        %calling solve_model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
        Pmax_values_nonstochastic(j,jj)=met_dict('Pmax');
    end
    disp(['row ' num2str(jj) ' of ' num2str(N)]);
end
    
[highest_P, best_k_index] = max(Pmax_values_nonstochastic);

%% calculating new max power

% recalculating max power based on the doctored k_optimal values

Power = zeros(N);
for i = 1:length(best_k_index)
    spring = linear_spring(k_values(best_k_index(i)),0,Inf);
    PvsF1D = zeros(size(Fmax_values));
    for jj = 1: length(Fmax_values)
        loading_motor = linear_motor(Fmax_values(jj),10,0.005,1);
        
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
        PvsF1D(jj)=met_dict('Pmax');
    end
    for j = 1:length(sigma_values)
        mu = Fmax_values(i);
        sigma = sigma_values(j)*mu;
        if (sigma < (Fmax_values(2)-Fmax_values(1)))
            Power(j,i) = PvsF1D(i);
        else
            gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
            Power(j,i) = trapz(Fmax_values,gaus(Fmax_values).*PvsF1D);
        end
    end
end
    
%% plotting best_k

hmap = figure;

% plotting heatmap
imagesc(xrange,yrange,repmat(k_values(best_k_index),[N,1]));

% setting axes labels
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','latex')
xlabel(xname,'Interpreter', 'Latex');
ylabel(yname, 'Interpreter', 'Latex');

% colorbar settings
c = colorbar;
c.Label.String = 'k';
set(c,'TickLabelInterpreter','latex')
c.Label.Interpreter="latex";

%% plotting highest_p

hmap2 = figure;

% plotting heatmap
imagesc(xrange,yrange,Power);

% setting axes labels
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','latex')
xlabel(xname,'Interpreter', 'Latex');
ylabel(yname, 'Interpreter', 'Latex');

% colorbar settings
c = colorbar;
c.Label.String = 'highest P with best k';
set(c,'TickLabelInterpreter','latex')
c.Label.Interpreter="latex";