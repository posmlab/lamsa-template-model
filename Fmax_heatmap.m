%% Fmax vs. sigma vs k_optimal/highest power

% this script plots k_optimal and Pmax as a function of Fmax and sigma

% Plot 1
% x-axis: Fmax
% y-axis: sigma
% z-axis: k_optimal

% Plot 2
% x-axis: Fmax
% y-axis: sigma
% z-axis: Pmax when using k_opt

% Make sure to ignore the data after x=15 because the 1D plot starts
% curving down there, so it is misleading (see lmFmax vs. Pmax for more 
% info). Basically, if your x range is 0-40, then only consider 0-20.

%% initializing variables

clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines resolution of heatmap
N=10; 

% setting x axis on the plot 
xname = 'loading motor Fmax';
xrange = [0 40];
Fmax_values = linspace(xrange(1),xrange(2),N);

% setting y axis value on plot 
yname = 'standard deviation';
yrange = [0 1];
simga_values = linspace(yrange(1),yrange(2),N);

% range of k values to loop through to find the best one
krange = [0 1E4];
k_values = linspace(krange(1),krange(2),N);

% component struct initialization
load = load_mass(0.01,0,1);
unlatching_motor = linear_motor(0.25,1,0.005,1);

% play with the min and max latching distance!
latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);
% latch = rounded_latch(0.005,0.003,0,0,0.003,0.0035,0);

%% calculations 

% initializing output matrices
highest_P = zeros(N);
best_k = zeros(N);

for jj = 1:length(k_values)
    % looping through k values to find the best one for each Fmax sigma
    % combination
    spring = linear_spring(k_values(jj),0,Inf);
    
    % initializing Pmax array for when sigma = 0
    Pmax_values_nonstochastic = zeros(size(Fmax_values));
    for j = 1:N 
        % looping through Fmax values
        loading_motor = linear_motor(Fmax_values(j),10,0.005,1);

        %calling solve_model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
        Pmax_values_nonstochastic(j)=met_dict('Pmax');
    end

    % initializing Pmax matrix for a non-zero sigma
    Pmax_values_sigma=zeros(N);

    % setting the first row of the Pmax values to the sigma = 0 values we already calcualted
    Pmax_values_sigma(1,:) = Pmax_values_nonstochastic;
    
    for i = 2:N % looping through sigmas
        Pmax_sigma = zeros(size(Pmax_values_nonstochastic));
        for j = 1:N % looping through Fmax values      
            mu = Fmax_values(j);
            sigma = simga_values(i)*mu;
            if (sigma < (Fmax_values(2)-Fmax_values(1)))
                Pmax_sigma(j) = Pmax_values_nonstochastic(j);
            else
                gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
                Pmax_sigma(j) = trapz(Fmax_values,gaus(Fmax_values).*Pmax_values_nonstochastic);
            end
            [maxval, index] = max(Pmax_sigma);
            Pmax_values_sigma(i,j) = maxval;
        end
    end

    % checking if this k value is better than previous k values
    for x = 1:N
        for y = 1:N
            if (Pmax_values_sigma(x,y) > highest_P(x,y))
                % saving the Pmax gotten from best_k
                highest_P(x,y) = Pmax_values_sigma(x,y);
                % saves the best k value
                best_k(x,y) = k_values(jj);
            end
        end
    end
    disp(['row ' num2str(jj) ' of ' num2str(N)]);
end

%% plotting best_k

hmap = figure;

% plotting heatmap
imagesc(xrange,yrange,best_k);

% setting axes
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','latex')
xlabel(xname,'Interpreter', 'Latex');
ylabel(yname, 'Interpreter', 'Latex');

% colorbar settings
c = colorbar;
c.Label.String = 'k';
set(c,'TickLabelInterpreter','latex')
c.Label.Interpreter = 'latex';

%% plotting highest_P

hmap2 = figure;

% plotting heatmap
imagesc(xrange,yrange,highest_P);

% setting axes
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','latex')
xlabel(xname,'Interpreter', 'Latex');
ylabel(yname, 'Interpreter', 'Latex');

% colorbar settings
c = colorbar;
c.Label.String = 'highest P with best k';
set(c,'TickLabelInterpreter','latex')
c.Label.Interpreter="latex";
