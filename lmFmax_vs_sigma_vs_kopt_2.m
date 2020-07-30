%% Fmax vs. sigma vs k_optimal/highest power 2

% this script does the same thing as Fmax vs. sigma vs k_optimal/highest
% power, except instead of finding a k_optimal for each Fmax vs. sigma, we
% find the k_optimal for Fmax vs. sigma = 0, and set that same k_optimal
% for sigma at that Fmax.

%% initializing variables

clearvars
addpath(genpath(fullfile(pwd,'..')));

% determines resolution of heatplots
N=10; 

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
% latch = rounded_latch(0.005,0.003,0,0,0,Inf,0);
latch = rounded_latch(0.005,0.003,0,0,0.003,0.0035,0);

%% calculating

% initializing output matrices
highest_P = zeros(N);
best_k = zeros(1,N);

for jj = 1:length(k_values)
    % looping through k values to find the best one
    spring = linear_spring(k_values(jj),0,Inf);
    
    Pmax_values_nonstochastic = zeros(size(Fmax_values));
    
    for j = 1:N %x var
        loading_motor = linear_motor(Fmax_values(j),10,0.005,1);

        %calling solve_model
        [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

        met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
        Pmax_values_nonstochastic(j)=met_dict('Pmax');
    end
    
    Pmax_values_sigma = zeros(N);
    Pmax_values_sigma(1,:) = Pmax_values_nonstochastic;
    for i = 2:N 
        Pmax_sigma = zeros(size(Pmax_values_nonstochastic));
        for j = 1:N       
            mu = Fmax_values(j);
            sigma = sigma_values(i)*mu;
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

    for x = 1:N
        if (Pmax_values_sigma(1,x) > highest_P(1,x))
            highest_P(1,x) = Pmax_values_sigma(1,x);
            best_k(x) = k_values(jj);
        end
    end
    
    disp(['row ' num2str(jj) ' of ' num2str(N)]);
end

%% calculating new max power

% recalculating max power based on the doctored k_optimal values

highest_P_sigma0 = zeros(size(Fmax_values));
for j = 1:N 
    spring = linear_spring(best_k(j),0,Inf);
    loading_motor = linear_motor(Fmax_values(j),10,0.005,1);

    %calling solve_model
    [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

    met_dict=get_metrics(sol,transition_times,load,{'Pmax'});
    highest_P_sigma0(j)=met_dict('Pmax');
end

new_highest_P = zeros(N);
new_highest_P(1,:) = highest_P_sigma0;
for i = 2:N %y var
    Pmax_sigma = zeros(size(highest_P_sigma0));
    for j = 1:N %x var        
        mu = Fmax_values(j);
        sigma = sigma_values(i)*mu;
        if (sigma < (Fmax_values(2)-Fmax_values(1)))
            Pmax_sigma(j) = highest_P_sigma0(j);
        else
            gaus = @(k)(1/(sigma*sqrt(2*pi)))*exp(-(((k-mu).^2)/(2*sigma.^2)));
            Pmax_sigma(j) = trapz(Fmax_values,gaus(Fmax_values).*highest_P_sigma0);
        end
        [maxval, index] = max(Pmax_sigma);
        new_highest_P(i,j) = maxval;
    end
end
    
%% plotting best_k

hmap = figure;

% plotting heatmap
imagesc(xrange,yrange,repmat(best_k,[N,1]));

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
imagesc(xrange,yrange,new_highest_P);

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