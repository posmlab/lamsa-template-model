%Vary F_max and m_eff and plot performance metrics
close all
clearvars
tic
debug = false;
N=50;



xname = 'Fmax';
xrange = [-1 3];
Fmaxs=logspace(xrange(1),xrange(2),N);

yname = 'vmax';
yrange = [-3 3];
v_maxs=logspace(yrange(1),yrange(2),N);

metrics = {'tto','vto','Pmax','ymax','tL','KEmax'};

L = 10E-3;
rho = 10;
sigma_f = 10E6;
Fmax = 20;
d = 5E-3;
v_max=5.0000;
motor_in=@(t,x) (Fmax*(1-x(2)/v_max)) .* (abs(x(1))<=d); %Linear F-v motor 

m=1E-3;
m_s=1E-4;
m_eff = m + m_s/3;
load_time_constraint=Inf;
F_spring_max=1E4;
k=1;
spring=@(t,x) -k*x(1).*(abs(k*x(1))<F_spring_max);


yL=@(x) Latch.max_width*(1-sqrt(1-x^2/Latch.max_width^2));
syms x;
yL_prime = diff(yL(x));
yL_doubleprime = diff(yL(x),2);
Latch.y_L = {yL, matlabFunction(yL_prime), matlabFunction(yL_doubleprime)};
Latch.mass=1E2;
Latch.v_0=0;
Latch.max_width=2E-1;

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
         motor_out=@(t,x) (Fmaxs(j)*(1-x(2)/v_maxs(i))) .* (abs(x(1))<=d); %Linear F-v motor 
        [sol,transition_times]=solve_model(motor_in,motor_out,spring,m_eff,Latch);
        %[sol,transition_times]=solve_model(motor_in,motor_out,spring,m_eff,m_L,R/v_0L,v_0L,y_L);
        if (debug)
            figure(h1)
            plot(sol(:,1),sol(:,2),'.');
            hold on;
            figure(h2)
            plot(sol(:,1),sol(:,3),'.');
            hold on;
            ginput(1)
        end
        met_dict=get_metrics(sol,transition_times,m,metrics);
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