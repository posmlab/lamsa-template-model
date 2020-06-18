
function [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring)
%Solve set of differential equations for loading, unlatching, and launching
%   phases of LAMSA motion
LARGE_NUM = 1E10; % subtract a large number in fzero function to trick fzero into identifying points where motor or spring suddnely go to 0
m_eff = load.mass + spring.mass/3;


%% Loading phase: Fs vs. Fin
% Assume quasistatic loading, end position is when Fs=Fin

% finding an order-of-magnitude initial guess for fzero call from motor properties
y_list = logspace(-20,20,40); % sweep through 40 orders of magnitude
F_list = zeros(size(y_list));
for i = 1:length(y_list)
    F_list(i) = loading_motor.Force(0,[y_list(i) 0]);
end
y_guess_motor = -y_list(find(F_list>0,1,'last'));

% initial guess based on initial spring stiffness
y_guess_spring = loading_motor.Force(0,[0 0])/((spring.Force(0,10*eps)-spring.Force(0,0))/(10*eps));

% use fzero to find when Fs=Fin does this work for exponential spring?
y_guess = max([y_guess_motor, y_guess_spring]);
options =  {};% optimset('Display','iter');
[y0,~,exitflag]=fzero(@(y) (loading_motor.Force(0,[y 0])-spring.Force(0,[y 0])) - LARGE_NUM*((~loading_motor.Force(0,[y 0]))||(~spring.Force(0,[y 0])))+LARGE_NUM*(y>0),y_guess,options);
if (exitflag<0)
    error('fzero failed');
end

%% Unlatching phase: Fs vs Flatch


F_friction = latch.coeff_fric*spring.Force(0,[y0,0]);
if unlatching_motor.max_force <= F_friction
    warning('latch cannot overcome friction force');
    sol = [0,0,0];
    transition_times = [0,0];
    return
end 


[inst_check,~,~]=unlatching_end(0,[0,latch.v_0],m_eff,y0,latch,spring,unlatching_motor);
if inst_check>0 
    unlatch_opts=odeset('Events',@(t,y) unlatching_end(t,y,m_eff,y0,latch,spring,unlatching_motor),'RelTol',1E-7,'AbsTol',1E-10);
    ode=@(t,y) unlatching_ode(t,y,m_eff,y0,latch,spring,unlatching_motor);
    
    a_0L = unlatching_motor.max_force / latch.mass;
    
    if (a_0L ~= 0)
        % calculate t_L_guess using quadratic formula
        % and the following kinematic equation: R = (1/2)a*t^2 + v_0*t  
        t_L_guess = (((-1*latch.v_0) + sqrt((latch.v_0)^2  + (2*a_0L*latch.max_width)))/(a_0L));
    elseif (latch.v_0 ~= 0 )
        t_L_guess = latch.max_width/latch.v_0;        
    else
        warning("The latch's initial velocity and acceleration are both zero.")
        sol = [0,0,0]
        transition_times = [0,0]
        return
    end
    
    tspan=linspace(0,t_L_guess,1000);%[0,t_L_guess];
    [t_unlatch,x_unlatch]=ode45(ode,tspan,[0 latch.v_0],unlatch_opts);
    % This ODE is for the latch x-coordinate, but we want the y-coordinate, so
    % convert
    y_unlatch=zeros(size(x_unlatch));
    for i=1:length(y_unlatch)
        y_unlatch(i,1)=latch.y_L{1}(x_unlatch(i,1))+y0;
        y_unlatch(i,2)=x_unlatch(i,2)*latch.y_L{2}(x_unlatch(i,1));
    end
    if (imag(y_unlatch))
        disp(y_unlatch);
        error('y unlatch imaginary');
    end
else % instantaneous unlatching
    y_unlatch=[y0,0]; %May cause a repeated time step and give NaNs on differentiation
    t_unlatch=[0];
end
t_unlatch = real(t_unlatch);
y_unlatch = real(y_unlatch);


%% Ballistic phase:Fs only
%guess launch times by treating the spring as ideal-ish and getting the
%   frequency

stiffness = abs( (spring.Force(0,y_unlatch(end,:)) - spring.Force(0,y_unlatch(end,:)+(10*eps))) / (10*eps)); 
%stiffness=abs(spring.Force(0,y_unlatch(end,:))/y_unlatch(end,1)); %Here be divide by 0 errors, probably
nat_freq=sqrt(stiffness/m_eff);
t_launch_guess=pi/nat_freq;
launch_opts=odeset('Events',@(t,y) launching_end(t,y,spring));
ode=@(t,y) launching_ode(t,y,m_eff,spring);
tspan=linspace(0,t_launch_guess,1E3);
y0=y_unlatch(end,:)';
[t_launch,y_launch]=ode45(ode,tspan,y0,launch_opts);

%% Stitch together solutions
transition_times=[t_unlatch(end),t_unlatch(end)+t_launch(end)];
T=[t_unlatch;t_unlatch(end)+t_launch];
Y=[y_unlatch;y_launch];
sol=[T Y];


%% Establishing Parameters for .json output
% params = struct();
% 
% params.m_L = m;
% 
% params.m_eff = m_eff;
% 
% params.v0 = Latch.v_0;
% 
% func_struct = functions(motor_in);
% params.loading_motor = rmfield(func_struct,{'file','type','within_file_path'});
% 
% func_struct = functions(spring);
% params.spring = rmfield(func_struct,{'file','type','within_file_path'});
% 
% func_struct = functions(yL);
% params.latch_shape = rmfield(func_struct,{'file','type','within_file_path'});
% 
% %% Writing .json output
% pretty = prettyjson(['parameters: ' jsonencode(params)]);
% timestamp = datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd_HH-mm-ss');
% fileName = sprintf("parameters_%s.json", timestamp);
% 
% fileID = fopen(fileName, 'w');
% if fileID == -1, error('Cannot create JSON file'); end
% fwrite(fileID, pretty, 'char');
% fclose(fileID);
% end
% 
% %% save solution data to csv and json files
% 
% outputDirectory = "output";
% if ~isdir(outputDirectory)
%     mkdir(outputDirectory)
% end
% 
% %replace spaces with underscores
% dateString = string(datetime);
% cleanDateString = regexprep(dateString, " ", "_");
% cleanDateString = regexprep(cleanDateString, ":", "_");
% 
% 
% csvFilePath = outputDirectory + "/output--" + cleanDateString + ".csv";
% 
% headers = ["Time", "y", "ydot", "x", "xdot", "normal force on latch x", ...
%     "normal force on load y", "frictional force on latch x", ...
%     "frictional force on load y", "spring force", ...
%     "unlatching motor force into"];
% 
% writematrix(headers, csvFilePath);
% writematrix(sol, csvFilePath, 'WriteMode', 'append');
% 
% end
% 
