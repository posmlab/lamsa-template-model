
function [sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
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
    F_list(i) = loading_motor.Force(Inf,[y_list(i) 0]);
end
y_guess_motor = -y_list(find(F_list>0,1,'last'));

% initial guess based on initial spring stiffness
y_guess_spring = loading_motor.Force(Inf,[0 0])/((spring.Force(0,-100*eps)-spring.Force(0,0))/(100*eps));

% use fzero to find when Fs=Fin does this work for exponential spring?
y_guess = max([y_guess_motor, y_guess_spring]);
options =  {};% optimset('Display','iter');
[y0,~,exitflag]=fzero(@(y) (loading_motor.Force(Inf,[y 0])-spring.Force(0,[y 0])) - LARGE_NUM*((~loading_motor.Force(Inf,[y 0]))||(~spring.Force(0,[y 0])))+LARGE_NUM*(y>0),y_guess,options);
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
inst_check = 1;
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
    
    tspan=linspace(0,t_L_guess*1000,1000000);%[0,t_L_guess]; Note 2nd and 3rd inputs are multiplied by 1000 for temporary fix
    [t_unlatch,x_unlatch]=ode45(ode,tspan,[0 latch.v_0],unlatch_opts);
    % This ODE is for the latch x-coordinate, but we want the y-coordinate, so
    % convert
    y_unlatch=zeros(size(x_unlatch));
    X = [x_unlatch(:,1), x_unlatch(:,2)];
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


%% Solving for Normal Force in the unlatching phase
for i=1:size(X, 1)% For derivation of this equation for F_n see Overleaf doc with LaMSA derivation
    num1 = (latch.mass*spring.Force(t_unlatch(i), y_unlatch(i, :))) - ...
        (m_eff*latch.y_L{3}(X(i,1))*(X(i,2)^2)*latch.mass) - ...
        (unlatching_motor.Force(t_unlatch(i), X(i,:))*m_eff*latch.y_L{2}(X(i,1)));
    rad = 1 + ((latch.y_L{2}(X(i,1)))^2);
    num2 = sqrt(rad);
    den1 = m_eff*latch.y_L{2}(X(i,1))*(latch.y_L{2}(X(i,1)) - latch.coeff_fric);
    den2 = latch.mass*(1+latch.coeff_fric*latch.y_L{2}(X(i,1)));
    F_n(i) =(num1*num2)/(den1 + den2);%filling in the F_n vector until unlatch time
end
F_n = F_n';%switching to a column vector so we can add it to sol
%disp(F_n)

%% Components of Normal Force And Frictional Force
% Currently not working, some trig or possibly t_L issues

% Defining the geometric definitions of sine and cosine

for i=1:size(X, 1)
    den = sqrt(1 + (latch.y_L{2}(X(i, 1))^2));
    sin_comp(i) = (latch.y_L{2}(X(i, 1)))/den;
    cos_comp(i) = 1/den;
end

sin_comp = sin_comp';
cos_comp = cos_comp';

%Calculating Normal and Frictional Force Components
for i=1:size(X, 1);
    F_nx(i) = F_n(i) .* sin_comp(i);
    F_ny(i) = F_n(i) .* cos_comp(i);
    F_fx(i) = F_nx(i) * latch.coeff_fric;
    F_fy(i) = F_ny(i) * latch.coeff_fric;
end
F_nx = F_nx';
F_ny = F_ny';
F_fx = F_fx';
F_fy = F_fy';

F_comp = [F_nx F_ny F_fx F_fy];


%% Ballistic phase:Fs only
%guess launch times by treating the spring as ideal-ish and getting the
%   frequency
stiffness = abs(( spring.Force(0,y_unlatch(end,:))-spring.Force(0,y_unlatch(end,:)+(100*eps)) ) / (100*eps)); %Here be divide by 0 errors, probably
%stiffness=abs(spring.Force(0,y_unlatch(end,:))/y_unlatch(end,1)); %Here be divide by 0 errors, probably
nat_freq=sqrt(stiffness/m_eff);
t_launch_guess=pi/nat_freq;
launch_opts=odeset('Events',@(t,y) launching_end(t,y,spring));
ode=@(t,y) launching_ode(t,y,m_eff,spring);
tspan=linspace(0,t_launch_guess,1E3);
y0=y_unlatch(end,:)';
[t_launch,y_launch]=ode45(ode,tspan,y0,launch_opts);

% Solve latch dynamics during Ballistic Phase
%     Currently assuming instaneous stopping of latch

%% Stitch together solutions
transition_times=[t_unlatch(end),t_unlatch(end)+t_launch(end)];
T=[t_unlatch;t_unlatch(end)+t_launch];
Y=[y_unlatch;y_launch];
fill = repmat([X(end, 1), 0],length(Y) - length(X), 1);%Makes X the right size to fit in sol
%fill2 = repmat([0],length(Y) - length(X), 1);%Makes F_n the right size to fit in sol
%F_N = [F_n;fill2];
xFinal = [X;fill];
%F_compFinal = [F_comp; fill];
for i = 1:size(T)
    fSpring(i) = spring.Force(T(i), [Y(i,1), Y(i,2)]);%fill out the fSpring vector to add to sol
    fUnlatchingMotor(i) = unlatching_motor.Force(T(i), [xFinal(i,1), xFinal(i,2)]);
end
fSpring = fSpring';
fUnlatchingMotor = fUnlatchingMotor';
% add zeros to the end of F_comp because F_comp consists of 
% forces during the unlatching phase, and the other vectors
% include forces during the ballistic phase. 
% Adding zeros makes this matrix the right size for appending 
% to the rest of the sol.
F_comp = [F_comp; zeros(size(T,1)-size(F_comp,1),4)];

% stitch together various numbers 
% for one big matrix to write to csv file 
sol=[T Y xFinal F_comp fSpring fUnlatchingMotor];


%% Establishing Parameters for .json output
% params = struct();
% 
% params.m_L = latch.mass;
% 
% params.m_eff = m_eff;
% 
% params.v0 = latch.v_0;
% 
% func_struct = functions(loading_motor.Force);
% params.loading_motor = rmfield(func_struct,{'file','type','within_file_path'});
% 
% func_struct = functions(spring.Force);
% params.spring = rmfield(func_struct,{'file','type','within_file_path'});
% 
% func_struct = functions(unlatching_motor.Force);
% params.unlatching_force = rmfield(func_struct,{'file','type','within_file_path'});
% 
% func_struct = functions(latch.y_L{1});
% params.latch_shape = rmfield(func_struct,{'file','type','within_file_path'});
% 
% params.unlatch_time = transition_times(1);
% 
% params.launch_time = transition_times(2);
% 
% %% Making output directory
% if ~isdir(outputDirectory)%checks for and possibly creates output directory
%     mkdir(outputDirectory)
% end
% 
% timeStampString = string(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd_HH-mm-ss-SSS'));
% 
% %% Writing .json output
% pretty = prettyjson(['parameters: ' jsonencode(params)]);%makes the parameter json file readable
% 
% fileName = outputDirectory + sprintf("/parameters_%s.json", timeStampString);%ensures the file is in the output directory
% 
% fileID = fopen(fileName, 'w');
% 
% fwrite(fileID, pretty, 'char');
% fclose(fileID);
% 
% %% save solution data to csv and json files
% 
% 
% 
% %replace spaces with underscores
% dateString = string(datetime);
% cleanDateString = regexprep(dateString, " ", "_");
% cleanDateString = regexprep(cleanDateString, ":", "_");%creates file friendly output name
% 
% 
% csvFilePath = outputDirectory + "/raw_data--" + timeStampString + ".csv";%same as above for json file, ensures location
% 
% headers = ["Time", "y", "ydot", "x", "xdot", "normal force on latch x", ...
%     "normal force on load y", "frictional force on latch x", ...
%     "frictional force on load y", "spring force", ...
%     "unlatching motor force"];
% 
% writematrix(headers, csvFilePath);%creates headers for output file
% writematrix(sol, csvFilePath, 'WriteMode', 'append');% addes actual data to csv file

end

