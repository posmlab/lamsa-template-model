%% LoebMuscleMotor class definition
% Mostly taken from Tsianos and Loeb 2013
% Looks at Af and general system as seen in Brown, Cheng, and Loeb 1999 (could be
% modified for slow-twitch)
% 
%arguments in required order
%     L_0 - rest length of the muscle (in m)
%     Fmax - muscle maximum tetanic force (in N), 
%     Vmax - maximum velocity at which the fascicle can lengthen in m/s
%     range_of_motion - range of motors motion
%     L_r1 - rest length of passive element 1 in L_0
%     L_r2 - rest length of passive element 2 in L_0
%     R - recruitment factor, active or recruited fraction of muscle 
%           cross-sectional area, 0<=R<=1, unitless (optional)
%     actually not sure about f units
%     f - stimulus frequency in pps(that's the same as Hertz)(optional)
% min # arguments = 8

classdef LoebMuscleMotor < Motor
    
    methods (Static)
        % first row contains parameter names
        % second row contains default values for the loading motor
        % third row contains default values for the unlatching motor
        function parameters = parameters()
            parameters = ["L_0" "Fmax" "Vmax" "range_of_motion" "L_r1" "L_r2"...
                          "R" "f";
                "2.45e-6" "10" "7.39" "10" "3.185e-6" "1.936e-6" "0.5" "120";
                "2.45e-6" "10" "7.39" "10"  "3.185e-6" "1.936e-6" "0.5" "120";
                "0" "0" "0" "0" " 0" "0" "0" "0";
                "Inf" "Inf" "Inf" "Inf" "Inf" "1" "Inf" "Inf"];
        end
    end
    
    methods
        function obj = LoebMuscleMotor(L_0,Fmax,Vmax,range_of_motion, L_r1, L_r2, varargin)
            % optional parameters
            varargin_param_names = {'R', 'f'};
            varargin_default_values = {0.5, 120};

            % check and assign optional parameters
            if (nargin < 6)
                error('Loeb muscle motor requires at least 6 arguments.');
            end
            if (length(varargin)>length(varargin_param_names))
                error('Too many input parameters');
            end
            for i=1:length(varargin)
                eval([varargin_param_names{i} '=varargin{i};'])
            end
            for i=(length(varargin)+1):length(varargin_param_names)
                eval([varargin_param_names{i} '=varargin_default_values{i};'])
            end
            if (Vmax == 0)
                Vmax = Inf;
                warning("Vmax argument must be nonzero, setting to Inf")
            end
            % model
            %model from Brown and Loeb 1999 and Tsianos and Loeb 2013
            %nomralize inputs
            
            %constants for fast-twitch=> change these to inputs mostly
            a_f=0.56;
            n_f0=2.1;
            n_f1=3.3;
            %for fast-twitch fibers
            c_1=23;
            k_1=0.046;
            
            %
            c_2=-0.02;
            k_2=-21;
            
            omega=0.75;
            B=1.55;
            rho=2.12;
            
            c_v0=-5.7;
            c_v1=9.18;
            a_v0=-1.53;
            a_v1=0;
            a_v2=0;
            b_v=0.69;
            n=0.01;
        
            %changing variables:

            L=@(t,x) x(1)+L_0;
            V=@(t,x) x(2);
            
            %Big equations

            n_f= @(t,x)n_f0+n_f1*(1/(L(t,x))-1);
            % use from 1999 paper because I don't know input frequencies, unitless, gives activation frequency, maybe just ask for it directly, 0<+Af<=1 maybe use *((1-exp(-(f/(a_f*n_f(t,x)))^n_f(t,x)))>=1)*((1-exp(-(f/(a_f*n_f(t,x)))^n_f(t,x)))<=0)
            Af=@(t,x) 1-exp(-(f/(a_f*n_f(t,x)))^n_f(t,x));
           
            F_PE1=@(t,x) c_1*k_1*log(exp((L(t,x)-L_r1)/k_1)+1)+n*V(t,x);
            %F_PE2<=0 because it only acts contrary to motion(and then only
            %some of the time) only acts at lenghs shorter than 0.75L_0
            F_PE2=@(t,x)c_2*(exp(k_2*(L(t,x)-L_r2))-1);
            %0<=FL<=1, unitless
            FL=@(t,x)exp(-abs((L(t,x)^B-1)/omega)^rho);
            %0<=FV<=1, unitless
            FV=@(t,x)(V(t,x)<=0)*((Vmax-V(t,x))/(Vmax+(c_v0+c_v1*L(t,x))*V(t,x)))+(V(t,x)>0)*(b_v-(a_v0+a_v1*L(t,x)+a_v2*L(t,x)^2)*V(t,x))/(b_v+V(t,x));
            
            F_PE=@(t,x)F_PE1(t,x)+R*Af(t,x)*F_PE2(t,x);
            F_CE=@(t,x)R*Af(t,x)*FL(t,x)*FV(t,x);
    %final properties
            range=range_of_motion;
            velocity=@(t,x)V(t,x);
            max_force = Fmax;
            Force=@(t,x)F_CE(t,x)*Fmax + F_PE(t,x)*Fmax;
            % call parent constructor
            obj = obj@Motor(max_force, range, velocity, Force);
        end  
    end
end
%% Citations
% Rosario MV, Sutton GP, Patek SN, Sawicki GS. 2016 Muscle�spring dynamics in time-limited, elastic movements.
%   Proc. R. Soc. B 283: 20161561. http://dx.doi.org/10.1098/rspb.2016.1561

