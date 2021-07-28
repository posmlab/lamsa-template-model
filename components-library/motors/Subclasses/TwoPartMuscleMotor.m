%% TwoPartMuscleMotor class definition
%as seen in:
% Guenther, M, S Schmitt, and V Wank. 2007. 'High-Frequency Oscillations as
% a Consequence of Neglected Serial Damping in Hill-Type Muscle Models.'
% Biological Cybernetics 97 (1) (July): 63-79.
% doi:10.1007/s00422-007-0160-6.
%
% arguments in required order:

%     L_0 - rest length of the contractile element
%     v_motor_max - maximum velocity at which the motor can travel
%     F_max - F_max in [N] for Extensor
%     l_CEopt - optimal length of CE in [m] for Extensor
%     L_PEE0 - rest length of PEE normalized to optimal lenght of CE (Guenther et al., 2007)

%optional
%     L_initial - initial stretch on the contractile element  
%       min # arguments = 5

classdef TwoPartMuscleMotor < Motor
    
    methods (Static)
        % first row contains parameter names
        % second row contains default values for the loading motor
        % third row contains default values for the unlatching motor
        function parameters = parameters()
            parameters = [ "L_0" "v_motor_max" "F_max" "l_CEopt"  "L_PEE0"...
                "L_initial";
                 "0.0828" "1" "5" "0.092" "0.9" "0";
                 "0.0828" "1" "5" "0.092" "0.9" "0";
                 "0" "0" "0" "0" "0" "0";
                 "Inf" "Inf" "Inf" "Inf" "Inf" "Inf"];
        end
    end
    
    methods
        function obj = TwoPartMuscleMotor( L_0, v_motor_max, F_max, l_CEopt, L_PEE0, varargin)
            
            % optional parameters
            varargin_param_names = {'L_initial'};
            varargin_default_values = {0};
           % 
             
%             % check and assign optional parameters
            if (nargin < 4)
                error('Two Part muscle motor requires at le%ast 4 arguments.');
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
            
            if (v_motor_max == 0)
                v_motor_max = Inf;
                warning("v_max argument must be nonzero, setting to Inf")
            end
            % model
            % constants
            v_PEE=2.5;%    v_PEE -  exponent of F_PEE (Moerl et al., 2012)
            F_PEE=2.0;%    F_PEE - force of PEE if l_CE is stretched to deltaWlimb_des (Moerl et al., 2012)
            %DeltaW_limb_des - width of normalized bell curve in descending branch of force-length relationship
            DeltaW_limb_des=0.35;
            % DeltaW_limb_asc - width of normalized bell curve in ascending branch of force length relationship
            DeltaW_limb_asc=0.35;
            %A_rel0 - parameter for contraction dynamics: maximum value of A_rel
            A_rel0=0.25;
            %    B_rel0 -  parameter for contraction dynmacis: maximum value of B_rel
            B_rel0=2.25;
            %    S_eccentric - relation between F(v) slopes at v_CE=0 (van Soest & Bobbert, 1993)
            S_eccentric=2;
            %    F_eccentric - factor by which the force can exceed F_isom for large eccentric velocities (van Soest & Bobbert, 1993)
            F_eccentric=1.5;
            % exponent for descending branch (Moerl et al., 2012)
            v_CElimb_des = 1.5; 
            % exponent for ascending branch (Moerl et al., 2012)
            v_CElimb_asc = 3.0;
            
            %changing variables
            l_CE=@(t,x) L_0-x(1);
            dot_l_CE=@(t,x) x(2);
                %muscles activity, 0<=q<=1
            r_activation=200;
            q=@(t,x)min(r_activation*t,1);
            %force calculations
            %chill initial stuff
            l_PEE0 = L_PEE0*l_CEopt;       % rest length of PEE (Guenther et al., 2007)
           
            K_PEE= F_PEE*( F_max/ ( l_CEopt*(DeltaW_limb_des+1-L_PEE0) )^v_PEE );%I think I have a units issue here??? This is the same equation so...
            % Isometric force (Force length relation)
            %Guenther et al. 2007
          %descending branch or ascending branch
            F_isom = @(t,x) exp( - ( abs( ((l_CE(t,x)/l_CEopt)-1)/DeltaW_limb_des ) )^v_CElimb_des )*(l_CE(t,x) >= l_CEopt)+ exp( -( abs( ((l_CE(t,x)/l_CEopt)-1)/DeltaW_limb_asc ) )^v_CElimb_asc )*(l_CE(t,x)<l_CEopt) ;

                                        % factor of non-linearity in F_PEE (Guenther et al., 2007)
            % Force of the parallel elastic element
            F_PE =@(t,x) K_PEE*(l_CE(t,x)-l_PEE0)^(v_PEE)*(l_CE(t,x) >= l_PEE0);
       
            
            % Hill Parameters concentric contraction
            A_rel=@(t,x)1*(l_CE(t,x)<l_CEopt)+ (l_CE(t,x)>=1)*(F_isom(t,x));
            
            A_rel = @(t,x) (A_rel(t,x) * A_rel0*1/4*(1+3*q(t,x)))*(dot_l_CE(t,x)<=0)+(dot_l_CE(t,x) > 0)*(-F_eccentric*q(t,x)*F_isom(t,x));
            
            B_rel = @(t,x)B_rel0*1*1/7*(3+4*q(t,x));

            B_rel=@(t,x) B_rel(t,x)*(dot_l_CE(t,x)<=0)+(dot_l_CE(t,x) > 0)*((q(t,x)*F_isom(t,x)*(1-F_eccentric)/(q(t,x)*F_isom(t,x)+A_rel(t,x))*B_rel(t,x)/S_eccentric));
            
            % Contractile element force (isometric)
            % I think this might just be to initialize F_CE = F_max*q*F_isom;
            F_CE = @(t,x) F_max*(( (q(t,x)*F_isom(t,x)+A_rel(t,x)) / (1 - dot_l_CE(t,x)/(l_CEopt*B_rel(t,x)) ) )-A_rel(t,x));
            %needed for constructor
            Force=@(t,x)F_CE(t,x)+F_PE(t,x);
            max_force = F_max;
            
            range=0.3*L_0-L_initial;
            %it needs to mean stretch minus 0.7 rest length. What's initial stretch?
        %     motor.range=muscle_length;
            velocity=@(t,x) dot_l_CE;
            
            % call parent constructor
            obj = obj@Motor(max_force, range, velocity, Force);
        end  
    end
end
%% Citations
% Rosario MV, Sutton GP, Patek SN, Sawicki GS. 2016 Muscle�spring dynamics in time-limited, elastic movements.
%   Proc. R. Soc. B 283: 20161561. http://dx.doi.org/10.1098/rspb.2016.1561
%Guenther, M, S Schmitt, and V Wank. 2007. 'High-Frequency Oscillations as
% a Consequence of Neglected Serial Damping in Hill-Type Muscle Models.'
% Biological Cybernetics 97 (1) (July): 63-79.
% doi:10.1007/s00422-007-0160-6.

