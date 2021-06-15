function [value,isterminal,direction]=direct_actuation_end(t,y,motor)
%End condition for loading
value=motor.Force(t,y);
isterminal=1;
direction=0;
end