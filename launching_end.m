function [value,isterminal,direction] = launching_end(t,y,spring)
%End condition for loading
value=spring.Force(t,y);
isterminal=1;
direction=0;
end