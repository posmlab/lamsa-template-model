function [value,isterminal,direction] = launching_end(t,y,F_s)
%End condition for loading
value=F_s(t,y);
isterminal=1;
direction=0;
end