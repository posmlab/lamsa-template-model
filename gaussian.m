%% gaussian

% This generates a gaussian bell curve! Use for figures, or whatever you 
% want.

%% code

figure
x = linspace(0,5E4,100);
y = gaussmf(x,[0.1*1E4 1E4]);
plot(x,y)