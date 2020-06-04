# MATLAB LaMSA Model

LaMSA toy model implemented in MATLAB. 

Run heatmap_demo.m for an example of usage.

get_metrics.m is a helper function to extract kinematic metrics (e.g. vmax, amax) from the solution (t,x,v)


Code to numerically integrate LaMSA model is primarily in function solve_model.m

Helper functions to solve_model.m:
launching_ode.m
launching_end.m
unlatching_ode.
unlatching_end.m



LaMSA Code May 2020 Updates
---------------------------

### To do

- [X] put axis labels on heatmap function
- [ ] perform sensitivity analysis with code Andres sent
- [X] test unlatching force 
- [ ] add in friction
- [ ] generate plots for paper showing trade-offs
- [ ] generate LaMSA zone plots
- [ ] write up model (including why y_L' = $\tan \theta$)

### May 28
+ improved initial guess for fzero function call to determine quasi-static loading in load phase solve_model.m
+ y_L is now a 3 element cell array: y_L{1} is latch shape function, y_L{2} is first derivative of latch shape function, y_L{3} is the second derivative of latch shape function
+ verified consistency with Science paper results for rounded latch
+ made tolerances smaller for ode45 call in unlatching phase of solve_model.m to prevent ode integration from over-stepping end of latch (which lead to complex valued results). Created a warning if complex valued results appear
