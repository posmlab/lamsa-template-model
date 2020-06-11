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


### To do

- [ ] perform sensitivity analysis with code Andres sent
- [ ] add in friction
- [ ] generate plots for paper showing trade-offs
- [ ] generate LaMSA zone plots
- [ ] write up model (including why y_L' = $\tan \theta$)


## How to Use this Code

Message from Andres Cook:

Congratulations, LaMSA person! You’ve volunteered for the coolest PoSM lab project. Of course, you’re probably just here because something stopped working, which is probably because of something I did, so…. Oops. Hopefully this documentation will at least be a good starting point, but if not, flag down Prof. Ilton and/or email me at acook@hmc.edu (or Facebook message me tbh)

-Andres Cook ‘21

### How to Call the Model
The main function to run one simulation is solve_model.m, which has the following argument syntax:

``` matlab
[sol,transition_times]=solve_model(F_in,F_out,F_s,m_eff,m_L,t_L_guess,v_0L,y_L)
```

##### Input

Name	|            Type           	| Description
----- | --------------------------- | -----------------
F_in	| Function (num, 1x3) -> num	| Loading motor force
F_out	| Function (num, 1x3) -> num	| Unlatching motor force
F_s	  | Function (num, 1x3) -> num	| Spring force
m_eff |           	num	            |   Load mass
m_L	  |             num	            |  Latch mass
t_L_guess |	        num	            | Upper bound on unlatching time
v_0L  |	            num	            | Initial latch removal velocity
y_L	  | 1x3 cell array of functions | y_L{1}	Latch shape y(x) <br> y_L{2} Latch slope y'(x) <br> y_L{3} Latch concavity y''(x)

##### Output

Name              |	Type      |	Description
---------------   | --------  | ---------------------------
transition_times	|  1x2	    | Unlatching and launch times
sol	              |  many x3	| Time, position, and velocity column vectors

##### Model Structure
The main function solve_model is generally structured as follows
+ Loading: use a root-finding function to determine where the total force on the load during loading is 0. This can be tricky, and often throws errors
+ Unlatching
  + Create an anonymous function wrapping unlatching_ode.m into a function that only takes time and state variables (position, velocity)
  + Perform a similar step for unlatching_ode.m
  + Set up the ODE solver to stop running when the value of unlatching_end changes sign (when the normal force becomes less than zero)
  + Using ode45, solve the unlatching ODE
  + The unlatching ODE is for the horizontal position of the latch, so use y_L to convert from x to y.
+ Launch
  + Essentially the same as loading, but using launching_ode.m and launching_end.m
+ Concatenate the solution arrays into one long one
