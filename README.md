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
writeInfoToFile.m
/components-library


### To do

- [ ] perform sensitivity analysis with code Andres sent
- [ ] add in friction
- [ ] generate plots for paper showing trade-offs
- [ ] generate LaMSA zone plots
- [ ] write up model (including why y_L' = $\tan \theta$)
- [ ] add recommended/Default parameters to GUI


## How to Use this Code

Message from Andres Cook:

Congratulations, LaMSA person! You’ve volunteered for the coolest PoSM lab project. Of course, you’re probably just here because something stopped working, which is probably because of something I did, so…. Oops. Hopefully this documentation will at least be a good starting point, but if not, flag down Prof. Ilton and/or email me at acook@hmc.edu (or Facebook message me tbh)

-Andres Cook ‘21

P.S. from Jackson Castro:

We've added some new features! As well as a ton of more ouputs to sol, we have also incorporated friction and reorganized a lot of stuff in the new /components-library. Further updates to this directory will allow for even easier calls to different types of system setups. As of this moment in time, the simulation seems to be working pretty well, but if for some reason something is giving an odd output, feel free to contact me at jcastro@g.hmc.edu or on Facebook messenger. Best of luck with the simulations!

-Jackson Castro '22


### /components-library
The /components-library directory stores several structs used to call solve_model. In this directory, there are the subdirectories of /springs, /motors, /load-masses, and /latches. Within these subdirectories are structs for varying types of each of these parameters, such as specific latch shapes, varying motor and spring types, as well as a basic one for the load mass. These structs are then called in whatever specified combination in solve_model and heatmap_demo. In heatmap_demo there already exists code for each of these structs that can be commented out or commented back in to account for the various scenarios.

### How to Call the Model
The main function to run one simulation is solve_model.m, which has the following argument syntax:

``` matlab
[sol,transition_times]=solve_model(loading_motor,unlatching_motor,load,latch,spring, outputDirectory)
```

##### Input

Name	           |         Type          | Description
---------------- | --------------------- | -----------------
loading_motor	   | struct with 5 fields	 | Loading motor Parameters
unlatching_motor | struct with 5 fields  | Unlatching motor Parameters
load             | struct with 1 field	 | Load mass
latch            | struct with 5 fields	 | Latch Parameters
spring           | struct with 2 fields  | mass and force function of spring
outputDirectory  |         string        | Name of the output subdirectory <br> to export to or create

##### Output

Name              |	Type      |	Description
---------------   | --------  | -------------------
transition_times	|  1x2	    | Unlatching and launch times
sol	              |  manyx11  | t, y, ydot, x, xdot, <br> normal and frictional components, <br> spring force, and unlatching force <br> column vectors

### Output Files
Currently, solve_model will create two output files every time it is run. One is a .json file that contains the parameters used for this specific call as well as the transition_times. It also outputs a .csv file with all of the data found in sol. Both of these files can be found in the outputDirectory specified in the call, and both the file names and the outputDirectory are organized with a timestamp. However, when heatmap_demo is run, a single output directory is created at the time heatmap_demo is called, and then all of the output from the looping solve_model calls are stored in this directory.

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
+ Write output files using writeInfoToFile.m
  + Create a .json parameter file and .csv output file to be stored in outputDirectory
