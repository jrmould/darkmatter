This is a fortran program, compilable with f77. It calls PGPLOT https://sites.astro.caltech.edu/~tjp/pgplot/
But these calls could be removed, since the output array of masses is sent to fort.93, from where it could be plotted
using other popular packages.

We start with a 1000 cubed map of zeros. "The map."
Two numbers are solicited, the number of time steps, and the degradation of spatial resolution, e.g. 10.
This allows fast running for testing modifications.
The initial radius is the inverse of the number of time steps/4.
(1) We call random numbers x0,y0,z0 for the first PBH.
Does its footprint land on any occupied area of the map?
If so, go back to 1. Record it as a failure. If not, place it on the map by setting its footprint = 1.
Record the radius, advance one time step.
Now, because our simulation is in co-moving coordinates, we actually shrink all previous PBH radii by
as many timesteps as have elapsed since they were formed.
Return to (1).
When the number of failures has exceeded the number of timesteps by a factor of 100, quit and make plots.

The simulation is entirely scale free. No seconds. No centimeters. No masses. 
Only after the mass of Phoebe is set equal to that peak, do all the MLT dimensions appear.