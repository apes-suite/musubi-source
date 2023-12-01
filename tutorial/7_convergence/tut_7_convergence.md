Setting up a Convergence Criterium {#tut_7convergence}
========

Navigate: \ref tut_6initial "&larr; Definition of Initial Conditions"
| \ref mus_tutorials_overview "Overview"
| \ref tut_00multilevel "Multilevel Simulations &rarr;"

In this tutorial, we show how to use a convergence criterium to
stop the code before reaching the maximum time.
Convergence sensors are actually objects part of the *time_control* table with an additional 
convergence table.

\Note See [tem_convergence_load](@ref tem_convergence_module::tem_convergence_load) for a set of options.


Let's consider the channel test case of before again.
If we now want to add a convergence sensor to the simulation, we must create a tracking object.
Here we do that by manually adding a table entry to the already existing tracking table.
Let's call this tracking object convergence and use the pressure as a convergence criterium.

\snippet testcases/channel/musubi.lua convergence

So we see that we are actually just adding a tracking object with its spatial extent 
with the 'shape' table 
and the definition of when to extract data with the 'time' table.

The tracking format has to be set to 'convergence', and the convergence table must be present.
In the convergence table you now specify under which condition a simulation is converged.
We already chose the quantity to use as the pressure above.
The norm to evaluate the convergence defines how to compute the converged condition on a set of available data.
The length of this set is defined by nvals. This basically means over how many points in time the convergence is evaluated.
The condition under which convergence is achieved is done with the condition table. You specify a threshold and the operator.
You can choose if you want an absolute or a relative metric by specifying absolute = true / false.

Each point in time when this tracking object is active (i.e. the conditions of the time definition table are met), 
a value is added to the convergence data set. 
The current convergence quantity is then compared to the norm of the convergence data set.
If the condition is met, convergence is achieved.
See [Evaluate residual](@ref tem_convergence_module::evaluate_residual).

In our case, the current pressure value is compared to the average over fifty points in time.
If the difference is smaller or equal to 1.e-5, the convergence is set to be achieved.

Once convergence is achieved, this is communicated to all other processes in the next interval, when the 
total density is computed. 
Then, the simulation is terminated after writing out all tracking, restart and output data.


## Automatic and Manual Stopping of Musubi ##

The code might have to stop in other cases. One case is when invalid numbers are encountered.
When the total density within the simulation domain is computed, Musubi checks if invalid numbers are 
found. If so, the result is invalid, as the code has just been crashing, and the simulation can be stopped.
Nevertheless, all tracking and restart files are written.

Sometimes computing resources are restricted to a certain time interval only.
It should therefore be ensured that Musubi terminates cleanly upon reaching such a time limit.
You can define maximum wall clock limits, upon which Musubi then stops.
In the time object, you just have to give the maximum number of seconds as wall_max = ...

\snippet testcases/channel/musubi.lua time settings

If you manually want to terminate Musubi, you can create a file in the musubi directory during runtime.
It must be named 'stop' and does not have to have any content.
Upon next check interval, Musubi checks for the existence of the file. 
If encountered, Musubi terminates cleanly.

Next chapter: \ref tut_00multilevel "Multilevel Simulations &rarr;"
