title: Abort Criteria
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Initial Conditions](tut_06_initial.md)
| [Overview](index.md)
| [Source Terms &rarr;](tut_08_source.md)

# Abort Criteria

In this tutorial, we show how to use a convergence criterium to
stop the code before reaching the maximum time.
Convergence sensors are actually objects part of the `sim_control` table with
an additional `abort_criteria` table.

```lua
sim_control = {
  time_control = {
    ...
  },
  abort_criteria = {
    ...
    convergence = {
      variable      = {...},
      shape         = {...},
      time_control  = {...},
      reduction     = {...},
      norm          =  ... , -- 'simple' or 'average'
      nvals         =  ... ,
      absolute      =  ... , -- 'true' or 'false'
      condition     = {...}, -- nConditions == nVariables
    }
  }
}
```

@Note See [[tem_convergence_load]] for a set of options and few examples.

Let's consider the channel test case of before again.
If we now want to add a convergence sensor to the simulation, we add the table
`abort_criteria` to `sim_control`.

```lua
--! [Simulation Control]
sim_control = {
  time_control = {
    max = { sim = tmax_phy },
    interval = { sim = interval_phy },
    clock = 3600 --s
   },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
    -- Abort if steady state is reached with condition defined in convergence
    -- table
    steady_state = true,
    -- Convergence condition for steady state check
    convergence = {
      -- Variables to check
      variable = {'pressure_phy'},
      -- Check only point in middle of domain
      shape = {
        kind = 'canoND',
        object = {{origin = {0.,0.,0.} }}
      },
      -- Reduce variables values in space to single average value
      reduction = {'average'},
      -- How often to do convergence?
      time_control = {
        min = tmax_phy/2.0,  -- Start convergence check after half sim. time
        max = tmax_phy,      -- DO convergence until end of simulation
        interval = 10*dt     -- Do convernce check every 10*dt [s]
      },
      norm='average',
      nvals = 50,
      absolute = false,
      -- Condition to statisfy to every variable
      condition = {
        { threshold = 1.e-5, operator = '<=' },
      }
    },
  }
}
--! [Simulation control]
```

In the convergence table you now specify under which condition a simulation is
converged. We already chose the quantity to use as the pressure above. The norm
to evaluate the convergence defines how to compute the converged condition on a
set of available data. The length of this set is defined by `nvals`. This
basically means over how many points in time the convergence is evaluated.
The condition under which convergence is achieved is done with the `condition`
table. You specify a threshold and the operator. You can choose if you want an
`absolute` or a `relative` metric by specifying `absolute = true / false`.

Each point in time when this object is active (i.e. the conditions of the time
definition table are met), a value is added to the convergence data set. The
current convergence quantity is then compared to the norm of the convergence
data set. If the condition is met, convergence is achieved.

> See residual is evaluated within [[tem_convergence_evaluate]]
> with function [[tem_convergence_module:evaluate_residual]].

In our case, the current pressure value is compared to the average over fifty
points in time. If the difference is smaller or equal to 1.e-5, the convergence
is set to be achieved.

Once convergence is achieved, this is communicated to all other processes in
the next interval, when the total density is computed.  Then, the simulation is
terminated after writing out all tracking, restart and output data.


## Automatic and Manual Stopping of Musubi ##

The code might have to stop in other cases. One case is when invalid numbers
are encountered. When the total density within the simulation domain is
computed, *Musubi* checks if invalid numbers are found. If so, the result is
invalid, as the code has just been crashing, and the simulation can be stopped.
Nevertheless, all tracking and restart files are written.

Sometimes computing resources are restricted to a certain time interval only.
It should therefore be ensured that *Musubi* terminates cleanly upon reaching
such a time limit. You can define maximum wall clock limits, upon which
*Musubi* then stops. In the time object, you just have to give the maximum
number of seconds as:

```lua
sim_control = {
  time_control = {
    max = tEnd,
    interval = interval,
    clock = 3600 --s
  },
  abort_criteria = {
    ...
  }
}
```

If you manually want to terminate *Musubi*, you can create a file in the musubi
directory during runtime. It must be named `stop` and does not have to have any
content. Upon next check interval, *Musubi* checks for the existence of the
file. If encountered, *Musubi* terminates cleanly.

Next chapter: [Source Terms &rarr;](tut_08_source.md)
