title: Restart
@warning WORK IN PROGRESS @endwarning

navigate: [&larr; Boundary Conditions](tut_04_boundaries.html)
| [Overview](index.html)
| [Initional Conditions &rarr;](tut_06_initial.html)

# Restart

In this tutorial we explain how to restart a simulation from a certain time
point and to make a simulation restartable from the last timestep of the
simulation run.

## How to Make it Possible to Restart the Simulation?

You are preparing your simulation with writing the musubi.lua file and you have
not yet run it.
To make a restart possible you have to do the following:

1.  At first, you have to create the restart folder in your simulation directory
with `mkdir restart`.
2.  Add the following code to your `musubi.lua` file somewhere at the end:
```lua
restart = {
  write = 'restart/',
  time_control = {
    min = { iter = rest_start },
    max = { iter = tmax_phy },
    interval = { iter = interval_phy }
  }
}
```

After that, you are able to test if the restart table is working. For every
timestep, *Musubi* writes each a `.lua` file and a `.lsb` file to `./restart/`.
For example, here is a lua file for a random timestep:

```lua
 binary_name = {
    'restart/channel_9.021E-03.lsb'
}
 solver_configFile = 'musubi.lua'
 mesh = './mesh/'
 weights = ''
 time_point = {
    sim =    9.021097956087902E-03,
    iter = 5,
    clock =  103.734125999999996E-03
}
 nElems = 2048
 nDofs = 1
 solver = 'Musubi_v2.0'
 varsys = {
    systemname = 'fluid_incompressible',
    variable = {
        {
            name = 'pdf',
            ncomponents = 19,
            state_varpos = { 1, 2, 3, 4, 5, 6, 7, 8,
                9, 10, 11, 12, 13, 14, 15, 16,
                17, 18, 19 }
        }
    },
    nScalars = 19,
    nStateVars = 1,
    nAuxScalars = 4,
    nAuxVars = 2
}
```
It refers itself to the binary file which is the `.lsb` file.

Now it is possible to restart your simulation from each timestep that is
available in the restart folder.

## Restart your Simulation from any Timestep

If you wish to restart from this random timestep, you will have to change the
restart table in `musubi.lua` to:

```lua
restart = {
  read = 'restart/channel_header_9.021E-03.lua',
  write = 'restart/',
  time_control = {
    min = { iter = rest_start },
    max = { iter = tmax_phy },
    interval = { iter = interval_phy }
  }
}
```

@note `tEnd` and `interval` are variables that are defined at the beginning of
the `musubi.lua` file of the channel testcase. You are free to use different
values in the `time_control` table. As you can see, you have to refer to the
**lua** and not to the **lsb** file in the restart table.

## Restart your Simulation right after the Last Timestep

If you wish to restart your simulation from the last timestep when data was
written to the restart folder, you will have to change:

```lua
read = 'restart/channel_header_9.021E-03.lua',
```
to
```lua
read = 'restart/channel_lastHeader.lua',
```

*Musubi* uses a pretty easy syntax in this case for restart files.
The filename is "simulation_name"+"lastHeader.lua".

Now you have learned how to restart a simulation from any time point you wish
to. This is very useful while you make use of expensive ressources like
supercomputers or while it takes a long time to simulate to time x and
you wish to go on from time x to time y.

Next chapter: [Initial Conditions &rarr;](tut_06_initial.html)
