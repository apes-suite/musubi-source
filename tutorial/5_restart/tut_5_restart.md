Restarting Simulations {#tut_5restart}
========

Navigate: \ref tut_4boundaries "&larr; Definition of Boundary Conditions"
| \ref mus_tutorials_overview "Overview"
| \ref tut_6initial "Definition of Initional Conditions &rarr;"

In this tutorial we explain how to restart a simulation from a certain time point and to make a simulation
restartable from the last time step of the simulation run.

## How to make it possible to restart the simulation?

You are preparing your simulation with writing the musubi.lua file and you have not yet run it.
To make a restart possible you have to do the following:

1.  At first you have to create the restart folder in your simulation directory with `mkdir restart`.
2.  Add the following code to your musubi.lua file somewhere at the end:
    \code
    restart = {
      write = 'restart/',
      time_control = { 
        min = {iter = tEnd/4},
        max = {iter = tEnd}, 
        interval = {iter = interval}
      }
    }
    \endcode
After that, you are able to test if the restart table is working. For every time step, Musubi writes 
each a .lua file and a .lsb file to `./restart/`. For example, here is a lua file for a random time step:
\code
binary_name = {
    'restart/lbm_incompbgkacoustic_106.250E-03.lsb' 
}
 mesh = './mesh/'
 time_point = {
    sim =  106.250000000000053E-03,
    iter = 34,
    clock =  409.535169601440430E-03 
}
 nElems = 2048
 nDofs = 1
 solver = 'Musubi_v1.0'
 varsys = {
    {
        systemname = 'lbm_incompbgkacoustic',
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
        nStateVars = 1 
    } 
}
\endcode
It refers itself to the binary file which is the .lsb file.

Now it is possible to restart your simulation from each time step that is available in the restart folder.

## Restart your simulation from any time step

If you wish to restart from this random time step, you will have to change the restart table in musubi.lua to:
\code
restart = {
  read = 'restart/channel_header106.250E-03.lua',
  write = 'restart/',
  time_control = {
    min = {iter = tEnd/4},
    max = {iter = tEnd}, 
    interval = {iter = interval}
  }
}
\endcode

\note `tEnd` and `interval` are variables that are defined at the beginning of the musubi.lua file of the
channel testcase. You are free to use different values in the time_control table. 
As you can see, you have to refer to the lua and not to the lsb file in the restart table.

## Restart your simulation right after the last time step

If you wish to restart your simulation from the last time step data was written to the restart folder,
you will have to change 
\code
read = 'restart/channel_header106.250E-03.lua',
\endcode
to
\code
read = 'restart/channel_lastHeader.lua',
\endcode

Musubi uses a pretty easy syntax in this case for restart files. 
The filename is "simulation_name"+"lastHeader.lua".

Now you have learned how to restart a simulation from any time point you wish to. This is very useful while
you make use of expensive ressources like supercomputers or while it takes a long time to simulate to time x and
you wish to go on from time x to time y.

Next chapter: \ref tut_6initial "Definition of Initial Conditions &rarr;" 
