title: Musubi-Harvesting

The former *Harvester* as a post-processing tool is not existing anymore.
Instead, there is a new implementation of its functions inside *Musubi* itself.
This one is called *Musubi-Harvesting* or shortly mus_harvesting.

This is why a few things in the handling of post-processing have changed. Inside
the *Musubi* main folder is a sample script that can be used as a template. One
can find it normally here: `~/apes/musubi/harvester.lua`.

In the following, the usage of this tool will be explained.

At first, one can define the name of the simulation that will be a part of the
filename of the produced files.

```lua
simulation_name = 'channel2D'
```

This basic information is read with the subroutine [[tem_load_general]].
It is stored inside `params%solver%simName`. The subroutine is called from
[[mus_hvs_config_load]].

After that, one provide the restart file that one would like to visualize:

```lua
restart = {read = 'restart/channel2D_lastHeader.lua'}
```
This character is stored inside [[tem_restartControl_type]] as 'readFileName'.

Then, one defines the verbosity of logging. [[tem_logging_load_primary]] is
looking for the 'level'. Default is `level = 1` in this case.

```lua
logging = { level = 1 }
```

Then one defines what shall be post-processed with mus_harvesting.
Therefore the `tracking` table is needed.

```lua
tracking = {
  { -- opens sub table for one tracking element
```

To give the resulting files a useful name one can use the `label`. Moreover one
can provide a special output `folder`. For both there are default values.

```lua
label = 'useful_name',
folder = 'here_is_your_output/',
```

But one must provide the output format of the files. There, one can choose
between `vtk`, `ascii` and `ascii_spatial`. To visualize the data with
Paraview one needs `vtk`:

```lua
output = { format = 'vtk' },
```

Then one can post-process every variable that is defined in the variable
system as a default variable or 'special' ones like a difference that is
used in the `variable` table of the *Musubi* config file. For more information
one can have a look at [this page](|temurl|/page/features/variables/index.html).

```lua
variable = {'treeid','solidsphere','density','momentum'},
```

It is possible to define a region from where one can get the data. Therefore
the `shape` table is used. In this example one gets the data for the whole
domain.
Another example of defining the shape can be found in [[tem_load_shape_single]].

```lua
shape = { kind = 'all' },
} -- close tracking sub table
} -- close tracking table
```

One should not forget to close the tracking table. But after that the required
setting is done. One can run the post-processor from the location of the config
file:

```sh
~/apes/musubi/build/mus_harvesting harvester.lua
```

After that, one can visualise the output with Paraview. There is a tutorial for
this software that one can find easily at Youtube.

The next part will be the tutorial for creating a dataset over all time steps.

Therefore one needs two files, the config file template `harvest_series.template`
and the config file for the python script `series.config`.

The template file looks quite similar to the previously introduced
`harvester.lua`.
The first step is to make the use of data from the *Musubi* config file possible.
With `require` one can provide the lua config file. Remember that the filename
without the lua extension is needed here.

Next, one needs the restart input files. Here one makes use of the `$!...!$`
writing that indicates a placeholder that is configured inside the python
script. `$!file!$` is the placeholder for the input files.

After that, one defines the tracking in the same way like mentioned before.
But for the output `folder` one can use the placeholder `$!out!$`.

```lua
require 'musubi'

restart.read = '$!file!$'
tracking = {
  {
    label = 'useful_name',
    variable = {'density'},
    shape = { kind = 'all' }, -- same as: shape = { kind = 'global' }
    folder = '$!out!$',
    output = { format = 'vtk' },
  }
}
```

After that, one configures the python input file that is called `series.config`.
In this case one must define the input files, the output folder, the path to
post-processor executable, the location of the template file, and the path to
the lua executable.

```python
# input files:
files: *_header_*.lua # for all header files
# Lua path:
lua: ~/apes/musubi/build/treelm/aotus/lua
# post processor:
harvester: ~/apes/musubi/build/mus_harvesting
# template file:
template: harvest_series.template
# output folder:
out: ./output/
```

One can run this script with:

```sh
~/apes/musubi/treelm/peons/harvest_series.py -c series.config
```

Another example of how to configure the harvester.lua file:

```lua
variable = {
  {
    name = 'velX',
    ncomponents = 1,
    vartype = 'operation',
    operation = {kind='extract',
                 input_varname = {'velocity'},
                 input_varindex = {1}
    }
  }
},
tracking = {
  {
    label = 'vtk',
    variable = {'velX'},
    shape = {kind = 'all'},
    output = {format = 'vtk'}
  }
}
```
