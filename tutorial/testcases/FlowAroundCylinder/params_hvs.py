## harvester input file template
## the filename for read and folder name for output are
## defined by harvest_series.py so they are fixed to
## ${filename}$ and ${folder}$
template: harvester.template

## input restart file which replaces the string
## ${filename}$ in the harvester.template.
## Can be single file or multiple files
files: restart/Re100/*.lua

## output folder where replaces the string 
out: harvest/

## harvester exectuable
harvester: ~/apes/musubi/build/debug/mus_harvesting

## command to run in parallel
run: mpirun -n 4
