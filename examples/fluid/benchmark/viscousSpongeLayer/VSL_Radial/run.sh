#!/bin/bash

# path to seeder executable
seeder_path=~/apes/seeder/build/seeder
# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking tracking_vtk mesh restart

# Create directories for Seeder and Musubi output
mkdir mesh tracking tracking_vtk restart

# Run Seeder
$seeder_path seeder.lua

# Run Musubi
mpirun --oversubscribe -np 8 $musubi_path musubi.lua

# plot 
python plot_track.py
