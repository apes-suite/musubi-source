#!/bin/bash

# path to seeder executable
seeder_path=~/apes/seeder/build/seeder
# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories and database
rm -rf mesh tracking harvest restart output *.db

# Remove old log files
rm timing.res runMusubi.log harvest_series.lua

# Create directories for Seeder and Musubi output
mkdir mesh harvest tracking restart output

# Run printParams.lua to print informations to screen
lua printParams.lua

# Run Seeder
$seeder_path seeder.lua

# Run Musubi
mpirun --oversubscribe -np 8 $musubi_path musubi.lua | tee runMusubi.log

# Generate mesh vtu file
~/apes/seeder/build/sdr_harvesting sdr_harvester.lua

# Create 2D Plots using Gleaner
python plot_track.py
