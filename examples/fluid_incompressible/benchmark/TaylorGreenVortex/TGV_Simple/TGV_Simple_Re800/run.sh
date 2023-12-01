#!/bin/bash

# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking restart vtkfiles *.db

# Create directories for Seeder and Musubi output
mkdir tracking restart vtkfiles

# Run Musubi
mpirun --oversubscribe -np 12 $musubi_path musubi.lua

# Create 3D Plots using Gleaner
python plot_track.py
# List the created plots
echo "List of created plots:"
ls *.png
