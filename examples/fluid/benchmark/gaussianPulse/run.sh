#!/bin/bash
#
# Run from with apes-pyenv activated

# Remove old directories
rm -rf tracking restart

# Create directories for Seeder and Musubi output
mkdir tracking restart

# Run Musubi
mpirun --oversubscribe -np 2 musubi

# Remove database if existent
rm *.db
# Create 3D Plots using Gleaner
python plot_track.py
# List the created plots
echo "List of created plots:"
ls *.png
