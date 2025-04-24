#!/bin/bash

# path to seeder executable
seeder_path=~/apes/seeder/build/seeder
# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

rm *.txt
rm *.dat
rm *.res
rm *.png
rm -rf tracking
rm -rf trajectories
rm -rf mesh

mkdir tracking
mkdir trajectories
mkdir -p mesh

$seeder_path seeder.lua

mpirun -np 2 $musubi_path musubi.lua

mv particle*.dat trajectories
python3 plot_trajectories.py

