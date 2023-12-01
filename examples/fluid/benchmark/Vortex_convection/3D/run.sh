#!/bin/bash

# path to musubi executable
musubi_path=~/apes/musubi/build/musubi
seeder_path=~/apes/seeder/build/seeder

# Remove old directories
rm -rf tracking restart mesh

# Create directories for Seeder and Musubi output
mkdir  tracking restart mesh

# Run Seeder
mpirun --oversubscribe -np 12 $seeder_path seeder.lua

# Run Musubi
mpirun --oversubscribe -np 12 $musubi_path musubi.lua

