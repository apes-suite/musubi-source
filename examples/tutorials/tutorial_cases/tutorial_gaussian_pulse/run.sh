#!/bin/bash

# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking restart

# Create directories for Seeder and Musubi output
mkdir tracking restart 

# Run Musubi
$musubi_path musubi.lua
