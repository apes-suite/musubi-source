#! /bin/bash

for F0 in 0.1 0.5 1 5 10 50 100
do
        # Write the value of F0 to args.lua
        echo "F0 = $F0" > args.lua
        # Run Musubi in the python virtual environment
        mpirun --use-hwthread-cpus -np 30 musubi
done

# Run the post-processing script to generate contour plots
# numpy and matplotlib are required
python plot_contour.py