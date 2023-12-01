mkdir mesh
mkdir harvest
mkdir tracking
#mkdir tracking/Re20
mkdir tracking/Re100

rm tracking/*
rm tracking/Re20/*
rm tracking/Re100/*
rm mesh/*
rm harvest/*

module add gcc

# generate mesh
~/apes/seeder/build/seeder seeder.lua

# run simulation
~/apes/musubi/build/musubi musubi.lua | tee musubi.log

# generate VTK files
~/apes/harvester/build/harvester harvester.lua
~/apes/harvester/scripts/harvest_series.py --config params_hvs.py

module rm gcc

# plot the parameters
~/apes/harvester/gleaner/gleaner.py params_plot.py

# compute maximum values for cD and cL
python3 get_coeff_max.py

ls -l filedata/
ls -l harvest/

module add gcc
