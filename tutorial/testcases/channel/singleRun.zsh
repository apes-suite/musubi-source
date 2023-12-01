#!/bin/zsh
#

# delete old files
rm mesh/*
rm tracking/channel*
rm harvest/*
rm restart/*
rm ./*.vtu
rm timing.res
rm runMusubi.log
rm harvest_series.lua

# generate mesh
$seeder

# run simulation
$musubi | tee runMusubi.log

# mv tracking/*hvsXY*lastHeader.lua

# generate vtk files
~/apes/musubi/treelm/peons/harvest_series.py -c series.config

# generate mesh vtu file
~/apes/seeder/build/sdr_harvesting sdr_harvester.lua
