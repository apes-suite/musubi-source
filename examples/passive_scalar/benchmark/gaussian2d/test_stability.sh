#! /bin/bash

for tau in $(seq 0.5001 0.0002 0.505)
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'first'" >> arg_given.lua

    musubi
    lua calculateD.lua >> stability.res
done

for tau in $(seq 0.507 0.002 0.8)
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'first'" >> arg_given.lua

    musubi
    lua calculateD.lua >> stability.res
done

for tau in $(seq 1 1 5)
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'first'" >> arg_given.lua

    musubi
    lua calculateD.lua >> stability.res
done
