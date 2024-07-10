#! /bin/bash

for collision in "bgk" "trt"
do
for tau in 0.51 0.55 0.6 0.7 0.8 0.9 1 2 3 4 5
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'second'" >> arg_given.lua
    echo "collision = '$collision'" >> arg_given.lua
    echo "sigma0 = 40" >> arg_given.lua

    musubi
    lua calculateD.lua >> trt_comparison.res
done


for tau in 0.51 0.8 2 5
do
    for u in 0.001 0.01 0.1
    do
        echo "tau = $tau" > arg_given.lua
        echo "u_field = $u" >> arg_given.lua
        echo "order = 'second'" >> arg_given.lua
        echo "collision = '$collision'" >> arg_given.lua
        echo "sigma0 = 40" >> arg_given.lua

        musubi
        lua calculateD.lua >> trt_comparison.res
    done
done

echo "\n here begins trt result" >> trt_comparison.res
done