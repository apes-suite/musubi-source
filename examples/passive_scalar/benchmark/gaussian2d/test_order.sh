#! /bin/bash

for tau in 0.501 0.503 0.505 0.51 0.55 0.8 2 5
do
    for u in 0.001 0.005 0.01 0.05 0.1
    do
        for order in 'first' 'second'
        do
            echo "tau = $tau" > arg_given.lua
            echo "u_field = $u" >> arg_given.lua
            echo "order = '$order'" >> arg_given.lua

            musubi
            lua calculateD.lua >> order.res
        done
    done
done