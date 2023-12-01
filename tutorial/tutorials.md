Tutorials {#mus_tutorials_overview}
========

## Tutorials

The Tutorials show you how to install, configure and set up Musubi. You are guided through all required steps from generating meshes, configuring the solver and post-processing the results.
You can choose a special topic or you can go through the tutorials step by step. 
 
- [Tutorial 00 - Prerequisites](tut_00_prerequisites.html)       
    What you need before you can start using Musubi?
    
- [Tutorial 01 - Configuration](tut_01_mus_config.html)  
    Basic configuration of Musubi simulations
    
- [Tutorial 02 - Toolchain](tut_02_mus_toolchain.html)  
    Using the complete APES simulation tool-chain. Generation of Meshes, Running the solver and post-processing the results
    
- [Tutorial 03 - Tracking](tut_03_tracking.html)  
    Tracking flow quantities during the simulation
    
- [Tutorial 04 - Boundary Conditions](tut_04_boundaries.html)  
    Setting Boundary Conditions
    
- [Tutorial 05 - Restart](tut_05_restart.html)  
    Restarting Simulations
    
- [Tutorial 06 - Initial Colditions](tut_06_initial.html)  
    Setting Initial Conditions ( Equilibrium and Non-equilibrium )
    
- [Tutorial 07 - Stopping criteria](tut_07_convergence.html)  
    Conditions for convergence and stopping the code in general

- [Tutorial 08 - Source terms](tut_08_source.html)

- [Tutorial 09 - Multilevel Simulations](tut_09_mus_multilevel.html)


## How to use the Tutorial ## {#howto}

You should go through this tutorial in a linear fashion.
Chapter 1 deals with the installation of Musubi. 
The following chapters dive deeper into the most important features and 
should give a good impression of what can be done with Musubi, and how.

After installation of Musubi, you will find a `tutorial` folder inside 
your `musubi` folder. 
This folder contains complete configuration files for each chapter, so 
you have working examples of everything we are going to build in this 
tutorial.

## Overview ## {#tut_overview}

- [What you need before you can start using Musubi](@ref tut_0prerequisites)
- [Basic configuration of Musubi simulations](@ref tut_1mus_config)
- [Using the complete APES simulation tool-chain. Generation of Meshes, Running the solver and post-processing the results](@ref tut_1mus_config)
  + [Tracking flow quantities during the simulation](@ref tut_3tracking)
  + [Setting Boundary Conditions](@ref tut_4boundaries)
  + [Setting Initial Conditions ( Equilibrium and Non-equilibrium )](@ref tut_6initial)
  + [Conditions for convergence and stopping the code in general](@ref tut_7convergence)
  + [Restarting Simulations](@ref tut_5restart)
- [Multilevel Simulations](@ref tut_00multilevel)


## Required Knowledge ## {#requ_know}

We assume that you know how to install programs on your computer 
and how to use a console. 
Also, we assume that you have some basic idea of the Lattice Boltzmann
 Method (see [here](http://www.scholarpedia.org/article/Lattice_Boltzmann_Method)
or [here](http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods)
for a brief description of it).
- Find requirements for building, running and visualizing  under
  \ref tem_requirements

## Test cases  

We provide several test cases to practice your skills with the Apes suite. We are working at the tutorials at the moment.

- [Flow around a cylinder](@ref tut_FlowAroundCyl2D)
