Musubi
======

Musubi is a MPI parallel Lattice-Boltzmann solver.
It utilizes treelm to represent meshes and allows for local refinements.
Treelm is incorporated as a subrepository.

Compilation is done via waf and achieved by:

    export FC=mpif90
    ./waf configure build

Note that MPI is required to compile Musubi.

License
-------

Musubi is licensed under the terms of the 2-clause BSD license reproduced below.
This means that Musubi is free software and can be used, reproduced, modified,
distributed and redistributed also for commercial purposes under the conditions
of the BSD license.
The only requirement is that some credit to the authors is given by putting this
copyright notice somewhere in your project.
See individual source files for copyright holders.

According to good scientific practice, publications on results achieved in whole
or in part due to Musubi should cite at least one paper presenting the Musubi
software.

An appropriate reference could be:
@article{Hasert:94TZ3_GF,
author = {Hasert, Manuel and Zimny, Simon and Masilamani, Kannan and Qi, Jiaxing and Klimach, Harald and Bernsdorf, J{\"o}rg and Roller, Sabine},
title = {{Complex Fluid Simulations with the Parallel Tree-based Lattice Boltzmann Solver Musubi}},
journal = {J. Comp. Sci.},
year = {submitted 2013},
month = mar
}
