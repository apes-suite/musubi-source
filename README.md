Musubi
======

Musubi is a MPI parallel Lattice-Boltzmann solver.
It utilizes treelm to represent meshes and allows for local refinements.

This repository only contains the sources of musubi, which can not be
built on its own.
Supporting libraries like [treelm](https://github.com/apes-suite/tem-source)
and the [build infrastructure](https://github.com/apes-suite/apes-bin) are
gathered in the [Musubi repository](https://github.com/apes-suite/musubi).

Please refer to that container repository for instruction on how to build Musubi.


Organization in Git Submodules
------------------------------

This repository provides the sources of Musubi, but to compile these sources
we need some dependencies and the build infrastructure.
These are gathered in a git supermodule besides the sources of this repository.
To ease the work with this setup, we provide a little script called "request".
The idea is to work on this source repository for the most part just, as if
there is no super repository.
Then when there is some pull request to be create and to share the changes,
you simply run `./request`, which takes care of dealing with the super
repository and creates pull requests in both projects accordingly.

Subsequently the request script can be used to update those pull requests as
needed. Have a look into the request script itself for details.

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
