Musubi
======

Musubi is an MPI parallel Lattice-Boltzmann solver.
It utilizes the treelm library to represent meshes
and allows for local refinements.

This repository only contains the sources of Musubi and can not be
built on its own.
Supporting libraries like [treelm](https://github.com/apes-suite/tem-source)
and the [build infrastructure](https://github.com/apes-suite/apes-bin) are
gathered in the [Musubi repository](https://github.com/apes-suite/musubi).

Please refer to that container repository and the
[documentation](https://apes-suite.org/musubi/page/)
for instructions on how to build Musubi.


Organization in Git Submodules
------------------------------

This repository provides the sources of Musubi, but to compile these sources
we need some dependencies and the build infrastructure.
These are gathered in a git supermodule besides the sources of this repository.
To ease the work with this setup, we provide a little script called "request"
using the github CLI.

The idea is to work on this source repository for the most part just, as if
there is no super repository.
Then when there is some pull request to be create and to share the changes,
you simply run `../bin/request`, which takes care of dealing with the super
repository and creates pull requests in both projects accordingly.

Subsequently the request script can be used to update those pull requests as
needed.
Have a look into the request script itself for details.
