title: Prerequisites
@warning WORK IN PROGRESS @endwarning

Navigate: [Overview](index.md)
| [Configuration &rarr;](tut_01_mus_config.md)

# Prerequisites
Download, build and run *Musubi*.

@note We assume that you are using a UNIX-like system.
If you run Windows, some commands might be different.

## Download ## {#tut_download}

We use *git* with submodules for revision control.
You need to have it installed on your system in order to download *Musubi*.

Create a directory for everything that happens in these tutorials
(we will call this directory `apes`), but you can use a different name.
Inside this directory, clone the *Musubi* repository from
`https://github.com/apes-suite/musubi.git` by running
```sh
git clone --recurse-submodules https://github.com/apes-suite/musubi.git
```
in your console.
If this worked, you have an up-to-date copy of the *Musubi* source code,
which we will compile now.

If you already cloned *Musubi* you can update *Musubi* via:
```sh
git pull --recurse-submodules
```

## Build ## {#tut_build}

We use the *waf* build system, you can learn more about it
[from its website](https://waf.io/).
Also, you need *MPI* installed on your system, see for example the
[OpenMPI website](http://www.open-mpi.org/) for instructions.
Finally, you need to set environment variables `FC` and `CC` in order to
assign the correct Fortran and C compilers to *waf*.
The compilers should point to the MPI wrappers from your MPI installation.
Typically these are `mpif90` (Fortran) and `mpicc` (C).
Set them with the bash commands
```sh
export CC=mpicc; export FC=mpif90
```

If you updated *Musubi* you need yo clean the old build and the coco
preprocessor. If you cloned *Musubi* there is no need to do so.
```sh
bin/waf cleanall
```

Once you have done all this, navigate to your *Musubi* directory and use the
command
```sh
bin/waf configure
```
to configure the compilation.

We are now ready to compile *Musubi*. Run
```sh
bin/waf build
```
to get a *Musubi* executable including *Mus_Harvester* for post-processing in
the *build* subdirectory. To compile musubi only add the following argument:
```sh
--target="musubi"
```
If the compilation finishes without errors, you have *Musubi* ready to run your
first test case!

## Run ## {#tut_run}

To check your *Musubi*-installation, navigate to your `musubi` directory and
create a required tracking directory for the output with
```sh
mkdir tracking
```

Execute *Musubi* with
``` sh
./build/musubi
```
which should result in loads of output, ending with a message similar to
``` sh
Done with Musubi in [s]    1.864378E-01
```
indicating that everything worked fine.
If you get any errors up to here, read the instructions again and follow
them carefully. You need to have this running before you proceed.

## Mesh Generation and Post-Processing ## {#tut_seedharvest}

Most of the tutorials will require creation of a mesh and some
post-processing of the results. The corresponding tools in our toolchain
"APES tool chain" are *Seeder* for mesh-generation, for post-processing
*Seeder-Harvesting* and *Musubi-Harvesting*. The installation procedure for
*Seeder* is very similar to *Musubi*. Again, navigate to your `apes` directory
and run

```sh
git clone https://github.com/apes-suite/seeder.git
```
to get a fresh copy of *Seeder*. Compile it by running
```sh
cd seeder
bin/waf configure build
```
and fix any errors before you proceed.

@note You can make your life easier by adding `apes/seeder/build`,
`apes/seeder/build/sdr_harvesting`, `apes/musubi/build` and
`apes/musubi/build/mus_harvesting` to your path,
for example by editing  `~/.profile` (MacOS X) or `~/.bashrc` (Unix with
bash) or whatever it is on your system.
In the following tutorials, we assume that you have done just that.
If you have not, you must add the correct paths to the command any time
you try to call *Musubi*, *Musubi-Harvesting*, *Seeder-Harvesting* or *Seeder*.

Once you are done with all that, we can start defining our first simulation.

## Troubleshooting

Once you get errors running *Musubi*, it is possible that something is wrong
with your code version. In order to get more detailed information concerning the
errors you can type this command inside the `musubi` directory which is
`apes/musubi/` as default:

```sh
bin/waf distclean configure debug
```
If you run your simulation once again, you will get more information about the
files that cause errors.

> Now, you have to run *Musubi* from a different directory which is
> `/build/debug/musubi` instead of `/build/musubi`.


Next chapter:
[Configuration &rarr;](tut_01_mus_config.md)
