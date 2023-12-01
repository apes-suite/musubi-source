Download, build and run Musubi {#tut_0prerequisites}
========

Navigate: \ref mus_tutorials_overview "Overview"
| \ref tut_1mus_config "Set up your first simulation with Musubi &rarr;"

\note We assume that you are using a UNIX-like system.
If you run Windows, some commands might be different.

## Download ## {#tut_download}

We use *Mercurial (hg)* for revision control. You need to have it
installed on your system in order to download Musubi.
Follow [these instructions](http://mercurial.selenic.com/
"Mercurial website") if you haven't, yet.
Then, create a directory for everything that happens in these tutorials
(we will call this directory `apes`, but you can use a different name.
Inside this directory, clone the Musubi repository from
`https://geb.inf.tu-dresden.de/hg/musubi` by running
\code
hg clone https://geb.inf.tu-dresden.de/hg/musubi
\endcode
in your console.
If this worked, you have an up-to-date copy of the musubi source code,
which we will compile now.

## Build ## {#tut_build}

We use the *waf* build system, you can learn more about it
[from its website](https://code.google.com/p/waf/ "waf website").
Also, you need *MPI* installed on your system, see for example the
[OpenMPI website](http://www.open-mpi.org/) for instructions.
Finally, you need to set environment variables `FC` and `CC` in order to
assign the correct Fortran and C compilers to *waf*.
The compilers should point to the MPI wrappers from your MPI installation.
Typically these are `mpif90` (Fortran) and `mpicc` (C).
Set them with the bash commands
\code
export CC=mpicc; export FC=mpif90
\endcode
Once you have done all this, navigate to your `musubi` directory and use the
command
\code
./waf configure
\endcode
to configure the compilation.

We are now ready to compile Musubi. Run
\code
./waf build
\endcode
to get a Musubi executable in the *build* subdirectory.
If the compilation finishes without errors, you have Musubi ready to run your
first test case!

## Run ## {#tut_run}

To check your Musubi-installation, navigate to your `musubi` directory and
create a required tracking directory for the output with
\code
mkdir tracking
\endcode

Execute Musubi with
\code
./build/musubi
\endcode
which should result in loads of output, ending with a message similar to
\verbatim
Done with Musubi in [s]    1.864378E-01
\endverbatim
indicating that everything worked fine.
If you get any errors up to here, read the instructions again and follow
them carefully. You need to have this running before you proceed.

## Mesh Generation and Post-processing ## {#tut_seedharvest}

Most of the tutorials will require creation of a mesh and some
post-processing of the results.
The corresponding tools in the \ref toolchain "APES tool chain" are *Seeder*
for mesh-generation and *Harvester* for post-processing.
The installation procedure for them is very similar to Musubi.
Again, navigate to your `apes` directory and run
\code
hg clone https://geb.inf.tu-dresden.de/hg/seeder
hg clone https://geb.inf.tu-dresden.de/hg/harvester
\endcode
to get fresh copies of Seeder and Harvester. Compile them by running
\code
cd seeder
./waf configure build
cd ../harvester
./waf configure build
\endcode
and fix any errors before you proceed.

\note You can make your life easier by adding `apes/seeder/build`,
`apes/harvester/build` and `apes/musubi/build` to your path,
for example by editing  `~/.profile` (MacOS X) or `~/.bashrc` (Unix with
bash) or whatever it is on your system.
In the following tutorials, we assume that you have done just that.
If you have not, you must add the correct paths to the command any time
you try to call Musubi, Harvester or Seeder.

Once you are done with all that, we can start defining our first simulation.

# Troubleshooting

Once you get errors running Musubi, it is possible that something is wrong with your code version. In order to get 
more detail information concerning the errors you can type this command inside the musubi directory which is `apes/musubi/`
as default:
\code
./waf distclean configure debug
\endcode
If you run your simulation once again, you will get more information about the files that cause errors.
> Now, you have to run musubi from a different directory which is `/build/debug/musubi` instead of `/build/musubi`.

Next chapter:
\ref tut_1mus_config "Set up your first simulation with Musubi &rarr;"
