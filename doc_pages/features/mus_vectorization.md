title: Vectorization
author: Jiaxing Qi
date: Feb 26, 2016

# Vectorization

## Introduction

Today’s CPUs are highly parallel processors with different levels of parallelism.
We find parallelism everywhere from the parallel execution units in a CPU core,
up to the SIMD (Single Instruction, Multiple Data) instruction set and the
parallel execution of multiple threads. The use of the Intel SIMD instruction
set (e.g. SSE, AVX), which is an extension to the x86 architecture, is called
vectorization.

An easy way to vectorize a Fortran code is to let compiler perform
auto-vectorization. Afterwards, one can examine the compiler report and make
sure that the compiler really does the correct vectorization.

In *Musubi*, one must use the Structure of Array (SOA) data layout and the
computer kernel must be implemented in a Block Loop pattern. Such an example can
be found in [[mus_d3q19_module:mus_advRel_kFluid_rBGK_vStd_lD3Q19]].

Modern CPU often equipped with a cache structure. The data transfer between CPU
and memory has to go through the cache. Normally, when the CPU wants to write
data to the memory, if the data is not in the cache at that moment, the CPU will
first load the data from memory to cache, then modify its value and finally
transfer data from cache to memory back. During this writing process, data were 
transfer through memory interface twice.

**Non-temporal stores** (also called "**streaming stores**")
do not require a prior cache line read for ownership (RFO) but write to memory
"directly". In this way the writing process requires data transfer only **once**!
Under such situation, the memory transfer requirement per element in LBM becomes
\( 8 \times 19 \times 2 + 4 \times 19 = 380 Bytes \), for the D3Q19 stencil.

## How to Enable Streaming Stores

On machines with Intel CPU, one needs to use Intel Fortran compiler and fulfill
the following key points:

1.  Instruct the compiler to align the state array which can not reside with a
    derived data type.

        real(kind=rk), allocatable :: state(:,:)
        !IBM* align(32, state)
        !dir$ attributes align : 32 :: state

2.  Padding is needed to insure the size of the array to be a factor of 4.

        remainder = mod( me%nElems_local, 4 )
        me%nSize  = me%nElems_local + mod(4 - remainder, 4)
        allocate( state( me%nSize, 2 ) )

3.  Use compiler directive before the loop to perform streaming store.

        !DIR$ IVDEP
        !DIR$ VECTOR aligned nontemporal (outState)
        do iLink = 1, ...
          outState(iLink) = ...
        end do

Finally, one has to check the compiler report to make sure that the compiler
really generate streaming stores.

    LOOP BEGIN at /zhome/academic/HLRS/pri/iprjiaqi/apes/musubi/build/source/compute/mus_compute_d3q19_module.f90(213,7)
      remark #15388: vectorization support: reference outstate has aligned access
      remark #15388: vectorization support: reference f000 has aligned access
      remark #15388: vectorization support: reference usqn_o1 has aligned access
      remark #15412: vectorization support: streaming store was generated for outstate
      remark #15300: LOOP WAS VECTORIZED remark #15448: unmasked aligned unit stride loads: 2
      remark #15449: unmasked aligned unit stride stores: 1
      remark #15467: unmasked aligned streaming stores: 1
      remark #15475: --- begin vector loop cost summary ---
      remark #15476: scalar loop cost: 10
      remark #15477: vector loop cost: 2.250
      remark #15478: estimated potential speedup: 4.430
      remark #15479: lightweight vector operations: 9
      remark #15488: --- end vector loop cost summary ---
    LOOP END

    LOOP BEGIN at /zhome/academic/HLRS/pri/iprjiaqi/apes/musubi/build/source/compute/mus_compute_d3q19_module.f90(223,7)
      remark #15388: vectorization support: reference u_x has aligned access
      remark #15388: vectorization support: reference u_y has aligned access
      remark #15388: vectorization support: reference rho has aligned access
      remark #15388: vectorization support: reference outstate has aligned access
      remark #15388: vectorization support: reference f110 has aligned access
      remark #15388: vectorization support: reference usqn_o1 has aligned access
      remark #15412: vectorization support: streaming store was generated for outstate
      remark #15300: LOOP WAS VECTORIZED remark #15448: unmasked aligned unit stride loads: 5
      remark #15449: unmasked aligned unit stride stores: 1
      remark #15467: unmasked aligned streaming stores: 1
      remark #15475: --- begin vector loop cost summary ---
      remark #15476: scalar loop cost: 28
      remark #15477: vector loop cost: 6.750
      remark #15478: estimated potential speedup: 4.130
      remark #15479: lightweight vector operations: 27
      remark #15488: --- end vector loop cost summary ---
    LOOP END
