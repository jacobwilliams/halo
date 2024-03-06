![NumDiff](media/logo.jpg)
============

## Description

Halo orbit solver with modern Fortran

This application can be used to generate long-duration Earth-Moon halo orbits in the ephemeris model. A work in progress.

## Build

The application is built using the [Fortran Package Manager](https://github.com/fortran-lang/fpm) (FPM). FPM will download all the required dependencies and compile the program, and can also be used to run the tests.

```
fpm build --profile release
```

To run an example:

```
fpm run --profile release -- examples/example.json
```

OpenMP can also be employed to speed up the process.
Another example, using the Intel Fortran Compiler with OpenMP:

```
fpm run --profile release --compiler ifort --flag "-qopenmp" -- examples/example.json
```

Or with Gfortran:

```
fpm run --profile release --flag "-fopenmp" -- examples/example.json
```

### Using qr_mumps

The [qr_mumps](https://gitlab.com/qr_mumps/qr_mumps) solver can also be used, which will run much faster than the default LAPACK solver for very large problems. For this, the code must be compiled with the `WITH_QRMUMPS` preprocessor directive. See the `build-mac.sh` and `run.sh` scripts for details.

### Building On Windows

Work in progress....

Example:

```
"<intel install dir>\windows\bin\ifortvars.bat" intel64
"<intel install dir>\windows\mkl\bin\mklvars.bat" intel64

set OMP_NUM_THREADS=12

fpm run --profile release --compiler ifort --flag "/fpp /Qmkl:parallel /Qopenmp" -- examples/example.json
```

Result is:
```
LINK : fatal error LNK1181: cannot open input file 'lapack.lib'
```

It's because `fmin` and `nlsolver-fortran` both have explicit `lapack` and `blas` dependencies specified in their `fpm.toml` files (`link = ["lapack", "blas"]`). If i delete those it compiles.

### Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/halo/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### See also

 * J. Williams, [Near Rectilinear Halo Orbits](https://degenerateconic.com/near-rectilinear-halo-orbits.html), July 31, 2022 [degenerateconic.com]
 * [qr_mumps](https://gitlab.com/qr_mumps/qr_mumps) -- a software package for the solution of sparse, linear
  systems on multicore computers.

### Reference

 * J. Williams, D. E. Lee, R. J. Whitley, K. A. Bokelmann, D. C. Davis, C. F. Berry, "[Targeting Cislunar Near Rectilinear Halo Orbits for Human Space Exploration](https://www.researchgate.net/publication/322526659_Targeting_Cislunar_Near_Rectilinear_Halo_Orbits_for_Human_Space_Exploration)", 27th AAS/AIAA Space Flight Mechanics Meeting, Feb. 2017
 7. J. Williams, R. D. Falck, and I. B. Beekman, "[Application of Modern Fortran to Spacecraft Trajectory Design and Optimization](https://ntrs.nasa.gov/api/citations/20180000413/downloads/20180000413.pdf)", 2018 Space Flight Mechanics Meeting, 8–12 January 2018, AIAA 2018-145.
 * E. M. Z. Spreen, "[Trajectory Design and Targeting For Applications to the Exploration Program in Cislunar Space](https://hammer.purdue.edu/articles/thesis/Trajectory_Design_and_Targeting_For_Applications_to_the_Exploration_Program_in_Cislunar_Space/14445717)", PhD Dissertation, Purdue University School of Aeronautics and Astronautics, May 2021