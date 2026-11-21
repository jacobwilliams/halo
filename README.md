![halo](media/logo.jpg)
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
fpm run halo_solver --profile release -- examples/example.json
```

OpenMP can also be employed to speed up the process.
Another example, using the Intel Fortran Compiler with OpenMP:

```
fpm run halo_solver --profile release --compiler ifort --flag "-qopenmp" -- examples/example.json
```

Or with Gfortran:

```
fpm run halo_solver --profile release --flag "-fopenmp" -- examples/example.json
```

### Jacobian options

The solver's Jacobian is computed by finite differences by default (`"jacobian_mode": 1`). With `"jacobian_mode": 2`, it is computed analytically instead: each segment is propagated once with its state transition matrix (using Pines' second partials for the lunar gravity field, in `src/astro_partials_module.f90`), and the epoch partials are computed in closed form. This is faster and more accurate. It requires `"use_splined_ephemeris": true`. Optional `"stm_rtol"` and `"stm_atol"` set the integrator tolerances for the state transition matrix (default: `rtol` and `atol`). See the technical note in `docs/` for details and a comparison.

The tests `partials_test` (the partial derivatives) and `jacobian_test` (both Jacobian methods, compared) are run from the root directory:

```
fpm test partials_test --profile release
fpm test jacobian_test --profile release --flag "-fopenmp" -- examples/example_sparse.json
```

`text_buffer_test` checks the number formatting used to write the trajectory JSON file (every number must read back exactly):

```
fpm test text_buffer_test --profile release
```

`ephemeris_spline_test` times the setup of the splined ephemeris and of the MOON_PA frame (and of the whole problem) for the given numbers of revs (default 100 and 1120), with a config file's ephemeris options:

```
fpm test ephemeris_spline_test --profile release --flag "-fopenmp" -- examples/example_sparse.json 100 1120
```

### Epoch variables

By default (`"epoch_variable_mode": 1`), each node's epoch variable is its epoch (days from the reference epoch), scaled by its own initial magnitude. The scale grows along the trajectory, so the solver's minimum-norm steps change the later epochs most, and the orbit's phase (its periapsis passage times) can drift by many hours from the initial guess. With `"epoch_variable_mode": 2`, the epoch variables are the offsets from the initial guess, all scaled by the same `"epoch_scale"` (days, default 1.0). This keeps the phase closer to the initial guess (for the 100-rev baseline: periapsis times within 4.4 h of the guess, instead of 18 h). A smaller scale holds the phase more tightly, but too small a scale doesn't converge (0.01 days diverges, and 1 day stalls for 1120 revs: use about 10 to 100 days there).

With `"xscale_mode": 2`, all the variables' scales are instead set from their sensitivities: each scale is divided by the norm of its column of the Jacobian at the initial guess, so a unit change of any scaled variable changes the constraints by the same amount (`"epoch_scale_factor"` multiplies the epoch scales: >1 makes the epochs cheaper to change). This is experimental: for the 100-rev baseline it converges as well as the default (periapsis times within 5.6 h of the guess), but it doesn't converge for 1120 revs (with any `epoch_scale_factor` from 1 to 100), where the default scaling does.

### GUI

`halo_gui` is a GUI for the solver ([qt-fortran](https://github.com/jacobwilliams/qt-fortran) for the window, [vkframes](https://github.com/jacobwilliams/vkframes) for the 3D view). Choose a config file (default `examples/example_sparse.json`) and click **Run solver**. The solver runs inside the GUI (in a worker thread), its output is shown in the window, and the resulting trajectory is plotted in the Earth-Moon rotating or inertial frame (Moon-centered, colored by the distance from the Moon). With **Plot each iteration** checked, each solver iteration's trajectory is plotted as it solves (sent straight to the 3D view, not through files). **Center** puts the camera's center (what it rotates about) on the central body or on the spacecraft (the camera then moves with it). **Stop** stops the solver at its next iteration. **File > Open Trajectory** plots an existing `traj_*.json` file.

The GUI is a separate fpm package (`gui/`), so the solver above still builds and runs with plain fpm, without Qt or Vulkan. Qt, VulkanSceneGraph and the C++ compiler come from the [pixi](https://pixi.sh) environment, and `tools/fpm_gui.py` gives fpm the C++ settings both libraries need (macOS needs version 26 or later, for the Vulkan driver):

```
pixi run gui-build                                     # build the GUI (release, OpenMP)
pixi run gui                                           # run it
pixi run gui -- examples/example.json --iterations --run   # check "Plot each iteration", and solve at startup
```

`pixi run gui-build-quad` and `pixi run gui-quad` do the same with the solver in quad precision (`-DREAL128`, much slower). The 3D view stays in double precision.

### Using qr_mumps

The [qr_mumps](https://gitlab.com/qr_mumps/qr_mumps) solver can also be used, which will run much faster than the default LAPACK solver for very large problems. For this, the code must be compiled with the `WITH_QRMUMPS` preprocessor directive. See the `build.sh` and `run.sh` scripts for details.

### Building On Windows

Work in progress....

Example:

```
"<intel install dir>\windows\bin\ifortvars.bat" intel64
"<intel install dir>\windows\mkl\bin\mklvars.bat" intel64

set OMP_NUM_THREADS=12

fpm run halo_solver --profile release --compiler ifort --flag "/fpp /Qmkl:parallel /Qopenmp" -- examples/example.json
```

Result is:
```
LINK : fatal error LNK1181: cannot open input file 'lapack.lib'
```

It's because `fmin` and `nlsolver-fortran` both have explicit `lapack` and `blas` dependencies specified in their `fpm.toml` files (`link = ["lapack", "blas"]`). If i delete those it compiles.

### Other notes

Another source of CR3BP halos that can be used as initial guesses for this tool can be found here: https://ssd.jpl.nasa.gov/tools/periodic_orbits.html

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