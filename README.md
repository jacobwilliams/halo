# halo

Halo orbit solver with modern Fortran

## Build

```
fpm build --profile release
```

To run an example:

```
fpm run --profile release -- example2.json
```

Another example, using the Intel Fortran Compiler with OpenMP:

```
fpm run --profile release --compiler ifort --flag -qopenmp -- example2.json
```

### See also

 * J. Williams, [Near Rectilinear Halo Orbits](https://degenerateconic.com/near-rectilinear-halo-orbits.html) (degenerateconic.com)
### Reference

 * J. Williams, D. E. Lee, R. Whitley, C. Berry, "[Targeting Cislunar Near Rectilinear Halo Orbits for Human Space Exploration](https://www.researchgate.net/publication/322526659_Targeting_Cislunar_Near_Rectilinear_Halo_Orbits_for_Human_Space_Exploration)", 27th AAS/AIAA Space Flight Mechanics Meeting, Feb. 2017
 7. J. Williams, R. D. Falck, and I. B. Beekman, "[Application of Modern Fortran to Spacecraft Trajectory Design and Optimization](https://ntrs.nasa.gov/api/citations/20180000413/downloads/20180000413.pdf)", 2018 Space Flight Mechanics Meeting, 8–12 January 2018, AIAA 2018-145.
 * E. M. Z. Spreen, "[Trajectory Design and Targeting For Applications to the Exploration Program in Cislunar Space](https://hammer.purdue.edu/articles/thesis/Trajectory_Design_and_Targeting_For_Applications_to_the_Exploration_Program_in_Cislunar_Space/14445717)", PhD Dissertation, Purdue University School of Aeronautics and Astronautics, May 2021