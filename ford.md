project: halo
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/halo
summary: Halo orbit solver with modern Fortran
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            json_module:https://jacobwilliams.github.io/json-fortran/
            ddeabm_module:https://jacobwilliams.github.io/ddeabm/
            pyplot_module:https://jacobwilliams.github.io/pyplot-fortran/
            bspline_module:https://jacobwilliams.github.io/bspline-fortran/
            nlesolver_module:https://jacobwilliams.github.io/nlesolver-fortran/
            numerical_differentiation_module:https://jacobwilliams.github.io/NumDiff/
            fortran_astrodynamics_toolkit:https://jacobwilliams.github.io/Fortran-Astrodynamics-Toolkit/

{!README.md!}