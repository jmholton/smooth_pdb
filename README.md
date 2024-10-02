# smooth_pdbs

turn a series of pdb files into a new series where x,y,z occ and B factor vary smoothly with a specified reaction coordinate

## Motivation

If you are doing time-resolved, multi-temperature, radiation damage or other types of studies when
 you are trying to probe how your structure changes as you gradually adjust your desired independent variable,
 then it would be nice if the dependent variables changed gradually too.  However, if you independently 
 auto-build into each data set, you will find the resulting models have a lot of random differences. Not just
 "noise" in the xyz coordinates, but the refined B factors and occupancies as well. What you really want, is 
 a kind of "joint" refinement, where information from all your data sets influences all your models and 
 carries information from neighboring points into each structure. You also want naming conventions to be
 consistent across the whole series. 
 
 One way to do this is to start at the starting point, build and refine that as best you can, and then 
 drop that model into the next data set as a starting point for a new round of building and refinement. 
 However, this will be time-consuming and prone to mistakes or other discontinuous events.
 
 The purpose of these scripts is to force your series of models into one that is consistent and smoothly-
 varying with the coordinate of interest. 

## Description

 These scripts use the so-called "natural smoothing spline" feature available in the gnuplot package.
a method for taking the constant N out of a constant NVT simulation. That is, you can safely add,
 remove or teleport waters in a simulation that is alredy flying. No need to run LEAP again, or to 
re-heat the system from absolute zero. Instead, you can keep your velocities and sprinkle in a few
more waters into any vacuum bubbles that have appeared. You may also strip out a few 
waters at a time, or even teleport a water from a place where there is
too much electron density to a place that needs more density, using a standard crystallographic
Fo-Fc difference map as a guide. <br>



## Getting Started

### Dependencies

* These are linux c-shell scripts, so you will need a working `/bin/tcsh`
* `gnuplot` to do the actual interpolation

### Installing

* git clone this repo
* copy the two files into somewhere in your shell `$PATH`, and make them executable:
```
    chmod u+x *.com
```
Yes, I know the extension says `*.com`, but these are not Windows executables. The use of `.com` to denote shell scripts pre-dates Windows.

### Executing program

#### one thing at a time
The first program you may want to run is:
```
smooth_pdbs_atom.com refined_00?.pdb atom=123 weight=1 smult=1
```
This will serve as a good test to see if everything is working. Success is indicated by the output of some files called smooth_vs_X.txt, which is the smoothed version of the X-coordinate of atom number 123 in all the files.


There are indeed other ways to do all these things, but this program is needed in the `$PATH` so it can be called by the above scripts.


## Help

Let me know if you have any questions or find any bugs.  In general, debugging information is provided by adding the command-line option: `debug=1`, and possibly `tempfile=temp` to change the default temporary file prefix.
```
rst2pdb_runme.com amber.rst7 debug=1
```

