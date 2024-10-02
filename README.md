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

 These scripts use the so-called "natural smoothing spline" feature available in the gnuplot package to take a series of related PDB files (say different time points) and create a new series of PDB files with arbitrary
 spacing and xyz coordinates, occupancies, and B factors derived from a smooth curve fit to the original files. It runs in parallel to speed things up.
<br>



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
smooth_pdbs_atom.com refined_00?.pdb atom=123
```
This will serve as a good test to see if everything is working. Success is indicated by the output of some files called smooth_vs_X.txt, which is the smoothed version of the X-coordinate of atom number 123 in all the files. The parallel parent script below uses these text files to make the new PDBs. Other command-line options:
- *.pdb    list of PDB files to smooth over
- atom     ordinal number of atom in first pdb file to extract and smooth its coordinates
- weight   increase data weight in smoothing function, higher numbers make result less smooth
- in_states  specify "x" axis for smoothing as comma-separated list or start-end:step range. default extract from filenames
- out_states specify "x" axis for output files.  default: same as inputs



#### the big run
The first program you may want to run is:
```
smooth_pdbs_multi.com refined_00?.pdb
```
This will take every atom in the first pdb file, extract copies with the same name from all the other files, make a smooth version of all the coordinates, occupancies and B factors, and then create a new series of PDB files containing those smoothed parameters.  If you re-refine these against the same data you will find your results are much more consistent across the series of data sets. Default is to make output files with the same spacings as the inputs, but you can also specify your own input and output coordinate series.
Full list of command-line options is available by running with no arguments, as well as here:
- *.pdb    list of PDB files to smooth over
- ref=     single PDB file containing all possible atoms
- weight   increase data weight in smoothing function, higher numbers make result less smooth
- in_states  specify "x" axis for smoothing as comma-separated list or start-end:step range. default extract from filenames
- out_states specify "x" axis for output files.  default: same as inputs



## Help

Let me know if you have any questions or find any bugs.  In general, debugging information is provided by adding the command-line option: `debug=1`, and possibly `tempfile=temp` to change the default temporary file prefix.
```
rst2pdb_runme.com amber.rst7 debug=1
```

