# randUTV

## Authors

* Gregorio Quintana-Orti,
  Depto. de Ingenieria y Ciencia de Computadores,
  Universitat Jaume I,
  12.071 Castellon, Spain.

* Amelia Simo,
  Depto. de Matematicas,
  Universitat Jaume I,
  12.071 Castellon, Spain.

## Correspondence

Please send correspondence about the code to 
Gregorio Quintana-Ortí: <gquintan@uji.es>

Please send correspondence about the paper to
Amelia Simo: <simo@uji.es>

## License

(( New 3-clause BSD. ))
(( See file License.txt for more details. ))

Free for non-commercial purposes.

## Disclaimer

This code is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED. 

## Description

This code repository contains several implementations to compute 
kernel regression of 3D shapes.
Moreover, some related auxiliary code is supplied to read and process 
datasets.

We will appreciate feedback from the community on the use of this code.

## Citing this work

We ask those who benefit from this work to cite the following article:

```
@ARTICLE{(( To be fixed )),
   author = {xxx },
    title = "{yyy}",
  journal = {ArXiv e-prints},
archivePrefix = "arXiv",
   eprint = {xxx1703.00998},
 primaryClass = "math.NA",
 keywords = {xxx - Mathematics - Numerical Analysis},
     year = 2018,
    month = jul,
   adsurl = {http://adsabs.harvard.edu/abs/2017arXiv170300998M},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Details of the code

We offer two variants of the code to compute the kernel regression:

* Basic code:
  It is written in the R programming language.
  It requires the `mvtnorm` R package.

* Accelerated code:
  It is written in the R programming language and the C programming language.
  It requires the `mvtnorm` R package and a high-performance `BLAS`
  (Basic Linear Algebra Subroutines) library (it is usually included in 
  the R programming frameworks).

In both cases, the `scatterplot3d`, `shapes`, and `rgl` R packages can be very 
useful to visualize the obtained shapes.

Both variants can be found in the `mdl_kernel_regression_in_shape_space.E` file.

### Details of the code writtten in the R programming language:

The code supplied contains the following main method:

```
# =============================================================================
compute_kernel_regression_of_varset = 
    function( dataset, numSteps, varset ) {
#
# Purpose:
#   It computes a regression in the shape space of the set of variables 
#   in "varset" by using the "dataset" received. The iterative process is 
#   limited to a maximum of "numSteps" iterations.
#
# Method's arguments:
#   dataset:   Dataset to be employed in the regression.
#   numSteps:  Maximum number of steps (iterations) in the iterative process.
#   varset:    Set of variables that define the shape to be regressed.
#
```

### Details of the accelerated code:

The code supplied contains the following main method:

```
# =============================================================================
accel_compute_kernel_regression_of_varset = 
    function( dataset, numSteps, varset ) {
#
# Purpose:
#   It computes a regression in the shape space of the set of variables 
#   in "varset" by using the "dataset" received. The iterative process is 
#   limited to a maximum of "numSteps" iterations.
#   This method has been accelerated by performing the most costly parts of the
#   computation in C with BLAS libraries.
# Restriction:
#   To use this method, the C code must be compiled and loaded with "dyn.load".
#
# Method's arguments:
#   dataset:   Dataset to be employed in the regression.
#   numSteps:  Maximum number of steps (iterations) in the iterative process.
#   varset:    Set of variables that define the shape to be regressed.
#
```

### Details of the dataset format

A method to load (read) a dataset from storage is supplied.
The dataset must be stored in two files:

* `name_coor.csv`: File containing coordinates.
  Every line stores information about one object in the dataset.
  Every line in the file contains the line number and then all the 
  landmarks of the object. For every landmark, the x, y, and z coordinates 
  must be stored.
  
* `name_vars.csv`: File containing variables.
  Every line stores information about one object in the dataset.
  Every line in the file contains the line number and then the values 
  of the variables for this object.
  Since the number of columns in this file is usually small, 
  this file requires a header line.

Both files must be stored in the `Data` folder.

As an example, we show the two files that describe a very basic 
dataset con two objects: two tetrahedrons.

The contents of the `tetras1_coor.csv` file is the following:

```
1;1;-0.5774;-0.4082;-1;-0.5774;-0.4082;0;1.1547;-0.4082;0;0;1.2247
2;1;-0.5774;-1.633;-1;-0.5774;-1.633;0;1.1547;-1.633;0;0;4.899
```

The contents of the `tetras1_vars.csv` file is the following:

```
lineNumber;varX;varY;varZ
1;1;1;1
2;1;1;4
```

These two files describe two tetrahedrons: 

* The first one has edge lengths approximately equal to 2 and its base 
is parallel to the XY plane.
It has the following landmarks:
(1, -0.5774, -0.4082),
(-1, -0.5774, -0.4082),
(0, 1.1547, -0.4082), and
(0, 0, 1.2247).

* The second one is the first one with the Z coordinates scaled by 4.
It has the following landmarks:
(1, -0.5774, -1.633),
(-1, -0.5774, -1.633),
(0, 1.1547, -1.633), and
(0, 0, 4.899).


