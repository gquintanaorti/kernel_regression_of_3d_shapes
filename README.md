# Kernel Regression of 3D Shapes

## Authors

* Gregorio Quintana-Orti,
  Depto. de Ingeniería y Ciencia de Computadores,
  Universitat Jaume I,
  12.071 Castellón, Spain.

* Amelia Simó,
  Depto. de Matemáticas,
  Universitat Jaume I,
  12.071 Castellón, Spain.

## Correspondence

Please send correspondence about the code to 
Gregorio Quintana-Ortí: <gquintan@icc.uji.es>

Please send correspondence about the paper to
Amelia Simó: <simo@mat.uji.es>

## License

<!---
(( New 3-clause BSD. ))
(( See file License.txt for more details. ))
-->

Free for non-commercial purposes.

## Disclaimer

This code is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED. 

## Description

This code repository contains several implementations to compute 
kernel regression of 3D shapes in the shape space.
Some related auxiliary code is supplied to read and process datasets.
Besides, more code is given to perform some advanced processing.

We will appreciate feedback from the community on the use of this code.

## Citing this work

We ask those who benefit from this work to cite the following article:

**(( paper or report ))**

<!---

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
-->

## Details of the code

Two variants of the main code to compute the kernel regression
of a 3D shape in the shape space are offered:

* **Basic code**:
  It is written in the R programming language.
  It requires the `mvtnorm` R package.

* **Accelerated code**:
  It is written in both the R and the C programming languages.
  This code uses R code and also some C code to accelerate the most 
  compute-intensive part of the process.
  It requires the `mvtnorm` R package and a high-performance `BLAS`
  (Basic Linear Algebra Subroutines) library, 
  which is usually included in the R programming frameworks.

In both cases, the following R packages can be very useful
to visualize and save the obtained shapes:
`scatterplot3d`, `shapes`, and `rgl`.

Both variants can be found in the `mdl_kernel_regression_in_shape_space.R` file.

The code supplied do not require any installation if you just want 
to use the basic code.
If you want to use the accelerated code, a basic installation is required.
It is described below.

### Details of the basic code:

The header of the main method that implements the kernel regression
is the following one:

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

The header of the main method that implements the kernel regression
is the following one:

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

The C code is embedded by using the `.Call` interface.

To be able to use this accelerated code, you must do the following:

1. Go to the `C_codes` directory.
2. Clean the directory with the `cleanDirectory` command.
3. Compile it with the `c` or the `compile` command.
   These are just scripts that compile the code 
   with the `R CMD SHLIB` command. 
4. Copy the `compute_shape_in_C.so` file (a dynamic library containing the
   compiled code) into the directory containing the R source code.
5. Open the R framework.
6. Load the compiled C code (the dynamic library) 
   with the `dyn.load( "compute_shape_in_C.so" )` command.

### Details of the dataset format:

A method to load (read) a dataset 
from secondary storage (two files) is supplied.
The dataset must be stored in two files:

* `name_coor.csv`: File containing coordinates.
  Every line stores information about one object in the dataset.
  Every line in the file contains the line number and then all the 
  landmarks of the object. 
  For every landmark, the `x`, `y`, and `z` coordinates must be given.
  This file does not require a header line 
  since the number of columns in this file is usually large.
  
* `name_vars.csv`: File containing variables.
  Every line stores information about one object in the dataset.
  Every line in the file contains the line number and then the values 
  of the variables for this object.
  This file requires a header line 
  since the number of columns in this file is usually small. 

Both files must be stored in the `Data` directory.

If those two files are stored in the `Data` directory, 
the data contained in both of them can be loaded into main memory 
with the following command:
```
ds = read_dataset( "name" )
```

Once the data is in main memory, the code supplied can process it.

#### An example of a basic dataset:

As an example, we show two files that describe a very simple
dataset with two basic objects: two tetrahedrons.

The contents of the `tetras1_coor.csv` file is the following:
```
1;1;-0.5774;-0.4082;-1;-0.5774;-0.4082;0;1.1547;-0.4082;0;0;1.2247
2;1;-0.5774;-1.633;-1;-0.5774;-1.633;0;1.1547;-1.633;0;0;4.899
```

The contents of the `tetras1_vars.csv` file is the following
(note the header in the first line):
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

The file with variables contain two sets of variables that define
the two tetrahedrons. 
The variables inside it define the scaling along each axis.

The data contained in both files can be loaded into main memory 
with the following command:
```
ds = read_dataset( "tetras1" )
```

