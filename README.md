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

The license for this code is:
GNU Affero General Public License v3.0.

Read the `LICENSE.txt` file for more details.

The `mdl_my_new_preshape.R` module of this repository
contains some modifications of several routines of the 
noteworthy `shapes` R package by Ian L. Dryden.
These methods have been rewritten to accelerate them 
on medium and large datasets.
License GPL-2 should be applied to this code too.

## Disclaimer

This code is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED. 

## Description

This code repository contains several implementations to compute 
the kernel regression of 3D shapes in the shape space.
Some related auxiliary code is supplied to read and process datasets,
and to perform some advanced processing.

The code provided is written in the R programming language.
Some other code written in the C programming language is provided to
optionally accelerate the most compute-intensive part.

We will appreciate feedback from the community on the use of this code.

## Citing this work

We ask those who benefit from this work to cite the following article:

```
To be updated soon.
Contact the authors if you want to reference this work and 
no paper is shown yet.
```

<!---
**(( paper or report ))**

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

This code repository contains three directories:

* Directory `R_source_codes`:
  It contains the main code.
  You need to download the contents of this directory to be able to work.

* Directory `C_source_codes`:
  It contains some auxiliary code written in the C programming language 
  to accelerate some parts of the R source codes.
  The downloading and installation of this directory is optional.

* Directory `Data`:
  It includes a basic sample dataset called `tetras1`.
  Besides, it contains the datasets employed in the paper.
  It includes two dataset of houses employed in the simulation study:
  Dataset of houses with an error of 0.01, and
  dataset of houses with an error of 0.05.
  These datasets could be different when generated again because of the
  values of the initial seeds for the random number generators.

Two variants of the main code to compute the kernel regression
of a 3D shape in the shape space are offered:

* **Basic code**:
  It is written in the R programming language.
  It requires the `mvtnorm` R package.

* **Accelerated code**:
  It is written in both the R and the C programming languages.
  This code uses R code to perform the basic work, and then 
  some C code to accelerate the most compute-intensive part of the process.
  It requires the `mvtnorm` R package and a high-performance `BLAS`
  (Basic Linear Algebra Subroutines) library, 
  which is usually included in the R programming frameworks.

In both cases, the following R packages can be very useful
to visualize and save the obtained shapes:
`scatterplot3d` and `rgl`.

Both variants can be found 
in the `mdl_kernel_regression_in_shape_space.R` file.

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

The C code is embedded inside the R code by using the `.Call` interface.
It works fine in Linux machines, but we have experienced some problems 
with this interface on some Windows machines.

To be able to use this accelerated code, you must do the following:

1. Go to the `C_source_codes` directory.
2. Clean the directory with the `cleanDirectory` command.
3. Compile it with the `c` command or the `compile` command.
   These are just scripts that compile the code 
   with the `R CMD SHLIB` command from R.
4. Copy the `compute_shape_in_C.so` file (a dynamic library containing the
   compiled code) into the directory containing the R source code.
5. Open the R framework.
6. Load the compiled C code (the dynamic library) 
   with the `dyn.load( "compute_shape_in_C.so" )` command.

### Details of the dataset format:

A method to load (read) a dataset 
from secondary storage is supplied.
The dataset must be partitioned and stored in two files:

* File `name_coor.csv`: 
  This file contains all the landmarks (coordinates) of the objects.
  Every line stores information about one object in the dataset.
  Every line in the file contains the line number and then all the 
  landmarks of the object. 
  For every landmark, the `x`, `y`, and `z` coordinates must be given.
  This file does not require a header line 
  since the number of columns in this file is usually large.
  
* File `name_vars.csv`: 
  This file contains the variables for the objects.
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

As an example, we show the contents of two files 
that describe a very simple dataset with two basic objects: 
two tetrahedrons.
Next we show the contents of the two files:

* File `tetras1_coor.csv`:
  The first column is the line number, and the rest of the columns 
  contain the landmarks of the tetrahedrons.
  The contents is the following:
  ```
  1;1;-0.5774;-0.4082;-1;-0.5774;-0.4082;0;1.1547;-0.4082;0;0;1.2247
  2;1;-0.5774;-1.633;-1;-0.5774;-1.633;0;1.1547;-1.633;0;0;4.899
  ```

* File `tetras1_vars.csv`:
  The first column is the line number, and the rest of the columns 
  contain a set of variables for each object.
  In this case, the variables inside it define the scaling along each axis.
  The contents is the following
  (note the header in the first line):
  ```
  lineNumber;varX;varY;varZ
  1;1;1;1
  2;1;1;4
  ```

These two files combined describe two tetrahedrons: 

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

The data contained in both files can be loaded into main memory 
with the following command:
```
ds = read_dataset( "tetras1" )
```

Once loaded the data in the files into the data structure `ds`,
it can be used as an argument to those methods that require a dataset.

### Summary of the R programming modules:

The code supplied in this repository includes the following modules
with R code:

* `mdl_load_all_modules.R` :
  It removes all the code and 
  then reloads all the modules containing R source code.

* `mdl_read_dataset.R` :
  It reads the two files described above and load the data into a new dataset
  (a new data structure) that can be used by many methods.

* `mdl_dataset_utils.R` :
  It includes several useful functions to process datasets:
  print the dataset dimensions, 
  randomly reorder the dataset,
  reduce the number of objects,
  reduce the resolution (number of landmarks),
  etc.

* `mdl_generate_new_dataset_with_houses.R` :
  It builds and stores a new dataset with houses. 
  The arguments allow to change some parameters.
  The source code can be easily modified to change some other parameters.

* `mdl_kernel_regression_in_shape_space.R` :
  It contains the two main methods to compute kernel regressions of 3D shapes
  in the shape space.
  Furthermore, it contains several auxiliary methods 
  used by these two main methods.

* `mdl_my_new_preshape.R` :
  It contains several modifications of some functions provided by the
  `shapes` R package by Ian L. Dryden. These new implementations 
  have been optimized for performance and are much faster for 
  medium and large datasets than the original codes. 
  For instance, on a dataset with about 3000 landmarks, the new codes
  are about between 100 and 650 times faster, depending on the computer
  and version of the R interpreter.
  This file might not be necessary in the future
  because the `shapes` R package will include the accelerated functions 
  in it in about late summer 2019.

* `mdl_cross_validation.R` :
  It contains code to perform a cross validation.

* `mdl_distances.R` :
  It contains some code to compute Riemannian distances.

* `mdl_tests.R` :
  It contains a test to compare the basic code and the accelerated code.

### A basic example of how to use the code:

Next we show how to use the code to predict the shape of a house:

```
source( "mdl_load_all_modules.R" ) 
ds = read_dataset( "houses_e001" )
varset = c( 10, 10, 15 )
shape = compute_kernel_regression_of_varset( ds, 1000, varset )
plot3d( shape, aspect = F )
```

The lines in the above code sample perform the following:

* The first line loads the codes to run the basic prediction.
* The second line loads the dataset provided with the houses.
* The third line defines a set of variables (varset) for the prediction.
* The fourth line computes the prediction of the previous varset 
  with up to 1000 iterations.
* The fifth line visualizes the predicted object.
  To run this line, 
  you need to install at least `scatterplot3d` and `rgl` R packages.

