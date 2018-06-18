# FEHD
Frequency Extracted Hierarchical Decomposition

This folder contains two implementations of FEHD. One is in FORTRAN, and this is for large systems (in terms of components). The other is in MATLAB. The MATLAB version is mainly for exposition. For analysis on data with a lot of components (eg. the number of principal components one chooses to represent their data) then we would highly recommend using the FORTRAN version 


-------------------------------------------------------------------------------------------------------------------
The FORTRAN version of FEHD.

We use gfortran from the GNU compiler collection to compile, as this compiler is freely available for any system. 

Files:
exampleDRIVER.f - The file the user will define parameters and data files
ARspectrum.f - Computes the spectra using the AR model
makeAR.f - Computes AR model from data
pca2.f - Computes the principal components
FEHD.f - The FEHD algorithm
invert.f - inverts a matrix, used to solve the least squares problem that forms the AR model
steepestDescent.f - The minimization routine
compGradient.f - Computes the gradient
grangerIntSINGLE.f - Single thread computation of integral of GC values over
    frequency
calcDet.f - Quickly computes a determinant.

The code also depends on the LAPACK and BLAS libraries. 

On the command line, one can compile as 

gfortran exampleDRIVER.f ARspectrum.f makeAR.f pca2.f FEHD.f invert.f 
steepestDescent.f compGradient.f grangerIntSINGLE.f calcDet.f -L/usr/lib 
-llapack -lblas -fopenmp

(you might have to run it more than once to get all of the modules written).

Then 

./a.out

A file will be created that contains the transformation matrix.

To use the method on your data, open exampleDRIVER.f and edit where indicated.
The data will have channels as columns and time points as rows. If there are
multiple (time-separated) segments (epochs), there are no line breaks or 
anything. The program uses the parameters given by the user to break up
the file into epochs.


----------------------------------------------
The MATLAB implementation

